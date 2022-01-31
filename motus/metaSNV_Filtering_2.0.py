#!/usr/bin/env python

from __future__ import division
import os
import sys
import argparse
import glob
import shutil
from multiprocessing import Pool
from functools import partial

basedir = os.path.dirname(os.path.abspath(__file__))

# ======================================================================================================================
# Parse Commandline Arguments


def get_arguments():
    """
    Get commandline arguments and return namespace
    """
    # Initialize Parser
    parser = argparse.ArgumentParser(prog='metaSNV_filtering.py', description='metaSNV filtering step',
                                     epilog='''Note:''', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Not Shown:
    parser.add_argument('--version', action='version', version='%(prog)s 2.0', help=argparse.SUPPRESS)
    parser.add_argument("--debug", action="store_true", help=argparse.SUPPRESS)

    # REQUIRED  arguments:
    parser.add_argument('projdir', help='project name', metavar='Proj')

    # OPTIONAL arguments:
    parser.add_argument('-b', metavar='FLOAT', type=float, default=40.0,
                        help="Coverage breadth: minimal horizontal genome coverage percentage per sample per species")
    parser.add_argument('-d', metavar='FLOAT', type=float, default=5.0,
                        help="Coverage depth: minimal average vertical genome coverage per sample per species")
    parser.add_argument('-m', metavar='INT', type=int, help="Minimum number of samples per species", default=2)
    parser.add_argument('-c', metavar='FLOAT', type=float,
                        help="FILTERING STEP II:"
                             "minimum coverage per position per sample per species", default=5.0)
    parser.add_argument('-p', metavar='FLOAT', type=float,
                        help="FILTERING STEP II:"
                             "required proportion of informative samples (coverage non-zero) per position",
                        default=0.50)
    parser.add_argument('--ind', action='store_true', help="Compute individual SNVs")
    parser.add_argument('--n_threads',metavar=': Number of Processes',default=1,type=int, help="Number of jobs to run simmultaneously.")


    return parser.parse_args()


# ======================================================================================================================
# Basic functions


def file_check():
    """Check if required files exist (True / False)"""
    args.projdir = args.projdir.rstrip('/')
    args.coverage_file = args.projdir + '/' + args.projdir.split('/')[-1] + '.all_cov.tab'
    args.percentage_file = args.projdir + '/' + args.projdir.split('/')[-1] + '.all_perc.tab'
    args.all_samples = args.projdir + '/' + 'all_samples'

    print("Checking for necessary input files...")
    if os.path.isfile(args.coverage_file) and os.path.isfile(args.percentage_file):
        print("found: '{}' \nfound:'{}'".format(args.coverage_file, args.percentage_file))
    else:
        sys.exit(
            "\nERROR: No such file '{}',\nERROR: No such file '{}'".format(args.coverage_file, args.percentage_file))

    if os.path.isfile(args.all_samples):
        print("found: '{}'\n".format(args.all_samples))
    else:
        sys.exit("\nERROR: No such file '{}'".format(args.all_samples))


def print_arguments():
    # Print Defined Thresholds:
    print("Options:")
    if args.b:
        print("threshold: percentage covered (breadth) {}".format(args.b))
    if args.d:
        print("threshold: average coverage (depth) {}".format(args.d))
    if args.m:
        print("threshold: Min. number samples_of_interest per taxid_of_interest {}".format(args.m))
    if args.c:
        print("threshold: Min. position coverage per sample within samples_of_interest {}".format(args.c))
    if args.p:
        print("threshold: Min. proportion of covered samples in samples_of_interest {}".format(args.p))
    if args.ind:
        print("Compute indiv SNVs : {}".format(args.ind))
    if args.n_threads:
        print("Number of parallel processes : {}".format(args.n_threads) )
    print("")


# ======================================================================================================================
# FILTER I: Count SAMPLES_OF_INTEREST per TAXON_OF_INTEREST
# Example: "Which taxon has at least m SoI?"
# Coverage Conditions:
#   Sample of Interest:
#       1. breadth > b (default 40 %)
#       2. depth   >  d (default 5 X)
#   Min. number Samples:
#       3. min samples/TaxID > m (default 2)


def relevant_taxa(args):
    """function that goes through the coverage files and determines taxa and samples of interest"""

    samples_of_interest = {}
    with open(args.coverage_file, 'r') as COV, open(args.percentage_file, 'r') as PER:
        header_cov = COV.readline().split()
        header_per = PER.readline().split()
        COV.readline()  # skip second row in COV_FILE
        PER.readline()  # skip second row in PER_FILE

        # Check if header match:
        if header_cov != header_per:
            sys.exit("ERROR: Coverage file headers do not match!")  # Exit with error message

        # Read taxon by taxon (line by line) check coverage conditions (thresholds)
        for cov, perc in zip(COV, PER):  # Line_cov: TaxID \t cov_valueS1 \t cov_value_S2[...]\n (no white spaces)
            cstring = cov.split()
            pstring = perc.split()
            cov_taxID = cstring.pop(0)
            perc_taxID = pstring.pop(0)
            coverage = list(map(float, cstring))  # convert list_of_strings to list_of_floats (TaxID=STR!)
            percentage = list(map(float, pstring))  # convert list_of_strings to list_of_floats
            sample_count = 0
            sample_names = []

            if cov_taxID != perc_taxID:  # check if taxIDs match!
                sys.exit("ERROR: TaxIDs in the coverage files are not in the same order!")

            for c, p in zip(coverage, percentage):  # two arrays with floats (taxID pop[removed])
                sample_count += 1

                if c >= args.d and p >= args.b:
                    sample_names.append(header_cov[sample_count - 1])  # starting at 0 in array

                if sample_count == len(header_cov) and len(sample_names) >= args.m:
                    samples_of_interest[cov_taxID] = sample_names

    return {'SoI': samples_of_interest, 'h': header_cov}  # return dict()
    COV.close()
    PER.close()


# ======================================================================================================================
# FILTER II: ACQUIRE SNPs WITH SUFFICIENT OCCURRENCE WITHIN SAMPLES_OF_INTEREST
#   SNP Conditions (Default):
#       1. Position covered by at least c (5) reads
#       2. Position present in at least proportion p (50 %) of the accepted samples_of_interest

def filter_two(species, args, snp_files, outdir):
    """position wise filtering"""

    snp_taxID = '_'

    # read <all_samples> for snp_file header:
    all_samples = open(args.all_samples, 'r')
    snp_header = all_samples.read().splitlines()
    snp_header = [i.split('/')[-1] for i in snp_header]  # get name /trim/off/path/to/sample.name.bam

    for best_split_x in snp_files:
        with open(best_split_x, 'r') as file:
            for snp_line in file:  # position wise loop
                snp_taxID = snp_line.split()[0].split('.')[0]  # Name of Genome change from . to ]

                # Species filter:
                if snp_taxID != species:  # Check if Genome is of interest
                    continue  # Taxon is not relevant, NEXT!
                else:
                    # Sample filter: only load samples with enough coverage
                    sample_list = samples_of_interest[snp_taxID]  # Load Sample List
                    # !!! Sample order, get indices - INDICES based on ORDER in COV/PERC file !!!
                    sample_indices = []
                    for name in sample_list:
                        sample_indices.append(snp_header.index(name))

                    # Position filter:
                    # Positions with sufficient coverage (c) and proportion (p) in samples of interests (SoIs).
                    whole_line = snp_line.split()
                    site_coverage = list(map(int, whole_line[4].split('|')))  # Site coverages as list of ints
                    nr_good = 0
                    for index in sample_indices:
                        if site_coverage[index] < args.c or site_coverage[index] == 0:
                            continue  # Not enough coverage depth at position in sample, no increment
                        else:
                            nr_good += 1

                    # Filter: Position incidence with sufficient coverage:
                    if float(nr_good) / len(sample_indices) < args.p:  # if snp_incidence < x % drop SNP
                        continue  # mainly uninformative position, SNP incidence < x %, SNP dropped!

                    else:
                        # Calculate SNP allele frequency:
                        # If file do not exist yet, create it:
                        if not os.path.isfile(outdir + '/' + '%s.filtered.freq' % snp_taxID):
                            outfile = open(outdir + '/' + '%s.filtered.freq' % snp_taxID, 'w')
                            print("Generating: {}".format(outdir + '/' + '%s.filtered.freq' % snp_taxID))
                            outfile.write('\t' + "\t".join(sample_list) + '\n')

                        # Loop through alternative alleles [5](comma separated):
                        line_id = ":".join(whole_line[:4])  # line_id composed of CHROM:REFGENE:POS:REFBASE
                        reference_base = whole_line[3]
                        alt_bases_totalcov = []  # total coverage per alt allele
                        # VCF format

                        # Loop Start:
                        for snp in whole_line[5].split(','):

                            xS = snp.split('|')
                            snp_coverage = list(map(float, xS[3:]))  # coverage string (unfiltered!)

                            # Sanity check:
                            if len(site_coverage) != len(snp_coverage):
                                print("ERROR: SNP FILE {} is corrupted".format(best_split_x))
                                sys.exit("ERROR: Site coverage and SNP coverage string have uneven length!")

                            # Compute allele frequency tables:
                            alt_base = snp.split('|')[1]  # alternative base

                            # Frequency Computation
                            total_reads = 0
                            snp_frq = []  # frequencies for SNPs (pos > cX in at least p% of the SoIs)
                            for index in sample_indices:

                                # prevent division by zero!
                                if site_coverage[index] >= args.c and site_coverage[index] != 0:
                                    snp_frq.append(snp_coverage[index] / site_coverage[index])
                                else:
                                    snp_frq.append(-1)

                            # Write Output Allele Frequencies (Default)
                            outfile.write(
                                ":".join(snp_line.split()[:4]) + '>' + alt_base + ':' + xS[2] + '\t' + "\t".join(
                                    str(x) for x in snp_frq) + '\n')
    if 'outfile' in locals():
        print("closing: {}".format(species))
        outfile.close()


# ======================================================================================================================
# Script

if __name__ == "__main__":

    args = get_arguments()
    print_arguments()
    file_check()

    # ==========================================
    # Filtering I - Determine Taxa of Interest:
    # ==========================================

    samples_of_interest = relevant_taxa(args)['SoI']
    header_cov = relevant_taxa(args)['h']

    print(samples_of_interest.keys())

    # =========================================
    # Filtering II - Position wise filtering
    # =========================================

    pars_toprint = '-m{}-d{}-b{}-c{}-p{}'.format(int(args.m), int(args.d), int(args.b), int(args.c), float(args.p))
    filt_folder = args.projdir + '/filtered' + pars_toprint + '/'

    if not os.path.exists(filt_folder):
        os.makedirs(filt_folder)
        os.makedirs(filt_folder + '/pop/')
    else:
        shutil.rmtree(filt_folder)
        os.makedirs(filt_folder)
        os.makedirs(filt_folder + '/pop/')

    p = Pool(processes=args.n_threads)
    partial_Div = partial(filter_two,
                          args=args,
                          snp_files=glob.glob(args.projdir + '/snpCaller/called*'),
                          outdir=filt_folder + '/pop')
    p.map(partial_Div, samples_of_interest.keys())
    p.close()
    p.join()

    if args.ind:
        if not os.path.exists(filt_folder + '/ind/'):
            os.makedirs(filt_folder + '/ind/')
        p = Pool(processes=args.n_threads)
        partial_Div = partial(filter_two,
                              args=args,
                              snp_files=glob.glob(args.projdir + '/snpCaller/indiv*'),
                              outdir=filt_folder + '/ind')
        p.map(partial_Div, samples_of_interest.keys())
        p.close()
        p.join()
