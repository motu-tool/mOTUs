#!/usr/bin/env python

from __future__ import division
import os
import sys
import argparse
import glob
from multiprocessing import Pool
from functools import partial

try:
    import numpy as np
except ImportError:
    sys.stderr.write("Numpy is necessary to run this script.\n")
    sys.exit(1)
try:
    import pandas as pd
except ImportError:
    sys.stderr.write("Pandas is necessary to run this script.\n")
    sys.exit(1)


basedir = os.path.dirname(os.path.abspath(__file__))


############################################################
### Parse Commandline Arguments

def get_arguments():
    '''
    Get commandline arguments and return namespace
    '''
    ## Initialize Parser
    parser = argparse.ArgumentParser(prog='metaSNV_DistDiv.py', description='metaSNV distance and diversity computation', epilog='''Note:''', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Not Showns:
    parser.add_argument('--version', action='version', version='%(prog)s 2.0', help=argparse.SUPPRESS)
    parser.add_argument("--debug", action="store_true", help=argparse.SUPPRESS)

    # REQUIRED  arguments:
    parser.add_argument('--filt', metavar=': Filtered frequency files', help="Folder containing /*.filtered.freq", required=True)

    # OPTIONAL  arguments:
    parser.add_argument('--dist',action='store_true', help="Compute distances")
    parser.add_argument('--div',action='store_true', help="Compute Diversity and FST")
    parser.add_argument('--divNS',action='store_true', help="Computing piN and piS")
    parser.add_argument('--matched',action='store_true', help="Computing on matched positions only")
    parser.add_argument('--n_threads',metavar=': Number of Processes',default=1,type=int, help="Number of jobs to run simmultaneously.")


    return parser.parse_args()

############################################################
### Basic functions

def file_check():
    ''' Check if required files exist (True / False)'''

    args.projdir = '/'.join(args.filt.rstrip('/').split('/')[:-1])
    args.pars = args.filt.rstrip('/').split('/')[-1].strip('filtered')
    args.coverage_file = args.projdir+'/'+args.projdir.split('/')[-1]+'.all_cov.tab'
    args.percentage_file = args.projdir+'/'+args.projdir.split('/')[-1]+'.all_perc.tab'
    args.bedfile = args.projdir+'/'+'bed_header'

    print("Checking for necessary input files...")
    if os.path.isfile(args.coverage_file) and os.path.isfile(args.percentage_file) and os.path.isfile(args.bedfile):
        print("found: '{}' \nfound:'{}' \nfound:'{}'".format(args.coverage_file, args.percentage_file, args.bedfile))
    else:
        sys.exit("\nERROR: No such file '{}',\nERROR: No such file '{}',\nERROR: No such file '{}'".format(args.coverage_file, args.percentage_file, args.bedfile))


def print_arguments():
    ## Print Defined Arguments:
    print("")

    print("Checking required arguments:")
    if args.filt:
        print("Filtered folder : {}".format(args.filt) )
    print("Checking optional arguments:")
    if args.dist:
        print("Computing distances : {}".format(args.dist) )
    if args.div:
        print("Computing diversity and FST : {}".format(args.div) )
    if args.divNS:
        print("Computing N and S diversities : {}".format(args.divNS) )
    if args.matched:
        print("Matching positions (present in 90% of the samples) : {}".format(args.matched) )
    if args.n_threads:
        print("Number of parallel processes : {}".format(args.n_threads) )
    print("")



############################################################
### Distances

def l1nonans(d1,d2):
    return np.abs(d1 - d2).mean()

def alleledist(d1,d2, threshold=.6):
    return (np.abs(d1 - d2) > threshold).mean()

def computeDist(filt_file):
    ''' Compute distances per species '''
    species = filt_file.split('/')[-1].replace('.freq','')
    data = pd.read_table(filt_file, index_col=0, na_values=['-1']).T

    dist = [[l1nonans(data.iloc[i], data.iloc[j]) for i in range(len(data))] for j in range(len(data))]
    dist = pd.DataFrame(dist, index=data.index, columns=data.index)
    dist.to_csv(outdir + '/' + '%s.mann.dist' % species, sep='\t')

    dist = [[alleledist(data.iloc[i], data.iloc[j]) for i in range(len(data))] for j in range(len(data))]
    dist = pd.DataFrame(dist, index=data.index, columns=data.index)
    dist.to_csv(outdir + '/' + '%s.allele.dist' %species, sep='\t')

def computeAllDist(args):

    print("Computing distances")

    allFreq = glob.glob(args.filt + '/*.freq')

    p = Pool(processes=args.n_threads)
    p.map(computeDist, allFreq)
    p.close()
    p.join()



############################################################
### Pairwise Diversity

def compute_diversity(sample1, sample2):

    '''Pairwise nucleotide diversity'''

    sample1nd = sample1.reset_index().drop_duplicates(subset='index', keep=False).set_index('index')
    sample2nd = sample2.reset_index().drop_duplicates(subset='index', keep=False).set_index('index')
    sample2nd = sample2nd.reindex(index=sample1nd.index)
    s1 = sample1nd.values
    s2 = sample2nd.values
    valid = ~(np.isnan(s1) | np.isnan(s2))
    s1 = s1[valid]
    s2 = s2[valid]
    s1 = np.vstack([s1, 1 - s1])
    s2 = np.vstack([s2, 1 - s2])
    dist_nd = (s1[0]*s2[1]+s1[1]*s2[0]).sum()

    def position_diversity(x):
        out = np.outer(x.s1.values, x.s2.values)
        return np.nansum(out) - np.nansum(out.diagonal())

    sample1d = sample1.ix[sample1.index[sample1.index.duplicated()]]
    sample2d = sample2.ix[sample2.index[sample2.index.duplicated()]]

    if not len(sample1d) or not len(sample2d):
        # No duplicates
        return dist_nd

    both = pd.DataFrame({'s1': sample1d, 's2': sample2d})
    both = both.reset_index()
    both = pd.concat([both,(1. - both.groupby('index').sum()).reset_index()])
    dist_d = both.groupby('index', group_keys=False).apply(position_diversity).sum()

    return dist_d + dist_nd

############################################################
### Per Species Diversity

def computeDiv(filt_file, horizontal_coverage, vertical_coverage, bedfile_tab, matched):
    '''Per species computation'''

    species = int(filt_file.split('/')[-1].split('.')[0])
    data = pd.read_table(filt_file, index_col=0, na_values=['-1'])

    pre_index = [i.split(':') for i in list(data.index)]
    # Setting index for each position
    data = data.set_index(pd.Index([item[0] + ':' + item[1] + ':' + item[2] for item in pre_index]))
    data = data.sort_index()

    ########
    ## If matched, filter for 'common' positions :
    if matched:
        def filt_proportion(x):
            if len(x) == 2:
                x = x.iloc[1]
            n = np.count_nonzero(np.isnan(x))
            return n>(len(x)*(0.1))

        index_drop = [index for index in data.index if filt_proportion(data.loc[index])]
        data = data.drop(index_drop)

    ########
    ## Number of bases observed :
    genome_length = bedfile_tab.loc[str(species), 2].sum()
    # Genome length corrected for horizontal coverage
    correction_coverage = [[(min(horizontal_coverage.loc[species, i], horizontal_coverage.loc[species, j]) * genome_length) / 100 for i in data.columns] for j in data.columns]

    ########
    ## Vertical coverage in pi within : AvgCov / (AvgCov - 1)
    for i in data.columns:
        j = list(data.columns).index(i)
        correction_within = vertical_coverage.loc[species,i] / (vertical_coverage.loc[species,i] - 1)
        correction_coverage[j][j] = correction_coverage[j][j] / correction_within

    div = [[compute_diversity(data.iloc[:, i], data.iloc[:, j]) / correction_coverage[j][i] for i in range(j + 1)] for j in range(len(data.columns))]
    FST = [[(1-(div[i][i]+div[j][j])/(2*div[j][i])) for i in range(j + 1)] for j in range(len(div))]

    div = pd.DataFrame(div, index=data.columns, columns=data.columns)

    FST = pd.DataFrame(FST, index=data.columns, columns=data.columns)

    div.to_csv(outdir + '/' + '%s.diversity' % species, sep='\t')
    FST.to_csv(outdir + '/' + '%s.FST' % species, sep='\t')


############################################################
### Per Species N & S Diversity

def computeDivNS(filt_file, horizontal_coverage, vertical_coverage, bedfile_tab, matched):
    '''Per species computation'''

    species = int(filt_file.split('/')[-1].split('.')[0])
    data = pd.read_table(filt_file, index_col=0, na_values=['-1'])

    pre_index = [i.split(':') for i in list(data.index)]
    # Setting index for each position
    index1 = [item[0] + ':' + item[1] + ':' + item[2] for item in pre_index]
    # Setting index for Non-synonymous vs Synonymous
    index2 = [item[4].split('[')[0] for item in pre_index]
    data = data.set_index(pd.MultiIndex.from_arrays([index1, index2], names=['index', 'significance']))
    data = data.sort_index()

    data_N = data.xs('N', level='significance')
    data_S = data.xs('S', level='significance')

    ########
    ## If matched, filter for 'common' positions :
    if matched:
        def filt_proportion(x):
            if len(x) == 2:
                x = x.iloc[1]
            n = np.count_nonzero(np.isnan(x))
            return n>(len(x)*(0.1))

        index_drop = [index for index in data_N.index if filt_proportion(data_N.loc[index])]
        data_N = data_N.drop(index_drop)

        index_drop = [index for index in data_S.index if filt_proportion(data_S.loc[index])]
        data_S = data_S.drop(index_drop)

    ########
    ## Number of bases observed :
    genome_length = bedfile_tab.loc[str(species), 2].sum()
    # Genome length corrected for horizontal coverage
    correction_coverage = [[(min(horizontal_coverage.loc[species, i], horizontal_coverage.loc[species, j]) * genome_length) / 100 for i in data.columns] for j in data.columns]

    ########
    ## Vertical coverage in pi within : AvgCov / (AvgCov - 1)
    for i in data.columns:
        j = list(data.columns).index(i)
        correction_within = vertical_coverage.loc[species,i] / (vertical_coverage.loc[species,i] - 1)
        correction_coverage[j][j] = correction_coverage[j][j] / correction_within

    div_N = [[compute_diversity(data_N.iloc[:, i], data_N.iloc[:, j]) / correction_coverage[j][i] for i in range(j + 1)] for j in range(len(data_N.columns))]
    div_N = pd.DataFrame(div_N, index=data_N.columns, columns=data_N.columns)
    div_N.to_csv(outdir + '/' + '%s.N_diversity' % species, sep='\t')

    div_S = [[compute_diversity(data_S.iloc[:, i], data_S.iloc[:, j]) / correction_coverage[j][i] for i in range(j + 1)] for j in range(len(data_S.columns))]
    div_S = pd.DataFrame(div_S, index=data_S.columns, columns=data_S.columns)
    div_S.to_csv(outdir + '/' + '%s.S_diversity' % species, sep='\t')




############################################################
### Compute Diversity for all Species

def computeAllDiv(args):

    '''Computing diversities & FST'''

    print("Computing diversities & FST")

    # Load external info : Coverage, genomes size, genes size
    horizontal_coverage = pd.read_table(args.percentage_file, skiprows=[1], index_col=0)
    vertical_coverage = pd.read_table(args.coverage_file, skiprows=[1], index_col=0)

    bedfile_tab = pd.read_table(args.bedfile, index_col=0, header=None)
    bed_index = [i.split('.')[0] for i in list(bedfile_tab.index)]
    bedfile_tab = bedfile_tab.set_index(pd.Index(bed_index))

    #All filtered.freq files in input folder
    allFreq = glob.glob(args.filt + '/*.freq')

    if args.div:
        p = Pool(processes=args.n_threads)
        partial_Div = partial(computeDiv, horizontal_coverage=horizontal_coverage, vertical_coverage=vertical_coverage, bedfile_tab=bedfile_tab, matched=args.matched)
        p.map(partial_Div, allFreq)
        p.close()
        p.join()

    if args.divNS:
        p = Pool(processes=args.n_threads)
        partial_DivNS = partial(computeDivNS, horizontal_coverage=horizontal_coverage, vertical_coverage=vertical_coverage, bedfile_tab=bedfile_tab, matched=args.matched)
        p.map(partial_DivNS, allFreq)
        p.close()
        p.join()


############################################################
### Script

if __name__ == "__main__":

    #####################

    args = get_arguments()
    print_arguments()
    file_check()

    if args.debug:
        print_arguments()

    #####################

    if args.matched:
        outdir = args.projdir + '/distances' + args.pars + '.matched_pos/'
    else:
        outdir = args.projdir + '/distances' + args.pars +'/'

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    #####################

    if args.dist:
        computeAllDist(args)

    if args.div or args.divNS:
        computeAllDiv(args)
