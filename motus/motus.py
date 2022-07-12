#!/usr/bin/env python

# ============================================================================ #
# motus - a tool for marker gene-based OTU (mOTU) profiling of metagenomes
#
# Authors: Alessio Milanese (milanese@embl.de),
#          Daniel R. Mende (danielrmende@gmail.com),
#          Georg Zeller (zeller@embl.de),
#          Shinichi Sunagawa(ssunagawa@ethz.ch)
#
# Type "motus" for usage help
#
#  LICENSE:
#    motus - a tool for marker gene-based OTU (mOTU) profiling
#    Copyright (C) 2018  A. Milanese, D. R. Mende, G. Zeller & S. Sunagawa
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ============================================================================ #


from __future__ import division
import os
import sys
import argparse
import shlex
import shutil
import time
import subprocess
import glob
import multiprocessing
import tempfile
import errno
import motus.map_genes_to_mOTUs as map_motu
import motus.runBWA as runbwa
import motus.runBWA_for_snv as runbwa_snv
import motus.map_mOTUs_to_LGs as map_lgs
import motus.msamtools_python as filter_sam
import motus.PEfiltering as PEfiltering
import motus.motu_utilities as motu_utilities
import motus.print_CAMI as print_CAMI
import motus.append as append
import motus.downloadDB
import motus.convert_long_reads as convert_long_reads

use_color = True
for i in range(len(sys.argv)):
    arg_this = sys.argv[i]
    if arg_this == "--no_colour":
        use_color = False

if not use_color:
    import motus.UTIL_log as log
else:
    import motus.UTIL_log_col as log







# check if we need to download the database ------------------------------------
for arg in sys.argv:
    if arg == "downloadDB":
        motus.downloadDB.main()
        sys.exit(0)








# position of the script -------------------------------------------------------
path_mOTUs = os.path.realpath(__file__)

path_array = path_mOTUs.split("/")
relative_path = "/".join(path_array[0:-1])

relative_path = relative_path + "/"
motus_binary_folder = "/".join(path_array[0:-1]) + "/"


# check if another database is specified ---------------------------------------
# it will be -db argument
count_db = 0
pos_db_info = 0
DATABASE = ""
for i in range(len(sys.argv)):
    a = sys.argv[i]
    if a == "-db":
        count_db = count_db + 1
        pos_db_info = i+1
# error if the -db appers too many times
if count_db > 1:
    sys.stderr.write("motus error: argument -db repeated multiple times\n")
    sys.exit(1)
# check which DB to use
if count_db == 1:
    if len(sys.argv) > pos_db_info:
        DATABASE = sys.argv[pos_db_info]
        # check if setup.py has been ran already ---------------------------------------
        if not(os.path.isdir(DATABASE)):
            sys.stderr.write("[E::main] Error: Cannot find '-db' directory "+DATABASE+"\n\n")
            sys.exit(1)
    else:
        sys.stderr.write("motus error: argument -db: expected one argument\n")
        sys.exit(1)
else: # default DB
    DATABASE = relative_path+'db_mOTU'


# remove "/" from the end of the database string
def strip_end(text, suffix):
    if not text.endswith(suffix):
        return text
    return strip_end(text[:-1], suffix)
DATABASE = strip_end(DATABASE, "/")
DATABASE_prefix = DATABASE.split("/")[-1]
# we expect all the files in the database to have the same prefix



# check if the database was downloaded already ---------------------------------
if not(os.path.isdir(DATABASE)): # here we check db_mOTU exactly, because want to see if it was installed
    sys.stderr.write("Error: database has not been downloaded. Run 'motus downloadDB' before using the motus profiler\n")
    sys.exit(1)

#-------------------------------------------------------------------------------
# tool version
path_info_version = DATABASE + "/" + DATABASE_prefix + "_versions"
try:
    location = open(path_info_version,'r')
    versions = dict()
    for line in location:
        l = line.rstrip().split('\t')
        if l[0] != "#":
            versions[l[0]] = l[1]
    location.close()
except:
    sys.stderr.write("Error loading file: "+path_info_version+"\nTry to download again the motus profiler\n")
    sys.exit(1)

version_tool = versions["motus"]
git_commit_id = "# git tag version "+version_tool






# ------------------------------------------------------------------------------
#       print the help informations
# ------------------------------------------------------------------------------
class CapitalisedHelpFormatter(argparse.HelpFormatter):
    def add_usage(self, usage, actions, groups, prefix=None):
        if prefix is None:
            prefix = ''
        return super(CapitalisedHelpFormatter, self).add_usage(usage, actions, groups, prefix)


# ------------------------------------------------------------------------------
# function for multiple cores for bwa
def run_bwa_multiple_cores(forward_reads, reverse_reads, single_reads, reference, threads, output, bamOutput, msam_script, technology, verbose, profile_mode, lane_id, result,default_min_len_align_length_map_tax):
    result1 = runbwa.runBWAmapping( forward_reads, reverse_reads, single_reads, reference, threads, output, bamOutput, msam_script, technology, verbose, profile_mode, lane_id,default_min_len_align_length_map_tax,log)
    for i in result1:
        result.append(i)

def prepare_output_bwa(output,outPutIsBam,header_file,all_sam_lines,str_end_header,verbose,str_info_min_len,str_perc_id,str_min_perc_query):

    if not(is_tool("samtools")) and outPutIsBam:
        log.print_warning("samtools is not in the path. Saving .sam format instead of .bam format")
        outPutIsBam = False

    # load the header
    try:
        h_file = open(header_file,'r')
    except:
        log.print_error("error loading file: "+header_file+"\nTry to download again the motus profiler")

    # create the temp sam file with the result
    log.print_message("")
    log.print_message("Create sam file")
    try:
        if output != "" or (output == "" and outPutIsBam):
            sam_temp_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
            os.chmod(sam_temp_file.name, 0o644)
        else:
            sam_temp_file = sys.stdout

        if verbose>4: log.print_message("  \n[map_db] saving intermediate sam file in "+sam_temp_file.name)

        sam_temp_file.write(str_end_header)
        sam_temp_file.write(str_info_min_len)
        sam_temp_file.write(str_perc_id)
        sam_temp_file.write(str_min_perc_query)

        for i in h_file:
            sam_temp_file.write(i)
        for i in all_sam_lines:
            sam_temp_file.write(i)
    except:
        log.print_error("failed to save intermediate sam file")

    # close the sam file
    if output == "" and (not outPutIsBam):
        return 0
        # it went to stdout, so we can close already
    else:
        try:
            sam_temp_file.flush()
            os.fsync(sam_temp_file.fileno())
            sam_temp_file.close()
        except:
            log.print_error("failed to save intermediate sam file")

    #IF WE RETURN SAM
    if not(outPutIsBam):
        # move the temp file to the final destination
        try:
            #os.rename(bam_temp_file.name,args.profile_bam_file) # atomic operation
            shutil.move(sam_temp_file.name,output) #It is not atomic if the files are on different filsystems.
        except:
            if verbose>4: log.print_error("Error when copying intermediate sam to the final destination")
            log.print_error("failed to save the sam file",exit = False)
            log.print_error("You can find the file here:\n"+sam_temp_file.name)
        return 0


    #IF WE RETURN BAM
    log.print_message("Convert sam file to bam file")
    try:
        if output != "":
            bam_temp_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
            os.chmod(bam_temp_file.name, 0o644)
        else:
            bam_temp_file = sys.stdout
        convertCMD = "samtools view -b -Sh "+ sam_temp_file.name
        convert_popenCMD = shlex.split(convertCMD)
        convert_cmd = subprocess.Popen(convert_popenCMD,stdout=bam_temp_file,)

        stdout_s,stderr_s = convert_cmd.communicate()
        if convert_cmd.returncode:
            log.print_error("failed to save intermediate bam file",exit = False)
            log.print_error(stderr_s.decode('ascii'))
    except:
        log.print_error("Error while converting sam to bam file")

    # move the temp file to the final destination
    if output != "":
        try:
            #os.rename(bam_temp_file.name,args.profile_bam_file) # atomic operation
            shutil.move(bam_temp_file.name,output) #It is not atomic if the files are on different filsystems.
        except:
            if verbose>4: log.print_error("Error when copying intermediate bam to the final destination\n")
            log.print_error("failed to save intermediate bam file\nYou can find the file here:\n"+bam_temp_file.name,exit = False)
            # remove the temporary sam file
            os.remove(sam_temp_file.name)
            sys.exit(1) # EXIT

    # remove the temporary sam file
    os.remove(sam_temp_file.name)
    return 0


# ------------------------------------------------------------------------------
# function to check if a specific tool exists
def is_tool(name):
    try:
        devnull = open(os.devnull)
        subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    return True

def is_tool_and_return0(name):
    try:
        devnull = open(os.devnull)
        popenCMD = shlex.split(name)
        child = subprocess.Popen(popenCMD, stdout=devnull)
        streamdata = child.communicate()
        rc = child.wait()
        if rc == 0:
            return True
        else:
            return False
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    return True

# ------------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------------
def main(argv=None):

    motu_call = "python "+(" ".join(sys.argv))

    parser = argparse.ArgumentParser(usage=log.msg(version_tool), formatter_class=CapitalisedHelpFormatter,add_help=False)
    #parser = argparse.ArgumentParser(description='This program calculates mOTU-LGs and specI abundances for one sample', add_help = True)
    parser.add_argument('command', action="store", default=None, help='mode to use the mOTU tool',choices=['profile','map_tax','calc_mgc','calc_motu','merge','map_snv','snv_call','util','downloadDB','prep_long'])
    parser.add_argument('-f', action="store", default=None,dest='forwardReads', help='name of input file for reads in forward orientation, fastq formatted, can be gzipped')
    parser.add_argument('-r', action="store", default=None,dest='reverseReads', help='name of input file for reads in reverse orientation, fastq formatted, can be gzipped')
    parser.add_argument('-s', action="store", default=None,dest='singleReads', help='name of input file for reads without mate, fastq formatted, can be gzipped')
    parser.add_argument('-o', action="store", dest='output', default=None, help='name of output file')
    parser.add_argument('-t', type=int, action="store", dest='threads', default=None, help='Number of threads to be used.')
    parser.add_argument('-e', action='store_true', default=None, dest='onlySpecI', help='Set if you want to profile only specI (mOTU-LGs will go in -1)')
    parser.add_argument('-b', action="store_true", default=None, dest='outPutIsBam', help='Specify if the final output should be a bam formatted file')
    parser.add_argument('-y', action="store", dest='type_output', default=None, help='type of output that you want to print',choices=['base.coverage', 'insert.raw_counts', 'insert.scaled_counts'])
    parser.add_argument('-n', action="store", dest='sampleName', default=None, help='sample name for the current mapping')
    parser.add_argument('-i', action="store", dest='listInputFiles', default=None, help='name of input file(s); sam or bam formatted files. If it is a list: insert all files separated by a comma')
    parser.add_argument('-a', action="store", dest='append_profiles', default=None,help='Environments to add to the profiles.')
    parser.add_argument('-m', action="store", dest='motu_read_counts_file', default=None, help='name of input file; for profiling. It is the file of mOTU read count')
    parser.add_argument('-v', action='store', type=int, default=None, dest='verbose', help='Verbose levels')
    parser.add_argument('-d', action="store", default=None, dest='directory_append', help='directory from where to take the files to append')
    parser.add_argument('-B', action='store_true', default=None, dest='BIOM_output', help='print output in BIOM format')
    parser.add_argument('-C', action='store', default=None, dest='CAMI_output', help='print output in CAMI format',choices=["precision", "recall", "parenthesis"])
    parser.add_argument('-c', action='store_true', default=None, dest='print_rel_ab', help='print result as count instead of relative abundance')
    parser.add_argument('-A', action='store_true', default=None, dest='print_all_taxa', help='print all taxa together')
    parser.add_argument('-u', action='store_true', default=None, dest='print_long_names', help='print names of the species in long format')
    parser.add_argument('-q', action='store_true', default=None, dest='print_full_rank', help='print the full rank of the taxonomy')
    parser.add_argument('-p', action='store_true', default=None, dest='print_NCBI_id', help='print NCBI id')
    parser.add_argument('-K', action='store_true', default=None, dest='keep_metaSNV_dirs', help='keep all the dir produced from metaSNV')
    parser.add_argument('-k', action="store", default=None, dest='taxonomic_level', help='Taxonomic level for the profiling',choices=['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'mOTU'])
    parser.add_argument('-I', action="store", default=None, dest='profile_bam_file', help='name of the bam file to save the intermediate bam file during profiling')
    parser.add_argument('-M', action="store", default=None, dest='profile_mOTU_reads_file', help='name of the output file when profiling, relative to the mOTU read counts')
    parser.add_argument('-g', type=int, action="store", dest='map_lgs_cutoff', default=None, help='cutoff when profiling for the number of genes')
    parser.add_argument('-l', type=int, action="store", dest='min_len_align_length', default=None, help='Minimum length of the alignment')
    #parser.add_argument('-L', type=float, action="store", dest='min_len_align_perc', default=None, help='Minimum length of the alignment, as percentage of the average read length')
    # develope, but to be checked
    parser.add_argument('-CC', type=int, action="store", dest='number_of_cores', default=None, help='Number of cores to use for bwa')
    # db has to be deleted
    #parser.add_argument('-db', action="store", default=None, dest='db', help='database of marker genes',choices=['nr', 'aug.cen', 'cen'])
    parser.add_argument('-db', action="store", default=None, dest='db', help='database of marker genes')
    parser.add_argument('-not_renormalise_cami', action="store_true", default=None, dest='not_renormalise_cami', help='do_not renormalise the cami result')
    parser.add_argument('-save_sam_lines', action="store", dest='save_sam_lines', default=None, help='name of output file')
    parser.add_argument('-load_sam_lines', action="store", dest='load_sam_lines', default=None, help='name of output file')
    parser.add_argument('-min_perc_id', action="store", default=None, dest='min_perc_id', help='minimum percentage of identity when filtering - choose between 97 and 100')
    parser.add_argument('-min_clip_length', action="store", default=None, dest='min_clip_length', help='min. length of alignment when clipped')
    parser.add_argument('-min_perc_align', action="store", default=None, dest='min_perc_align', help='min. percent of the query that must be aligned, between 0 and 100')
    parser.add_argument('-fb', metavar='FLOAT', type=float, default=None,
                        help="Coverage breadth: minimal horizontal genome coverage percentage per sample per species")
    parser.add_argument('-fd', metavar='FLOAT', type=float, default=None,
                        help="Coverage depth: minimal average vertical genome coverage per sample per species")
    parser.add_argument('-fm', metavar='INT', type=int, help="Minimum number of samples per species", default=None)
    parser.add_argument('-fc', metavar='FLOAT', type=float,
                        help="FILTERING STEP II:"
                             "minimum coverage per position per sample per species", default=None)
    parser.add_argument('-fp', metavar='FLOAT', type=float,
                        help="FILTERING STEP II:"
                             "required proportion of informative samples (coverage non-zero) per position",
                        default=None)

    parser.add_argument('-sl', metavar='INT', type=int, help="splitting length for the long reads", default=None)
    parser.add_argument('-ml', metavar='INT', type=int, help="minimum length for the reads", default=None)
    parser.add_argument('-no_gz',action='store_true', default=None, help='do not zip the output')

    parser.add_argument('--version', action='version', version='%(prog)s {0} on python {1}'.format(version_tool, sys.version.split()[0]))
    parser.add_argument('--test', action='store_true', default=None, dest='test_motu', help='test that motus has all the dependencies and is working correctly')
    parser.add_argument('--split_cami_file', action="store", default=None,dest='cami_file_to_split', help='split a gzipped CAMI file into two fzipped for and rev files')
    parser.add_argument('--remove_strain_from_cami_profile', action="store_true", default=None,dest='remove_strain_cami', help='in the cami header, do not print |strain')
    parser.add_argument('--remove_cami_comments', action="store_true", default=None,dest='remove_cami_comments', help='remove comments from the beginning of a cami profile')
    parser.add_argument('--no_colour', action='store_true', default=None, dest='no_colour', help='Do not use colours for the interface')


    args = parser.parse_args()


    # run test.py --------------------------------------------------------------
    if args.test_motu:
        popenCMD = shlex.split("python "+motus_binary_folder+"test.py")
        child = subprocess.Popen(popenCMD)
        child.communicate()
        rc = child.wait()
        return(rc)

    # hidden command: split fastq file from CAMI -------------------------------
    if args.command == 'util':
        if args.cami_file_to_split:
            motu_utilities.split_cami_file(args.cami_file_to_split,args.verbose,log)


    # print menus ----------------------------------------------------------------
    if (args.forwardReads is None) and (args.reverseReads is None) and (args.singleReads is None) and (args.output is None) and (args.threads is None) and (args.onlySpecI is None) and (args.outPutIsBam is None) and (args.keep_metaSNV_dirs is None):
        if (args.type_output is None) and (args.sampleName is None) and (args.listInputFiles is None) and (args.motu_read_counts_file is None) and (args.verbose is None) and (args.directory_append is None) and (args.CAMI_output is None):
            if (args.BIOM_output is None) and (args.taxonomic_level is None) and (args.profile_bam_file is None) and (args.profile_mOTU_reads_file is None) and (args.map_lgs_cutoff is None) and (args.number_of_cores is None):
                if (args.print_rel_ab is None) and (args.print_NCBI_id is None) and (args.print_long_names is None) and (args.fb is None) and (args.fd is None) and (args.fm is None) and (args.fc is None) and (args.fp is None) and (args.print_all_taxa is None):
                    if args.command == 'profile': log.print_menu_profile()
                    if args.command == 'map_snv': log.print_menu_map_snv()
                    if args.command == 'map_tax': log.print_menu_bwa()
                    if args.command == 'calc_mgc': log.print_menu_map_genes()
                    if args.command == 'calc_motu': log.print_menu_map_lgs()
                    if args.command == 'merge': log.print_menu_append()
                    if args.command == 'snv_call': log.print_menu_snv_call()
                    if args.command == 'prep_long': log.print_menu_prep_long()
                    sys.exit(1)
    # set default for args.verbose
    if (args.verbose is None): args.verbose = 3

    # print parameters that are ignored -------------------------------------------
    if args.command == 'profile':
        if ((args.outPutIsBam is not None) and (args.verbose>=2)): log.print_warning("-b ignored")
        if ((args.directory_append is not None) and (args.verbose>=2)): log.print_warning("-d ignored")
        if ((args.fb is not None) and (args.verbose>=2)): log.print_warning("-fb ignored")
        if ((args.fd is not None) and (args.verbose>=2)): log.print_warning("-fd ignored")
        if ((args.fm is not None) and (args.verbose>=2)): log.print_warning("-fm ignored")
        if ((args.fp is not None) and (args.verbose>=2)): log.print_warning("-fp ignored")
        if ((args.fc is not None) and (args.verbose>=2)): log.print_warning("-fc ignored")
        if ((args.keep_metaSNV_dirs is not None) and (args.verbose>=2)): log.print_warning("-K ignored")
    if args.command == 'map_tax':
        if ((args.keep_metaSNV_dirs is not None) and (args.verbose>=2)): log.print_warning("-K ignored")
        if ((args.directory_append is not None) and (args.verbose>=2)): log.print_warning("-d ignored")
        if ((args.onlySpecI is not None) and (args.verbose>=2)): log.print_warning("-e ignored")
        if ((args.type_output is not None) and (args.verbose>=2)): log.print_warning("-y ignored")
        if ((args.sampleName is not None) and (args.verbose>=2)): log.print_warning("-n ignored")
        if ((args.listInputFiles is not None) and (args.verbose>=2)): log.print_warning("-i ignored")
        if ((args.motu_read_counts_file is not None) and (args.verbose>=2)): log.print_warning("-m ignored")
        if ((args.BIOM_output is not None) and (args.verbose>=2)): log.print_warning("-B ignored")
        if ((args.CAMI_output is not None) and (args.verbose>=2)): log.print_warning("-C ignored")
        if ((args.print_rel_ab is not None) and (args.verbose>=2)): log.print_warning("-c ignored")
        if ((args.print_full_rank is not None) and (args.verbose>=2)): log.print_warning("-q ignored")
        if ((args.print_NCBI_id is not None) and (args.verbose>=2)): log.print_warning("-p ignored")
        if ((args.print_long_names is not None) and (args.verbose>=2)): log.print_warning("-u ignored")
        if ((args.taxonomic_level is not None) and (args.verbose>=2)): log.print_warning("-k ignored")
        if ((args.profile_bam_file is not None) and (args.verbose>=2)): log.print_warning("-I ignored")
        if ((args.profile_mOTU_reads_file is not None) and (args.verbose>=2)): log.print_warning("-M ignored")
        if ((args.map_lgs_cutoff is not None) and (args.verbose>=2)): log.print_warning("-g ignored")
        if ((args.fb is not None) and (args.verbose>=2)): log.print_warning("-fb ignored")
        if ((args.fd is not None) and (args.verbose>=2)): log.print_warning("-fd ignored")
        if ((args.fm is not None) and (args.verbose>=2)): log.print_warning("-fm ignored")
        if ((args.fp is not None) and (args.verbose>=2)): log.print_warning("-fp ignored")
        if ((args.fc is not None) and (args.verbose>=2)): log.print_warning("-fc ignored")
    if args.command == 'merge':
        if ((args.keep_metaSNV_dirs is not None) and (args.verbose>=2)): log.print_warning("-K ignored")
        if ((args.CAMI_output is not None) and (args.verbose>=2)): log.print_warning("-C ignored")
        if ((args.onlySpecI is not None) and (args.verbose>=2)): log.print_warning("-e ignored")
        if ((args.type_output is not None) and (args.verbose>=2)): log.print_warning("-y ignored")
        if ((args.sampleName is not None) and (args.verbose>=2)): log.print_warning("-n ignored")
        if ((args.motu_read_counts_file is not None) and (args.verbose>=2)): log.print_warning("-m ignored")
        if ((args.print_rel_ab is not None) and (args.verbose>=2)): log.print_warning("-c ignored")
        if ((args.print_full_rank is not None) and (args.verbose>=2)): log.print_warning("-q ignored")
        if ((args.print_NCBI_id is not None) and (args.verbose>=2)): log.print_warning("-p ignored")
        if ((args.print_long_names is not None) and (args.verbose>=2)): log.print_warning("-u ignored")
        if ((args.taxonomic_level is not None) and (args.verbose>=2)): log.print_warning("-k ignored")
        if ((args.profile_bam_file is not None) and (args.verbose>=2)): log.print_warning("-I ignored")
        if ((args.profile_mOTU_reads_file is not None) and (args.verbose>=2)): log.print_warning("-M ignored")
        if ((args.map_lgs_cutoff is not None) and (args.verbose>=2)): log.print_warning("-g ignored")
        if ((args.forwardReads is not None) and (args.verbose>=2)): log.print_warning("-f ignored")
        if ((args.reverseReads is not None) and (args.verbose>=2)): log.print_warning("-r ignored")
        if ((args.singleReads is not None) and (args.verbose>=2)): log.print_warning("-s ignored")
        if ((args.threads is not None) and (args.verbose>=2)): log.print_warning("-t ignored")
        if ((args.outPutIsBam is not None) and (args.verbose>=2)): log.print_warning("-b ignored")
        if ((args.min_len_align_length is not None) and (args.verbose>=2)): log.print_warning("-l ignored")
        if ((args.fb is not None) and (args.verbose>=2)): log.print_warning("-fb ignored")
        if ((args.fd is not None) and (args.verbose>=2)): log.print_warning("-fd ignored")
        if ((args.fm is not None) and (args.verbose>=2)): log.print_warning("-fm ignored")
        if ((args.fp is not None) and (args.verbose>=2)): log.print_warning("-fp ignored")
        if ((args.fc is not None) and (args.verbose>=2)): log.print_warning("-fc ignored")
    if args.command == 'map_snv':
        if ((args.keep_metaSNV_dirs is not None) and (args.verbose>=2)): log.print_warning("-K ignored")
        if ((args.directory_append is not None) and (args.verbose>=2)): log.print_warning("-d ignored")
        if ((args.onlySpecI is not None) and (args.verbose>=2)): log.print_warning("-e ignored")
        if ((args.type_output is not None) and (args.verbose>=2)): log.print_warning("-y ignored")
        if ((args.sampleName is not None) and (args.verbose>=2)): log.print_warning("-n ignored")
        if ((args.listInputFiles is not None) and (args.verbose>=2)): log.print_warning("-i ignored")
        if ((args.motu_read_counts_file is not None) and (args.verbose>=2)): log.print_warning("-m ignored")
        if ((args.BIOM_output is not None) and (args.verbose>=2)): log.print_warning("-B ignored")
        if ((args.CAMI_output is not None) and (args.verbose>=2)): log.print_warning("-C ignored")
        if ((args.print_rel_ab is not None) and (args.verbose>=2)): log.print_warning("-c ignored")
        if ((args.print_full_rank is not None) and (args.verbose>=2)): log.print_warning("-q ignored")
        if ((args.print_NCBI_id is not None) and (args.verbose>=2)): log.print_warning("-p ignored")
        if ((args.print_long_names is not None) and (args.verbose>=2)): log.print_warning("-u ignored")
        if ((args.taxonomic_level is not None) and (args.verbose>=2)): log.print_warning("-k ignored")
        if ((args.profile_bam_file is not None) and (args.verbose>=2)): log.print_warning("-I ignored")
        if ((args.profile_mOTU_reads_file is not None) and (args.verbose>=2)): log.print_warning("-M ignored")
        if ((args.map_lgs_cutoff is not None) and (args.verbose>=2)): log.print_warning("-g ignored")
        if ((args.outPutIsBam is not None) and (args.verbose>=2)): log.print_warning("-b ignored")
        if ((args.fb is not None) and (args.verbose>=2)): log.print_warning("-fb ignored")
        if ((args.fd is not None) and (args.verbose>=2)): log.print_warning("-fd ignored")
        if ((args.fm is not None) and (args.verbose>=2)): log.print_warning("-fm ignored")
        if ((args.fp is not None) and (args.verbose>=2)): log.print_warning("-fp ignored")
        if ((args.fc is not None) and (args.verbose>=2)): log.print_warning("-fc ignored")
    if args.command == 'calc_mgc':
        if ((args.keep_metaSNV_dirs is not None) and (args.verbose>=2)): log.print_warning("-K ignored")
        if ((args.outPutIsBam is not None) and (args.verbose>=2)): log.print_warning("-b ignored")
        if ((args.directory_append is not None) and (args.verbose>=2)): log.print_warning("-d ignored")
        if ((args.forwardReads is not None) and (args.verbose>=2)): log.print_warning("-f ignored")
        if ((args.reverseReads is not None) and (args.verbose>=2)): log.print_warning("-r ignored")
        if ((args.singleReads is not None) and (args.verbose>=2)): log.print_warning("-s ignored")
        if ((args.threads is not None) and (args.verbose>=2)): log.print_warning("-t ignored")
        if ((args.onlySpecI is not None) and (args.verbose>=2)): log.print_warning("-e ignored")
        if ((args.motu_read_counts_file is not None) and (args.verbose>=2)): log.print_warning("-m ignored")
        if ((args.BIOM_output is not None) and (args.verbose>=2)): log.print_warning("-B ignored")
        if ((args.CAMI_output is not None) and (args.verbose>=2)): log.print_warning("-C ignored")
        if ((args.print_rel_ab is not None) and (args.verbose>=2)): log.print_warning("-c ignored")
        if ((args.print_full_rank is not None) and (args.verbose>=2)): log.print_warning("-q ignored")
        if ((args.print_NCBI_id is not None) and (args.verbose>=2)): log.print_warning("-p ignored")
        if ((args.print_long_names is not None) and (args.verbose>=2)): log.print_warning("-u ignored")
        if ((args.taxonomic_level is not None) and (args.verbose>=2)): log.print_warning("-k ignored")
        if ((args.profile_bam_file is not None) and (args.verbose>=2)): log.print_warning("-I ignored")
        if ((args.profile_mOTU_reads_file is not None) and (args.verbose>=2)): log.print_warning("-M ignored")
        if ((args.map_lgs_cutoff is not None) and (args.verbose>=2)): log.print_warning("-g ignored")
        if ((args.fb is not None) and (args.verbose>=2)): log.print_warning("-fb ignored")
        if ((args.fd is not None) and (args.verbose>=2)): log.print_warning("-fd ignored")
        if ((args.fm is not None) and (args.verbose>=2)): log.print_warning("-fm ignored")
        if ((args.fp is not None) and (args.verbose>=2)): log.print_warning("-fp ignored")
        if ((args.fc is not None) and (args.verbose>=2)): log.print_warning("-fc ignored")
    if args.command == 'calc_motu':
        if ((args.keep_metaSNV_dirs is not None) and (args.verbose>=2)): log.print_warning("-K ignored")
        if ((args.outPutIsBam is not None) and (args.verbose>=2)): log.print_warning("-b ignored")
        if ((args.directory_append is not None) and (args.verbose>=2)): log.print_warning("-d ignored")
        if ((args.profile_bam_file is not None) and (args.verbose>=2)): log.print_warning("-I ignored")
        if ((args.profile_mOTU_reads_file is not None) and (args.verbose>=2)): log.print_warning("-M ignored")
        if ((args.min_len_align_length is not None) and (args.verbose>=2)): log.print_warning("-l ignored")
        if ((args.forwardReads is not None) and (args.verbose>=2)): log.print_warning("-f ignored")
        if ((args.reverseReads is not None) and (args.verbose>=2)): log.print_warning("-r ignored")
        if ((args.singleReads is not None) and (args.verbose>=2)): log.print_warning("-s ignored")
        if ((args.threads is not None) and (args.verbose>=2)): log.print_warning("-t ignored")
        if ((args.type_output is not None) and (args.verbose>=2)): log.print_warning("-y ignored")
        if ((args.motu_read_counts_file is not None) and (args.verbose>=2)): log.print_warning("-m ignored")
        if ((args.fb is not None) and (args.verbose>=2)): log.print_warning("-fb ignored")
        if ((args.fd is not None) and (args.verbose>=2)): log.print_warning("-fd ignored")
        if ((args.fm is not None) and (args.verbose>=2)): log.print_warning("-fm ignored")
        if ((args.fp is not None) and (args.verbose>=2)): log.print_warning("-fp ignored")
        if ((args.fc is not None) and (args.verbose>=2)): log.print_warning("-fc ignored")
    if args.command == 'snv_call':
        if ((args.onlySpecI is not None) and (args.verbose>=2)): log.print_warning("-e ignored")
        if ((args.type_output is not None) and (args.verbose>=2)): log.print_warning("-y ignored")
        if ((args.sampleName is not None) and (args.verbose>=2)): log.print_warning("-n ignored")
        if ((args.motu_read_counts_file is not None) and (args.verbose>=2)): log.print_warning("-m ignored")
        if ((args.print_rel_ab is not None) and (args.verbose>=2)): log.print_warning("-c ignored")
        if ((args.print_full_rank is not None) and (args.verbose>=2)): log.print_warning("-q ignored")
        if ((args.print_NCBI_id is not None) and (args.verbose>=2)): log.print_warning("-p ignored")
        if ((args.print_long_names is not None) and (args.verbose>=2)): log.print_warning("-u ignored")
        if ((args.taxonomic_level is not None) and (args.verbose>=2)): log.print_warning("-k ignored")
        if ((args.profile_bam_file is not None) and (args.verbose>=2)): log.print_warning("-I ignored")
        if ((args.profile_mOTU_reads_file is not None) and (args.verbose>=2)): log.print_warning("-M ignored")
        if ((args.map_lgs_cutoff is not None) and (args.verbose>=2)): log.print_warning("-g ignored")
        if ((args.forwardReads is not None) and (args.verbose>=2)): log.print_warning("-f ignored")
        if ((args.reverseReads is not None) and (args.verbose>=2)): log.print_warning("-r ignored")
        if ((args.singleReads is not None) and (args.verbose>=2)): log.print_warning("-s ignored")
        if ((args.outPutIsBam is not None) and (args.verbose>=2)): log.print_warning("-b ignored")
        if ((args.min_len_align_length is not None) and (args.verbose>=2)): log.print_warning("-l ignored")
        if ((args.BIOM_output is not None) and (args.verbose>=2)): log.print_warning("-B ignored")
        if ((args.CAMI_output is not None) and (args.verbose>=2)): log.print_warning("-C ignored")
        if ((args.listInputFiles is not None) and (args.verbose>=2)): log.print_warning("-i ignored")



    # set the default
    if (args.threads is None): args.threads = 1
    if (args.verbose is None): args.verbose = 3
    if (args.type_output is None): args.type_output = 'insert.scaled_counts'
    if (args.map_lgs_cutoff is None): args.map_lgs_cutoff = 3
    if (args.number_of_cores is None): args.number_of_cores = 1
    if (args.taxonomic_level is None): args.taxonomic_level = "mOTU"
    if (args.sampleName is None): args.sampleName = "unnamed sample"
    if (args.outPutIsBam is None): args.outPutIsBam = False
    if (args.onlySpecI is None): args.onlySpecI = False
    if (args.print_full_rank is None): args.print_full_rank = False
    if (args.BIOM_output is None): args.BIOM_output = False
    if (args.CAMI_output is None): args.CAMI_output = "no_CAMI"
    if (args.print_rel_ab is None): args.print_rel_ab = False
    if (args.print_all_taxa is None): args.print_all_taxa = False
    if (args.print_NCBI_id is None): args.print_NCBI_id = False
    if (args.print_long_names is None): args.print_long_names = False
    if (args.output is None): args.output = ""
    if (args.keep_metaSNV_dirs is None): args.keep_metaSNV_dirs = False

    # default for snv_call
    if (args.fb is None): args.fb = 80.0
    if (args.fd is None): args.fd = 5.0
    if (args.fm is None): args.fm = 2
    if (args.fc is None): args.fc = 5.0
    if (args.fp is None): args.fp = 0.9

    # for cami renormalise
    if (args.not_renormalise_cami is None):
        # if -not_renormalise_cami is not added, then we do normalise - DEFAULT: normalise
        renormalise_cami = True
    else:
        # if -not_renormalise_cami is added, then we do not normalise
        renormalise_cami = False

    # check min length alignment
    default_min_len_align_length = 75 # min length align that we use to filter for profiling
    default_min_len_align_length_map_tax = 75 # min length align that we use in map_tax

    if (args.command == 'map_tax') and (args.min_len_align_length is not None): # map_tax default is 45, unless we set -l
        default_min_len_align_length_map_tax = args.min_len_align_length

    if (args.min_len_align_length is None): args.min_len_align_length = default_min_len_align_length
    if (args.min_len_align_length < 20):
        log.print_error("-l should be greater than 20")
    if (args.min_len_align_length < 45) and (args.verbose>=2):
        log.print_warning("minimum suggested value for -l is 45")
    if (args.min_len_align_length < default_min_len_align_length_map_tax):
        default_min_len_align_length_map_tax = args.min_len_align_length
    if args.command == 'map_snv': # if we use map_snv then we filter with the value that was inserted
        default_min_len_align_length_map_tax = args.min_len_align_length
    info_parameters_p = "-l "+str(args.min_len_align_length)



    # set default for parameters that are not displayed
    if (args.min_perc_id is None): args.min_perc_id = 97
    if (args.min_perc_align is None): args.min_perc_align = 45
    if (args.min_clip_length is None): args.min_clip_length = 60
    args.min_perc_id = float(args.min_perc_id)
    args.min_perc_align = float(args.min_perc_align)
    args.min_clip_length = float(args.min_clip_length)
    # check parameters that are not displayed
    if (args.min_perc_id < 97 or args.min_perc_id >100):
        log.print_error("-min_perc_id should be between 97 and 100")
    if (args.min_perc_align < 45 or args.min_perc_align >100):
        log.print_error("-min_perc_align should be between 45 and 100")
    if (args.min_clip_length < 60):
        log.print_error("-min_clip_length should be greater than 60")

    # ---------------------- check general parameters --------------------------
    if args.verbose < 1:
        log.print_error("verbose level (-v) is less than 1")
    if args.threads < 1:
        log.print_error("number of threads (-t) is less than 1")
    if args.number_of_cores < 1:
        log.print_error("number of cores (-CC) is less than 1")
    if args.map_lgs_cutoff<1 or args.map_lgs_cutoff>10:
        log.print_error("invalid number of marker genes cutoff (-g): "+str(args.map_lgs_cutoff)+" (possible values: from 1 to 10)")

    # check if CAMI output interfer with other parameters ----------------------
    if (args.CAMI_output != "no_CAMI") and args.BIOM_output:
        log.print_error("cannot print result with both -B and -C")
    # check other output options
    if (args.CAMI_output != "no_CAMI") and args.print_all_taxa:
        log.print_error("cannot print result with both -A and -C")
    if args.BIOM_output and args.print_all_taxa:
        log.print_error("cannot print result with both -A and -B")









    # --------------------------------------------------------------------------
    #                          COMMAND TO SPLIT LONG READS
    if args.command == "prep_long":
        # set splitting length and minimum length
        if args.sl is None:
            args.sl = 300
        if args.ml is None:
            args.ml = 50
        # set if gzip the result or not
        gzipped_result = True
        if args.no_gz is not None:
            if args.no_gz:
                gzipped_result = False
        # check that input and output are provided
        if args.listInputFiles is None:
            log.print_error("Missing input file (-i)")
        if args.output is None:
            log.print_error("Missing output file (-o)")

        if args.verbose>2:
            log.print_log("Split long read file into shorter reads")
            log.print_message_execution("The long reads in the file are split into short reads of length "+str(args.sl))

        convert_long_reads.convert_long_reads(args.listInputFiles, args.output, split_len = args.sl, min_len= args.ml, quality = "D", gz_out = gzipped_result, verbose = args.verbose, logger=log)


    # --------------------------------------------------------------------------
    #                                  SNVCALL
    # Begin edit hans 8.3.2018
    # I leave all settings here. They might need to be put to a different location later.
    # --------------------------------------------------------------------------
    if args.command == 'snv_call':
        initial_start_time = time.time()
        if args.verbose>2:
            log.print_log("Calculate SNVs")
            log.print_message_execution("Distances of SNV profiles between samples are calculated with metaSNV")
            log.print_message_execution("(Costea et al. metaSNV: A tool for metagenomic strain level analysis. PLoS One, 2017)\n")

        if args.directory_append is None:
            log.print_error("-d is missing")
        if args.output == "":
            log.print_error("-o is missing")
        # set threads to 1
        args.threads = 1

    # --------------------------------------------------------------------------
    #                       Settings/Params
    # --------------------------------------------------------------------------


        reference = DATABASE+"/"+DATABASE_prefix+"_DB_CEN.fasta"
        annotation = DATABASE+"/"+DATABASE_prefix+"_DB_CEN.fasta.annotations"
        if not os.path.isfile(annotation):
            error_text = "The annotation file "+annotation+" does not exist. Create annotation file with the get.annotations.py command from the metaSNV package"
            log.print_error(error_text)
        bamdir = args.directory_append
        fb = args.fb
        fd = args.fd
        fm = args.fm
        fp = args.fp
        fc = args.fc
        paramsdir = '-m{}-d{}-b{}-c{}-p{}'.format(int(args.fm), int(args.fd), int(args.fb), int(args.fc), float(args.fp))
        # threads = args.threads
        sample_file = os.path.abspath(os.path.join(bamdir, 'sample_list'))
        outdir = os.path.abspath(args.output)
        if os.path.exists(outdir):
            log.print_error("The output directory already exists: "+outdir+" Existing. Delete and restart")
        msnv_log_file = os.path.abspath(os.path.join(outdir, 'metaSNV.log'))
        if (args.verbose>=2): log.print_message("Log file saved in: " + msnv_log_file)
    # --------------------------------------------------------------------------
    #                           create sample_list file
    # --------------------------------------------------------------------------




        if args.verbose>=4:
            log.print_message("[main] Creating sample list file for metaSNV: %s\n" % (sample_file))

        bamfiles = [os.path.abspath(p) for p in glob.glob(os.path.join(bamdir, '*bam'))]

        if len(bamfiles) == 0:
            log.print_error("No BAM files found in the input folder")
        with open(sample_file, 'w') as handle:
            for bamfile in bamfiles:
                handle.write(bamfile + '\n')

    # --------------------------------------------------------------------------
    #                           metaSNV
    # --------------------------------------------------------------------------
        if args.verbose>=4:
            log.print_message("[main] Executing metaSNV main routine\n")


        metaSNV_command = 'metaSNV.py %s %s %s --db_ann %s --threads %d --n_splits %d' % (outdir, sample_file, reference, annotation, args.threads, args.threads)
        if args.verbose>=4:
            log.print_message("[main] Command: %s\n" % metaSNV_command)
        metaSNV_popenCMD = shlex.split(metaSNV_command)
        p = subprocess.Popen(metaSNV_popenCMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, stderr) = p.communicate()
        mlf = open(msnv_log_file, 'w')
        mlf.write(metaSNV_command + '\n')
        mlf.write(stdout.decode('ascii') + '\n')
        mlf.write(stderr.decode('ascii') + '\n')
        if p.returncode:
            log.print_error("metaSNV finished with an error", exit = False)
            log.print_error(str(stderr), exit = False)
            log.print_error(str(stdout), exit = False)
            mlf.close()
            sys.exit(1) # EXIT
    # --------------------------------------------------------------------------
    #                                  FILTER
    # --------------------------------------------------------------------------
        if args.verbose>=4:
            log.print_message("[main] Executing metaSNV filter routine\n")
            # params = '-m 5 -d 10 -b 80 -p 0.9'
        metaSNVFilter_command = 'python '+motus_binary_folder+'metaSNV_Filtering_2.0.py %s -m %d -d %f -b %f -p %f -c %f --n_threads %d' %(outdir, fm, fd, fb, fp, fc, args.threads)
        if args.verbose>=4:
            log.print_message("[main] Command: %s\n" % metaSNVFilter_command)
        metaSNVFilter_popenCMD = shlex.split(metaSNVFilter_command)
        p = subprocess.Popen(metaSNVFilter_popenCMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, stderr) = p.communicate()
        mlf.write(metaSNVFilter_command + '\n')
        mlf.write(stdout.decode('ascii') + '\n')
        mlf.write(stderr.decode('ascii') + '\n')
        if p.returncode:
            log.print_error("metaSNV finished with an error", exit = False)
            log.print_error(str(stderr), exit = False)
            log.print_error(str(stdout), exit = False)
            mlf.close()
            sys.exit(1) # EXIT
    # --------------------------------------------------------------------------
    #                                  REMOVE PADDING
    # --------------------------------------------------------------------------
        if args.verbose>=4:
            log.print_message("[main] Executing metaSNV remove padded routine\n")
        metaSNVRMPadded_command = 'bash '+motus_binary_folder+'motus.remove.padded.sh %s/filtered%s/' % (outdir, paramsdir)
        if args.verbose>=4:
            log.print_message("[main] Command: %s\n" % metaSNVRMPadded_command)
        metaSNVRMPadded_popenCMD = shlex.split(metaSNVRMPadded_command)
        p = subprocess.Popen(metaSNVRMPadded_popenCMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, stderr) = p.communicate()
        mlf.write(metaSNVRMPadded_command + '\n')
        mlf.write(stdout.decode('ascii') + '\n')
        mlf.write(stderr.decode('ascii') + '\n')
        if p.returncode:
            log.print_error("metaSNVRMPadded finished with an error", exit = False)
            log.print_error(str(stderr), exit = False)
            log.print_error(str(stdout), exit = False)
            mlf.close()
            sys.exit(1) # EXIT
    # --------------------------------------------------------------------------
    #                                  DISTANCE
    # --------------------------------------------------------------------------
        if args.verbose>=4:
            log.print_message("[main] Executing metaSNV distance routine\n")
        metaSNVDist_command = 'python '+motus_binary_folder+'metaSNV_DistDiv.py --filt %s/filtered%s --dist --n_threads %d' % (outdir, paramsdir, args.threads)
        if args.verbose>=4:
            log.print_message("[main] Command: %s\n" % metaSNVDist_command)
        metaSNVDist_popenCMD = shlex.split(metaSNVDist_command)
        p = subprocess.Popen(metaSNVDist_popenCMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, stderr) = p.communicate()
        mlf.write(metaSNVDist_command + '\n')
        mlf.write(stdout.decode('ascii') + '\n')
        mlf.write(stderr.decode('ascii') + '\n')
        if p.returncode:
            log.print_error("metaSNVDist finished with an error", exit = False)
            log.print_error(str(stderr), exit = False)
            log.print_error(str(stdout), exit = False)
            mlf.close()
            sys.exit(1) # EXIT
    # --------------------------------------------------------------------------
    #                          KEEP ONLY INTERESTING FILES
    # --------------------------------------------------------------------------
        # remove files
        if not(args.keep_metaSNV_dirs):
            os.remove(os.path.abspath(os.path.join(outdir, 'all_samples')))
            os.remove(os.path.abspath(os.path.join(outdir, 'bed_header')))
            shutil.rmtree(os.path.abspath(os.path.join(outdir, 'bestsplits')))
            shutil.rmtree(os.path.abspath(os.path.join(outdir, 'cov')))
            shutil.rmtree(os.path.abspath(os.path.join(outdir, 'distances')))
            shutil.rmtree(os.path.abspath(os.path.join(outdir, 'filtered')))
            shutil.rmtree(os.path.abspath(os.path.join(outdir, 'snpCaller')))


        # check if the 'distances' dir is empty
        dist_files = [os.path.abspath(p) for p in glob.glob(os.path.join(outdir, 'dist*/*dist'))]
        if len(dist_files) == 0 and args.verbose >= 2:
            log.print_warning("The result is empty")


        # end
        if args.verbose>2:
            log.print_message("")
            log.print_log("Finished computation")
            log.print_log("Total time: "+ str("{0:.2f}".format(time.time() - initial_start_time)) + " s")


        mlf.close()
        sys.exit(0)
    # --------------------------------------------------------------------------
    #                                  END snv_call
    # --------------------------------------------------------------------------











    #set database of sequences
    reference = DATABASE+"/"+DATABASE_prefix+"_DB_NR.fasta" # must keep this
    header_file = DATABASE+"/"+DATABASE_prefix+"_bam_header_NR" # must keep this
    args_db = 'nr' # must keep this

    # if we use snv, then we use centroids
    if args.command == 'map_snv':
        reference = DATABASE+"/"+DATABASE_prefix+"_DB_CEN.fasta"
        header_file = DATABASE+"/"+DATABASE_prefix+"_bam_header_CEN"
        args_db == 'cen'

    # --------------------------------------------------------------------------
    #                      VERSION of DB,scripts,taxonomy
    # --------------------------------------------------------------------------

    # db -----------
    version_db = args_db+versions[args_db]
    if (args.verbose >= 5): log.print_message("[main] Map reads to db "+version_db+"\n")
    # scripts
    version_map_lgs = "calc_motu "+str(versions["map_mOTUs_to_LGs"])
    version_map_lgs = version_map_lgs + " -k "+args.taxonomic_level
    version_map_lgs = version_map_lgs + " -C "+args.CAMI_output
    version_map_lgs = version_map_lgs + " -g "+str(args.map_lgs_cutoff)
    if args.onlySpecI: version_map_lgs = version_map_lgs + " -e"
    if args.print_rel_ab: version_map_lgs = version_map_lgs + " -c"
    if args.print_NCBI_id: version_map_lgs = version_map_lgs + " -p"
    if args.print_long_names: version_map_lgs = version_map_lgs + " -u"
    version_map_lgs = version_map_lgs + " | taxonomy: ref_mOTU_"+str(versions["specI_tax"])
    version_map_lgs = version_map_lgs + " meta_mOTU_"+str(versions["mOTULG_tax"])

    version_append = "# motus version "+str(version_tool)+" | merge "+str(versions["append"])


    # --------------------------------------------------------------------------
    #                                     PROFILE
    # --------------------------------------------------------------------------
    if args.command == 'profile' or args.command == 'map_snv':
        # ---------------- check input -----------------------------------------
        if args.verbose>2:
            initial_start_time = time.time()
            time_after_bwa = time.time() # for now we initialize it, so that if we skip bwa, then we use the initial time
            time_after_map_genes = time.time() # for now we initialize it, so that if we skip calc_mgc, then we use the initial time
            if (args.command == 'profile' and (args.listInputFiles is None and args.motu_read_counts_file is None)) or args.command == 'map_snv':
                log.print_log("Map reads to the marker gene database")
                log.print_message_execution("Reads are aligned (by BWA) to marker gene sequences in the reference database")
                if args.command == 'map_snv':
                    log.print_message_execution("for SNV profiling")

                log.print_message_execution("Reads are filtered for:")
                log.print_message_execution(" - similarity to the reference gene (>97%), and")
                log.print_message_execution(" - length of the alignment (at least "+str(args.min_len_align_length)+" nucleotides, -l option)")

        # check that there is at least one file with reads
        if (args.forwardReads is None) and (args.reverseReads is None) and (args.singleReads is None) and (args.listInputFiles is None) and (args.motu_read_counts_file is None):
            log.print_error("input is missing. Please provide -f,-r or -s or -i or -m")
        # check that for and rev reads are present togehter
        if ((args.forwardReads is not None) and (args.reverseReads is None) and (args.listInputFiles is None) and (args.motu_read_counts_file is None)):
            log.print_error("reverse reads (-r) missing")
        if ((args.forwardReads is None) and (args.reverseReads is not None) and (args.listInputFiles is None) and (args.motu_read_counts_file is None)):
            log.print_error("forward reads (-f) missing")
        # if we have -i bam/sam file, then we ignore the reads input
        if (args.listInputFiles is not None and (args.motu_read_counts_file is None)):
            if (args.verbose >= 2) and (args.forwardReads is not None): log.print_warning("-f FILE ignored since there is -i")
            if (args.verbose >= 2) and (args.reverseReads is not None): log.print_warning("-r FILE ignored since there is -i")
            if (args.verbose >= 2) and (args.singleReads is not None): log.print_warning("-s FILE ignored since there is -i")
        # if we have -m, then we ignore other inputs
        if (args.motu_read_counts_file is not None):
            if (args.verbose >= 2) and (args.forwardReads is not None): log.print_warning("-f FILE ignored since there is -m")
            if (args.verbose >= 2) and (args.reverseReads is not None): log.print_warning("-r FILE ignored since there is -m")
            if (args.verbose >= 2) and (args.singleReads is not None): log.print_warning("-s FILE ignored since there is -m")
            if (args.verbose >= 2) and (args.listInputFiles is not None): log.print_warning("-i FILE ignored since there is -m")

        if (args.listInputFiles is None and args.motu_read_counts_file is None): ## run BWA
            # ---------------------- divide the lanes ------------------------------
            singles = list()
            forw = list()
            reve = list()
            if (args.singleReads is not None): singles = args.singleReads.split(",")
            if (args.forwardReads is not None): forw = args.forwardReads.split(",")
            if (args.reverseReads is not None): reve = args.reverseReads.split(",")

            number_of_lanes = max(len(singles),len(forw),len(reve))
            if args.verbose > 2: log.print_message("\n   Number of detected lanes: "+str(number_of_lanes)+"\n")

            # ----check input: check number of files for forw and rev --------------
            if (args.forwardReads is not None):
                if len(forw) != len(reve):
                    log.print_error("number of files for forward reads ("+str(len(forw))+") is different from number of files for reverse reads ("+str(len(reve))+")")

            # ----check if use multicores and prepere data for multicores ----------
            if args.number_of_cores > 1:
                multiple_cores = True
            else:
                multiple_cores = False

            # check that we have the same number
            # of cores as the number of lanes
            if args.number_of_cores > number_of_lanes:
                if args.verbose >= 2:
                     log.print_warning("We use only "+str(number_of_lanes)+" out of "+str(args.number_of_cores)+" cores")
            if args.number_of_cores < number_of_lanes and args.number_of_cores != 1:
                if args.verbose >= 2:
                     log.print_warning("multiple cores computation is implemented only if the number of cores is equal to the number of lanes ("+str(number_of_lanes)+"). Set -CC "+str(number_of_lanes))
                multiple_cores = False
            #prepare data
            if multiple_cores:
                processes = list()
                results_bwa = list()

            # ------------ execute bwa ---------------------------------------------
            all_sam_lines = list()
            avg_length_reads = list()
            # check total number of reads:
            if args.verbose > 4:
                list_files_check = list()
                if (args.reverseReads is not None):
                    list_files_check = list_files_check + forw
                if (args.singleReads is not None):
                    list_files_check = list_files_check + singles
                motu_utilities.print_n_reads(list_files_check,args.verbose,log)


            for i in range(number_of_lanes):
                if args.verbose>2: log.print_message_time("Run bwa on lane "+str(i+1))

                forward_reads = ""
                reverse_reads = ""
                single_reads = ""
                if (args.reverseReads is not None):
                    if len(forw)> i:
                        forward_reads = forw[i]
                        reverse_reads = reve[i]
                if (args.singleReads is not None):
                    if len(singles)> i:
                        single_reads = singles[i]

                lane_id = "lane"+str(i)
                output = ""
                bamOutput = False
                msam_script = motus_binary_folder+"msamtools_python.py"
                technology = "" # should we implement this?

                profile_mode = True

                # check that the files are fastq and get the average reads length
                if forward_reads != "":
                    avg_length_reads.append(motu_utilities.is_fastq(forward_reads,args.verbose,log))
                if reverse_reads != "":
                    avg_length_reads.append(motu_utilities.is_fastq(reverse_reads,args.verbose,log))
                if single_reads != "":
                    avg_length_reads.append(motu_utilities.is_fastq(single_reads,args.verbose,log))



                if multiple_cores:
                    manager = multiprocessing.Manager()
                    result = manager.list()
                    processes.append(multiprocessing.Process(target=run_bwa_multiple_cores, args=(forward_reads, reverse_reads, single_reads, reference, args.threads, output, bamOutput, msam_script, technology, args.verbose, profile_mode, lane_id, result,default_min_len_align_length_map_tax)) )
                    results_bwa.append(result)
                else:
                    if args.command == 'profile':
                        sam_lines_list = runbwa.runBWAmapping( forward_reads, reverse_reads, single_reads, reference, args.threads, output, bamOutput, msam_script, technology, args.verbose, profile_mode, lane_id,default_min_len_align_length_map_tax,log)
                    else: # if it is map_snv
                        if (args.load_sam_lines is None):
                            sam_lines_list = runbwa_snv.runBWAmapping( forward_reads, reverse_reads, single_reads, reference, args.threads, output, bamOutput, msam_script, technology, args.verbose, profile_mode, lane_id,default_min_len_align_length_map_tax,log)
                        else:
                            sam_lines_list = list()
                    all_sam_lines = all_sam_lines + sam_lines_list

            if multiple_cores:
                for i in range(max(len(singles),len(forw),len(reve))):
                    processes[i].start()

                for i in range(max(len(singles),len(forw),len(reve))):
                    processes[i].join()
                for i in range(max(len(singles),len(forw),len(reve))):
                    all_sam_lines = all_sam_lines + list(results_bwa[i])

            ####### ONLY SNVs   ====================================================
            # if we compute the SNVs, then we have to keep only uniq hits and sort
            if args.command == 'map_snv':

                # quick and dirty save and load all_sam_lines
                if (args.load_sam_lines is not None):
                    temp_f = open(args.load_sam_lines,"r")
                    for jj in temp_f:
                        all_sam_lines.append(jj)
                    temp_f.close()
                if (args.save_sam_lines is not None):
                    temp_f = open(args.save_sam_lines,"w")
                    for jj in all_sam_lines:
                        temp_f.write(jj)
                    temp_f.close()




                start_time_metaSNV = time.time()
                if args.verbose>3:
                    log.print_message("[main] -- After calling bwa there are "+str(len(all_sam_lines))+" sam lines\n")
                if args.verbose>2:
                    log.print_message("")
                    log.print_message("Process the bam file for the snv calling")

                # First: PE filtering
                # coordinates for the padding of the genes are different for centroids
                dictReference2geneLocation = DATABASE+"/"+DATABASE_prefix+"_padding_coordinates_CEN.tsv"
                if args.verbose>3:
                    log.print_message("[main] -- File with gene padding coordinates: "+dictReference2geneLocation)

                all_sam_lines_pe_filtered = PEfiltering.parseBWA_SAMoutput(all_sam_lines, dictReference2geneLocation)
                if args.verbose>3:
                    log.print_message("[main] -- After calling PE filtering there are "+str(len(all_sam_lines_pe_filtered))+" sam lines")

                # Second: keep only uniq mappers
                sam_lines_uniques = dict()
                list_remove = set()
                for sam_line in all_sam_lines_pe_filtered:
                    read_id = sam_line.split('\t')[0]
                    if not(read_id in sam_lines_uniques):
                        sam_lines_uniques[read_id] = sam_line
                    else:
                        sam_lines_uniques[read_id] = "discard"
                        list_remove.add(read_id)

                for id_rem in list_remove:
                    del sam_lines_uniques[id_rem]

                if args.verbose>2:
                    time_after_unique = time.time()
                    log.print_message("Select only reads that map uniquely")

                if args.verbose>3:
                    log.print_message("[main] -- After selecting reads that map uniquely there are "+str(len(all_sam_lines_pe_filtered))+" sam lines\n")



                # Third: samtool sort
                # samtool sort needs a bam file as input
                if args.verbose>4: log.print_message("   [map_snv] saving temporary sam file")
                temp_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
                try:
                    h_file = open(header_file,'r')
                except:
                    log.print_error("Error loading file: "+header_file+"\n [map_db] Try to download again the motus profiler")

                try:
                    for i in h_file:
                        temp_file.write(i)
                    for i in sam_lines_uniques:
                        temp_file.write(sam_lines_uniques[i])
                    temp_file.flush()
                    os.fsync(temp_file.fileno())
                    temp_file.close()
                except:
                    log.print_error("error saving the intermediate sam file, used for sorting")

                if args.verbose>4: log.print_message("   [man_snv] saving intermediate sam file in "+temp_file.name + "\n")

                # take stdout of samtools view and go through sort
                if args.output == "":
                    outfile_bam = sys.stdout # if the output is not specified, then we choose stdout
                else:
                    #outfile_bam = open(args.output, "w")
                    outfile_bam = tempfile.NamedTemporaryFile(delete=False, mode="w")
                    os.chmod(outfile_bam.name, 0o644)



                log.print_message("Sort bam file")

                convertCMD = "samtools view -b -Sh "+temp_file.name + " |  samtools sort - "
                if args.verbose>4: log.print_message("call for sorting: samtools view -b -Sh "+temp_file.name + " |  samtools sort - ")
                try:
                    if is_tool("samtools"):
                        convert_cmd = subprocess.call(convertCMD, stdout=outfile_bam, shell=True) # NOTE: even if samtools fails for "dyld: Library not loaded: @rpathdyld: Library not loaded: @rpath/libcrypto.1.0.0.dylib", it doesnt return an error
                    else:
                        log.print_error("samtools is not in the path")
                except:
                    log.print_error("couldn't save the sorted file")

                # remove temporary file1
                os.remove(temp_file.name) #the file used before sorting
                # move temporary file2 to the final destination
                if args.output != "":
                    try:
                        outfile_bam.flush()
                        os.fsync(outfile_bam.fileno())
                        outfile_bam.close()
                    except:
                        log.print_error("failed to save the final bam/sam file")
                    try:
                        #os.rename(outfile_bam.name,args.output) # atomic operation
                        shutil.move(outfile_bam.name,args.output) #It is not atomic if the files are on different filsystems.
                    except:
                        log.print_error("failed to save the final bam/sam file",exit = False)
                        log.print_error("you can find the file here:\n"+outfile_bam.name)

                # end map_SNV
                if args.verbose>2:
                    log.print_message("")
                    log.print_log("Finished computation")
                    log.print_log("Total time: "+ str("{0:.2f}".format(time.time() - initial_start_time)) + " s")


                return 0 #map_snv is finished, so we can exit




            ####### end ONLY SNVs   ================================================

            # calculate average length of the reads
            reads_avg_length = int(sum(avg_length_reads) / float(len(avg_length_reads)))


            first_script_header = "map_tax "+versions["runBWA"]+" | gene database: "+version_db+" | "+str(reads_avg_length)

            if args.profile_bam_file is not None: # save bam file -----------------------
                if is_tool("samtools"):
                    if args.verbose>2:
                        log.print_message("")
                        log.print_message("Save bam file produced by bwa (provided with -I)")

                    start_time_save_bam = time.time()
                    error_save_intermediate_bam_file = False

                    # load the header
                    try:
                        h_file = open(header_file,'r')
                    except:
                        log.print_error("error loading file: "+header_file+"\nTry to download again the motus profiler")

                    # create the temp sam file with the result
                    try:
                        sam_temp_file = tempfile.NamedTemporaryFile(delete=False, mode="w")

                        if args.verbose>4: log.print_message("  \n[map_db] saving intermediate sam file in "+sam_temp_file.name)

                        str_end_header = "@PG\tID:bwa\tPN:map_tax "+versions["runBWA"]+" | gene database: "+version_db+" | "+str(reads_avg_length)+"\n"
                        sam_temp_file.write(str_end_header)

                        str_info_min_len = "@CO\tmin_len_alignment "+str(default_min_len_align_length_map_tax)+"\n"
                        sam_temp_file.write(str_info_min_len)

                        str_perc_id = "@CO\tmin_perc_id 97\n"
                        sam_temp_file.write(str_perc_id)

                        str_min_perc_query = "@CO\tmin_perc_query 45\n"
                        sam_temp_file.write(str_min_perc_query)

                        for i in h_file:
                            sam_temp_file.write(i)
                        for i in all_sam_lines:
                            sam_temp_file.write(i)
                    except:
                        if args.verbose>1: log.print_warning("failed to save intermediate sam file")
                        if args.verbose>4: log.print_warning("Error when saving the intermediate sam file")
                        error_save_intermediate_bam_file = True

                    # close the sam file
                    try:
                        sam_temp_file.flush()
                        os.fsync(sam_temp_file.fileno())
                        sam_temp_file.close()
                    except:
                        if args.verbose>4: log.print_warning("Error when closing sam file")
                        if not error_save_intermediate_bam_file:
                            if args.verbose>1: log.print_warning("failed to save intermediate sam file")

                    # convert to bam
                    created_bam_file = False
                    try:
                        bam_temp_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
                        os.chmod(bam_temp_file.name, 0o644)
                        convertCMD = "samtools view -b -Sh "+ sam_temp_file.name
                        convert_popenCMD = shlex.split(convertCMD)
                        convert_cmd = subprocess.Popen(convert_popenCMD,stdout=bam_temp_file,)

                        stdout_s,stderr_s = convert_cmd.communicate()
                        if convert_cmd.returncode:
                            if not error_save_intermediate_bam_file:
                                if args.verbose>1: log.print_warning("failed to save intermediate bam file")
                                if args.verbose>1: log.print_warning(stderr_s.decode('ascii'))
                        created_bam_file = True
                    except:
                        created_bam_file = False
                        if args.verbose>4: log.print_warning("Error when converting to bam")
                        if not error_save_intermediate_bam_file:
                            if args.verbose>1: log.print_warning("failed to save intermediate bam file")

                    # move the temp file to the final destination
                    if created_bam_file:
                        try:
                            #os.rename(bam_temp_file.name,args.profile_bam_file) # atomic operation
                            shutil.move(bam_temp_file.name,args.profile_bam_file) #It is not atomic if the files are on different filsystems.
                        except:
                            if args.verbose>4: log.print_warning("Error when copying intermediate bam to the final destination")
                            if not error_save_intermediate_bam_file:
                                if args.verbose>1: log.print_warning("failed to save intermediate bam file")
                                if args.verbose>1: log.print_warning("you can find the file here:\n"+bam_temp_file.name)

                    # remove the temporary sam file
                    os.remove(sam_temp_file.name)


                else:
                    if args.verbose>1: log.print_warning("samtools is not in the path. The intemediate bam file (-I) was not saved")

            if args.verbose>2:
                time_after_bwa = time.time()
                log.print_message("")
                log.print_message("Total time for mapping reads: " + str("{0:.2f}".format(time_after_bwa - initial_start_time))+" s")
                log.print_message("")


#-------# map genes ------------------------------------------------------------
        if args.motu_read_counts_file is None:

            if args.verbose>2:
                log.print_log("Calculate marker gene cluster (MGC) abundance")
                log.print_message_execution("Abundance of MGCs are quantified by the number of inserts aligning to their")
                log.print_message_execution("member genes")

                log.print_message_execution("Inserts can either map to one MGC ('Unique mappers'), or map to many")
                log.print_message_execution("MGCs from different species ('Multiple mappers')\n")

            if (args.listInputFiles is not None):
                if args.verbose>2:
                    log.print_message("Load sam/bam files and inizialize the parametrs for calc_mgc\n")


            # prepare inputs
            listInputFiles = ['unused list']
            database_prefix = DATABASE_prefix
            database_dir = DATABASE
            sample_name = "trial"
            if (args.sampleName is not None): sample_name = args.sampleName
            multThreshold = 3
            winnerThreshold = 0.95
            loserThreshold = 0.01

            output = ""
            msam_script = motus_binary_folder+"msamtools_python.py"
            return_dictionary = True
            profile_mode = True

            if (args.listInputFiles is not None):
                listInputFiles = args.listInputFiles.split(",")
                profile_mode = False
                all_sam_lines_input_map_motu = ""
                #check what was the filter (-l) used during map_tax
                for kk in listInputFiles:
                    filter_l_temp = motu_utilities.read_filter_len_from_bam_file(kk,log)
                    if args.verbose>5 and filter_l_temp is None: log.print_message("filter from bwa not found")
                    if filter_l_temp is not None:
                        if args.verbose>5: log.print_message("Map_tax used -l "+str(filter_l_temp))
                        if args.min_len_align_length < filter_l_temp:
                            if args.verbose>1: log.print_warning("The bam input file was filtered with -l "+str(filter_l_temp)+". If you want to set a lower value for -l, you have to run profile again from the fastq files")



                # we have to check what is the average length of the reads
                all_length_avg = list()
                for kk in listInputFiles:
                    l_temp = motu_utilities.read_length_from_bam_file(kk,log)
                    if l_temp is None:
                        if args.verbose>1: log.print_warning("Cannot find the average read length for the file: "+kk)
                    else:
                        all_length_avg.append(l_temp)
                # if we are inside this if, then reads_avg_length is not present and we create it here
                if len(all_length_avg) != 0:
                    reads_avg_length = int(sum(all_length_avg) / float(len(all_length_avg)))
                else:
                    if args.verbose>1: log.print_warning("no file with information of the length of the reads")
                    reads_avg_length = "unknown" # if we dont know the average length of the file, we set it to 100


            # choose the proper value for min_len_align -------------------
            min_len_align = args.min_len_align_length
            if reads_avg_length != "unknown":
                if reads_avg_length < min_len_align:
                    if args.verbose>1: log.print_warning("Average read length ("+str(reads_avg_length)+") is lower than the -l filter ("+str(min_len_align)+"). We suggest to decrease the value of -l")

            if args.verbose>3: log.print_warning("[main] Minimum alignment length: "+str(min_len_align)+" (average read length: "+str(reads_avg_length)+")")

            # set min clipped length
            if args.verbose>4: log.print_warning("[main] Selecting the clipped length:...")
            #minClippedAlignLength = max(args.min_clip_length,min_len_align)
            minClippedAlignLength = min_len_align
            if args.verbose>4: log.print_warning(str(minClippedAlignLength)+"\n")



            # we have to filter the sam lines if they are not given as input -- note that the filtering (inside the second script) is done only for the one that are loaded, and not during bwa
            if (args.listInputFiles is None):
                all_sam_lines_input_map_motu = filter_sam.run_all_lines((args.min_perc_id/100),min_len_align,args.min_perc_align,all_sam_lines)


            version_information_map_read,mOTU_counts = map_motu.run_mOTUs_v2_mapping(listInputFiles, database_dir, database_prefix, sample_name, multThreshold, winnerThreshold, loserThreshold, minClippedAlignLength, output, msam_script, args.type_output,args.verbose,profile_mode,all_sam_lines_input_map_motu,return_dictionary, args.min_perc_id,min_len_align,args.min_perc_align,log)

            # header for mgc table ----------------
            mgc_table_header = "# "
            if version_information_map_read != "no_info":
                mgc_table_header = mgc_table_header + version_information_map_read + " | "
            else:
                mgc_table_header = mgc_table_header + first_script_header + " | "
            # add info of the parameters
            mgc_table_header = mgc_table_header + "calc_mgc "+versions["map_genes_to_mOTUs"]+" -y "+args.type_output+" "+info_parameters_p


            #save the mOTU read count, actually the mgc table
            if args.profile_mOTU_reads_file is not None:
                start_time_save_motu = time.time()

                try:
                    mgc_temp_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
                    os.chmod(mgc_temp_file.name, 0o644)
                    mgc_temp_file.write(mgc_table_header+"\n")
                    mgc_temp_file.write(sample_name+"\n")
                    for m, value in sorted(mOTU_counts.items()):
                        mgc_temp_file.write("{0}\t{1:.10f}\n".format(m, value))
                except:
                    if args.verbose>1: log.print_warning("failed to save the intermediate mgc table")


                try:
                    mgc_temp_file.flush()
                    os.fsync(mgc_temp_file.fileno())
                    mgc_temp_file.close()
                except:
                    if args.verbose>1: log.print_warning("failed to save the intermediate mgc table\n")
                try:
                    #os.rename(mgc_temp_file.name,args.profile_mOTU_reads_file) # atomic operation
                    shutil.move(mgc_temp_file.name,args.profile_mOTU_reads_file) #It is not atomic if the files are on different filsystems.
                except:
                    if args.verbose>1: log.print_warning("failed to save the intermediate mgc table")
                    if args.verbose>1: log.print_warning("you can find the file here:\n"+mgc_temp_file.name)

                if args.verbose>2:
                    log.print_message("")
                    log.print_message("Save mgc abundance table (provided with -M)")

            if args.verbose>2:
                time_after_map_genes = time.time()
                log.print_message("")
                log.print_message("Total time to calculate the MGCs: " + str("{0:.2f}".format(time_after_map_genes - time_after_bwa))+" s")
                log.print_message("")


#-------# create LGs -----------------------------------------------------------
        if args.verbose>2:
            log.print_log("Generate mOTU profile")
            log.print_message_execution("Each mOTU is composed of 6 to 10 MGCs")
            if args.type_output == 'base.coverage':
                log.print_message_execution("The final mOTU base coverage is obtained as the median of its MGC base coverage values\n")
            else:
                log.print_message_execution("The final mOTU insert count is obtained as the median of its MGC read counts\n")

            log.print_message("At least "+str(args.map_lgs_cutoff)+" (-g) MGCs need to be detected to consider a mOTU as being present\n")

        if args.motu_read_counts_file is not None:# we load the mOTU read count table
            if args.verbose>2: log.print_message("Load the mgc table and compute the profile\n")

        database_dir = DATABASE+"/"
        LGs_map = database_dir+DATABASE_prefix+"_MAP_MGCs_to_mOTUs.tsv"
        LGs_map_l = database_dir+DATABASE_prefix+"_MAP_MGCs_to_mOTUs_in-line.tsv"
        specI_taxonomy = database_dir+DATABASE_prefix+"_taxonomy_ref-mOTUs.tsv"
        mOTULG_taxonomy = database_dir+DATABASE_prefix+"_taxonomy_meta-mOTUs.tsv"
        short_name_file = database_dir+DATABASE_prefix+"_taxonomy_ref-mOTUs_short_names.tsv"
        CAMI_file = database_dir+DATABASE_prefix+"_taxonomy_CAMI.tsv"
        infile=""
        sample_name = "trial"
        if (args.sampleName is not None): sample_name = args.sampleName

        profile_mode = True

        if args.motu_read_counts_file is not None: # we load the mOTU read count table
            profile_mode = False
            mOTU_counts = ""
            mgc_table_header = ""
            infile = args.motu_read_counts_file
            if not os.path.isfile(infile):
                log.print_error(infile+': No such file')

        if args.verbose>2: log.print_message("Create taxonomy profile")

        if args.CAMI_output != "no_CAMI":
            print_CAMI.calculate_abundance(infile, LGs_map, LGs_map_l, args.output, args.map_lgs_cutoff, args.onlySpecI, sample_name, profile_mode,mOTU_counts,mgc_table_header,version_map_lgs,version_tool,args.verbose,motu_call,git_commit_id,version_tool,CAMI_file,args.CAMI_output,renormalise_cami,args.remove_strain_cami,args.remove_cami_comments,log)
        else:
            if args.print_all_taxa:
                map_lgs.calculate_abundance_all(infile, LGs_map, LGs_map_l, specI_taxonomy, mOTULG_taxonomy, args.output, args.map_lgs_cutoff, args.onlySpecI, sample_name, args.taxonomic_level, args.BIOM_output, profile_mode,mOTU_counts, args.print_NCBI_id, args.print_rel_ab,mgc_table_header,version_map_lgs,version_tool,args.verbose,motu_call,git_commit_id,args.print_full_rank,args.print_long_names,short_name_file,version_tool,log)
            else:
                map_lgs.calculate_abundance(infile, LGs_map, LGs_map_l, specI_taxonomy, mOTULG_taxonomy, args.output, args.map_lgs_cutoff, args.onlySpecI, sample_name, args.taxonomic_level, args.BIOM_output, profile_mode,mOTU_counts, args.print_NCBI_id, args.print_rel_ab,mgc_table_header,version_map_lgs,version_tool,args.verbose,motu_call,git_commit_id,args.print_full_rank,args.print_long_names,short_name_file,version_tool,log)



        # end
        if args.verbose>2:
            log.print_message("")
            log.print_log("Finished computation")
            log.print_log("Total time: "+ str("{0:.2f}".format(time.time() - initial_start_time)) + " s")











    # --------------------------------------------------------------------------
    #                                     BWA
    # --------------------------------------------------------------------------
    if args.command == 'map_tax':

        if args.verbose>2:
            log.print_log("Map reads to the marker gene database")
            log.print_message_execution("Reads are aligned (by BWA) to marker gene sequences in the reference database")
            log.print_message_execution("Reads are filtered for:")
            log.print_message_execution(" - similarity to the reference gene (>97%), and")
            log.print_message_execution(" - length of the alignment (at least "+str(args.min_len_align_length)+" nucleotides, -l option)")

        if args.verbose>2:
            initial_start_time = time.time()
            time_after_bwa = time.time() # for now we initialize it, so that if we skip bwa, then we use the initial time
            time_after_map_genes = time.time() # for now we initialize it, so that if we skip calc_mgc, then we use the initial time

        # check that there is at least one file with reads
        if (args.forwardReads is None) and (args.reverseReads is None) and (args.singleReads is None):
            log.print_error("input is missing. Please provide -f,-r or -s")
        # check that for and rev reads are present togehter
        if ((args.forwardReads is not None) and (args.reverseReads is None) and (args.listInputFiles is None) and (args.motu_read_counts_file is None)):
            log.print_error("reverse reads (-r) missing")
        if ((args.forwardReads is None) and (args.reverseReads is not None) and (args.listInputFiles is None) and (args.motu_read_counts_file is None)):
            log.print_error("forward reads (-f) missing")

        # ---------------------- divide the lanes ------------------------------
        singles = list()
        forw = list()
        reve = list()
        if (args.singleReads is not None): singles = args.singleReads.split(",")
        if (args.forwardReads is not None): forw = args.forwardReads.split(",")
        if (args.reverseReads is not None): reve = args.reverseReads.split(",")

        number_of_lanes = max(len(singles),len(forw),len(reve))
        if args.verbose > 2: log.print_message("\n   Number of detected lanes: "+str(number_of_lanes)+"\n")

        # ----check input: check number of files for forw and rev --------------
        if (args.forwardReads is not None):
            if len(forw) != len(reve):
                log.print_error("number of files for forward reads ("+str(len(forw))+") is different from number of files for reverse reads ("+str(len(reve))+")")

        # ----check if use multicores and prepere data for multicores ----------
        if args.number_of_cores > 1:
            multiple_cores = True
        else:
            multiple_cores = False

        # check that we have the same number
        #of cores as the number of lanes
        if args.number_of_cores > number_of_lanes:
            if args.verbose >= 2: log.print_warning("We use only "+str(number_of_lanes)+" out of "+str(args.number_of_cores)+" cores")
        if args.number_of_cores < number_of_lanes and args.number_of_cores != 1:
            if args.verbose >= 2: log.print_warning("multiple cores computation is implemented only if the number of cores is equal to the number of lanes ("+str(number_of_lanes)+"). Set -CC "+str(number_of_lanes))
            multiple_cores = False
        #prepare data
        if multiple_cores:
            processes = list()
            results_bwa = list()

        # ------------ execute bwa ---------------------------------------------
        all_sam_lines = list()
        avg_length_reads = list()
        for i in range(number_of_lanes):
            if args.verbose>2: log.print_message_time("Run bwa on lane "+str(i+1))

            forward_reads = ""
            reverse_reads = ""
            single_reads = ""
            if (args.reverseReads is not None):
                if len(forw)> i:
                    forward_reads = forw[i]
                    reverse_reads = reve[i]
            if (args.singleReads is not None):
                if len(singles)> i:
                    single_reads = singles[i]

            lane_id = "lane"+str(i)
            output = ""
            bamOutput = False
            msam_script = motus_binary_folder+"msamtools_python.py"
            technology = "" # should we implement this?

            profile_mode = True

            # check that the files are fastq and get the average reads length
            if forward_reads != "":
                avg_length_reads.append(motu_utilities.is_fastq(forward_reads,args.verbose,log))
            if reverse_reads != "":
                avg_length_reads.append(motu_utilities.is_fastq(reverse_reads,args.verbose,log))
            if single_reads != "":
                avg_length_reads.append(motu_utilities.is_fastq(single_reads,args.verbose,log))

            if multiple_cores:
                manager = multiprocessing.Manager()
                result = manager.list()
                processes.append(multiprocessing.Process(target=run_bwa_multiple_cores, args=(forward_reads, reverse_reads, single_reads, reference, args.threads, output, bamOutput, msam_script, technology, args.verbose, profile_mode, lane_id, result,default_min_len_align_length_map_tax)) )
                results_bwa.append(result)
            else:
                sam_lines_list = runbwa.runBWAmapping( forward_reads, reverse_reads, single_reads, reference, args.threads, output, bamOutput, msam_script, technology, args.verbose, profile_mode, lane_id,default_min_len_align_length_map_tax,log)
                all_sam_lines = all_sam_lines + sam_lines_list

        if multiple_cores:
            for i in range(max(len(singles),len(forw),len(reve))):
                processes[i].start()

            for i in range(max(len(singles),len(forw),len(reve))):
                processes[i].join()
            for i in range(max(len(singles),len(forw),len(reve))):
                all_sam_lines = all_sam_lines + list(results_bwa[i])



        # print output -------------------------------------------------------

        # calculate average length of the reads
        reads_avg_length = int(sum(avg_length_reads) / float(len(avg_length_reads)))

        str_end_header = "@PG\tID:bwa\tPN:map_tax "+versions["runBWA"]+" | gene database: "+version_db+" | "+str(reads_avg_length)+"\n"

        str_info_min_len = "@CO\tmin_len_alignment "+str(default_min_len_align_length_map_tax)+"\n"
        str_perc_id = "@CO\tmin_perc_id 97\n"
        str_min_perc_query = "@CO\tmin_perc_query 45\n"

        prepare_output_bwa(args.output,args.outPutIsBam,header_file,all_sam_lines,str_end_header,args.verbose,str_info_min_len,str_perc_id,str_min_perc_query)

        # end
        if args.verbose>2:
            log.print_message("")
            log.print_log("Finished computation")
            log.print_log("Total time: "+ str("{0:.2f}".format(time.time() - initial_start_time)) + " s")














    # --------------------------------------------------------------------------
    #                                MAP mOTU
    # --------------------------------------------------------------------------
    if args.command == 'calc_mgc':
        if args.verbose>2:
            log.print_log("Load mapped reads")
            log.print_message("Load SAM lines created with map_tax\n")

        initial_start_time = time.time()

        if (args.listInputFiles is None):
            log.print_error("Please provide at least one input file")

        listInputFiles = args.listInputFiles.split(",")
        database_prefix = DATABASE_prefix
        database_dir = DATABASE
        sample_name = "trial"
        if (args.sampleName is not None): sample_name = args.sampleName
        multThreshold = 3
        winnerThreshold = 0.95
        loserThreshold = 0.01
        msam_script = motus_binary_folder+"msamtools_python.py"
        return_dictionary = True # the result is printed inside the function map_motu.run_mOTUs_v2_mapping
        profile_mode = False
        all_sam_lines = ""
        output = ""

        #check what was the filter (-l) used during map_tax
        for kk in listInputFiles:
            filter_l_temp = motu_utilities.read_filter_len_from_bam_file(kk,log)
            if args.verbose>5 and filter_l_temp is None: log.print_message("filter from bwa not found")
            if filter_l_temp is not None:
                if args.verbose>5: log.print_message("Map_tax used -l "+str(filter_l_temp))
                if args.min_len_align_length < filter_l_temp:
                    if args.verbose>1: log.print_warning("the bam input file was filtered with -l "+str(filter_l_temp)+". If you want to set a lower value for -l, you have to run profile again from the fastq files")


        # we have to check what is the average length of the reads
        all_length_avg = list()
        for kk in listInputFiles:
            l_temp = motu_utilities.read_length_from_bam_file(kk,log)
            if l_temp is None:
                if args.verbose>1: log.print_warning("Cannot find the average read length for the file: "+kk)
            else:
                all_length_avg.append(l_temp)
        # if we are inside this if, then reads_avg_length is not present and we create it here
        if len(all_length_avg) != 0:
            reads_avg_length = int(sum(all_length_avg) / float(len(all_length_avg)))
        else:
            if args.verbose>1: log.print_warning("no file with information of the length of the reads")
            reads_avg_length = "unknown" # if we dont know the average length of the file, we set it to 100
        # choose the proper valur for min_len_align -------------------
        min_len_align = args.min_len_align_length
        if reads_avg_length != "unknown":
            if reads_avg_length < min_len_align:
                if args.verbose>1: log.print_warning("Average read length ("+str(reads_avg_length)+") is lower than the -l filter ("+str(min_len_align)+"). We suggest to decrease the value of -l")

        if args.verbose>2: log.print_message("Minimum alignment length: "+str(min_len_align)+" (average read length: "+str(reads_avg_length)+")\n")

        # set min clipped length
        if args.verbose>4: log.print_message("Selecting the clipped length:...")
        #minClippedAlignLength = max(args.min_clip_length,min_len_align)
        minClippedAlignLength = min_len_align
        if args.verbose>4: log.print_message(str(minClippedAlignLength)+"\n")


        if args.verbose>2:
            log.print_log("Calculate marker gene cluster (MGC) abundance")
            log.print_message_execution("Abundance of MGCs are quantified by the number of inserts aligning to their")
            log.print_message_execution("member genes")

            log.print_message_execution("Inserts can either map to one MGC ('Unique mappers'), or map to many")
            log.print_message_execution("MGCs from different species ('Multiple mappers')\n")

        version_information_map_read,mOTU_counts = map_motu.run_mOTUs_v2_mapping(listInputFiles, database_dir, database_prefix, sample_name, multThreshold, winnerThreshold, loserThreshold, minClippedAlignLength, output, msam_script, args.type_output,args.verbose,profile_mode,all_sam_lines,return_dictionary, args.min_perc_id,min_len_align,args.min_perc_align,log)

        # header for mgc table ----------------
        mgc_table_header = "# "
        if version_information_map_read != "no_info":
            mgc_table_header = mgc_table_header + version_information_map_read + " | "
        # add info of the parameters
        mgc_table_header = mgc_table_header + "calc_mgc "+versions["map_genes_to_mOTUs"]+" -y "+args.type_output+" "+info_parameters_p


        #save the mOTU read count, actually the mgc table
        if args.output != "":
            outfile = tempfile.NamedTemporaryFile(delete=False, mode="w")
            os.chmod(outfile.name, 0o644)
            #outfile = open(args.output, "w")
            if args.verbose>2:
                log.print_message("")
                log.print_message("Save mgc abundance table")
        else:
            outfile = sys.stdout

        try:
            outfile.write(mgc_table_header+"\n")
            outfile.write(sample_name+"\n")
            for i in mOTU_counts:
                outfile.write(i+"\t"+str(mOTU_counts[i])+"\n")
        except:
            log.print_error("failed to save the mgc table")

        if args.output != "":
            try:
                outfile.flush()
                os.fsync(outfile.fileno())
                outfile.close()
            except:
                log.print_error("failed to save the intermediate mgc table")
            try:
                #os.rename(outfile.name,args.output) # atomic operation
                shutil.move(outfile.name,args.output) #It is not atomic if the files are on different filsystems.
            except:
                log.print_error("failed to save the intermediate mgc table\n       You can find the file here:\n"+outfile.name)


        # end
        if args.verbose>2:
            log.print_message("")
            log.print_log("Finished computation")
            log.print_log("Total time: "+ str("{0:.2f}".format(time.time() - initial_start_time)) + " s")










    # --------------------------------------------------------------------------
    #                                 MAP LGs
    # --------------------------------------------------------------------------
    if args.command == 'calc_motu':
        initial_start_time = time.time()
        if args.verbose>2:
            log.print_log("Generate mOTU profile")
            log.print_message_execution("Each mOTU is composed of 6 to 10 MGCs")
            if args.type_output == 'base.coverage':
                log.print_message_execution("The final mOTU base coverage is obtained as the median of its MGC base coverage values\n")
            else:
                log.print_message_execution("The final mOTU insert count is obtained as the median of its MGC read counts\n")

            log.print_message("At least "+str(args.map_lgs_cutoff)+" (-g) MGCs need to be detected to consider a mOTU as being present\n")

        if (args.listInputFiles is None):
            log.print_error("Please provide at least one input file")

        database_dir = DATABASE+"/"
        LGs_map = database_dir+DATABASE_prefix+"_MAP_MGCs_to_mOTUs.tsv"
        LGs_map_l = database_dir+DATABASE_prefix+"_MAP_MGCs_to_mOTUs_in-line.tsv"
        specI_taxonomy = database_dir+DATABASE_prefix+"_taxonomy_ref-mOTUs.tsv"
        mOTULG_taxonomy = database_dir+DATABASE_prefix+"_taxonomy_meta-mOTUs.tsv"
        short_name_file = database_dir+DATABASE_prefix+"_taxonomy_ref-mOTUs_short_names.tsv"
        CAMI_file = database_dir+DATABASE_prefix+"_taxonomy_CAMI.tsv"
        sample_name = "trial"
        if (args.sampleName is not None): sample_name = args.sampleName

        profile_mode = False
        mOTU_counts = ""
        infile = args.listInputFiles
        mgc_table_header=""
        if not os.path.isfile(infile):
            log.print_error(infile+': No such file')

        if args.verbose>2: log.print_message("Create taxonomy profile")

        if args.CAMI_output != "no_CAMI":
            print_CAMI.calculate_abundance(infile, LGs_map, LGs_map_l, args.output, args.map_lgs_cutoff, args.onlySpecI, sample_name, profile_mode,mOTU_counts,mgc_table_header,version_map_lgs,version_tool,args.verbose,motu_call,git_commit_id,version_tool,CAMI_file,args.CAMI_output,renormalise_cami,args.remove_strain_cami,args.remove_cami_comments,log)
        else:
            if args.print_all_taxa:
                map_lgs.calculate_abundance_all(infile, LGs_map, LGs_map_l, specI_taxonomy,mOTULG_taxonomy, args.output, args.map_lgs_cutoff, args.onlySpecI, sample_name, args.taxonomic_level, args.BIOM_output, profile_mode,mOTU_counts, args.print_NCBI_id, args.print_rel_ab,mgc_table_header,version_map_lgs,version_tool,args.verbose,motu_call,git_commit_id,args.print_full_rank,args.print_long_names,short_name_file,version_tool,log)
            else:
                map_lgs.calculate_abundance(infile, LGs_map, LGs_map_l, specI_taxonomy,mOTULG_taxonomy, args.output, args.map_lgs_cutoff, args.onlySpecI, sample_name, args.taxonomic_level, args.BIOM_output, profile_mode,mOTU_counts, args.print_NCBI_id, args.print_rel_ab,mgc_table_header,version_map_lgs,version_tool,args.verbose,motu_call,git_commit_id,args.print_full_rank,args.print_long_names,short_name_file,version_tool,log)


        # end
        if args.verbose>2:
            log.print_message("")
            log.print_log("Finished computation")
            log.print_log("Total time: "+ str("{0:.2f}".format(time.time() - initial_start_time)) + " s")









    # --------------------------------------------------------------------------
    #                                  APPEND
    # --------------------------------------------------------------------------
    if args.command == 'merge':
        initial_start_time = time.time()
        if args.verbose>2:
            log.print_log("Merge taxonomic profiles")
        # check that there is at least one input
        if (args.directory_append is None) and (args.listInputFiles is None):
            log.print_error("input option missing")

        # check that we dont have both inputs
        if (args.directory_append is not None) and (args.listInputFiles is not None):
            log.print_error("too many inputs (-i and -d)")

        if (args.listInputFiles is not None):
            if args.verbose>2:
                log.print_message_execution("Merge multiple mOTU profiles into a taxonomy table, in which rows are mOTUs and")
                log.print_message_execution("columns are samples")
                log.print_message_execution("")

        if args.directory_append is not None:
            if args.verbose>2: log.print_message("Merge all the files in the input directory")
            if (args.directory_append[-1]!="/"):
                args.directory_append = args.directory_append + "/"
        environments_to_merge = []
        public_profiles = DATABASE+"/public_profiles/mOTUs.profiles.gz"
        public_profiles_envo = DATABASE+"/public_profiles/mOTUs.profiles_environments.gz"
        if (args.append_profiles is not None):
            allowed_envo = ['all', 'air', 'bioreactor', 'bee', 'cat', 'cattle', 'chicken', 'dog', 'fish', 'freshwater', 'human', 'marine', 'mouse', 'pig', 'sheep', 'soil', 'termite', 'wastewater']

            for envo in args.append_profiles.split(','):
                if envo not in allowed_envo:
                    log.print_error("Unknown environment " + envo)
                if envo == 'all':
                    environments_to_merge = ['all']
                    break
                else:
                    environments_to_merge.append(envo)


        append.append_profilings(args.directory_append, args.listInputFiles, args.output, args.verbose, args.BIOM_output,version_append,motu_call,version_tool, environments_to_merge, public_profiles, public_profiles_envo, log)

        if args.verbose>2:
            log.print_message("")
            log.print_log("Finished computation")
            log.print_log("Total time: "+ str("{0:.2f}".format(time.time() - initial_start_time))+" s")


    return 0        # success

#-------------------------------- run main -------------------------------------
if __name__ == '__main__':
    status = main()
    sys.exit(status)
