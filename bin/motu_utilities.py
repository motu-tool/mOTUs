#!/usr/bin/env python

# ============================================================================ #
# motu_utilities.py: some useful functions
#
# ============================================================================ #

from itertools import islice
from gzip import GzipFile
from bz2 import BZ2File
import sys
import os
import shlex
import subprocess

# ------------------------------------------------------------------------------
# function to check if a specific tool exists
# ------------------------------------------------------------------------------
def is_tool(name):
    try:
        devnull = open(os.devnull)
        subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            return False
    return True

# ------------------------------------------------------------------------------
# function to load a sam/bam file
# ------------------------------------------------------------------------------
def readSAMfile(strSAMfilename):

    try:
        from subprocess import DEVNULL
    except ImportError:
        DEVNULL = open(os.devnull, 'wb')

    if not os.path.isfile(strSAMfilename):
        sys.stderr.write("[E::main] Error: "+strSAMfilename+': No such file.\n')
        sys.exit(1)

    #we can understand if it is SAM or BAM reading the first line:
    # in the new version of samtools it is detected automatically if we have a bam or sam file
    isBam = False
    try:
        location = open(strSAMfilename,'r')
    except:
        sys.stderr.write("[E::main] Error: failed to load "+strSAMfilename+'\n')
        sys.exit(1)
    try:
        first_line = location.readline()
    except UnicodeDecodeError:
        # is bam
        isBam = True
    location.close()
    # here we keep it for compatibility with old versions of samtools
    samtoolsCMD = "samtools view -Sh -F 4 " + strSAMfilename
    if (isBam):
        samtoolsCMD = "samtools view -h -F 4 " + strSAMfilename

    #sys.stderr.write(samtoolsCMD+"\n")
    samtools_popenCMD = shlex.split(samtoolsCMD)

    if is_tool("samtools"):
        # execute samtools view 2 times, first to check errors and second to save output. Cannot find better solution
        samtools_cmd_t = subprocess.Popen(samtools_popenCMD,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        stdout_s,stderr_s = samtools_cmd_t.communicate()
        stderr_samtool = list()
        for k in stderr_s:
            stderr_samtool.append(k)
        if len(stderr_samtool) != 0:
            sys.stderr.write("[E::calc_mgc] Error from samtools view for file "+strSAMfilename+":\n")
            subprocess.call(samtools_popenCMD,stdout=DEVNULL)
            sys.exit(1)
        samtools_cmd = subprocess.Popen(samtools_popenCMD,stdout=subprocess.PIPE,)
    else:
        sys.stderr.write("[E::calc_mgc] Error: samtools is not in the path. Cannot load the files.\n")
        sys.exit(1)

    output = samtools_cmd.stdout
    return(output)


# ------------------------------------------------------------------------------
# print error in the fasta file
# ------------------------------------------------------------------------------
def print_error_fasta(message,cont,fastq_file):
    sys.stderr.write("[E::main] Error: file "+fastq_file+" is not a fastq file\n")
    sys.stderr.write("[E::main] Line "+str(cont)+": "+message+"\n")
    sys.exit(1)

# ------------------------------------------------------------------------------
# function to check if the input is a fastq file and check the average length of
# the first 10,000 reads
# ------------------------------------------------------------------------------
def is_fastq(fastq_file,verbose):
    n_lines_head = 10000

    if not os.path.isfile(fastq_file):
        sys.stderr.write("[E::main] Error: "+fastq_file+': No such file.\n')
        sys.exit(1)

    try:
        zippedInput = False
        if (fastq_file.endswith(".gz")):
            zippedInput = True
            with GzipFile(fastq_file) as filef:
                head = list(islice(filef, n_lines_head))
        elif (fastq_file.endswith(".bz2")):
            zippedInput = True
            with BZ2File(fastq_file) as filef:
                head = list(islice(filef, n_lines_head))
        else:
            try:
                with open(fastq_file) as filef:
                    head = list(islice(filef, n_lines_head))
            except:
                sys.stderr.write("[E::main] Error loading file: "+fastq_file+"\n")
                sys.exit(1)
    except:
        sys.stderr.write("[E::main] Error. Cannot load the file "+fastq_file+"\n")
        sys.stderr.write("[E::main] Supported file are zipped .gz or .bz2, or plain text\n")
        sys.exit(1)


    # check min length of fastq file:
    if len(head) < 4:
        sys.stderr.write("[E::main] Error: file "+fastq_file+" is not a fastq file\n")
        sys.exit(1)

    # go through the lines in head
    all_lengths = list()
    cont_line = 1
    cont = 0
    for line in head:
        cont = cont + 1
        if zippedInput:
            line = line.decode('ascii')

        if cont_line == 1:
            if line[0] != "@": print_error_fasta("Line does not start with @",cont,fastq_file)

        if cont_line == 2:
            l1 = len(line)

        if cont_line == 3:
            if line[0] != "+": print_error_fasta("Line does not start with +",cont,fastq_file)

        if cont_line == 4:
            l2 = len(line)
            if (l1 != l2):
                print_error_fasta("Length of quality line is different from length of the read",cont,fastq_file)
            else:
                all_lengths.append(l1)

        if cont_line != 4:
            cont_line = cont_line + 1
        else:
            cont_line = 1

    if len(all_lengths) < 1:
        sys.stderr.write("[E::main] Error: file "+fastq_file+" is not a fastq file\n")
        sys.exit(1)

    if (all(x == all_lengths[0] for x in all_lengths)) and len(all_lengths) > 1:
        if verbose>=2: sys.stderr.write("[W::main] Warning: file "+fastq_file+". The length of the first "+str(len(all_lengths))+" reads is "+str(all_lengths[0])+". It is suggested to quality control the reads before profiling\n")


    avg_length = sum(all_lengths) / float(len(all_lengths))
    return int(avg_length)

# ------------------------------------------------------------------------------
# read the length from a bam/sam file
# ------------------------------------------------------------------------------
def read_length_from_bam_file(SAM_BAM_file):
    samLines = readSAMfile(SAM_BAM_file)
    try:
        for strSamline in samLines:
            strSamline = strSamline.decode('ascii')
            if (strSamline.startswith('@PG')):
                version_information_map_read = strSamline[14:].rstrip()
                all_info = version_information_map_read.split(" | ")
                if len(all_info) != 3: return None
                avg_length = float(all_info[2])
                return avg_length
            else:
                continue
    except:
        return None
    return None



# ------------------------------------------------------------------------------
# find the -l filter that was used to filter in map_tax (default is 45) from a bam/sam file
# ------------------------------------------------------------------------------
def read_filter_len_from_bam_file (SAM_BAM_file):
    samLines = readSAMfile(SAM_BAM_file)
    try:
        for strSamline in samLines:
            strSamline = strSamline.decode('ascii')
            if (strSamline.startswith('@CO min_len_alignment')):
                all_info = strSamline.rstrip().split(" ")
                length_filter = float(all_info[2])
                return length_filter
            else:
                continue
    except:
        return None
    return None
