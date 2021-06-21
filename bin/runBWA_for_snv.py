#!/usr/bin/env python

# ============================================================================ #
# runBWA.py: first step of the mOTU-LGs profiling
#
# ============================================================================ #

from __future__ import division
import os
import sys
import argparse
import shlex
import time
import subprocess
import re
import errno

log = ""


# ------------------------------------------------------------------------------
# function to check if a specific tool exists
# ------------------------------------------------------------------------------
def is_tool(name):
    try:
        devnull = open(os.devnull)
        subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    return True

# ------------------------------------------------------------------------------
# run bwa on a file that contains reads that are single end
# ------------------------------------------------------------------------------
def runBWA_singleEnd(strFilteredReadFile, reference, msamPercID, msamminLength, threads, technology, msam_script, msamOverlap,verbose):
    if verbose >= 5: log.print_message("bwa: values. msamPercID: "+str(msamPercID)+" msamminLength: "+str(msamminLength)+" msamOverlap: "+str(msamOverlap))
    try:
        from subprocess import DEVNULL
    except ImportError:
        DEVNULL = open(os.devnull, 'wb')

    techFlag = ""
    if (technology and technology != "illumina" and (technology == "pacbio" or technology == "ont2d" or technology == "intractg" )):
        techFlag = " -x " + technology

    if (threads):
        threadsFlag = " -t " + str(threads)
    else:
        threadsFlag = " -t 1"

    zippedInput = False
    # bwa can handle .gz files
    #if (strFilteredReadFile.endswith(".gz")):
    #    unzipCMD = "gunzip -c " + strFilteredReadFile
    #    zippedInput = True
    #    if not(is_tool("gunzip")):
    #        sys.stderr.write("[E::map_db] Error: gunzip is not installed. Cannot unzip the files\n")
    #        sys.exit(1)
    if (strFilteredReadFile.endswith(".bz2")):
        unzipCMD = "bunzip2 -c " + strFilteredReadFile
        zippedInput = True
        if not(is_tool("bunzip2")):
            log.print_error("bunzip2 is not installed. Cannot unzip the files")

    # check that bwa is in the path
    if not(is_tool("bwa")):
        log.print_error("BWA is not in the path. Cannot map the reads")


    # run bwa
    try:
        if (zippedInput):
            bwaCMD = "bwa mem -Y -v 1 -a" + techFlag + threadsFlag + " " + reference + " -"
        else:
            bwaCMD = "bwa mem -Y -v 1 -a" + techFlag + threadsFlag + " " + reference + " " + strFilteredReadFile

        if verbose >= 5: log.print_message("bwa call: "+bwaCMD)

        if (zippedInput):
            unzip_popenCMD = shlex.split(unzipCMD)
            unzip_cmd = subprocess.Popen(unzip_popenCMD,stdout=subprocess.PIPE,)

            bwa_popenCMD = shlex.split(bwaCMD)
            bwa_cmd = subprocess.Popen(bwa_popenCMD,stdin=unzip_cmd.stdout,stdout=subprocess.PIPE,stderr=DEVNULL)
        else:
            bwa_popenCMD = shlex.split(bwaCMD)
            bwa_cmd = subprocess.Popen(bwa_popenCMD,stdout=subprocess.PIPE,stderr=DEVNULL)

        min_perc_id=msamPercID
        min_length_align=msamminLength
        min_perc_cover=msamOverlap

        if verbose >= 5: log.print_message("Filter in bwa: MIN_PERC_ID:"+str(min_perc_id)+" MIN_LENGTH_ALIGN: "+str(min_length_align)+" MIN_PERC_COVER: "+str(min_perc_cover))

        # for the metaSNV we need to have the fasta/q sequence
        # we can find it from the previous reads
        prev_read_id = ""
        prev_fasta_seq = "A"
        prev_fastaq_seq = "A"
        for line in bwa_cmd.stdout:
            line = line.decode('ascii')

            if line[0]!="@":
                # insert the fasta/q lines if necessary
                all_vals = line.split("\t")
                if all_vals[9] == "*": # it means that it is a suboptimal alignment for bwa
                    if all_vals[0] == prev_read_id:
                        all_vals[9] = prev_fasta_seq
                        all_vals[10] = prev_fastaq_seq
                        line = "\t".join(all_vals)
                prev_read_id = all_vals[0]
                prev_fasta_seq = all_vals[9]
                prev_fastaq_seq = all_vals[10]



            #filter lines
            if line[0]!="@": # header
                arr = line.split("\t")
                if not((arr[1] == '4') or (arr[2] == '*') or (arr[5] == '*')):
                    len_seq = 0

                    min_query_al = 0

                    tott = 0

                    flag2 = False
                    tot = 0
                    for n,i in (re.findall('(\d+)([IDMSHX=])', arr[5])):
                        tott = tott + int(n)
                        if i == 'M' or i == 'D' or i == 'X' or i == '=':
                            tot = tot + int(n) # tot is the number of nucleotide that align
                        if i != 'I' and i != 'S' and i!= 'H':
                            len_seq = len_seq + int(n)
                        if i != 'H' and i != 'S' and i!= 'D':
                            min_query_al = min_query_al + int(n)

                    if tot >= float(min_length_align):
                        flag2 = True

                    flag1 = False
                    if ((int(arr[11].split(":")[2])/float(len_seq)) < (1-float(min_perc_id))): # TODO: here we assume that the NM value is in position 11, it should be checked
                        flag1 = True

                    # min. percent of the query that must be aligned, between 0 and 100 (required)
                    flag3 = False
                    if (float(min_query_al)/float(tott)) >= (float(min_perc_cover)/100 ):
                        flag3 = True

                    if flag1 and flag2 and flag3:
                        yield line

        #check that bzip finished correctly
        if (zippedInput):
            unzip_cmd.stdout.close()
            return_code = unzip_cmd.wait()
            if return_code:
                log.print_error("bunzip2 failed")

        # chack that bwa finished correctly
        bwa_cmd.stdout.close()
        return_code = bwa_cmd.wait()
        if return_code:
            log.print_error("bwa failed")


    except:
        log.print_error("Cannot call bwa on the file "+strFilteredReadFile)


# ------------------------------------------------------------------------------
# run the bwa mapping considering paired end
# ------------------------------------------------------------------------------
def runBWAmapping(forwardReads, reverseReads, singleReads, reference, threads, output, bamOutput, msam_script, technology, verbose, profile_mode, lane_id,msamminLength_from_motus,log_):

    # set up log
    global log
    log = log_
    # ----------------------

    # parameters for msamtools are fixed
    msamPercID = 0.97
    msamminLength = msamminLength_from_motus
    msamOverlap = 45

    ## check that the files exists
    if (forwardReads):
        if not os.path.isfile(forwardReads):
            log.print_error(forwardReads+': No such file')
        if not os.path.isfile(reverseReads):
            log.print_error(reverseReads+': No such file')
    if (singleReads):
        if not os.path.isfile(singleReads):
            log.print_error(singleReads+': No such file')

    # files
    sam_lines = list()
    mapped_sam_lines = list()
    sam_header = list()

    # computation for and rev --------------------------------------------------
    if (forwardReads):
        if verbose>2: start_time = time.time()

        #forward -----
        output_fwd = runBWA_singleEnd(forwardReads, reference, msamPercID, msamminLength, threads, technology, msam_script, msamOverlap,verbose)
        for line_b in output_fwd:
            line = line_b#.decode('ascii') # convert from binary to str
            if (not line.startswith("@")):
                orientation = "."+lane_id+'/1'
                line = line.replace("\t", orientation+"\t", 1)
                mapped_sam_lines.append(line)
            else:
                sam_header.append(line)

        if verbose>2: log.print_message("  map forward reads: " + str("{0:.2f}".format(time.time() - start_time))+" sec")

        # reverse -----
        if verbose>2: start_time = time.time()

        output_rev = runBWA_singleEnd(reverseReads, reference, msamPercID, msamminLength, threads, technology, msam_script, msamOverlap,verbose)

        for line_b in output_rev:
            line = line_b#.decode('ascii') # convert from binary to str
            if (not line.startswith("@")):
                orientation = "."+lane_id+'/2'
                line = line.replace("\t", orientation+"\t", 1)
                mapped_sam_lines.append(line)

        if verbose>2: log.print_message("  map reverse reads: " + str("{0:.2f}".format(time.time() - start_time))+" sec")

        # sort for and rev
        if verbose>2: start_time = time.time()
        mapped_sam_lines.sort()
        if verbose>2: log.print_message("  sort alignments: " + str("{0:.2f}".format(time.time() - start_time))+" sec")


    # computation single -------------------------------------------------------
    if (singleReads):
        if verbose>2: start_time = time.time()

        output_single = runBWA_singleEnd(singleReads, reference, msamPercID, msamminLength, threads, technology, msam_script, msamOverlap,verbose)
        for line_b in output_single:
            line = line_b#.decode('ascii') # convert from binary to str
            if (not line.startswith("@")):
                orientation = "."+lane_id+'.single'
                line = line.replace("\t", orientation+"\t", 1)
                mapped_sam_lines.append(line) #single ends are just appended (no need to sort)
            else:
                if (forwardReads==""): # if the header has not been printed already, then we print the header
                    sam_header.append(line)

        if verbose>2: log.print_message("  map single reads: " + str("{0:.2f}".format(time.time() - start_time))+" sec")

    # if we are running this as profile mode, then we return the list of the sam lines
    # without the header.
    if profile_mode:
        return mapped_sam_lines

    ## we always use profile mode, hence this part that follow of printing output is never used ##
    # print output -------------------------------------------------------------
    # concatenate header with sorted mapped reads
    sam_lines = sam_header + mapped_sam_lines

    if output != "":
        outfile = open(output, "w")
    else:
        outfile = sys.stdout

    if (bamOutput):
        #check that samtools is in the path
        if not(is_tool("samtools")):
            log.print_error("samtools is not in the path. Cannot save the file")

        convertCMD = "samtools view -b -S -"
        convert_popenCMD = shlex.split(convertCMD)
        convert_cmd = subprocess.Popen(convert_popenCMD,stdin=subprocess.PIPE,stdout=outfile,)
        for i in sam_lines:
            convert_cmd.stdin.write(i.encode())
    else:
        for i in sam_lines:
            outfile.write(i)
