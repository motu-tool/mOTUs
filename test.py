#!/usr/bin/env python

# ============================================================================ #
# test.py: test the motus profiler
#
# Author: Alessio Milanese (milanese@embl.de)
#
# Main steps:
#    * Test the tools and their versions
#    * Check the taxonomy profiling
#
# ============================================================================ #

import os
import sys
import tempfile
import subprocess
import shlex
import errno

# position of the script -------------------------------------------------------
path_mOTUs = os.path.realpath(__file__)
path_array = path_mOTUs.split("/")
relative_path = "/".join(path_array[0:-1])
relative_path = relative_path + "/"

# ------------------------------------------------------------------------------
# function to check if a specific tool exists
# ------------------------------------------------------------------------------
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

def is_tool(name):
    try:
        devnull = open(os.devnull)
        subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    return True


# ------------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------------
def main(argv=None):
    sys.stderr.write(" ------------------------------------------------------------------------------\n")
    sys.stderr.write("|                               TEST MOTUS TOOL                                |\n")
    sys.stderr.write(" ------------------------------------------------------------------------------\n")

    error_found = False

    # check if setup.py has been ran already -----------------------------------
    sys.stderr.write("\n1-- ran setup.py: ")
    if "--skip-setup-check" in sys.argv:
        sys.stderr.write("SKIPPED. As requested due to --skip-setup-check\n\n")
    else:
        # check that db_mOTU exists
        if not(os.path.isdir(relative_path+'db_mOTU')):
            sys.stderr.write("ERROR. Run setup.py\n\n")
            sys.exit(1)
        # check that temp dir exists
        if not(os.path.isdir(relative_path+'db_mOTU/db_mOTU_test')):
            sys.stderr.write("ERROR. db_mOTU/test directory is missing. Run setup.py\n\n")
            sys.exit(1)
        # check that files inside temp exists
        if not(os.path.isfile(relative_path+'db_mOTU/db_mOTU_test/test1_single.fastq') and os.path.isfile(relative_path+'db_mOTU/db_mOTU_test/test1.motus')):
            sys.stderr.write("ERROR. test files are missing. Run setup.py\n\n")
            sys.exit(1)
        sys.stderr.write("done\n\n")


    sys.stderr.write("2-- Tools and versions:\n")
    # check python version -----------------------------------------------------
    sys.stderr.write("- python:   ")
    python_version = sys.version_info
    if(python_version >= (2,7,0)):
        if(python_version >= (3,0,0)):
            sys.stderr.write(" correct\n")
        else:
            sys.stderr.write(" WARNING: python2 is not supported\n")
    else:
        sys.stderr.write(" ERROR: found v "+str(python_version[0])+"."+str(python_version[1])+"."+str(python_version[2])+". Required version 2.7 or 3.0 (or higher)\n")
        error_found = True


    # check bwa ----------------------------------------------------------------
    sys.stderr.write("- bwa:      ") # with bwa we don't check the return code, because the normal bwa return 1 if you dont add anything
    if is_tool("bwa"):
        sys.stderr.write(" correct\n")
        #TODO: maybe check version? at least 0.7.15-r1140
    else:
        sys.stderr.write(" WARNING. BWA is not in the path\n\n")
        error_found = True

    # check samtools -----------------------------------------------------------
    sys.stderr.write("- samtools: ")
    if is_tool("samtools"):
        if is_tool_and_return0("samtools --version"):
            sys.stderr.write(" correct\n")
        else:
            sys.stderr.write(" WARNING. Samtools is not correctly installed\n")
            error_found = True
    else:
        sys.stderr.write(" WARNING. Samtools is not in the path\n\n")
        error_found = True

    # check metaSNV ------------------------------------------------------------
    sys.stderr.write("- metaSNV:  ")
    if is_tool("metaSNV.py"):
        sys.stderr.write(" correct\n")
        #TODO: maybe check version? at least metaSNV v1.0.3
    else:
        sys.stderr.write(" WARNING. metaSNV is not in the path\n\n")
        error_found = True


    #===========================================================================
    #==================== test metagenomic profiling ===========================
    # Run motus on a test file
    sys.stderr.write("\n3-- Taxonomy profiling test:")

    error_flag = False

    if "--skip-setup-check" in sys.argv:
        sys.stderr.write("SKIPPED. Due to --skip-setup-check\n\n")
    else:
        test_file = relative_path+'db_mOTU/db_mOTU_test/test1_single.fastq'
        ground_truth_file = relative_path+'db_mOTU/db_mOTU_test/test1.motus'
        temp_file_profile = tempfile.NamedTemporaryFile(delete=False, mode="w")

        sys.stderr.write("\n- Run motus (-v 1, only error messages):\n")

        motus_command = "python "+relative_path+"motus profile -s "+test_file+" -g 1 -o "+temp_file_profile.name+" -v 1 -c -y insert.raw_counts"
        process = subprocess.Popen(motus_command.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()

        sys.stderr.write("- end motus call\n")
        sys.stderr.write("\nCheck resulting profile: ")
        profiled = open(temp_file_profile.name,"r")
        ground_truth = open(ground_truth_file,"r")

        gt = list()
        pr = list()
        for i in profiled:
            pr.append(i)
        for i in ground_truth:
            gt.append(i)
        profiled.close()
        ground_truth.close()

        # check values
        if len(pr) != len(gt):
            error_flag = True
            sys.stderr.write("ERROR. profiled sample is not correct\n\n")
        else:
            for i in range(1,len(pr)):
                if i != 0 and i != 1 and i != 2: # do not check the first 3 lines
                    gt1 = gt[i].split("\t")[0]
                    gt2 = int(float(gt[i].rstrip().split("\t")[1]))
                    pr1 = pr[i].split("\t")[0]
                    pr2 = int(float(pr[i].rstrip().split("\t")[1]))
                    if (gt1 != pr1):
                        error_flag = True
                        sys.stderr.write("ERROR. Different values in line "+str(i+1)+":\n"+gt1+"\n"+pr1+" \n\n")
                    if (gt2 != pr2):
                        error_flag = True
                        sys.stderr.write("ERROR. Different values in line "+str(i+1)+":\n"+str(gt2)+"\n"+str(pr2)+" \n\n")

        if not(error_flag):
            sys.stderr.write("correct\n\n")

        # remove temp file
        os.remove(temp_file_profile.name)


    if (error_found or error_flag):
        return 1
    else:
        return 0        # success

#-------------------------------- run main -------------------------------------
if __name__ == '__main__':
    status = main()
    sys.exit(status)
