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
import shutil
import subprocess
import hashlib

try:
	import requests
except:
	sys.stderr.write("Error: request library is not installed. Run:\npipenv install requests\n(check http://docs.python-requests.org/en/master/user/install/)")
	sys.exit(1)

# position of the script -------------------------------------------------------
path_mOTUs = os.path.realpath(__file__)
path_array = path_mOTUs.split("/")
relative_path = "/".join(path_array[0:-1])
relative_path = relative_path + "/"

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
# MAIN
# ------------------------------------------------------------------------------
def main(argv=None):
	sys.stderr.write(" ------------------------------------------------------------------------------\n")
	sys.stderr.write("|                               TEST MOTUS TOOL                                |\n")
	sys.stderr.write(" ------------------------------------------------------------------------------\n")

	# check if setup.py has been ran already -----------------------------------
	sys.stderr.write(" ran setup.py: ")
	if os.path.isdir(relative_path+'db_mOTU'):
	    sys.stderr.write("correct\n\n")
	else:
		sys.stderr.write("ERROR. Run setup.py\n\n")


	sys.stderr.write("-- Tools and versions\n")
	# check python version -----------------------------------------------------
	sys.stderr.write(" python:   ")
	python_version = sys.version_info
	if(python_version >= (2,7,0)):
		sys.stderr.write(" correct\n")
	else:
		sys.stderr.write(" ERROR: found v "+str(python_version[0])+"."+str(python_version[1])+"."+str(python_version[2])+". Required version 2.7 or 3.0 (or higher)\n")


	# check bwa ----------------------------------------------------------------
	sys.stderr.write(" bwa:      ")
	if is_tool("bwa"):
		sys.stderr.write(" correct\n")
		#TODO: maybe check version? at least 0.7.15-r1140
	else:
		sys.stderr.write(" WARNING\n         bwa is not in the path\n\n")

	# check samtools -----------------------------------------------------------
	sys.stderr.write(" samtools: ")
	if is_tool("samtools"):
		sys.stderr.write(" correct\n")
		#TODO: maybe check version? at least 1.5
	else:
		sys.stderr.write(" WARNING\n         samtools is not in the path\n\n")

	# check metaSNV ------------------------------------------------------------
	sys.stderr.write(" metaSNV:  ")
	if is_tool("metaSNV.py"):
		sys.stderr.write(" correct\n")
		#TODO: maybe check version? at least metaSNV v1.0.3
	else:
		sys.stderr.write(" WARNING\n         metaSNV is not in the path\n\n")








	return 0		# success







#-------------------------------- run main -------------------------------------
if __name__ == '__main__':
	status = main()
	sys.exit(status)
