#!/usr/bin/env python

# ============================================================================ #
# setup.py: prepare the mOTU profiler
#
# This will be better implemented as a package, for now it just creates a couple
# of files after cloning the directory
# ============================================================================ #

import os
import sys
import tempfile
import shutil
import subprocess



# position of the script -------------------------------------------------------
path_mOTUs = os.path.realpath(__file__)
path_array = path_mOTUs.split("/")
relative_path = "/".join(path_array[0:-1])
relative_path = relative_path + "/"

# ------------------------------------------------------------------------------
# function to check if a specific tool exists
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
	# path of the file that will contain the git commit id--------------------------
	sys.stdout.write("Creating file with the information of the git repository\n")
	path_version_git = relative_path + "db_mOTU/version_git_commit"

	if not(is_tool("git")):
		sys.stderr.write("Error: git is not in the path\n")
		sys.exit(1)

	gitCMD = "git --git-dir="+relative_path+".git log"
	try:
		git_cmd = subprocess.Popen(gitCMD,stdout=subprocess.PIPE,shell = True)
		stdout_s,stderr_s = git_cmd.communicate()
		if git_cmd.returncode:
			sys.stderr.write("Error while parsing git output\n")
			sys.exit(1)
		else:
			list_stdout = (stdout_s.decode('ascii')).split("\n")
		commit_id_dir = (list_stdout[0].split(" "))[1]
	except:
		sys.stderr.write("Error while parsing git output\n")
		sys.exit(1)


	# write file to tempfile
	try:
		outfile = tempfile.NamedTemporaryFile(delete=False, mode = "w")
		outfile.write("#\tgit commit id\n")
		outfile.write("git_commit_id\t"+commit_id_dir+"\n")
		outfile.write("#\ttag version\n")

		# ------------------------------------------------------tag version number ----------------
		outfile.write("git_tag_version\t0.5\n")

		# close file
		outfile.flush()
		os.fsync(outfile.fileno())
		outfile.close()
	except:
		sys.stderr.write("Error while saving the file\n")
		sys.exit(1)

	# move temp file to the final destination
	if os.path.isfile(path_version_git):
		sys.stderr.write("Warning. overwriting the old file\n")

	try:
		shutil.move(outfile.name,path_version_git)
	except:
		sys.stderr.write("Error while saving the file\n")
		sys.exit(1)


	return 0		# success

#-------------------------------- run main -------------------------------------
if __name__ == '__main__':
	status = main()
	sys.exit(status)
