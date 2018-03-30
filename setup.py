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
	sys.stdout.write("\n--- INSTALL MOTUS v2 ---\n")
	# path of the file that will contain the git commit id--------------------------
	sys.stdout.write("Download the compressed motus database\n")

	link = "https://oc.embl.de/index.php/s/wrz9YfKfrNyYCaY/download"
	db_name = relative_path+"db_mOTU.tar.gz"
	with open(db_name, "wb") as f:
		response = requests.get(link, stream=True)
		total_length = response.headers.get('content-length')

		if total_length is None: # no content length header
			f.write(response.content)
		else:
			dl = 0
			total_length = int(total_length)
			for data in response.iter_content(chunk_size=4096):
				dl += len(data)
				f.write(data)
				done = int(50 * dl / total_length)
				sys.stdout.write("\r[%s%s]" % ('=' * done, ' ' * (50-done)) )
				sys.stdout.flush()
		sys.stdout.write("\n")

	sys.stdout.write("\nExtract files from the archive:\n")
	extract_cmd = "tar -zxvf "+db_name+" -C "+relative_path
	try:
		process = subprocess.Popen(extract_cmd.split(), stdout=subprocess.PIPE)
		output, error = process.communicate()
	except:
		sys.stderr.write("Error: failed to extract files\n")
		sys.exit(1)
	if process.returncode:
		sys.stderr.write("Error: failed to extract files\n")
		sys.exit(1)




	return 0		# success

#-------------------------------- run main -------------------------------------
if __name__ == '__main__':
	status = main()
	sys.exit(status)
