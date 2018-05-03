#!/usr/bin/env python

# ============================================================================ #
# setup.py: prepare the motus tool after cloning from github
#
# Author: Alessio Milanese (milanese@embl.de)
#
# Main steps:
#    * Download the mOTUs database
#    * Create a file with the version information
#
# ============================================================================ #

motus_version = "0.6"
link_db = "https://oc.embl.de/index.php/s/znSp2EH6ca407hk/download"
md5_db = "a031a0e3b9402db1653bc7a447470e6c"

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



def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()



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
	sys.stderr.write(" ------------------------------------------------------------------------------\n")
	sys.stderr.write("|                              SETUP MOTUS TOOL                                |\n")
	sys.stderr.write(" ------------------------------------------------------------------------------\n")
	# download the files -------------------------------------------------------
	sys.stdout.write("Download the compressed motus database\n")

	db_name = relative_path+"db_mOTU.tar.gz"
	with open(db_name, "wb") as f:
		response = requests.get(link_db, stream=True)
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

	# check md5 ----------------------------------------------------------------
	sys.stdout.write("\nCheck md5: ")
	current_md5 = md5(db_name)

	if current_md5 == md5_db:
		sys.stdout.write("MD5 verified\n")
	else:
		sys.stdout.write("MD5 verification failed!\n")
		os.remove(db_name)
		sys.exit(1)


	# extract files ------------------------------------------------------------
	sys.stdout.write("Extract files from the archive...")
	extract_cmd = "tar -zxvf "+db_name+" -C "+relative_path
	try:
		FNULL = open(os.devnull, 'w')
		process = subprocess.Popen(extract_cmd.split(),stderr=FNULL)
		output, error = process.communicate()
	except:
		sys.stderr.write("Error: failed to extract files\n")
		sys.exit(1)
	if process.returncode:
		sys.stderr.write("Error: failed to extract files\n")
		sys.exit(1)
	else:
		sys.stderr.write("done\n")

	# --- remove db file
	sys.stdout.write("Remove zipped file...")
	os.remove(db_name)
	sys.stdout.write("done\n")



	# --------------- add file with version informations -----------------------
	sys.stdout.write("Add version file...")
	path_versions = relative_path + "db_mOTU/versions"
	try:
		outfile = tempfile.NamedTemporaryFile(delete=False, mode = "w")
		outfile.write("motus\t"+motus_version+"\n")
		outfile.write("#\tdatabase\n")
		outfile.write("nr\t"+motus_version+"\n")
		outfile.write("cen\t"+motus_version+"\n")
		outfile.write("#\tscripts\n")
		outfile.write("append\t"+motus_version+"\n")
		outfile.write("map_genes_to_mOTUs\t"+motus_version+"\n")
		outfile.write("map_mOTUs_to_LGs\t"+motus_version+"\n")
		outfile.write("runBWA\t"+motus_version+"\n")
		outfile.write("#\ttaxonomy\n")
		outfile.write("specI_tax\t"+motus_version+"\n")
		outfile.write("mOTULG_tax\t"+motus_version+"\n")

		# close file
		outfile.flush()
		os.fsync(outfile.fileno())
		outfile.close()
	except:
		sys.stderr.write("\nError while saving the file\n")
		sys.exit(1)

	try:
		shutil.move(outfile.name,path_versions)
	except:
		sys.stderr.write("\nError while saving the file\n")
		sys.exit(1)

	sys.stdout.write("done\n\n")




	return 0		# success

#-------------------------------- run main -------------------------------------
if __name__ == '__main__':
	status = main()
	sys.exit(status)
