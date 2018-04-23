#!/usr/bin/env python

# ============================================================================ #
# setup.py: prepare the mOTU profiler
# ============================================================================ #

motus_version = "0.6"

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
	# download the files -------------------------------------------------------
	sys.stdout.write("Download the compressed motus database\n")

	link = "https://oc.embl.de/index.php/s/gqM8uOn3OZfRzB5/download"
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

	# extract files ------------------------------------------------------------
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

	# --- remove db file
	sys.stdout.write("\nRemove zipped file...")
	os.remove(db_name)
	sys.stdout.write("done\n")



	# --------------- add file with version informations -----------------------
	sys.stdout.write("\nAdd version file...")
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

	sys.stdout.write("done\n")




	return 0		# success

#-------------------------------- run main -------------------------------------
if __name__ == '__main__':
	status = main()
	sys.exit(status)
