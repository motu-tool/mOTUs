import os
import sys
import tempfile
import shutil
import subprocess
import hashlib
import time
import urllib.request
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

motus_version = "3.1.0"
link_db = "https://zenodo.org/record/7778108/files/db_mOTU_v3.1.0.tar.gz"
md5_db = "f841c36150025af837f7a9a358c9a3c3"
DOI_db = "10.5281/zenodo.7778108"



# function to print progress bar for python 3
def reporthook(count, block_size, total_size):
    global start_time
    if count == 0:
        start_time = time.time()
        return
    duration = time.time() - start_time
    progress_size = int(count * block_size)
    speed = int(progress_size / (1024 * duration))
    percent = int(count * block_size * 100 / total_size)
    sys.stdout.write("\r %d%%, %d MB, %d KB/s, %d seconds passed" %
                    (percent, progress_size / (1024 * 1024), speed, duration))
    sys.stdout.flush()

def save_f(url, filename):
    if "--no-download-progress" in sys.argv:
        urllib.request.urlretrieve(url, filename)
    else:
        urllib.request.urlretrieve(url, filename, reporthook)


# function to check md5
def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()



# position of the script -------------------------------------------------------
path_mOTUs = os.path.realpath(__file__) # RECODE
path_array = path_mOTUs.split("/")
relative_path = "/".join(path_array[0:-1])
relative_path = relative_path + "/"

# ------------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------------
def main(argv=None):
    sys.stderr.write(" ------------------------------------------------------------------------------\n")
    sys.stderr.write("|                              SETUP MOTUS TOOL                                |\n")
    sys.stderr.write(" ------------------------------------------------------------------------------\n")
    # download the files -------------------------------------------------------
    path_versions = relative_path + "db_mOTU/db_mOTU_versions"
    if os.path.isfile(path_versions) and "--force-redownload" not in sys.argv:
        sys.stderr.write("Database already downloaded. Not doing anything.\n"
                         "Use --force-redownload to download again.\n")
        sys.exit(0)

    sys.stderr.write("Download the compressed motus database (~3.5Gb)\n")
    db_name = relative_path+"db_mOTU.tar.gz"



    save_f(link_db, db_name)
    sys.stderr.write("\n")

    # check md5 ----------------------------------------------------------------
    sys.stderr.write("\nCheck md5: ")
    current_md5 = md5(db_name)
    sys.stderr.flush()

    if current_md5 == md5_db:
        sys.stderr.write("MD5 verified\n")
    else:
        sys.stderr.write("MD5 verification failed!\n")
        os.remove(db_name)
        sys.exit(1)
    sys.stderr.flush()


    # extract files ------------------------------------------------------------
    sys.stderr.write("Extract files from the archive...")
    sys.stderr.flush()
    extract_cmd = "tar -zxvf "+db_name+" -C "+relative_path
    try:
        FNULL = open(os.devnull, 'w')
        process = subprocess.Popen(extract_cmd.split(),stderr=FNULL,stdout=FNULL)
        output, error = process.communicate()
    except:
        sys.stderr.write("Error: failed to extract files\n")
        sys.exit(1)
    if process.returncode:
        sys.stderr.write("Error: failed to extract files\n")
        sys.exit(1)
    else:
        sys.stderr.write("done\n")
    sys.stderr.flush()

    # --- remove db file
    sys.stderr.write("Remove zipped file...")
    sys.stderr.flush()

    os.remove(db_name)
    sys.stderr.write("done\n")
    sys.stderr.flush()



    # --------------- add file with version informations -----------------------
    sys.stderr.write("Add version file...")
    sys.stderr.flush()
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

    sys.stderr.write("done\n\n")
    sys.stderr.flush()




    return 0        # success

#-------------------------------- run main -------------------------------------
if __name__ == '__main__':
    status = main()
    sys.exit(status)
