#!/usr/bin/env python

# ============================================================================ #
# append.py: append profiles from different samples creating a table
#
# ============================================================================ #

from __future__ import division
import os
import sys
import argparse
import tempfile
import shutil
import datetime
import re

#function that detect the python version
def python_version():
    if(sys.version_info >= (3,0,0)):
        return(3)
    else:
        return(2)

def cArray(size1,size2):
    c = [[0. for i in range(size1)] for j in range(size2)]
    return c

# ------------------------------------------------------------------------------
def append_profilings(directory, list_files, output, verbose, BIOM_output,version_append,motu_call,version_tool):
    #--------------------------- save files ------------------------------------
    if directory is None and list_files is None:
        sys.stderr.write("[E::merge] Error: both -d and -i are empty")
        sys.exit(1)

    if list_files is not None:
        list_files = list_files.split(",")
        first_file = list_files[0]

    if directory is not None:
        try:
            list_files = os.listdir(directory)
        except:
            sys.stderr.write("[E::merge] Error: failed to open directory: "+directory+"\n")
            sys.exit(1)
        list_files = sorted(list_files)
        if verbose>2: sys.stderr.write(" [merge] Number of detected files: " +str(len(list_files))+"\n")
        first_file = directory+list_files[0]

    # first file to get informations
    if verbose>=5: sys.stderr.write("[merge] Opening file: "+first_file+"\n")
    try:
        location = open(first_file,'r')
    except:
        sys.stderr.write("[E::merge] Error: failed to read "+first_file+"\n")
        sys.exit(1)

    try:
        header_execution = location.readline()
        if header_execution[0:5] != "# git":
            if verbose > 5: sys.stderr.write("[E::merge] Error reading the first file - first line\n")
            sys.stderr.write("[E::merge] Error. truncated file: "+first_file+"\n")
            sys.exit(1)

        header_call = location.readline()
        if header_call[0:6] != "# call":
            sys.stderr.write("[E::merge] Error. truncated file: "+first_file+"\n")
            sys.exit(1)

        header_sample_name = location.readline().split('\t')
        len_info = len(header_sample_name) - 1
        header_ref = "\t".join(header_sample_name[0:len_info])

        taxa_id = list()
        for line in location:
            l = line.rstrip().split('\t')
            name = "\t".join(l[0:len_info])
            taxa_id.append(name)
        location.close()
    except:
        if verbose > 5: sys.stderr.write("[E::merge] Error reading the first file\n")
        sys.stderr.write("[E::merge] Error: failed to parse "+first_file+"\n")
        sys.exit(1)

    # create array
    array_c = cArray(len(list_files),len(taxa_id))

    #for all the files
    headers = ""
    cont_files = 0
    all_info_version = dict()
    for i in list_files:
        if directory is not None:
            file_name = directory+i
        else:
            file_name = i

        try:
            location = open(file_name,'r')
        except:
            sys.stderr.write("[E::merge] Error. failed to read "+file_name+"\n")
            sys.exit(1)

        # check the header/headers
        try:
            header_execution = location.readline().rstrip()
            if header_execution[0:5] != "# git":
                sys.stderr.write("[E::merge] Error. truncated file: "+file_name+"\n")
                sys.exit(1)
            else:
                all_info_version[cont_files] = header_execution

            header_call = location.readline()
            if header_call[0:6] != "# call":
                sys.stderr.write("[E::merge] Error. truncated file: "+file_name+"\n")
                sys.exit(1)

            # third header is the name
            header_3 = location.readline().rstrip().split('\t')
            if len(header_3) != len_info + 1:
                sys.stderr.write("[E::merge] Error: Different number of columns in file "+file_name+"\n")
                sys.exit(1)
            else:
                header_sample_name = header_3[len(header_3)-1]

        except:
            if verbose > 5: sys.stderr.write("[E::merge] Error parsing file\n")
            sys.stderr.write("[E::merge] Error: failed to parse "+file_name+"\n")
            sys.exit(1)

        headers = headers+header_sample_name+"\t"
        cont = 0
        for line in location:
            l = line.rstrip().split('\t')
            if len(l) != len_info + 1:
                sys.stderr.write("[E::merge] Error: Inconsistent number of columns in file "+file_name+"\n")
                sys.exit(1)

            # find the name
            name = "\t".join(l[0:len_info])

            if name == taxa_id[cont]:
                array_c[cont][cont_files] = l[len_info]
                cont += 1
            else:
                sys.stderr.write("[E::merge] Error: The taxa ids are different:\n")
                sys.stderr.write(i+": "+name+" -- "+taxa_id[cont]+"\n")
                sys.exit(1)
        cont_files += 1
        location.close()

    headers = headers[0:-1]

    # ---------------------- check the information of the version --------------
    if len(all_info_version) != 0:
        info_version_uniq = list(set(all_info_version.values()))
        if len(info_version_uniq) != 1:
            if verbose>=2:
                sys.stderr.write("[W::merge] Warning: The profiles that you are merging were analysed with different parameters:\n")
                for u in info_version_uniq:
                    sys.stderr.write(u+"\n")



    #----------------------------- append files --------------------------------
    if output != "":
        outfile = tempfile.NamedTemporaryFile(delete=False, mode="w")
        os.chmod(outfile.name, 0o644)
        #outfile = open(output, "w")
    else:
        outfile = sys.stdout

    if not BIOM_output: # header for normal output -----------------------------
        # first header
        if len(all_info_version) != 0:
            outfile.write(version_append+" | info merged profiles: ")
            for k in info_version_uniq:
                outfile.write(k)
                outfile.write(" ")
            outfile.write("\n")

        # second header
        outfile.write("# call: "+motu_call+"\n")

        # third header
        outfile.write(header_ref+"\t"+headers+"\n")
    else: # header for BIOM output ---------------------------------------------
        now = datetime.datetime.now()
        outfile.write("{\n    \"id\": \"None\",\n")
        outfile.write("    \"format\": \"Biological Observation Matrix 1.0.0\",\n")
        outfile.write("    \"format_url\": \"http://biom-format.org\",\n")
        outfile.write("    \"type\": \"OTU table\",\n")
        outfile.write("    \"generated_by\": \"motus v"+version_tool+"\",\n")
        outfile.write("    \"date\": \""+now.strftime("%Y-%m-%dT%H:%M:00")+"\",\n")
        outfile.write("    \"rows\":[\n")



    # VALUES - not BIOM format -------------------------------------------------
    if not BIOM_output:
        for i, taxa in enumerate(taxa_id):
            values = "\t".join(array_c[i])
            outfile.write(taxa + "\t" + values + "\n")

    # VALUES - BIOM format -----------------------------------------------------
    if BIOM_output:
        # rows informations
        n_col_header = len(taxa_id[1].split("\t"))
        for taxa in taxa_id[0:-1]:
            tt = taxa.split("\t")
            if n_col_header == 1:
                name = "null"
            if n_col_header == 2:
                name = "{\"NCBI_id\":\""+tt[1]+"\"}"
            if n_col_header == 3:
                name = "{\"name\":\""+tt[1]+"\",\n"
                name = name + "                                       \"NCBI_id\":\""+tt[2]+"\"}"

            outfile.write("            {\"id\":\""+tt[0]+"\", \"metadata\":"+name+"},\n")
        tt = taxa_id[-1].split("\t")
        if n_col_header == 1:
            name = "null"
        if n_col_header == 2:
            name = "{\"NCBI_id\":\""+tt[1]+"\"}"
        if n_col_header == 3:
            name = "{\"name\":\""+tt[1]+"\",\n"
            name = name + "                                       \"NCBI_id\":\""+tt[2]+"\"}"

        outfile.write("            {\"id\":\""+tt[0]+"\", \"metadata\":"+name+"}\n")

        # columns info --
        outfile.write("        ],\n")
        outfile.write("    \"columns\": [\n")
        headers_d = headers.split("\t")
        for k in headers_d[0:-1]:
            outfile.write('            {"id":"'+k+'", "metadata":null},\n')
        outfile.write('            {"id":"'+headers_d[-1]+'", "metadata":null}\n')
        outfile.write('        ],\n')


        outfile.write('    "matrix_type": "dense",\n')
        if re.search("\.", array_c[0][0]):
            outfile.write("    \"matrix_element_type\": \"float\",\n")
        else:
            outfile.write("    \"matrix_element_type\": \"int\",\n")
        outfile.write('    "shape": ['+str(len(taxa_id))+','+str(len(headers_d))+'],\n')

        # data ------------
        outfile.write('    "data":  [')
        values = ",".join(array_c[0])
        outfile.write('['+values+'],\n')
        for i in range(1,(len(taxa_id)-1)):
            values = ",".join(array_c[i])
            outfile.write('              ['+values+'],\n')

        values = ",".join(array_c[len(taxa_id)-1])
        outfile.write('              ['+values+']]\n')
        outfile.write('}')




    if output != "":
        if verbose>2: sys.stderr.write(" [merge] (Saving the the merged profiles)\n")
        try:
            outfile.flush()
            os.fsync(outfile.fileno())
            outfile.close()
        except:
            sys.stderr.write("[E::main] Error: failed to save the merged profiles\n")
            sys.exit(1)
        try:
            #os.rename(outfile.name,output) # atomic operation
            shutil.move(outfile.name,output) #It is not atomic if the files are on different filsystems.
        except:
            sys.stderr.write("[E::merge] Error: failed to save the merged profiles\n")
            sys.stderr.write("[E::merge] you can find the file here:\n"+outfile.name+"\n")
            sys.exit(1)




# ------------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------------
def main(argv=None):

    #----------------------------- input parameters ----------------------------
    parser = argparse.ArgumentParser(description='This program calculates mOTU abundances for one sample', add_help=True)
    parser.add_argument('--directory', '-d', action="store", default=None,dest='directory', help='append all the files in the directory')
    parser.add_argument('--output', '-o', action="store", dest='output', default="", help='name of the output file, if not specified is stdout')
    parser.add_argument('-v', action='store', type=int, default=3, dest='verbose', help='Verbose levels')
    parser.add_argument('-i', action="store", dest='listInputFiles', default="", help='name of input file(s); sam or bam formatted files. If it is a list: insert all files separated by a comma')
    parser.add_argument('-B', action='store_true', default=False, dest='BIOM_output', help='print output in BIOM format')
    args = parser.parse_args()

    #-------------------------------- check input ------------------------------
    # check that there is at least one file with reads
    if (args.directory==""):
        sys.stderr.write("[E::merge] Error: insert the directory value.\n")
        sys.exit(1)

    if args.directory is not None:
        if (args.directory[-1]!="/"):
            args.directory = args.directory + "/"

    # call the function --------------------------------------------------------
    append_profilings(args.directory, args.listInputFiles, args.output, args.verbose, args.BIOM_output)


    return 0        # success

#-------------------------------- run main -------------------------------------
if __name__ == '__main__':
    status = main()
    sys.exit(status)
