#!/usr/bin/env python

from __future__ import division
import os
import sys
import argparse
import re
import tempfile
import shutil
import datetime

def save_file_to_dict(file_r,col_key,col_value,header,divide_with_tab,skip_not_profilable):
    # skip_not_profilable, if it's True, we dont save in the dictionary if col_key = "not_profilable"
    try:
        location = open(file_r,'r')
        if header:
            header_res = location.readline().rstrip()

        res_dict = dict()
        for line in location:
            l = line.rstrip().split('\t')
            if not(l[col_key] == "not_profilable" and skip_not_profilable is True):
                if divide_with_tab:
                    # divide first word from the rest with a tab
                    ncbi_id = l[col_value].split(" ")[0]
                    value_correct = l[col_value].split(" ")[1:]
                    res_dict[l[col_key]] = ncbi_id+"\t"+(" ".join(value_correct))
                else:
                    res_dict[l[col_key]] = l[col_value]
        location.close()
    except:
        sys.stderr.write("[E::calc_motu] Error loading file: "+file_r+"\n[E::calc_motu] Try to download again the motus profiler\n\n")
        sys.exit(1)

    if header:
        return header_res,res_dict
    else:
        return res_dict



def save_file_to_dict_full_rank(file_r,col_key,col_value,header,skip_not_profilable):
    # skip_not_profilable, if it's True, we dont save in the dictionary if col_key = "not_profilable"
    initial_val = 1
    if skip_not_profilable: initial_val = 2
    beginning_tax = ["k__","p__","c__","o__","f__","g__","s__"]
    # col_key should be zero
    try:
        location = open(file_r,'r')
        if header:
            header_res = location.readline().rstrip()

        res_dict = dict()
        for line in location:
            # for every line we have to take the full rank
            l = line.rstrip().split('\t')
            ncbi_id_temp = list()
            consensus_name_temp = list()
            for k in range(initial_val,col_value+1):
                # divide string in the spaces
                temp_name = l[k].split(' ')
                # ncbi id is the first
                ncbi_id_temp.append(temp_name[0])
                # consensus name is the rest
                tempor_name = " ".join(temp_name[1:])
                tempor_name = beginning_tax[k-initial_val] + tempor_name
                consensus_name_temp.append(tempor_name)

            ncbi_id_join = "|".join(ncbi_id_temp)
            cons_name_join = "|".join(consensus_name_temp)

            if not(l[col_key] == "not_profilable" and skip_not_profilable is True):
                res_dict[l[col_key]] = ncbi_id_join + "\t" + cons_name_join

        location.close()
    except:
        sys.stderr.write("[E::calc_motu] Error loading file: "+file_r+"\n[E::calc_motu] Try to download again the motus profiler\n\n")
        sys.exit(1)

    if header:
        return header_res,res_dict
    else:
        return res_dict



def save_file_to_dict_two_headers(file_r,col_key,col_value,header,remove_first_value):
    try:
        location = open(file_r,'r')
        if header:
            header_res1 = location.readline().rstrip()
            header_res2 = location.readline().rstrip()
    except:
        sys.stderr.write("[E::calc_motu] Error loading file: "+file_r+"\n")
        sys.exit(1)

    res_dict = dict()
    for line in location:
        try:
            l = line.rstrip().split('\t')
            if remove_first_value:
                # remove the first value l[col_value]
                value_correct = l[col_value].split(" ")[1:]
                res_dict[l[col_key]] = " ".join(value_correct)
            else:
                res_dict[l[col_key]] = l[col_value]
        except:
            sys.stderr.write("[E::calc_motu] Error with input file: "+file_r+"\n")
            sys.stderr.write("[E::calc_motu] truncated file\n")
            sys.exit(1)

    location.close()
    if header:
        return header_res1,header_res2,res_dict
    else:
        return res_dict

def calculate_abundance(infile, LGs_map, LGs_map_l, specI_taxonomy, mOTULG_taxonomy, output, cutoff, onlySpecI, sampleName, taxonomic_level, BIOM_output, profile_mode,input_for_profile, print_NCBI_id, print_rel_ab,mgc_table_header,version_map_lgs,motu_version_tool,verbose,motu_call,git_commit_id,print_full_rank,print_full_name,short_names_file,version_tool):

    if print_full_rank:
        print_full_name = True
        # if we print all the rank of the taxonomy, then we print also the full
        # name (and not the short names)
    if BIOM_output:
        print_full_name = True

    # load data ----------------------------------------------------------------
    if print_full_rank:
        # load the taxonomy for the specI - first always map at the species level
        taxonomy_header_s, taxonomy_s = save_file_to_dict_full_rank(specI_taxonomy,1,8,True,True)
        # load the taxonomy for the mOTU_LGs
        taxonomy_header_m, taxonomy_m = save_file_to_dict_full_rank(mOTULG_taxonomy,0,7,True,False)
    else:
        # load the taxonomy for the specI - first always map at the species level
        taxonomy_header_s, taxonomy_s = save_file_to_dict(specI_taxonomy,1,8,True,True,True)
        # load the taxonomy for the mOTU_LGs
        taxonomy_header_m, taxonomy_m = save_file_to_dict(mOTULG_taxonomy,0,7,True,True,False)

    # laod short names
    shortNames_header, shortNames = save_file_to_dict(short_names_file,0,1,True,True,True)

    # load the mOTU read counts (output from map_genes_to_mOTUs.py)
    if profile_mode:
        mOTUs_ab = input_for_profile
        sample_id_header = sampleName
    else:
        info_computation_so_far,sample_id_header, mOTUs_ab = save_file_to_dict_two_headers(infile,0,1,True,False)
        if sampleName != "":
            sample_id_header = sampleName
        mgc_table_header = info_computation_so_far

    # check that the header from the mgc_table is correct
    if len(mgc_table_header.split(" | ")) != 4:
        if verbose>1: sys.stderr.write(" [W::calc_motu] Warning: couldn't find any header in the mgc table\n")
        mgc_table_header = "# map_tax unknown | gene database: unknown | 100 | calc_mgc unknown"


    ## map_lgs header --
    map_lgs_header = mgc_table_header+" | " + version_map_lgs
    map_lgs_header = map_lgs_header[1:]
    map_lgs_header = "# motus version "+str(motu_version_tool)+" |"+map_lgs_header

    # open the map from mOTU to LGs
    mOTUs_LGs = save_file_to_dict(LGs_map,0,1,False,False,False)

    # open the map from mOTU to LGs that are in a line
    try:
        location = open(LGs_map_l,'r')
        mOTUs_LGs_l = dict()

        list_LGs = list() # we want to preserve the order. Hence, we create the list here

        for line in location:
            l = line.rstrip().split('\t')
            mOTUs_LGs_l[l[0]] = l[1]
            list_LGs.append(l[0])
        location.close()
    except:
        sys.stderr.write("[E::calc_motu] Error loading file: "+LGs_map_l+"\n[E::calc_motu] Try to download again the motus profiler\n\n")
        sys.exit(1)

    # check that the mgc_table is correct --------------------------------------
    error_flag_mgc_table = False
    all_wrong = True
    for k in mOTUs_ab:
        if k in mOTUs_LGs:
            all_wrong = False
        else:
            error_flag_mgc_table = True


    if len(mOTUs_ab) == 0:
        all_wrong = False
        if verbose>1: sys.stderr.write(" [W::calc_motu] Warning: The mgc table is empty\n")


    if all_wrong:
        if profile_mode:
            sys.stderr.write("[E::calc_motu] Error: the mgc table does not contain information of the mgc\n")
        else:
            sys.stderr.write("[E::calc_motu] Error in file "+infile+":\n[E::calc_motu] the mgc table does not contain information of the mgc\n")
        sys.exit(1)

    if (not all_wrong) and (error_flag_mgc_table):
        for k in mOTUs_ab:
            if not (k in mOTUs_LGs):
                sys.stderr.write(" [W::calc_motu] Warning: \'"+k+"\' not a mgc. Ignore the line\n")



    # create some useful variables ---------------------------------------------
    list_mOTUs = list(set(mOTUs_LGs.keys()))
    counts_mOTUs = dict()
    rel_ab_LGs = dict()

    #calculate -----------------------------------------------------------------
    for i in list_mOTUs:
        counts_mOTUs[i] = 0
    for i in mOTUs_ab:
        counts_mOTUs[i] = mOTUs_ab[i]

    for j in list_LGs: # for every LG
        genes_l = mOTUs_LGs_l[j] # we find the mOTUs that compose the LG
        genes_list = genes_l.split(';')
        counts_mOTUs_j = [counts_mOTUs[x] for x in genes_list] # vector that represents the read counts of the mOTUs of the LGs
        counts_mOTUs_j = [float(numeric_string) for numeric_string in counts_mOTUs_j] # transform from string to float

        list_diff_zero = list()
        cog_type = list()
        for i in range(len(counts_mOTUs_j)):
            if counts_mOTUs_j[i]>0:
                list_diff_zero.append(counts_mOTUs_j[i]) # find the one that are different from zero
                cog_type.append(genes_list[i].split('.')[0]) # and save the COG type

        rel_ab_LGs[j] = 0
        if j != '-1':
            if len(list_diff_zero) >= cutoff:
                if len(list_diff_zero) == 1:
                    rel_ab_LGs[j] = float(sum(list_diff_zero))
                else:
                    list_diff_zero.sort()
                    if (len(list_diff_zero) % 2) == 1:
                        #take median
                        pos_median = int(len(list_diff_zero)/2 - 0.5)
                        rel_ab_LGs[j] = float(list_diff_zero[pos_median])
                    else:
                        #take average of two medians
                        pos_median1 = int(len(list_diff_zero)/2)
                        pos_median2 = int(len(list_diff_zero)/2 - 1)
                        rel_ab_LGs[j] = float(list_diff_zero[pos_median1] + list_diff_zero[pos_median2]) / 2
        else: # what to do with -1: take the average of the different COGs values
            all_cog_type = list(set(cog_type))
            count_cogs_m_1 = dict()
            # set to zero to start
            for ll in all_cog_type:
                count_cogs_m_1[ll] = 0
            # add values to the dicrionary
            for ll in range(len(cog_type)):
                count_cogs_m_1[cog_type[ll]] = count_cogs_m_1[cog_type[ll]] + list_diff_zero[ll]
            # calculate average for every COG
            mean_vals_cogs = list()
            for ll in all_cog_type:
                if count_cogs_m_1[ll] != 0:
                    mean_vals_cogs.append(count_cogs_m_1[ll])

            if sum(mean_vals_cogs) != 0:
                mean_vals_cogs.sort()
                if (len(mean_vals_cogs) % 2) == 1:
                    pos_median = int(len(mean_vals_cogs)/2 - 0.5)
                else:
                    pos_median = int(len(mean_vals_cogs)/2 - 1)
                rel_ab_LGs[j] = mean_vals_cogs[pos_median]
            else:
                rel_ab_LGs[j] = 0


    # divide by sum
    rel_ab_LGs_rel = dict()
    s = sum(rel_ab_LGs.values())
    if not print_rel_ab:
        if s != 0:
            for j in list_LGs:
                rel_ab_LGs_rel[j] = float(rel_ab_LGs[j])/s
        else:
            for j in list_LGs:
                rel_ab_LGs_rel[j] = 0
            if verbose>1: sys.stderr.write(" [W::calc_motu] Warning: The relative abundance is 0 for all the mOTUs\n")
    else: # if we dont print the rel. ab.
        base_coverage_flag = False
        if re.search("base.coverage", mgc_table_header):
            base_coverage_flag = True

        for j in list_LGs:
            if base_coverage_flag:
                rel_ab_LGs_rel[j] = float(rel_ab_LGs[j])
            else:# if we are using insert_* then we round the counts
                rel_ab_LGs_rel[j] = round(rel_ab_LGs[j])

    # keep only specI
    if onlySpecI:
        rel_ab_LGs_rel_temp = dict(rel_ab_LGs_rel)
        rel_ab_LGs_rel = dict()
        value_minus1 = 0
        for j in list_LGs:
            if j == '-1':
                value_minus1 = value_minus1 + rel_ab_LGs_rel_temp[j]
            else:
                type_c = j.split("_")[0]
                if (type_c == 'meta'):
                    value_minus1 = value_minus1 + rel_ab_LGs_rel_temp[j]
                else:
                    rel_ab_LGs_rel[j] = rel_ab_LGs_rel_temp[j]
        rel_ab_LGs_rel['-1'] = value_minus1

    # general print
    if output != "":
        #outfile = open(output, "w")
        outfile = tempfile.NamedTemporaryFile(delete=False, mode="w")
        os.chmod(outfile.name, 0o644)
    else:
        outfile = sys.stdout

    # remove the 4th value of the header because is the avg. read length
    list_parameters = map_lgs_header.split(' | ')
    del list_parameters[3]
    map_lgs_header = ' | '.join(list_parameters)

    if not BIOM_output:
        outfile.write(git_commit_id+" | "+map_lgs_header[1:]+"\n")
        outfile.write("# call: "+motu_call+"\n")

    # third header
    if taxonomic_level == "mOTU" and print_full_name:
        header_complete = taxonomic_level+"\tconsensus_taxonomy\t"
    else:
        header_complete = "consensus_taxonomy\t"
    if print_NCBI_id: header_complete = header_complete + "NCBI_tax_id\t"
    header_complete += sample_id_header+ "\n"

    if not BIOM_output:
        outfile.write("#"+header_complete)


    # print header for BIOM format (hardcoded) =================================
    if BIOM_output:
        now = datetime.datetime.now()
        outfile.write("{\n    \"id\": \""+sample_id_header+"\",\n")
        outfile.write("    \"format\": \"Biological Observation Matrix 1.0.0\",\n")
        outfile.write("    \"format_url\": \"http://biom-format.org\",\n")
        outfile.write("    \"type\": \"OTU table\",\n")
        outfile.write("    \"generated_by\": \"motus v"+version_tool+"\",\n")
        outfile.write("    \"date\": \""+now.strftime("%Y-%m-%dT%H:%M:00")+"\",\n")
        outfile.write("    \"rows\":[\n")
        list_rows = list()
        list_vals = list()
        first_val = rel_ab_LGs_rel["-1"]


    # print result FOR SPECIES LEVEL ===========================================
    # print full name------------------------
    if taxonomic_level == "mOTU" and print_full_name:
        # preapre data
        list_LGs_print = list(list_LGs)
        if onlySpecI:
            list_LGs_print = list()
            for j in list_LGs:
                if j != '-1':
                    type_c = j.split("_")[0]
                    if (type_c != 'meta'):
                        list_LGs_print.append(j)
            list_LGs_print.append("-1")

        for j in list_LGs_print:
            if ((j in taxonomy_s) or (j in taxonomy_m)):
                if j in taxonomy_s: all_val = taxonomy_s[j].split("\t")
                if j in taxonomy_m: all_val = taxonomy_m[j].split("\t")
                # line to print
                if not BIOM_output:
                    name = j+"\t" # mOTU id
                    name = name + all_val[1]# consensus_name
                    if print_NCBI_id: name = name + "\t"+all_val[0] # NCBI_tax_id
                    name = name + "\t" + str(rel_ab_LGs_rel[j]) +"\n" # value
                    outfile.write(name)
                else:
                    name = "{\"name\":\""+all_val[1]+"\",\n"
                    name = name + "                                       \"NCBI_id\":\""+all_val[0]+"\"}"
                    list_rows.append("            {\"id\":\""+j+"\", \"metadata\":"+name+"},\n")
                    list_vals.append("["+str(rel_ab_LGs_rel[j])+"],\n")
            elif j == "-1": # -1
                name = "-1\t-1"
                if print_NCBI_id: name = name + "\tNA"
                name = name + "\t" +str(rel_ab_LGs_rel[j])+"\n"
                if not BIOM_output:
                    outfile.write(name)
                else:
                    list_rows.append("            {\"id\":\"-1\", \"metadata\":{\"name\":\"unknown\",\n                                    \"NCBI_id\":\"NA\"}}\n")
                    list_vals.append("["+str(rel_ab_LGs_rel[j])+"]]\n")

            else: # if it not in anyone (it should not happen)
                if verbose>1: sys.stderr.write(" [W::calc_motu] Warning: find mOTU "+j+" that is not present in the taxonomy\n")

    # print short name - deafult ----------------
    if taxonomic_level == "mOTU" and (not (print_full_name)):
        # preapre data
        list_LGs_print = list(list_LGs)
        if onlySpecI:
            list_LGs_print = list()
            for j in list_LGs:
                if j != '-1':
                    type_c = j.split("_")[0]
                    if (type_c != 'meta'):
                        list_LGs_print.append(j)
            list_LGs_print.append("-1")

        for j in list_LGs_print:
            if ((j in taxonomy_s) or (j in taxonomy_m)):
                if j in taxonomy_s:
                    all_val = taxonomy_s[j].split("\t")
                    name = shortNames[j].split("\t")[1]
                if j in taxonomy_m:
                    all_val = taxonomy_m[j].split("\t")
                    name = all_val[1] # the short name for meta-mOTUs is the normal name
                # line to print
                name = name + " ["+j+"]"
                if print_NCBI_id: name = name + "\t"+all_val[0] # NCBI_tax_id
                name = name + "\t" + str(rel_ab_LGs_rel[j]) +"\n" # value
                outfile.write(name)
            elif j == "-1": # -1
                name = "-1"
                if print_NCBI_id: name = name + "\tNA"
                name = name + "\t" +str(rel_ab_LGs_rel[j])+"\n"
                outfile.write(name)
            else: # if it not in anyone (it should not happen)
                if verbose>1: sys.stderr.write(" [W::calc_motu] Warning: find mOTU "+j+" that is not present in the taxonomy\n")


    # print result FOR NOT SPECIES LEVEL =======================================
    if taxonomic_level != "mOTU":
        # choose the right taxonomy level
        taxonomic_levels = ["mOTULG_cluster","kingdom","phylum","class","order","family","genus","species"]
        pos = taxonomic_levels.index(taxonomic_level)

        # load the taxonomy for specific taxonomic level
        if print_full_rank:
            taxonomy_header_s, taxonomy_s_2 = save_file_to_dict_full_rank(specI_taxonomy,1,pos+1,True,True)
            taxonomy_header_m, taxonomy_m_2 = save_file_to_dict_full_rank(mOTULG_taxonomy,0,pos,True,False)
        else:
            taxonomy_header_s, taxonomy_s_2 = save_file_to_dict(specI_taxonomy,1,pos+1,True,True,True)
            taxonomy_header_m, taxonomy_m_2 = save_file_to_dict(mOTULG_taxonomy,0,pos,True,True,False)


        # create list of unique values ---> it is based on the string: NCBI_id + consensus_name
        if onlySpecI:
            list_taxon = list(set(taxonomy_s_2.values()))
        else:
            list_taxon = list(set(list(taxonomy_s_2.values())+list(taxonomy_m_2.values())))

        list_taxon.sort()
        list_taxon.append("-1")
        rel_abundance_taxon = dict()
        for i in list_taxon:
            rel_abundance_taxon[i] = 0
        for i in rel_ab_LGs_rel:
            if i != "-1":
                if i in taxonomy_s_2:
                    rel_abundance_taxon[taxonomy_s_2[i]] = rel_abundance_taxon[taxonomy_s_2[i]] + rel_ab_LGs_rel[i]
                if i in taxonomy_m_2:
                    rel_abundance_taxon[taxonomy_m_2[i]] = rel_abundance_taxon[taxonomy_m_2[i]] + rel_ab_LGs_rel[i]
            else:
                rel_abundance_taxon[i] = rel_ab_LGs_rel[i]

        # print
        for i in list_taxon:
            if i != "-1":
                all_val = i.split("\t")
                # prepare line to print
                name = all_val[1] # consensus_name
                if print_NCBI_id: name = name + "\t"+all_val[0] # NCBI_tax_id
                name = name + "\t" + str(rel_abundance_taxon[i]) +"\n" # value
                if not BIOM_output:
                    outfile.write(name)
                else:
                    name = "{\"NCBI_id\":\""+all_val[0]+"\"}"
                    list_rows.append("            {\"id\":\""+all_val[1]+"\", \"metadata\":"+name+"},\n")
                    list_vals.append("["+str(rel_abundance_taxon[i])+"],\n")

            else:
                if not BIOM_output:
                    if print_NCBI_id:
                        name = "-1\tNA\t"+str(rel_abundance_taxon[i])+"\n"
                    else:
                        name = "-1\t"+str(rel_abundance_taxon[i])+"\n"
                    outfile.write(name)
                else:
                    list_rows.append("            {\"id\":\"-1\", \"metadata\":{\"NCBI_id\":\"NA\"}}\n")
                    list_vals.append("["+str(rel_abundance_taxon[i])+"]]\n")

    ####### PRINT BIOM FORMAT ##################################################
    if BIOM_output:
        for i in list_rows:
            outfile.write(i)
        outfile.write("        ],\n")
        outfile.write("    \"columns\": [\n")
        outfile.write("            {\"id\":\""+sample_id_header+"\", \"metadata\":null}\n")
        outfile.write("        ],\n")
        outfile.write("    \"matrix_type\": \"dense\",\n")
        if (type(first_val) is int):
            outfile.write("    \"matrix_element_type\": \"int\",\n")
        else:
            outfile.write("    \"matrix_element_type\": \"float\",\n")
        outfile.write("    \"shape\": ["+str(len(list_rows))+",1],\n")
        outfile.write("    \"data\":  [")
        outfile.write(list_vals[0])
        for i in list_vals[1:]:
            outfile.write("              "+i)
        outfile.write("}")


    # move the temp file to the final destination ------------------------------
    if output != "":
        try:
            outfile.flush()
            os.fsync(outfile.fileno())
            outfile.close()
        except:
            sys.stderr.write("[E::main] Error: failed to save the profile\n")
            sys.exit(1)
        try:
            #os.rename(outfile.name,output) # atomic operation
            shutil.move(outfile.name,output) #It is not atomic if the files are on different filsystems.
        except:
            sys.stderr.write("[E::main] Error: failed to save the profile\n")
            sys.stderr.write("[E::main] you can find the file here:\n"+outfile.name+"\n")
            sys.exit(1)



########################
def main(argv=None):
    if(not argv):
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser(description='This script sum the mOTUs reads and calulates the LGs abundances ', add_help=True)
    parser.add_argument('infile', action="store", help='Infile represents a vector of relative abundace of the mOTUs')
    parser.add_argument('--outfile','-o', action="store", default="", dest='outfile', help='outfile for resulting table. If it is not set, then it will be print in standard output')
    parser.add_argument('--databaseDir', action='store', default="", dest='databaseDir', help='directory of the database')
    parser.add_argument('--onlySpecI', '-s', action='store_true', default=False, dest='onlySpecI', help='Set if you want to profile only specI (mOTU-LGs will go in -1)')
    parser.add_argument('--specI_ID', '-i', action='store_true', default=False, dest='specI_ID', help='Visualize the specI ID instead of the species name')
    parser.add_argument('--cutoff', '-c', action='store', dest='cutoff', default=2, type=int, help="minimum number of genes different from zero for deciding that a LG is present.")
    parser.add_argument('--sampleName', '-sn', action="store", dest='sampleName', default="", help='sample name for the current mapping')
    parser.add_argument('--taxonomic_level', action="store", default="species", dest='taxonomic_level', help='Taxonomic level for the profiling')
    parser.add_argument('--BIOM_output', action="store_true", default=False, dest='BIOM_output', help='print the result in BIOM format')
    parser.add_argument('-c', action='store_true', default=False, dest='print_rel_ab', help='print result as counts instead of relative abundance')
    parser.add_argument('-e', action='store_true', default=False, dest='print_NCBI_id', help='print NCBI id')
    parser.add_argument('--version', '-v', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()

    LGs_map = args.databaseDir+"mOTU-LG.map.tsv"
    LGs_map_l = args.databaseDir+"mOTU-LG.map.line.tsv"
    specI_taxonomy = args.databaseDir+"specI.taxonomy"
    mOTULG_taxonomy = args.databaseDir+"mOTULG.taxonomy"

    if (args.cutoff > 10) or (args.cutoff < 1):
        sys.stderr.write("[E::calc_motu] Error: cutoff should be between 1 and 10.\n")
        sys.exit(1)

    profile_mode = False # when using motu profile, this is set to True
    input_for_profile = ""

    calculate_abundance(args.infile, LGs_map, LGs_map_l, specI_taxonomy, mOTULG_taxonomy, args.outfile, args.cutoff, args.onlySpecI, args.sampleName, args.taxonomic_level, args.BIOM_output, profile_mode,input_for_profile, args.print_NCBI_id, args.print_rel_ab)

    return 0        # success



if __name__ == '__main__':
    status = main()
    sys.exit(status)
