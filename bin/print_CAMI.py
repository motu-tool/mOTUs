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

def calculate_abundance(infile, LGs_map, LGs_map_l, output, cutoff, onlySpecI, sampleName, profile_mode,input_for_profile, mgc_table_header,version_map_lgs,motu_version_tool,verbose,motu_call,git_commit_id,version_tool,CAMI_file,type_output,renormalise_cami,remove_strain_cami,remove_cami_comments):

    # load the CAMI annotation for the mOTUs
    CAMI_annotation = dict()
    try:
        location = open(CAMI_file,'r')
        # date is first line
        date_NCBI_dump = location.readline().rstrip()
        for line in location:
            vals = line.rstrip().split("\t")
            CAMI_annotation[vals[0]] = vals[1:]
    except:
        sys.stderr.write("[E::calc_motu] Error loading CAMI taxonomy file\n")
        sys.exit(1)


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
        if j != 'unassigned':
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
        else: # what to do with -1(unassigned): take the average of the different COGs values
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
    rel_ab_is_rounded = False
    rel_ab_LGs_rel = dict()
    s = sum(rel_ab_LGs.values())
    if s != 0:
        for j in list_LGs:
            rel_ab_LGs_rel[j] = (float(rel_ab_LGs[j])/s)*100
    else:
        for j in list_LGs:
            rel_ab_LGs_rel[j] = 0
        if verbose>1: sys.stderr.write(" [W::calc_motu] Warning: The relative abundance is 0 for all the mOTUs\n")

    # keep only specI
    if onlySpecI:
        rel_ab_LGs_rel_temp = dict(rel_ab_LGs_rel)
        rel_ab_LGs_rel = dict()
        value_minus1 = 0
        for j in list_LGs:
            if j == 'unassigned':
                value_minus1 = value_minus1 + rel_ab_LGs_rel_temp[j]
            else:
                type_c = j.split("_")[0]
                if (type_c == 'meta'):
                    value_minus1 = value_minus1 + rel_ab_LGs_rel_temp[j]
                else:
                    rel_ab_LGs_rel[j] = rel_ab_LGs_rel_temp[j]
        rel_ab_LGs_rel['unassigned'] = value_minus1

    # select only the mOTUs that are different from zero -----------------------
    mOTUs_to_print = dict()
    for i in rel_ab_LGs_rel:
        if rel_ab_LGs_rel[i] > 0:
            mOTUs_to_print[i] = rel_ab_LGs_rel[i]

    # --------------------------------------------------------------------------
    # merge at the different taxonomy levels -----------------------------------
    lines_to_print = dict()
    percentage_to_print = dict()
    for k in mOTUs_to_print:
        if k != "unassigned":
            for base in range(0,27,4):
                current_tax_id = CAMI_annotation[k][base+0]
                if current_tax_id != "NA":
                    if current_tax_id in lines_to_print:
                        percentage_to_print[current_tax_id] = percentage_to_print[current_tax_id] + mOTUs_to_print[k]
                    else:
                        percentage_to_print[current_tax_id] = mOTUs_to_print[k]
                        lines_to_print[current_tax_id] = CAMI_annotation[k][base+0]+"\t"+CAMI_annotation[k][base+1]+"\t"+CAMI_annotation[k][base+2]+"\t"+CAMI_annotation[k][base+3]+"\t"


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

    # print header
    if not remove_cami_comments:
        outfile.write("# Taxonomic Profiling Output\n")
        outfile.write(git_commit_id+" | "+map_lgs_header[1:]+"\n")
        outfile.write("# call: "+motu_call+"\n\n")


    outfile.write("@SampleID: "+sample_id_header+"\n")
    outfile.write("@Version:0.9.3\n")
    if remove_strain_cami == True:
        outfile.write("@Ranks:superkingdom|phylum|class|order|family|genus|species\n")
    else:
        outfile.write("@Ranks:superkingdom|phylum|class|order|family|genus|species|strain\n")
    outfile.write("@TaxonomyID: "+date_NCBI_dump+"\n")
    outfile.write("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n")

    # print the values
    tax_levels = ["superkingdom","phylum","class","order","family","genus","species"]
    # if recall, we have to check that there are no duplicates =================
    # otherwise this can happen:
    # 1328	species	2|1239|91061|186826|1300|1301|1328	Bacteria|Firmicutes|Bacilli|Lactobacillales|Streptococcaceae|Streptococcus|Streptococcus anginosus	3,285633475
    # 1328	species	2|1239|91061|186826|1300|1301|(1328/1338)	Bacteria|Firmicutes|Bacilli|Lactobacillales|Streptococcaceae|Streptococcus|(Streptococcus anginosus/Streptococcus intermedius)	14,64201763
    # 1338	species	2|1239|91061|186826|1300|1301|(1328/1338)	Bacteria|Firmicutes|Bacilli|Lactobacillales|Streptococcaceae|Streptococcus|(Streptococcus anginosus/Streptococcus intermedius)	14,64201763
    # 1338	species	2|1239|91061|186826|1300|1301|(1338/76860)	Bacteria|Firmicutes|Bacilli|Lactobacillales|Streptococcaceae|Streptococcus|(Streptococcus intermedius/Streptococcus constellatus)	6,063508333
    if type_output == "recall":
        for i in tax_levels:
            all_lines = list()
            for taxa in lines_to_print:
                if lines_to_print[taxa].split("\t")[1] == i:
                    # repleace "NA" with ""
                    lines_to_print[taxa] = lines_to_print[taxa].replace("NA", "")
                    if not (lines_to_print[taxa].startswith("(")): # print the one without parenthesis normally
                        all_lines.append(lines_to_print[taxa]+str(percentage_to_print[taxa])+"\n")
                    else: # for the one with parenthesis, we split the relative abundance equally
                        line_values = lines_to_print[taxa].split("\t")
                        all_taxa = line_values[0][1:-1].split("/")
                        splitted_perc = percentage_to_print[taxa]/len(all_taxa)
                        for t in all_taxa:
                            all_lines.append(t+"\t"+("\t".join(line_values[1:]))+str(splitted_perc)+"\n")
            # now find all unique tax ids (pos 1)
            merged_lines = dict()
            for line in all_lines:
                t_id = line.split("\t")[0]
                if t_id in merged_lines:
                    vals = line.rstrip().split("\t")
                    vals_dict = merged_lines[t_id].rstrip().split("\t")
                    vals_dict[4] = str( float(vals_dict[4]) + float(vals[4]) )
                    merged_lines[t_id] = "\t".join(vals_dict) + "\n"
                else:
                    merged_lines[t_id] = line
            # renormalise if needed ----
            if renormalise_cami:
                sum_A = 0
                for ll in merged_lines:
                    sum_A = sum_A + float(merged_lines[ll].rstrip().split("\t")[4])
                if sum_A != 0:
                    for line in merged_lines:
                        vals_dict = merged_lines[line].rstrip().split("\t")
                        vals_dict[4] = str( float(vals_dict[4]) / sum_A * 100 )
                        merged_lines[line] = "\t".join(vals_dict) + "\n"

            # print all unique lines
            for ll in merged_lines:
                outfile.write(merged_lines[ll])


    # if no recall: ============================================================
    if type_output != "recall":
        for i in tax_levels:
            merged_lines = dict()
            for taxa in lines_to_print:
                if lines_to_print[taxa].split("\t")[1] == i:
                    # repleace "NA" with ""
                    lines_to_print[taxa] = lines_to_print[taxa].replace("NA", "")
                    # decide how to print it
                    if type_output == "parenthesis": # we print the result with parenthesis - this is not tecnically in CAMI format
                        dummy_v = lines_to_print[taxa]+str(percentage_to_print[taxa])+"\n"
                        merged_lines[dummy_v] = dummy_v
                    if type_output == "precision":
                        if not (lines_to_print[taxa].startswith("(")): # print only the ones without parenthesis
                            dummy_v = lines_to_print[taxa]+str(percentage_to_print[taxa])+"\n"
                            merged_lines[dummy_v] = dummy_v
            # renormalise if needed ----
            if renormalise_cami:
                sum_A = 0
                for ll in merged_lines:
                    sum_A = sum_A + float(merged_lines[ll].rstrip().split("\t")[4])
                if sum_A != 0:
                    for line in merged_lines:
                        vals_dict = merged_lines[line].rstrip().split("\t")
                        vals_dict[4] = str( float(vals_dict[4]) / sum_A * 100 )
                        merged_lines[line] = "\t".join(vals_dict) + "\n"

            # print all unique lines
            for ll in merged_lines:
                outfile.write(merged_lines[ll])




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

    sys.exit(0)
