#!/usr/bin/env python

import sys
import re

# to execute for 1 line
def check_one_line(line,min_perc_id,min_length_align,min_perc_cover,res_stdout):
    if len(line) == 0:
        if not res_stdout:
            return "None"
    if line[0]=="@": # reprint the header
        if res_stdout:
            sys.stdout.write(line)
        else:
            return line
    else:
        arr = line.split("\t")
        if not((arr[1] == '4') or (arr[2] == '*') or (arr[5] == '*')):
            len_seq = 0

            min_query_al = 0

            tott = 0

            flag2 = False
            tot = 0
            for n,i in (re.findall('(\d+)([IDMSHX=])', arr[5])):
                tott = tott + int(n)
                if i == 'M' or i == 'D' or i == 'X' or i == '=':
                    tot = tot + int(n)
                if i != 'I' and i != 'S' and i!= 'H':
                    len_seq = len_seq + int(n)
                if i != 'H' and i != 'S' and i!= 'D':
                    min_query_al = min_query_al + int(n)

            if tot >= float(min_length_align):
                flag2 = True

            flag1 = False
            if ((int(arr[11].split(":")[2])/float(len_seq)) < (1-float(min_perc_id))): # TODO: here we assume that the NM value is in position 11, it should be checked
                flag1 = True

            # min. percent of the query that must be aligned, between 0 and 100 (required)
            flag3 = False
            if (float(min_query_al)/float(tott)) >= (float(min_perc_cover)/100 ):
                flag3 = True



            if flag1 and flag2 and flag3:
                if res_stdout:
                    sys.stdout.write(line)
                else:
                    return line
            else:
                if not res_stdout:
                    return "None"
    if not res_stdout:
        return "None"





# run for all lines
def run_all_lines(min_perc_id,min_length_align,min_perc_cover,input_t):
    if input_t == "None":
        for line in sys.stdin:
            check_one_line(line,min_perc_id,min_length_align,min_perc_cover,True)
    else:
        res = list() # if we dont use stdin, then we dont use also stdout
        for line in input_t:
            rr = check_one_line(line,min_perc_id,min_length_align,min_perc_cover,False)
            if rr != "None":
                res.append(rr)
        return res


def main(argv=None):
    #min. percent identity of alignment, between 0 and 1; requires NM field to be present (required)
    min_perc_id = sys.argv[1]
    #min. length of alignment (required)
    min_length_align = sys.argv[2]
    #min. percent of the query that must be aligned, between 0 and 100 (required)
    min_perc_cover = sys.argv[3]
    #type of input - if it is null then we choose stdin, otherwise it is the input itself
    input_t = sys.argv[4]

    run_all_lines(min_perc_id,min_length_align,min_perc_cover,input_t)



if __name__ == '__main__':
	status = main()
	sys.exit(status)
