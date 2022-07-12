#!/usr/bin/env python

# ============================================================================ #
# convert_long_reads.py: Take long reads and split them into short reads that
#                        can be analysed by mOTUs
#
# ============================================================================ #


import gzip
import sys

possible_nucleotides = ["A","C","G","T","U","W","S","M","K","R","Y","B","D","H","V","N"]

def split_read(read_name, read_seq, split_len, min_len, quality):
    count = 0
    # remove first character of the read name (either ">" or "@")
    read_name = read_name[1:]
    # check
    ll = len(read_seq)
    if ll < min_len:
        return ""

    # from now on we have long enough reads
    res = ""
    if ll < split_len:
        # between 100 and 300 (or split_len)
        count = count + 1
        res = res + "@"+read_name+"_"+str(count) + "\n"
        res = res + read_seq + "\n"
        res = res + "+\n"
        res = res + quality*ll + "\n"
        return res
    else:
        # we need to split the reads
        n_splits = round(ll/split_len)
        length_split = round(ll/n_splits)
        chunks = [read_seq[i:i+length_split] for i in range(0, len(read_seq), length_split)]
        for i in chunks:
            count = count + 1
            res = res + "@"+read_name+"_"+str(count) + "\n"
            res = res + i + "\n"
            res = res + "+\n"
            res = res + quality*len(i) + "\n"
        return res


def convert_long_reads(path_original, path_converted, split_len = 300, min_len= 50, quality = "D", gz_out = True, verbose = 3, log_ = None):
    # - path_original, file with the long reads (can be fasta or fastq) (can be
    #                  gzipped)
    # - path_converted, where to save the converted reads
    # - split_len, length to which spiltting the reads (if they are shorter than
    #              this they will stay the same)
    # - min_len, minimum length of the reads (if shorter they are removed)
    # - quality, quality that we use for the output fastq files, default is "D"
    #            which is 35
    # - gz_out, it's a BOOL that says if the output should be gzipped or not

    # set up log
    global log
    log = log_
    # ----------------------

    # First we check that the file exist ---------------------------------------
    # and it is a fasta or fastq
    try:
        is_gz = False
        if path_original.endswith(".gz"):
            is_gz = True
            f = gzip.open(path_original,'rb')
            line1 = f.readline().rstrip().decode(encoding='UTF-8')
            line2 = f.readline().rstrip().decode(encoding='UTF-8')
            line3 = f.readline().rstrip().decode(encoding='UTF-8')
        else:
            f = open(path_original,'r')
            line1 = f.readline().rstrip()
            line2 = f.readline().rstrip()
            line3 = f.readline().rstrip()
        f.close()
        # check if it is fasta or fastq
        file_type = "unknown"
        if line1.startswith(">"):
            if line2[0] in possible_nucleotides:
                file_type = "fasta"
        if line1.startswith("@"):
            if line3.startswith("+"):
                file_type = "fastq"
    except Exception as e:
        log.print_error("Cannot open the input file", exit = False)
        if verbose > 2:
            sys.stderr.write(str(e)+"\n")
        sys.exit(1)

    if file_type == "unknown":
        log.print_error("Cannot recognise file type (not fastq or fasta)")

    # Transform the reads ------------------------------------------------------
    # open file
    if is_gz:
        f = gzip.open(path_original,'rb')
    else:
        f = open(path_original,'r')
    # open file to save
    if gz_out:
        w = gzip.open(path_converted,'wb')
    else:
        w = open(path_converted,'w')

    # go through fastq file
    if file_type == "fastq":
        for line_b in f:
            # read the 4 lines
            head = line_b.rstrip()
            seq = f.readline().rstrip()
            f.readline()
            f.readline()
            if is_gz:
                splitted_reads = split_read(head.decode(encoding='UTF-8'), seq.decode(encoding='UTF-8'), split_len, min_len, quality)
            else:
                splitted_reads = split_read(head, seq, split_len, min_len, quality)
            # save to file
            if gz_out:
                w.write(splitted_reads.encode())
            else:
                w.write(splitted_reads)

    if file_type == "fasta":
        head = ""
        for line_b in f:
            if is_gz:
                line = line_b.decode(encoding='UTF-8')
            else:
                line = line_b
            #
            if line.startswith(">"):
                if head != "":
                    splitted_reads = split_read(head, seq, split_len, min_len, quality)
                    # save to file
                    if gz_out:
                        w.write(splitted_reads.encode())
                    else:
                        w.write(splitted_reads)
                # make new head
                head = line.rstrip()
                seq = ""
            else:
                seq = seq + line.rstrip()
        # print the last read
        splitted_reads = split_read(head, seq, split_len, min_len, quality)
        # save to file
        if gz_out:
            w.write(splitted_reads.encode())
        else:
            w.write(splitted_reads)


    f.close()
    w.close()
