import sys
import datetime

def get_timestamp():
    return '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())


def print_log(text):
    sys.stderr.write(get_timestamp()+" "+text+"\n")
    sys.stderr.flush()

def print_message(text):
    sys.stderr.write("   "+text+"\n")
    sys.stderr.flush()

def print_message_execution(text):
    sys.stderr.write("   "+text+"\n")
    sys.stderr.flush()

def print_message_time(text):
    sys.stderr.write("  "+get_timestamp()+" "+text+"\n")
    sys.stderr.flush()

def print_error(text, exit = True):
    sys.stderr.write("   "+"Error: "+text+"\n")
    sys.stderr.flush()
    if exit:
        sys.exit(1)

def print_warning(text):
    sys.stderr.write("   "+"Warning: "+text+"\n")
    sys.stderr.flush()









################################################################################
################################################################################
# print menus
def msg(version_tool):
    str_msg = '''
\00
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
Version: '''+version_tool+'''
Reference: Milanese et al. Microbial abundance, activity and population genomic profiling with mOTUs2. Nature Communications (2019). doi: 10.1038/s41467-019-08844-4

Usage: motus <command> [options]

Command:
 -- Taxonomic profiling
      profile     Perform taxonomic profiling (map_tax + calc_mgc + calc_motu) in a single step
      merge       Merge several taxonomic profiling results into one table

      map_tax     Map reads to the marker gene database
      calc_mgc    Calculate marker gene cluster (MGC) abundance
      calc_motu   Summarize MGC abundances into a mOTU abundance table

 -- SNV calling
      map_snv     Map reads to the marker gene database for SNV calling
      snv_call    Generate SNV profiles (using metaSNV)

Type motus <command> to print the help menu for a specific command
        '''
    return str_msg






# ------------------------------------------------------------------------------
def print_menu_profile():
    sys.stderr.write("\n")
    sys.stderr.write("Usage: motus profile [options]\n\n")
    sys.stderr.write("Input options:\n")
    sys.stderr.write("   -f  FILE[,FILE]  input file(s) for reads in forward orientation, fastq formatted\n")
    sys.stderr.write("   -r  FILE[,FILE]  input file(s) for reads in reverse orientation, fastq formatted\n")
    sys.stderr.write("   -s  FILE[,FILE]  input file(s) for reads without mate, fastq formatted\n")
    sys.stderr.write("   -n  STR          sample name\n")
    sys.stderr.write("   -i  FILE[,FILE]  provide SAM or BAM input file(s)\n")
    sys.stderr.write("   -m  FILE         provide a mgc reads count file\n")
    sys.stderr.write("   -db DIR          provide a database directory\n\n")
    sys.stderr.write("Output options:\n")
    sys.stderr.write("   -o  FILE         output file name [stdout]\n")
    sys.stderr.write("   -I  FILE         save the result of BWA in BAM format (intermediate step) [None]\n")
    sys.stderr.write("   -M  FILE         save the mgc reads count (intermediate step) [None]\n")
    sys.stderr.write("   -e               profile only reference species (ref_mOTUs)\n")
    sys.stderr.write("   -c               print result as counts instead of relative abundances\n")
    sys.stderr.write("   -p               print NCBI id\n")
    sys.stderr.write("   -u               print the full name of the species\n")
    sys.stderr.write("   -q               print the full rank taxonomy\n")
    sys.stderr.write("   -B               print result in BIOM format\n")
    sys.stderr.write("   -C  STR          print result in CAMI format (BioBoxes format 0.9.1)\n")
    sys.stderr.write("                    Values: [precision, recall, parenthesis]\n")
    sys.stderr.write("   -A               print all taxonomy levels together\n")
    sys.stderr.write("   -k  STR          taxonomic level [mOTU]\n")
    sys.stderr.write("                    Values: [kingdom, phylum, class, order, family, genus, mOTU]\n\n")
    sys.stderr.write("Algorithm options:\n")
    sys.stderr.write("   -g  INT          number of marker genes cutoff: 1=higher recall, 10=higher precision [3]\n")
    sys.stderr.write("   -l  INT          min. length of alignment for the reads (number of nucleotides) [75]\n")
    sys.stderr.write("   -t  INT          number of threads [1]\n")
    sys.stderr.write("   -v  INT          verbose level: 1=error, 2=warning, 3=message, 4+=debugging [3]\n")
    sys.stderr.write("   -y  STR          type of read counts [insert.scaled_counts]\n")
    sys.stderr.write("                    Values: [base.coverage, insert.raw_counts, insert.scaled_counts]\n\n")

# ------------------------------------------------------------------------------
def print_menu_map_snv():
    sys.stderr.write("\n")
    sys.stderr.write("Usage: motus map_snv [options]\n\n")
    sys.stderr.write("Input options:\n")
    sys.stderr.write("   -f  FILE[,FILE]  input file(s) for reads in forward orientation, fastq formatted\n")
    sys.stderr.write("   -r  FILE[,FILE]  input file(s) for reads in reverse orientation, fastq formatted\n")
    sys.stderr.write("   -s  FILE[,FILE]  input file(s) for reads without mate, fastq formatted\n")
    sys.stderr.write("   -db DIR          provide a database directory\n\n")
    sys.stderr.write("Output options:\n")
    sys.stderr.write("   -o  FILE         output BAM file name [stdout]\n\n")
    sys.stderr.write("Algorithm options:\n")
    sys.stderr.write("   -l  INT          min. length of alignment for the reads (number of nucleotides) [75]\n")
    #sys.stderr.write("   -CC INT          number of cores for bwa (max one per lane) [1]\n")
    sys.stderr.write("   -t  INT          number of threads [1]\n")
    sys.stderr.write("   -v  INT          verbose level: 1=error, 2=warning, 3=message, 4+=debugging [3]\n\n")

# ------------------------------------------------------------------------------
def print_menu_snv_call():
    sys.stderr.write("\n")
    sys.stderr.write("Usage: motus snv_call -d Directory -o Directory [options]\n\n")
    sys.stderr.write("Input options:\n")
    sys.stderr.write("   -d  DIR     Call metaSNV on all BAM files in the directory. [Mandatory]\n")
    sys.stderr.write("   -fb FLOAT   Coverage breadth: minimal horizontal genome coverage percentage per sample per species. Default=80.0\n")
    sys.stderr.write("   -fd FLOAT   Coverage depth: minimal average vertical genome coverage per sample per species. Default=5.0\n")
    sys.stderr.write("   -fm INT     Minimum number of samples per species. Default=2\n")
    sys.stderr.write("   -fp FLOAT   FILTERING STEP II: Required proportion of informative samples (coverage non-zero) per position. Default=0.90\n")
    sys.stderr.write("   -fc FLOAT   FILTERING STEP II: Minimum coverage per position per sample per species. Default=5.0\n")
    sys.stderr.write("   -db DIR     Provide a database directory\n\n")
    sys.stderr.write("Output options:\n")
    sys.stderr.write("   -o  DIR     Output directory. Will fail if already exists. [Mandatory]\n")
    sys.stderr.write("   -K          Keep all the directories produced by metaSNV. Default is to remove cov, distances, filtered, snpCaller\n\n")
    sys.stderr.write("Algorithm options:\n")
    sys.stderr.write("   -v INT      Verbose level: 1=error, 2=warning, 3=message, 4+=debugging. Default=3\n\n")

# ------------------------------------------------------------------------------
def print_menu_bwa():
    sys.stderr.write("\n")
    sys.stderr.write("Usage: motus map_tax [options]\n\n")
    sys.stderr.write("Input options:\n")
    sys.stderr.write("   -f  FILE[,FILE]  input file(s) for reads in forward orientation, fastq formatted\n")
    sys.stderr.write("   -r  FILE[,FILE]  input file(s) for reads in reverse orientation, fastq formatted\n")
    sys.stderr.write("   -s  FILE[,FILE]  input file(s) for reads without mate, fastq formatted\n")
    sys.stderr.write("   -db DIR          provide a database directory\n\n")
    sys.stderr.write("Output options:\n")
    sys.stderr.write("   -o  FILE         output file name [stdout]\n")
    sys.stderr.write("   -b               save the result of BWA in BAM format\n\n")
    sys.stderr.write("Algorithm options:\n")
    sys.stderr.write("   -l  INT          min. length of alignment for the reads (number of nucleotides) [75]\n")
    #sys.stderr.write("   -CC INT          number of cores for bwa (max one per lane) [1]\n")
    sys.stderr.write("   -t  INT          number of threads [1]\n")
    sys.stderr.write("   -v  INT          verbose level: 1=error, 2=warning, 3=message, 4+=debugging [3]\n\n")

# ------------------------------------------------------------------------------
def print_menu_map_genes():
    sys.stderr.write("\n")
    sys.stderr.write("Usage: motus calc_mgc [options]\n\n")
    sys.stderr.write("Input options:\n")
    sys.stderr.write("   -n  STR          sample name\n")
    sys.stderr.write("   -i  FILE[,FILE]  provide a SAM or BAM input file (or list of files) output of motus map_tax\n")
    sys.stderr.write("   -db DIR          provide a database directory\n\n")
    sys.stderr.write("Output options:\n")
    sys.stderr.write("   -o  FILE         output file name [stdout]\n\n")
    sys.stderr.write("Algorithm options:\n")
    sys.stderr.write("   -l  INT          min. length of alignment for the reads (number of nucleotides) [75]\n")
    #sys.stderr.write("   -L FLOAT        min. length of alignment for the reads (percentage of the average read length) [None]\n")
    sys.stderr.write("   -v  INT          verbose level: 1=error, 2=warning, 3=message, 4+=debugging [3]\n")
    sys.stderr.write("   -y  STR          type of read counts [insert.scaled_counts]\n")
    sys.stderr.write("                    Values: [base.coverage, insert.raw_counts, insert.scaled_counts]\n\n")

# ------------------------------------------------------------------------------
def print_menu_map_lgs():
    sys.stderr.write("\n")
    sys.stderr.write("Usage: motus calc_motu [options]\n\n")
    sys.stderr.write("Input options:\n")
    sys.stderr.write("   -n  STR   sample name\n")
    sys.stderr.write("   -i  FILE  provide the mgc abundance table (output of motus calc_mgc)\n")
    sys.stderr.write("   -db DIR   provide a database directory\n\n")
    sys.stderr.write("Output options:\n")
    sys.stderr.write("   -o  FILE  output file name [stdout]\n")
    sys.stderr.write("   -e        profile only reference species (ref_mOTUs)\n")
    sys.stderr.write("   -B        print result in BIOM format\n")
    sys.stderr.write("   -C  STR   print result in CAMI format (BioBoxes format 0.9.1)\n")
    sys.stderr.write("             Values: [precision, recall, parenthesis]\n")
    sys.stderr.write("   -A        print all taxonomy levels together\n")
    sys.stderr.write("   -c        print result as counts instead of relative abundances\n")
    sys.stderr.write("   -p        print NCBI id\n")
    sys.stderr.write("   -u        print the full name of the species\n")
    sys.stderr.write("   -q        print the full rank taxonomy\n")
    sys.stderr.write("   -k  STR   taxonomic level [mOTU]\n")
    sys.stderr.write("             Values: [kingdom, phylum, class, order, family, genus, mOTU]\n\n")
    sys.stderr.write("Algorithm options:\n")
    sys.stderr.write("   -g  INT   number of marker genes cutoff: 1=higher recall, 10=higher precision [3]\n")
    sys.stderr.write("   -v  INT   verbose level: 1=error, 2=warning, 3=message, 4+=debugging [3]\n\n")

# ------------------------------------------------------------------------------
def print_menu_append():
    sys.stderr.write("\n")
    sys.stderr.write("Usage: motus merge [options]\n\n")
    sys.stderr.write("Input options:\n")
    sys.stderr.write("   -i FILE[,FILE] list of mOTU profiles to merge (comma separated)\n")
    sys.stderr.write("   -d DIR         merge all the files in the directory DIR\n")
    sys.stderr.write("   -a STR[,STR]   add pre-computed profiles from different environmental samples\n                  Values: [all, air, bioreactor, bee, cat,\n                  cattle, chicken, dog, fish, freshwater, human,\n                  marine, mouse, pig, sheep, soil, termite, wastewater]\n\n")
    sys.stderr.write("Output options:\n")
    sys.stderr.write("   -o FILE        output file name [stdout]\n")
    sys.stderr.write("   -B             print result in BIOM format\n\n")
    sys.stderr.write("Algorithm options:\n")
    sys.stderr.write("   -v INT         verbose level: 1=error, 2=warning, 3=message, 4+=debugging [3]\n\n")
