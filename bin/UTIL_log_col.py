import sys
import datetime

END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
RED = '\033[31m'
GREEN = '\033[32m'
MAGENTA = '\033[35m'
YELLOW = '\033[93m'
CYAN = "\033[36m"
BLUE = "\033[34m"
DIM = '\033[2m'


def colour(text, text_colour):
    bold_text = 'bold' in text_colour
    text_colour = text_colour.replace('bold', '')
    underline_text = 'underline' in text_colour
    text_colour = text_colour.replace('underline', '')
    text_colour = text_colour.replace('_', '')
    text_colour = text_colour.replace(' ', '')
    text_colour = text_colour.lower()
    if 'red' in text_colour:
        coloured_text = RED
    elif 'green' in text_colour:
        coloured_text = GREEN
    elif 'yellow' in text_colour:
        coloured_text = YELLOW
    elif 'blue' in text_colour:
        coloured_text = BLUE
    elif 'cyan' in text_colour:
        coloured_text = CYAN
    elif 'magenta' in text_colour:
        coloured_text = MAGENTA
    elif 'dim' in text_colour:
        coloured_text = DIM
    else:
        coloured_text = ''
    if bold_text:
        coloured_text += BOLD
    if underline_text:
        coloured_text += UNDERLINE
    if not coloured_text:
        return text
    coloured_text += text + END_FORMATTING
    return coloured_text



def get_timestamp():
    return '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())


def print_log(text):
    sys.stderr.write(colour(get_timestamp(),"bold_blue")+" "+colour(text,"bold_blue")+"\n")
    sys.stderr.flush()

def print_message(text):
    sys.stderr.write("   "+text+"\n")
    sys.stderr.flush()

def print_message_time(text):
    sys.stderr.write("   "+colour(get_timestamp(),"bold")+" "+text+"\n")
    sys.stderr.flush()

def print_message_execution(text):
    sys.stderr.write("   "+colour(text,"cyan")+"\n")
    sys.stderr.flush()

def print_error(text, exit = True):
    sys.stderr.write(colour("Error: ","red_bold")+colour(text,"red")+"\n")
    sys.stderr.flush()
    if exit:
        sys.exit(1)

def print_warning(text):
    sys.stderr.write("   "+colour("Warning: ","magenta_bold")+colour(text,"magenta")+"\n")
    sys.stderr.flush()









################################################################################
################################################################################
# print menus
def msg(version_tool):
    str_msg = '''
\00
'''+colour("Program:","blue_bold")+''' motus - a tool for marker gene-based OTU (mOTU) profiling
'''+colour("Version: ","blue_bold")+version_tool+'''
'''+colour("Reference:","blue_bold")+''' Milanese et al. Microbial abundance, activity and population genomic profiling with mOTUs2. Nature Communications (2019). doi: 10.1038/s41467-019-08844-4

'''+colour("Usage:","blue_bold")+''' motus <command> [options]

'''+colour("Command:","blue_bold")+'''
'''+colour("-- Taxonomic profiling","bold")+'''
'''+colour("      profile","cyan")+'''     Perform taxonomic profiling (map_tax + calc_mgc + calc_motu) in a single step
'''+colour("      merge","cyan")+'''       Merge several taxonomic profiling results into one table

'''+colour("      map_tax","cyan")+'''     Map reads to the marker gene database
'''+colour("      calc_mgc","cyan")+'''    Calculate marker gene cluster (MGC) abundance
'''+colour("      calc_motu","cyan")+'''   Summarize MGC abundances into a mOTU abundance table

'''+colour("-- SNV calling","bold")+'''
'''+colour("      map_snv","cyan")+'''     Map reads to the marker gene database for SNV calling
'''+colour("      snv_call","cyan")+'''    Generate SNV profiles (using metaSNV)

Type motus <command> to print the help menu for a specific command
        '''
    return str_msg




# ------------------------------------------------------------------------------
def print_menu_profile():
    sys.stderr.write("\n")
    sys.stderr.write(colour("Usage:","blue_bold")+" motus profile [options]\n\n")
    sys.stderr.write(colour("Input options:\n","bold"))
    sys.stderr.write("   "+colour("-f ","blue_bold")+" FILE[,FILE]  input file(s) for reads in forward orientation, fastq formatted\n")
    sys.stderr.write("   "+colour("-r ","blue_bold")+" FILE[,FILE]  input file(s) for reads in reverse orientation, fastq formatted\n")
    sys.stderr.write("   "+colour("-s ","blue_bold")+" FILE[,FILE]  input file(s) for reads without mate, fastq formatted\n")
    sys.stderr.write("   "+colour("-n ","blue_bold")+" STR          sample name "+colour("['unnamed sample']","magenta")+"\n")
    sys.stderr.write("   "+colour("-i ","blue_bold")+" FILE[,FILE]  provide SAM or BAM input file(s)\n")
    sys.stderr.write("   "+colour("-m ","blue_bold")+" FILE         provide a mgc reads count file\n")
    sys.stderr.write("   "+colour("-db ","blue_bold")+"DIR          provide a database directory\n\n")
    sys.stderr.write(colour("Output options:\n","bold"))
    sys.stderr.write("   "+colour("-o ","blue_bold")+" FILE         output file name "+colour("[stdout]","magenta")+"\n")
    sys.stderr.write("   "+colour("-I ","blue_bold")+" FILE         save the result of BWA in BAM format (intermediate step)\n")
    sys.stderr.write("   "+colour("-M ","blue_bold")+" FILE         save the mgc reads count (intermediate step)\n")
    sys.stderr.write("   "+colour("-e ","blue_bold")+"              profile only reference species (ref_mOTUs)\n")
    sys.stderr.write("   "+colour("-c ","blue_bold")+"              print result as counts instead of relative abundances\n")
    sys.stderr.write("   "+colour("-p ","blue_bold")+"              print NCBI id\n")
    sys.stderr.write("   "+colour("-u ","blue_bold")+"              print the full name of the species\n")
    sys.stderr.write("   "+colour("-q ","blue_bold")+"              print the full rank taxonomy\n")
    sys.stderr.write("   "+colour("-B ","blue_bold")+"              print result in BIOM format\n")
    sys.stderr.write("   "+colour("-C ","blue_bold")+" STR          print result in CAMI format (BioBoxes format 0.9.1)\n")
    sys.stderr.write("                    Values: [precision, recall, parenthesis]\n")
    sys.stderr.write("   "+colour("-A ","blue_bold")+"              print all taxonomy levels together\n")
    sys.stderr.write("   "+colour("-k ","blue_bold")+" STR          taxonomic level "+colour("[mOTU]","magenta")+"\n")
    sys.stderr.write("                    Values: [kingdom, phylum, class, order, family, genus, mOTU]\n\n")
    sys.stderr.write(colour("Algorithm options:\n","bold"))
    sys.stderr.write("   "+colour("-g ","blue_bold")+" INT          number of marker genes cutoff: 1=higher recall, 10=higher precision "+colour("[3]","magenta")+"\n")
    sys.stderr.write("   "+colour("-l ","blue_bold")+" INT          min. length of alignment for the reads (number of nucleotides) "+colour("[75]","magenta")+"\n")
    sys.stderr.write("   "+colour("-t ","blue_bold")+" INT          number of threads "+colour("[1]","magenta")+"\n")
    sys.stderr.write("   "+colour("-v ","blue_bold")+" INT          verbose level: 1=error, 2=warning, 3=message, 4+=debugging "+colour("[3]","magenta")+"\n")
    sys.stderr.write("   "+colour("-y ","blue_bold")+" STR          type of read counts "+colour("[insert.scaled_counts]","magenta")+"\n")
    sys.stderr.write("                    Values: [base.coverage, insert.raw_counts, insert.scaled_counts]\n\n")

# ------------------------------------------------------------------------------
def print_menu_map_snv():
    sys.stderr.write("\n")
    sys.stderr.write(colour("Usage:","blue_bold")+" motus map_snv [options]\n\n")
    sys.stderr.write(colour("Input options:\n","bold"))
    sys.stderr.write("   "+colour("-f ","blue_bold")+" FILE[,FILE]  input file(s) for reads in forward orientation, fastq formatted\n")
    sys.stderr.write("   "+colour("-r ","blue_bold")+" FILE[,FILE]  input file(s) for reads in reverse orientation, fastq formatted\n")
    sys.stderr.write("   "+colour("-s ","blue_bold")+" FILE[,FILE]  input file(s) for reads without mate, fastq formatted\n")
    sys.stderr.write("   "+colour("-db ","blue_bold")+" DIR          provide a database directory\n\n")
    sys.stderr.write(colour("Output options:\n","bold"))
    sys.stderr.write("   "+colour("-o ","blue_bold")+" FILE         output BAM file name "+colour("[stdout]","magenta")+"\n\n")
    sys.stderr.write(colour("Algorithm options:\n","bold"))
    sys.stderr.write("   "+colour("-l ","blue_bold")+" INT          min. length of alignment for the reads (number of nucleotides) "+colour("[75]","magenta")+"\n")
    #sys.stderr.write("   -CC INT          number of cores for bwa (max one per lane) [1]\n")
    sys.stderr.write("   "+colour("-t ","blue_bold")+" INT          number of threads "+colour("[1]","magenta")+"\n")
    sys.stderr.write("   "+colour("-v ","blue_bold")+" INT          verbose level: 1=error, 2=warning, 3=message, 4+=debugging "+colour("[3]","magenta")+"\n\n")

# ------------------------------------------------------------------------------
def print_menu_snv_call():
    sys.stderr.write("\n")
    sys.stderr.write(colour("Usage:","blue_bold")+" motus snv_call "+colour("-d","blue_bold")+" <DIR> "+colour("-o","blue_bold")+" <DIR> [options]\n\n")
    sys.stderr.write(colour("Input options:\n","bold"))
    sys.stderr.write("   "+colour("-d ","blue_bold")+" DIR     Call metaSNV on all BAM files in the directory\n")
    sys.stderr.write("   "+colour("-fb ","blue_bold")+"FLOAT   Coverage breadth: minimal horizontal genome coverage percentage per sample per species "+colour("[80.0]","magenta")+"\n")
    sys.stderr.write("   "+colour("-fd ","blue_bold")+"FLOAT   Coverage depth: minimal average vertical genome coverage per sample per species "+colour("[5.0]","magenta")+"\n")
    sys.stderr.write("   "+colour("-fm ","blue_bold")+"INT     Minimum number of samples per species "+colour("[2]","magenta")+"\n")
    sys.stderr.write("   "+colour("-fp ","blue_bold")+"FLOAT   FILTERING STEP II: Required proportion of informative samples (coverage non-zero) per position "+colour("[0.90]","magenta")+"\n")
    sys.stderr.write("   "+colour("-fc ","blue_bold")+"FLOAT   FILTERING STEP II: Minimum coverage per position per sample per species "+colour("[5.0]","magenta")+"\n")
    sys.stderr.write("   "+colour("-db ","blue_bold")+"DIR     Provide a database directory\n\n")
    sys.stderr.write(colour("Output options:\n","bold"))
    sys.stderr.write("   "+colour("-o ","blue_bold")+" DIR     Output directory. Will fail if already exists\n")
    sys.stderr.write("   "+colour("-K ","blue_bold")+"         Keep all the directories produced by metaSNV. Default is to remove cov, distances, filtered, snpCaller\n\n")
    sys.stderr.write(colour("Algorithm options:\n","bold"))
    sys.stderr.write("   "+colour("-v ","blue_bold")+" INT      Verbose level: 1=error, 2=warning, 3=message, 4+=debugging "+colour("[3]","magenta")+"\n\n")

# ------------------------------------------------------------------------------
def print_menu_bwa():
    sys.stderr.write("\n")
    sys.stderr.write(colour("Usage:","blue_bold")+" motus map_tax [options]\n\n")
    sys.stderr.write(colour("Input options:\n","bold"))
    sys.stderr.write("   "+colour("-f ","blue_bold")+"  FILE[,FILE]  input file(s) for reads in forward orientation, fastq formatted\n")
    sys.stderr.write("   "+colour("-r ","blue_bold")+"  FILE[,FILE]  input file(s) for reads in reverse orientation, fastq formatted\n")
    sys.stderr.write("   "+colour("-s ","blue_bold")+"  FILE[,FILE]  input file(s) for reads without mate, fastq formatted\n")
    sys.stderr.write("   "+colour("-db ","blue_bold")+" DIR          provide a database directory\n\n")
    sys.stderr.write(colour("Output options:\n","bold"))
    sys.stderr.write("   "+colour("-o ","blue_bold")+"  FILE         output file name "+colour("[stdout]","magenta")+"\n")
    sys.stderr.write("   "+colour("-b ","blue_bold")+"               save the result of BWA in BAM format\n\n")
    sys.stderr.write(colour("Algorithm options:\n","bold"))
    sys.stderr.write("   "+colour("-l ","blue_bold")+"  INT          min. length of alignment for the reads (number of nucleotides) "+colour("[75]","magenta")+"\n")
    #sys.stderr.write("   -CC INT          number of cores for bwa (max one per lane) [1]\n")
    sys.stderr.write("   "+colour("-t ","blue_bold")+"  INT          number of threads "+colour("[1]","magenta")+"\n")
    sys.stderr.write("   "+colour("-v ","blue_bold")+"  INT          verbose level: 1=error, 2=warning, 3=message, 4+=debugging "+colour("[3]","magenta")+"\n\n")

# ------------------------------------------------------------------------------
def print_menu_map_genes():
    sys.stderr.write("\n")
    sys.stderr.write(colour("Usage:","blue_bold")+" motus calc_mgc [options]\n\n")
    sys.stderr.write(colour("Input options:\n","bold"))
    sys.stderr.write("   "+colour("-n ","blue_bold")+" STR          sample name\n")
    sys.stderr.write("   "+colour("-i ","blue_bold")+" FILE[,FILE]  provide a SAM or BAM input file (or list of files) output of motus map_tax\n")
    sys.stderr.write("   "+colour("-db ","blue_bold")+"DIR          provide a database directory\n\n")
    sys.stderr.write(colour("Output options:\n","bold"))
    sys.stderr.write("   "+colour("-o ","blue_bold")+" FILE         output file name "+colour("[stdout]","magenta")+"\n\n")
    sys.stderr.write(colour("Algorithm options:\n","bold"))
    sys.stderr.write("   "+colour("-l ","blue_bold")+" INT          min. length of alignment for the reads (number of nucleotides) "+colour("[75]","magenta")+"\n")
    #sys.stderr.write("   -L FLOAT        min. length of alignment for the reads (percentage of the average read length) [None]\n")
    sys.stderr.write("   "+colour("-v ","blue_bold")+" INT          verbose level: 1=error, 2=warning, 3=message, 4+=debugging "+colour("[3]","magenta")+"\n")
    sys.stderr.write("   "+colour("-y ","blue_bold")+" STR          type of read counts "+colour("[insert.scaled_counts]","magenta")+"\n")
    sys.stderr.write("                    Values: [base.coverage, insert.raw_counts, insert.scaled_counts]\n\n")

# ------------------------------------------------------------------------------
def print_menu_map_lgs():
    sys.stderr.write("\n")
    sys.stderr.write(colour("Usage:","blue_bold")+" motus calc_motu [options]\n\n")
    sys.stderr.write(colour("Input options:\n","bold"))
    sys.stderr.write("   "+colour("-n ","blue_bold")+" STR   sample name\n")
    sys.stderr.write("   "+colour("-i ","blue_bold")+" FILE  provide the mgc abundance table (output of motus calc_mgc)\n")
    sys.stderr.write("   "+colour("-db ","blue_bold")+"DIR   provide a database directory\n\n")
    sys.stderr.write(colour("Output options:\n","bold"))
    sys.stderr.write("   "+colour("-o ","blue_bold")+" FILE  output file name "+colour("[stdout]","magenta")+"\n")
    sys.stderr.write("   "+colour("-e ","blue_bold")+"       profile only reference species (ref_mOTUs)\n")
    sys.stderr.write("   "+colour("-B ","blue_bold")+"       print result in BIOM format\n")
    sys.stderr.write("   "+colour("-C ","blue_bold")+" STR   print result in CAMI format (BioBoxes format 0.9.1)\n")
    sys.stderr.write("             Values: [precision, recall, parenthesis]\n")
    sys.stderr.write("   "+colour("-A ","blue_bold")+"       print all taxonomy levels together\n")
    sys.stderr.write("   "+colour("-c ","blue_bold")+"       print result as counts instead of relative abundances\n")
    sys.stderr.write("   "+colour("-p ","blue_bold")+"       print NCBI id\n")
    sys.stderr.write("   "+colour("-u ","blue_bold")+"       print the full name of the species\n")
    sys.stderr.write("   "+colour("-q ","blue_bold")+"       print the full rank taxonomy\n")
    sys.stderr.write("   "+colour("-k ","blue_bold")+" STR   taxonomic level "+colour("[mOTU]","magenta")+"\n")
    sys.stderr.write("             Values: [kingdom, phylum, class, order, family, genus, mOTU]\n\n")
    sys.stderr.write(colour("Algorithm options:\n","bold"))
    sys.stderr.write("   "+colour("-g ","blue_bold")+"  INT   number of marker genes cutoff: 1=higher recall, 10=higher precision "+colour("[3]","magenta")+"\n")
    sys.stderr.write("   "+colour("-v ","blue_bold")+"  INT   verbose level: 1=error, 2=warning, 3=message, 4+=debugging "+colour("[3]","magenta")+"\n\n")

# ------------------------------------------------------------------------------
def print_menu_append():
    sys.stderr.write("\n")
    sys.stderr.write(colour("Usage:","blue_bold")+" motus merge [options]\n\n")
    sys.stderr.write(colour("Input options:\n","bold"))
    sys.stderr.write("   "+colour("-i ","blue_bold")+"FILE[,FILE] list of mOTU profiles to merge (comma separated)\n")
    sys.stderr.write("   "+colour("-d ","blue_bold")+"DIR         merge all the files in the directory DIR\n")
    sys.stderr.write("   "+colour("-a ","blue_bold")+"STR[,STR]   add pre-computed profiles from different environmental samples\n                  Values: [all, air, bioreactor, bee, cat,\n                  cattle, chicken, dog, fish, freshwater, human,\n                  marine, mouse, pig, sheep, soil, termite, wastewater]\n\n")
    sys.stderr.write(colour("Output options:\n","bold"))
    sys.stderr.write("   "+colour("-o ","blue_bold")+"FILE        output file name "+colour("[stdout]","magenta")+"\n")
    sys.stderr.write("   "+colour("-B ","blue_bold")+"             print result in BIOM format\n\n")
    sys.stderr.write(colour("Algorithm options:\n","bold"))
    sys.stderr.write("   "+colour("-v ","blue_bold")+"INT         verbose level: 1=error, 2=warning, 3=message, 4+=debugging "+colour("[3]","magenta")+"\n\n")
