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
'''+colour("      profile","cyan")+'''     Perform a taxonomic profiling (map_tax + calc_mgc + calc_motu)
'''+colour("      merge","cyan")+'''       Merge different taxonomic profiles to create a table

'''+colour("      map_tax","cyan")+'''     Map reads to the marker gene database
'''+colour("      calc_mgc","cyan")+'''    Aggregate reads from the same marker gene cluster (mgc)
'''+colour("      calc_motu","cyan")+'''   From a mgc abundance table, produce the mOTU abundance table

'''+colour("-- SNV calling","bold")+'''
'''+colour("      map_snv","cyan")+'''     Map reads to create bam/sam file for snv calling
'''+colour("      snv_call","cyan")+'''    SNV calling using metaSNV

Type motus <command> to print the help for a specific command
        '''
    return str_msg




# ------------------------------------------------------------------------------
def print_menu_profile():
    sys.stderr.write("\n")
    sys.stderr.write(colour("Usage:","blue_bold")+" motus profile [options]\n\n")
    sys.stderr.write(colour("Input options:\n","bold"))
    sys.stderr.write("   "+colour("-f ","blue")+" FILE[,FILE]  input file(s) for reads in forward orientation, fastq formatted\n")
    sys.stderr.write("   "+colour("-r ","blue")+" FILE[,FILE]  input file(s) for reads in reverse orientation, fastq formatted\n")
    sys.stderr.write("   "+colour("-s ","blue")+" FILE[,FILE]  input file(s) for reads without mate, fastq formatted\n")
    sys.stderr.write("   "+colour("-n ","blue")+" STR          sample name "+colour("['unnamed sample']","magenta")+"\n")
    sys.stderr.write("   "+colour("-i ","blue")+" FILE[,FILE]  provide sam or bam input file(s)\n")
    sys.stderr.write("   "+colour("-m ","blue")+" FILE         provide a mgc reads count file\n")
    sys.stderr.write("   "+colour("-db ","blue")+"DIR          provide a database directory\n\n")
    sys.stderr.write(colour("Output options:\n","bold"))
    sys.stderr.write("   "+colour("-o ","blue")+" FILE         output file name "+colour("[stdout]","magenta")+"\n")
    sys.stderr.write("   "+colour("-I ","blue")+" FILE         save the result of bwa in bam format (intermediate step)\n")
    sys.stderr.write("   "+colour("-M ","blue")+" FILE         save the mgc reads count (intermediate step)\n")
    sys.stderr.write("   "+colour("-e ","blue")+"              profile only reference species (ref_mOTUs)\n")
    sys.stderr.write("   "+colour("-c ","blue")+"              print result as counts instead of relative abundances\n")
    sys.stderr.write("   "+colour("-p ","blue")+"              print NCBI id\n")
    sys.stderr.write("   "+colour("-u ","blue")+"              print the full name of the species\n")
    sys.stderr.write("   "+colour("-q ","blue")+"              print the full rank taxonomy\n")
    sys.stderr.write("   "+colour("-B ","blue")+"              print result in BIOM format\n")
    sys.stderr.write("   "+colour("-C ","blue")+" STR          print result in CAMI format (BioBoxes format 0.9.1)\n")
    sys.stderr.write("                    Values: [precision, recall, parenthesis]\n")
    sys.stderr.write("   "+colour("-A ","blue")+"              print all taxonomy levels together\n")
    sys.stderr.write("   "+colour("-k ","blue")+" STR          taxonomic level "+colour("[mOTU]","magenta")+"\n")
    sys.stderr.write("                    Values: [kingdom, phylum, class, order, family, genus, mOTU]\n\n")
    sys.stderr.write(colour("Algorithm options:\n","bold"))
    sys.stderr.write("   "+colour("-g ","blue")+" INT          number of marker genes cutoff: 1=higher recall, 10=higher precision "+colour("[3]","magenta")+"\n")
    sys.stderr.write("   "+colour("-l ","blue")+" INT          min. length of alignment for the reads (number of nucleotides) "+colour("[75]","magenta")+"\n")
    sys.stderr.write("   "+colour("-t ","blue")+" INT          number of threads "+colour("[1]","magenta")+"\n")
    sys.stderr.write("   "+colour("-v ","blue")+" INT          verbose level: 1=error, 2=warning, 3=message, 4+=debugging "+colour("[3]","magenta")+"\n")
    sys.stderr.write("   "+colour("-y ","blue")+" STR          type of read counts "+colour("[insert.scaled_counts]","magenta")+"\n")
    sys.stderr.write("                    Values: [base.coverage, insert.raw_counts, insert.scaled_counts]\n\n")

# ------------------------------------------------------------------------------
def print_menu_map_snv():
    sys.stderr.write("\n")
    sys.stderr.write(colour("Usage:","blue_bold")+" motus map_snv [options]\n\n")
    sys.stderr.write(colour("Input options:\n","bold"))
    sys.stderr.write("   "+colour("-f ","blue")+" FILE[,FILE]  input file(s) for reads in forward orientation, fastq formatted\n")
    sys.stderr.write("   "+colour("-r ","blue")+" FILE[,FILE]  input file(s) for reads in reverse orientation, fastq formatted\n")
    sys.stderr.write("   "+colour("-s ","blue")+" FILE[,FILE]  input file(s) for reads without mate, fastq formatted\n")
    sys.stderr.write("   "+colour("-db ","blue")+"DIR          provide a database directory\n\n")
    sys.stderr.write(colour("Output options:\n","bold"))
    sys.stderr.write("   "+colour("-o ","blue")+" FILE         output bam file name "+colour("[stdout]","magenta")+"\n\n")
    sys.stderr.write(colour("Algorithm options:\n","bold"))
    sys.stderr.write("   "+colour("-l ","blue")+" INT          min. length of alignment for the reads (number of nucleotides) "+colour("[75]","magenta")+"\n")
    #sys.stderr.write("   -CC INT          number of cores for bwa (max one per lane) [1]\n")
    sys.stderr.write("   "+colour("-t ","blue")+" INT          number of threads "+colour("[1]","magenta")+"\n")
    sys.stderr.write("   "+colour("-v ","blue")+" INT          verbose level: 1=error, 2=warning, 3=message, 4+=debugging "+colour("[3]","magenta")+"\n\n")

# ------------------------------------------------------------------------------
def print_menu_snv_call():
    sys.stderr.write("\n")
    sys.stderr.write(colour("Usage:","blue_bold")+" motus snv_call "+colour("-d ","blue")+"Directory "+colour("-o ","blue")+"Directory [options]\n\n")
    sys.stderr.write(colour("Input options:\n","bold"))
    sys.stderr.write("   "+colour("-d ","blue")+" DIR     Call metaSNV on all bam files in the directory\n")
    sys.stderr.write("   "+colour("-fb ","blue")+"FLOAT   Coverage breadth: minimal horizontal genome coverage percentage per sample per species "+colour("[80.0]","magenta")+"\n")
    sys.stderr.write("   "+colour("-fd ","blue")+"FLOAT   Coverage depth: minimal average vertical genome coverage per sample per species "+colour("[5.0]","magenta")+"\n")
    sys.stderr.write("   "+colour("-fm ","blue")+"INT     Minimum number of samples per species "+colour("[2]","magenta")+"\n")
    sys.stderr.write("   "+colour("-fp ","blue")+"FLOAT   FILTERING STEP II: Required proportion of informative samples (coverage non-zero) per position "+colour("[0.90]","magenta")+"\n")
    sys.stderr.write("   "+colour("-fc ","blue")+"FLOAT   FILTERING STEP II: Minimum coverage per position per sample per species "+colour("[5.0]","magenta")+"\n")
    sys.stderr.write("   "+colour("-db ","blue")+"DIR     Provide a database directory\n\n")
    sys.stderr.write(colour("Output options:\n","bold"))
    sys.stderr.write("   "+colour("-o ","blue")+" DIR     Output directory. Will fail if already exists\n")
    sys.stderr.write("   "+colour("-K ","blue")+"         Keep all the directories produced by metaSNV. Default is to remove cov, distances, filtered, snpCaller\n\n")
    sys.stderr.write(colour("Algorithm options:\n","bold"))
    sys.stderr.write("   "+colour("-v ","blue")+"INT      Verbose level: 1=error, 2=warning, 3=message, 4+=debugging "+colour("[3]","magenta")+"\n\n")

# ------------------------------------------------------------------------------
def print_menu_bwa():
    sys.stderr.write("\n")
    sys.stderr.write(colour("Usage:","blue_bold")+" motus map_tax [options]\n\n")
    sys.stderr.write(colour("Input options:\n","bold"))
    sys.stderr.write("   "+colour("-f ","blue")+" FILE[,FILE]  input file(s) for reads in forward orientation, fastq formatted\n")
    sys.stderr.write("   "+colour("-r ","blue")+" FILE[,FILE]  input file(s) for reads in reverse orientation, fastq formatted\n")
    sys.stderr.write("   "+colour("-s ","blue")+" FILE[,FILE]  input file(s) for reads without mate, fastq formatted\n")
    sys.stderr.write("   "+colour("-db ","blue")+"DIR          provide a database directory\n\n")
    sys.stderr.write(colour("Output options:\n","bold"))
    sys.stderr.write("   "+colour("-o ","blue")+" FILE         output file name "+colour("[stdout]","magenta")+"\n")
    sys.stderr.write("   "+colour("-b ","blue")+"              save the result of bwa in bam format\n\n")
    sys.stderr.write(colour("Algorithm options:\n","bold"))
    sys.stderr.write("   "+colour("-l ","blue")+" INT          min. length of alignment for the reads (number of nucleotides) "+colour("[75]","magenta")+"\n")
    #sys.stderr.write("   -CC INT          number of cores for bwa (max one per lane) [1]\n")
    sys.stderr.write("   "+colour("-t ","blue")+" INT          number of threads "+colour("[1]","magenta")+"\n")
    sys.stderr.write("   "+colour("-v ","blue")+" INT          verbose level: 1=error, 2=warning, 3=message, 4+=debugging "+colour("[3]","magenta")+"\n\n")

# ------------------------------------------------------------------------------
def print_menu_map_genes():
    sys.stderr.write("\n")
    sys.stderr.write(colour("Usage:","blue_bold")+" motus calc_mgc [options]\n\n")
    sys.stderr.write(colour("Input options:\n","bold"))
    sys.stderr.write("   "+colour("-n ","blue")+" STR          sample name\n")
    sys.stderr.write("   "+colour("-i ","blue")+" FILE[,FILE]  provide a sam or bam input file (or list of files)\n")
    sys.stderr.write("   "+colour("-db ","blue")+"DIR          provide a database directory\n\n")
    sys.stderr.write(colour("Output options:\n","bold"))
    sys.stderr.write("   "+colour("-o ","blue")+" FILE         output file name "+colour("[stdout]","magenta")+"\n\n")
    sys.stderr.write(colour("Algorithm options:\n","bold"))
    sys.stderr.write("   "+colour("-l ","blue")+" INT          min. length of alignment for the reads (number of nucleotides) "+colour("[75]","magenta")+"\n")
    #sys.stderr.write("   -L FLOAT        min. length of alignment for the reads (percentage of the average read length) [None]\n")
    sys.stderr.write("   "+colour("-v ","blue")+" INT          verbose level: 1=error, 2=warning, 3=message, 4+=debugging "+colour("[3]","magenta")+"\n")
    sys.stderr.write("   "+colour("-y ","blue")+" STR          type of read counts "+colour("[insert.scaled_counts]","magenta")+"\n")
    sys.stderr.write("                    Values: [base.coverage, insert.raw_counts, insert.scaled_counts]\n\n")

# ------------------------------------------------------------------------------
def print_menu_map_lgs():
    sys.stderr.write("\n")
    sys.stderr.write(colour("Usage:","blue_bold")+" motus calc_motu [options]\n\n")
    sys.stderr.write(colour("Input options:\n","bold"))
    sys.stderr.write("   "+colour("-n ","blue")+" STR   sample name\n")
    sys.stderr.write("   "+colour("-i ","blue")+" FILE  provide the mgc abundance table\n")
    sys.stderr.write("   "+colour("-db ","blue")+"DIR   provide a database directory\n\n")
    sys.stderr.write(colour("Output options:\n","bold"))
    sys.stderr.write("   "+colour("-o ","blue")+" FILE  output file name "+colour("[stdout]","magenta")+"\n")
    sys.stderr.write("   "+colour("-e ","blue")+"       profile only reference species (ref_mOTUs)\n")
    sys.stderr.write("   "+colour("-B ","blue")+"       print result in BIOM format\n")
    sys.stderr.write("   "+colour("-C ","blue")+" STR   print result in CAMI format (BioBoxes format 0.9.1)\n")
    sys.stderr.write("             Values: [precision, recall, parenthesis]\n")
    sys.stderr.write("   "+colour("-A ","blue")+"       print all taxonomy levels together\n")
    sys.stderr.write("   "+colour("-c ","blue")+"       print result as counts instead of relative abundances\n")
    sys.stderr.write("   "+colour("-p ","blue")+"       print NCBI id\n")
    sys.stderr.write("   "+colour("-u ","blue")+"       print the full name of the species\n")
    sys.stderr.write("   "+colour("-q ","blue")+"       print the full rank taxonomy\n")
    sys.stderr.write("   "+colour("-k ","blue")+" STR   taxonomic level "+colour("[mOTU]","magenta")+"\n")
    sys.stderr.write("             Values: [kingdom, phylum, class, order, family, genus, mOTU]\n\n")
    sys.stderr.write(colour("Algorithm options:\n","bold"))
    sys.stderr.write("   "+colour("-g ","blue")+" INT   number of marker genes cutoff: 1=higher recall, 10=higher precision "+colour("[3]","magenta")+"\n")
    sys.stderr.write("   "+colour("-v ","blue")+" INT   verbose level: 1=error, 2=warning, 3=message, 4+=debugging "+colour("[3]","magenta")+"\n\n")

# ------------------------------------------------------------------------------
def print_menu_append():
    sys.stderr.write("\n")
    sys.stderr.write(colour("Usage:","blue_bold")+" motus merge [options]\n\n")
    sys.stderr.write(colour("Input options:\n","bold"))
    sys.stderr.write("   "+colour("-i ","blue")+"STR[,STR]  list of files (comma separated)\n")
    sys.stderr.write("   "+colour("-d ","blue")+"DIR        merge all the files in the directory\n")
    sys.stderr.write("   "+colour("-a ","blue")+"STR[,STR]  Add public profiles from different environments. [all, air, bioreactor, bee, cat, \n\t\t cattle, chicken, dog, fish, freshwater, human, \n\t\t marine, mouse, pig, sheep, soil, termite, wastewater]\n\n")
    sys.stderr.write(colour("Output options:\n","bold"))
    sys.stderr.write("   "+colour("-o ","blue")+"FILE       output file name "+colour("[stdout]","magenta")+"\n")
    sys.stderr.write("   "+colour("-B ","blue")+"           print result in BIOM format\n\n")
    sys.stderr.write(colour("Algorithm options:\n","bold"))
    sys.stderr.write("   "+colour("-v ","blue")+"INT        verbose level: 1=error, 2=warning, 3=message, 4+=debugging "+colour("[3]","magenta")+"\n\n")