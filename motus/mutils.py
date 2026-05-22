import logging
import sys
import pathlib
import shutil
from motus.mentities import MOTUS_PARAMETERS
from motus.mentities import MOTUS_DB
import Bio.SeqIO.FastaIO as FastaIO
import Bio.SeqIO.QualityIO as QualityIO
import gzip

SAM_ID_FLAG = 'mOTUs4'
MOTUS_VERSION = '4.1.0'
DEFAULT_MOTUS_MGDB_PARENT_LOCATION = pathlib.Path(__file__).resolve().parent
DEFAULT_MOTUS_MGDB_LOCATION = DEFAULT_MOTUS_MGDB_PARENT_LOCATION.joinpath('db_mOTU')
DEFAULT_MOTUS_MGDB_LOCATION_MARKER = DEFAULT_MOTUS_MGDB_LOCATION.joinpath('db_mOTU.downloaded')
DEFAULT_MOTUS_ANNODB_LOCATION = DEFAULT_MOTUS_MGDB_LOCATION.joinpath('mOTUsv4.1.annotation.db')
DEFAULT_MOTUS_ANNODB_LOCATION_MARKER = DEFAULT_MOTUS_MGDB_LOCATION.joinpath('mOTUsv4.1.annotation.db.downloaded')


MOTUS_MGDB_REMOTE_LOCATION_41_toy = 'https://zenodo.org/records/20322003/files/db_mOTU.tar.gz'
MOTUS_MGDB_REMOTE_LOCATION_41_toy_version = '4.1-toy'
MOTUS_MGDB_REMOTE_LOCATION_40 = 'https://zenodo.org/records/17668622/files/db_mOTU.tar.gz'
MOTUS_MGDB_REMOTE_LOCATION_40_version = '4.0'
MOTUS_MGDB_REMOTE_LOCATION_41 = 'https://zenodo.org/records/20322482/files/db_mOTU.tar.gz'
MOTUS_MGDB_REMOTE_LOCATION_41_version = '4.1'

MOTUS_MGDB_REMOTE_LOCATION = MOTUS_MGDB_REMOTE_LOCATION_41

MOTUS_GENOME_REMOTE_PREFIX = 'https://sunagawalab.ethz.ch/share/MOTUS/database/4.0/data/genomes/'


MOTUS_ANNODB_REMOTE_LOCATION_40 = 'https://zenodo.org/records/17669279/files/mOTUsv4.0.annotation.db' 
MOTUS_ANNODB_REMOTE_LOCATION_41 = 'https://zenodo.org/records/20343612/files/mOTUsv4.1.annotation.db' 
MOTUS_ANNODB_REMOTE_LOCATION_40_version = '4.0' 
MOTUS_ANNODB_REMOTE_LOCATION_41_version = '4.1' 
MOTUS_ANNODB_REMOTE_LOCATION = MOTUS_ANNODB_REMOTE_LOCATION_41




def normalize_header(header: str) -> str:
    if header.endswith('/1') or header.endswith('/2'):
        return header[:-2]
    return header


def shutdown(exitcode: int) -> None:
    """
    Securily shutdown the mOTU tool.
    Args:
        exitcode: The exitcode with which mOTU should shutdown.

    Returns:
        None
    """
    logging.info(f'mOTU tool shutting down with exitcode {exitcode}')
    sys.exit(exitcode)


def startup() -> None:
    """A method to group all functions that should be
    executed during startup of the mOTU tool
    """

    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=logging.INFO, datefmt='%Y-%m-%d,%H:%M:%S')

    logging.info(f'mOTU tool starting - {SAM_ID_FLAG}:{MOTUS_VERSION}')
    if shutil.which('bwa') is None:
        logging.error('bwa is not installed or not on PATH. Please install bwa 0.7.18.')
        shutdown(1)
    if shutil.which('vsearch') is None:
        logging.warning('vsearch is not installed or not on PATH. The classify command will not work.')

def cite_text() -> str:
    """Returns the mOTUs4 citations

    Returns:
        str: mOTUs4 citation text
    """    
    tmp = '''
    References:
        Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand
        taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022).
        doi: https://doi.org/10.1186/s40168-022-01410-z

        Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible
        genomic context to taxonomic profiling of microbial communities. Nucleic Acids Research (2025).
        doi: https://doi.org/10.1093/nar/gkae1004
    '''
    return tmp


def check_validity_of_mgc_header(header_line: str) -> None:
    """check of the header line matches the parameters used in the current
    call. Example:
    #tool_version=4.0.2     database_version=4.0    min_alignment_length=110

    what is checked?
    1. tool_version
    2. database_version

    what is not checked
    1. min_alignment_length: This method only makes sense for calc_motus and there
        alignment length is not used. So it is ignored here

    This method will kill the current job if parameters don\'t match

    Args:
        header_line (str): The header line of an MGC file
    """    

    if not header_line or not header_line.startswith('#'):
        logging.error(f'Header line: {header_line} doesn\'t look like a valid header')
        shutdown(1)
    header_line_splits = header_line.strip().split('\t')

    tool = header_line_splits[0].strip().split('#tool_version=')[-1]
    version = header_line_splits[1].strip().split('database_version=')[-1]
    min_aln_length = int(header_line_splits[2].strip().split('min_alignment_length=')[-1])
    if tool != MOTUS_DB.get_tool_version():
        logging.error('MGC file was created with a different version of the tool.')
        logging.error(f'Tool version in MGC file: {tool}')
        logging.error(f'This tool version: {MOTUS_DB.get_tool_version()}')
        shutdown(1)

    if version != MOTUS_DB.get_database_version():
        logging.error('MGC file was created with a different version of the database.')
        logging.error(f'Database version in MGC file: {version}')
        logging.error(f'This database version: {MOTUS_DB.get_database_version()}')
        shutdown(1)
    MOTUS_PARAMETERS.set_minimal_alignment_length(min_aln_length)


def create_mgc_header_line():

    header_line = f'#tool_version={MOTUS_DB.get_tool_version()}\tdatabase_version={MOTUS_DB.get_database_version()}\tmin_alignment_length={MOTUS_PARAMETERS.get_minimal_alignment_length()}'
    return header_line


def is_gzipped(filename: pathlib.Path):
    """
    Check if a file is gzipped by inspecting its magic number.
    """
    with open(filename, 'rb') as f:
        magic = f.read(2)
    return magic == b'\x1f\x8b'


def yield_reads(reads_file: pathlib.Path):
    allowed_file_fq_endings = ['fq.gz', 'fq', 'fastq', 'fastq.gz']
    allowed_file_fa_endings = ['fa', 'fa.gz', 'fasta', 'fasta.gz', 'fna', 'fna.gz']
    is_fq = False
    is_fa = False
    is_gz = False
    if is_gzipped(reads_file):
        is_gz = True

    for allowed_file_fa_ending in allowed_file_fa_endings:
        if str(reads_file).endswith(allowed_file_fa_ending):
            is_fa = True
    for allowed_file_fq_ending in allowed_file_fq_endings:
        if str(reads_file).endswith(allowed_file_fq_ending):
            is_fq = True

    if is_gz:
        of = gzip.open(reads_file, 'rt')
    else:
        of = open(reads_file, 'r')

    if is_fa:
        for (header, sequence) in FastaIO.SimpleFastaParser(of):
            yield (header.strip().split()[0], sequence)
    elif is_fq:
        for header, sequence, qual in QualityIO.FastqGeneralIterator(of):
            yield (header.strip().split()[0], sequence)
    else:
        logging.error(f'Unknown file format: {reads_file}. Expecting a fasta or fastq file, can be gzipped.')
        shutdown(1)
    of.close()
