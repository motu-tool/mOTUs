import logging
import sys
import pathlib
from mentities import MOTUS_PARAMETERS
from mentities import MOTUS_DB

SAM_ID_FLAG = 'mOTUs4'
MOTUS_VERSION = '4.0.2'
DEFAULT_MOTUS_MGDB_PARENT_LOCATION = pathlib.Path(__file__).resolve().parent
DEFAULT_MOTUS_MGDB_LOCATION = DEFAULT_MOTUS_MGDB_PARENT_LOCATION.joinpath('db_mOTU')
DEFAULT_MOTUS_MGDB_LOCATION_MARKER = DEFAULT_MOTUS_MGDB_LOCATION.joinpath('db_mOTU.downloaded')
MOTUS_MGDB_REMOTE_LOCATION = 'https://sunagawalab.ethz.ch/share/MOTUS/database/4.0/data/mOTUS-MGDB/current/db_mOTU.tar.gz'
MOTUS_GENOME_REMOTE_PREFIX = 'https://sunagawalab.ethz.ch/share/MOTUS/database/4.0/data/genomes/'


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
    """
    A method to group all functions that should be
    executed during startup of the mOTU tool.
    Returns:
        None
    """
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=logging.INFO, datefmt='%Y-%m-%d,%H:%M:%S')

    #TODO TEST if bwa is installed and working
    logging.info(f'mOTU tool starting - {SAM_ID_FLAG}:{MOTUS_VERSION}')

def cite_text() -> str:
    tmp = '''
    References:
    
    Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand 
    taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022). 
    doi: https://doi.org/10.1186/s40168-022-01410-z

    Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible 
    genomic context to taxonomic profiling of microbial communities. Nuclic Acids Research (2025). 
    doi: https://doi.org/10.1093/nar/gkae1004
    '''
    return tmp


def check_validity_of_mgc_header(header_line: str) -> None:
    '''check of the header line matches the parameters used in the current
    call. Example:
    #tool_version=4.0.2     database_version=4.0    min_alignment_length=110

    what is checked?
    1. tool_version
    2. database_version

    what is not checked
    1. min_alignment_length: This method only makes sense for calc_motus and there
        alignment length is not used. So it is ignored here

    This method will kill the current job if parameters don\'t match

    Returns:
        None

    '''

    if not header_line or not header_line.startswith('#'):
        logging.error(f'Header line: {header_line} doesn\'t look like a valid header')
        shutdown()
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