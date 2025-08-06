import logging
import sys
import pathlib


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