import pysam
import Bio
import logging
from typing import List, Dict


class MotusDB():
    """
    A class to keep all relevant database information such as:
    - MG - MGC - MOTU
    - Taxonomy per mOTU
    - Version
    """

    database_version: str = None
    mg_2_mgc: Dict[str, str] = {}
    mg_2_mglength: Dict[str, int] = {}
    mgc_2_motu: Dict[str, str] = {}
    motu_2_taxonomy: Dict[str, str] = {}
    index_location: str = None


    def __int__(self, mOTUsdb_folder: str) -> None:
        """
        loads the contents of the mOTUs database
        Following files are expected:
        1. mOTUs.version --> holds the version of the database
        2. mOTUsNR.fasta.gz --> Marker gene sequences in gzipped fasta file
        3. mOTUsNR.fasta.gz.* --> the BWA index
        4. mOTUs.MG.metadata.tsv --> MG MGC MOTU LENGTH
        5. mOTUs.MOTU.metadata.tsv --> MOTU TAX_GTDB TAX_NCBI

        :param mOTUsdb_folder:
        :return: None
        """
        x = 0
        # TODO







def profile(motusdb: MotusDB, forward_files: List[str], reverse_files: List[str], unpaired_files: List[str], motus_file: str, bam_file: str =None, mgc_file: str = None, samplename: str = 'unnamed sample', threads: int = 1, count_mode: str = 'insert.scaled_counts', minlength: int = 45, mg_cutoff: int = 3) -> None:
    """
    TODO summarize the 3 methods (map_tax, calc_mgc, calc_motu)
    :param forward_files: List of forward read files in fasta or fastq format, optionally gzipped. Has to match the reverse files.
    :param reverse_files: List of reverse read files in fasta or fastq format, optionally gzipped. Has to match the forward files.
    :param unpaired_files: List of single/merged read files in fasta or fastq format, optionally gzipped.
    :param motus_file: The output file for the mOTUs profile
    :param bam_file: [Optional] Location of intermediate bam file. In case of None, a temporary file will be created in /tmp/. Default=[None]
    :param mgc_file: [Optional] Location of intermediate mgc file. In case of None, a temporary file will be created in /tmp/. Default=[None]
    :param samplename: [Optional] Name of the sample used in mgc and mOTUs file. Default='unnamed sample'
    :param threads: [Optional] Number of threads used for the alignment of reads against the mOTUs database
    :param count_mode: [Optional] Mode of counting inserts/bases. insert.scaled_counts, insert.raw_counts, base.coverage. Default=[insert.scaled_counts]
    :param minlength: [Optional] minimal length of alignment. Default=[45]
    :param mg_cutoff: [Optional] minimal number of MGCs that require to have abundance>0 for a mOTU to be counted as present. Default=[3]
    :return: None
    """
    return None
def map_tax(motusdb: MotusDB, forward_files: List[str], reverse_files: List[str], unpaired_files: List[str], bam_file: str = None, threads: int = 1, minlength: int = 45) -> None:
    """
    Takes a list of forward/reverse/unpaired read files and aligns them against the mOTUs database using the number of specified threads.
    Alignments will be filtered by 97% identity and the defined minimal alignment length. The resulting alignments will be stored in the
    sorted BAM file which is either specified as a parameter or as a temporary file.

    :param forward_files: List of forward read files in fasta or fastq format, optionally gzipped. Has to match the reverse files.
    :param reverse_files: List of reverse read files in fasta or fastq format, optionally gzipped. Has to match the forward files.
    :param unpaired_files: List of single/merged read files in fasta or fastq format, optionally gzipped.
    :param bam_file: [Optional] Location of intermediate bam file. In case of None, a temporary file will be created in /tmp/. Default=[None]
    :param threads: [Optional] Number of threads used for the alignment of reads against the mOTUs database
    :param minlength: [Optional] minimal length of alignment. Default=[45]
    :return: None
    """
    return None
def calc_mgc(motusdb: MotusDB, bam_file: str, mgc_file: str = None, samplename: str = 'unnamed sample', count_mode: str = 'insert.scaled_counts') -> None:
    """
    Takes the BAM file created in the map_tax method and assigned individual alignments to marker genes and then to marker gene clusters.
    Details for the default mode - insert_scaled:
    1. Read the name sorted alignments in the BAM file and pair by insert (R1/R2/S). Always have only one insert in memory.
    2. Find the best alignment per insert using combined paired alignment score.
    3. Check if this is a unique mapper or a multimapper. A multimapper that maps only against MG from the same MGC count as unique mapper
    4. Distribute the unique mappers to individual MGs
    5. Distribute the multimappers to MGs based on the fractional abundance of the unique mappers in those MGs
    6. Normalise abundance by MG/MGC --> TODO write the exact method
    7. Group abundance by MGC and write to file


    :param bam_file: Location of intermediate bam file.
    :param mgc_file: [Optional] Location of intermediate mgc file. In case of None, a temporary file will be created in /tmp/. Default=[None]
    :param count_mode: [Optional] Mode of counting inserts/bases. insert.scaled_counts, insert.raw_counts, base.coverage. Default=[insert.scaled_counts
    :param samplename: [Optional] Name of the sample used in mgc and mOTUs file. Default='unnamed sample'
    :return:
    """
    return None
def calc_motu(motusdb: MotusDB, motus_file: str, mgc_file: str, samplename: str = 'unnamed sample', mg_cutoff: int = 3) -> None:
    """
    Takes the MGC file produced by calc_mgc and produces a mOTUs profile file.


    :param motus_file: The output file for the mOTUs profile
    :param mgc_file: Location of intermediate mgc file.
    :param samplename: [Optional] Name of the sample used in mgc and mOTUs file. Default='unnamed sample'
    :param mg_cutoff: [Optional] minimal number of MGCs that require to have abundance>0 for a mOTU to be counted as present. Default=[3]
    :return:
    """
    return None
def merge_profiles(merged_motus_file: str, motus_files: List[str]) -> None:
    """
    Takes a list of mOTUs profiles created with the same version of mOTUs and the same parameters
    and merges them into a single profile
    :param merged_motus_file: The output file with the merged mOTUs profiles
    :param motus_files:  The mOTUs files in default profile/calc_motu format to merge profiles.

    :return:
    """
    return None
if __name__ == '__main__':
    x = 0
