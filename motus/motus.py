import os
import statistics

import pysam
import Bio.SeqIO.FastaIO as FastaIO
import Bio.SeqIO.QualityIO as QualityIO
import logging
import pathlib
import csv
import subprocess
import gzip
import collections
import random
import argparse
import sys
from typing import List, Dict, Set, Tuple, Generator




"""
Terminology

markergeneheader = mgh = an instance of an markergene
markergene = mg = one of the 10 mOTUs markergenes
markergenecluster = mgc = a set of mgh that come from the same mOTU and markergene
motu = Top level unit, species level cluster

"""
R1IDENTIFIER = '1'
R2IDENTIFIER = '2'
SIDENTIFIER = 'S'

motusfiles = None
motusdb = None

MOTUS_VERSION = '4.0.0'


Mgc_values = collections.namedtuple("Mgc_values", "insert_raw insert_norm insert_scaled base_raw base_norm")
def check_call(command: str) -> None:
    """
    Simple wrapper to execute check_call and catch exceptions
    :param command:
    :return:
    """

    returncode = 1
    try:
        returncode = subprocess.check_call(command, shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        logging.error('Command {} failed with message:\t{}'.format(e.cmd, e.stderr))
        shutdown(returncode)

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

    #TODO TEST if bwa and samtools are installed and working
    logging.info('mOTU tool starting')


class MotusParameters:

    _forward_files: List[pathlib.Path] = []
    _reverse_files: List[pathlib.Path] = []
    _unpaired_files: List[pathlib.Path] = []
    _read_files_were_checked: bool= False
    _alignment_file: pathlib.Path = None
    _temp_alignment_file: pathlib.Path = None

    _mgc_file: pathlib.Path = None
    _motu_file: pathlib.Path = 'UNNAMED_SAMPLE'
    _samplename: str = None
    _min_alignment_length: int = 0
    _threads: int = 1
    _is_strict_db_mode = True

    _count_mode: str = 'INSERT_SCALED'

    '''
    Count modes:
    insert_raw (not actual mode but used for preproc):
        The sum of all inserts that map against a markergene.
        Multimapper inserts are counted fractional
    insert_norm:
        length normalised insert counts:
        for each markergene do:
        mg(insert_raw)/len(mg) / sum(foreach mg: mg(insert_raw)/len(mg))
    insert_scaled:
        add scaling factor to have values above 1
        mg(insert_scaled) = tot_inserts * mg(insert_norm)
        
    base_raw (not actual mode but used for preproc):
        The sum of all bases that map against a markergene
        Multimapper inserts are counted fractional
    base_norm:
        for each markergene do:
        mg(base_raw)/len(mg) / sum(foreach mg: mg(base_raw)/len(mg))
    base_scaled:
        add scaling factor to have values above 1
        mg(base_scaled) = tot_bases * mg(base_norm)
    '''

    count_mode_insert_raw_mode: str = 'INSERT_RAW'
    count_mode_insert_norm_mode: str = 'INSERT_NORM'
    count_mode_insert_scaled_mode: str = 'INSERT_SCALED'
    count_mode_base_raw_mode: str = 'BASE_RAW'
    count_mode_base_norm_mode: str = 'BASE_NORM'
    #count_mode_base_scaled_mode: str = 'base_scaled'

    _min_mgcs: str = 3
    _report_mode = 'counts'

    def is_strict_db_mode(self):
        return self._is_strict_db_mode

    def set_report_mode_rel_abundance(self):
        self._report_mode = 'relab'

    def set_minimal_number_of_mgcs(self, min_mgcs: int) -> None:
        self._min_mgcs = min_mgcs
    def set_count_mode(self, count_mode: str) -> None:
        self._count_mode = count_mode

    def set_minimal_alignment_length(self, minimal_alignment_length: int):
        if minimal_alignment_length < 30:
            logging.error('Minimal alignment length is below aligner threshold. Pick a larger value. Quitting ...')
            shutdown(1)
        if minimal_alignment_length > 150:
            logging.warning('Minimal alignment length set to above average read length of metagenomic sequencing data.')
        self._min_alignment_length = int(minimal_alignment_length)

    def get_minimal_alignment_length(self) -> int:
        return self._min_alignment_length

    def set_threads(self, threads: int):
        if threads < 1:
            logging.error('Threads have to be at least 1')
            shutdown()
        if threads > os.cpu_count():
            logging.warning('Number of threads exceeds the total number of CPU cores.')
        self._threads = int(threads)

    def get_threads(self) -> int:
        return self._threads

    def set_sample_name(self, samplename: str):
        if len(samplename) == 0:
            logging.error('Sample name cannot be empty. Quitting')
            shutdown(1)
        self._samplename = samplename

    def get_sample_name(self) -> str:
        return self._samplename
    def get_count_mode(self) -> str:
        return self._count_mode

    def get_min_mgcs(self) -> int:
        return self._min_mgcs


    def get_mgc_file(self) -> pathlib.Path:
        self._mgc_file.parent.mkdir(exist_ok=True, parents=True)
        return self._mgc_file

    def set_mgc_file(self, mgc_file: pathlib.Path, required_to_exist=True) -> None:
        self._mgc_file = mgc_file
        if required_to_exist:
            if not mgc_file.exists():
                logging.error(f'MGC file {mgc_file} does not exist. Shutting down ...')
                shutdown(1)

    def get_motu_file(self) -> pathlib.Path:
        self._motu_file.parent.mkdir(exist_ok=True, parents=True)
        return self._motu_file


    def set_motu_file(self, motu_file: pathlib.Path, required_to_exist=True) -> None:
        self._motu_file = motu_file
        if required_to_exist:
            if not motu_file.exists():
                logging.error(f'mOTU file {motu_file} does not exist. Shutting down ...')
                shutdown(1)

    def set_alignment_file(self, alignment_file: pathlib.Path, required_to_exist=True) -> None:
        self._alignment_file = alignment_file
        self._temp_alignment_file = pathlib.Path(str(alignment_file) + '_tmp.bam')
        if required_to_exist:
            if not alignment_file.exists():
                logging.error(f'Alignment file {alignment_file} does not exist. Shutting down ...')
                shutdown(1)

        if not str(alignment_file).endswith('.bam'):
            logging.error(f'Alignment file {alignment_file} is/will be a BAM formatted file. Please set file suffix accordingly. Shutting down ...')
            shutdown(1)


    def get_read_files(self) -> List[Tuple[pathlib.Path, str]]:
        read_files = []
        for (r1_file, r2_file) in zip(self._forward_files, self._reverse_files):#, strict=True):
            read_files.append((r1_file, f'/{R1IDENTIFIER}'))
            read_files.append((r2_file, f'/{R2IDENTIFIER}'))
        for u_file in self._unpaired_files:
            read_files.append((u_file, f'/{SIDENTIFIER}'))
        return read_files

    def get_temporary_alignment_file(self) -> pathlib.Path:
        self._temp_alignment_file.parent.mkdir(exist_ok=True, parents=True)
        return self._temp_alignment_file

    def delete_temporary_alignment_file(self) -> None:
        self._temp_alignment_file.unlink(missing_ok=True)

    def get_alignment_file(self) -> pathlib.Path:
        self._alignment_file.parent.mkdir(exist_ok=True, parents=True)
        return self._alignment_file

    def get_first_1000_reads(self, reads_file: pathlib.Path) -> List[Tuple[str, str]]:
        """ Read the first thousand reads
        and check if the file endings are correct.

        Params:
            reads_file: The file with the short read sequencing data

        Returns:
            A list with the first 1000 reads of the file as tuples of
                header and sequence

        """
        # self._forward_files = forward_files
        # self._reverse_files = reverse_files
        # self._unpaired_files = unpaired_files

        allowed_file_fq_endings = ['fq.gz', 'fq', 'fastq', 'fastq.gz']
        allowed_file_fa_endings = ['fa', 'fa.gz', 'fasta', 'fasta.gz', 'fna', 'fna.gz']
        is_fq = False
        is_fa = False
        is_gz = False
        if str(reads_file).endswith('.gz'):
            is_gz = True
        for allowed_file_fa_ending in allowed_file_fa_endings:
            if str(reads_file).endswith(allowed_file_fa_ending):
                is_fa = True
        for allowed_file_fq_ending in allowed_file_fq_endings:
            if str(reads_file).endswith(allowed_file_fq_ending):
                is_fq = True

        reads = []
        if is_gz:
            of = gzip.open(reads_file, 'rt')
        else:
            of = open(reads_file, 'r')
        if is_fa:
            for (header, sequence) in FastaIO.SimpleFastaParser(of):
                if len(reads) >= 1000:
                    break
                reads.append((header.strip().split()[0], sequence))
        elif is_fq:
            for header, sequence, qual in QualityIO.FastqGeneralIterator(of):
                if len(reads) >= 1000:
                    break
                reads.append((header.strip().split()[0], sequence))
        else:
            logging.error(f'Unknown file format: {reads_file}')
            shutdown(1)
        of.close()
        return reads

    def set_read_files(self, forward_files: List[pathlib.Path],  reverse_files: List[pathlib.Path], unpaired_files: List[pathlib.Path],  check_files: bool=True) -> None:
        """ Define set of read files
        that we should align against the mOTUs
        database. This step can/will also check
        if the files exist and to check if foward
        and reverse read files have the same read
        headers

        Params:
            forward_files: A list of pathlike objects
                which are forward read files. Can be fasta
                or fastq. Can be gzipped or uncompressed
            reverse_files: A list of pathlike objects
                which are reverse read files. Can be fasta
                or fastq. Can be gzipped or uncompressed
            unpaired_files: A list of pathlike objects
                which are unpaired read files. Can be fasta
                or fastq. Can be gzipped or uncompressed
        """


        # check existence
        if check_files:
            files_that_dont_exist = []
            if len(forward_files + reverse_files + unpaired_files) == 0:
                logging.error('No input files defined with -f -r or -s. Quitting ...')
                shutdown(1)
            for f in forward_files + reverse_files + unpaired_files:
                if not f.exists():
                    files_that_dont_exist.append(f)
            if len(files_that_dont_exist) != 0:
                logging.error(f'Some read files dont exist: {files_that_dont_exist}')
                for f in files_that_dont_exist:
                    logging.error(f'\t{f}')
                shutdown(1)
            if len(set(forward_files + reverse_files + unpaired_files)) != len(forward_files + reverse_files + unpaired_files):
                logging.error(f'Duplicated read files. Please submit every file only once. Shutting down ...')
                shutdown(1)
            if len(forward_files) != len(reverse_files):
                logging.error('Unequal number of files submitted with -r and -f. Quitting ...')
                shutdown(1)
            for (r1_file, r2_file) in zip(forward_files, reverse_files): #, strict=True):
                r1_reads = self.get_first_1000_reads(r1_file)
                r2_reads = self.get_first_1000_reads(r2_file)
                r1_header = set([r[0] for r in r1_reads])
                r2_header = set([r[0] for r in r2_reads])
                if len(r1_header.symmetric_difference(r2_header)) != 0:
                    logging.error(f'Headers of reads are not identical. Shutting down ...')
                    logging.error(f'Differing read headers: {r1_header.symmetric_difference(r2_header)}')
                    logging.error(f'Differing read headers file 1: {r1_file}')
                    logging.error(f'Differing read headers file 1: {r2_file}')
                    shutdown(1)

            for u_file in unpaired_files:
                u_reads = self.get_first_1000_reads(u_file)

        self._forward_files = forward_files
        self._reverse_files = reverse_files
        self._unpaired_files = unpaired_files





class MotusDB:
    """
    A class to keep all relevant database information such as:
    - MG - MGC - MOTU
    - Taxonomy per mOTU
    - Version
    """

    database_version: str = None
    database_date: str = None
    mgh_2_mgc: Dict[str, str] = {}
    mgh_2_mglength: Dict[str, int] = {}
    mgc_2_motu: Dict[str, str] = {}
    motus: Set[str] = set()
    blocklist_mg = set()
    mgh_2_mg: Dict[str, str] = {}
    mgc_2_mg: Dict[str, str] = {}
    #motu_2_taxonomy: Dict[str, str] = {}
    index_location: pathlib.Path = None
    _motus_core_mgs = ['COG0012','COG0016','COG0018','COG0172','COG0215','COG0495','COG0525','COG0533','COG0541','COG0552']
    _unassigned_motu_name = None

    def __init__(self, mOTUsdb_folder: pathlib.Path) -> None:
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
        logging.info('Loading database ... ')
        versions_file = mOTUsdb_folder.joinpath('mOTUsv4.0.db').resolve()
        index_files = [mOTUsdb_folder.joinpath(f).resolve() for f in ['mOTUsv4.0.db.fna.gz', 'mOTUsv4.0.db.fna.gz.amb','mOTUsv4.0.db.fna.gz.ann','mOTUsv4.0.db.fna.gz.bwt','mOTUsv4.0.db.fna.gz.pac','mOTUsv4.0.db.fna.gz.sa']]
        mgs_file = mOTUsdb_folder.joinpath('mOTUsv4.0.map.tsv.gz').resolve()
        blocklist_file = mOTUsdb_folder.joinpath('mOTUsv4.0.db.blocklist.gz').resolve()
        with open(versions_file) as handle:
            self.database_version = handle.readline().strip().split()[-1]
            self.database_date = handle.readline().strip().split()[-1]
        self.index_location = index_files[0]
        for index_file in index_files + [mgs_file, blocklist_file]:
            if not index_file.exists():
                logging.error(f'Database file {index_file} is missing. Quitting mOTUs...')
                shutdown(1)
        with gzip.open(mgs_file, 'rt') as handle:
            for entry in  csv.DictReader(handle, delimiter='\t'):
                self.mgh_2_mgc[entry['MG']] = entry['MGC']
                self.mgh_2_mglength[entry['MG']] = int(entry['LENGTH'])
                self.mgc_2_motu[entry['MGC']] = entry['#MOTU']
                self.mgh_2_mg[entry['MG']] = entry['COG']
                self.motus.add(entry['#MOTU'])
                self.mgc_2_mg[entry['MGC']] = entry['COG']
                if 'unassigned' in entry['#MOTU']:
                    self._unassigned_motu_name = entry['#MOTU']
        with gzip.open(blocklist_file, 'rt') as handle:
            for line in handle:
                self.blocklist_mg.add(line.strip())

        logging.info(f'Loading database finished. Version {self.database_version} (version date: {self.database_date}) contains {len(self.motus)} mOTUs, {len(self.mgc_2_motu)} markergeneclusters and {len(self.mgh_2_mglength)} markergenes.')

    def is_mg_blocked(self, mg: str) -> bool:
        if mg in self.blocklist_mg:
            return True
        else:
            return False

    def get_full_version(self):
        return 'TOOL:' + MOTUS_VERSION + '_DB:' + self.database_version
    def get_full_sam_id(self):
        return 'mOTUs4'

    def get_mg_by_mgc(self, mgc):
        return self.mgc_2_mg[mgc]
    def is_unassigned_motu(self, motu):
        if not self._unassigned_motu_name:
            logging.error('The unassigned mOTU was not set. This indicates a corrupted database. Please re-download database. Quitting...')
            shutdown(1)
        if motu == self._unassigned_motu_name:
            return True
        else:
            return False
    def get_unassigned_motu(self):
        return self._unassigned_motu_name
    def get_motu_by_mgc(self, mgc):
        return self.mgc_2_motu[mgc]
    def get_bwa_index(self):
        return self.index_location

    def get_mgc_by_mg(self, mgh) -> str:
        return self.mgh_2_mgc[mgh]
    def get_length_by_mg(self, mgh) -> int:
        return self.mgh_2_mglength[mgh]
    def get_mg_by_mgh(self, mgh):
        return self.mgh_2_mg[mgh]

    def get_core_motus_mgs(self) -> List[str]:
        return self._motus_core_mgs







def map_tax() -> None:
    """
    Takes a list of forward/reverse/unpaired read files and aligns them against the mOTUs database.
    Alignments will be filtered by 97% identity and the defined minimal alignment length.
    The resulting alignments will be stored in the sorted BAM file which is either specified as
    a parameter or as a temporary file.


    Returns:
        None

    """
    logging.info('Starting mOTUs - map_tax routine - Alignment against the mOTUs database ... ')
    min_perc_id: float = 97.0
    threads: int = motusfiles.get_threads()
    minlength: int = motusfiles.get_minimal_alignment_length()

    temp_bam_file = motusfiles.get_temporary_alignment_file()
    temp_bam_file_handle = None



    total_reads = 0
    total_mapped_reads = 0

    for readsfile, orientation in motusfiles.get_read_files():
        total_reads_this_file: int = 0
        total_mapped_reads_this_file: Set[str] = set()
        logging.info(f'Aligning {readsfile}')
        command: str = f'bwa mem -a -t {threads} {motusdb.get_bwa_index()} {readsfile}'
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        in_bam_file_handle = pysam.AlignmentFile(process.stdout, 'rb')
        if not temp_bam_file_handle:

            alignmentfile_header = in_bam_file_handle.header.to_dict()
            pg_header = {}
            pg_header['CL'] = 'motus.py map_tax '
            pg_header['PN'] = 'motus.py'
            motus_version = motusdb.get_full_version()
            pg_header['VN'] = motus_version
            pg_header['ID'] = motusdb.get_full_sam_id()
            alignmentfile_header['PG'].append(pg_header)
            temp_bam_file_handle = pysam.AlignmentFile(temp_bam_file, "wb", header = alignmentfile_header)


        for record in in_bam_file_handle:

            if record.is_unmapped:
                total_reads_this_file += 1
                continue
            else:
                if not record.is_secondary and not record.is_supplementary:
                    total_reads_this_file += 1
                if motusdb.is_mg_blocked(record.reference_name):
                    continue
                alnlength: int = sum(record.get_cigar_stats()[0][0:3])
                if alnlength < minlength:
                    continue
                query_covered_bases: int = sum(record.get_cigar_stats()[0][0:2])
                query_length: int = record.infer_read_length()
                mismatches: int = record.get_tag('NM')
                percid: float = (alnlength - mismatches) / float(alnlength) * 100.0
                percid: float = round(percid, 2)
                if min_perc_id > percid:
                    continue
                qcov: float = query_covered_bases / float(query_length)
                record.set_tag('id', percid, 'f')
                record.set_tag('qc', qcov, 'f')
                record.set_tag('al', alnlength, 'i')

                total_mapped_reads_this_file.add(record.qname)
                record.qname = ''.join([record.qname, orientation])
                temp_bam_file_handle.write(record)
        logging.info(f'Finished alignment. Total reads: {total_reads_this_file}, Total aligned reads {len(total_mapped_reads_this_file)}, {round(len(total_mapped_reads_this_file) * 100.0 / total_reads_this_file, 4)}% aligned.')
        total_mapped_reads += len(total_mapped_reads_this_file)
        total_reads += total_reads_this_file
        in_bam_file_handle.close()
    process.stdout.close()
    return_code: int = process.wait()
    if return_code != 0:
        logging.error(f'BWA command failed with return code {return_code}')
        shutdown(1)
    logging.info(f'Finished all alignments. Total reads: {total_reads}, Total aligned reads {total_mapped_reads}, {round(total_mapped_reads * 100.0 / total_reads, 4)}%')
    temp_bam_file_handle.close()

    logging.info(f'Sorting BAM file')
    pysam.sort('-n', '-m', '1G', '-@', '1', '-o', str(motusfiles.get_alignment_file()),  str(motusfiles.get_temporary_alignment_file()))
    logging.info(f'Finished sorting BAM file')
    motusfiles.delete_temporary_alignment_file()
    logging.info('Finished mOTUs - map_tax routine - Alignment against the mOTUs database ...')
    return None

def _get_orientation_of_aligned_segment_by_name(alignment: pysam.AlignedSegment) -> Tuple[str, str]:
        splits = alignment.query_name.rsplit('/', 1)
        if len(splits) == 2:
            return splits[0], splits[1]
        else:
            return alignment.query_name, SIDENTIFIER




class BestAlignment:
    """
    Store information of each best alignment
    of an insert
    """
    _mg_2_blocks = None

    def __init__(self):
        self._mg_2_blocks = {}

    def append(self, mg, blocks):
        self._mg_2_blocks[mg] = blocks

    def isMultimapper(self):
        if len(self._mg_2_blocks) == 1:
            return False
        else:
            return True

    def get_mg_and_blocks(self):
        if self.isMultimapper():
            logging.error('This method doesnt work for multi mappers.')
            shutdown(1)
        for mg, blocks in self._mg_2_blocks.items():
            return mg, blocks

    def get_mgs_and_blocks(self):
        if not self.isMultimapper():
            logging.error('This method doesnt work for unique mappers.')
            shutdown(1)
        return self._mg_2_blocks


class InsertCounter:
    _unique_mappers = None
    _multi_mappers = None


    _mg_2_edge_corrected_raw_uniquemapper_insert_counts = {}
    _mg_2_edge_corrected_raw_uniquemapper_base_counts = {}
    _mg_2_edge_corrected_raw_multimapper_insert_counts = {}
    _mg_2_edge_corrected_raw_multimapper_base_counts = {}

    _mg_2_edge_corrected_raw_insert_counts = {}
    _mg_2_edge_corrected_raw_base_counts = {}
    _mg_2_edge_corrected_norm_insert_counts = {}
    _mg_2_edge_corrected_scaled_insert_counts = {}
    _mg_2_edge_corrected_norm_base_counts = {}
    _mg_2_edge_corrected_scaled_base_counts = {}

    def get_mg_insert_raw(self):
        return self._mg_2_edge_corrected_raw_insert_counts
    def get_mg_base_raw(self):
        return self._mg_2_edge_corrected_raw_base_counts
    def get_mg_insert_norm(self):
        return self._mg_2_edge_corrected_norm_insert_counts
    def get_mg_base_norm(self):
        return self._mg_2_edge_corrected_norm_base_counts
    def get_mg_insert_scaled(self):
        return self._mg_2_edge_corrected_scaled_insert_counts

    def __init__(self):
        self._unique_mappers = []
        self._multi_mappers = []

    def appendmapper(self, insert_name:str, bestAlignment: BestAlignment) -> None:
        if bestAlignment.isMultimapper():
            self._multi_mappers.append((insert_name, bestAlignment))
        else:
            self._unique_mappers.append((insert_name, bestAlignment))

    def get_unique_mapper_count(self):
        return len(self._unique_mappers)
    def get_multi_mapper_count(self):
        return len(self._multi_mappers)

    def correct_multi_mapper_edges(self, min_alignment_length: int):
        mg_2_alignments = collections.defaultdict(list)

        for insert_name, bestAlignment in self._multi_mappers:
            mg_2_blocks = bestAlignment.get_mgs_and_blocks()
            tot_weight = sum([self._mg_2_edge_corrected_raw_uniquemapper_insert_counts.get(mg, 0.0) for mg in mg_2_blocks.keys()])
            if tot_weight < 1.0:
                tot_weight = 1.0

            mg_2_weight = {mg: self._mg_2_edge_corrected_raw_uniquemapper_insert_counts.get(mg, 0.0) / tot_weight for mg in mg_2_blocks.keys()}
            for mg, alignment_blocks in mg_2_blocks.items():
                if mg_2_weight[mg] != 0.0:
                    mg_2_alignments[mg].append((alignment_blocks, mg_2_weight[mg]))

        mg_2_edge_corrected_insert_counts, mg_2_edge_corrected_base_counts = self._correct_edges(mg_2_alignments,min_alignment_length)
        self._mg_2_edge_corrected_raw_multimapper_insert_counts = mg_2_edge_corrected_insert_counts
        self._mg_2_edge_corrected_raw_multimapper_base_counts = mg_2_edge_corrected_base_counts



    def _correct_edges(self, mg_2_alignments, min_alignment_length):

        mg_2_trunc_insert_counts = collections.Counter()
        mg_2_untrunc_insert_counts = collections.Counter()
        mg_2_trunc_base_counts = collections.Counter()
        mg_2_untrunc_base_counts = collections.Counter()

        # alignments -> All alignments against a mg. One alignment represents multiple alignment blocks
        for mg, alignments in mg_2_alignments.items():
            first_allowed_base = min_alignment_length + 1
            mg_len = motusdb.get_length_by_mg(mg)
            last_allowed_base = mg_len - min_alignment_length - 1

            alignments_trunc = []

            for (alignment_blocks, weight) in alignments:

                aligned_bases_untrunc = 0
                aligned_bases_trunc = 0
                for (aln_start, aln_end) in alignment_blocks:
                    if aln_start > aln_end:
                        aln_start_tmp = aln_end
                        aln_end = aln_start
                        aln_start = aln_start_tmp
                    aligned_bases_untrunc += aln_end - aln_start
                    if aln_end < first_allowed_base:
                        continue
                    if aln_start > last_allowed_base:
                        continue
                    if aln_start < first_allowed_base:
                        aln_start = first_allowed_base
                    if aln_end > last_allowed_base:
                        aln_end = last_allowed_base

                    aligned_bases_trunc += aln_end - aln_start
                if aligned_bases_trunc != 0:
                    alignments_trunc.append(aligned_bases_trunc / aligned_bases_untrunc)
                mg_2_untrunc_base_counts[mg] += aligned_bases_untrunc * weight
                mg_2_trunc_base_counts[mg] += aligned_bases_trunc * weight

            mg_2_untrunc_insert_counts[mg] = len(alignments) * weight
            mg_2_trunc_insert_counts[mg] = len(alignments_trunc) * weight

        mg_2_edge_corrected_insert_counts = collections.Counter()
        mg_2_edge_corrected_base_counts = collections.Counter()

        for mg, trunc_insert_count in mg_2_trunc_insert_counts.items():
            mg_len = motusdb.get_length_by_mg(mg)
            mg_trunc_len = mg_len - 2 * min_alignment_length
            edge_corrected_insert_count = mg_len * (trunc_insert_count / mg_trunc_len)
            mg_2_edge_corrected_insert_counts[mg] = edge_corrected_insert_count

        for mg, trunc_base_count in mg_2_trunc_base_counts.items():
            mg_len = motusdb.get_length_by_mg(mg)
            mg_trunc_len = mg_len - 2 * min_alignment_length
            edge_corrected_base_count = mg_len * (trunc_base_count / mg_trunc_len)
            mg_2_edge_corrected_base_counts[mg] = edge_corrected_base_count
        return mg_2_edge_corrected_insert_counts, mg_2_edge_corrected_base_counts
    def correct_uniq_mapper_edges(self, min_alignment_length: int):
        """
        ======================================================================
        Correct the alignment abundances which are biased due to partial
        alignments at the left and right parts of the genes.

        Assume that a read aligns against the left side of a gene but overlaps
        only 5 bases (usually more, just as an example)

        read = =====================
        gene =              ==============================================

        This alignment will not be counted by default for two reasons:
        1. The alignment length is too short which means that motus will not
            report it. In mOTUs the default is 75
        2. The aligner doesn't report it as it is too short for the aligner
            to confidently report it. In BWA this is 30 bases

        This means that abundances at the edges of genes are incorrectly
        reported, short genes suffer proportinally more of this error.


        Edge correction is using abundance of the parts of the gene
        that is correctly reported and extrapolates it to the full length of
        the gene. E.g here a gene has length 30, and anything below 5 is
        reported wrongly:
                      ABUNDANCE
        FALSE         CORRECT      FALSE
        =====|====================|=====

             =====================
             Use this abundance
             Extrapolate to full length
        abundance(full_gene) = len(full_gene) * abundance(trunc_gene) / len(trunc_gene)

        """

        mg_2_alignments = collections.defaultdict(list)
        for insert_name, bestAlignment in self._unique_mappers:
            mg, alignment_blocks = bestAlignment.get_mg_and_blocks()
            mg_2_alignments[mg].append((alignment_blocks, 1.0))
        mg_2_edge_corrected_insert_counts, mg_2_edge_corrected_base_counts = self._correct_edges(mg_2_alignments, 30)

        self._mg_2_edge_corrected_raw_uniquemapper_insert_counts = mg_2_edge_corrected_insert_counts
        self._mg_2_edge_corrected_raw_uniquemapper_base_counts = mg_2_edge_corrected_base_counts


    def combined_raw_counts(self):
        mg_2_edge_corrected_raw_insert_counts = collections.Counter()
        for mg, count in self._mg_2_edge_corrected_raw_uniquemapper_insert_counts.items():
            #print(mg, count, self._mg_2_edge_corrected_raw_multimapper_insert_counts.get(mg, 0.0))
            mg_2_edge_corrected_raw_insert_counts[mg] = count
        for mg, count in self._mg_2_edge_corrected_raw_multimapper_insert_counts.items():
            mg_2_edge_corrected_raw_insert_counts[mg] += count
        mg_2_edge_corrected_raw_base_counts = collections.Counter()
        for mg, count in self._mg_2_edge_corrected_raw_uniquemapper_base_counts.items():
            mg_2_edge_corrected_raw_base_counts[mg] = count
        for mg, count in self._mg_2_edge_corrected_raw_multimapper_base_counts.items():
            mg_2_edge_corrected_raw_base_counts[mg] += count

        self._mg_2_edge_corrected_raw_insert_counts = mg_2_edge_corrected_raw_insert_counts
        self._mg_2_edge_corrected_raw_base_counts = mg_2_edge_corrected_raw_base_counts

    def _norm_and_scale_counts2(self, mg_2_raw_counts):
        tot_cnt = float(sum(mg_2_raw_counts.values()))
        denominator: float = sum([float(mgh_2_count[1]) / float(motusdb.get_length_by_mg(mgh_2_count[0])) for mgh_2_count in mg_2_raw_counts.items()])
        scaled_mg_2_counts = {}
        norm_mg_2_counts = {}
        for mg, count in mg_2_raw_counts.items():
            numerator = float(count) / float(motusdb.get_length_by_mg(mg))
            norm_count = numerator / denominator
            scaled_count = norm_count * tot_cnt
            scaled_mg_2_counts[mg] = scaled_count
            norm_mg_2_counts[mg] = norm_count

        return norm_mg_2_counts, scaled_mg_2_counts




    def norm_and_scale_counts(self):
        mg_2_edge_corrected_norm_insert_counts, mg_2_edge_corrected_scaled_insert_counts = self._norm_and_scale_counts2(self._mg_2_edge_corrected_raw_insert_counts)
        mg_2_edge_corrected_norm_base_counts, mg_2_edge_corrected_scaled_base_counts = self._norm_and_scale_counts2(self._mg_2_edge_corrected_raw_base_counts)
        self._mg_2_edge_corrected_norm_insert_counts = mg_2_edge_corrected_norm_insert_counts
        self._mg_2_edge_corrected_scaled_insert_counts = mg_2_edge_corrected_scaled_insert_counts
        self._mg_2_edge_corrected_norm_base_counts = mg_2_edge_corrected_norm_base_counts
        self._mg_2_edge_corrected_scaled_base_counts = mg_2_edge_corrected_scaled_base_counts


    def _get_alignment_blocks(self, alignments: List[pysam.AlignedSegment]) -> List[Tuple[int, int]]:
        '''
        Get the aligned positions from the alignments. Depending on the
        number of alignments (1 or 2) and the number of indels this number
        can be between 1 and n where n is the readlength. However, in most
        cases the result will be one Tuple for a singleton insert and
        two Tuples for a paired end insert

        Params:
            alignments: A list of pysam AlignedSegments

        Returns:
            A list of Tuples with aligned positions
        '''

        blocks = []
        for alignment in alignments:
            for block in alignment.get_blocks():
                blocks.append(block)
        blocks.sort(key=lambda a: a[0])
        return blocks


    def _filter_best_alignment(self, current_insert: Dict[str, List[pysam.AlignedSegment]],
                               random_mgc_resolver=True) -> BestAlignment:
        """Take a list of alignments against the mOTUs
        genes and find the mOTUs gene with the highest score.
        Will also report aligned bases

        Params:
            current_insert: All alignments associated to one insert.
                Can be paired end or singleton. Can have multiple
                alignments against different genes
            random_mgc_resolver: pick randomly if an insert aligns
                against multiple mgs from the same mgc
        Returns:
            A list of at least one mg this insert is assigned to
            together with the alignment positions
        """

        markergeneheader_2_orientation_2_alignments = collections.defaultdict(lambda: collections.defaultdict(list))
        for orientation, alignments in current_insert.items():
            for alignment in alignments:
                reference_name = alignment.reference_name
                markergeneheader_2_orientation_2_alignments[reference_name][orientation].append(alignment)
        markergeneheader_2_best_alignments = collections.defaultdict(list)
        for markergeneheader, orientation_2_alignments in markergeneheader_2_orientation_2_alignments.items():
            for orientation, alignments in orientation_2_alignments.items():
                best_alignment = alignments[0]
                if len(alignments) != 1:
                    max_score = max([alignment.get_tag('AS') for alignment in alignments])
                    best_alignment = [alignment for alignment in alignments if alignment.get_tag('AS') == max_score][0]
                markergeneheader_2_best_alignments[markergeneheader].append(best_alignment)
        all_alignment_scores = []
        for markergeneheader, best_alignments in markergeneheader_2_best_alignments.items():
            all_alignment_scores.append(sum([alignment.get_tag('AS') for alignment in best_alignments]))
        best_alignment_score = max(all_alignment_scores)
        markergeneheader_2_best_alignments2 = {}
        for markergeneheader, best_alignments in markergeneheader_2_best_alignments.items():
            if sum([alignment.get_tag('AS') for alignment in best_alignments]) == best_alignment_score:
                markergeneheader_2_best_alignments2[markergeneheader] = best_alignments
        best_mgs = BestAlignment()
        if len(markergeneheader_2_best_alignments2) > 1 and random_mgc_resolver:
            mgc_2_mg = collections.defaultdict(list)
            for mg, alns in markergeneheader_2_best_alignments2.items():
                mgc = motusdb.get_mgc_by_mg(mg)
                mgc_2_mg[mgc].append(mg)
            for mgc, mgs in mgc_2_mg.items():
                if len(mgs) == 1:
                    best_mgs.append(mgs[0], self._get_alignment_blocks(markergeneheader_2_best_alignments2[mgs[0]]))
                else:
                    picked_mg = random.choice(mgs)
                    best_mgs.append(picked_mg,
                                    self._get_alignment_blocks(markergeneheader_2_best_alignments2[picked_mg]))
        else:
            for mg in markergeneheader_2_best_alignments2.keys():
                best_mgs.append(mg, self._get_alignment_blocks(markergeneheader_2_best_alignments2[mg]))

        return best_mgs

    def count(self, bam_insert_iterator) -> None:
        for insert_name, alignments in bam_insert_iterator:
            self.appendmapper(insert_name, self._filter_best_alignment(alignments))

        logging.info('Finished reading alignment file ...')
        logging.info(f'Read {self.get_unique_mapper_count() + self.get_multi_mapper_count()} aligned inserts of which {round(100.0 * self.get_multi_mapper_count() / (self.get_unique_mapper_count() + self.get_multi_mapper_count()),2)}% are multimappers')
        self.correct_uniq_mapper_edges(motusfiles.get_minimal_alignment_length())
        self.correct_multi_mapper_edges(motusfiles.get_minimal_alignment_length())
        self.combined_raw_counts()
        self.norm_and_scale_counts()










class MGCCounter:
    """
    A class which takes care of
    reading, parsing and interpreting
    the inserts mapped against the mOTUs
    database.
    """


    def aggregate_mgc(self, mgh_2_scaled_counts, mgh_2_unscaled_counts):
        mgc_2_count = {}
        for mgh, count in mgh_2_scaled_counts.items():
            mgc = motusdb.get_mgc_by_mg(mgh)
            [scaled, unscaled] = mgc_2_count.get(mgc, [0.0, 0.0])
            scaled = scaled + count
            mgc_2_count[mgc] = [scaled, unscaled]

        for mgh, count in mgh_2_unscaled_counts.items():
            mgc = motusdb.get_mgc_by_mg(mgh)
            [scaled, unscaled] = mgc_2_count.get(mgc, [0.0, 0.0])
            unscaled = unscaled + count
            mgc_2_count[mgc] = [scaled, unscaled]
        return mgc_2_count


    def count(self) -> Dict[str, Mgc_values]:
        '''
        Entry Level method for this class
        Read the BAM file and counts abundances
        using different modes (insert_raw, insert_scaled,...)
        '''
        insertcounter = InsertCounter()
        logging.info('Reading alignment file ...')
        insertcounter.count(self._bam_insert_iterator())
        # now aggregate by MGC

        '''
        How to aggregate
        1. for each counting method (insert, base, norm, scaled)
            for each mg
                find mgc
                sum up value for mgc
        2. report
            for each counting method
                for each mgc
                one line with each counting method
        '''

        mgc_insert_raw = collections.defaultdict(lambda: 0.0)
        mgc_insert_norm = collections.defaultdict(lambda: 0.0)
        mgc_insert_scaled = collections.defaultdict(lambda: 0.0)
        mgc_base_raw = collections.defaultdict(lambda: 0.0)
        mgc_base_norm = collections.defaultdict(lambda: 0.0)
        all_mgcs = set()
        for (mg_data, mgc_data) in zip([insertcounter.get_mg_insert_raw(), insertcounter.get_mg_insert_norm(), insertcounter.get_mg_insert_scaled(), insertcounter.get_mg_base_raw(), insertcounter.get_mg_base_norm()], [mgc_insert_raw, mgc_insert_norm, mgc_insert_scaled, mgc_base_raw, mgc_base_norm]):
            for mg, abundance in mg_data.items():
                mgc = motusdb.get_mgc_by_mg(mg)
                mgc_data[mgc] = mgc_data[mgc] + abundance
                all_mgcs.add(mgc)

        mgc_2_all_counts = {}
        for mgc in all_mgcs:
            #Mgc_values = collections.namedtuple("Mgc_values", "insert_raw insert_norm insert_scaled base_raw base_norm")
            insert_raw = mgc_insert_raw[mgc]
            insert_norm = mgc_insert_norm[mgc]
            insert_scaled = mgc_insert_scaled[mgc]
            base_raw = mgc_base_raw[mgc]
            base_norm = mgc_base_norm[mgc]
            mgc_vals = Mgc_values(insert_raw=insert_raw, insert_norm=insert_norm, insert_scaled=insert_scaled, base_raw=base_raw, base_norm=base_norm)
            mgc_2_all_counts[mgc] = mgc_vals

        return mgc_2_all_counts













    def _bam_insert_iterator(self) -> Generator[Tuple[str, Dict[str, List[pysam.AlignedSegment]]], None, None]:
        """Reads through a sorted BAM file and
        finds the best alignment(s) per insert
        """

        alignments = pysam.AlignmentFile(motusfiles.get_alignment_file(), 'r')
        motus_version = [entry for entry in alignments.header.to_dict()['PG'] if entry['ID'] == motusdb.get_full_sam_id()]
        header_valid = False
        if len(motus_version) > 0:
            if motus_version[0]['VN'] == motusdb.get_full_version():
                header_valid = True
        if not header_valid:
            if motusfiles.is_strict_db_mode():
                logging.error('mOTUs tool/database have changed and bam file is invalid. Please profile with updated database. Quitting ...')
                shutdown(1)
            else:
                logging.warning('mOTUs tool/database have changed and bam file is invalid. Lenient mode enabled, will continue but results might be broken ...')

        try:
            alignment: pysam.AlignedSegment = next(alignments)
        except StopIteration:
            alignments.close()
            logging.info(f'The alignmentfile {motusfiles.get_alignment_file()} has no valid alignments. Quitting ...')
            shutdown(1)
            return
        current_name, orientation = _get_orientation_of_aligned_segment_by_name(alignment)
        current_insert = collections.defaultdict(list)
        orientations = set()
        orientations.add(orientation)
        current_insert[orientation].append(alignment)
        minlength: int = motusfiles.get_minimal_alignment_length()
        readname = None

        for alignment in alignments:
            if motusdb.is_mg_blocked(alignment.reference_name):
                continue
            alnlength: int = sum(alignment.get_cigar_stats()[0][0:3])
            if alnlength < minlength:
                continue
            (readname, orientation) = _get_orientation_of_aligned_segment_by_name(alignment)
            if readname == current_name:
                current_insert[orientation].append(alignment)
                orientations.add(orientation)
            else:
                if len(orientations) == 3 or (len(orientations) > 1 and SIDENTIFIER in orientations):
                    raise Exception('An alignment cannot be Paired End and Single End at the same time. Problematic insert: {}'.format(readname))

                yield current_name, current_insert
                current_name = readname
                current_insert = collections.defaultdict(list)
                orientations = set()
                current_insert[orientation].append(alignment)
                orientations.add(orientation)

        if len(orientations) == 3 or (len(orientations) > 1 and SIDENTIFIER in orientations):
            raise Exception('An alignment cannot be Paired End and Single End at the same time. Problematic insert: {}'.format(readname))
        yield current_name, current_insert
        alignments.close()




def calc_mgc() -> None:
    """
    Takes the BAM file created in the map_tax method and assigns individual alignments to marker genes and next to marker gene clusters.
    1. Read the name sorted alignments in the BAM file and pair by insert (1/2/S). Always have only one insert in memory.
    2. Find the best alignment per insert using additive paired alignment score.
    3. Check if this is a unique mapper or a multimapper. A multimapper that maps only against MG from the same MGC counts as unique mapper (and will be randomly assigned to a mg)
    4. Distribute the unique mappers to individual MGs
    5. Distribute the multimappers to MGs based on the fractional abundance of the unique mappers in those MGs
    6. Apply edge correction
    7. calculate using different count modes
    8. Group abundance by MGC and write to file


    Returns:
        None

    """
    logging.info('Starting mOTUs - calc_mgc routine - Calculating abundances per MGC ... ')
    mgc_counter = MGCCounter()
    mgc_2_counts = mgc_counter.count()

    with open(motusfiles.get_mgc_file(), 'w') as handle:
        header_line = motusdb.get_full_version()
        handle.write(f'#{header_line}\n')
        handle.write('MGC\tINSERT_RAW\tINSERT_NORM\tINSERT_SCALED\tBASE_RAW\tBASE_NORM\n')
        for mgc in sorted(mgc_2_counts.keys()):
            counts =  mgc_2_counts[mgc]
            handle.write(f'{mgc}\t{round(counts.insert_raw, 4):.4f}\t{round(counts.insert_norm, 10):.10f}\t{round(counts.insert_scaled, 4):.4f}\t{round(counts.base_raw, 4):.4f}\t{round(counts.base_norm, 10):.10f}\n')


    logging.info('Finished mOTUs - calc_mgc routine - Calculating abundances per MGC ... ')


    return None


def calc_motu() -> None:
    """
    Takes the MGC file produced by calc_mgc and produces a mOTUs profile file.


    Returns
        None
    """
    mgc_file = motusfiles.get_mgc_file()
    has_header = False
    with open(mgc_file) as handle:
        first_line = handle.readline().strip()
        if first_line.startswith('#'):
            has_header = True
        motus_version = first_line.replace('#', '')
        if motus_version == motusdb.get_full_version():
            header_valid = True
        if not header_valid:
            if motusfiles.is_strict_db_mode():
                logging.error('mOTUs tool/database have changed and bam file is invalid. Please profile with updated database. Quitting ...')
                shutdown(1)
            else:
                logging.warning('mOTUs tool/database have changed and bam file is invalid. Lenient mode enabled, will continue but results might be broken ...')



    mgc_2_count = {}
    count_mode = motusfiles.get_count_mode()
    with open(mgc_file) as handle:
        if has_header:
            handle.readline()
        for entry in csv.DictReader(handle, delimiter='\t'):
            mgc_2_count[entry['MGC']] = float(entry[count_mode])


    motu_2_mgccounts = collections.defaultdict(lambda: collections.defaultdict(lambda: 0.0))


    for mgc, count in mgc_2_count.items():
        motu = motusdb.get_motu_by_mgc(mgc)
        mg = motusdb.get_mg_by_mgc(mgc)
        motu_2_mgccounts[motu][mg] += count

    with open(motusfiles.get_motu_file(), 'w') as handle:
        handle.write(f'MOTU\t{motusfiles.get_sample_name()}\n')
        for motu in sorted(list(motu_2_mgccounts.keys())):
            counts = list(motu_2_mgccounts[motu].values())

            median_count = statistics.median(counts)



            count = '{number:.{digits}f}'.format(number=median_count, digits=8)
            if len(counts) >= motusfiles.get_min_mgcs() or motusdb.is_unassigned_motu(motu):
                handle.write(f'{motu}\t{count}\n')







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




class CapitalisedHelpFormatter(argparse.HelpFormatter):
    def add_usage(self, usage, actions, groups, prefix=None):
        if prefix is None:
            prefix = ''
        return super(CapitalisedHelpFormatter, self).add_usage(usage, actions, groups, prefix)


def parse_map_tax():
    parser = argparse.ArgumentParser(usage = '''Program: motus - a tool for marker gene-based OTU (mOTU) profiling
Version: 4.0.0
Reference: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand 
taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022). 
doi: https://doi.org/10.1186/s40168-022-01410-z
    
motus map_tax [options]

Input options:
   -f   FILE[ FILE]  input file(s) for reads in forward orientation, fastq(.gz)-formatted
   -r   FILE[ FILE]  input file(s) for reads in reverse orientation, fastq(.gz)-formatted
   -s   FILE[ FILE]  input file(s) for unpaired reads, fastq(.gz)-formatted


Output options:
   -o   FILE         output file name

Algorithm options:
   -l   INT          min length of the alignment (bp) [75]
   -t   INT          number of threads [1]
   -v   INT          verbosity level: 1=error, 2=warning, 3=message, 4+=debugging [3]
      ''', formatter_class=CapitalisedHelpFormatter,add_help=False)

    # Input options
    parser.add_argument("-f", nargs="+",default=[])  # input files(s) for reads in forward orientation, fastq(.gz)-formatted
    parser.add_argument("-r", nargs="+",default=[])  # input files(s) for reads in reverse orientation, fastq(.gz)-formatted
    parser.add_argument("-s", nargs="+", default=[])  # input files(s) for unpaired reads, fastq(.gz)-formatted
    #parser.add_argument("-db")  # provide a different database directory

    # Output options
    parser.add_argument("-o", required=True)  # output file name
    #parser.add_argument("-b", action="store_true")  # save the result in BAM format

    # ALgorithm options
    parser.add_argument("-l", type=int, default=75)  # min length of the alignment (bp) [75]
    parser.add_argument("-t", type=int, default=1)  # number of threads
    parser.add_argument("-v", type=int, default=1)  # verbodisty level:

    args = parser.parse_args(sys.argv[2:])

    # print usage and exit if no arguments are passed
    if sys.argv[2:] == []:
        parser.print_usage()
        shutdown(1)

    # convert string arguments into Pathlib objects
    forward_files = [pathlib.Path(el) for el in args.f]
    reverse_files = [pathlib.Path(el) for el in args.r]
    unpaired_files = [pathlib.Path(el) for el in args.s]
    alignment_file = pathlib.Path(args.o)
    startup()
    threads = args.t
    min_alignment_length = args.l


    global motusdb
    motusdb = MotusDB(db_folder)
    global motusfiles
    motusfiles = MotusParameters()

    motusfiles.set_read_files(forward_files, reverse_files, unpaired_files, check_files=True)
    motusfiles.set_alignment_file(alignment_file, required_to_exist=False)
    motusfiles.set_minimal_alignment_length(min_alignment_length)
    motusfiles.set_threads(threads)
    map_tax()




db_folder = pathlib.Path('/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/DEV/MOTUSv4/mOTUs4-dev/db_mOTU/')

def parse_profile():
    parser = argparse.ArgumentParser(usage = '''Program: motus - a tool for marker gene-based OTU (mOTU) profiling
Version: 4.0.0
Reference: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand 
taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022). 
doi: https://doi.org/10.1186/s40168-022-01410-z
    
motus profile [options]

Input options:
   -f  FILE[ FILE]  input file(s) for reads in forward orientation, fastq(.gz)-formatted
   -r  FILE[ FILE]  input file(s) for reads in reverse orientation, fastq(.gz)-formatted
   -s  FILE[ FILE]  input file(s) for unpaired reads, fastq(.gz)-formatted
   -n  STR          sample name ['unnamed sample']

Output options:
   -o  FILE         output file name [required]
   -c               print result as counts instead of relative abundances

Algorithm options:
   -g  INT          number of marker genes cutoff: 1=higher recall, 6=higher precision [3]
   -l  INT          min length of the alignment (bp) [75]
   -t  INT          number of threads [1]
   -v  INT          verbosity level: 1=error, 2=warning, 3=message, 4+=debugging [1]
   -y  STR          type of read counts [INSERT_SCALED]
                    Values: [INSERT_RAW, INSERT_NORM, INSERT_SCALED, BASE_RAW, BASE_NORM]
]''', formatter_class=CapitalisedHelpFormatter,add_help=False)

    # Input options
    parser.add_argument("-f", nargs="+", default=[])  # input file(s) for reads in forward direction
    parser.add_argument("-r", nargs="+", default=[])  # input file(s) for reads in reverse direction
    parser.add_argument("-s", nargs="+", default=[])  # input file(s) for unpaired reads
    parser.add_argument("-n", type=str, default='unnamed sample')  # sample name
    #parser.add_argument("-i", nargs="+")  # provide SAM or BAM input files (generated by motus map_tax)
    #parser.add_argument("-m")  # provide mgc reads count file (generated by motus calc_mgc)
    #parser.add_argument("-db")  # provide a different DB directory

    # Output options
    parser.add_argument("-o", required=True)  # output file name
    #parser.add_argument("-e", action="store_true")  # only species with reference genomes (ref-mOTUs)
    #parser.add_argument("-u", action="store_true")  # print the full name of the species
    parser.add_argument("-c", action="store_true")  # print result as counts instead of relative abundances
    #parser.add_argument("-p", action="store_true")  # print NCBI taxonomy identifiers
    #parser.add_argument("-B", action="store_true")  # print result in BIOM format
    #parser.add_argument("-C", type=str)  # print result in CAMI format (BioBoes format 0.9.1)
    #parser.add_argument("-q", action="store_true")  # print the full rank taxonomy
    #parser.add_argument("-A", action="store_true")  # print all taxonomic levels together
    #parser.add_argument("-k", type=str)  # taxonomic level [mOTU]

    # Algorithm options
    parser.add_argument("-g", type=int, default=3, choices=[1,2,3,4,5,6,7,8,9,10])  # number of marker genes cutoff
    parser.add_argument("-l", type=int, default=75)  # min length of the alignment (bp) [75]
    parser.add_argument("-t", type=int, default=1)  # number of thread [1]
    parser.add_argument("-v", type=int, default=1)  # verbosity level
    parser.add_argument("-y", type=str, default='INSERT_SCALED', choices=['INSERT_RAW', 'INSERT_NORM', 'INSERT_SCALED', 'BASE_RAW', 'BASE_NORM'])

    args = parser.parse_args(sys.argv[2:])

    # print usage and exit if no arguments are passed
    if sys.argv[2:] == []:
        parser.print_usage()
        shutdown(1)

    # converting string arguments to pathlib objects
    forward_files = [pathlib.Path(el) for el in args.f]
    reverse_files = [pathlib.Path(el) for el in args.r]
    unpaired_files = [pathlib.Path(el) for el in args.s]
    #mgcInputFile = pathlib.Path(args.m) if args.m != None else None
    #dbDir = pathlib.Path(args.db) if args.db != None else None
    motu_file = pathlib.Path(args.o)
    alignment_file = pathlib.Path(args.o + '.bam')
    mgc_file = pathlib.Path(args.o + '.mgc')
    startup()
    threads = args.t
    min_alignment_length = args.l
    samplename = args.n

    global motusdb
    motusdb = MotusDB(db_folder)
    global motusfiles
    motusfiles = MotusParameters()
    motusfiles.set_read_files(forward_files, reverse_files, unpaired_files, check_files=True)
    motusfiles.set_alignment_file(alignment_file, required_to_exist=False)
    motusfiles.set_mgc_file(mgc_file,required_to_exist=False)
    motusfiles.set_motu_file(motu_file, required_to_exist=False)
    motusfiles.set_sample_name(samplename)
    motusfiles.set_minimal_alignment_length(min_alignment_length)
    motusfiles.set_threads(threads)
    if not args.c:
        motusfiles.set_report_mode_rel_abundance()
    motusfiles.set_count_mode(args.y)
    motusfiles.set_minimal_number_of_mgcs(args.g)
    map_tax()

    calc_mgc()
    calc_motu()
    shutdown(0)


def parse_calc_mgc():
    parser = argparse.ArgumentParser(usage = '''Program: motus - a tool for marker gene-based OTU (mOTU) profiling
Version: 4.0.0
Reference: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand 
taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022). 
doi: https://doi.org/10.1186/s40168-022-01410-z
    
motus calc_mgc [options]

Input options:
   -i  FILE         provide the SAM or BAM input file (output of motus map_tax)

Output options:
   -o  FILE         output file name

Algorithm options:
   -l  INT          min length of the alignment (bp) [75]
   -v  INT          verbosity level: 1=error, 2=warning, 3=message, 4+=debugging [3]''', formatter_class=CapitalisedHelpFormatter,add_help=False)

    # Input options
    parser.add_argument("-i", type=str)  # provide a SAM or BAM input file (or list of files) output of motus map_tax

    # Output options
    parser.add_argument("-o", required=True)  # output file name [stdout]

    # Algorithm options
    parser.add_argument("-l", type=int, default=75)  # min length of the alignment (bp) [75]
    parser.add_argument("-v", type=int, default=1)  # verbosity level

    args = parser.parse_args(sys.argv[2:])

    # print usage and exit if no arguments are passed
    if sys.argv[2:] == []:
        parser.print_usage()
        shutdown(1)

    # converting string arguments into pathlib objects
    alignment_file = pathlib.Path(args.i)
    mgc_file = pathlib.Path(args.o)


    startup()
    min_alignment_length = args.l

    global motusdb
    motusdb = MotusDB(db_folder)
    global motusfiles
    motusfiles = MotusParameters()
    motusfiles.set_alignment_file(alignment_file, required_to_exist=True)
    motusfiles.set_mgc_file(mgc_file,required_to_exist=False)
    motusfiles.set_minimal_alignment_length(min_alignment_length)
    motusfiles.set_threads(1)
    calc_mgc()
    shutdown(0)


def parse_calc_motu():
    parser = argparse.ArgumentParser(usage = '''Program: motus - a tool for marker gene-based OTU (mOTU) profiling
Version: 4.0.0
Reference: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand 
taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022). 
doi: https://doi.org/10.1186/s40168-022-01410-z
    
motus calc_motu [options]

    Input options:
       -n  STR   sample name [unnamed sample]
       -i  FILE  provide the mgc abundance table (output of motus calc_mgc)
    
    Output options:
       -o  FILE  output file name 
       -c        print result as counts instead of relative abundances
    
    Algorithm options:
       -g   INT   number of marker genes cutoff: 1=higher recall, 6=higher precision [3]
       -v   INT   verbosity level: 1=error, 2=warning, 3=message, 4+=debugging [3]
       -y  STR    type of read counts [INSERT_SCALED]
                    Values: [INSERT_RAW, INSERT_NORM, INSERT_SCALED, BASE_RAW, BASE_NORM]
      
      ''', formatter_class=CapitalisedHelpFormatter,add_help=False)

    # Input options
    parser.add_argument("-n", type=str, default='unnamed sample')  # sample name
    parser.add_argument("-i", required=True)  # provide the mgc abundance table(output of motus calc_mgc)

    # Output options
    parser.add_argument("-o", required=True)  # output fil name [stdout]
    parser.add_argument("-c", action="store_true")  # print result as counts instead of realtive abundances
    parser.add_argument("-y", type=str, default='INSERT_SCALED', choices=['INSERT_RAW', 'INSERT_NORM', 'INSERT_SCALED', 'BASE_RAW', 'BASE_NORM'])
    parser.add_argument("-g", type=int, default=3,choices=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10])  # number of marker genes cutoff

    args = parser.parse_args(sys.argv[2:])

    # print usage and exit if no arguments are passed
    if sys.argv[2:] == []:
        parser.print_usage()
        shutdown(1)

    # converting string arguments into pathlib objects
    mgc_file = pathlib.Path(args.i)
    motu_file = pathlib.Path(args.o)
    startup()
    samplename = args.n

    global motusdb
    motusdb = MotusDB(db_folder)
    global motusfiles
    motusfiles = MotusParameters()
    motusfiles.set_mgc_file(mgc_file, required_to_exist=True)
    motusfiles.set_motu_file(motu_file, required_to_exist=False)
    motusfiles.set_sample_name(samplename)
    motusfiles.set_threads(1)
    motusfiles.set_count_mode(args.y)
    if not args.c:
        motusfiles.set_report_mode_rel_abundance()
    motusfiles.set_count_mode(args.y)
    motusfiles.set_minimal_number_of_mgcs(args.g)
    calc_motu()
    shutdown(0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage = '''Program: motus - a tool for marker gene-based OTU (mOTU) profiling
Version: 4.0.0
Reference: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand 
taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022). 
doi: https://doi.org/10.1186/s40168-022-01410-z
    
motus <command> [options]
    
    -- Taxonomic profiling
          profile     Perform taxonomic profiling (map_tax + calc_mgc + calc_motu) in a single step
          merge       Merge several taxonomic profiling results into one table

          map_tax     Map reads to the marker gene database
          calc_mgc    Calculate marker gene cluster (MGC) abundance
          calc_motu   Summarize MGC abundances into a mOTU profile

          prep_long   Prepare long reads to be profiled by mOTUs


    Type motus <command> to print the help menu for a specific command
    ''',formatter_class=CapitalisedHelpFormatter,add_help=False)

    parser.add_argument('command', choices=["profile", "merge", "map_tax", "calc_mgc", "calc_motu", "prep_long"])
    args: argparse.Namespace = parser.parse_args(sys.argv[1:2])
    if args.command == 'profile':
        parse_profile()
    elif args.command == 'merge':
        logging.error('Command merge not implemented yet')
    elif args.command == 'map_tax':
        parse_map_tax()
    elif args.command == 'calc_mgc':
        parse_calc_mgc()
    elif args.command == 'calc_motu':
        parse_calc_motu()
    elif args.command == 'prep_long':
        logging.error('Command prep_long not implemented yet')
    else:
        parser.print_usage()
        print(f'Unrecognized command {args}')
        shutdown(1)
    shutdown(0)



