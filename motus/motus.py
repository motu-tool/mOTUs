import pysam
import Bio.SeqIO.FastaIO as FastaIO
import Bio.SeqIO.QualityIO as QualityIO
import logging
import pathlib
import sys
import csv
import subprocess
import gzip
from typing import List, Dict, Set, Tuple





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


class MotusFiles():

    _forward_files: List[pathlib.Path] = []
    _reverse_files: List[pathlib.Path] = []
    _unpaired_files: List[pathlib.Path] = []
    _read_files_were_checked: bool= False
    _alignment_file: pathlib.Path = None
    _temp_alignment_file: pathlib.Path = None

    _mgc_file = None
    _motu_file = None

    def set_alignment_file(self, alignment_file: pathlib.Path, required_to_exist=True) -> None:
        self._alignment_file = alignment_file
        self._temp_alignment_file = pathlib.Path(str(alignment_file) + '_tmp.bam')
        if required_to_exist:
            if not alignment_file.exists():
                logging.error(f'Alignment file {alignment_file} does not exist. Shutting down ...')
                shutdown(1)

        if not str(alignment_file).endswith('.bam'):
            logging.error(f'Alignment file {alignment_file} is/will be a BAM formatted file. Please set file suffix accordingly. Shutting down ...')
            shutdown()


    def get_read_files(self):
        read_files = []
        for (r1_file, r2_file) in zip(self._forward_files, self._reverse_files, strict=True):
            read_files.append((r1_file, '/1'))
            read_files.append((r2_file, '/2'))
        for u_file in self._unpaired_files:
            read_files.append((u_file, '/S'))
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
        self._forward_files = forward_files
        self._reverse_files = reverse_files
        self._unpaired_files = unpaired_files

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
            for f in forward_files + reverse_files + unpaired_files:
                if not f.exists():
                    files_that_dont_exist.append(f)
            if len(files_that_dont_exist) != 0:
                logging.error(f'Some read files dont exist:')
                for f in files_that_dont_exist:
                    logging.error(f'\t{f}')
                shutdown(1)
            if len(set(forward_files + reverse_files + unpaired_files)) != len(forward_files + reverse_files + unpaired_files):
                logging.error(f'Duplicated read files. Please submit every file only once. Shutting down ...')
                shutdown(1)
            # check if correct file ending
            # check if reads are paired
            for (r1_file, r2_file) in zip(forward_files, reverse_files, strict=True):
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
    motus: Set[str] = set()
    #motu_2_taxonomy: Dict[str, str] = {}
    index_location: str = None


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
        versions_file = mOTUsdb_folder.joinpath('mOTUsv4.version').resolve()
        index_files = [mOTUsdb_folder.joinpath(f).resolve() for f in ['mOTUsv4.fna', 'mOTUsv4.fna.amb', 'mOTUsv4.fna.ann', 'mOTUsv4.fna.bwt', 'mOTUsv4.fna.pac', 'mOTUsv4.fna.sa']]
        mgs_file = mOTUsdb_folder.joinpath('mOTUsv4.mgs.tsv').resolve()

        with open(versions_file) as handle:
            self.database_version = handle.readline().strip()
        self.index_location = index_files[0]
        for index_file in index_files:
            if not index_file.exists():
                logging.error(f'Database file {index_file} is missing. Quitting mOTUs...')
                shutdown(1)
        with open(mgs_file) as handle:
            for entry in  csv.DictReader(handle, delimiter='\t'):
                self.mg_2_mgc[entry['MG']] = entry['MGC']
                self.mg_2_mglength[entry['MG']] = int(entry['LENGTH'])
                self.mgc_2_motu[entry['MGC']] = entry['#MOTU']
                self.motus.add(entry['#MOTU'])
        logging.info(f'Loading database finished. Version {self.database_version} contains {len(self.motus)} mOTUs, {len(self.mgc_2_motu)} markergeneclusters and {len(self.mg_2_mglength)} markergenes.')


    def get_bwa_index(self):
        return self.index_location







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


def map_tax(threads: int = 1, minlength: int = 45) -> None:
    """
    Takes a list of forward/reverse/unpaired read files and aligns them against the mOTUs database using the number of specified threads.
    Alignments will be filtered by 97% identity and the defined minimal alignment length. The resulting alignments will be stored in the
    sorted BAM file which is either specified as a parameter or as a temporary file.

    # :param forward_files: List of forward read files in fasta or fastq format, optionally gzipped. Has to match the reverse files.
    # :param reverse_files: List of reverse read files in fasta or fastq format, optionally gzipped. Has to match the forward files.
    # :param unpaired_files: List of single/merged read files in fasta or fastq format, optionally gzipped.
    # :param bam_file: [Optional] Location of intermediate bam file. In case of None, a temporary file will be created in /tmp/. Default=[None]
    :param threads: [Optional] Number of threads used for the alignment of reads against the mOTUs database
    :param minlength: [Optional] minimal length of alignment. Default=[45]
    :return: None
    """
    logging.info('Starting mOTUs - map_tax routine - Alignment against the mOTUs database ... ')
    min_perc_id = 97.0

    temp_bam_file = motusfiles.get_temporary_alignment_file()
    temp_bam_file_handle = None



    total_reads = 0
    total_mapped_reads = 0

    for readsfile, orientation in motusfiles.get_read_files():
        total_reads_this_file = 0
        total_mapped_reads_this_file = set()
        logging.info(f'Aligning {readsfile}')
        command = f'bwa mem -a -t {threads} {motusdb.get_bwa_index()} {readsfile}'
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        in_bam_file_handle = pysam.AlignmentFile(process.stdout, 'rb')
        if not temp_bam_file_handle:
            temp_bam_file_handle = pysam.AlignmentFile(temp_bam_file, "wb", template=in_bam_file_handle)

        for record in in_bam_file_handle:


            if record.is_unmapped:
                total_reads_this_file += 1
                continue
            else:
                if not record.is_secondary and not record.is_supplementary:
                    total_reads_this_file += 1
                alnlength = sum(record.get_cigar_stats()[0][0:3])
                if alnlength < minlength:
                    continue
                query_covered_bases = sum(record.get_cigar_stats()[0][0:2])
                query_length = record.infer_read_length()
                mismatches = record.get_tag('NM')
                percid = (alnlength - mismatches) / float(alnlength) * 100.0
                percid = round(percid, 2)
                if min_perc_id > percid:
                    continue
                qcov = query_covered_bases / float(query_length)
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
    return_code = process.wait()
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

motufiles = None
motusdb = None
if __name__ == '__main__':
    startup()
    db_folder = pathlib.Path('/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/DEV/MOTUSv4/speci/speci_workfolder/motus/15database_build/')
    forward_files = [pathlib.Path('/nfs/cds-shini.ethz.ch/exports/biol_micro_cds_gr_sunagawa/SequenceStorage/analysis/GENERAL-INTERNAL/WIRB19-1/METAG/WIRB19-1_SAMEA4817971_METAG/qc_v3/ERR2726527/WIRB19-1_SAMEA4817971_METAG_ERR2726527.1.fq.gz'), pathlib.Path('/nfs/cds-shini.ethz.ch/exports/biol_micro_cds_gr_sunagawa/SequenceStorage/analysis/GENERAL-INTERNAL/WIRB19-1/METAG/WIRB19-1_SAMEA4817971_METAG/qc_v3/ERR2726529/WIRB19-1_SAMEA4817971_METAG_ERR2726529.1.fq.gz')]
    reverse_files = [pathlib.Path('/nfs/cds-shini.ethz.ch/exports/biol_micro_cds_gr_sunagawa/SequenceStorage/analysis/GENERAL-INTERNAL/WIRB19-1/METAG/WIRB19-1_SAMEA4817971_METAG/qc_v3/ERR2726527/WIRB19-1_SAMEA4817971_METAG_ERR2726527.2.fq.gz'), pathlib.Path('/nfs/cds-shini.ethz.ch/exports/biol_micro_cds_gr_sunagawa/SequenceStorage/analysis/GENERAL-INTERNAL/WIRB19-1/METAG/WIRB19-1_SAMEA4817971_METAG/qc_v3/ERR2726529/WIRB19-1_SAMEA4817971_METAG_ERR2726529.2.fq.gz')]
    unpaired_files = [pathlib.Path('/nfs/cds-shini.ethz.ch/exports/biol_micro_cds_gr_sunagawa/SequenceStorage/analysis/GENERAL-INTERNAL/WIRB19-1/METAG/WIRB19-1_SAMEA4817971_METAG/qc_v3/ERR2726527/WIRB19-1_SAMEA4817971_METAG_ERR2726527.s.fq.gz'), pathlib.Path('/nfs/cds-shini.ethz.ch/exports/biol_micro_cds_gr_sunagawa/SequenceStorage/analysis/GENERAL-INTERNAL/WIRB19-1/METAG/WIRB19-1_SAMEA4817971_METAG/qc_v3/ERR2726529/WIRB19-1_SAMEA4817971_METAG_ERR2726529.s.fq.gz')]
    output_bam_file = pathlib.Path('/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/DEV/MOTUSv4/example.bam')
    threads = 32
    min_alignment_length = 45
    alignment_file = pathlib.Path('/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/DEV/MOTUSv4/wirb_test.bam')
    motusdb = MotusDB(db_folder)
    motusfiles = MotusFiles()
    motusfiles.set_read_files(forward_files, reverse_files, unpaired_files, check_files=True)
    motusfiles.set_alignment_file(alignment_file, required_to_exist=False)
    map_tax(threads, minlength=min_alignment_length)
    shutdown(0)