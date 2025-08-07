from typing import List, Tuple, Dict, Set
import pathlib
import logging
import os
import mutils
import Bio.SeqIO.FastaIO as FastaIO
import Bio.SeqIO.QualityIO as QualityIO
import gzip
import csv

R1IDENTIFIER = '1'
R2IDENTIFIER = '2'
SIDENTIFIER = 'S'


class MotusParameters:
    _forward_files: List[pathlib.Path] = []
    _reverse_files: List[pathlib.Path] = []
    _unpaired_files: List[pathlib.Path] = []
    _read_files_were_checked: bool = False
    _alignment_file: pathlib.Path = None
    _temp_alignment_file: pathlib.Path = None

    _mgc_file: pathlib.Path = None
    _motu_file: pathlib.Path = None
    _motu_file_relab: pathlib.Path = None
    _inserts_file: pathlib.Path = None
    _samplename: str = None
    _min_alignment_length: int = 0
    _threads: int = 1
    _is_strict_db_mode = True

    _count_mode: str = 'INSERT_SCALED'


    def __init__(self):
        x = 0

        # nothing to do
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

    _min_mgcs: str = 3

    def is_strict_db_mode(self):
        return self._is_strict_db_mode

    def enable_lenient_mode(self):
        self._is_strict_db_mode = False

    def set_minimal_number_of_mgcs(self, min_mgcs: int) -> None:
        self._min_mgcs = min_mgcs

    def set_count_mode(self, count_mode: str) -> None:
        self._count_mode = count_mode

    def get_count_type(self) -> str:
        if 'NORM' in self._count_mode:
            return 'float'
        else:
            return 'int'

    def set_minimal_alignment_length(self, minimal_alignment_length: int):
        if minimal_alignment_length < 30:
            logging.error('Minimal alignment length is below aligner threshold. Pick a larger value. Quitting ...')
            mutils.shutdown(1)
        if minimal_alignment_length > 150:
            logging.warning('Minimal alignment length set to above average read length of metagenomic sequencing data.')
        self._min_alignment_length = int(minimal_alignment_length)

    def get_minimal_alignment_length(self) -> int:
        return self._min_alignment_length

    def set_threads(self, threads: int):
        if threads < 1:
            logging.error('Threads have to be at least 1')
            mutils.shutdown(1)
        if threads > os.cpu_count():
            logging.warning('Number of threads exceeds the total number of CPU cores.')
        self._threads = int(threads)

    def get_threads(self) -> int:
        return self._threads

    def set_sample_name(self, samplename: str):
        if len(samplename) == 0:
            logging.error('Sample name cannot be empty. Quitting')
            mutils.shutdown(1)
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
                mutils.shutdown(1)

    def get_inserts_file(self) -> pathlib.Path:
        self._inserts_file.parent.mkdir(exist_ok=True, parents=True)
        return self._inserts_file

    def set_inserts_file(self, inserts_file: pathlib.Path, required_to_exist=True) -> None:
        self._inserts_file = inserts_file
        if required_to_exist:
            if not inserts_file.exists():
                logging.error(f'Inserts file {inserts_file} does not exist. Shutting down ...')
                mutils.shutdown(1)

    def get_motu_file(self) -> pathlib.Path:
        self._motu_file.parent.mkdir(exist_ok=True, parents=True)
        return self._motu_file

    def get_motu_file_relab(self) -> pathlib.Path:
        self._motu_file_relab.parent.mkdir(exist_ok=True, parents=True)
        return self._motu_file_relab

    def set_motu_file(self, motu_file: pathlib.Path, required_to_exist=True) -> None:
        self._motu_file = motu_file
        self._motu_file_relab = pathlib.Path(str(motu_file) + '.relab')
        if required_to_exist:
            if not motu_file.exists():
                logging.error(f'mOTU file {motu_file} does not exist. Shutting down ...')
                mutils.shutdown(1)

    def set_alignment_file(self, alignment_file: pathlib.Path, required_to_exist=True) -> None:
        self._alignment_file = alignment_file
        self._temp_alignment_file = pathlib.Path(str(alignment_file) + '_tmp.bam')
        if required_to_exist:
            if not alignment_file.exists():
                logging.error(f'Alignment file {alignment_file} does not exist. Shutting down ...')
                mutils.shutdown(1)

        if not str(alignment_file).endswith('.bam'):
            logging.error(
                f'Alignment file {alignment_file} is/will be a BAM formatted file. Please set file suffix accordingly. Shutting down ...')
            mutils.shutdown(1)

    def get_read_files(self) -> List[Tuple[pathlib.Path, str]]:
        read_files = []
        for (r1_file, r2_file) in zip(self._forward_files, self._reverse_files):  # , strict=True):
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
            logging.error(f'Unknown file format: {reads_file}. Expecting a fasta or fastq file, can be gzipped.')
            mutils.shutdown(1)
        of.close()
        return reads

    def set_read_files(self, forward_files: List[pathlib.Path], reverse_files: List[pathlib.Path],
                       unpaired_files: List[pathlib.Path], check_files: bool = True) -> None:
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

        if check_files:
            files_that_dont_exist = []
            if len(forward_files + reverse_files + unpaired_files) == 0:
                logging.error('No input files defined with -f -r or -s. Quitting ...')
                mutils.shutdown(1)
            for f in forward_files + reverse_files + unpaired_files:
                if not f.exists():
                    files_that_dont_exist.append(f)
            if len(files_that_dont_exist) != 0:
                logging.error(f'Some read files dont exist: {files_that_dont_exist}')
                for f in files_that_dont_exist:
                    logging.error(f'\t{f}')
                mutils.shutdown(1)
            if len(set(forward_files + reverse_files + unpaired_files)) != len(
                    forward_files + reverse_files + unpaired_files):
                logging.error(f'Duplicated read files. Please submit every file only once. Shutting down ...')
                mutils.shutdown(1)
            if len(forward_files) != len(reverse_files):
                logging.error('Unequal number of files submitted with -r and -f. Quitting ...')
                mutils.shutdown(1)
            for (r1_file, r2_file) in zip(forward_files, reverse_files):  # , strict=True):
                r1_reads = self.get_first_1000_reads(r1_file)
                r2_reads = self.get_first_1000_reads(r2_file)
                r1_header = set([r[0] for r in r1_reads])
                r2_header = set([r[0] for r in r2_reads])
                if len(r1_header.symmetric_difference(r2_header)) != 0:
                    logging.error(f'Headers of reads are not identical. Shutting down ...')
                    logging.error(f'Differing read headers: {r1_header.symmetric_difference(r2_header)}')
                    logging.error(f'Differing read headers file 1: {r1_file}')
                    logging.error(f'Differing read headers file 2: {r2_file}')
                    mutils.shutdown(1)

            for u_file in unpaired_files:
                u_reads = self.get_first_1000_reads(u_file)

        self._forward_files = forward_files
        self._reverse_files = reverse_files
        self._unpaired_files = unpaired_files


class MotusDB:
    """
    The MotusDB class contains all relevant information that
    represents the current mOTUs marker gene database.
    E.g. (but not limited to)
    - MG - MGC - MOTU
    - Taxonomy per mOTU
    - Version
    - Member genome information
    """

    database_version: str = None
    database_date: str = None
    mgh_2_mgc: Dict[str, str] = {}
    mgh_2_mglength: Dict[str, int] = {}
    mgc_2_motu: Dict[str, str] = {}
    motus: Set[str] = set()
    blocklist_mg = set()
    motu_2_representative_genome_gtdb_tax = {}
    motu_2_mv_gtdb_tax = {}
    motu_2_representative = {}
    mgh_2_mg: Dict[str, str] = {}
    mgc_2_mg: Dict[str, str] = {}
    index_location: pathlib.Path = None
    _motus_core_mgs = ['COG0012','COG0016','COG0018','COG0172','COG0215','COG0495','COG0525','COG0533','COG0541','COG0552']
    _unassigned_motu_name = None

    motus_mv_taxonomy_file = None
    genome_metadata_file = None

    def __init__(self):
        x = 0

    def load_motus_db(self, mOTUsdb_folder: pathlib.Path, load=True) -> None:
        """Collect the mOTUs MGDB files, check their existence and,
        if the load parameter is set, load their contents into memory.

        Following files are expected:
        1. mOTUs.version --> holds the version of the database
        2. mOTUsNR.fasta.gz --> Marker gene sequences in gzipped fasta file
        3. mOTUsNR.fasta.gz.* --> the BWA index
        4. mOTUsv4.0.map.tsv.gz --> Link between mOTU, MGC and MG
        5. mOTUsv4.0.db.blocklist.gz --> contains MGs that should be removed from the alignment file
        6. mOTUsv4.0.gtdb.taxonomy.rep.tsv.gz --> The GTDB R220 annotation of mOTUs by their rep genome
        7. mOTUsv4.0.gtdb.taxonomy.80mv.tsv.gz --> The GTDB R220 annotation of mOTUs by their 80% majority vote
        8. mOTUsv4.0.genomes.tsv.gz --> Genome metadata
        9.

        :param mOTUsdb_folder:
        :return: None
        """
        logging.info('Loading database ... ')
        versions_file = mOTUsdb_folder.joinpath('mOTUsv4.0.db').resolve()
        index_files = [mOTUsdb_folder.joinpath(f).resolve() for f in ['mOTUsv4.0.db.fna.gz', 'mOTUsv4.0.db.fna.gz.amb', 'mOTUsv4.0.db.fna.gz.ann', 'mOTUsv4.0.db.fna.gz.bwt', 'mOTUsv4.0.db.fna.gz.pac', 'mOTUsv4.0.db.fna.gz.sa']]
        mgs_file = mOTUsdb_folder.joinpath('mOTUsv4.0.map.tsv.gz').resolve()
        blocklist_file = mOTUsdb_folder.joinpath('mOTUsv4.0.db.blocklist.gz').resolve()
        gtdb_taxonomy_file_reps = mOTUsdb_folder.joinpath('mOTUsv4.0.gtdb.taxonomy.rep.tsv.gz').resolve()
        gtdb_taxonomy_file_mv = mOTUsdb_folder.joinpath('mOTUsv4.0.gtdb.taxonomy.80mv.tsv.gz').resolve()
        genome_data_file = mOTUsdb_folder.joinpath('mOTUsv4.0.genomes.tsv.gz').resolve()

        with open(versions_file) as handle:
            self.database_version = handle.readline().strip().split()[-1]
            self.database_date = handle.readline().strip().split()[-1]
        self.index_location = index_files[0]
        for index_file in index_files + [mgs_file, blocklist_file, gtdb_taxonomy_file_reps, gtdb_taxonomy_file_mv, genome_data_file]:
            if not index_file.exists():
                logging.error(f'Database file {index_file} is missing. Quitting mOTUs...')
                mutils.shutdown(1)

        self.motus_mv_taxonomy_file = gtdb_taxonomy_file_mv
        self.genome_metadata_file = genome_data_file

        if load:
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
            with gzip.open(gtdb_taxonomy_file_reps, 'rt') as handle:
                #MOTU    GENOME  GTDBR220
                handle.readline()
                for line in handle:
                    [motu, representative, gtdb_taxonomy] = line.strip().split('\t')
                    self.motu_2_representative_genome_gtdb_tax[motu] = gtdb_taxonomy
                    self.motu_2_representative[motu] = representative
            with gzip.open(gtdb_taxonomy_file_mv, 'rt') as handle:
                #MOTU    GENOME  GTDBR220
                handle.readline()
                for line in handle:
                    [motu, gtdb_taxonomy] = line.strip().split('\t')
                    self.motu_2_mv_gtdb_tax[motu] = gtdb_taxonomy
            logging.info(f'Loading database finished. Version {self.database_version} (version date: {self.database_date}) contains {len(self.motus)} mOTUs, {len(self.mgc_2_motu)} markergeneclusters and {len(self.mgh_2_mglength)} markergenes.')

    def is_mg_blocked(self, mg: str) -> bool:
        """
        Checks if the markergene is in the blocklist
        and return True if that is the case

        Params:
            mg: a mOTUs markergene

        Return:
            True is markergene is in the blocklist. Does not check
            if markergene is in the mOTUs database
        """
        if mg in self.blocklist_mg:
            return True
        else:
            return False

    def get_full_version(self):
        """
        Get the full database version used in the
        header section of all mOTUs files.
        """

        return 'TOOL:' + mutils.MOTUS_VERSION + '_DB:' + self.database_version


    def get_tool_version(self):
        return mutils.MOTUS_VERSION
    def get_database_version(self):
        return self.database_version
    def get_full_sam_id(self):
        """
        Get the ID flag name for SAM/BAM header lines
        """

        return mutils.SAM_ID_FLAG

    def get_mg_by_mgc(self, mgc: str) -> str:
        """
        Get the COG of the markergenecluster

        Params:
            mgc: Name of the markergenecluster

        Returns:
            COG assoicated to markergenecluster

        """

        return self.mgc_2_mg[mgc]

    def is_unassigned_motu(self, motu: str) -> bool:
        """
        Checks whether the mOTU in the parameters
        is the unassigned mOTU

        Params:
            motu: a mOTUs

        Return:
            Whether the motu variable is the unassigned mOTU
        """

        if not self._unassigned_motu_name:
            logging.error('The unassigned mOTU was not set. This indicates a corrupted database. Please re-download database. Quitting...')
            mutils.shutdown(1)
        if motu == self._unassigned_motu_name:
            return True
        else:
            return False

    def get_unassigned_motu(self) -> str:
        """
        Get the name of the unassigned mOTU

        Returns:
            the name of the unassigned mOTU
        """

        if not self._unassigned_motu_name:
            logging.error('The unassigned mOTU was not set. This indicates a corrupted database. Please re-download database. Quitting...')
            mutils.shutdown(1)

        return self._unassigned_motu_name

    def get_motu_by_mgc(self, mgc: str) -> str:
        """
        Get the name of the mOTU associated with
        the mgc

        Params:
            mgc: The name of the markergenecluster

        Returns:
            The name of the associated mOTU
        """

        return self.mgc_2_motu[mgc]

    def get_bwa_index(self) -> pathlib.Path:
        """
        Get the location of the bwa index

        Returns:
            The Location of the bwa index
        """

        return self.index_location

    def get_mgc_by_mg(self, mgh: str) -> str:
        """
        Report markergenecluster by a markergene

        Params:
            mgh: the name of the markergene

        Returns:
            the associated markergenecluster
        """

        return self.mgh_2_mgc[mgh]

    def get_length_by_mg(self, mgh: str) -> int:
        """
        Get the length of the markergeneheader

        Params:
            mgh: a markergeneheader

        Returns:
            The length of the markergeneheader
        """

        return self.mgh_2_mglength[mgh]

    def get_mg_by_mgh(self, mgh: str) -> str:
        """
        Get the markergene (COG) associated with the
        markergeneheader

        Params:
            mgh: the markergeneheader

        Return:
            The COG associated with the markergeneheader
        """

        return self.mgh_2_mg[mgh]

    def get_core_motus_mgs(self) -> List[str]:
        """
        Return a list of the (currently) 10 mOTUs
        markergenes

        Return:
            List of 10 mOTUs markergenes
        """

        return self._motus_core_mgs





MOTUS_PARAMETERS = MotusParameters()
MOTUS_DB = MotusDB()