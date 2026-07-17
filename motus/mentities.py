from typing import List, Tuple, Dict, Set, Self
import pathlib
import logging
import os
from motus import mutils
import Bio.SeqIO.FastaIO as FastaIO
import Bio.SeqIO.QualityIO as QualityIO
import gzip
import polars as pl
import collections
import statistics

R1IDENTIFIER = '1'
R2IDENTIFIER = '2'
SIDENTIFIER = 'S'


class GenomeLocator:

    def __init__(self, genome_metadata_file: pathlib.Path) -> None:
        self._valid_genomes = set()
        self._representative_genomes = set()

        df = pl.scan_csv(genome_metadata_file, separator="\t", has_header=True, infer_schema_length=0).select(["GENOME", 'MOTU4_STATUS']).collect()
        for row in df.iter_rows():
            [genome, motu4_status] = row
            self._valid_genomes.add(genome)
            if 'representative' in motu4_status:
                self._representative_genomes.add(genome)

    def is_represenative_genome(self, genome: str) -> bool:
        if genome not in self._valid_genomes:
            logging.error(f'Genome "{genome}" does not exist. Quitting')
            mutils.shutdown(1)
        return genome in self._representative_genomes

    def get_genome_path(self, genome: str, file_type: str = 'genome') -> str:
        if genome not in self._valid_genomes:
            logging.error(f'Genome "{genome}" does not exist. Quitting')
            mutils.shutdown(1)
        return f'{mutils.MOTUS_GENOME_API_BASE_URL}/{genome}/download?file_type={file_type}'



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

    _write_relative_abundances = False

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

    _min_mgcs: int = 3

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

    def set_write_relabundances(self):
        self._write_relative_abundances = True

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
                       unpaired_files: List[pathlib.Path], check_files: bool = True, skip_pair_check: bool = False) -> None:
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
                if skip_pair_check:
                    continue
                r1_reads = self.get_first_1000_reads(r1_file)
                r2_reads = self.get_first_1000_reads(r2_file)
                read1 = r1_reads[0][0]
                read1_norm = mutils.normalize_header(read1)
                if read1 != read1_norm:
                    logging.warning('Detected /1 /2 on the end of readheaders. Removing those for downstream analysis')
                r1_header = set(mutils.normalize_header(r[0]) for r in r1_reads)
                r2_header = set(mutils.normalize_header(r[0]) for r in r2_reads)
                if len(r1_header.symmetric_difference(r2_header)) != 0:
                    logging.error(f'Read headers do not match between paired files. Shutting down ...')
                    logging.error(f'Differing read headers: {r1_header.symmetric_difference(r2_header)}')
                    logging.error(f'Forward file: {r1_file}')
                    logging.error(f'Reverse file: {r2_file}')
                    logging.error(f'If your reads are unsorted or contain singletons, use --skip-pair-check to bypass this validation.')
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
    _has_unassigned_motu = False
    _mOTUsdb_folder: pathlib.Path = None

    motus_mv_taxonomy_file = None
    genome_metadata_file = None

    def __init__(self):
        x = 0




    def load_motus_db(self, mOTUsdb_folder: pathlib.Path, load=True) -> None:
        """Collect the mOTUs MGDB files, check their existence and,
        if the load parameter is set, load their contents into memory.

        The pattern of the file determines the version

        e.g 
        mOTUsv4.0.db --> version 4.0 of the database
        mOTUsv4.1-toy.db --> version 4.1-toy
        mOTUsv4.1.db --> version 4.0

        version shortened to XXX

        Following files are expected:
        1. mOTUsvXXX.db --> information on database, including version
        2. mOTUsvXXX.db.fna.gz --> Marker gene sequences in gzipped fasta file
        3. mOTUsvXXX.db.fna.gz.* --> the BWA index
        4. mOTUsvXXX.map.tsv.gz --> Link between mOTU, MGC and MG
        5. mOTUsvXXX.db.blocklist.gz --> contains MGs that should be removed from the alignment file
        6. mOTUsvXXX.gtdb.taxonomy.rep.tsv.gz --> The GTDB R220 annotation of mOTUs by their rep genome
        7. mOTUsvXXX.gtdb.taxonomy.80mv.tsv.gz --> The GTDB R220 annotation of mOTUs by their 80% majority vote
        8. mOTUsvXXX.genomes.tsv.gz --> Genome metadata

        9. mOTUsvXXX.annotation.db --> the mOTUs annotation sql database


        :param mOTUsdb_folder:
        :return: None
        """
        logging.info(f'Loading database from {mOTUsdb_folder} ...')

        db_marker = mOTUsdb_folder / 'db_mOTU.downloaded'
        if not db_marker.exists():
            logging.error('mOTUs marker gene database not downloaded. Download database with "motus downloadMGDB"')
            mutils.shutdown(1)

        # Discover the database version from the mOTUsv*.db file on disk,
        # excluding mOTUsv*.annotation.db which lives in the same folder.
        db_files = sorted(
            f for f in mOTUsdb_folder.glob('mOTUsv*.db')
            if not f.name.endswith('.annotation.db')
        )
        if len(db_files) == 0:
            logging.error(f'No mOTUsv*.db file found in {mOTUsdb_folder}. Is the database complete?')
            mutils.shutdown(1)
        if len(db_files) > 1:
            logging.error(f'Multiple mOTUsv*.db files found in {mOTUsdb_folder}: {[f.name for f in db_files]}. Cannot determine version.')
            mutils.shutdown(1)
        versions_file = db_files[0].resolve()
        db_prefix_from_filename = versions_file.stem[len('mOTUsv'):]  # e.g. '4.0' or '4.1-toy'

        self._mOTUsdb_folder = mOTUsdb_folder

        with open(versions_file) as handle:
            self.database_version = handle.readline().strip().split()[-1]
            self.database_date = handle.readline().strip().split()[-1]

        if self.database_version != db_prefix_from_filename:
            logging.error(f'Database version mismatch: filename "{versions_file.name}" implies version "{db_prefix_from_filename}" but the file reports version "{self.database_version}". Quitting ...')
            mutils.shutdown(1)
        logging.info(f'mOTUs database version: {self.database_version}')
        v = f'mOTUsv{self.database_version}'
        index_files = [mOTUsdb_folder.joinpath(f).resolve() for f in [
            f'{v}.db.fna.gz', f'{v}.db.fna.gz.amb', f'{v}.db.fna.gz.ann',
            f'{v}.db.fna.gz.bwt', f'{v}.db.fna.gz.pac', f'{v}.db.fna.gz.sa'
        ]]
        mgs_file = mOTUsdb_folder.joinpath(f'{v}.map.tsv.gz').resolve()
        blocklist_file = mOTUsdb_folder.joinpath(f'{v}.db.blocklist.gz').resolve()
        gtdb_taxonomy_file_reps = mOTUsdb_folder.joinpath(f'{v}.gtdb.taxonomy.rep.tsv.gz').resolve()
        gtdb_taxonomy_file_mv = mOTUsdb_folder.joinpath(f'{v}.gtdb.taxonomy.80mv.tsv.gz').resolve()
        genome_data_file = mOTUsdb_folder.joinpath(f'{v}.genomes.tsv.gz').resolve()
        self.index_location = index_files[0]
        for index_file in index_files + [mgs_file, blocklist_file, gtdb_taxonomy_file_reps, gtdb_taxonomy_file_mv, genome_data_file]:
            if not index_file.exists():
                logging.error(f'Database file {index_file} is missing. Quitting mOTUs...')
                mutils.shutdown(1)

        self.motus_mv_taxonomy_file = gtdb_taxonomy_file_mv
        self.genome_metadata_file = genome_data_file

        if load:
            df = pl.scan_csv(mgs_file, separator="\t", has_header=True, infer_schema_length=0,
                             schema_overrides={'LENGTH': pl.UInt32}).select(
                ['#MOTU', 'MGC', 'COG', 'MG', 'LENGTH']).collect()
            mg_col    = df['MG'].to_list()
            mgc_col   = df['MGC'].to_list()
            cog_col   = df['COG'].to_list()
            motu_col  = df['#MOTU'].to_list()
            len_col   = df['LENGTH'].to_list()
            self.mgh_2_mgc      = dict(zip(mg_col, mgc_col))
            self.mgh_2_mglength = dict(zip(mg_col, len_col))
            self.mgc_2_motu     = dict(zip(mgc_col, motu_col))
            self.mgh_2_mg       = dict(zip(mg_col, cog_col))
            self.mgc_2_mg       = dict(zip(mgc_col, cog_col))
            self.motus          = set(motu_col)
            unassigned = df.filter(pl.col('#MOTU').str.contains('unassigned'))['#MOTU']
            if len(unassigned) > 0:
                self._unassigned_motu_name = unassigned[0]
                self._has_unassigned_motu = True

            with gzip.open(blocklist_file, 'rt') as handle:
                for line in handle:
                    self.blocklist_mg.add(line.strip())

            df_reps = pl.scan_csv(gtdb_taxonomy_file_reps, separator="\t", has_header=True,
                                  infer_schema_length=0).collect()
            df_reps = df_reps.rename({df_reps.columns[0]: df_reps.columns[0].lstrip('#')})
            gtdb_col_reps = 'GTDB' if 'GTDB' in df_reps.columns else 'GTDBR220'
            self.motu_2_representative_genome_gtdb_tax = dict(zip(
                df_reps['MOTU'].to_list(), df_reps[gtdb_col_reps].to_list()))
            self.motu_2_representative = dict(zip(
                df_reps['MOTU'].to_list(), df_reps['GENOME'].to_list()))

            df_mv = pl.scan_csv(gtdb_taxonomy_file_mv, separator="\t", has_header=True,
                                infer_schema_length=0).collect()
            df_mv = df_mv.rename({df_mv.columns[0]: df_mv.columns[0].lstrip('#')})
            gtdb_col_mv = 'GTDB' if 'GTDB' in df_mv.columns else 'GTDBR220'
            self.motu_2_mv_gtdb_tax = dict(zip(
                df_mv['MOTU'].to_list(), df_mv[gtdb_col_mv].to_list()))
            logging.info(f'Loading database finished. Version {self.database_version} (version date: {self.database_date}) contains {len(self.motus)} mOTUs, {len(self.mgc_2_motu)} markergeneclusters and {len(self.mgh_2_mglength)} markergenes.')


    def get_motu_2_median_mgc_gene_length(self) -> Dict[str, Dict[str, int]]:
        """Get for every MGC the median
        gene length and store by motu_2_mgc dict
        Caution: Ignores the unassigned mOTU

        Returns:
            Dict[str, Dict[str, int]]: motu --> mgcs --> median_length
        """

        mgc_2_lengths = collections.defaultdict(list)
        for mgh, mgc in self.mgh_2_mgc.items():
            if 'unassigned' in mgc:
                continue
            mgh_length = self.mgh_2_mglength[mgh]
            mgc_2_lengths[mgc].append(mgh_length)
        motu_2_cog = collections.defaultdict(lambda: {})

        for mgc, lengths in mgc_2_lengths.items():
            cog = self.mgc_2_mg[mgc]
            motu = self.mgc_2_motu[mgc]
            motu_2_cog[motu][cog] = int(statistics.median(lengths))

        return motu_2_cog


    def get_mv_tax_for_motu(self, motu: str) -> str:
        if 'unassigned' in motu:
            return 'd__;p___A;c__;o__;f__;g__;s__'
        else:
            return self.motu_2_mv_gtdb_tax[motu]
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

    def is_valid_motu(self, motu: str):
        """ checks if the text in the motu variable is
        an actual mOTU in the current database

        Params:
            motu: Name of a mOTU
        Returns:
            if the motu is an actual motu
        """
        return motu in self.motus

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

        if not self._has_unassigned_motu:
            return False
        return motu == self._unassigned_motu_name

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

class SinglemOTUsFile:
    # metadata
    _count_mode = None
    _count_mode_options = ['INSERT_RAW', 'INSERT_NORM', 'INSERT_SCALED', 'BASE_RAW', 'BASE_NORM']
    _min_mgcs = None
    _min_mgcs_options = [1,2,3,4,5,6,7,8,9,10]
    _database_version = None
    _tool_version = None
    _min_alignment_length = None
    _value_type = None
    _value_type_options = ['counts', 'relative_abundances']
    # descriptors
    _samplename = None

    # data
    _motu_2_values = None



    def __init__(self, motu_2_values: Dict[str, float], min_alignment_length: int, min_mgcs: int, count_mode:str, database_version: str, tool_version: str, value_type: str, samplename: str) -> None:
        self._motu_2_values = {}

        if count_mode not in self._count_mode_options:
            logging.error(f'Unknown count mode: {count_mode}. Options: {self._count_mode_options}')
            mutils.shutdown(1)
        self._count_mode = count_mode

        if value_type not in self._value_type_options:
            logging.error(f'Unknown value type: {value_type}. Options: {self._value_type_options}')
            mutils.shutdown(1)
        self._value_type = value_type

        if min_mgcs not in self._min_mgcs_options:
            logging.error(f'Unknown number of MGCs: {min_mgcs}. Options: {self._min_mgcs_options}')
            mutils.shutdown(1)
        self._min_mgcs = min_mgcs

        if database_version != MOTUS_DB.get_database_version():
            logging.error(f'Invalid database version: {database_version}. Options: {MOTUS_DB.get_database_version()}')
            mutils.shutdown(1)
        self._database_version = database_version

        if tool_version != MOTUS_DB.get_tool_version():
            logging.error(f'Invalid tool version: {tool_version}. Options: {MOTUS_DB.get_tool_version()}')
            mutils.shutdown(1)
        self._tool_version = tool_version

        self._min_alignment_length = min_alignment_length

        self._samplename = samplename

        for motu, value in motu_2_values.items():
            if value != 0.0:
                self._motu_2_values[motu] = value
                if not MOTUS_DB.is_valid_motu(motu):
                    logging.error(f'Unknown mOTU {motu}')
                    mutils.shutdown(1)




    @staticmethod
    def read_motus_file(motus_file: pathlib.Path) -> 'SinglemOTUsFile':
        '''Read the input path into
        a SinglemOTUsFile object.
        Does also perform standard checks
        Params:
            motus_file: A file path to the single motus file
        Returns:
            A SinglemOTUsFile object

        '''

        if mutils.is_gzipped(motus_file):
            fo = gzip.open(motus_file, 'rt')
        else:
            fo = open(motus_file, 'r')
        # #tool_version=4.0.2     database_version=4.0    min_alignment_length=110        min_mgcs=3      count_mode=INSERT_NORM  value_type=counts
        header = fo.readline()
        [tool_version, database_version, min_alignment_length, min_mgcs, count_mode, value_type] = header.strip().split('\t')
        tool_version = tool_version.replace('#tool_version=', '')
        database_version = database_version.replace('database_version=', '')
        min_alignment_length = int(min_alignment_length.replace('min_alignment_length=', ''))
        min_mgcs = int(min_mgcs.replace('min_mgcs=', ''))
        count_mode = count_mode.replace('count_mode=', '')
        value_type = value_type.replace('value_type=', '')

        header2 = fo.readline()
        [motu, taxonomy, samplename] = header2.strip().split('\t')

        motu_2_values = {}
        for line in fo:
            [motu, _, value] = line.strip().split('\t')
            motu_2_values[motu] = float(value)
        fo.close()

        smf = SinglemOTUsFile(motu_2_values, min_alignment_length, min_mgcs, count_mode, database_version, tool_version, value_type, samplename)
        return smf



    def get_motus_file_header(self):
        tmp = f'#tool_version={self._tool_version}\tdatabase_version={self._database_version}\tmin_alignment_length={self._min_alignment_length}\t'
        tmp = tmp + f'min_mgcs={self._min_mgcs}\tcount_mode={self._count_mode}\tvalue_type={self._value_type}'
        return tmp


    def write_to_file(self, output_file: pathlib.Path) -> None:
        '''write the contents of this object to a
        the output file
        '''
        sorted_motus = sorted(self._motu_2_values.keys())
        header_line = self.get_motus_file_header()
        with open(output_file, 'w') as outhandle:
            outhandle.write(f'{header_line}\n')
            outhandle.write(f'mOTU\tTaxonomy\t{self._samplename}\n')
            for motu in sorted_motus:
                tax = MOTUS_DB.get_mv_tax_for_motu(motu)
                value = self._motu_2_values[motu]
                if self._value_type == 'relative_abundances':
                    report_value = '{number:.{digits}f}'.format(number=value, digits=8)
                else:
                    if 'NORM' in self._count_mode:
                        report_value = '{number:.{digits}f}'.format(number=value, digits=8)
                    else:
                        report_value = round(value)
                outhandle.write(f'{motu}\t{tax}\t{report_value}\n')






    def get_relative_abundances(self) -> Self:
        '''Uses the data in this object to create a new
        or existing SingleMotus file object which has
        relative abundances instead of counts as values

        '''
        if self._value_type == 'relative_abundances':
            return self

        motu_2_count = self._motu_2_values
        motu_2_relab = {}
        total = float(sum(motu_2_count.values()))
        if total == 0.0:
            return self
        for motu, count in motu_2_count.items():
            relab = float(count) / total
            motu_2_relab[motu] = relab

        smf = SinglemOTUsFile(motu_2_relab, self._min_alignment_length, self._min_mgcs, self._count_mode, self._database_version, self._tool_version, 'relative_abundances', self._samplename)
        return smf


class MergedmOTUsFile:

    def __init__(self, motus_files: List[pathlib.Path]) -> None:
        self._singlemotusfiles = {}
        # check that there are >1 motus files
        if not motus_files or len(motus_files) < 2:
            logging.error(f'Provide at least 2 mOTUs files for merging.')
            mutils.shutdown(1)
        # load the individual motus files
        singlemotusfiles = []
        for motus_file in motus_files:
            smf = SinglemOTUsFile.read_motus_file(motus_file)
            singlemotusfiles.append(smf)

        # check that they're mergable
        # 1. headers have to be the same

        fields = {
            "count_mode": set([smf._count_mode for smf in singlemotusfiles]),
            "min_mgcs": set([smf._min_mgcs for smf in singlemotusfiles]),
            "database_version": set([smf._database_version for smf in singlemotusfiles]),
            "tool_version": set([smf._tool_version for smf in singlemotusfiles]),
            "min_alignment_length": set([smf._min_alignment_length for smf in singlemotusfiles]),
            "value_type": set([smf._value_type for smf in singlemotusfiles]),
        }
        samplenames = set([smf._samplename for smf in singlemotusfiles])

        killnow = False
        for name, values in fields.items():
            if len(values) > 1:
                logging.error(f"Different {name} detected: {values}")
                killnow = True
        # 2. different sample names
        if len(samplenames) != len(singlemotusfiles):
            logging.error(f"Duplicated sample names detected.")
            killnow = True
        if killnow:
            mutils.shutdown(1)
        # Then merge into this object
        for smf in singlemotusfiles:
            self._singlemotusfiles[smf._samplename] = smf

    def write_to_file(self, output_file: pathlib.Path) -> None:

        sorted_motus = set()
        random_smf = None
        for _, smf in self._singlemotusfiles.items():
            random_smf = smf

            for k, v in smf._motu_2_values.items():
                sorted_motus.add(k)


        sorted_motus = sorted(list(sorted_motus))
        sorted_samples = sorted(self._singlemotusfiles.keys())

        with open(output_file, 'w') as outhandle:
            outhandle.write(f'{random_smf.get_motus_file_header()}\n')
            xx = "\t".join(sorted_samples)
            outhandle.write(f'mOTU\tTaxonomy\t{xx}\n')
            for motu in sorted_motus:
                tax = MOTUS_DB.get_mv_tax_for_motu(motu)
                values = []
                for samplename in sorted_samples:
                    values.append(self._singlemotusfiles[samplename]._motu_2_values.get(motu, 0.0))
                if random_smf._value_type == 'relative_abundances':
                    report_values = ['{number:.{digits}f}'.format(number=v, digits=8) for v in values]
                else:
                    if 'NORM' in random_smf._count_mode:
                        report_values = ['{number:.{digits}f}'.format(number=v, digits=8) for v in values]
                    else:
                        report_values = [round(v) for v in values]
                tmp = [motu, tax] + report_values
                tmp = '\t'.join([str(t) for t in tmp])
                outhandle.write(tmp + '\n')


MOTUS_PARAMETERS = MotusParameters()
MOTUS_DB = MotusDB() 