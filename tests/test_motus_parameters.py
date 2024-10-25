import unittest
from motus import motus
from unittest.mock import patch, mock_open
import pathlib
import os


class TestMotusParametersSetGet(unittest.TestCase):

    def setUp(self):
        self.obj = motus.MotusParameters()

        # set initial values to ensure a known state
        self.obj = motus.MotusParameters()
        self.obj._min_alignment_length = 100
        self.obj._threads = 4
        self.obj._samplename = "sample_01"
        self.obj._is_strict_db_mode = True
        self.obj._min_mgcs = 3
        self.obj._count_mode = 'raw'
        self.obj._mgc_file = pathlib.Path("mgc_file")
        self.obj._inserts_file = pathlib.Path("inserts_file")
        self.obj._motu_file = pathlib.Path("motu_file")
        self.obj._motu_file_rel_ab = pathlib.Path("motu_file.relab")
        self.obj._alignment_file = pathlib.Path("alignment_file.bam")
        self.obj._temp_alignment_file = pathlib.Path("alignment_file_tmp.bam")
        self.obj._forward_files = [pathlib.Path("forward_file1"), pathlib.Path("forward_file2")]
        self.obj._reverse_files = [pathlib.Path("reverse_file1"), pathlib.Path("reverse_file2")]
        self.obj._unpaired_files = [pathlib.Path("unpaired_file")]

    def test_is_strict_db_mode(self):
        # test strict mode is initially True
        self.assertTrue(self.obj.is_strict_db_mode())

    def test_enable_lenient_mode(self):
        # test that enable_lenient_mode sets _is_strict_db_mode to False
        self.obj.enable_lenient_mode()
        self.assertFalse(self.obj.is_strict_db_mode())

    def test_set_minimal_number_of_mgcs(self):
        # test that set_minimal_number_of_mgcs sets the correct value
        self.obj.set_minimal_number_of_mgcs(10)
        self.assertEqual(self.obj._min_mgcs, 10)

        self.obj.set_minimal_number_of_mgcs(5)
        self.assertEqual(self.obj._min_mgcs, 5)

    def test_set_count_mode(self):
        # test that set_count_mode sets the correct count mode
        self.obj.set_count_mode('raw')
        self.assertEqual(self.obj._count_mode, 'raw')

        self.obj.set_count_mode('NORM')
        self.assertEqual(self.obj._count_mode, 'NORM')

    def test_get_count_type_norm(self):
        # test that set_count_mode sets the correct count mode
        self.obj._count_mode = "NORM"

        self.assertEqual(self.obj.get_count_type(), 'float')

    def test_get_count_type_other(self):
        # test that set_count_mode sets the correct count mode
        self.obj._count_mode = "OTHER"

        self.assertEqual(self.obj.get_count_type(), 'int')

    def test_set_minimal_alignment_length(self):
        # test that set_minimal_alignment_length
        minimal_alignment_length = 50
        self.obj.set_minimal_alignment_length(minimal_alignment_length)

        self.assertEqual(self.obj._min_alignment_length, 50)

    def test_set_minimal_alignment_length_below(self):
        # test that set_minimal_alignment_length fails for values below 30
        minimal_alignment_length = 29

        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.set_minimal_alignment_length(minimal_alignment_length)
            # check that log message is accurate
            self.assertEqual(log_capture.output, [
                'ERROR:root:Minimal alignment length is below aligner threshold. Pick a larger value. Quitting ...'])

    def test_set_minimal_alignment_length_above(self):
        # test that set_minimal_alignment_length fails for values above 152
        minimal_alignment_length = 152

        with self.assertLogs('root', level='INFO') as log_capture:
            self.obj.set_minimal_alignment_length(minimal_alignment_length)
            # check that log message is accurate
            self.assertEqual(log_capture.output, [
                'WARNING:root:Minimal alignment length set to above average read length of metagenomic sequencing data.'])

    def test_get_minimal_alignment_length(self):
        # test that minimal alignment length is set correctly
        self.assertEqual(self.obj.get_minimal_alignment_length(), 100)

    def test_set_threads(self):
        # test valid thread setting
        self.obj.set_threads(4)
        self.assertEqual(self.obj.get_threads(), 4)

    def test_set_too_many_threads(self):
        # test setting more threads than CPU cores
        with self.assertLogs('root', level='INFO') as log_capture:
            self.obj.set_threads(os.cpu_count() + 1)
            self.assertEqual(log_capture.output,
                             ['WARNING:root:Number of threads exceeds the total number of CPU cores.'])

    def test_set_invalid_number_of_threads(self):
        # test setting invalid number of threads
        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.set_threads(0)
            self.assertEqual(log_capture.output, ['ERROR:root:Threads have to be at least 1'])

    def test_set_sample_name(self):
        # test valid sample name
        self.obj.set_sample_name("new_sample")
        self.assertEqual(self.obj.get_sample_name(), "new_sample")

    def test_set_sample_name_empty(self):
        # test invalid (empty) sample name
        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.set_sample_name("")
            self.assertEqual(log_capture.output, ['ERROR:root:Sample name cannot be empty. Quitting'])

    def test_get_count_mode(self):
        # test get_count_mode returns expected value
        self.assertEqual(self.obj.get_count_mode(), "raw")

    def test_get_min_mgcs(self):
        # test get_min_mgcs returns expected value
        self.assertEqual(self.obj.get_min_mgcs(), 3)

    @patch("pathlib.Path.mkdir")
    def test_get_mgc_file(self, mock_mkdir):
        # call get_mgc_file
        result = self.obj.get_mgc_file()

        # assert that the mkdir method was called with correct parameters
        mock_mkdir.assert_called_once_with(exist_ok=True, parents=True)

        # assert that the method returns the correct Path object
        self.assertEqual(result, self.obj._mgc_file)

    @patch('pathlib.Path.exists', return_value=True)
    def test_set_mgc_file_exists(self, mock_exists):
        # test valid MGC file
        self.obj.set_mgc_file(pathlib.Path(mock_exists), required_to_exist=True)
        self.assertEqual(self.obj._mgc_file, pathlib.Path(mock_exists))

    def test_set_mgc_file_does_not_exist(self):
        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.set_mgc_file(pathlib.Path("/invalid/path/to/mgc_file"), required_to_exist=True)
            self.assertEqual(log_capture.output,
                             ['ERROR:root:MGC file /invalid/path/to/mgc_file does not exist. Shutting down ...'])

    @patch("pathlib.Path.mkdir")
    def test_get_inserts_file(self, mock_mkdir):
        # call get_inserts_file
        result = self.obj.get_inserts_file()

        # assert that the mkdir method was called
        mock_mkdir.assert_called_once_with(exist_ok=True, parents=True)

        # assert that the method returns the correct Path object
        self.assertEqual(result, self.obj._inserts_file)

    @patch('pathlib.Path.exists', return_value=True)
    def test_set_inserts_file_exists(self, mock_exists):
        # test valid MGC file
        self.obj.set_inserts_file(pathlib.Path(mock_exists), required_to_exist=True)
        self.assertEqual(self.obj._inserts_file, pathlib.Path(mock_exists))

    def test_set_inserts_file_not_exist(self):
        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.set_inserts_file(pathlib.Path("/invalid/path/to/inserts_file"), required_to_exist=True)
            self.assertEqual(log_capture.output, [
                'ERROR:root:Inserts file /invalid/path/to/inserts_file does not exist. Shutting down ...'])

    @patch("pathlib.Path.mkdir")
    def test_get_motu_file(self, mock_mkdir):
        # call get_motu_file
        result = self.obj.get_motu_file()

        # assert that the mkdir method was called
        mock_mkdir.assert_called_once_with(exist_ok=True, parents=True)

        # assert that the method returns the correct Path object
        self.assertEqual(result, self.obj._motu_file)

    @patch("pathlib.Path.mkdir")
    def test_get_motu_file_rel_ab(self, mock_mkdir):
        # call get_motu_file_relab
        self.obj._motu_file_relab = pathlib.Path("motu_file.relab")
        result = self.obj.get_motu_file_relab()

        # assert that the mkdir method was called
        mock_mkdir.assert_called_once_with(exist_ok=True, parents=True)

        # assert that the method returns the correct Path object
        self.assertEqual(result, self.obj._motu_file_rel_ab)

    @patch('pathlib.Path.exists', return_value=True)
    def test_set_motu_file_exists(self, mock_exists):
        # call set_motu_file
        self.obj.set_motu_file(pathlib.Path(mock_exists), required_to_exist=True)

        # assert that the paths were correctly set
        self.assertEqual(self.obj._motu_file, pathlib.Path(mock_exists))
        self.assertEqual(self.obj._motu_file_rel_ab, pathlib.Path('motu_file.relab'))

    def test_set_motu_file_does_not_exist(self):
        # do not patch the pathlib.Path function, so it fails because the motu file doesn't exist
        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.set_motu_file(pathlib.Path("/invalid/path/to/mOTU_file"), required_to_exist=True)
            self.assertEqual(log_capture.output,
                             ['ERROR:root:mOTU file /invalid/path/to/mOTU_file does not exist. Shutting down ...'])

    @patch('pathlib.Path.exists', return_value=True)
    def test_set_alignment_file_exists(self, mock_exists):
        # call set_alignment_file
        file_name = pathlib.Path("alignment_file.bam")
        self.obj.set_alignment_file(pathlib.Path(file_name), required_to_exist=True)
        self.assertEqual(self.obj._alignment_file, pathlib.Path(file_name))

    @patch('pathlib.Path.exists', return_value=True)
    def test_set_alignment_file_invalid_suffix(self, mock_exists):
        # pass alignment file with invalid suffix, expect failure
        with self.assertLogs('root', level='ERROR') as log_capture:
            invalid_file = pathlib.Path("alignment_file.txt")
            with self.assertRaises(SystemExit):
                self.obj.set_alignment_file(invalid_file, required_to_exist=True)
            self.assertEqual(log_capture.output, [
                'ERROR:root:Alignment file alignment_file.txt is/will be a BAM formatted file. Please set file suffix '
                'accordingly. Shutting down ...'])

    def test_set_alignment_file_does_not_exist(self):
        # do not patch pathlib.Path.exists, so it will fail because the file doesn't exist
        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.set_alignment_file(pathlib.Path("/invalid/path/to/alignment_file"), required_to_exist=True)
            self.assertEqual(log_capture.output, [
                'ERROR:root:Alignment file /invalid/path/to/alignment_file does not exist. Shutting down ...'])

    def test_get_read_files(self):
        # test if get_read_files returns the correct values
        expected_read_files = [
            (pathlib.Path("forward_file1"), '/1'),
            (pathlib.Path("reverse_file1"), '/2'),
            (pathlib.Path("forward_file2"), '/1'),
            (pathlib.Path("reverse_file2"), '/2'),
            (pathlib.Path("unpaired_file"), '/S')
        ]
        self.assertEqual(self.obj.get_read_files(), expected_read_files)

    def test_get_temporary_alignment_file(self):
        # test get_temporary_alignment_file
        temp_file = self.obj.get_temporary_alignment_file()
        self.assertEqual(temp_file, pathlib.Path("alignment_file_tmp.bam"))

    @patch('pathlib.Path.unlink', return_value=None)
    def test_delete_temporary_alignment_file(self, mock_unlink):
        # test delete_temporary_alignment_file
        self.obj.delete_temporary_alignment_file()
        mock_unlink.assert_called_once_with(missing_ok=True)

    @patch("pathlib.Path.mkdir")
    def test_get_alignment_file(self, mock_mkdir):
        # test get_alignment_file
        alignment_file = self.obj.get_alignment_file()
        self.assertEqual(alignment_file, pathlib.Path("alignment_file.bam"))
        mock_mkdir.assert_called_once_with(exist_ok=True, parents=True)


class TestGetFirst1000Reads(unittest.TestCase):
    def setUp(self):
        # mock object MotusParameters for testing
        self.obj = motus.MotusParameters()

    @patch("gzip.open", new_callable=mock_open, read_data=">header\nATCG\n")
    @patch("Bio.SeqIO.FastaIO.SimpleFastaParser", return_value=[("header", "ATCG")] * 1001)
    def test_get_first_1000_reads_fa_gz(self, mock_fasta_parser, mock_gzip_open):
        # test get_first_1000_reads with fa.gz file
        reads_file = pathlib.Path("/fake/path/to/file.fa.gz")
        result = self.obj.get_first_1000_reads(reads_file)

        # ensure gzip.open is called and only first 1000 reads are returned
        mock_gzip_open.assert_called_once_with(reads_file, 'rt')
        self.assertEqual(len(result), 1000)
        self.assertEqual(result[0], ("header", "ATCG"))

    @patch("builtins.open", new_callable=mock_open, read_data=">header\nATCG\n")
    @patch("Bio.SeqIO.FastaIO.SimpleFastaParser", return_value=[("header", "ATCG")] * 1001)
    def test_get_first_1000_reads_fa(self, mock_fasta_parser, mock_open_file):
        # test get_first_1000_reads with fa file
        reads_file = pathlib.Path("/fake/path/to/file.fa")
        result = self.obj.get_first_1000_reads(reads_file)

        # ensure open is called and only first 1000 reads are returned
        mock_open_file.assert_called_once_with(reads_file, 'r')
        self.assertEqual(len(result), 1000)
        self.assertEqual(result[0], ("header", "ATCG"))

    @patch("gzip.open", new_callable=mock_open, read_data=">header\nATCG\n")
    @patch("Bio.SeqIO.FastaIO.SimpleFastaParser", return_value=[("header", "ATCG")] * 1001)
    def test_get_first_1000_reads_fna_gz(self, mock_fasta_parser, mock_gzip_open):
        # test get_first_1000_reads with fna.gz file
        reads_file = pathlib.Path("/fake/path/to/file.fna.gz")
        result = self.obj.get_first_1000_reads(reads_file)

        # ensure gzip.open is called and only first 1000 reads are returned
        mock_gzip_open.assert_called_once_with(reads_file, 'rt')
        self.assertEqual(len(result), 1000)
        self.assertEqual(result[0], ("header", "ATCG"))

    @patch("builtins.open", new_callable=mock_open, read_data=">header\nATCG\n")
    @patch("Bio.SeqIO.FastaIO.SimpleFastaParser", return_value=[("header", "ATCG")] * 1001)
    def test_get_first_1000_reads_fna(self, mock_fasta_parser, mock_open_file):
        # test get_first_1000_reads with fna file
        reads_file = pathlib.Path("/fake/path/to/file.fna")
        result = self.obj.get_first_1000_reads(reads_file)

        # ensure open is called and only first 1000 reads are returned
        mock_open_file.assert_called_once_with(reads_file, 'r')
        self.assertEqual(len(result), 1000)
        self.assertEqual(result[0], ("header", "ATCG"))

    @patch("gzip.open", new_callable=mock_open, read_data=">header\nATCG\n")
    @patch("Bio.SeqIO.FastaIO.SimpleFastaParser", return_value=[("header", "ATCG")] * 1001)
    def test_get_first_1000_reads_fasta_gz(self, mock_fasta_parser, mock_gzip_open):
        # test get_first_1000_reads with fasta.gz file
        reads_file = pathlib.Path("/fake/path/to/file.fasta.gz")
        result = self.obj.get_first_1000_reads(reads_file)

        # ensure gzip.open is called and only first 1000 reads are returned
        mock_gzip_open.assert_called_once_with(reads_file, 'rt')
        self.assertEqual(len(result), 1000)
        self.assertEqual(result[0], ("header", "ATCG"))

    @patch("builtins.open", new_callable=mock_open, read_data=">header\nATCG\n")
    @patch("Bio.SeqIO.FastaIO.SimpleFastaParser", return_value=[("header", "ATCG")] * 1001)
    def test_get_first_1000_reads_fasta(self, mock_fasta_parser, mock_open_file):
        # test get_first_1000_reads with fasta file
        reads_file = pathlib.Path("/fake/path/to/file.fasta")
        result = self.obj.get_first_1000_reads(reads_file)

        # ensure open is called and only first 1000 reads are returned
        mock_open_file.assert_called_once_with(reads_file, 'r')
        self.assertEqual(len(result), 1000)
        self.assertEqual(result[0], ("header", "ATCG"))

    @patch("gzip.open", new_callable=mock_open, read_data="@header\nATCG\n+\n!!!!\n")
    @patch("Bio.SeqIO.QualityIO.FastqGeneralIterator", return_value=[("header", "ATCG", "!!!!")] * 1001)
    def test_get_first_1000_reads_fq_gz(self, mock_fastq_iterator, mock_gzip_open):
        # test get_first_1000_reads with fq.gz file
        reads_file = pathlib.Path("/fake/path/to/file.fq.gz")
        result = self.obj.get_first_1000_reads(reads_file)

        # ensure gzip.open is called and only first 1000 reads are returned
        mock_gzip_open.assert_called_once_with(reads_file, 'rt')
        self.assertEqual(len(result), 1000)
        self.assertEqual(result[0], ("header", "ATCG"))

    @patch("builtins.open", new_callable=mock_open, read_data="@header\nATCG\n+\n!!!!\n")
    @patch("Bio.SeqIO.QualityIO.FastqGeneralIterator", return_value=[("header", "ATCG", "!!!!")] * 1001)
    def test_get_first_1000_reads_fq(self, mock_fastq_iterator, mock_open_file):
        # test get_first_1000_reads with fq file
        reads_file = pathlib.Path("/fake/path/to/file.fq")
        result = self.obj.get_first_1000_reads(reads_file)

        # ensure open is called and only first 1000 reads are returned
        mock_open_file.assert_called_once_with(reads_file, 'r')
        self.assertEqual(len(result), 1000)
        self.assertEqual(result[0], ("header", "ATCG"))

    @patch("gzip.open", new_callable=mock_open, read_data="@header\nATCG\n+\n!!!!\n")
    @patch("Bio.SeqIO.QualityIO.FastqGeneralIterator", return_value=[("header", "ATCG", "!!!!")] * 1001)
    def test_get_first_1000_reads_fastq_gz(self, mock_fastq_iterator, mock_gzip_open):
        # test get_first_1000_reads with fastq.gz file
        reads_file = pathlib.Path("/fake/path/to/file.fastq.gz")
        result = self.obj.get_first_1000_reads(reads_file)

        # Ensure gzip.open is called and only first 1000 reads are returned
        mock_gzip_open.assert_called_once_with(reads_file, 'rt')
        self.assertEqual(len(result), 1000)
        self.assertEqual(result[0], ("header", "ATCG"))

    @patch("builtins.open", new_callable=mock_open, read_data="@header\nATCG\n+\n!!!!\n")
    @patch("Bio.SeqIO.QualityIO.FastqGeneralIterator", return_value=[("header", "ATCG", "!!!!")] * 1001)
    def test_get_first_1000_reads_fastq(self, mock_fastq_iterator, mock_open_file):
        # test get_first_1000_reads with fastq file
        reads_file = pathlib.Path("/fake/path/to/file.fastq")
        result = self.obj.get_first_1000_reads(reads_file)

        # ensure open is called and only first 1000 reads are returned
        mock_open_file.assert_called_once_with(reads_file, 'r')
        self.assertEqual(len(result), 1000)
        self.assertEqual(result[0], ("header", "ATCG"))

    @patch("builtins.open", new_callable=mock_open)
    def test_unknown_file_format(self, mock_open_file):
        # test failure when passing an unknown file format to get_first_1000_reads
        reads_file = pathlib.Path("/fake/path/to/file.unknown")
        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.get_first_1000_reads(reads_file)
            self.assertEqual(log_capture.output, [f'ERROR:root:Unknown file format: {reads_file}. Expecting a '
 'fasta or fastq file, can be gzipped.'])


class TestSetReadFiles(unittest.TestCase):
    def setUp(self):
        # mock object MotusParameters for testing
        self.obj = motus.MotusParameters()

    def test_no_input_files(self):
        # empty lists for forward, reverse, and unpaired files
        forward_files = []
        reverse_files = []
        unpaired_files = []

        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.set_read_files(forward_files, reverse_files, unpaired_files)
            self.assertEqual(log_capture.output, ['ERROR:root:No input files defined with -f -r or -s. Quitting ...'])

    @patch("pathlib.Path.exists", return_value=False)
    def test_files_dont_exist(self, mock_exists):
        # define some mock file paths
        forward_files = [pathlib.Path("/fake/path/to/forward_1.fq")]
        reverse_files = [pathlib.Path("/fake/path/to/reverse_2.fq")]
        unpaired_files = []

        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.set_read_files(forward_files, reverse_files, unpaired_files)
            self.assertEqual(log_capture.output, ['ERROR:root:Some read files dont exist: '
                                                  "[PosixPath('/fake/path/to/forward_1.fq'), "
                                                  "PosixPath('/fake/path/to/reverse_2.fq')]",
                                                  'ERROR:root:\t/fake/path/to/forward_1.fq',
                                                  'ERROR:root:\t/fake/path/to/reverse_2.fq'])

    @patch("pathlib.Path.exists", return_value=True)
    def test_duplicate_files(self, mock_exists):
        # duplicate files in forward and reverse lists
        forward_files = [pathlib.Path("/fake/path/to/file.fq")]
        reverse_files = [pathlib.Path("/fake/path/to/file.fq")]
        unpaired_files = []

        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.set_read_files(forward_files, reverse_files, unpaired_files)
            self.assertEqual(log_capture.output, ['ERROR:root:Duplicated read files. Please submit every file only '
                                                  'once. Shutting down ...'])

    @patch("pathlib.Path.exists", return_value=True)
    def test_unequal_forward_reverse_files(self, mock_exists):
        # unequal number of forward and reverse files
        forward_files = [pathlib.Path("/fake/path/to/forward_1.fq")]
        reverse_files = [pathlib.Path("/fake/path/to/reverse_1.fq"), pathlib.Path("/fake/path/to/reverse_2.fq")]
        unpaired_files = []

        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.set_read_files(forward_files, reverse_files, unpaired_files)
            self.assertEqual(log_capture.output, ['ERROR:root:Unequal number of files submitted with -r and -f. '
                                                  'Quitting ...'])

    @patch("pathlib.Path.exists", return_value=True)
    @patch("builtins.open")
    def test_mismatched_read_headers(self, mocked_open, mocked_exists):
        # define the content of the two files
        mock_file1 = mock_open(read_data=">header1\nATCG\n")
        mock_file2 = mock_open(read_data=">header2\nATCG\n")

        # create a side effect that returns different mock files based on the file being opened
        def mock_file_selector(file, mode='r'):
            if "forward_1.fa" in str(file):
                return mock_file1.return_value
            elif "reverse_2.fa" in str(file):
                return mock_file2.return_value

        # set the side effect for open to simulate reading different files
        mocked_open.side_effect = mock_file_selector

        # mismatched headers should cause a SystemExit
        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.set_read_files([pathlib.Path("forward_1.fa")], [pathlib.Path("reverse_2.fa")],
                                        [pathlib.Path("unpaired_file.fa")])
            self.assertIn(log_capture.output, [['ERROR:root:Headers of reads are not identical. Shutting down ...',
                                                "ERROR:root:Differing read headers: {'header1', 'header2'}",
                                                'ERROR:root:Differing read headers file 1: forward_1.fa',
                                                'ERROR:root:Differing read headers file 2: reverse_2.fa'],
                                               ['ERROR:root:Headers of reads are not identical. Shutting down ...',
                                                "ERROR:root:Differing read headers: {'header2', 'header1'}",
                                                'ERROR:root:Differing read headers file 1: forward_1.fa',
                                                'ERROR:root:Differing read headers file 2: reverse_2.fa']])

    @patch("pathlib.Path.exists", return_value=True)
    @patch("builtins.open")
    # mock content of both forward, reverse and unpaired file (will be identical)
    @patch("Bio.SeqIO.QualityIO.FastqGeneralIterator", return_value=[("header", "ATCG", "!!!!")] * 1001)
    def test_valid_files(self, mock_generator, mocked_open, mock_exists):
        # call the method with valid files
        self.obj.set_read_files([pathlib.Path("/fake/path/to/forward_1.fq")],
                                [pathlib.Path("/fake/path/to/reverse_1.fq")],
                                [pathlib.Path("/fake/path/to/unpaired_1.fq")])

        # mock file paths
        forward_files = [pathlib.Path("/fake/path/to/forward_1.fq")]
        reverse_files = [pathlib.Path("/fake/path/to/reverse_1.fq")]
        unpaired_files = [pathlib.Path("/fake/path/to/unpaired_1.fq")]

        # ensure files were set correctly in the object
        self.assertEqual(self.obj._forward_files, forward_files)
        self.assertEqual(self.obj._reverse_files, reverse_files)
        self.assertEqual(self.obj._unpaired_files, unpaired_files)


if __name__ == '__main__':
    unittest.main()
