import unittest
from motus import motus
from unittest.mock import patch, mock_open
import pathlib
import os

class TestDBModeClass(unittest.TestCase):

    def setUp(self):
        self.obj = motus.MotusParameters()

        # Set initial values to ensure a known state
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
        # Test strict mode is initially True
        self.assertTrue(self.obj.is_strict_db_mode())

    def test_enable_lenient_mode(self):
        # Test that enable_lenient_mode sets _is_strict_db_mode to False
        self.obj.enable_lenient_mode()
        self.assertFalse(self.obj.is_strict_db_mode())

    def test_set_minimal_number_of_mgcs(self):
        # Test that set_minimal_number_of_mgcs sets the correct value
        self.obj.set_minimal_number_of_mgcs(10)
        self.assertEqual(self.obj._min_mgcs, 10)

        self.obj.set_minimal_number_of_mgcs(5)
        self.assertEqual(self.obj._min_mgcs, 5)

    def test_set_count_mode(self):
        # Test that set_count_mode sets the correct count mode
        self.obj.set_count_mode('raw')
        self.assertEqual(self.obj._count_mode, 'raw')

        self.obj.set_count_mode('normalized')
        self.assertEqual(self.obj._count_mode, 'normalized')

    def test_get_count_type(self):
        # Test that set_count_mode sets the correct count mode
        self.obj._count_mode = "NORM"

        self.assertEqual(self.obj.get_count_type(), 'float')

    def test_get_count_type(self):
        # Test that set_count_mode sets the correct count mode
        self.obj._count_mode = "OTHER"

        self.assertEqual(self.obj.get_count_type(), 'int')

    def test_set_minimal_alignment_length(self):
        # Test that set_count_mode sets the correct count mode
        minimal_alignment_length = 50
        self.obj.set_minimal_alignment_length(minimal_alignment_length)

        self.assertEqual(self.obj._min_alignment_length, 50)

    def test_set_minimal_alignment_length_below(self):
        # Test that set_count_mode sets the correct count mode
        minimal_alignment_length = 29

        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.set_minimal_alignment_length(minimal_alignment_length)
            # check that log message is accurate
            self.assertEqual(log_capture.output, ['ERROR:root:Minimal alignment length is below aligner threshold. Pick a larger value. Quitting ...'])

    def test_set_minimal_alignment_length_above(self):
        # Test that set_count_mode sets the correct count mode
        minimal_alignment_length = 152

        with self.assertLogs('root', level='INFO') as log_capture:
            self.obj.set_minimal_alignment_length(minimal_alignment_length)
            # check that log message is accurate
            self.assertEqual(log_capture.output, ['WARNING:root:Minimal alignment length set to above average read length of metagenomic sequencing data.'])

    def test_get_minimal_alignment_length(self):
        self.assertEqual(self.obj.get_minimal_alignment_length(), 100)

    def test_set_threads(self):
        # Test valid thread setting
        self.obj.set_threads(4)
        self.assertEqual(self.obj.get_threads(), 4)

        # Test setting more threads than CPU cores
        with self.assertLogs('root', level='INFO') as log_capture:
            #with self.assertRaises(Warning):
            self.obj.set_threads(os.cpu_count() + 1)
            self.assertEqual(log_capture.output, ['WARNING:root:Number of threads exceeds the total number of CPU cores.'])

        # Test setting invalid number of threads
        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.set_threads(0)
            self.assertEqual(log_capture.output, ['ERROR:root:Threads have to be at least 1'])


    def test_set_sample_name(self):
        # Test valid sample name
        self.obj.set_sample_name("new_sample")
        self.assertEqual(self.obj.get_sample_name(), "new_sample")

        # Test invalid (empty) sample name
        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.set_sample_name("")
            self.assertEqual(log_capture.output, ['ERROR:root:Sample name cannot be empty. Quitting'])


    def test_get_count_mode(self):
        self.assertEqual(self.obj.get_count_mode(), "raw")

    def test_get_min_mgcs(self):
        self.assertEqual(self.obj.get_min_mgcs(), 3)

    def test_set_mgc_file_does_not_exist(self):
        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.set_mgc_file(pathlib.Path("/invalid/path/to/mgc_file"), required_to_exist=True)
            self.assertEqual(log_capture.output, ['ERROR:root:MGC file /invalid/path/to/mgc_file does not exist. Shutting down ...'])


    @patch('pathlib.Path.exists', return_value=True)
    def test_set_mgc_file_exists(self, mock_exists):
        # Test valid MGC file
        self.obj.set_mgc_file(pathlib.Path("/valid/path/to/mgc_file"), required_to_exist=True)
        self.assertEqual(self.obj._mgc_file, pathlib.Path("/valid/path/to/mgc_file"))

    def test_set_alignment_file_invalid_suffix(self):
        with self.assertLogs('root', level='ERROR') as log_capture:
            invalid_file = pathlib.Path("alignment_file.txt")
            with self.assertRaises(SystemExit):
                self.obj.set_alignment_file(invalid_file, required_to_exist=True)
            self.assertEqual(log_capture.output, ['ERROR:root:Alignment file alignment_file.txt is/will be a BAM formatted file. Please set file suffix accordingly. Shutting down ...'])

    def test_get_read_files(self):
        expected_read_files = [
            (pathlib.Path("forward_file1"), '/1'),
            (pathlib.Path("reverse_file1"), '/2'),
            (pathlib.Path("forward_file2"), '/1'),
            (pathlib.Path("reverse_file2"), '/2'),
            (pathlib.Path("unpaired_file"), '/S')
        ]
        self.assertEqual(self.obj.get_read_files(), expected_read_files)

    def test_get_temporary_alignment_file(self):
        temp_file = self.obj.get_temporary_alignment_file()
        self.assertEqual(temp_file, pathlib.Path("alignment_file_tmp.bam"))

    @patch('pathlib.Path.unlink', return_value=None)
    def test_delete_temporary_alignment_file(self, mock_unlink):
        self.obj.delete_temporary_alignment_file()
        mock_unlink.assert_called_once_with(missing_ok=True)

    def test_get_alignment_file(self):
        alignment_file = self.obj.get_alignment_file()
        self.assertEqual(alignment_file, pathlib.Path("alignment_file.bam"))


if __name__ == '__main__':
    unittest.main()

