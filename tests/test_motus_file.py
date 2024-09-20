import unittest
from unittest.mock import patch, mock_open
from motus import motus
import pathlib

class TestSetmOTUCounts(unittest.TestCase):

    def setUp(self):
        self.obj = motus.MotusFile()

    def test_set_mOTU_counts_int(self):
        data_int = {
            "sample1": {"motu1": 1, "motu2": 2},
            "sample2": {"motu1": 100, "motu3": 3},
        }

        self.obj.set_mOTU_counts(data_int, "int")
        # Check if counts are correctly rounded
        expected_counts = {'sample1': {'motu1': 1, 'motu2': 2}, 'sample2': {'motu1': 100, 'motu3': 3}}
        self.assertEqual(self.obj._samplename_2_motus_2_counts, expected_counts)

        # Check if relative abundances are calculated correctly
        expected_relab = {'sample1': {'motu1': 1 / 3, 'motu2': 2 / 3}, 'sample2': {'motu1': 100 / 103, 'motu3': 3 / 103}}
        self.assertEqual(self.obj._samplename_2_motus_2_relab, expected_relab)

    def test_set_mOTU_counts_float(self):
        # Test with count_type = 'float'
        data = {
            "sample1": {"motu1": 1.8, "motu2": 2.5},
            "sample2": {"motu1": 1.2, "motu3": 3.6},
        }

        self.obj.set_mOTU_counts(data, "float")

        # Check if counts remain as floats
        expected_counts = {
            "sample1": {"motu1": 1.8, "motu2": 2.5},
            "sample2": {"motu1": 1.2, "motu3": 3.6},
        }
        self.assertEqual(self.obj._samplename_2_motus_2_counts, expected_counts)

        # Check if relative abundances are calculated correctly
        expected_relab = {
            "sample1": {"motu1": 1.8 / 4.3, "motu2": 2.5 / 4.3},
            "sample2": {"motu1": 1.2 / 4.8, "motu3": 3.6 / 4.8},
        }
        self.assertEqual(self.obj._samplename_2_motus_2_relab, expected_relab)

class TestGetmOTUsFileHeader(unittest.TestCase):

    def setUp(self):
        self.obj = motus.MotusFile()
        self.obj._min_mgcs = 10
        self.obj._count_mode = 'raw'
        self.obj._full_version = 'v3.0'

    def test_valid_header_with_counts(self):
        # Test for a valid header with counts mode (relabundance=False)
        expected_header = "#v3.0\treport_mode=counts\tcount_mode=raw\tmin_mgcs=10"
        header = self.obj.get_mOTUs_file_header(relabundance=False)
        self.assertEqual(header, expected_header)

    def test_valid_header_with_relative_abundance(self):
        # Test for a valid header with relative_abundance mode (relabundance=True)
        expected_header = "#v3.0\treport_mode=relative_abundance\tcount_mode=raw\tmin_mgcs=10"
        header = self.obj.get_mOTUs_file_header(relabundance=True)
        self.assertEqual(header, expected_header)

    @patch('logging.error')
    def test_missing_min_mgcs_logs_error(self, mock_log):
        # Test when _min_mgcs is missing (None) and check logging
        self.obj._min_mgcs = None
        self.obj.get_mOTUs_file_header()
        mock_log.assert_called_with("min_mgcs parameter not set. Can't create header. Quitting...")

    @patch('logging.error')
    def test_missing_full_version_logs_error(self, mock_log):
        # Test when _full_version is missing (None) and check logging
        self.obj._full_version = None
        self.obj.get_mOTUs_file_header()
        mock_log.assert_called_with("full_version parameter not set. Can't create header. Quitting...")

    @patch('logging.error')
    def test_missing_count_mode_logs_error(self, mock_log):
        # Test when _count_mode is missing (None) and check logging
        self.obj._count_mode = None
        self.obj.get_mOTUs_file_header()
        mock_log.assert_called_with("count_mode parameter not set. Can't create header. Quitting...")

# Mocking a simple MotusFile class for testing
class MotusFile:
    def __init__(self, full_version, count_mode, min_mgcs, counts=None, relab=None):
        self._full_version = full_version
        self._count_mode = count_mode
        self._min_mgcs = min_mgcs
        self._samplename_2_motus_2_counts = counts or {}
        self._samplename_2_motus_2_relab = relab or {}

class TestMergeProfiles(unittest.TestCase):

    def setUp(self):
        self.obj = motus.MotusFile()


    def test_no_motus_files(self):
        # Test with no motus files provided
        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.merge_profiles([])
            # check that log message is accurate
            self.assertEqual(log_capture.output, ['ERROR:root:No MotusFiles found to merge. Quitting ...'])


    def test_incompatible_versions(self):
        # Test with incompatible versions
        motus_files = [
            MotusFile(full_version='v1', count_mode='raw', min_mgcs=10),
            MotusFile(full_version='v2', count_mode='raw', min_mgcs=10)
        ]
        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.merge_profiles(motus_files)
            # check that log message is accurate
            # accept both orders of versions
            string = "ERROR:root:Incompatible versions in profiles that should be merged. versions = {'v1', 'v2'} ERROR:root:Incompatible versions in profiles that should be merged. versions = {'v2', 'v1'}"
            self.assertIn(log_capture.output[0], string)



    def test_incompatible_count_modes(self):
        # Test with incompatible count modes
        motus_files = [
            MotusFile(full_version='v1', count_mode='raw', min_mgcs=10),
            MotusFile(full_version='v1', count_mode='norm', min_mgcs=10)
        ]
        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.merge_profiles(motus_files)
            # check that log message is accurate
            string = "ERROR:root:Incompatible count modes in profiles that should be merged. count modes = {'raw', 'norm'} ERROR:root:Incompatible count modes in profiles that should be merged. count modes = {'norm', 'raw'}"
            self.assertIn(log_capture.output[0], string)


    def test_incompatible_min_mgcs(self):
        # Test with incompatible min_mgcs
        motus_files = [
            MotusFile(full_version='v1', count_mode='raw', min_mgcs=10),
            MotusFile(full_version='v1', count_mode='raw', min_mgcs=5)
        ]
        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.merge_profiles(motus_files)
            # check that log message is accurate
            self.assertEqual(log_capture.output, ['ERROR:root:Incompatible min mgcs in profiles that should be merged. min mgcs ' '= {10, 5}'])


    def test_mixed_count_modes(self):
        # Test when some profiles have counts and others have relative abundances
        motus_files = [
            MotusFile(full_version='v1', count_mode='raw', min_mgcs=10, counts={"sample1": {"motu1": 10}}),
            MotusFile(full_version='v1', count_mode='raw', min_mgcs=10, relab={"sample2": {"motu1": 0.5}})
        ]
        with self.assertLogs('root', level='INFO') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.merge_profiles(motus_files)
            # check that log message is accurate
            self.assertEqual(log_capture.output, ['INFO:root:Profile files are mixed. Some are reported as counts, some as relative abundances. Quitting ...', 'INFO:root:mOTU tool shutting down with exitcode 1'])


    def test_duplicate_sample_names(self):
        # Test when sample names are duplicated across profiles
        motus_files = [
            MotusFile(full_version='v1', count_mode='raw', min_mgcs=10, relab={"sample1": {"motu1": 0.5}}),
            MotusFile(full_version='v1', count_mode='raw', min_mgcs=10, relab={"sample1": {"motu2": 0.3}})
        ]
        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.merge_profiles(motus_files)
            # check that log message is accurate
            self.assertEqual(log_capture.output, ['ERROR:root:Samplename duplicated: sample1. Quitting ...'])


    def test_successful_merge_counts(self):
        # Test successful merge of compatible profiles
        motus_files = [
            MotusFile(full_version='v1', count_mode='raw', min_mgcs=10, counts={"sample1": {"motu1": 10}}),
            MotusFile(full_version='v1', count_mode='raw', min_mgcs=10, counts={"sample2": {"motu1": 1}})
        ]
        self.obj.merge_profiles(motus_files)

        # Check that the data was merged correctly
        self.assertEqual(self.obj._motus_with_abundance, ['motu1'])
        self.assertEqual(self.obj._samplename_2_motus_2_counts['sample1']['motu1'], 10)
        self.assertEqual(self.obj._samplename_2_motus_2_counts['sample2']['motu1'], 1)
        self.assertEqual(self.obj._full_version, 'v1')
        self.assertEqual(self.obj._count_mode, 'raw')
        self.assertEqual(self.obj._min_mgcs, 10)

    def test_successful_merge_relab(self):
        # Test successful merge of relative abundance profiles
        motus_files = [
            MotusFile(full_version='v1', count_mode='raw', min_mgcs=10, relab={"sample1": {"motu1": 0.5}}),
            MotusFile(full_version='v1', count_mode='raw', min_mgcs=10, relab={"sample2": {"motu1": 2}})
        ]
        self.obj.merge_profiles(motus_files)

        # Check that the data was merged correctly
        self.assertEqual(self.obj._motus_with_abundance, ['motu1'])
        self.assertEqual(self.obj._samplename_2_motus_2_relab['sample1']['motu1'], 0.5)
        self.assertEqual(self.obj._samplename_2_motus_2_relab['sample2']['motu1'], 2)
        self.assertEqual(self.obj._full_version, 'v1')
        self.assertEqual(self.obj._count_mode, 'raw')
        self.assertEqual(self.obj._min_mgcs, 10)


class TestReadMotusFile(unittest.TestCase):

    def setUp(self):
        # Initialize the object that contains the read_mOTUs_file method
        self.obj = motus.MotusFile()

    @patch('builtins.open', new_callable=mock_open,
           read_data='#TOOL1.0\treport_mode=counts\tcount_mode=raw\tmin_mgcs=10\nMOTU\tsample1'
                     '\nmOTU1\t74')
    def test_valid_mOTUs_file_4_fields(self, mock_file):
        # Test a valid mOTUs file with 4 fields in the header
        self.obj.read_mOTUs_file(mock_file)

        self.assertEqual(self.obj._count_mode, 'raw')
        self.assertEqual(self.obj._min_mgcs, 10)
        self.assertEqual(self.obj._full_version, 'TOOL1.0')

# TODO write test for 7 columns
    """    @patch('builtins.open', new_callable=mock_open,
           read_data='#TOOL1.0\treport_mode=counts\tcount_mode=raw\tmin_mgcs=10\ttaxonomy=species\taggregated=0\tlevel=species')
    def test_valid_mOTUs_file_7_fields(self, mock_file):
        # Test a valid mOTUs file with 7 fields in the header
        self.obj.read_mOTUs_file(mock_file)

        self.assertEqual(self.obj._count_mode, 'raw')
        self.assertEqual(self.obj._min_mgcs, 10)
        self.assertEqual(self.obj._taxonomy, 'species')
        self.assertEqual(self.obj._taxonomy_level, 'species')
        self.assertEqual(self.obj._full_version, 'TOOL1.0')
        self.assertEqual(self.obj._samplename_2_motus_2_counts['sample1']['motu1'], 10)"""

    @patch('builtins.open', new_callable=mock_open,
           read_data='TOOL1.0\tcounts\tmin_mgcs=10\nMOTU\tsample1\tsample2\nmotu1\t10\t20\n')
    def test_malformed_header(self, mock_file):
        # Test a malformed header (too few fields)
        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.read_mOTUs_file(mock_file)
            # check that log message is accurate
            self.assertEqual(log_capture.output, ['ERROR:root:The header of this mOTUs file looks malformed. Please check. '  'Quitting ...', 'ERROR:root:TOOL1.0\tcounts\tmin_mgcs=10'])


    @patch('builtins.open', new_callable=mock_open,
           read_data='#TOOL1.0\treport_mode=counts\tmin_mgcs=10\nMOTU\tsample1'
                     '\nmOTU1\t74')
    def test_missing_count_mode(self, mock_file):
        # Test missing min_mgcs field in the header
        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.read_mOTUs_file(mock_file)
            # check that log message is accurate
            self.assertEqual(log_capture.output, ['ERROR:root:The header of this mOTUs file looks malformed. Please check. ' 'Quitting ...', 'ERROR:root:#TOOL1.0\treport_mode=counts\tmin_mgcs=10'])

    @patch('builtins.open', new_callable=mock_open,
           read_data='#TOOL1.0\treport_mode=counts\tcount_mode=raw\nMOTU\tsample1\tsample2\nmotu1\t10\t20\n')
    def test_missing_min_mgcs(self, mock_file):
        # Test missing min_mgcs field in the header
        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.read_mOTUs_file(mock_file)
            # check that log message is accurate
            self.assertEqual(log_capture.output, ['ERROR:root:The header of this mOTUs file looks malformed. Please check. ' 'Quitting ...', 'ERROR:root:#TOOL1.0\treport_mode=counts\tcount_mode=raw'])

    @patch('builtins.open', new_callable=mock_open,
           read_data='#TOOL1.0\treport_mode=counts\tcount_mode=raw\tmin_mgcs=10\taggregated=False\tlevel=species')
    def test_missing_taxonomy(self, mock_file):
        # Test missing min_mgcs field in the header
        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.read_mOTUs_file(mock_file)
            # check that log message is accurate
            self.assertEqual(log_capture.output, ['ERROR:root:The header of this mOTUs file looks malformed. Please check. ' 'Quitting ...',  'ERROR:root:#TOOL1.0\treport_mode=counts\tcount_mode=raw\tmin_mgcs=10\taggregated=False\tlevel=species'])

    @patch('builtins.open', new_callable=mock_open,
           read_data='#TOOL1.0\treport_mode=counts\tcount_mode=raw\tmin_mgcs=10\ttaxonomy=species\taggregated=True\t\nMOTU\tsample1\tsample2\nmotu1\t10\t20\n')
    def test_missing_level(self, mock_file):
        # Test an invalid report mode
        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.read_mOTUs_file(mock_file)
            # check that log message is accurate
            self.assertEqual(log_capture.output, ['ERROR:root:The header of this mOTUs file looks malformed. Please check. Quitting ...', 'ERROR:root:#TOOL1.0\treport_mode=counts\tcount_mode=raw\tmin_mgcs=10\ttaxonomy=species\taggregated=True'])


    @patch('builtins.open', new_callable=mock_open,
           read_data='#TOOL1.0\treport_mode=invalid_mode\tcount_mode=raw\tmin_mgcs=10\nMOTU\tsample1\tsample2\nmotu1\t10\t20\n')
    def test_invalid_report_mode(self, mock_file):
        # Test an invalid report mode
        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj.read_mOTUs_file(mock_file)
            # check that log message is accurate
            self.assertEqual(log_capture.output, ['ERROR:root:Report mode can only be counts or relative_abundance but is invalid_mode. Quitting ...'])


    @patch('builtins.open', new_callable=mock_open,
           read_data='#TOOL1.0\treport_mode=counts\tcount_mode=raw\tmin_mgcs=10\ttaxonomy=species\taggregated=True\tlevel=species\nMOTU\tsample1\tsample2\nmotu1\t10\t20\n')
    def test_aggregated_table(self, mock_file):
        # Test an aggregated table that leads to a shutdown
        with self.assertLogs('root', level='ERROR') as log_capture:
            mock_file = pathlib.Path("test_motus_file.tsv")
            with self.assertRaises(SystemExit):
                self.obj.read_mOTUs_file(mock_file)
            # check that log message is accurate
            self.assertEqual(log_capture.output, ['ERROR:root:The mOTUs profile (test_motus_file.tsv) has values aggregated at non-mOTU level (species). This table is an endproduct and cannot be used in mOTUs anymore. Quitting ...'])


class TestWriteMOTUsFile(unittest.TestCase):

    def setUp(self):
        # Mocking the internal attributes of the object
        self.obj = motus.MotusFile()
        self.obj._samplename_2_motus_2_counts = {
            'sample1': {'motu1': 10, 'motu2': 5},
            'sample2': {'motu1': 7, 'motu2': 3}
        }
        self.obj._samplename_2_motus_2_relab = {
            'sample1': {'motu1': 0.5, 'motu2': 0.25},
            'sample2': {'motu1': 0.35, 'motu2': 0.15}
        }
        self.obj._motus_with_abundance = ['motu1', 'motu2']
        self.obj._count_mode = "RAW"
        self.obj._min_mgcs=5
        self.obj._full_version = "TOOL:4.0.0_DB:4.0"


    @patch('builtins.open', new_callable=mock_open)
    def test_write_mOTUs_file_relabundance_false(self, mock_file):
        # Test with relabundance=False (normal counts)
        self.obj.write_mOTUs_file('output.tsv', relabundance=False)

        # Check the headers written to the file
        mock_file().write.assert_any_call("#TOOL:4.0.0_DB:4.0\treport_mode=counts\tcount_mode=RAW\tmin_mgcs=5\n")
        mock_file().write.assert_any_call("MOTU\tsample1\tsample2\n")

        # Check the data written for motu1 and motu2
        mock_file().write.assert_any_call("motu1\t10\t7\n")
        mock_file().write.assert_any_call("motu2\t5\t3\n")

    @patch('builtins.open', new_callable=mock_open)
    def test_write_mOTUs_file_relabundance_true(self, mock_file):
        # Test with relabundance=True (relative abundance)
        self.obj.write_mOTUs_file('output.tsv', relabundance=True)

        # Check the contents written to the file
        mock_file().write.assert_any_call("#TOOL:4.0.0_DB:4.0\treport_mode=relative_abundance\tcount_mode=RAW\tmin_mgcs=5\n")
        mock_file().write.assert_any_call("MOTU\tsample1\tsample2\n")  # The sample header

        # Check the data written for motu1 and motu2 with 8 decimal places
        mock_file().write.assert_any_call("motu1\t0.50000000\t0.35000000\n")
        mock_file().write.assert_any_call("motu2\t0.25000000\t0.15000000\n")

    @patch('builtins.open', new_callable=mock_open)
    def test_write_mOTUs_file_with_normalized_counts(self, mock_file):
        # Test with normalized counts, checking that integers are still written if count mode doesn't contain "NORM"
        self.obj._count_mode = "NORM"
        self.obj.write_mOTUs_file('output.tsv', relabundance=False)

        # Check the data written for motu1 and motu2 (should remain as integers)
        mock_file().write.assert_any_call("motu1\t10\t7\n")
        mock_file().write.assert_any_call("motu2\t5\t3\n")


if __name__ == '__main__':
    unittest.main()
