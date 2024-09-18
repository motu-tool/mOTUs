import unittest
from unittest.mock import patch, Mock
import collections
from motus import motus


class TestMotusFile(unittest.TestCase):
    def setup(self):
        self.my_class = motus.MotusFile()
    def test_has_no_counts(self):
        self.assertFalse(self.my_class.has_counts())
    def test_has_counts(self):
        self.my_class._samplename_2_motus_2_counts = {'sample1': {'motus1': 10}}
        self.assertTrue(self.my_class.has_counts())  # add assertion here


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
        with self.assertRaises(SystemExit):
            self.obj.merge_profiles([])

    def test_incompatible_versions(self):
        # Test with incompatible versions
        motus_files = [
            MotusFile(full_version='v1', count_mode='raw', min_mgcs=10),
            MotusFile(full_version='v2', count_mode='raw', min_mgcs=10)
        ]
        with self.assertRaises(SystemExit):
            self.obj.merge_profiles(motus_files)


    def test_incompatible_count_modes(self):
        # Test with incompatible count modes
        motus_files = [
            MotusFile(full_version='v1', count_mode='raw', min_mgcs=10),
            MotusFile(full_version='v1', count_mode='norm', min_mgcs=10)
        ]
        with self.assertRaises(SystemExit):
            self.obj.merge_profiles(motus_files)

    def test_incompatible_min_mgcs(self):
        # Test with incompatible min_mgcs
        motus_files = [
            MotusFile(full_version='v1', count_mode='raw', min_mgcs=10),
            MotusFile(full_version='v1', count_mode='raw', min_mgcs=5)
        ]
        with self.assertRaises(SystemExit):
            self.obj.merge_profiles(motus_files)

    def test_mixed_count_modes(self):
        # Test when some profiles have counts and others have relative abundances
        motus_files = [
            MotusFile(full_version='v1', count_mode='raw', min_mgcs=10, counts={"sample1": {"motu1": 10}}),
            MotusFile(full_version='v1', count_mode='raw', min_mgcs=10, relab={"sample2": {"motu1": 0.5}})
        ]
        with self.assertRaises(SystemExit):
            self.obj.merge_profiles(motus_files)

    def test_duplicate_sample_names(self):
        # Test when sample names are duplicated across profiles
        motus_files = [
            MotusFile(full_version='v1', count_mode='raw', min_mgcs=10, relab={"sample1": {"motu1": 0.5}}),
            MotusFile(full_version='v1', count_mode='raw', min_mgcs=10, relab={"sample1": {"motu2": 0.3}})
        ]
        with self.assertRaises(SystemExit):
            self.obj.merge_profiles(motus_files)

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



if __name__ == '__main__':
    unittest.main()
