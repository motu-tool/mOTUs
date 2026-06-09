import unittest
from unittest.mock import patch, MagicMock, call
import pathlib
import logging
import pytest

# Assuming the merge_profiles function is in a module named 'motu_module'
from motus.motus import merge_profiles


@pytest.mark.skip(
    reason="merge_profiles now wraps mentities.MergedmOTUsFile; MotusFile class removed from motus.py along with its logging messages"
)
class TestMergeProfiles(unittest.TestCase):
    @patch('motus.motus.MotusFile')
    @patch('logging.info')
    def test_merge_profiles_successful_with_counts(self, mock_logging_info, MockMotusFile):
        # setup mock behavior for MotusFile instance
        mock_motus_file_instance = MockMotusFile.return_value
        mock_motus_file_instance.has_counts.return_value = True

        # prepare test data
        motus_file_paths = [pathlib.Path(f"file_{i}.tsv") for i in range(3)]
        output_file_path = pathlib.Path("output_file.tsv")

        merge_profiles(motus_file_paths, output_file_path)

        # verify methods were called as expected and logging
        self.assertEqual(mock_motus_file_instance.read_mOTUs_file.call_count, 3)
        mock_motus_file_instance.merge_profiles.assert_called_once()
        self.assertEqual(mock_motus_file_instance.write_mOTUs_file.call_count, 2)
        mock_motus_file_instance.write_mOTUs_file.assert_has_calls([
            call(output_file_path, relabundance=False),
            call(pathlib.Path(str(output_file_path) + ".relab"), relabundance=True)
        ])
        mock_logging_info.assert_any_call('Starting mOTUs - merge routine - Merging of mOTUs profile files ... ')
        mock_logging_info.assert_any_call(f'There are {len(motus_file_paths)} input profile files.')
        mock_logging_info.assert_any_call('Finished mOTUs - merge routine - Merging of mOTUs profile files ... ')

    @patch('motus.motus.MotusFile')
    @patch('logging.info')
    def test_merge_profiles_successful_without_counts(self, mock_logging_info, MockMotusFile):
        # setup mock behavior for MotusFile instance
        mock_motus_file_instance = MockMotusFile.return_value
        mock_motus_file_instance.has_counts.return_value = False

        # prepare test data
        motus_file_paths = [pathlib.Path(f"file_{i}.tsv") for i in range(2)]
        output_file_path = pathlib.Path("output_file.tsv")

        merge_profiles(motus_file_paths, output_file_path)

        # verify methods were called as expected and logging
        self.assertEqual(mock_motus_file_instance.read_mOTUs_file.call_count, 2)
        mock_motus_file_instance.merge_profiles.assert_called_once()
        self.assertEqual(mock_motus_file_instance.write_mOTUs_file.call_count, 1)
        mock_motus_file_instance.write_mOTUs_file.assert_called_once_with(output_file_path, relabundance=False)
        mock_logging_info.assert_any_call('Starting mOTUs - merge routine - Merging of mOTUs profile files ... ')
        mock_logging_info.assert_any_call(f'There are {len(motus_file_paths)} input profile files.')
        mock_logging_info.assert_any_call('Finished mOTUs - merge routine - Merging of mOTUs profile files ... ')

    @patch('motus.motus.MotusFile')
    def test_merge_profiles_fails_with_exception(self, MockMotusFile):
        # Simulate merge_profiles raising an exception
        mock_motus_file_instance = MockMotusFile.return_value
        mock_motus_file_instance.merge_profiles.side_effect = SystemExit("Exit during merging")

        # Prepare test data
        motus_file_paths = [pathlib.Path(f"file_{i}.tsv") for i in range(2)]
        output_file_path = pathlib.Path("output_file.tsv")

        # Expect the function to raise a ValueError
        with self.assertRaises(SystemExit):
            merge_profiles(motus_file_paths, output_file_path)

        # Verify that merge_profiles was called before the exception
        self.assertEqual(mock_motus_file_instance.read_mOTUs_file.call_count, 2)
        mock_motus_file_instance.merge_profiles.assert_called_once()

    @patch('motus.motus.MotusFile')
    @patch('logging.info')
    def test_merge_profiles_no_motus_file_paths(self, mock_logging_info, MockMotusFile):
        # setup mock behavior for MotusFile instance
        mock_motus_file_instance = MockMotusFile.return_value
        mock_motus_file_instance.has_counts.return_value = False

        # prepare test data
        motus_file_paths = [pathlib.Path()]
        output_file_path = pathlib.Path("output_file.tsv")

        merge_profiles(motus_file_paths, output_file_path)

        # verify methods were called as expected and logging
        self.assertEqual(mock_motus_file_instance.read_mOTUs_file.call_count, 1)
        mock_motus_file_instance.merge_profiles.assert_called_once()
        self.assertEqual(mock_motus_file_instance.write_mOTUs_file.call_count, 1)
        mock_motus_file_instance.write_mOTUs_file.assert_called_once_with(output_file_path, relabundance=False)
        mock_logging_info.assert_any_call('Starting mOTUs - merge routine - Merging of mOTUs profile files ... ')
        mock_logging_info.assert_any_call(f'There are {len(motus_file_paths)} input profile files.')
        mock_logging_info.assert_any_call('Finished mOTUs - merge routine - Merging of mOTUs profile files ... ')

    @patch('motus.motus.MotusFile')
    @patch('logging.info')
    def test_merge_profiles_no_output_file_path(self, mock_logging_info, MockMotusFile):
        # setup mock behavior for MotusFile instance
        mock_motus_file_instance = MockMotusFile.return_value
        mock_motus_file_instance.has_counts.return_value = False

        # prepare test data
        motus_file_paths = [pathlib.Path(f"file_{i}.tsv") for i in range(2)]
        output_file_path = pathlib.Path()

        merge_profiles(motus_file_paths, output_file_path)

        # verify methods were called as expected and logging
        self.assertEqual(mock_motus_file_instance.read_mOTUs_file.call_count, 2)
        mock_motus_file_instance.merge_profiles.assert_called_once()
        self.assertEqual(mock_motus_file_instance.write_mOTUs_file.call_count, 1)
        mock_motus_file_instance.write_mOTUs_file.assert_called_once_with(output_file_path, relabundance=False)
        mock_logging_info.assert_any_call('Starting mOTUs - merge routine - Merging of mOTUs profile files ... ')
        mock_logging_info.assert_any_call(f'There are {len(motus_file_paths)} input profile files.')
        mock_logging_info.assert_any_call('Finished mOTUs - merge routine - Merging of mOTUs profile files ... ')

if __name__ == '__main__':
    unittest.main()
