import io
import unittest
from unittest.mock import patch, MagicMock
import sys
from motus.motus import parse_merge


class TestParseMerge(unittest.TestCase):
    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    @patch('logging.info')
    def test_parse_merge_minimum_parameters(self, mock_logging_info, mock_exists, mock_gzip_open, mock_open):
        # pass required parameters
        sys.argv = ["motus", "merge", "-o", "output_file", "-i", "my_motus_profile1", "my_motus_profile2"]
        mock_open.return_value.__enter__.return_value = MagicMock(readline=MagicMock(side_effect=[
            "#TOOL:4.0.0_DB:4.0", "report_mode=relative_abundance", "count_mode=RAW", "min_mgcs=5"
        ]))

        with self.assertRaises(SystemExit):
            parse_merge()
            mock_logging_info.assert_any_call('Starting mOTUs - merge routine - Merging of mOTUs profile files ... ')
            mock_logging_info.assert_any_call('There are 2 input profile files. ')

    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    @patch('logging.info')
    def test_parse_merge_empty_motus_profile(self, mock_logging_info, mock_exists, mock_gzip_open, mock_open):
        # pass empty motus profile
        sys.argv = ["motus", "merge", "-o", "output_file", "-i", "my_motus_profile1"]
        mock_open.return_value.__enter__.return_value = MagicMock(readline=MagicMock(side_effect=[
            "version: 4.0", "date: 2024-01-01"
        ]))

        with self.assertRaises(SystemExit):
            parse_merge()
            mock_logging_info.assert_any_call('Starting mOTUs - merge routine - Merging of mOTUs profile files ... ')
            mock_logging_info.assert_any_call('There are 2 input profile files. ')

    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    @patch('logging.info')
    def test_parse_merge_no_motus_profile(self, mock_logging_info, mock_exists, mock_gzip_open, mock_open):
        # pass no motus profile
        sys.argv = ["motus", "merge", "-o", "output_file", "-i"]
        mock_open.return_value.__enter__.return_value = MagicMock(readline=MagicMock(side_effect=[
            "version: 4.0", "date: 2024-01-01"
        ]))

        with self.assertRaises(SystemExit):
            sys.stderr = io.StringIO()
            parse_merge()
            mock_logging_info.assert_any_call('Starting mOTUs - merge routine - Merging of mOTUs profile files ... ')
            mock_logging_info.assert_any_call('There are 2 input profile files. ')
        self.assertIn('error: argument -i: expected at least one argument', sys.stderr.getvalue())

    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    @patch('logging.info')
    def test_parse_merge_nothing_passed(self, mock_logging_info, mock_exists, mock_gzip_open, mock_open):
        # pass no parameters
        sys.argv = ["motus", "merge"]
        mock_open.return_value.__enter__.return_value = MagicMock(readline=MagicMock(side_effect=[
            "version: 4.0", "date: 2024-01-01"
        ]))

        with self.assertRaises(SystemExit):
            sys.stderr = io.StringIO()
            parse_merge()
            mock_logging_info.assert_any_call('Starting mOTUs - merge routine - Merging of mOTUs profile files ... ')
            mock_logging_info.assert_any_call('There are 2 input profile files. ')
        self.assertIn('error: the following arguments are required: -i, -o', sys.stderr.getvalue())

if __name__ == '__main__':
    unittest.main()