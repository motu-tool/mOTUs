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
    @patch('motus.mentities.MOTUS_DB.load_motus_db')
    @patch('motus.motus.merge_profiles')
    def test_parse_merge_minimum_parameters(self, mock_merge, mock_load_db, mock_logging_info, mock_exists, mock_gzip_open, mock_open):
        # pass required parameters; merge_profiles and load_motus_db are mocked so no real files are needed
        sys.argv = ["motus", "merge", "-o", "output_file", "-i", "my_motus_profile1", "my_motus_profile2"]
        mock_open.return_value.__enter__.return_value = MagicMock(readline=MagicMock(side_effect=[
            "#TOOL:4.0.0_DB:4.0", "report_mode=relative_abundance", "count_mode=RAW", "min_mgcs=5"
        ]))

        with self.assertRaises(SystemExit):
            parse_merge()

    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    @patch('motus.mentities.MOTUS_DB.load_motus_db')
    def test_parse_merge_empty_motus_profile(self, mock_load_db, mock_exists, mock_gzip_open, mock_open):
        # pass a single file; parse_merge reads it as a list file, producing an empty list,
        # which causes MergedmOTUsFile to call shutdown(1) -> SystemExit
        sys.argv = ["motus", "merge", "-o", "output_file", "-i", "my_motus_profile1"]
        mock_open.return_value.__enter__.return_value = MagicMock(readline=MagicMock(side_effect=[
            "version: 4.0", "date: 2024-01-01"
        ]))

        with self.assertRaises(SystemExit):
            parse_merge()

    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_parse_merge_no_motus_profile(self, mock_exists, mock_gzip_open, mock_open):
        # pass no motus profile
        # argparse now includes the long name in error messages: -i/--input-files
        sys.argv = ["motus", "merge", "-o", "output_file", "-i"]
        mock_open.return_value.__enter__.return_value = MagicMock(readline=MagicMock(side_effect=[
            "version: 4.0", "date: 2024-01-01"
        ]))

        with self.assertRaises(SystemExit):
            sys.stderr = io.StringIO()
            parse_merge()
        self.assertIn('error: argument -i/--input-files: expected at least one argument', sys.stderr.getvalue())

    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_parse_merge_nothing_passed(self, mock_exists, mock_gzip_open, mock_open):
        # pass no parameters
        # argparse now includes the long names in error messages: -i/--input-files, -o/--output-file
        sys.argv = ["motus", "merge"]
        mock_open.return_value.__enter__.return_value = MagicMock(readline=MagicMock(side_effect=[
            "version: 4.0", "date: 2024-01-01"
        ]))

        with self.assertRaises(SystemExit):
            sys.stderr = io.StringIO()
            parse_merge()
        self.assertIn('error: the following arguments are required: -i/--input-files, -o/--output-file', sys.stderr.getvalue())

if __name__ == '__main__':
    unittest.main()
