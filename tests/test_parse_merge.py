import builtins
import io
import pathlib
import unittest
from unittest.mock import patch, MagicMock
import sys
from motus.motus import parse_merge

TOY_DB = pathlib.Path(__file__).parent / 'data' / 'motus4.1-toy-db'

_real_open = builtins.open


class TestParseMerge(unittest.TestCase):
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    @patch('motus.motus.merge_profiles')
    def test_parse_merge_minimum_parameters(self, mock_merge, mock_exists, mock_gzip_open):
        # pass required parameters; merge_profiles is mocked so no real profile files are needed.
        # With 2 files parse_merge skips the list-file branch so no builtins.open mock is needed.
        sys.argv = ["motus", "merge", "-o", "output_file", "-i", "my_motus_profile1", "my_motus_profile2",
                    "-db", str(TOY_DB)]

        with self.assertRaises(SystemExit):
            parse_merge()

    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_parse_merge_empty_motus_profile(self, mock_exists, mock_gzip_open, mock_open):
        # pass a single file; parse_merge reads it as a list file, producing an empty list,
        # which causes MergedmOTUsFile to call shutdown(1) -> SystemExit
        sys.argv = ["motus", "merge", "-o", "output_file", "-i", "my_motus_profile1",
                    "-db", str(TOY_DB)]

        # Fake the list file as empty; let the toy DB version file go through real open.
        def selective_open(file, mode='r', *args, **kwargs):
            if 'my_motus_profile1' in str(file):
                m = MagicMock()
                m.__enter__ = lambda s: s
                m.__exit__ = MagicMock(return_value=False)
                m.__iter__ = lambda s: iter([])
                return m
            return _real_open(file, mode, *args, **kwargs)

        mock_open.side_effect = selective_open

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
        # pass no parameters, expecting help to be shown
        sys.argv = ["motus", "merge"]

        with self.assertRaises(SystemExit):
            parse_merge()

if __name__ == '__main__':
    unittest.main()
