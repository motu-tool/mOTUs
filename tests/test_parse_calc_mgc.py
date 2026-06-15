import io
import pathlib
import unittest
from unittest.mock import patch, MagicMock
import sys
from motus.motus import parse_calc_mgc

TOY_DB = pathlib.Path(__file__).parent / 'data' / 'motus4.1-toy-db'


class TestParseCalcMGC(unittest.TestCase):
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_parse_calc_mgc_minimum_parameters(self, mock_exists, mock_gzip_open):
        sys.argv = ["motus", "calc_mgc", "-o", "output_file", "-i", "my_file.bam", "-db", str(TOY_DB)]

        with self.assertRaises(FileNotFoundError):
            parse_calc_mgc()

    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_parse_calc_mgc_all_parameters(self, mock_exists, mock_gzip_open):
        # pass all possible parameters; -v (verbosity) was removed from the CLI in v4.1
        sys.argv = ["motus", "calc_mgc", "-i", "my_files.bam", "-o", "output_file", "-l", "80",
                    "-db", str(TOY_DB)]

        with self.assertRaises(FileNotFoundError):
            parse_calc_mgc()

    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_parse_calc_mgc_missing_parameters(self, mock_exists, mock_gzip_open, mock_open):
        # pass parameters, some are missing; -v (verbosity) was removed from the CLI in v4.1
        # argparse now includes long argument names in error messages (e.g. -i/--input-file)
        sys.argv = ["motus", "calc_mgc", "-l", "80"]

        with self.assertRaises(SystemExit):
                sys.stderr = io.StringIO()
                parse_calc_mgc()
        self.assertIn('error: the following arguments are required: -i/--input-file, -o/--output-file', sys.stderr.getvalue())

    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_parse_calc_mgc_no_parameters(self, mock_exists, mock_gzip_open, mock_open):
        # pass no parameters, expecting help to be shown
        sys.argv = ["motus", "calc_mgc"]

        with self.assertRaises(SystemExit):
            parse_calc_mgc()


if __name__ == '__main__':
    unittest.main()
