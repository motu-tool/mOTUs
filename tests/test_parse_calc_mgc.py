import io
import unittest
from unittest.mock import patch, MagicMock
import sys
from motus.motus import parse_calc_mgc


class TestParseCalcMGC(unittest.TestCase):
    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_parse_calc_mgc_minimum_parameters(self, mock_exists, mock_gzip_open, mock_open):
        # pass required parameters
        sys.argv = ["motus", "calc_mgc", "-o", "output_file", "-i", "my_file.bam"]
        mock_open.return_value.__enter__.return_value = MagicMock(readline=MagicMock(side_effect=[
            "version: 4.0", "date: 2024-01-01"
        ]))

        with self.assertRaises(FileNotFoundError):
            parse_calc_mgc()

    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_parse_calc_mgc_all_parameters(self, mock_exists, mock_gzip_open, mock_open):
        # pass all possible parameters
        sys.argv = ["motus", "calc_mgc", "-i", "my_files.bam", "-o", "output_file", "-l", "80", "-v", "4"]

        with self.assertRaises(FileNotFoundError):
            parse_calc_mgc()

    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_parse_calc_mgc_missing_parameters(self, mock_exists, mock_gzip_open, mock_open):
        # pass parameters, some are missing
        sys.argv = ["motus", "calc_mgc", "-l", "80", "-v", "4"]

        with self.assertRaises(SystemExit):
                sys.stderr = io.StringIO()
                parse_calc_mgc()
        self.assertIn('error: the following arguments are required: -i, -o', sys.stderr.getvalue())

    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_parse_calc_mgc_no_parameters(self, mock_exists, mock_gzip_open, mock_open):
        # pass no parameters, expecting failure
        sys.argv = ["motus", "calc_mgc"]

        with self.assertRaises(SystemExit):
            sys.stderr = io.StringIO()
            parse_calc_mgc()
        self.assertIn('error: the following arguments are required: -i, -o', sys.stderr.getvalue())


if __name__ == '__main__':
    unittest.main()