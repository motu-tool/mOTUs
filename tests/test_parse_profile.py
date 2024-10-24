import io
import unittest
from unittest.mock import patch, MagicMock
import sys
from motus.motus import parse_profile


class TestParseProfile(unittest.TestCase):
    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_parse_profile_minimum_parameters(self, mock_exists, mock_gzip_open, mock_open):
        # pass required parameters
        sys.argv = ["motus", "profile", "-o", "output_file", "-f", "input_fasta_forward.fasta", "-r", "input_fasta_reverse.fasta", "-s", "input_fasta_unpaired.fasta"]
        mock_open.return_value.__enter__.return_value = MagicMock(readline=MagicMock(side_effect=[
            "version: 4.0", "date: 2024-01-01"
        ]))

        with self.assertRaises(ValueError):
            parse_profile()

    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_parse_profile_all_parameters(self, mock_exists, mock_gzip_open, mock_open):
        # pass all possible parameters
        sys.argv = ["motus", "profile", "-o", "output_file", "-f", "input_fasta_forward.fasta", "-r", "input_fasta_reverse.fasta", "-s", "input_fasta_unpaired.fasta", "-n", "SAMPLE-1", "-g", "10", "-l", "80", "-t", "4", "-v", "4", "-y", "INSERT_NORM"]
        mock_open.return_value.__enter__.return_value = MagicMock(readline=MagicMock(side_effect=[
            "version: 4.0", "date: 2024-01-01"
        ]))

        with self.assertRaises(ValueError):
            parse_profile()

    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_parse_profile_missing_parameters(self, mock_exists, mock_gzip_open, mock_open):
        # pass parameters, some are missing
        sys.argv = ["motus", "profile", "-o", "output_file", "-n", "SAMPLE-1", "-g", "10", "-l", "80", "-t", "4", "-v", "4", "-y", "INSERT_NORM"]
        mock_open.return_value.__enter__.return_value = MagicMock(readline=MagicMock(side_effect=[
            "version: 4.0", "date: 2024-01-01"
        ]))

        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                parse_profile()
            self.assertIn(
                'ERROR:root:No input files defined with -f -r or -s. Quitting ...',
                log_capture.output)

    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_parse_profile_no_parameters(self, mock_exists, mock_gzip_open, mock_open):
        # pass no parameters, expecting failure
        sys.argv = ["motus", "profile"]
        mock_open.return_value.__enter__.return_value = MagicMock(readline=MagicMock(side_effect=[
            "version: 4.0", "date: 2024-01-01"
        ]))

        with self.assertRaises(SystemExit):
            sys.stderr = io.StringIO()
            parse_profile()
        self.assertIn('error: the following arguments are required: -o', sys.stderr.getvalue())


if __name__ == '__main__':
    unittest.main()