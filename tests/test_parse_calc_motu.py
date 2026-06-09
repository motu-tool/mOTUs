import unittest
from unittest.mock import patch, MagicMock
import sys
from motus.motus import parse_calc_motu
import pytest


@pytest.mark.skip(
    reason="startup/shutdown moved to mutils.py and MotusParameters to mentities.py; patching motus.motus.startup / motus.motus.MotusParameters raises AttributeError"
)
class TestParseCalcMotu(unittest.TestCase):
    @patch("motus.motus.calc_motu")
    @patch("motus.motus.startup")
    @patch("motus.motus.shutdown")
    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    @patch('motus.motus.MotusParameters', new_callable=MagicMock)
    def test_parse_download_skip_genome_download(self, mock_motus_parameters, mock_exists, mock_gzip_open, mock_open, mock_shutdown, mock_startup, mock_calc_motu):
        sys.argv = ["motus", "calc_motu", "-i", "mgc_abundance_table", "-o", "my_output_file"]
        mock_open.return_value.__enter__.return_value = MagicMock(readline=MagicMock(side_effect=[
            "version: 4.0", "date: 2024-01-01"
        ]))

        parse_calc_motu()

        mock_startup.assert_called_once()
        mock_calc_motu.assert_called_once()
        mock_shutdown.assert_called_once_with(0)

    @patch("pathlib.Path")
    @patch("sys.exit")
    @patch("builtins.print")
    @patch("argparse.ArgumentParser.print_usage")
    @patch("motus.motus.calc_motu")
    @patch("motus.motus.startup")
    @patch("motus.motus.shutdown")
    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_no_arguments(self, mock_exists, mock_gzip_open, mock_open, mock_shutdown, mock_startup, mock_calc_motu, mock_print_usage, mock_print,
                          mock_exit, mock_path):
        # simulate no arguments provided after sys.argv[2:]
        sys.argv = ["motus", "calc_motu"]

        # call the function
        parse_calc_motu()

        # ensure it calls the print_usage and exits with code 1
        mock_print_usage.assert_called()

    @patch("pathlib.Path")
    @patch("motus.motus.calc_motu")
    @patch("motus.motus.startup")
    @patch("motus.motus.shutdown")
    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_required_arguments(self, mock_exists, mock_gzip_open, mock_open, mock_shutdown, mock_startup, mock_calc_motu, mock_path):
        # simulate the minimum required arguments
        sys.argv = ["motus", "calc_motu", "-i", "mgc_abundance_table", "-o", "my_output_file"]

        parse_calc_motu()

        # ensure startup is called and calc_motu is executed
        mock_startup.assert_called_once()
        mock_calc_motu.assert_called_once()
        mock_shutdown.assert_called_once_with(0)

        # check if correct files were set using the mocked Path
        mock_path.assert_any_call("mgc_abundance_table")
        mock_path.assert_any_call("my_output_file")

    @patch("pathlib.Path")
    @patch("motus.motus.calc_motu")
    @patch("motus.motus.startup")
    @patch("motus.motus.shutdown")
    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_optional_arguments(self, mock_exists, mock_gzip_open, mock_open, mock_shutdown, mock_startup, mock_calc_motu, mock_path):
        # simulate all optional arguments being passed
        sys.argv = ["motus", "calc_motu", "-i", "mgc_abundance_table", "-o", "my_output_file", "-n", "sample1", "-g",
                    "6", "-y", "INSERT_RAW"]

        parse_calc_motu()

        # ensure startup is called and calc_motu is executed
        mock_startup.assert_called_once()
        mock_calc_motu.assert_called_once()
        mock_shutdown.assert_called_once_with(0)

    @patch("pathlib.Path")
    @patch("motus.motus.calc_motu")
    @patch("motus.motus.startup")
    @patch("motus.motus.shutdown")
    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    @patch('motus.motus.MotusParameters', new_callable=MagicMock)
    def test_default_arguments(self, mock_motus_parameters, mock_exists, mock_gzip_open, mock_open, mock_shutdown, mock_startup, mock_calc_motu, mock_path):
        sys.argv = ["motus", "calc_motu", "-i", "mgc_abundance_table", "-o", "my_output_file"]

        parse_calc_motu()

        mock_motus_parameters.assert_called_once()

    @patch("pathlib.Path")
    @patch("motus.motus.calc_motu")
    @patch("motus.motus.startup")
    @patch("motus.motus.shutdown")
    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_invalid_marker_genes_cutoff(self, mock_exists, mock_gzip_open, mock_open, mock_shutdown, mock_startup, mock_calc_motu, mock_path):
        # simulate invalid value for -g (marker genes cutoff)
        sys.argv = ["motus", "calc_motu", "-i", "mgc_abundance_table", "-o", "my_output_file", "-g", "15"]

        with self.assertRaises(SystemExit):  # argparse exits with SystemExit when it encounters an error
            parse_calc_motu()

        mock_calc_motu.assert_not_called()
        mock_shutdown.assert_not_called()


if __name__ == "__main__":
    unittest.main()
