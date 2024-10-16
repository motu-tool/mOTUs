import io
import unittest
from unittest.mock import patch, MagicMock
import sys
from motus.motus import parse_download


class TestParseDownload(unittest.TestCase):
    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    @patch('logging.info')
    def test_parse_download_skip_genome_download(self, mock_logging_info, mock_exists, mock_gzip_open, mock_open):
        key_word = "Enterococcus"
        sys.argv = ["motus", "download", "-s", "metadata_file", "-w", key_word, "-l"]

        with self.assertRaises(SystemExit):
            parse_download()

        mock_logging_info.assert_any_call(f"Searching for keyword: {key_word}.")
        mock_logging_info.assert_any_call("Finished writing genome information to metadata_file")

    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    @patch('logging.info')
    def test_parse_download_skip_genome_download_pass_folder(self, mock_logging_info, mock_exists, mock_gzip_open, mock_open):
        key_word = "Enterococcus"
        output_folder = "my_output_folder"
        sys.argv = ["motus", "download", "-s", "metadata_file", "-w", key_word, "-l", "-o", output_folder]

        with self.assertRaises(SystemExit):
            parse_download()

        mock_logging_info.assert_any_call(f"Searching for keyword: {key_word}.")
        mock_logging_info.assert_any_call("Finished writing genome information to metadata_file")

    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    @patch('logging.info')
    @patch('os.mkdir')
    def test_parse_download_skip_with_genome_download(self, mock_mkdir, mock_logging_info, mock_exists, mock_gzip_open, mock_open):
        key_word = "Enterococcus"
        output_folder = "my_output_folder"
        sys.argv = ["motus", "download", "-s", "metadata_file", "-w", key_word, "-o", output_folder]

        with self.assertRaises(SystemExit):
            parse_download()

        mock_logging_info.assert_any_call(f"Searching for keyword: {key_word}.")
        mock_logging_info.assert_any_call("Finished writing genome information to metadata_file")
        mock_logging_info.assert_any_call(f"Downloading genomes to {output_folder}")

    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    @patch('logging.info')
    @patch('os.mkdir')
    def test_parse_download_only_representative_genomes(self, mock_mkdir, mock_logging_info, mock_exists, mock_gzip_open, mock_open):
        key_word = "Enterococcus"
        output_folder = "my_output_folder"
        sys.argv = ["motus", "download", "-s", "metadata_file", "-w", key_word, "-r", "-o", output_folder]

        with self.assertRaises(SystemExit):
            parse_download()
        mock_logging_info.assert_any_call(f"Searching for keyword: {key_word}.")
        mock_logging_info.assert_any_call("Finished writing genome information to metadata_file")
        mock_logging_info.assert_any_call(f"Downloading genomes to {output_folder}")

    def test_parse_download_no_parameters(self):
        sys.argv = ["motus", "download"]

        with self.assertRaises(SystemExit):
            sys.stderr = io.StringIO()
            parse_download()
        self.assertIn('error: the following arguments are required: -s, -w', sys.stderr.getvalue())