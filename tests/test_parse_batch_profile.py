import builtins
import io
import pathlib
import unittest
from unittest.mock import patch, MagicMock
import sys
from motus.motus import parse_batch_profile

TOY_DB = pathlib.Path(__file__).parent / 'data' / 'motus4.1-toy-db'

_real_open = builtins.open


def mock_version_and_tsv_file(tsv_file_content, version_file_content, mock_open):
    """Set up mock_open so that 'input_file.tsv' returns fake TSV content and all
    other file opens (e.g. the toy DB version file) fall through to the real open."""
    tsv_file_mock = MagicMock()
    tsv_file_mock.__enter__.return_value = tsv_file_mock
    tsv_file_mock.__iter__.return_value = iter(tsv_file_content.splitlines())

    def mock_file_selector(file, mode='r', *args, **kwargs):
        if "input_file.tsv" in str(file):
            return tsv_file_mock
        return _real_open(file, mode, *args, **kwargs)

    mock_open.side_effect = mock_file_selector
    return None, tsv_file_mock, mock_open.side_effect


class TestParseBatchProfile(unittest.TestCase):

    def test_no_input_file_provided(self):
        # test for case where no input files are provided, expect an error
        sys.argv = ["motus", "batch_profile", "-g", "1"]
        with self.assertRaises(SystemExit):
            sys.stderr = io.StringIO()
            parse_batch_profile()
        self.assertIn('error: the following arguments are required: -f', sys.stderr.getvalue())

    def test_no_arguments_provided(self):
        # test for case where no arguments are provided, expecting help to be shown
        sys.argv = ["motus", "batch_profile"]
        with self.assertRaises(SystemExit):
            parse_batch_profile()

    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    @patch('builtins.open')
    def test_all_arguments_provided(self, mock_open, mock_exists, mock_gzip_open):
        # pass all arguments possible to the function, since the input file is empty, it will not run all steps
        # -v (verbosity) was removed from the CLI in v4.1
        sys.argv = ["motus", "batch_profile", "-f", "input_file.tsv", "-l", "5", "-t", "80", "-y",
                    "INSERT_NORM", "-db", str(TOY_DB)]

        # Fake empty TSV; let all other opens (including the toy DB version file) go through.
        def selective_open(file, mode='r', *args, **kwargs):
            if 'input_file.tsv' in str(file):
                m = MagicMock()
                m.__enter__ = lambda s: s
                m.__exit__ = MagicMock(return_value=False)
                m.__iter__ = lambda s: iter([])
                return m
            return _real_open(file, mode, *args, **kwargs)

        mock_open.side_effect = selective_open

        with self.assertLogs('root', level='INFO') as log_capture:
            with self.assertRaises(SystemExit):
                parse_batch_profile()
            self.assertIn('INFO:root:mOTU tool shutting down with exitcode 0', log_capture.output)

    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_map_file_parsing(self, mock_exists, mock_gzip_open, mock_open):
        # -v (verbosity) was removed from the CLI in v4.1
        sys.argv = ["motus", "batch_profile", "-f", "input_file.tsv", "-l", "5", "-t", "80", "-y",
                    "INSERT_NORM", "-db", str(TOY_DB)]

        # mock the tsv file content; DB version file goes through real open
        tsv_file_content = "SAMPLE-1\tmy_file1.bam\nSAMPLE-2\tmy_file2.bam\n"
        version_file_content = "Version: 4.0.0\nDate: 2024-10-10\n"

        _, tsv_file_mock, mock_open.side_effect = mock_version_and_tsv_file(tsv_file_content,
                                                                             version_file_content, mock_open)

        with self.assertRaises(FileNotFoundError):
            parse_batch_profile()

    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_map_file_parsing_no_bam_ending(self, mock_exists, mock_gzip_open, mock_open):
        # pass wrong file ending; -v (verbosity) was removed from the CLI in v4.1
        sys.argv = ["motus", "batch_profile", "-f", "input_file.tsv", "-l", "5", "-t", "80", "-y",
                    "INSERT_NORM"]

        tsv_file_content = "SAMPLE-1\tmy_file1\nSAMPLE-2\tmy_file2.bam\n"
        version_file_content = "Version: 4.0.0\nDate: 2024-10-10\n"

        _, tsv_file_mock, mock_open.side_effect = mock_version_and_tsv_file(tsv_file_content,
                                                                             version_file_content, mock_open)

        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                parse_batch_profile()
            self.assertIn(
                'ERROR:root:Submitted BAM file my_file1 does not end with .bam. Probably malformed file. Quitting ...',
                log_capture.output)

    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    def test_non_existent_bam_file(self, mock_gzip_open, mock_open):
        # -v (verbosity) was removed from the CLI in v4.1
        sys.argv = ["motus", "batch_profile", "-f", "input_file.tsv", "-l", "5", "-t", "80", "-y",
                    "INSERT_NORM"]

        tsv_file_content = "SAMPLE-1\tmy_file1.bam\nSAMPLE-2\tmy_file2.bam\n"
        version_file_content = "Version: 4.0.0\nDate: 2024-10-10\n"

        _, tsv_file_mock, mock_open.side_effect = mock_version_and_tsv_file(tsv_file_content,
                                                                             version_file_content, mock_open)

        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                parse_batch_profile()
            self.assertIn('ERROR:root:Submitted BAM file my_file1.bam does not exist. Quitting ...', log_capture.output)

    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_duplicated_sample_names(self, mock_exists, mock_gzip_open, mock_open):
        # -v (verbosity) was removed from the CLI in v4.1
        sys.argv = ["motus", "batch_profile", "-f", "input_file.tsv", "-l", "5", "-t", "80", "-y",
                    "INSERT_NORM"]

        # duplicated sample name
        tsv_file_content = "SAMPLE-1\tmy_file1.bam\nSAMPLE-1\tmy_file2.bam\nSAMPLE-3\tmy_file2.bam\n"
        version_file_content = "Version: 4.0.0\nDate: 2024-10-10\n"

        _, tsv_file_mock, mock_open.side_effect = mock_version_and_tsv_file(tsv_file_content,
                                                                             version_file_content, mock_open)

        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                parse_batch_profile()
            self.assertIn('ERROR:root:Submitted samplename SAMPLE-1 duplicated. Quitting ...', log_capture.output)


if __name__ == '__main__':
    unittest.main()
