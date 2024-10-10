import io
import unittest
from unittest.mock import patch, mock_open, MagicMock
import sys
from motus.motus import parse_batch_profile



class TestParseBatchProfile(unittest.TestCase):

    def mock_version_and_tsv_file(self, tsv_file_content, version_file_content):
        tsv_file_mock = MagicMock()
        tsv_file_mock.__enter__.return_value = tsv_file_mock
        tsv_file_mock.__iter__.return_value = tsv_file_content.splitlines().__iter__()  # Simulates 'for line in file'
        version_file_mock = MagicMock()
        version_file_mock.__enter__.return_value = version_file_mock
        version_file_mock.readline.side_effect = version_file_content.splitlines()

        # Define side_effects for different file paths
        def mock_file_selector(file, mode='r'):
            if "input_file.tsv" in str(file):
                return tsv_file_mock
            if "mOTUsv4.0.db" in str(file):
                return version_file_mock

        # Set the side effect for open to simulate reading different files
        mock_open.side_effect = mock_file_selector

        return version_file_mock, tsv_file_mock, mock_open.side_effect

    def test_no_input_file_provided(self):
        # Test for case where no input files are provided, expect an error
        sys.argv = ["motus", "batch_profile", "-g", "1"]
        with self.assertRaises(SystemExit):
            sys.stderr = io.StringIO()
            parse_batch_profile()
        self.assertIn('error: the following arguments are required: -f', sys.stderr.getvalue())

    def test_no_arguments_provided(self):
        # Test for case where no input files are provided, expect an error
        sys.argv = ["motus", "batch_profile"]
        with self.assertRaises(SystemExit):
            sys.stderr = io.StringIO()
            parse_batch_profile()
        self.assertIn('error: the following arguments are required: -f', sys.stderr.getvalue())

    @patch('builtins.open', new_callable=MagicMock)
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_all_arguments_provided(self, mock_exists, mock_gzip_open, mock_open):
        # pass all arguments possible to the function, since the input file is empty, it will not run all steps
        sys.argv = ["motus", "batch_profile", "-f", "input_file.tsv", "-l", "5", "-t", "80", "-v", "6", "-y", "INSERT_NORM"]

        mock_open.return_value.__enter__.return_value = MagicMock(readline=MagicMock(side_effect=[
            "version: 4.0", "date: 2024-01-01"
        ]))

        with self.assertLogs('root', level='INFO') as log_capture:
            with self.assertRaises(SystemExit):
                parse_batch_profile()
            self.assertIn('INFO:root:mOTU tool shutting down with exitcode 0', log_capture.output)


    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_map_file_parsing(self, mock_exists, mock_gzip_open, mock_open):
        # Simulate the command-line arguments passed to the function
        sys.argv = ["motus", "batch_profile", "-f", "input_file.tsv", "-l", "5", "-t", "80", "-v", "6", "-y",
                    "INSERT_NORM"]

        # Mock the tsv file and version file content
        tsv_file_content = "SAMPLE-1\tmy_file1.bam\nSAMPLE-2\tmy_file2.bam\n"
        version_file_content = "Version: 4.0.0\nDate: 2024-10-10\n"

        # Create a mock for the 'input_file.tsv'
        version_file_mock, tsv_file_mock, mock_open.side_effect = self.mock_version_and_tsv_file(tsv_file_content, version_file_content)

        #with self.assertLogs('root', level='ERROR') as log_capture:
        with self.assertRaises(FileNotFoundError):
            parse_batch_profile()
            #self.assertIn("[Errno 2] No such file or directory: 'my_file1.bam'", log_capture.output)

    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_map_file_parsing_no_bam_ending(self, mock_exists, mock_gzip_open, mock_open):
        # Simulate the command-line arguments passed to the function
        sys.argv = ["motus", "batch_profile", "-f", "input_file.tsv", "-l", "5", "-t", "80", "-v", "6", "-y",
                    "INSERT_NORM"]

        # Mock the tsv file and version file content
        tsv_file_content = "SAMPLE-1\tmy_file1\nSAMPLE-2\tmy_file2.bam\n"
        version_file_content = "Version: 4.0.0\nDate: 2024-10-10\n"

        version_file_mock, tsv_file_mock, mock_open.side_effect = self.mock_version_and_tsv_file(tsv_file_content, version_file_content)

        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                parse_batch_profile()
            self.assertIn('ERROR:root:Submitted BAM file my_file1 does not end with .bam. Probably malformed file. Quitting ...', log_capture.output)

    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    def test_non_existent_bam_file(self, mock_gzip_open, mock_open):
        # Simulate the command-line arguments passed to the function
        sys.argv = ["motus", "batch_profile", "-f", "input_file.tsv", "-l", "5", "-t", "80", "-v", "6", "-y",
                    "INSERT_NORM"]

        # Mock the tsv file and version file content
        tsv_file_content = "SAMPLE-1\tmy_file1.bam\nSAMPLE-2\tmy_file2.bam\n"
        version_file_content = "Version: 4.0.0\nDate: 2024-10-10\n"

        version_file_mock, tsv_file_mock, mock_open.side_effect = self.mock_version_and_tsv_file(tsv_file_content, version_file_content)

        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                parse_batch_profile()
            self.assertIn('ERROR:root:Submitted BAM file my_file1.bam does not exist. Quitting ...', log_capture.output)

    @patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_duplicated_sample_names(self, mock_exists, mock_gzip_open, mock_open):
        # Simulate the command-line arguments passed to the function
        sys.argv = ["motus", "batch_profile", "-f", "input_file.tsv", "-l", "5", "-t", "80", "-v", "6", "-y",
                    "INSERT_NORM"]

        # Mock the tsv file and version file content
        tsv_file_content = "SAMPLE-1\tmy_file1.bam\nSAMPLE-1\tmy_file2.bam\nSAMPLE-3\tmy_file2.bam\n"
        version_file_content = "Version: 4.0.0\nDate: 2024-10-10\n"

        version_file_mock, tsv_file_mock, mock_open.side_effect = self.mock_version_and_tsv_file(tsv_file_content, version_file_content)

        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                parse_batch_profile()
            self.assertIn('ERROR:root:Submitted samplename SAMPLE-1 duplicated. Quitting ...', log_capture.output)


    # if 1965 changes to "if os.path.samefile(obf,bamfile)", this will work --> then we can test the remainder of the function
    """@patch('builtins.open')
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    @patch('os.path.samefile', return_value=True)
    def test_duplicated_bam_files(self, mock_samefile, mock_exists, mock_gzip_open, mock_open):
        # Simulate the command-line arguments passed to the function
        sys.argv = ["motus", "batch_profile", "-f", "input_file.tsv", "-l", "5", "-t", "80", "-v", "6", "-y",
                    "INSERT_NORM"]

        # Mock the tsv file and version file content
        tsv_file_content = "SAMPLE-1\tmy_file1.bam\nSAMPLE-2\tmy_file1.bam\n"
        version_file_content = "Version: 4.0.0\nDate: 2024-10-10\n"

        # Create a mock for the 'input_file.tsv'
        version_file_mock, tsv_file_mock, mock_open.side_effect = self.mock_version_and_tsv_file(tsv_file_content, version_file_content)

        with self.assertLogs('root', level='INFO') as log_capture:
            with self.assertRaises(SystemExit):
                parse_batch_profile()
            self.assertIn('ERROR:root:Submitted BAM file my_file1.bam', log_capture.output)"""

if __name__ == '__main__':
    unittest.main()
