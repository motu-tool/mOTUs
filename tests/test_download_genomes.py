import unittest
from unittest.mock import patch, MagicMock, mock_open, call
from motus.motus import download_genomes
import pathlib


class TestDownloadGenomes(unittest.TestCase):

    @patch('motus.motus.urllib.request.urlretrieve')
    @patch('motus.motus.open', new_callable=mock_open)
    @patch('motus.motus.MotusSearchDB')
    @patch('motus.motus.pathlib.Path.mkdir')
    @patch('logging.info')
    def test_download_genomes_success(self, mock_logging_info, mock_mkdir, MockMotusSearchDB,
                                      mock_file, mock_urlretrieve):
        # set up the mock MotusSearchDB instance
        mock_db_instance = MockMotusSearchDB.return_value
        mock_db_instance.search_for_genomes.return_value = ['genome1', 'genome2']
        mock_db_instance.get_genome_path.side_effect = lambda genome: f"/path/to/{genome}"
        mock_db_instance.get_genome_motu.side_effect = lambda genome: f"motu_{genome}"
        mock_db_instance.get_genome_tax.side_effect = lambda \
                genome: "Domain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies"

        # prepare test data
        keyword = "test_keyword"
        output_folder = pathlib.Path("/output/folder")
        output_file = pathlib.Path("/output/file.txt")

        # call the function
        download_genomes(keyword, mock_db_instance, output_folder, output_file,
                         download_representative_genomes_only=True)

        # verify the genome search was performed
        mock_db_instance.search_for_genomes.assert_called_once_with(keyword, only_representatives=True)

        # verify genome information was written to file
        mock_file().write.assert_any_call('GENOME\tMOTU\tPATH\tDOMAIN\tPHYLUM\tCLASS\tORDER\tFAMILY\tGENUS\tSPECIES\n')
        mock_file().write.assert_any_call(
            'genome1\tmotu_genome1\t/path/to/genome1\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n')
        mock_file().write.assert_any_call(
            'genome2\tmotu_genome2\t/path/to/genome2\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n')

        # verify directory creation and downloads
        mock_mkdir.assert_called_once_with(exist_ok=True, parents=True)
        mock_urlretrieve.assert_has_calls([
            call("/path/to/genome1", "/output/folder/genome1"),
            call("/path/to/genome2", "/output/folder/genome2")
        ])

        # verify logging messages
        mock_logging_info.assert_any_call(f'Searching for keyword: {keyword}.')
        mock_logging_info.assert_any_call(f'Found: 2 hits.')
        mock_logging_info.assert_any_call(f'Finished writing genome information to {output_file}')
        mock_logging_info.assert_any_call(f'Downloading genome (1 / 2) genome1 to /output/folder/genome1')
        mock_logging_info.assert_any_call(f'Downloading genome (2 / 2) genome2 to /output/folder/genome2')
        mock_logging_info.assert_any_call(f'Finished downloading genomes')

    @patch('motus.motus.urllib.request.urlretrieve')
    @patch('motus.motus.open', new_callable=mock_open)
    @patch('motus.motus.MotusSearchDB')
    @patch('motus.motus.pathlib.Path.mkdir')
    @patch('logging.info')
    def test_download_genomes_success_all_genomes(self, mock_logging_info, mock_mkdir, MockMotusSearchDB,
                                      mock_file, mock_urlretrieve):
        # set up the mock MotusSearchDB instance
        mock_db_instance = MockMotusSearchDB.return_value
        mock_db_instance.search_for_genomes.return_value = ['genome1', 'genome2']
        mock_db_instance.get_genome_path.side_effect = lambda genome: f"/path/to/{genome}"
        mock_db_instance.get_genome_motu.side_effect = lambda genome: f"motu_{genome}"
        mock_db_instance.get_genome_tax.side_effect = lambda \
                genome: "Domain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies"

        # prepare test data
        keyword = "test_keyword"
        output_folder = pathlib.Path("/output/folder")
        output_file = pathlib.Path("/output/file.txt")

        # call the function
        download_genomes(keyword, mock_db_instance, output_folder, output_file,
                         download_representative_genomes_only=False)

        # verify the genome search was performed
        mock_db_instance.search_for_genomes.assert_called_once_with(keyword, only_representatives=False)

        # verify genome information was written to file
        mock_file().write.assert_any_call('GENOME\tMOTU\tPATH\tDOMAIN\tPHYLUM\tCLASS\tORDER\tFAMILY\tGENUS\tSPECIES\n')
        mock_file().write.assert_any_call(
            'genome1\tmotu_genome1\t/path/to/genome1\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n')
        mock_file().write.assert_any_call(
            'genome2\tmotu_genome2\t/path/to/genome2\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n')

        # verify directory creation and downloads
        mock_mkdir.assert_called_once_with(exist_ok=True, parents=True)
        mock_urlretrieve.assert_has_calls([
            call("/path/to/genome1", "/output/folder/genome1"),
            call("/path/to/genome2", "/output/folder/genome2")
        ])

        # verify logging messages
        mock_logging_info.assert_any_call(f'Searching for keyword: {keyword}.')
        mock_logging_info.assert_any_call(f'Found: 2 hits.')
        mock_logging_info.assert_any_call(f'Finished writing genome information to {output_file}')
        mock_logging_info.assert_any_call(f'Downloading genome (1 / 2) genome1 to /output/folder/genome1')
        mock_logging_info.assert_any_call(f'Downloading genome (2 / 2) genome2 to /output/folder/genome2')
        mock_logging_info.assert_any_call(f'Finished downloading genomes')

    @patch('motus.motus.MotusSearchDB')
    @patch('motus.motus.open', new_callable=mock_open)
    @patch('logging.info')
    def test_no_genomes_found(self, mock_logging_info, mock_file, MockMotusSearchDB):
        # set up the mock MotusSearchDB instance
        mock_db_instance = MockMotusSearchDB.return_value
        mock_db_instance.search_for_genomes.return_value = []

        # prepare test data
        keyword = "test_keyword"
        output_folder = None
        output_file = pathlib.Path("/output/file.txt")

        # call the function
        download_genomes(keyword, mock_db_instance, output_folder, output_file,
                         download_representative_genomes_only=True)

        # verify that the genome search was performed
        mock_db_instance.search_for_genomes.assert_called_once_with(keyword, only_representatives=True)

        # ensure no download occurs
        mock_file().write.assert_any_call('GENOME\tMOTU\tPATH\tDOMAIN\tPHYLUM\tCLASS\tORDER\tFAMILY\tGENUS\tSPECIES\n')

        # verify logging
        mock_logging_info.assert_any_call(f'Searching for keyword: {keyword}.')
        mock_logging_info.assert_any_call(f'Found: 0 hits.')
        mock_logging_info.assert_any_call(f'Finished writing genome information to {output_file}')

    @patch('motus.motus.urllib.request.urlretrieve')
    @patch('motus.motus.MotusSearchDB')
    @patch('logging.error')
    @patch('builtins.open')
    def test_output_folder_is_file(self, mock_open, mock_logging_error, MockMotusSearchDB, mock_urlretrieve):
        # setup mock MotusSearchDB instance
        mock_db_instance = MockMotusSearchDB.return_value
        mock_db_instance.search_for_genomes.return_value = ['genome1']

        # prepare test data
        keyword = "test_keyword"
        output_folder = MagicMock()
        output_folder.is_file.return_value = True  # Simulate output_folder being a file
        output_file = pathlib.Path("/output/file.txt")

        # call the function and expect it to raise a SystemExit
        with self.assertRaises(SystemExit):
            download_genomes(keyword, mock_db_instance, output_folder, output_file,
                             download_representative_genomes_only=True)

        # verify that the error was logged
        mock_logging_error.assert_called_once_with(
            'Output Path exists and is file. Cannot download genomes to this location')

    @patch('motus.motus.MotusSearchDB')
    @patch('motus.motus.open', new_callable=mock_open)
    @patch('logging.info')
    def test_only_write_genomes_no_download(self, mock_logging_info, mock_file, MockMotusSearchDB):
        # set up mock MotusSearchDB instance
        mock_db_instance = MockMotusSearchDB.return_value
        mock_db_instance.search_for_genomes.return_value = ['genome1']
        mock_db_instance.get_genome_path.side_effect = lambda genome: f"/path/to/{genome}"
        mock_db_instance.get_genome_motu.side_effect = lambda genome: f"motu_{genome}"
        mock_db_instance.get_genome_tax.side_effect = lambda \
                genome: "Domain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies"

        # prepare test data
        keyword = "test_keyword"
        output_folder = None  # No download, just write to file
        output_file = pathlib.Path("/output/file.txt")

        # call the function
        download_genomes(keyword, mock_db_instance, output_folder, output_file,
                         download_representative_genomes_only=True)

        # verify genome search
        mock_db_instance.search_for_genomes.assert_called_once_with(keyword, only_representatives=True)

        # verify file writing
        mock_file().write.assert_any_call('GENOME\tMOTU\tPATH\tDOMAIN\tPHYLUM\tCLASS\tORDER\tFAMILY\tGENUS\tSPECIES\n')
        mock_file().write.assert_any_call(
            'genome1\tmotu_genome1\t/path/to/genome1\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n')

        # ensure no directory creation or download happens
        mock_logging_info.assert_any_call(f'Finished writing genome information to {output_file}')
        self.assertNotIn("downloading", mock_logging_info)
        self.assertNotIn("Downloading", mock_logging_info)


if __name__ == '__main__':
    unittest.main()
