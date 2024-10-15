from motus.motus import parse_downloadDB
import unittest
from unittest.mock import patch, MagicMock
import sys


class TestMotusDownloadDB(unittest.TestCase):

    @patch('motus.motus.DEFAULT_MOTUS_MGDB_LOCATION_MARKER', new_callable=MagicMock)
    @patch('urllib.request.urlretrieve')
    @patch('tarfile.open')
    @patch('shutil.rmtree')
    @patch('logging.info')
    def test_download_database_no_force(self, mock_logging_info, mock_rmtree, mock_tarfile, mock_urlretrieve,
                                        mock_marker):
        sys.argv = ['motus', 'downloadDB']

        # simulate that the marker file exists, which means the database is already downloaded
        mock_marker.exists.return_value = True

        # run code and check that system exit is raised and correct logs are written
        with self.assertRaises(SystemExit):
            parse_downloadDB()
        mock_logging_info.assert_any_call('Database already downloaded and -f not set. All good.')
        mock_logging_info.assert_any_call('mOTU tool shutting down with exitcode 0')

    @patch('motus.motus.DEFAULT_MOTUS_MGDB_LOCATION_MARKER', new_callable=MagicMock)
    @patch('motus.motus.DEFAULT_MOTUS_MGDB_LOCATION', new_callable=MagicMock)
    @patch('motus.motus.DEFAULT_MOTUS_MGDB_PARENT_LOCATION', new_callable=MagicMock)
    @patch('motus.motus.MOTUS_MGDB_REMOTE_LOCATION', new_callable=MagicMock)
    @patch('urllib.request.urlretrieve')
    @patch('tarfile.open')
    @patch('shutil.rmtree')
    @patch('logging.info')
    def test_download_database_force_flag(self, mock_logging_info, mock_rmtree, mock_tarfile, mock_urlretrieve,
                                          mock_remote_location, mock_parent_location, mock_location,
                                          mock_marker):
        # simulate that the marker file exists and force flag is set
        mock_marker.exists.return_value = True
        mock_location.exists.return_value = True

        # run the code
        sys.argv = ['motus', 'downloadDB', '-f']
        with self.assertRaises(SystemExit):
            parse_downloadDB()

        # force flag is set, database exists, so it should delete and re-download
        mock_logging_info.assert_any_call(
            'Database already downloaded and -f set. Will delete current database and download again.')
        mock_logging_info.assert_any_call(
            'Finished untaring mOTUs marker gene database.')
        mock_rmtree.assert_called_with(mock_location)
        mock_urlretrieve.assert_called_once_with(mock_remote_location,
                                                 str(mock_parent_location.joinpath('db_mOTU.tar.gz')))
        mock_marker.touch.assert_called_once()
        mock_tarfile.assert_called_once()

    @patch('motus.motus.DEFAULT_MOTUS_MGDB_LOCATION_MARKER', new_callable=MagicMock)
    @patch('motus.motus.DEFAULT_MOTUS_MGDB_PARENT_LOCATION', new_callable=MagicMock)
    @patch('motus.motus.MOTUS_MGDB_REMOTE_LOCATION', new_callable=MagicMock)
    @patch('sys.argv', new=['motus', 'downloadDB'])
    @patch('urllib.request.urlretrieve')  # Mock downloading the database
    @patch('tarfile.open')  # Mock tarfile for extracting
    @patch('shutil.rmtree')  # Mock shutil.rmtree for deleting directories
    @patch('logging.info')  # Mock logging info
    def test_download_database_first_time(self, mock_logging_info, mock_rmtree, mock_tarfile, mock_urlretrieve,
                                          mock_remote_location, mock_parent_location, mock_marker):
        # simulate that the marker file does not exist, meaning database has not been downloaded yet
        mock_marker.exists.return_value = False

        # run the code
        with self.assertRaises(SystemExit):
            parse_downloadDB()

        # database doesn't exist, so it should download
        mock_logging_info.assert_any_call('Start downloading mOTUs marker gene database. ~6GB')
        mock_logging_info.assert_any_call('Finished untaring mOTUs marker gene database.')
        mock_urlretrieve.assert_called_once_with(mock_remote_location,
                                                 str(mock_parent_location.joinpath('db_mOTU.tar.gz')))
        mock_marker.touch.assert_called_once()
        mock_tarfile.assert_called_once()

    @patch('motus.motus.DEFAULT_MOTUS_MGDB_LOCATION_MARKER', new_callable=MagicMock)
    @patch('motus.motus.DEFAULT_MOTUS_MGDB_LOCATION', new_callable=MagicMock)
    @patch('motus.motus.DEFAULT_MOTUS_MGDB_PARENT_LOCATION', new_callable=MagicMock)
    @patch('motus.motus.MOTUS_MGDB_REMOTE_LOCATION', new_callable=MagicMock)
    @patch('urllib.request.urlretrieve')
    @patch('tarfile.open')
    @patch('shutil.rmtree')
    @patch('logging.info')
    def test_download_database_force_flag_MGDB_location_does_not_exist(self, mock_logging_info, mock_rmtree,
                                                                       mock_tarfile, mock_urlretrieve,
                                                                       mock_remote_location, mock_parent_location,
                                                                       mock_location,
                                                                       mock_marker):
        # simulate that the marker file exists and force flag is set
        mock_marker.exists.return_value = True
        mock_location.exists.return_value = False

        # run the code
        sys.argv = ['motus', 'downloadDB', '-f']
        with self.assertRaises(SystemExit):
            parse_downloadDB()

        # force flag is set, database exists, so it should delete and re-download
        mock_logging_info.assert_any_call(
            'Database already downloaded and -f set. Will delete current database and download again.')
        mock_logging_info.assert_any_call(
            'Finished untaring mOTUs marker gene database.')
        mock_rmtree.assert_called_with(mock_location)
        mock_urlretrieve.assert_called_once_with(mock_remote_location,
                                                 str(mock_parent_location.joinpath('db_mOTU.tar.gz')))
        mock_marker.touch.assert_called_once()
        mock_tarfile.assert_called_once()


if __name__ == '__main__':
    unittest.main()
