import unittest
from unittest.mock import MagicMock, patch
import pathlib
from motus import motus


class TestMotusDB(unittest.TestCase):
    @patch('builtins.open', new_callable=MagicMock)
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def setUp(self, mock_exists, mock_gzip_open, mock_open):
        # setup mock data for testing initialization (versions file)
        mock_open.return_value.__enter__.return_value = MagicMock(readline=MagicMock(side_effect=[
            "version: 4.0", "date: 2024-01-01"
        ]))

        # setup dummy folder path
        self.db_folder = pathlib.Path("/path/to/db")

        # mock gzip file reads
        mgc_data = "MG\tMGC\tLENGTH\t#MOTU\tCOG\nMG1\tMGC1\t100\tMOTU1\tCOG0012\n"
        mock_gzip_open.return_value.__enter__.return_value = MagicMock(read=MagicMock(return_value=mgc_data))
        self.motus_db = motus.MotusDB(self.db_folder, load=True)
        self.obj = motus.MotusDB

    def test_initialization(self):
        """Test that initialization works and loads the database version correctly."""
        self.assertEqual(self.motus_db.database_version, '4.0')
        self.assertEqual(self.motus_db.database_date, '2024-01-01')
        self.assertEqual(self.motus_db.index_location, pathlib.Path('/path/to/db/mOTUsv4.0.db.fna.gz'))

    @patch('builtins.open', new_callable=MagicMock)
    def test_index_file_does_not_exist(self, mock_open):
        """Test that it raises an error when the path to the mOTUs database does not exist."""
        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.obj(pathlib.Path("/path/does/not/exist"), load=True)
            string = "ERROR:root:Database file /path/does/not/exist/mOTUsv4.0.db.fna.gz is missing. Quitting mOTUs..."
            self.assertIn(log_capture.output[0], string)

    def test_is_mg_blocked(self):
        """Test that blocked marker genes are correctly identified."""
        self.motus_db.blocklist_mg.add('MG1')
        self.assertTrue(self.motus_db.is_mg_blocked('MG1'))
        self.assertFalse(self.motus_db.is_mg_blocked('MG2'))

    @patch('motus.motus.MotusDB.database_version', '4.0')
    def test_get_full_version(self):
        """Test the get_full_version method."""
        self.assertEqual(self.motus_db.get_full_version(), 'TOOL:4.0.0_DB:4.0')

    def test_get_full_sam_id(self):
        """Test the get_full_sam_id method."""
        self.assertEqual(self.motus_db.get_full_sam_id(), 'mOTUs4')

    def test_get_mg_by_mgc(self):
        """Test getting markergene by MGC."""
        self.motus_db.mgc_2_mg['MGC'] = 'COG0012'
        self.assertEqual(self.motus_db.get_mg_by_mgc('MGC'), 'COG0012')

    def test_is_unassigned_motu(self):
        """Test if a MOTU is unassigned."""
        self.motus_db._unassigned_motu_name = 'MOTU_unassigned'
        self.assertTrue(self.motus_db.is_unassigned_motu('MOTU_unassigned'))
        self.assertFalse(self.motus_db.is_unassigned_motu('MOTU1'))

    def test_get_unassigned_motu(self):
        """Test getting the unassigned MOTU name."""
        self.motus_db._unassigned_motu_name = 'MOTU_unassigned'
        self.assertEqual(self.motus_db.get_unassigned_motu(), 'MOTU_unassigned')

    def test_get_motu_by_mgc(self):
        """Test getting MOTU by MGC."""
        self.motus_db.mgc_2_motu['MGC'] = 'MOTU1'
        self.assertEqual(self.motus_db.get_motu_by_mgc('MGC'), 'MOTU1')

    def test_get_bwa_index(self):
        """Test getting the BWA index location."""
        self.motus_db.index_location = "/path/to/index"
        self.assertEqual(self.motus_db.get_bwa_index(), "/path/to/index")

    def test_get_mgc_by_mg(self):
        """Test getting MGC by MG."""
        self.motus_db.mgh_2_mgc['MG1'] = 'MGC1'
        self.assertEqual(self.motus_db.get_mgc_by_mg('MG1'), 'MGC1')

    def test_get_length_by_mg(self):
        """Test getting length of MG."""
        self.motus_db.mgh_2_mglength['MG1'] = 100
        self.assertEqual(self.motus_db.get_length_by_mg('MG1'), 100)

    def test_get_core_motus_mgs(self):
        """Test getting core MOTU MGs."""
        self.assertEqual(self.motus_db.get_core_motus_mgs(),
                         ['COG0012', 'COG0016', 'COG0018', 'COG0172', 'COG0215', 'COG0495', 'COG0525', 'COG0533',
                          'COG0541', 'COG0552'])


if __name__ == '__main__':
    unittest.main()
