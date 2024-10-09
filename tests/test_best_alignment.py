import unittest
from unittest.mock import patch
from motus.motus import BestAlignment
from typing import List, Tuple


# Assuming the BestAlignment class is defined in a module called alignment_module
# from alignment_module import BestAlignment

class TestBestAlignment(unittest.TestCase):

    def setUp(self):
        """Set up a new instance of BestAlignment for each test."""
        self.ba = BestAlignment()

    def test_initialization(self):
        """Test that BestAlignment initializes with an empty dictionary."""
        self.assertEqual(self.ba._mg_2_blocks, {})

    def test_append_single_entry(self):
        """Test appending a single marker gene and its blocks."""
        mg = 'marker_gene_1'
        blocks = [(100, 200), (300, 400)]

        self.ba.append(mg, blocks)
        self.assertIn(mg, self.ba._mg_2_blocks)
        self.assertEqual(self.ba._mg_2_blocks[mg], blocks)

    def test_append_multiple_entries(self):
        """Test appending multiple marker genes and their blocks."""
        mg1 = 'marker_gene_1'
        blocks1 = [(100, 200), (300, 400)]
        mg2 = 'marker_gene_2'
        blocks2 = [(500, 600), (700, 800)]

        self.ba.append(mg1, blocks1)
        self.ba.append(mg2, blocks2)

        self.assertIn(mg1, self.ba._mg_2_blocks)
        self.assertIn(mg2, self.ba._mg_2_blocks)
        self.assertEqual(self.ba._mg_2_blocks[mg1], blocks1)
        self.assertEqual(self.ba._mg_2_blocks[mg2], blocks2)

    def test_is_multimapper_false(self):
        """Test isMultimapper returns False for a single marker gene."""
        self.ba.append('marker_gene_1', [(100, 200)])

        self.assertFalse(self.ba.isMultimapper())

    def test_is_multimapper_true(self):
        """Test isMultimapper returns True for multiple marker genes."""
        self.ba.append('marker_gene_1', [(100, 200)])
        self.ba.append('marker_gene_2', [(300, 400)])

        self.assertTrue(self.ba.isMultimapper())

    def test_get_mg_and_blocks_unique(self):
        """Test get_mg_and_blocks returns correct data for unique mappers."""
        mg = 'marker_gene_1'
        blocks = [(100, 200), (300, 400)]
        self.ba.append(mg, blocks)

        result = self.ba.get_mg_and_blocks()
        self.assertEqual(result, (mg, blocks))

        with self.assertNoLogs('root', level='ERROR'):
            self.ba.get_mg_and_blocks()


    def test_get_mg_and_blocks_multimapper(self):
        """Test get_mg_and_blocks fails for multimappers."""
        self.ba.append('marker_gene_1', [(100, 200)])
        self.ba.append('marker_gene_2', [(300, 400)])

        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.ba.get_mg_and_blocks()
            self.assertEqual(log_capture.output,
                             ['ERROR:root:This method doesnt work for multi mappers.'])


    def test_get_mgs_and_blocks_multimapper(self):
        """Test get_mgs_and_blocks returns correct data for multimappers."""
        mg1 = 'marker_gene_1'
        blocks1 = [(100, 200)]
        mg2 = 'marker_gene_2'
        blocks2 = [(300, 400)]

        self.ba.append(mg1, blocks1)
        self.ba.append(mg2, blocks2)

        with self.assertNoLogs('root', level='ERROR'):
            result = self.ba.get_mgs_and_blocks()
            expected = {
                mg1: blocks1,
                mg2: blocks2
            }

            self.assertEqual(result, expected)


    def test_get_mgs_and_blocks_unique(self):
        """Test get_mgs_and_blocks fails for unique mappers."""
        mg = 'marker_gene_1'
        blocks = [(100, 200)]
        self.ba.append(mg, blocks)

        with self.assertLogs('root', level='ERROR') as log_capture:
            with self.assertRaises(SystemExit):
                self.ba.get_mgs_and_blocks()
                self.assertEqual(log_capture.output,
                             ['ERROR:root:This method doesnt work for unique mappers.'])


if __name__ == '__main__':
    unittest.main()
