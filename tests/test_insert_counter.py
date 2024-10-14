from motus.motus import InsertCounter, BestAlignment
import unittest
from unittest.mock import MagicMock
import collections


class TestInsertCounter(unittest.TestCase):

    def setUp(self):
        self.insert_counter = InsertCounter()
        self.ba = BestAlignment()

    def test_initial_state(self):
        """Test that the initial state of the InsertCounter is as expected."""
        self.assertEqual(self.insert_counter.get_unique_mapper_count(), 0)
        self.assertEqual(self.insert_counter.get_multi_mapper_count(), 0)

    def test_append_mapper_unique(self):
        """Test appending a unique mapper."""
        self.ba.append('MG1', [(100, 200)])
        self.insert_counter.appendmapper('insert1', self.ba)

        self.assertEqual(self.insert_counter.get_unique_mapper_count(), 1)
        self.assertEqual(self.insert_counter.get_multi_mapper_count(), 0)

    def test_append_mapper_multi_mapper(self):
        """Test appending a multi mapper."""
        self.ba.append('MG1', [(100, 200)])
        self.ba.append('MG2', [(300, 400)])
        self.insert_counter.appendmapper('insert1', self.ba)

        self.assertEqual(self.insert_counter.get_unique_mapper_count(), 0)
        self.assertEqual(self.insert_counter.get_multi_mapper_count(), 1)

    """def test_correct_uniq_mapper_edges(self):
        #Test correcting unique mapper edges.
        # Setup mock for file writer
        inserts_file_writer = StringIO()

        # Mock behavior of BestAlignment
        self.ba.append('MG1', [(0, 100)])

        # Add a unique mapper
        self.insert_counter.appendmapper('insert1', self.ba)

        # Perform the edge correction
        self.insert_counter.correct_uniq_mapper_edges(inserts_file_writer, 30)

        # Validate results
        expected_output = 'insert1\tMG1\t1.0000\n'
        self.assertEqual(inserts_file_writer.getvalue(), expected_output)"""

    """@patch('motus.motus.MotusDB.get_length_by_mg')
    @patch('motus.motus.InsertCounter._get_edge_corrected_raw_uniquemapper_insert_counts_per_mgc', return_value = collections.Counter("[(0, 30)]"))
    @patch('motus.motus.MotusDB.get_mgc_by_mg', side_effect=lambda mg: 'mgc1')
    def test_correct_multi_mapper_edges(self, mock_get_length, mock_edge_corrected, mock_lambda):
        #Test correcting multi mapper edges.
        inserts_file_writer = StringIO()
        self.ba.append('MG1', [(0, 30)])
        self.ba.append('MG2', [(0, 40)])
        self.insert_counter.appendmapper('insert1', self.ba)

        # Mock lengths of MGs
        mock_get_length.side_effect = lambda mg: 100 if mg in ['MG1', 'MG2'] else 0

        # Prepare unique mapper counts
        self.insert_counter._mg_2_edge_corrected_raw_uniquemapper_insert_counts = {
            'MG1': 10,
            'MG2': 20
        }

        # Perform the edge correction
        self.insert_counter.correct_multi_mapper_edges(inserts_file_writer, 30)

        expected_output = 'insert1\tMG1\t0.3333\ninsert1\tMG2\t0.6667\n'
        self.assertIn(expected_output, inserts_file_writer.getvalue())"""

    def test_combined_raw_counts(self):
        """Test the combined raw counts function."""
        self.insert_counter._mg_2_edge_corrected_raw_uniquemapper_insert_counts = {'MG1': 10}
        self.insert_counter._mg_2_edge_corrected_raw_multimapper_insert_counts = {'MG1': 5}

        self.insert_counter.combined_raw_counts()

        self.assertEqual(self.insert_counter._mg_2_edge_corrected_raw_insert_counts['MG1'], 15)

    """@patch('motus.motus.MotusDB.get_length_by_mg')
    def test_norm_and_scale_counts(self, mock_get_length):
        #Test normalization and scaling of counts.
        # Setup mock lengths
        #mock_get_length.side_effect = lambda mg: 100 if mg == 'MG1' else 200

        self.insert_counter._mg_2_edge_corrected_raw_insert_counts = {'MG1': 10}
        self.insert_counter.norm_and_scale_counts()

        # Validate normalized and scaled counts
        expected_norm = {'MG1': 0.1}
        expected_scaled = {'MG1': 10.0}

        self.assertEqual(self.insert_counter._mg_2_edge_corrected_norm_insert_counts, expected_norm)
        self.assertEqual(self.insert_counter._mg_2_edge_corrected_scaled_insert_counts, expected_scaled)"""

    def test_combined_raw_counts_empty(self):
        """Test when both unique and multimapper counts are empty"""
        self.insert_counter._mg_2_edge_corrected_raw_uniquemapper_insert_counts = collections.Counter()
        self.insert_counter._mg_2_edge_corrected_raw_multimapper_insert_counts = collections.Counter()
        self.insert_counter._mg_2_edge_corrected_raw_uniquemapper_base_counts = collections.Counter()
        self.insert_counter._mg_2_edge_corrected_raw_multimapper_base_counts = collections.Counter()

        self.insert_counter.combined_raw_counts()

        # both combined counts should be empty
        self.assertEqual(self.insert_counter._mg_2_edge_corrected_raw_insert_counts, collections.Counter())
        self.assertEqual(self.insert_counter._mg_2_edge_corrected_raw_base_counts, collections.Counter())

    def test_combined_raw_counts_unique_only(self):
        """Test when only unique mapper counts are present"""
        self.insert_counter._mg_2_edge_corrected_raw_uniquemapper_insert_counts = collections.Counter(
            {'MG1': 100, 'MG2': 150})
        self.insert_counter._mg_2_edge_corrected_raw_multimapper_insert_counts = collections.Counter()  # Empty
        self.insert_counter._mg_2_edge_corrected_raw_uniquemapper_base_counts = collections.Counter(
            {'MG1': 500, 'MG2': 600})
        self.insert_counter._mg_2_edge_corrected_raw_multimapper_base_counts = collections.Counter()  # Empty

        self.insert_counter.combined_raw_counts()

        # combined counts should equal the unique counts
        self.assertEqual(self.insert_counter._mg_2_edge_corrected_raw_insert_counts, {'MG1': 100, 'MG2': 150})
        self.assertEqual(self.insert_counter._mg_2_edge_corrected_raw_base_counts, {'MG1': 500, 'MG2': 600})

    def test_combined_raw_counts_multimapper_only(self):
        """Test when only multimapper counts are present"""
        self.insert_counter._mg_2_edge_corrected_raw_uniquemapper_insert_counts = collections.Counter()  # Empty
        self.insert_counter._mg_2_edge_corrected_raw_multimapper_insert_counts = collections.Counter(
            {'MG1': 80, 'MG2': 120})
        self.insert_counter._mg_2_edge_corrected_raw_uniquemapper_base_counts = collections.Counter()  # Empty
        self.insert_counter._mg_2_edge_corrected_raw_multimapper_base_counts = collections.Counter(
            {'MG1': 400, 'MG2': 500})

        self.insert_counter.combined_raw_counts()

        # combined counts should equal the multimapper counts
        self.assertEqual(self.insert_counter._mg_2_edge_corrected_raw_insert_counts, {'MG1': 80, 'MG2': 120})
        self.assertEqual(self.insert_counter._mg_2_edge_corrected_raw_base_counts, {'MG1': 400, 'MG2': 500})

    def test_combined_raw_counts_multiple_mgs(self):
        """Test with multiple `mg` values"""
        self.insert_counter._mg_2_edge_corrected_raw_uniquemapper_insert_counts = collections.Counter(
            {'MG1': 100, 'MG2': 150, 'MG3': 200})
        self.insert_counter._mg_2_edge_corrected_raw_multimapper_insert_counts = collections.Counter(
            {'MG1': 50, 'MG2': 75, 'MG3': 100})
        self.insert_counter._mg_2_edge_corrected_raw_uniquemapper_base_counts = collections.Counter(
            {'MG1': 500, 'MG2': 600, 'MG3': 700})
        self.insert_counter._mg_2_edge_corrected_raw_multimapper_base_counts = collections.Counter(
            {'MG1': 250, 'MG2': 300, 'MG3': 350})

        self.insert_counter.combined_raw_counts()

        # check combined counts for multiple `mg`s
        self.assertEqual(self.insert_counter._mg_2_edge_corrected_raw_insert_counts,
                         {'MG1': 150, 'MG2': 225, 'MG3': 300})
        self.assertEqual(self.insert_counter._mg_2_edge_corrected_raw_base_counts,
                         {'MG1': 750, 'MG2': 900, 'MG3': 1050})

    def test_get_alignment_blocks_empty_input(self):
        """Test when the input is an empty list"""
        alignments = []
        result = self.insert_counter._get_alignment_blocks(alignments)
        self.assertEqual(result, [])

    def test_get_alignment_blocks_single_block(self):
        """Test with a single alignment having one block"""
        mock_alignment = MagicMock()
        mock_alignment.get_blocks.return_value = [(100, 150)]  # Mock one alignment block

        alignments = [mock_alignment]
        result = self.insert_counter._get_alignment_blocks(alignments)

        # expect one block returned as is
        self.assertEqual(result, [(100, 150)])

    def test_get_alignment_blocks_multiple_blocks(self):
        """Test with a single alignment having multiple blocks (simulating an indel)"""
        mock_alignment = MagicMock()
        mock_alignment.get_blocks.return_value = [(100, 150), (200, 250)]  # Mock multiple blocks

        alignments = [mock_alignment]
        result = self.insert_counter._get_alignment_blocks(alignments)

        # expect multiple blocks returned as is, in sorted order
        self.assertEqual(result, [(100, 150), (200, 250)])

    def test_get_alignment_blocks_multiple_alignments(self):
        """Test with multiple alignments having multiple blocks"""
        mock_alignment1 = MagicMock()
        mock_alignment1.get_blocks.return_value = [(100, 150)]  # First alignment block
        mock_alignment2 = MagicMock()
        mock_alignment2.get_blocks.return_value = [(200, 250), (300, 350)]  # Second alignment with two blocks

        alignments = [mock_alignment1, mock_alignment2]
        result = self.insert_counter._get_alignment_blocks(alignments)

        # expect all blocks sorted by their start position
        self.assertEqual(result, [(100, 150), (200, 250), (300, 350)])

    def test_get_alignment_blocks_unsorted_blocks(self):
        """Test with multiple blocks in unsorted order"""
        mock_alignment1 = MagicMock()
        mock_alignment1.get_blocks.return_value = [(300, 350)]  # This block starts later
        mock_alignment2 = MagicMock()
        mock_alignment2.get_blocks.return_value = [(100, 150), (200, 250)]  # These blocks start earlier

        alignments = [mock_alignment1, mock_alignment2]
        result = self.insert_counter._get_alignment_blocks(alignments)

        # result should be sorted by start position, so the result should be:
        self.assertEqual(result, [(100, 150), (200, 250), (300, 350)])


if __name__ == '__main__':
    unittest.main()
