import unittest
from unittest.mock import Mock
from motus.motus import _get_orientation_of_aligned_segment_by_name
import pysam

SIDENTIFIER = 'S'


class TestGetOrientationOfAlignedSegmentByName(unittest.TestCase):

    def setUp(self):
        # Define a common setup for the AlignedSegment mock
        self.mock_alignment = Mock(spec=pysam.AlignedSegment)

    def test_standard_case(self):
        # Test the case where query_name contains a '/' to separate insert and orientation.
        # one '/' separator
        self.mock_alignment.query_name = 'insert_name/orientation'

        result = _get_orientation_of_aligned_segment_by_name(self.mock_alignment)
        expected = ('insert_name', 'orientation')

        self.assertEqual(result, expected)

    def test_edge_case_no_separator(self):
        # Test the case where query_name does not contain a '/'.
        # without a '/' separator
        self.mock_alignment.query_name = 'insert_name'

        result = _get_orientation_of_aligned_segment_by_name(self.mock_alignment)
        expected = ('insert_name', SIDENTIFIER)  # Assuming SIDENTIFIER is a predefined constant

        self.assertEqual(result, expected)

    def test_multiple_slashes_in_query_name(self):
        # Test the case where query_name contains multiple '/'.
        # multiple '/' characters
        self.mock_alignment.query_name = 'some/insert_name/with/multiple/slashes/orientation'

        # split on the last '/'
        result = _get_orientation_of_aligned_segment_by_name(self.mock_alignment)
        expected = ('some/insert_name/with/multiple/slashes', 'orientation')

        self.assertEqual(result, expected)

    def test_empty_query_name(self):
        # Test the case where query_name is empty.
        # empty string
        self.mock_alignment.query_name = ''

        result = _get_orientation_of_aligned_segment_by_name(self.mock_alignment)
        expected = ('', SIDENTIFIER)

        self.assertEqual(result, expected)


if __name__ == '__main__':
    unittest.main()
