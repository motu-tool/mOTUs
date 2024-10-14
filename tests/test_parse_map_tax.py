import io
import unittest
from unittest import mock
from unittest.mock import patch, MagicMock
import sys
import pathlib
from motus.motus import parse_map_tax
from io import StringIO
import pysam
from collections import OrderedDict


class PysamFakeBam:
    def __init__(self, header, reads):
        """
        a mock object that mimics the pysam.AlignmentFile object
        :param pysam.AlignmentHeader header: header of the mock sam file
        :param List[pysam.AlignedSegment] reads: reads of the mock sam file
        """
        self.header = header
        self.reads = reads

    def __iter__(self):
        return iter(self.reads)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return self

    def close(self):
        return self


def mock_bam_header(contig_list):
    """
    making a mock pysam.AlignmentHeader object
    Example::
        contigs = [("chr1", 10), ("chr2", 20)]
        mock_header = mock_bam_header(contigs)
    :param List[Tuple[str, int]] contig_list: a list of tuples of (contig name, contig length)
    :return: a pysam.AlignmentHeader object
    :rtype: pysam.AlignmentHeader
    """
    header_dict = OrderedDict(
        [
            ("SQ", [dict(SN=contig[0], LN=contig[1]) for contig in contig_list]),
            ("PG", [  # Program information
                {
                    "ID": "my_program",  # Unique ID for the program
                    "PN": "My Program",  # Program name
                    "VN": "1.0",  # Version
                    "CL": "python script.py"  # Command line used to run the program
                }
            ])
        ]
    )
    return pysam.AlignmentHeader.from_dict(header_dict)


def mock_alignment(
        header,
        reference_name,
        query_name,
        query_sequence,
        reference_start,
        cigar,
        flag,
        mapping_quality,
        next_reference_name=None,
        next_reference_start=None,
        is_unmapped=True,
):
    """
    making a mock pysam.AlignedSegment object
    :param pysam.AlignmentHeader header: a pysam alignment header object (can be created by mock_bam_header)
    :param str reference_name: reference name
    :param str query_name: query name
    :param str query_sequence: query sequence
    :param int reference_start: reference start
    :param list cigar: cigar
    :param int flag: flag
    :param int mapping_quality: mapping quality
    :param str next_reference_name: reference name for the paired end alignment mapped
    :param int next_reference_start: reference start of the paired end alignment
    :param bool is_unmapped: whether mapped or not
    """
    alignment = pysam.AlignedSegment(header)
    alignment.reference_name = reference_name
    alignment.query_name = query_name
    alignment.query_sequence = query_sequence
    alignment.reference_start = reference_start
    alignment.cigar = cigar
    alignment.flag = flag
    alignment.mapping_quality = mapping_quality
    alignment.is_unmapped = is_unmapped
    if next_reference_name is not None and next_reference_start is not None and next_reference_start > 0:
        alignment.next_reference_name = next_reference_name
        alignment.next_reference_start = next_reference_start
    return alignment


class TestParseMapTax(unittest.TestCase):

    @patch('builtins.open', new_callable=MagicMock)
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    @patch('subprocess.Popen.wait', new_callable=MagicMock, return_value=0)
    @patch('pysam.sort', return_value=None)
    def run_parse_map_tax(self, mock_pysam_sort, mock_popen, mock_exists, mock_gzip_open, mock_open):
        with (mock.patch('pathlib.Path', return_value=pathlib.Path("/fakepath")),
              mock.patch('motus.motus.MotusParameters.set_read_files') as mock_set_read_files,
              mock.patch('motus.motus.MotusParameters.set_alignment_file') as mock_set_alignment_file,
              mock.patch(
                  'motus.motus.MotusParameters.set_minimal_alignment_length') as mock_set_minimal_alignment_length,
              mock.patch('motus.motus.MotusParameters.set_threads') as mock_set_threads,
              mock.patch('motus.motus.MotusParameters._temp_alignment_file', return_value="temporary.bam"),
              mock.patch('motus.motus.MotusParameters.get_read_files', return_value=[
                  (pathlib.Path("forward_file1"), '/1'),
                  (pathlib.Path("reverse_file1"), '/2'),
                  (pathlib.Path("forward_file2"), '/1'),
                  (pathlib.Path("reverse_file2"), '/2'),
                  (pathlib.Path("unpaired_file"), '/S')
              ]),
              mock.patch('motus.motus.MotusParameters.get_alignment_file',
                         return_value="aligment.bam") as mock_get_alignment_file,
              mock.patch('pysam.AlignmentFile') as pysam_bam,
              mock.patch('sys.stdout', new=StringIO())):
            header = mock_bam_header([('chr1', 100)])  # mock a 100 bp chr1 contig
            in_alignment = mock_alignment(
                header=header,
                reference_name='chr1',
                query_name='aln1',
                query_sequence='ACTGAGAGACGAGAGTT',
                reference_start=10,
                cigar=[(0, 17)],
                flag=0,
                mapping_quality=30,
                is_unmapped=True,
            )
            mock_in_bam = PysamFakeBam(header,
                                       [in_alignment])  # the mock in bam iterator will return our mock alignment
            pysam_bam.return_value = mock_in_bam
            mock_open.return_value.__enter__.return_value = MagicMock(readline=MagicMock(side_effect=[
                "version: 4.0", "date: 2024-01-01"
            ]))

            parse_map_tax()

            mock_set_read_files.assert_called_once()
            mock_get_alignment_file.assert_called_once()
            mock_set_alignment_file.assert_called_once_with(pathlib.Path("output.txt"), required_to_exist=False)
            mock_set_minimal_alignment_length.assert_called_once_with(75)  # Default length is 75
            mock_set_threads.assert_called_once_with(1)  # Default thread count is 1

    def test_default_values(self):
        sys.argv = ["motus", "map_tax", "-f", "forward.fastq", "-o", "output.txt"]
        self.run_parse_map_tax()

    def test_multiple_input_files(self):
        sys.argv = ["motus", "map_tax", "-f", "forward1.fastq", "forward2.fastq", "forward3.fastq", "-r",
                    "reverse1.fastq", "reverse2.fastq", "reverse3.fastq", "-s", "unpaired1.fastq",
                    "unpaired2.fastq", "-o", "output.txt"]
        self.run_parse_map_tax()

    def test_missing_output_argument(self):
        # simulate command-line args without the required output file (-o)
        sys.argv = ["motus", "map_tax", "-f", "forward.fastq"]

        with self.assertRaises(SystemExit):
            sys.stderr = io.StringIO()
            parse_map_tax()
            sys.stderr = sys.__stderr__
        self.assertIn('error: the following arguments are required: -o', sys.stderr.getvalue().strip())

    @patch('builtins.open', new_callable=MagicMock)
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_forward_reverse_unpaired_input(self, mock_exists, mock_gzip_open, mock_open):
        sys.argv = ["motus", "map_tax", "-f", "forward.fastq", "-r", "reverse.fastq", "-s", "unpaired.fastq",
                    "-o", "output.txt"]

        with mock.patch('pathlib.Path', return_value=pathlib.Path("/fakepath")), \
                mock.patch('motus.motus.MotusParameters.set_read_files') as mock_set_read_files, \
                mock.patch('motus.motus.MotusParameters.set_alignment_file') as mock_set_alignment_file, \
                mock.patch('motus.motus.MotusParameters.set_minimal_alignment_length'), \
                mock.patch('motus.motus.MotusParameters.set_threads'), \
                mock.patch('motus.motus.map_tax'):
            parse_map_tax()

            # ensure that the correct file paths are passed to set_read_files
            expected_forward = [pathlib.Path("forward.fastq")]
            expected_reverse = [pathlib.Path("reverse.fastq")]
            expected_unpaired = [pathlib.Path("unpaired.fastq")]
            mock_set_read_files.assert_called_once_with(expected_forward, expected_reverse, expected_unpaired,
                                                        check_files=True)

    @patch('builtins.open', new_callable=MagicMock)
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_custom_min_length_threads_verbosity(self, mock_exists, mock_gzip_open, mock_open):
        # test with custom alignment length, thread count, and verbosity level
        sys.argv = ["motus", "map_tax", "-f", "forward.fastq", "-o", "output.txt", "-l", "100", "-t", "4", "-v", "2"]

        with mock.patch('pathlib.Path', return_value=pathlib.Path("/fakepath")), \
                mock.patch('motus.motus.MotusParameters.set_read_files'), \
                mock.patch('motus.motus.MotusParameters.set_alignment_file'), \
                mock.patch(
                    'motus.motus.MotusParameters.set_minimal_alignment_length') as mock_set_minimal_alignment_length, \
                mock.patch('motus.motus.MotusParameters.set_threads') as mock_set_threads, \
                mock.patch('motus.motus.map_tax'):
            mock_open.return_value.__enter__.return_value = MagicMock(readline=MagicMock(side_effect=[
                "version: 4.0", "date: 2024-01-01"
            ]))
            parse_map_tax()

            mock_set_minimal_alignment_length.assert_called_once_with(100)  # Custom alignment length
            mock_set_threads.assert_called_once_with(4)  # Custom thread count

    @patch('builtins.open', new_callable=MagicMock)
    @patch('gzip.open', new_callable=MagicMock)
    @patch('pathlib.Path.exists', new_callable=MagicMock)
    def test_no_input_files_provided(self, mock_exists, mock_gzip_open, mock_open):
        # test for case where no input files are provided, expect an error
        sys.argv = ["motus", "map_tax", "-o", "output.txt"]
        with self.assertLogs('root', level='ERROR') as log_capture:
            with mock.patch('pathlib.Path', return_value=pathlib.Path("/fakepath")), \
                    self.assertRaises(SystemExit):
                parse_map_tax()
                self.assertEqual(log_capture.output, ['ERROR: No input files defined with -f -r or -s. Quitting ...'])


if __name__ == '__main__':
    unittest.main()
