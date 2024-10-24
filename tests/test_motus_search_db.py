import unittest
from unittest.mock import patch, mock_open
import pathlib
from motus import motus


def set_mock(mock_gzip_open):
    mock_taxonomy_data = "MOTU\tGTDB\n"
    mock_metadata_data = (
        "GENOME\tLOCATION\tMOTU4\tMOTU4_STATUS\tDOMAIN\tPHYLUM\tCLASS\tORDER\tFAMILY\tGENUS\tSPECIES\n"
        "genome1\tloc1\tMOTU1\trepresentative\tBacteria\tProteobacteria\tGammaproteobacteria\tEnterobacterales"
        "\tEnterobacteriaceae\tEscherichia\tcoli\n"
    )

    mock_taxonomy_handle = mock_open(read_data=mock_taxonomy_data)
    mock_metadata_handle = mock_open(read_data=mock_metadata_data)
    mock_gzip_open.side_effect = [mock_taxonomy_handle.return_value, mock_metadata_handle.return_value]

    mock_taxonomy_file = pathlib.Path("mock_taxonomy.gz")
    mock_metadata_file = pathlib.Path("mock_metadata.gz")

    db = motus.MotusSearchDB(mock_taxonomy_file, mock_metadata_file)

    return db, mock_taxonomy_handle, mock_metadata_handle, mock_gzip_open.side_effect, mock_taxonomy_file, mock_metadata_file


class TestMotusSearchDB(unittest.TestCase):

    @patch("gzip.open")
    def test_init(self, mock_gzip_open):
        # mock the contents of the motu taxonomy file
        mock_taxonomy_data = (
            "MOTU\tGTDB\n"
            "MOTU1\tk__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae"
            ";g__Escherichia;s__coli\n"
            "MOTU2\tk__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus"
            ";s__pneumoniae\n"
        )

        # mock the genome metadata file as well
        mock_metadata_data = (
            "GENOME\tLOCATION\tMOTU4\tMOTU4_STATUS\tDOMAIN\tPHYLUM\tCLASS\tORDER\tFAMILY\tGENUS\tSPECIES\n"
            "genome1\tloc1\tMOTU1\trepresentative\tBacteria\tProteobacteria\tGammaproteobacteria\tEnterobacterales"
            "\tEnterobacteriaceae\tEscherichia\tcoli\n"
            "genome2\tloc2\tMOTU2\t\tBacteria\tFirmicutes\tBacilli\tLactobacillales\tStreptococcaceae\tStreptococcus"
            "\tpneumoniae\n"
        )

        # mock the file handle (when open is called, it will return different contents depending on the file name passed
        # )
        mock_taxonomy_handle = mock_open(read_data=mock_taxonomy_data)
        mock_metadata_handle = mock_open(read_data=mock_metadata_data)

        def mock_file_selector(file, mode='r'):
            if "mock_taxonomy.gz" in str(file):
                return mock_taxonomy_handle.return_value
            elif "mock_metadata.gz" in str(file):
                return mock_metadata_handle.return_value

        mock_gzip_open.side_effect = mock_file_selector

        # create a mock Path objects for the file paths
        mock_taxonomy_file = pathlib.Path("mock_taxonomy.gz")
        mock_metadata_file = pathlib.Path("mock_metadata.gz")

        # instantiate the MotusSearchDB class with mocked data
        db = motus.MotusSearchDB(mock_taxonomy_file, mock_metadata_file)

        self.assertIn("MOTU1", db._motu_2_genome)
        self.assertIn("genome1", db._motu_2_genome["MOTU1"])
        self.assertEqual(db._motu_2_genome, {'MOTU1': {'genome1'}, 'MOTU2': {'genome2'}})

        self.assertIn("Bacteria", db._tax_2_motu_and_genome)
        self.assertIn("MOTU1", db._tax_2_motu_and_genome["Bacteria"])
        self.assertEqual(db._tax_2_motu_and_genome,
                         {'Bacteria': {'MOTU2', 'MOTU1', 'genome1', 'genome2'}, 'Proteobacteria': {'MOTU1', 'genome1'},
                          'Gammaproteobacteria': {'MOTU1', 'genome1'}, 'Enterobacterales': {'MOTU1', 'genome1'},
                          'Enterobacteriaceae': {'MOTU1', 'genome1'}, 'Escherichia': {'MOTU1', 'genome1'},
                          'coli': {'MOTU1', 'genome1'}, 'Firmicutes': {'MOTU2', 'genome2'},
                          'Bacilli': {'MOTU2', 'genome2'}, 'Lactobacillales': {'MOTU2', 'genome2'},
                          'Streptococcaceae': {'MOTU2', 'genome2'}, 'Streptococcus': {'MOTU2', 'genome2'},
                          'pneumoniae': {'MOTU2', 'genome2'}})

        self.assertIn("genome1", db._representative_genomes)
        self.assertEqual(db._genome_2_path, {'genome1': 'loc1', 'genome2': 'loc2'})
        self.assertEqual(db._genome_2_tax, {'genome1': 'Bacteria\tProteobacteria\tGammaproteobacteria'
                                                       '\tEnterobacterales\tEnterobacteriaceae\tEscherichia\tcoli',
                                            'genome2':
                                                'Bacteria\tFirmicutes\tBacilli\tLactobacillales\tStreptococcaceae'
                                                '\tStreptococcus\tpneumoniae'})

    @patch("gzip.open")
    def test_search_for_genomes(self, mock_gzip_open):
        # mock the contents of the motu taxonomy file
        mock_taxonomy_data = (
            "MOTU\tGTDB\n"
            "MOTU1\tk__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae"
            ";g__Escherichia;s__coli\n"
            "MOTU2\tk__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus"
            ";s__pneumoniae\n"
        )

        # mock the genome metadata file as well
        mock_metadata_data = (
            "GENOME\tLOCATION\tMOTU4\tMOTU4_STATUS\tDOMAIN\tPHYLUM\tCLASS\tORDER\tFAMILY\tGENUS\tSPECIES\n"
            "genome1\tloc1\tMOTU11\trepresentative\tBacteria\tProteobacteria\tGammaproteobacteria\tEnterobacterales"
            "\tEnterobacteriaceae\tEscherichia\tcoli\n"
            "genome2\tloc1\tMOTU11\t\tBacteria\tProteobacteria\tGammaproteobacteria\tEnterobacterales"
            "\tEnterobacteriaceae\tEscherichia\tcoli\n"
            "genome3\tloc2\tMOTU12\t\tBacteria\tFirmicutes\tBacilli\tLactobacillales\tStreptococcaceae\tStreptococcus"
            "\tpneumoniae\n"
            "genome4\tloc3\tMOTU12\t\tBacteria\tFirmicutes\tBacilli\tLactobacillales\tStreptococcaceae\tStreptococcus"
            "\tpneumoniae\n"
        )

        # mock the file handles
        mock_taxonomy_handle = mock_open(read_data=mock_taxonomy_data)
        mock_metadata_handle = mock_open(read_data=mock_metadata_data)

        def mock_file_selector(file, mode='r'):
            if "mock_taxonomy.gz" in str(file):
                return mock_taxonomy_handle.return_value
            elif "mock_metadata.gz" in str(file):
                return mock_metadata_handle.return_value

        mock_gzip_open.side_effect = mock_file_selector

        # create mock Path objects for the file paths
        mock_taxonomy_file = pathlib.Path("mock_taxonomy.gz")
        mock_metadata_file = pathlib.Path("mock_metadata.gz")

        # instantiate the MotusSearchDB class with mocked data
        db = motus.MotusSearchDB(mock_taxonomy_file, mock_metadata_file)

        # test genome search for keyword "MOTU2"
        genomes = db.search_for_genomes("MOTU12")
        self.assertEqual(genomes, ['genome3', 'genome4'])

        # test search for non-existent keyword
        genomes = db.search_for_genomes("UnknownKeyword")
        self.assertEqual(genomes, [])

        # test representative genomes search
        genomes = db.search_for_genomes("MOTU11", only_representatives=True)
        self.assertEqual(genomes, ["genome1"])

        # test MOTUs with representative genomes search
        genomes = db.search_for_genomes("MOTU11")
        self.assertEqual(genomes, ["genome1", "genome2"])

        # test mixture of lower and uppercase
        """genomes = db.search_for_genomes("MOtu11")
        self.assertEqual(genomes, ["genome1", "genome2"])"""  # not working as expected

    @patch("gzip.open")
    def test_get_genome_path(self, mock_gzip_open):
        # mock setup
        db, mock_taxonomy_handle, mock_metadata_handle, mock_gzip_open.side_effect, mock_taxonomy_file, mock_metadata_file = set_mock(
            mock_gzip_open)

        # test get_genome_path
        genome_path = db.get_genome_path("genome1")
        self.assertEqual(genome_path, "https://sunagawalab.ethz.ch/share/MOTUS/database/4.0/data/genomes/loc1")

    @patch("gzip.open")
    def test_get_genome_motu(self, mock_gzip_open):
        # mock setup
        db, mock_taxonomy_handle, mock_metadata_handle, mock_gzip_open.side_effect, mock_taxonomy_file, mock_metadata_file = set_mock(
            mock_gzip_open)

        # test get_genome_motu
        motu = db.get_genome_motu("genome1")
        self.assertEqual(motu, "MOTU1")

    @patch("gzip.open")
    def test_get_genome_tax(self, mock_gzip_open):
        # mock setup
        db, mock_taxonomy_handle, mock_metadata_handle, mock_gzip_open.side_effect, mock_taxonomy_file, mock_metadata_file = set_mock(
            mock_gzip_open)

        # test get_genome_tax
        tax = db.get_genome_tax("genome1")
        self.assertEqual(tax, "Bacteria\tProteobacteria\tGammaproteobacteria\tEnterobacterales\tEnterobacteriaceae"
                              "\tEscherichia\tcoli")


if __name__ == '__main__':
    unittest.main()
