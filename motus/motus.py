#!/usr/bin/env python


# ============================================================================ #
# motus - a tool for marker gene-based OTU (mOTU) profiling of metagenomes
#
# Authors: Hans-Joachim Ruscheweyh (hansr@ethz.ch),
#          Lilith Feer,
#          Marija Dmitrijeva,
#          Kang Li,
#          Florian Ruscheweyh,
#          Daniel R. Mende
#          Georg Zeller,
#          Shinichi Sunagawa
#
# Type "motus" for usage help
#
# Copyright (c) ${2024} ${SunagawaLab}.
#
# This file is part of ${projectname}
# (see ${https://motus-tool.org/}).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# ============================================================================ #




import statistics
import pysam
import Bio.SeqIO.FastaIO as FastaIO
import logging
import pathlib
import csv
import subprocess
import gzip
import collections
import random
import argparse
import sys
from typing import List, Dict, Set, Tuple, Generator, TextIO
import urllib.request
import tarfile
import shutil
from enum import Enum
import mutils
from mentities import MOTUS_PARAMETERS
from mentities import MOTUS_DB
import mentities


__author__ = ('Hans-Joachim Ruscheweyh (hansr@ethz.ch), '
              'Lilith Feer, '
              'Marija Dmitrijeva, '
              'Kang Li, '
              'Florian Ruscheweyh'
              'Daniel Mende, '
              'Georg Zeller, '
              'Shinichi Sunagawa')
__version__ = mutils.MOTUS_VERSION
__date__ = '27 June 2024'
__license__ = "GPL - v3"
__maintainer__ = "Hans-Joachim Ruscheweyh"


"""
Terminology

markergeneheader = mgh = an instance of an markergene
markergene = mg = one of the 10 mOTUs markergenes
markergenecluster = mgc = a set of mgh that come from the same mOTU and markergene
motu = Top level unit, species level cluster
"""






Mgc_values = collections.namedtuple("Mgc_values", "insert_raw insert_norm insert_scaled base_raw base_norm")


class TaxonomicRank(Enum):
    DOMAIN = 'domain'
    PHYLUM = 'phylum'
    CLASS = 'class'
    ORDER = 'order'
    FAMILY = 'family'
    GENUS = 'genus'
    SPECIES = 'species'



class SinglemOTUsFile:
    _count_mode = None
    _min_mgcs = None
    _full_version = None
    _taxonomy = None
    _taxonomy_level = None # TODO remove. we dont aggregate in motus
    _counts_aggregated_by_taxonomy = None # TODO remove. we dont aggregate in motus
    _motus_with_abundance = None
    _samplename_2_motus_2_counts = None
    _samplename_2_motus_2_relab = None
    _motu_2_taxonomy = None # TODO not sure this will be needed





class MotusFile:
    _count_mode = None
    _min_mgcs = None
    _full_version = None
    _taxonomy = None
    _taxonomy_level = None # TODO remove. we dont aggregate in motus
    _counts_aggregated_by_taxonomy = None # TODO remove. we dont aggregate in motus
    _motus_with_abundance = None
    _samplename_2_motus_2_counts = None
    _samplename_2_motus_2_relab = None
    _motu_2_taxonomy = None # TODO not sure this will be needed

    def has_counts(self):
        if self._samplename_2_motus_2_counts:
            return True
        else:
            return False

    def has_motus_with_abundance(self):
        if len(self._motus_with_abundance) == 0:
            return False
        else:
            return True

    def set_mOTU_counts(self, samplename_2_motus_2_counts: Dict[str, float], count_type):
        self._motus_with_abundance = set()
        self._samplename_2_motus_2_counts = {}
        for samplename, motu_2_counts in samplename_2_motus_2_counts.items():
            self._samplename_2_motus_2_counts[samplename] = {}
            for motu, count in motu_2_counts.items():
                self._motus_with_abundance.add(motu)
                if count_type == 'int':
                    count = round(count)
                    self._samplename_2_motus_2_counts[samplename][motu] = count
                else:
                    self._samplename_2_motus_2_counts[samplename][motu] = count
        self._motus_with_abundance = sorted(list(self._motus_with_abundance))

        self._samplename_2_motus_2_relab = {}

        for samplename, motus_2_counts in samplename_2_motus_2_counts.items():
            self._samplename_2_motus_2_relab[samplename] = {}
            tot_abundance = sum(motus_2_counts.values())
            for motu, value in motus_2_counts.items():
                self._samplename_2_motus_2_relab[samplename][motu] = value / tot_abundance #'{number:.{digits}f}'.format(number=value / tot_abundance, digits=8)

    def set_sample_name(self, sample_name):
        self._samplename = sample_name

    def set_count_mode(self, count_mode):
        self._count_mode = count_mode

    def set_min_mgcs(self, min_mgcs):
        self._min_mgcs = min_mgcs

    def set_full_version(self, full_version):
        self._full_version = full_version

    def get_mOTUs_file_header(self, relabundance=False) -> str:

        relab = relabundance
        min_mgcs = self._min_mgcs
        count_mode = self._count_mode
        full_version = self._full_version
        if not min_mgcs:
            logging.error('min_mgcs parameter not set. Can\'t create header. Quitting...')
        if not full_version:
            logging.error('full_version parameter not set. Can\'t create header. Quitting...')
        if not count_mode:
            logging.error('count_mode parameter not set. Can\'t create header. Quitting...')

        if relab:
            relab = 'relative_abundance'
        else:
            relab = 'counts'
        header = f'#{full_version}\treport_mode={relab}\tcount_mode={count_mode}\tmin_mgcs={min_mgcs}'
        return header

    def merge_profiles(self, motus_files) -> None: # Typing --> This is a List[MotusFile]
        """ merge multiple profiles into one
        Will also check if profiles are compatible

        Params:
            motus_files: A list with MotusFile objects, each containing
                        one or multiple motus profiles

        Returns:
            None, will update the MotusFile object
        """

        if not motus_files or len(motus_files) == 0:
            logging.error(f'No MotusFiles found to merge. Quitting ...')
            mutils.shutdown(1)

        # check if versions are compatible
        versions = set([mf._full_version for mf in motus_files])
        if len(versions) != 1:
            logging.error(f'Incompatible versions in profiles that should be merged. versions = {versions}')
            mutils.shutdown(1)

        # check if count_mode is compatible
        count_modes = set([mf._count_mode for mf in motus_files])
        if len(count_modes) != 1:
            logging.error(f'Incompatible count modes in profiles that should be merged. count modes = {count_modes}')
            mutils.shutdown(1)

        # check if num_mgcs is compatible
        min_mgcs = set([mf._min_mgcs for mf in motus_files])
        if len(min_mgcs) != 1:
            logging.error(f'Incompatible min mgcs in profiles that should be merged. min mgcs = {min_mgcs}')
            mutils.shutdown(1)

        # check if all or none have counts --> either all or none have to have counts1
        # if none have counts --> report rel abundances, if all have counts --> report counts

        total_motus_files = len(motus_files)
        total_motus_files_w_counts = len([mf for mf in motus_files if mf._samplename_2_motus_2_counts])
        total_motus_files_w_relab = len([mf for mf in motus_files if mf._samplename_2_motus_2_relab])
        report_counts = False
        report_relab = False
        if total_motus_files_w_counts == 0:
            report_relab = True
        elif total_motus_files_w_counts != total_motus_files:
            logging.info('Profile files are mixed. Some are reported as counts, some as relative abundances. Quitting ...')
            mutils.shutdown(1)
        else:
            report_counts = True
            report_relab = True

        # check if all sample names are distinct
        samplenames = collections.Counter()
        if report_relab:
            for mf in motus_files:
                for samplename in mf._samplename_2_motus_2_relab.keys():
                    samplenames[samplename] += 1
        for samplename, samplename_count in samplenames.items():
            if samplename_count != 1:
                logging.error(f'Samplename duplicated: {samplename}. Quitting ...')
                mutils.shutdown(1)
        present_motus = set()
        merged_sample_2_motus_counts = None
        merged_sample_2_motus_relab = {}
        if report_counts:
            merged_sample_2_motus_counts = {}
            for mf in motus_files:
                for sample, motu_2_counts in mf._samplename_2_motus_2_counts.items():
                    merged_sample_2_motus_counts[sample] = {}
                    for motu, count in motu_2_counts.items():
                        present_motus.add(motu)
                        merged_sample_2_motus_counts[sample][motu] = count
        for mf in motus_files:
            for sample, motu_2_relab in mf._samplename_2_motus_2_relab.items():
                merged_sample_2_motus_relab[sample] = {}
                for motu, relab in motu_2_relab.items():
                    present_motus.add(motu)
                    merged_sample_2_motus_relab[sample][motu] = relab
        present_motus = sorted(list(present_motus))

        self._motus_with_abundance = present_motus
        self._samplename_2_motus_2_relab = merged_sample_2_motus_relab
        self._samplename_2_motus_2_counts = merged_sample_2_motus_counts
        self._count_mode = count_modes.pop()
        self._full_version = versions.pop()
        self._min_mgcs = min_mgcs.pop()

    def read_mOTUs_file(self, mOTUs_file: pathlib.Path) -> None:
        """
        Read the mOTUs_file into this object.
        """
        with open(mOTUs_file) as handle:
            header = handle.readline().strip()
            splits = header.split('\t')

            if len(splits) == 4:
                [full_version, report_mode, count_mode, min_mgcs] = splits
                taxonomy = None
                aggregated = False
                level = None
            elif len(splits) == 7:
                [full_version, report_mode, count_mode, min_mgcs, taxonomy, aggregated, level] = splits
            else:
                logging.error('The header of this mOTUs file looks malformed. Expected 7 columns. Please check. Quitting ...')
                logging.error(f'{header}')
                mutils.shutdown(1)

            if not full_version.startswith('#TOOL'):
                logging.error('The header of this mOTUs file looks malformed. #TOOL token is missing. Please check. Quitting ...')
                logging.error('Malformed header:')
                logging.error(f'{header}')
                mutils.shutdown(1)
            full_version = full_version[1:]

            if 'report_mode' not in report_mode:
                logging.error('The header of this mOTUs file looks malformed. report_mode token missing. Please check. Quitting ...')
                logging.error('Malformed header:')
                logging.error(f'{header}')
                mutils.shutdown(1)
            report_mode = report_mode.split('=')[1]

            if 'count_mode' not in count_mode:
                logging.error('The header of this mOTUs file looks malformed. count_mode token missing. Please check. Quitting ...')
                logging.error('Malformed header:')
                logging.error(f'{header}')
                mutils.shutdown(1)
            count_mode = count_mode.split('=')[1]

            if 'min_mgcs' not in min_mgcs:
                logging.error('The header of this mOTUs file looks malformed. min_mgcs token missing. Please check. Quitting ...')
                logging.error('Malformed header:')
                logging.error(f'{header}')
                mutils.shutdown(1)
            min_mgcs = int(min_mgcs.split('=')[1])

            if taxonomy:
                if 'taxonomy' not in taxonomy:
                    logging.error('The header of this mOTUs file looks malformed. taxonomy token missing. Please check. Quitting ...')
                    logging.error('Malformed header:')
                    logging.error(f'{header}')
                    mutils.shutdown(1)
                taxonomy = taxonomy.split('=')[1]

            if level:
                if 'level' not in level:
                    logging.error('The header of this mOTUs file looks malformed. level token missing. Please check. Quitting ...')
                    logging.error('Malformed header:')
                    logging.error(f'{header}')
                    mutils.shutdown(1)
                level = level.split('=')[1]
            if aggregated:
                if 'aggregated' not in aggregated:
                    logging.error('The header of this mOTUs file looks malformed. aggregated token missing. Please check. Quitting ...')
                    logging.error('Malformed header:')
                    logging.error(f'{header}')
                    mutils.shutdown(1)
                aggregated = aggregated.split('=')[1]
                if bool(aggregated):
                    logging.error(f'The mOTUs profile ({mOTUs_file}) has values aggregated at non-mOTU level ({level}). This table is an endproduct and cannot be used in mOTUs anymore. Quitting ...')
                    mutils.shutdown(1)

            self._count_mode = count_mode
            self._min_mgcs = min_mgcs
            self._full_version = full_version
            self._taxonomy = taxonomy
            self._taxonomy_level = level
            self._counts_aggregated_by_taxonomy = None
            self._motus_with_abundance = set()
            if taxonomy:
                self._motu_2_taxonomy = {}
                self._taxonomy = taxonomy
                self._taxonomy_level = level
                self._counts_aggregated_by_taxonomy = bool(aggregated)
            else:
                self._motu_2_taxonomy = None
                self._taxonomy = None
                self._taxonomy_level = None
                self._counts_aggregated_by_taxonomy = None

            sample_2_motu_2_cnts = {}
            for entry in csv.DictReader(handle, delimiter='\t'):
                motu = entry.pop('MOTU')
                self._motus_with_abundance.add(motu)
                if taxonomy:
                    taxstring = entry.pop('TAXONOMY')
                    self._motu_2_taxonomy[motu] = taxstring
                for sample, counts in entry.items():
                    if sample not in sample_2_motu_2_cnts:
                        sample_2_motu_2_cnts[sample] = {}
                    sample_2_motu_2_cnts[sample][motu] = float(counts)

            if report_mode == 'counts':
                self._samplename_2_motus_2_counts = sample_2_motu_2_cnts
                self._samplename_2_motus_2_relab = {}

                for samplename, motus_2_counts in sample_2_motu_2_cnts.items():
                    self._samplename_2_motus_2_relab[samplename] = {}
                    tot_abundance = sum(motus_2_counts.values())
                    for motu, value in motus_2_counts.items():
                        if int(value) == 0:
                            self._samplename_2_motus_2_relab[samplename][motu] = 0.0
                        else:
                            self._samplename_2_motus_2_relab[samplename][motu] = value / tot_abundance
            elif report_mode == 'relative_abundance':
                self._samplename_2_motus_2_counts = None
                self._samplename_2_motus_2_relab = sample_2_motu_2_cnts
            else:
                logging.error(f'Report mode can only be counts or relative_abundance but is {report_mode}. Quitting ...')
                mutils.shutdown(1)

            self._motus_with_abundance = sorted(list(self._motus_with_abundance))

    def write_mOTUs_file(self, filename, relabundance) -> None:
        # currently implemented without taxonomy
        with open(filename, 'w') as handle:
            sample_2_motu_2_report_vals = self._samplename_2_motus_2_counts
            if relabundance:
                sample_2_motu_2_report_vals = self._samplename_2_motus_2_relab
            handle.write(self.get_mOTUs_file_header(relabundance=relabundance) + '\n')
            samplenames = sorted(list(sample_2_motu_2_report_vals.keys()))
            tmp = '\t'.join(samplenames)
            handle.write(f'MOTU\t{tmp}\n')

            for motu in self._motus_with_abundance:
                tmp = [motu]
                for samplename in samplenames:
                    if relabundance:
                        abundance = '{number:.{digits}f}'.format(number=sample_2_motu_2_report_vals[samplename].get(motu, 0), digits=8)
                    else:
                        abundance = sample_2_motu_2_report_vals[samplename].get(motu, 0)
                        if not 'NORM' in self._count_mode:
                            abundance = int(abundance)
                    tmp.append(str(abundance))
                tmp = '\t'.join(tmp)
                handle.write(f'{tmp}\n')



class MotusSearchDB:
    _motu_2_genome = collections.defaultdict(set)
    _tax_2_motu_and_genome = collections.defaultdict(set)
    _genome_2_motu = {}
    _genome_2_path = {}
    _genome_2_tax = {}
    _representative_genomes = set()

    def __init__(self, motu_taxonomy_file: pathlib.Path, genome_metadata_file: pathlib.Path) -> None:

        logging.info('Initialising the mOTUs search database.')
        with gzip.open(motu_taxonomy_file, 'rt') as handle:
            handle.readline()
            for line in handle:
                [motu, gtdb] = line.strip().split('\t')
                [domain, phylum, class_rank, order, family, genus, species] = [x.split('__')[1] for x in gtdb.split(';')]
                self._tax_2_motu_and_genome[domain].add(motu)
                self._tax_2_motu_and_genome[phylum].add(motu)
                self._tax_2_motu_and_genome[class_rank].add(motu)
                self._tax_2_motu_and_genome[order].add(motu)
                self._tax_2_motu_and_genome[family].add(motu)
                self._tax_2_motu_and_genome[genus].add(motu)
                self._tax_2_motu_and_genome[species].add(motu)
        with gzip.open(genome_metadata_file, 'rt') as handle:
            for entry in csv.DictReader(handle, delimiter='\t'):
                genome = entry['GENOME']
                location = entry['LOCATION']
                motu = entry['MOTU4']
                self._genome_2_motu[genome] = motu
                self._motu_2_genome[motu].add(genome)
                self._genome_2_path[genome] = location
                if 'representative' in entry['MOTU4_STATUS']:
                    self._representative_genomes.add(genome)
                [domain, phylum, class_rank, order, family, genus, species] = [entry['DOMAIN'], entry['PHYLUM'], entry['CLASS'], entry['ORDER'], entry['FAMILY'], entry['GENUS'], entry['SPECIES']]
                self._genome_2_tax[genome] = '\t'.join([domain, phylum, class_rank, order, family, genus, species])

                self._tax_2_motu_and_genome[domain].add(genome)
                self._tax_2_motu_and_genome[phylum].add(genome)
                self._tax_2_motu_and_genome[class_rank].add(genome)
                self._tax_2_motu_and_genome[order].add(genome)
                self._tax_2_motu_and_genome[family].add(genome)
                self._tax_2_motu_and_genome[genus].add(genome)
                self._tax_2_motu_and_genome[species].add(genome)

        logging.info(f'Finished initialising the mOTUs search database. Found {len(self._motu_2_genome)} mOTUs, {len(self._genome_2_path)} genomes and {len(self._tax_2_motu_and_genome)} taxonomy search words.')

    def search_for_genomes(self, keyword, only_representatives = False) -> List[str]:
        """
        Search for genomes in the mOTUs database using an exact keyword
        Search is exact but tolerates upper/lowercase differences

        Params:
            keyword: a string of at least one word, has to match
                    exactly a mOTU, genome name or a taxon in GTDB
            only_representative: If True, only return representative
                    genomes
        Returns:
            a set with all genomes that have been found
        """
        genomes = set()
        for genome in self._motu_2_genome.get(keyword, []):
            genomes.add(genome)
        if keyword in self._genome_2_path:
            genomes.add(self._genome_2_path[keyword])
        for genome_or_motu in self._tax_2_motu_and_genome.get(keyword, []):
            if genome_or_motu in self._genome_2_path:
                genomes.add(genome_or_motu)
            for genome in self._motu_2_genome.get(genome_or_motu, []):
                genomes.add(genome)
        report_genomes = set()
        if only_representatives:
            for genome in genomes:
                if genome in self._representative_genomes:
                    report_genomes.add(genome)
        else:
            report_genomes = genomes
        report_genomes = sorted(list(report_genomes))
        return report_genomes
    def get_genome_path(self, genome:str):
        p = mutils.MOTUS_GENOME_REMOTE_PREFIX + self._genome_2_path[genome]
        return p
    def get_genome_motu(self, genome: str):
        return self._genome_2_motu[genome]
    def get_genome_tax(self, genome: str):
        return self._genome_2_tax[genome]

def map_tax() -> None:
    """
    Alignment of read files against the mOTUs MGDB:

    Takes a list of forward/reverse/unpaired read files and aligns them against the mOTUs database.
    Alignments will be filtered by 97% identity and the defined minimal alignment length.
    The resulting alignments will be stored in the sorted BAM file which is either specified as
    a parameter or as a temporary file.

    Returns:
        None

    """

    logging.info('Starting mOTUs - map_tax routine - Alignment against the mOTUs database ... ')
    min_perc_id: float = 97.0
    threads: int = MOTUS_PARAMETERS.get_threads()
    minlength: int = MOTUS_PARAMETERS.get_minimal_alignment_length()
    temp_bam_file = MOTUS_PARAMETERS.get_temporary_alignment_file()
    temp_bam_file_handle = None

    total_reads = 0
    total_mapped_reads = 0

    for readsfile, orientation in MOTUS_PARAMETERS.get_read_files():
        total_reads_this_file: int = 0
        total_mapped_reads_this_file: Set[str] = set()
        logging.info(f'Aligning {readsfile}')
        command: str = f'bwa mem -a -t {threads} {MOTUS_DB.get_bwa_index()} {readsfile}'
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        in_bam_file_handle = pysam.AlignmentFile(process.stdout, 'rb')

        if not temp_bam_file_handle:
            alignmentfile_header = in_bam_file_handle.header.to_dict()
            pg_header = {}
            pg_header['CL'] = 'motus.py map_tax '
            pg_header['PN'] = 'motus.py'
            motus_version = MOTUS_DB.get_full_version()
            pg_header['VN'] = motus_version
            pg_header['ID'] = MOTUS_DB.get_full_sam_id()
            alignmentfile_header['PG'].append(pg_header)
            temp_bam_file_handle = pysam.AlignmentFile(temp_bam_file, "wb", header = alignmentfile_header)

        for record in in_bam_file_handle:
            if record.is_unmapped:
                total_reads_this_file += 1
                continue
            else:
                if not record.is_secondary and not record.is_supplementary:
                    total_reads_this_file += 1
                if MOTUS_DB.is_mg_blocked(record.reference_name):
                    continue
                alnlength: int = sum(record.get_cigar_stats()[0][0:3])
                if alnlength < minlength:
                    continue
                query_covered_bases: int = sum(record.get_cigar_stats()[0][0:2])
                query_length: int = record.infer_read_length()
                mismatches: int = record.get_tag('NM')
                percid: float = (alnlength - mismatches) / float(alnlength) * 100.0
                percid: float = round(percid, 2)
                if min_perc_id > percid:
                    continue
                qcov: float = query_covered_bases / float(query_length)
                record.set_tag('id', percid, 'f')
                record.set_tag('qc', qcov, 'f')
                record.set_tag('al', alnlength, 'i')

                total_mapped_reads_this_file.add(record.qname)
                record.qname = ''.join([record.qname, orientation])
                temp_bam_file_handle.write(record)

        logging.info(f'Finished alignment. Total reads: {total_reads_this_file}, Total aligned reads {len(total_mapped_reads_this_file)}, {round(len(total_mapped_reads_this_file) * 100.0 / total_reads_this_file, 4)}% aligned.')
        total_mapped_reads += len(total_mapped_reads_this_file)
        total_reads += total_reads_this_file
        in_bam_file_handle.close()

    process.stdout.close()
    return_code: int = process.wait()
    if return_code != 0:
        logging.error(f'BWA command failed with return code {return_code}')
        mutils.shutdown(1)
    logging.info(f'Finished all alignments. Total reads: {total_reads}, Total aligned reads {total_mapped_reads}, {round(total_mapped_reads * 100.0 / total_reads, 4)}%')
    temp_bam_file_handle.close()

    logging.info(f'Sorting BAM file')
    pysam.sort('-n', '-m', '1G', '-@', '1', '-o', str(MOTUS_PARAMETERS.get_alignment_file()),  str(MOTUS_PARAMETERS.get_temporary_alignment_file()))
    logging.info(f'Finished sorting BAM file')

    MOTUS_PARAMETERS.delete_temporary_alignment_file()
    logging.info('Finished mOTUs - map_tax routine - Alignment against the mOTUs database ...')
    return None


def _get_orientation_of_aligned_segment_by_name(alignment: pysam.AlignedSegment) -> Tuple[str, str]:
    """
    Receives an alignment, checks the orientation of the alignment
    (forward, reverse, singleton) and return the name of the insert and the
    orientation

    Params:
        alignment: Alignment information of a read against a markergenesequence

    Return:
        A Tuple with (insert_name, orientation)
    """

    splits = alignment.query_name.rsplit('/', 1)
    if len(splits) == 2:
        return splits[0], splits[1]
    else:
        return alignment.query_name, mentities.SIDENTIFIER


class BestAlignment:
    """
    An object to keep track of the best alignment of an insert.
    """

    _mg_2_blocks = None

    def __init__(self):
        """
        Init this class
        """

        self._mg_2_blocks = {}

    def append(self, mg: str, blocks: List[(Tuple[int, int])]) -> None:
        """
        Set the alignment blocks

        Params:
            mg: the markergene an insert has aligned against
            blocks: the alignmentblocks as reported by pysam

        Returns:
            None
        """
        self._mg_2_blocks[mg] = blocks

    def isMultimapper(self) -> bool:
        """
        Checks whether this alignment
        is a multimapper

        Returns:
            True is this alignment is a multimapper
        """

        if len(self._mg_2_blocks) == 1:
            return False
        else:
            return True

    def get_mg_and_blocks(self) -> List[(Tuple[int, int])]:
        """
        Get the alignmentblocks on case this insert is a
        unique mapper. Fail if this insert is a multimapper

        Returns:
            markergeneheader and alignment blocks
        """

        if self.isMultimapper():
            logging.error('This method doesnt work for multi mappers.')
            mutils.shutdown(1)
        for mg, blocks in self._mg_2_blocks.items():
            return mg, blocks

    def get_mgs_and_blocks(self):
        """
        Get markergeneheaders and alignment blocks.
        Will fail if this insert is a unique mapper

        Returns:
            Dictionary with markergeneheaders to alignment blocks
        """

        if not self.isMultimapper():
            logging.error('This method doesnt work for unique mappers.')
            mutils.shutdown(1)
        return self._mg_2_blocks


class InsertCounter:
    _unique_mappers = None
    _multi_mappers = None

    _mg_2_edge_corrected_raw_uniquemapper_insert_counts = {}
    _mg_2_edge_corrected_raw_uniquemapper_base_counts = {}
    _mg_2_edge_corrected_raw_multimapper_insert_counts = {}
    _mg_2_edge_corrected_raw_multimapper_base_counts = {}

    _mg_2_edge_corrected_raw_insert_counts = {}
    _mg_2_edge_corrected_raw_base_counts = {}
    _mg_2_edge_corrected_norm_insert_counts = {}
    _mg_2_edge_corrected_scaled_insert_counts = {}
    _mg_2_edge_corrected_norm_base_counts = {}
    _mg_2_edge_corrected_scaled_base_counts = {}

    def get_mg_insert_raw(self):
        return self._mg_2_edge_corrected_raw_insert_counts

    def get_mg_base_raw(self):
        return self._mg_2_edge_corrected_raw_base_counts

    def get_mg_insert_norm(self):
        return self._mg_2_edge_corrected_norm_insert_counts

    def get_mg_base_norm(self):
        return self._mg_2_edge_corrected_norm_base_counts

    def get_mg_insert_scaled(self):
        return self._mg_2_edge_corrected_scaled_insert_counts

    def __init__(self):
        self._unique_mappers = []
        self._multi_mappers = []

    def appendmapper(self, insert_name:str, bestAlignment: BestAlignment) -> None:
        if bestAlignment.isMultimapper():
            self._multi_mappers.append((insert_name, bestAlignment))
        else:
            self._unique_mappers.append((insert_name, bestAlignment))

    def get_unique_mapper_count(self):
        return len(self._unique_mappers)

    def get_multi_mapper_count(self):
        return len(self._multi_mappers)


    def _get_edge_corrected_raw_uniquemapper_insert_counts_per_mgc(self):
        '''
        Summarize the raw unique mapper counts
        '''
        mgc_2_count = collections.Counter()
        for mg, count in self._mg_2_edge_corrected_raw_uniquemapper_insert_counts.items():
            mgc_2_count[MOTUS_DB.get_mgc_by_mg(mg)] += count
        return mgc_2_count

    def correct_multi_mapper_edges(self, insert_file_writer: TextIO, min_alignment_length: int):
        '''
        Routine to distribute multimapper based on fractional MGC abundance and to correct edges

        1. Collect MGC abundances of unique mappers
        2. Distribute multimappers:
            For each multimapper (only one mg per MGC)
            collect the total weight by summing up unique MGC
            abundances of all aligned MGs in this multimapper.
            Then create a fractional weight for each aligned MG.
            Then distribute abundances based on that weight.
            Give every MG the same weight in case all MGC
            abundances are 0.

        '''

        mg_2_alignments = collections.defaultdict(list)
        multimapper_with_no_prior_mgc_abundance = 0
        mgc_2_uniquemapper_insert_counts = self._get_edge_corrected_raw_uniquemapper_insert_counts_per_mgc()
        for insert_name, bestAlignment in self._multi_mappers:
            mg_2_blocks = bestAlignment.get_mgs_and_blocks()
            if False: # MG weighting
                tot_weight = sum([self._mg_2_edge_corrected_raw_uniquemapper_insert_counts.get(mg, 0.0) for mg in mg_2_blocks.keys()])
                if tot_weight < 1.0:
                    tot_weight = 1.0
                mg_2_weight = {mg: self._mg_2_edge_corrected_raw_uniquemapper_insert_counts.get(mg, 0.0) / tot_weight for mg in mg_2_blocks.keys()}
                if sum(mg_2_weight.values()) == 0.0: # this part is different in the MGC weighing part
                    mg_2_weight2 = {}
                    for mg in mg_2_weight:
                        mg_2_weight2[mg] = 1.0 / len(mg_2_weight)
                    mg_2_weight = mg_2_weight2
                for mg, alignment_blocks in mg_2_blocks.items():
                    if mg_2_weight[mg] != 0.0:
                        mg_2_alignments[mg].append((alignment_blocks, mg_2_weight[mg]))
            else: # MGC weighting
                tot_mgc_weight = sum([mgc_2_uniquemapper_insert_counts.get(MOTUS_DB.get_mgc_by_mg(mg), 0.0) for mg in mg_2_blocks.keys()])
                if tot_mgc_weight < 1.0:
                    tot_mgc_weight = 1.0
                mg_2_mgc_weight = {mg: mgc_2_uniquemapper_insert_counts.get(MOTUS_DB.get_mgc_by_mg(mg), 0.0) / tot_mgc_weight for mg in mg_2_blocks.keys()}
                if sum(mg_2_mgc_weight.values()) == 0.0:
                    multimapper_with_no_prior_mgc_abundance += 1
                    continue


                for mg, alignment_blocks in mg_2_blocks.items():
                    if mg_2_mgc_weight[mg] != 0.0:
                        insert_file_writer.write(f'{insert_name}\t{mg}\t{round(mg_2_mgc_weight[mg], 4):.4f}\n')
                        mg_2_alignments[mg].append((alignment_blocks, mg_2_mgc_weight[mg]))

        # TODO correct the minimum alignment length. It is set differently in the unique mapper routine

        logging.info(f'Processed {len(self._multi_mappers)} multimappers. {multimapper_with_no_prior_mgc_abundance} were discarded, {len(self._multi_mappers) - multimapper_with_no_prior_mgc_abundance} were used.')
        mg_2_edge_corrected_insert_counts, mg_2_edge_corrected_base_counts = self._correct_edges(mg_2_alignments,min_alignment_length)
        self._mg_2_edge_corrected_raw_multimapper_insert_counts = mg_2_edge_corrected_insert_counts
        self._mg_2_edge_corrected_raw_multimapper_base_counts = mg_2_edge_corrected_base_counts

    def _correct_edges(self, mg_2_alignments, min_alignment_length):
        """Corrects for missing alignments towards the gene edges.
        Assuming you have a read which is only partly overlapping with a
        markergene therefor the alignment is too short is filtered. This
        happens mostly at the edges of genes. We can correct for that by
        something we call inverse padding (formerly known as edge correction).
        We remove all aligned bases from the end regions of the gene and then
        extrapolate the abundance based on the median abundance of the rest of
        the gene

        """

        mg_2_trunc_insert_counts = collections.Counter()
        mg_2_untrunc_insert_counts = collections.Counter()
        mg_2_trunc_base_counts = collections.Counter()
        mg_2_untrunc_base_counts = collections.Counter()

        # alignments -> All alignments against a mg. One alignment represents multiple alignment blocks
        for mg, alignments in mg_2_alignments.items():
            first_allowed_base = min_alignment_length + 1
            mg_len = MOTUS_DB.get_length_by_mg(mg)
            last_allowed_base = mg_len - min_alignment_length - 1
            alignments_trunc = []
            for (alignment_blocks, weight) in alignments:
                aligned_bases_untrunc = 0
                aligned_bases_trunc = 0
                for (aln_start, aln_end) in alignment_blocks:
                    if aln_start > aln_end:
                        aln_start_tmp = aln_end
                        aln_end = aln_start
                        aln_start = aln_start_tmp
                    aligned_bases_untrunc += aln_end - aln_start
                    if aln_end < first_allowed_base:
                        continue
                    if aln_start > last_allowed_base:
                        continue
                    if aln_start < first_allowed_base:
                        aln_start = first_allowed_base
                    if aln_end > last_allowed_base:
                        aln_end = last_allowed_base

                    aligned_bases_trunc += aln_end - aln_start
                if aligned_bases_trunc != 0:
                    alignments_trunc.append(aligned_bases_trunc / aligned_bases_untrunc) # a value of 1 means that the alignment was not truncated. a value <1 means that the alignment was truncated
                mg_2_untrunc_base_counts[mg] += aligned_bases_untrunc * weight
                mg_2_trunc_base_counts[mg] += aligned_bases_trunc * weight

            mg_2_untrunc_insert_counts[mg] = len(alignments) * weight
            mg_2_trunc_insert_counts[mg] = len(alignments_trunc) * weight

        mg_2_edge_corrected_insert_counts = collections.Counter()
        mg_2_edge_corrected_base_counts = collections.Counter()

        for mg, trunc_insert_count in mg_2_trunc_insert_counts.items():
            mg_len = MOTUS_DB.get_length_by_mg(mg)
            mg_trunc_len = mg_len - 2 * min_alignment_length
            edge_corrected_insert_count = mg_len * (trunc_insert_count / mg_trunc_len)
            mg_2_edge_corrected_insert_counts[mg] = edge_corrected_insert_count

        for mg, trunc_base_count in mg_2_trunc_base_counts.items():
            mg_len = MOTUS_DB.get_length_by_mg(mg)
            mg_trunc_len = mg_len - 2 * min_alignment_length
            edge_corrected_base_count = mg_len * (trunc_base_count / mg_trunc_len)
            mg_2_edge_corrected_base_counts[mg] = edge_corrected_base_count
        return mg_2_edge_corrected_insert_counts, mg_2_edge_corrected_base_counts

    def correct_uniq_mapper_edges(self, inserts_file_writer: TextIO, min_alignment_length: int):
        """
        ======================================================================
        Correct the alignment abundances which are biased due to partial
        alignments at the left and right parts of the genes.

        Assume that a read aligns against the left side of a gene but overlaps
        only 5 bases (usually more, just as an example)

        read = =====================
        gene =              ==============================================

        This alignment will not be counted by default for two reasons:
        1. The alignment length is too short which means that motus will not
            report it. In mOTUs the default is 75
        2. The aligner doesn't report it as it is too short for the aligner
            to confidently report it. In BWA this is 30 bases

        This means that abundances at the edges of genes are incorrectly
        reported, short genes suffer proportinally more of this error.


        Edge correction is using abundance of the parts of the gene
        that is correctly reported and extrapolates it to the full length of
        the gene. E.g here a gene has length 30, and anything below 5 is
        reported wrongly:
                      ABUNDANCE
        FALSE         CORRECT      FALSE
        =====|====================|=====

             =====================
             Use this abundance
             Extrapolate to full length
        abundance(full_gene) = len(full_gene) * abundance(trunc_gene) / len(trunc_gene)

        """

        mg_2_alignments = collections.defaultdict(list)
        for insert_name, bestAlignment in self._unique_mappers:
            mg, alignment_blocks = bestAlignment.get_mg_and_blocks()
            inserts_file_writer.write(f'{insert_name}\t{mg}\t1.0000\n')
            mg_2_alignments[mg].append((alignment_blocks, 1.0))



        mg_2_edge_corrected_insert_counts, mg_2_edge_corrected_base_counts = self._correct_edges(mg_2_alignments, 30)

        self._mg_2_edge_corrected_raw_uniquemapper_insert_counts = mg_2_edge_corrected_insert_counts
        self._mg_2_edge_corrected_raw_uniquemapper_base_counts = mg_2_edge_corrected_base_counts

    def combined_raw_counts(self):
        mg_2_edge_corrected_raw_insert_counts = collections.Counter()
        for mg, count in self._mg_2_edge_corrected_raw_uniquemapper_insert_counts.items():
            #print(mg, count, self._mg_2_edge_corrected_raw_multimapper_insert_counts.get(mg, 0.0))
            mg_2_edge_corrected_raw_insert_counts[mg] = count
        for mg, count in self._mg_2_edge_corrected_raw_multimapper_insert_counts.items():
            mg_2_edge_corrected_raw_insert_counts[mg] += count
        mg_2_edge_corrected_raw_base_counts = collections.Counter()
        for mg, count in self._mg_2_edge_corrected_raw_uniquemapper_base_counts.items():
            mg_2_edge_corrected_raw_base_counts[mg] = count
        for mg, count in self._mg_2_edge_corrected_raw_multimapper_base_counts.items():
            mg_2_edge_corrected_raw_base_counts[mg] += count

        self._mg_2_edge_corrected_raw_insert_counts = mg_2_edge_corrected_raw_insert_counts
        self._mg_2_edge_corrected_raw_base_counts = mg_2_edge_corrected_raw_base_counts

    def _norm_and_scale_counts2(self, mg_2_raw_counts):
        """ Raw insert counts are being normalised (and scaled).

        Normalisation per MG:
        norm(mg_count) = mg_cnt / sum(mgn_cnt/len(mgn) + ... + mgm_cnt/len(mgm))
        scaled(mg_count) = norm(mg_count) / sum(mg_counts)

        """
        tot_cnt = float(sum(mg_2_raw_counts.values()))
        denominator: float = sum([float(mgh_2_count[1]) / float(MOTUS_DB.get_length_by_mg(mgh_2_count[0])) for mgh_2_count in mg_2_raw_counts.items()])
        scaled_mg_2_counts = {}
        norm_mg_2_counts = {}
        for mg, count in mg_2_raw_counts.items():
            numerator = float(count) / float(MOTUS_DB.get_length_by_mg(mg))
            norm_count = numerator / denominator
            scaled_count = norm_count * tot_cnt
            scaled_mg_2_counts[mg] = scaled_count
            norm_mg_2_counts[mg] = norm_count

        return norm_mg_2_counts, scaled_mg_2_counts

    def norm_and_scale_counts(self):
        mg_2_edge_corrected_norm_insert_counts, mg_2_edge_corrected_scaled_insert_counts = self._norm_and_scale_counts2(self._mg_2_edge_corrected_raw_insert_counts)
        mg_2_edge_corrected_norm_base_counts, mg_2_edge_corrected_scaled_base_counts = self._norm_and_scale_counts2(self._mg_2_edge_corrected_raw_base_counts)
        self._mg_2_edge_corrected_norm_insert_counts = mg_2_edge_corrected_norm_insert_counts
        self._mg_2_edge_corrected_scaled_insert_counts = mg_2_edge_corrected_scaled_insert_counts
        self._mg_2_edge_corrected_norm_base_counts = mg_2_edge_corrected_norm_base_counts
        self._mg_2_edge_corrected_scaled_base_counts = mg_2_edge_corrected_scaled_base_counts

    def _get_alignment_blocks(self, alignments: List[pysam.AlignedSegment]) -> List[Tuple[int, int]]:
        '''
        Get the aligned positions from the alignments. Depending on the
        number of alignments (1 or 2) and the number of indels this number
        can be between 1 and n where n is the readlength. However, in most
        cases the result will be one Tuple for a singleton insert and
        two Tuples for a paired end insert

        Params:
            alignments: A list of pysam AlignedSegments

        Returns:
            A list of Tuples with aligned positions
        '''

        blocks = []
        for alignment in alignments:
            for block in alignment.get_blocks():
                blocks.append(block)
        blocks.sort(key=lambda a: a[0])
        return blocks

    def _filter_best_alignment(self, current_insert: Dict[str, List[pysam.AlignedSegment]],
                               random_mgc_resolver=True) -> BestAlignment:
        """Take a list of alignments against the mOTUs
        genes and find the mOTUs gene with the highest score.
        Will also report aligned bases

        Params:
            current_insert: All alignments associated to one insert.
                Can be paired end or singleton. Can have multiple
                alignments against different genes
            random_mgc_resolver: pick randomly if an insert aligns
                against multiple mgs from the same mgc
        Returns:
            A list of at least one mg this insert is assigned to
            together with the alignment positions
        """

        markergeneheader_2_orientation_2_alignments = collections.defaultdict(lambda: collections.defaultdict(list))
        for orientation, alignments in current_insert.items():
            for alignment in alignments:
                reference_name = alignment.reference_name
                markergeneheader_2_orientation_2_alignments[reference_name][orientation].append(alignment)
        markergeneheader_2_best_alignments = collections.defaultdict(list)
        for markergeneheader, orientation_2_alignments in markergeneheader_2_orientation_2_alignments.items():
            for orientation, alignments in orientation_2_alignments.items():
                best_alignment = alignments[0]
                if len(alignments) != 1:
                    max_score = max([alignment.get_tag('AS') for alignment in alignments])
                    best_alignment = [alignment for alignment in alignments if alignment.get_tag('AS') == max_score][0]
                markergeneheader_2_best_alignments[markergeneheader].append(best_alignment)
        all_alignment_scores = []
        for markergeneheader, best_alignments in markergeneheader_2_best_alignments.items():
            all_alignment_scores.append(sum([alignment.get_tag('AS') for alignment in best_alignments]))
        best_alignment_score = max(all_alignment_scores)
        markergeneheader_2_best_alignments2 = {}
        for markergeneheader, best_alignments in markergeneheader_2_best_alignments.items():
            if sum([alignment.get_tag('AS') for alignment in best_alignments]) == best_alignment_score:
                markergeneheader_2_best_alignments2[markergeneheader] = best_alignments
        best_mgs = BestAlignment()
        if len(markergeneheader_2_best_alignments2) > 1 and random_mgc_resolver:
            mgc_2_mg = collections.defaultdict(list)
            for mg, alns in markergeneheader_2_best_alignments2.items():
                mgc = MOTUS_DB.get_mgc_by_mg(mg)
                mgc_2_mg[mgc].append(mg)
            for mgc, mgs in mgc_2_mg.items():
                if len(mgs) == 1:
                    best_mgs.append(mgs[0], self._get_alignment_blocks(markergeneheader_2_best_alignments2[mgs[0]]))
                else:
                    picked_mg = random.choice(mgs)
                    best_mgs.append(picked_mg,
                                    self._get_alignment_blocks(markergeneheader_2_best_alignments2[picked_mg]))
        else:
            for mg in markergeneheader_2_best_alignments2.keys():
                best_mgs.append(mg, self._get_alignment_blocks(markergeneheader_2_best_alignments2[mg]))

        return best_mgs

    def count(self) -> None:
        """Umbrella count method
        Reads the alignments and stores them based on insert name.
        Then counts unique mappers and distributes multimappers
        based on the abundances of mgcs

        """

        bam_insert_iterator = self._bam_insert_iterator()
        for insert_name, alignments in bam_insert_iterator:
            self.appendmapper(insert_name, self._filter_best_alignment(alignments))

        logging.info('Finished reading alignment file ...')
        logging.info(f'Read {self.get_unique_mapper_count() + self.get_multi_mapper_count()} aligned inserts of which {round(100.0 * self.get_multi_mapper_count() / (self.get_unique_mapper_count() + self.get_multi_mapper_count()),2)}% are multimappers')
        inserts_file_writer = gzip.open(MOTUS_PARAMETERS.get_inserts_file(), 'wt')
        self.correct_uniq_mapper_edges(inserts_file_writer, MOTUS_PARAMETERS.get_minimal_alignment_length())
        self.correct_multi_mapper_edges(inserts_file_writer, MOTUS_PARAMETERS.get_minimal_alignment_length())
        inserts_file_writer.close()
        self.combined_raw_counts()
        self.norm_and_scale_counts()



    def _bam_insert_iterator(self) -> Generator[Tuple[str, Dict[str, List[pysam.AlignedSegment]]], None, None]:
        """Reads through a sorted BAM file and
        finds the best alignment(s) per insert
        """
        alignments = pysam.AlignmentFile(MOTUS_PARAMETERS.get_alignment_file(), 'r')
        motus_version = [entry for entry in alignments.header.to_dict()['PG'] if entry['ID'] == MOTUS_DB.get_full_sam_id()]
        header_valid = False
        if len(motus_version) > 0:
            if motus_version[0]['VN'] == MOTUS_DB.get_full_version():
                header_valid = True
        if not header_valid:
            if MOTUS_PARAMETERS.is_strict_db_mode():
                logging.error('mOTUs tool/database have changed and bam file is invalid. Please profile with updated database. Quitting ...')
                mutils.shutdown(1)
            else:
                logging.warning('mOTUs tool/database have changed and bam file is invalid. Lenient mode enabled, will continue but results might be broken ...')

        try:
            alignment: pysam.AlignedSegment = next(alignments)
        except StopIteration:
            alignments.close()
            logging.info(f'The alignmentfile {MOTUS_PARAMETERS.get_alignment_file()} has no valid alignments. Quitting ...')
            mutils.shutdown(1)
            return
        current_name, orientation = _get_orientation_of_aligned_segment_by_name(alignment)
        current_insert = collections.defaultdict(list)
        orientations = set()
        orientations.add(orientation)
        current_insert[orientation].append(alignment)
        minlength: int = MOTUS_PARAMETERS.get_minimal_alignment_length()
        readname = None

        for alignment in alignments:
            if MOTUS_DB.is_mg_blocked(alignment.reference_name):
                continue
            alnlength: int = sum(alignment.get_cigar_stats()[0][0:3])
            if alnlength < minlength:
                continue
            (readname, orientation) = _get_orientation_of_aligned_segment_by_name(alignment)
            if readname == current_name:
                current_insert[orientation].append(alignment)
                orientations.add(orientation)
            else:
                if len(orientations) == 3 or (len(orientations) > 1 and mentities.SIDENTIFIER in orientations):
                    raise Exception('An alignment cannot be Paired End and Single End at the same time. Problematic insert: {}'.format(readname))

                yield current_name, current_insert
                current_name = readname
                current_insert = collections.defaultdict(list)
                orientations = set()
                current_insert[orientation].append(alignment)
                orientations.add(orientation)

        if len(orientations) == 3 or (len(orientations) > 1 and mentities.SIDENTIFIER in orientations):
            raise Exception('An alignment cannot be Paired End and Single End at the same time. Problematic insert: {}'.format(readname))
        yield current_name, current_insert
        alignments.close()


class MGCCounter:
    """
    A class which takes care of
    reading, parsing and interpreting
    the inserts mapped against the mOTUs
    database.
    """


    def aggregate_mgc(self, mgh_2_scaled_counts, mgh_2_unscaled_counts):
        """Aggregate counts by markergenes by markergeneclusters

        """
        mgc_2_count = {}
        for mgh, count in mgh_2_scaled_counts.items():
            mgc = MOTUS_DB.get_mgc_by_mg(mgh)
            [scaled, unscaled] = mgc_2_count.get(mgc, [0.0, 0.0])
            scaled = scaled + count
            mgc_2_count[mgc] = [scaled, unscaled]

        for mgh, count in mgh_2_unscaled_counts.items():
            mgc = MOTUS_DB.get_mgc_by_mg(mgh)
            [scaled, unscaled] = mgc_2_count.get(mgc, [0.0, 0.0])
            unscaled = unscaled + count
            mgc_2_count[mgc] = [scaled, unscaled]
        return mgc_2_count


    def count(self) -> Dict[str, Mgc_values]:
        '''
        Entry Level method for this class
        Read the BAM file and counts abundances
        using different modes (insert_raw, insert_scaled,...)
        '''
        insertcounter = InsertCounter()
        logging.info('Reading alignment file ...')
        insertcounter.count()
        # now aggregate by MGC

        '''
        How to aggregate
        1. for each counting method (insert, base, norm, scaled)
            for each mg
                find mgc
                sum up value for mgc
        2. report
            for each counting method
                for each mgc
                one line with each counting method
        '''

        mgc_insert_raw = collections.defaultdict(lambda: 0.0)
        mgc_insert_norm = collections.defaultdict(lambda: 0.0)
        mgc_insert_scaled = collections.defaultdict(lambda: 0.0)
        mgc_base_raw = collections.defaultdict(lambda: 0.0)
        mgc_base_norm = collections.defaultdict(lambda: 0.0)
        all_mgcs = set()
        for (mg_data, mgc_data) in zip([insertcounter.get_mg_insert_raw(), insertcounter.get_mg_insert_norm(), insertcounter.get_mg_insert_scaled(), insertcounter.get_mg_base_raw(), insertcounter.get_mg_base_norm()], [mgc_insert_raw, mgc_insert_norm, mgc_insert_scaled, mgc_base_raw, mgc_base_norm]):
            for mg, abundance in mg_data.items():
                mgc = MOTUS_DB.get_mgc_by_mg(mg)
                mgc_data[mgc] = mgc_data[mgc] + abundance
                all_mgcs.add(mgc)

        mgc_2_all_counts = {}
        for mgc in all_mgcs:
            #Mgc_values = collections.namedtuple("Mgc_values", "insert_raw insert_norm insert_scaled base_raw base_norm")
            insert_raw = mgc_insert_raw[mgc]
            insert_norm = mgc_insert_norm[mgc]
            insert_scaled = mgc_insert_scaled[mgc]
            base_raw = mgc_base_raw[mgc]
            base_norm = mgc_base_norm[mgc]
            mgc_vals = Mgc_values(insert_raw=insert_raw, insert_norm=insert_norm, insert_scaled=insert_scaled, base_raw=base_raw, base_norm=base_norm)
            mgc_2_all_counts[mgc] = mgc_vals

        return mgc_2_all_counts



def calc_mgc() -> None:
    """
    Takes the BAM file created in the map_tax method and assigns individual alignments to marker genes and next to marker gene clusters.
    1. Read the name sorted alignments in the BAM file and pair by insert (1/2/S). Always have only one insert in memory.
    2. Find the best alignment per insert using additive paired alignment score.
    3. Check if this is a unique mapper or a multimapper. A multimapper that maps only against MG from the same MGC counts as unique mapper (and will be randomly assigned to a mg)
    4. Distribute the unique mappers to individual MGs
    5. Distribute the multimappers to MGs based on the fractional abundance of the unique mappers in those MGs
    6. Apply edge correction
    7. calculate using different count modes
    8. Group abundance by MGC and write to file

    Returns:
        None

    """
    logging.info('Starting mOTUs - calc_mgc routine - Calculating abundances per MGC ... ')
    mgc_counter = MGCCounter()
    mgc_2_counts = mgc_counter.count()

    with open(MOTUS_PARAMETERS.get_mgc_file(), 'w') as handle:
        header_line = MOTUS_DB.get_full_version()
        handle.write(f'#{header_line}\n')
        handle.write('MGC\tINSERT_RAW\tINSERT_NORM\tINSERT_SCALED\tBASE_RAW\tBASE_NORM\n')
        for mgc in sorted(mgc_2_counts.keys()):
            counts = mgc_2_counts[mgc]
            handle.write(f'{mgc}\t{round(counts.insert_raw, 4):.4f}\t{round(counts.insert_norm, 10):.10f}\t{round(counts.insert_scaled, 4):.4f}\t{round(counts.base_raw, 4):.4f}\t{round(counts.base_norm, 10):.10f}\n')

    logging.info('Finished mOTUs - calc_mgc routine - Calculating abundances per MGC ... ')
    return None


def calc_motu() -> None:
    """
    Takes the MGC file produced by calc_mgc and produces a mOTUs profile file.


    Returns
        None
    """

    mgc_file = MOTUS_PARAMETERS.get_mgc_file()
    has_header = False
    with open(mgc_file) as handle:
        first_line = handle.readline().strip()
        if first_line.startswith('#'):
            has_header = True
        motus_version = first_line.replace('#', '')
        if motus_version == MOTUS_DB.get_full_version():
            header_valid = True
        if not header_valid:
            if MOTUS_PARAMETERS.is_strict_db_mode():
                logging.error('mOTUs tool/database have changed and bam file is invalid. Please profile with updated database. Quitting ...')
                mutils.shutdown(1)
            else:
                logging.warning('mOTUs tool/database have changed and bam file is invalid. Lenient mode enabled, will continue but results might be broken ...')

    mgc_2_count = {}
    count_mode = MOTUS_PARAMETERS.get_count_mode()
    with open(mgc_file) as handle:
        if has_header:
            handle.readline()
        for entry in csv.DictReader(handle, delimiter='\t'):
            mgc_2_count[entry['MGC']] = float(entry[count_mode])
    motu_2_mgccounts = collections.defaultdict(lambda: collections.defaultdict(lambda: 0.0))

    for mgc, count in mgc_2_count.items():
        motu = MOTUS_DB.get_motu_by_mgc(mgc)
        mg = MOTUS_DB.get_mg_by_mgc(mgc)
        motu_2_mgccounts[motu][mg] += count

    motu_counts = {}
    for motu in sorted(list(motu_2_mgccounts.keys())):
        counts = list(motu_2_mgccounts[motu].values())
        median_count = statistics.median(counts)
        if len(counts) >= MOTUS_PARAMETERS.get_min_mgcs() or MOTUS_DB.is_unassigned_motu(motu):
            motu_counts[motu] = median_count

    mOTU_file = MotusFile()
    sample_2_motus_counts = {}
    sample_2_motus_counts[MOTUS_PARAMETERS.get_sample_name()] = motu_counts
    mOTU_file.set_mOTU_counts(sample_2_motus_counts, MOTUS_PARAMETERS.get_count_type())
    mOTU_file.set_min_mgcs(MOTUS_PARAMETERS.get_min_mgcs())
    mOTU_file.set_full_version(MOTUS_DB.get_full_version())
    mOTU_file.set_count_mode(MOTUS_PARAMETERS.get_count_mode())
    mOTU_file.write_mOTUs_file(MOTUS_PARAMETERS.get_motu_file(), False)
    mOTU_file.write_mOTUs_file(MOTUS_PARAMETERS.get_motu_file_relab(), True)
    return None


class CapitalisedHelpFormatter(argparse.HelpFormatter):
    def add_usage(self, usage, actions, groups, prefix=None):
        if prefix is None:
            prefix = ''
        return super(CapitalisedHelpFormatter, self).add_usage(usage, actions, groups, prefix)


def parse_map_tax():
    parser = argparse.ArgumentParser(usage = f'''Program: motus - a tool for marker gene-based OTU (mOTU) profiling
Version: {mutils.MOTUS_VERSION}
Reference: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand 
taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022). 
doi: https://doi.org/10.1186/s40168-022-01410-z
    
motus map_tax [options]

Input options:
   -f   FILE[ FILE]  input file(s) for reads in forward orientation, fastq(.gz)-formatted
   -r   FILE[ FILE]  input file(s) for reads in reverse orientation, fastq(.gz)-formatted
   -s   FILE[ FILE]  input file(s) for unpaired reads, fastq(.gz)-formatted


Output options:
   -o   FILE         output file name

Algorithm options:
   -l   INT          min length of the alignment (bp) [75]
   -t   INT          number of threads [1]
   -v   INT          verbosity level: 1=error, 2=warning, 3=message, 4+=debugging [3]
      ''', formatter_class=CapitalisedHelpFormatter,add_help=False)

    # Input options
    parser.add_argument("-f", nargs="+",default=[])  # input files(s) for reads in forward orientation, fastq(.gz)-formatted
    parser.add_argument("-r", nargs="+",default=[])  # input files(s) for reads in reverse orientation, fastq(.gz)-formatted
    parser.add_argument("-s", nargs="+", default=[])  # input files(s) for unpaired reads, fastq(.gz)-formatted
    #parser.add_argument("-db")  # provide a different database directory

    # Output options
    parser.add_argument("-o", required=True)  # output file name
    #parser.add_argument("-b", action="store_true")  # save the result in BAM format

    # ALgorithm options
    parser.add_argument("-l", type=int, default=75)  # min length of the alignment (bp) [75]
    parser.add_argument("-t", type=int, default=1)  # number of threads
    parser.add_argument("-v", type=int, default=1)  # verbodisty level:

    args = parser.parse_args(sys.argv[2:])

    # print usage and exit if no arguments are passed
    if sys.argv[2:] == []:
        parser.print_usage()
        mutils.shutdown(1)

    # convert string arguments into Pathlib objects
    forward_files = [pathlib.Path(el) for el in args.f]
    reverse_files = [pathlib.Path(el) for el in args.r]
    unpaired_files = [pathlib.Path(el) for el in args.s]
    alignment_file = pathlib.Path(args.o)
    mutils.startup()
    threads = args.t
    min_alignment_length = args.l

    MOTUS_DB.load_motus_db(mutils.DEFAULT_MOTUS_MGDB_LOCATION)



    MOTUS_PARAMETERS.set_read_files(forward_files, reverse_files, unpaired_files, check_files=True)
    MOTUS_PARAMETERS.set_alignment_file(alignment_file, required_to_exist=False)
    MOTUS_PARAMETERS.set_minimal_alignment_length(min_alignment_length)
    MOTUS_PARAMETERS.set_threads(threads)
    map_tax()



def parse_batch_profile():
    '''
    Hidden routine which takes a list of bam files to parse, create motus and mgc files.

    Advantage --> needs to load the database only once
    '''

    parser = argparse.ArgumentParser(usage=f'''Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: {mutils.MOTUS_VERSION}
    Reference: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand 
    taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022). 
    doi: https://doi.org/10.1186/s40168-022-01410-z

    motus batch_profile [options]

    Input options:
       -f  FILE[ FILE]  input tsv file. First column=samplename, second column=bam file

    Algorithm options:
       -g  INT          number of marker genes cutoff: 1=higher recall, 6=higher precision [3]
       -l  INT          min length of the alignment (bp) [75]
       -v  INT          verbosity level: 1=error, 2=warning, 3=message, 4+=debugging [1]
       -y  STR          type of read counts [INSERT_SCALED]
                        Values: [INSERT_RAW, INSERT_NORM, INSERT_SCALED, BASE_RAW, BASE_NORM]
    ]''', formatter_class=CapitalisedHelpFormatter, add_help=False)

    # Input options
    parser.add_argument("-f", required=True)  # input file(s) for reads in forward direction

    parser.add_argument("-g", type=int, default=3,
                        choices=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10])  # number of marker genes cutoff
    parser.add_argument("-l", type=int, default=75)  # min length of the alignment (bp) [75]
    parser.add_argument("-t", type=int, default=1)  # number of thread [1]
    parser.add_argument("-v", type=int, default=1)  # verbosity level
    parser.add_argument("-y", type=str, default='INSERT_SCALED',
                        choices=['INSERT_RAW', 'INSERT_NORM', 'INSERT_SCALED', 'BASE_RAW', 'BASE_NORM'])

    args = parser.parse_args(sys.argv[2:])
    if sys.argv[2:] == []:
        parser.print_usage()
        mutils.shutdown(1)

    map_file = args.f
    samplename_2_files = {}
    with open(map_file) as handle:
        for line in handle:
            [samplename, bamfile] = line.strip().split('\t')
            bamfile = pathlib.Path(bamfile)
            if not bamfile.exists():
                logging.error(f'Submitted BAM file {bamfile} does not exist. Quitting ...')
                mutils.shutdown(1)
            if not str(bamfile).endswith('.bam'):
                logging.error(f'Submitted BAM file {bamfile} does not end with .bam. Probably malformed file. Quitting ...')
                mutils.shutdown(1)

            mgc_file = str(bamfile).rsplit('.bam', 1)[0] + '.mgc'
            inserts_file = str(bamfile).rsplit('.bam', 1)[0] + '.inserts.gz'
            motu_file = str(bamfile).rsplit('.bam', 1)[0]
            if samplename in samplename_2_files:
                logging.error(f'Submitted samplename {samplename} duplicated. Quitting ...')
                mutils.shutdown(1)
            for (obf,mf,mof,iof) in samplename_2_files.values(): #output_bam_file, mgc_file, motu_file, inserts_output_file
                if obf.samefile(bamfile):
                    logging.error(f'Submitted BAM file {bamfile} duplicated. Quitting ...')
                    mutils.shutdown(1)
            samplename_2_files[samplename] = (bamfile, pathlib.Path(mgc_file), pathlib.Path(motu_file), pathlib.Path(inserts_file))


    mutils.startup()
    threads = args.t
    min_alignment_length = args.l


    MOTUS_DB.load_motus_db(mutils.DEFAULT_MOTUS_MGDB_LOCATION)

    for samplename, (bamfile, mgc_file, motu_file, inserts_file) in samplename_2_files.items():


        MOTUS_PARAMETERS.set_alignment_file(bamfile)
        MOTUS_PARAMETERS.set_mgc_file(mgc_file, required_to_exist=False)
        MOTUS_PARAMETERS.set_inserts_file(inserts_file, required_to_exist=False)

        MOTUS_PARAMETERS.set_motu_file(motu_file, required_to_exist=False)
        MOTUS_PARAMETERS.set_sample_name(samplename)
        MOTUS_PARAMETERS.set_minimal_alignment_length(min_alignment_length)
        MOTUS_PARAMETERS.set_threads(threads)
        MOTUS_PARAMETERS.set_count_mode(args.y)
        MOTUS_PARAMETERS.set_minimal_number_of_mgcs(args.g)
        MOTUS_PARAMETERS.enable_lenient_mode()
        calc_mgc()
        calc_motu()
    mutils.shutdown(0)

def parse_profile():
    parser = argparse.ArgumentParser(usage = f'''Program: motus - a tool for marker gene-based OTU (mOTU) profiling
Version: {mutils.MOTUS_VERSION}
Reference: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand 
taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022). 
doi: https://doi.org/10.1186/s40168-022-01410-z
    
motus profile [options]

Input options:
   -f  FILE[ FILE]  input file(s) for reads in forward orientation, fastq(.gz)-formatted
   -r  FILE[ FILE]  input file(s) for reads in reverse orientation, fastq(.gz)-formatted
   -s  FILE[ FILE]  input file(s) for unpaired reads, fastq(.gz)-formatted
   -n  STR          sample name ['unnamed sample']

Output options:
   -o  FILE         output file name [required]

Algorithm options:
   -g  INT          number of marker genes cutoff: 1=higher recall, 6=higher precision, 10=maximum [3]
   -l  INT          min length of the alignment (bp) [75]
   -t  INT          number of threads [1]
   -v  INT          verbosity level: 1=error, 2=warning, 3=message, 4+=debugging [1]
   -y  STR          type of read counts [INSERT_SCALED]
                    Values: [INSERT_RAW, INSERT_NORM, INSERT_SCALED, BASE_RAW, BASE_NORM]
]''', formatter_class=CapitalisedHelpFormatter,add_help=False)

    # Input options
    parser.add_argument("-f", nargs="+", default=[])  # input file(s) for reads in forward direction
    parser.add_argument("-r", nargs="+", default=[])  # input file(s) for reads in reverse direction
    parser.add_argument("-s", nargs="+", default=[])  # input file(s) for unpaired reads
    parser.add_argument("-n", type=str, default='unnamed sample')  # sample name

    # Output options
    parser.add_argument("-o", required=True)
    parser.add_argument("-g", type=int, default=3, choices=[1,2,3,4,5,6,7,8,9,10])  # number of marker genes cutoff
    parser.add_argument("-l", type=int, default=75)  # min length of the alignment (bp) [75]
    parser.add_argument("-t", type=int, default=1)  # number of thread [1]
    parser.add_argument("-v", type=int, default=1)  # verbosity level
    parser.add_argument("-y", type=str, default='INSERT_SCALED', choices=['INSERT_RAW', 'INSERT_NORM', 'INSERT_SCALED', 'BASE_RAW', 'BASE_NORM'])

    args = parser.parse_args(sys.argv[2:])
    if sys.argv[2:] == []:
        parser.print_usage()
        mutils.shutdown(1)

    forward_files = [pathlib.Path(el) for el in args.f]
    reverse_files = [pathlib.Path(el) for el in args.r]
    unpaired_files = [pathlib.Path(el) for el in args.s]

    motu_file = pathlib.Path(args.o)
    alignment_file = pathlib.Path(args.o + '.bam')
    mgc_file = pathlib.Path(args.o + '.mgc')
    inserts_file = pathlib.Path(args.o + '.inserts.gz')
    mutils.startup()
    threads = args.t
    min_alignment_length = args.l
    samplename = args.n


    MOTUS_DB.load_motus_db(mutils.DEFAULT_MOTUS_MGDB_LOCATION)


    MOTUS_PARAMETERS.set_read_files(forward_files, reverse_files, unpaired_files, check_files=True)
    MOTUS_PARAMETERS.set_alignment_file(alignment_file, required_to_exist=False)
    MOTUS_PARAMETERS.set_mgc_file(mgc_file,required_to_exist=False)
    MOTUS_PARAMETERS.set_inserts_file(inserts_file, required_to_exist=False)
    MOTUS_PARAMETERS.set_motu_file(motu_file, required_to_exist=False)
    MOTUS_PARAMETERS.set_sample_name(samplename)
    MOTUS_PARAMETERS.set_minimal_alignment_length(min_alignment_length)
    MOTUS_PARAMETERS.set_threads(threads)

    MOTUS_PARAMETERS.set_count_mode(args.y)
    MOTUS_PARAMETERS.set_minimal_number_of_mgcs(args.g)

    map_tax()

    calc_mgc()
    calc_motu()
    mutils.shutdown(0)


def parse_calc_mgc():
    parser = argparse.ArgumentParser(usage = f'''Program: motus - a tool for marker gene-based OTU (mOTU) profiling
Version: {mutils.MOTUS_VERSION}
Reference: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand 
taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022). 
doi: https://doi.org/10.1186/s40168-022-01410-z
    
motus calc_mgc [options]

Input options:
   -i  FILE         provide the SAM or BAM input file (output of motus map_tax)

Output options:
   -o  FILE         output file name

Algorithm options:
   -l  INT          min length of the alignment (bp) [75]
   -v  INT          verbosity level: 1=error, 2=warning, 3=message, 4+=debugging [3]''', formatter_class=CapitalisedHelpFormatter,add_help=False)

    parser.add_argument("-i", type=str, required=True)  # provide a SAM or BAM input file (or list of files) output of motus map_tax
    parser.add_argument("-o", required=True)  # output file name [stdout]
    parser.add_argument("-l", type=int, default=75)  # min length of the alignment (bp) [75]
    parser.add_argument("-v", type=int, default=1)  # verbosity level

    args = parser.parse_args(sys.argv[2:])
    # print usage and exit if no arguments are passed
    if sys.argv[2:] == []:
        parser.print_usage()
        mutils.shutdown(1)

    alignment_file = pathlib.Path(args.i)
    mgc_file = pathlib.Path(args.o)
    inserts_file = pathlib.Path(args.o + '.inserts')


    mutils.startup()
    min_alignment_length = args.l

    MOTUS_DB.load_motus_db(mutils.DEFAULT_MOTUS_MGDB_LOCATION)

    MOTUS_PARAMETERS.set_alignment_file(alignment_file, required_to_exist=True)
    MOTUS_PARAMETERS.set_mgc_file(mgc_file,required_to_exist=False)
    MOTUS_PARAMETERS.set_inserts_file(inserts_file, required_to_exist=False)
    MOTUS_PARAMETERS.set_minimal_alignment_length(min_alignment_length)
    MOTUS_PARAMETERS.set_threads(1)
    calc_mgc()
    mutils.shutdown(0)

def assign_taxonomy(input_motu_file: pathlib.Path, output_motu_file: pathlib.Path, aggregate: bool, use_representative_taxonomy: bool, taxonomy_to_use: str, taxonomic_rank: TaxonomicRank) -> None:
    """
    Assign taxonomy to existing motusfile and
    potentially aggregate
    Params:
        input_motu_file: a file with motus abundances. can also be merged. Not allowed to be aggregate
        output_motu_file: the output file where results will be stored
        aggregate: should motus abundances be aggregated at the provided level
        use_representative_taxonomy: mOTUs supports 2 taxonomy types. representative and concensus.
            Setting this flag will disable consensus and enable representative taxonomy
        taxonomy_to_use: The taxonomy to use. Currently only GTDB supported
        taxonomic_rank: One of the standard ranks to report
    """
    mf = MotusFile()
    mf.read_mOTUs_file(input_motu_file)
    x = 0

def merge_profiles(motus_file_paths: List[pathlib.Path], output_motus_file_path: pathlib.Path) -> None:
    """Routine which merges mOTU profiles.

    - Will read the mOTUs profiles into Motu_file objects
    - Merge them into a new Motu_file object (only if parameters are the same, throw an exception instead)
    - Write the merged profile into the output file

    """

    #TODO check for empty motus files

    logging.info('Starting mOTUs - merge routine - Merging of mOTUs profile files ... ')
    motus_files = []
    logging.info(f'There are {len(motus_file_paths)} input profile files.')
    for cnt, motus_file_path in enumerate(motus_file_paths, 1):
        if cnt % 100 == 0:
            logging.info(f'{cnt} / {len(motus_file_paths)} processed')
        mf = MotusFile()
        mf.read_mOTUs_file(motus_file_path)
        if not mf.has_motus_with_abundance():
            logging.info(f'mOTUs4 file has no mOTUs with abundance. Will be ignored. {motus_file_path}')
            continue
        motus_files.append(mf)
    logging.info(f'{cnt} / {len(motus_file_paths)} processed')
    mf = MotusFile()
    mf.merge_profiles(motus_files)
    if mf.has_counts():
        mf.write_mOTUs_file(output_motus_file_path, relabundance=False)
        mf.write_mOTUs_file(pathlib.Path(str(output_motus_file_path) + '.relab'), relabundance=True)
    else:
        #TODO shouldnt the relabundance be True here?
        mf.write_mOTUs_file(output_motus_file_path, relabundance=False)
    logging.info('Finished mOTUs - merge routine - Merging of mOTUs profile files ... ')


def parse_merge():
    parser = argparse.ArgumentParser(usage=f'''Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: {mutils.MOTUS_VERSION}
    Reference: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand 
    taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022). 
    doi: https://doi.org/10.1186/s40168-022-01410-z

    motus merge [options]

        Input options:
           -i  FILE[ FILE]  A list of mOTUs profile files or a text file with one line
                            per mOTUs profile files to be merged
                            

        Output options:
           -o  FILE  output file name       

          ''', formatter_class=CapitalisedHelpFormatter, add_help=False)


    parser.add_argument("-i", nargs="+", required=True)

    parser.add_argument("-o", required=True)
    args = parser.parse_args(sys.argv[2:])

    if sys.argv[2:] == []:
        parser.print_usage()
        mutils.shutdown(1)

    mutils.startup()

    input_mOTUs_files = [pathlib.Path(el) for el in args.i]
    output_mOTUs_file = pathlib.Path(args.o)

    if len(input_mOTUs_files) == 1: # this means that you provided a file with one line per motus profile
        input_mOTUs_files_tmp = []
        with open(input_mOTUs_files[0]) as handle:
            for line in handle:
                input_mOTUs_files_tmp.append(pathlib.Path(line.strip()))
        input_mOTUs_files = input_mOTUs_files_tmp


    input_mOTUs_files = sorted(input_mOTUs_files)
    merge_profiles(input_mOTUs_files, output_mOTUs_file)

    mutils.shutdown(0)


def download_genomes(keyword: str, motusSearchDB: MotusSearchDB, output_folder: pathlib.Path, output_file: pathlib.Path, download_representative_genomes_only: bool) -> None:
    logging.info(f'Searching for keyword: {keyword}.')
    genomes_to_download = motusSearchDB.search_for_genomes(keyword, only_representatives=download_representative_genomes_only)
    logging.info(f'Found: {len(genomes_to_download)} hits.')

    logging.info(f'Found {len(genomes_to_download)} genomes. Writing genome information to {output_file}')
    with open(output_file, 'w') as handle:
        handle.write('GENOME\tMOTU\tPATH\tDOMAIN\tPHYLUM\tCLASS\tORDER\tFAMILY\tGENUS\tSPECIES\n')
        for genome in genomes_to_download:
            genome_path = motusSearchDB.get_genome_path(genome)
            genome_motu = motusSearchDB.get_genome_motu(genome)
            genome_tax = motusSearchDB.get_genome_tax(genome)
            handle.write(f'{genome}\t{genome_motu}\t{genome_path}\t{genome_tax}\n')
    logging.info(f'Finished writing genome information to {output_file}')
    if output_folder:
        if output_folder.is_file():
            logging.error('Output Path exists and is file. Cannot download genomes to this location')
            mutils.shutdown(1)
        logging.info(f'Downloading genomes to {output_folder}')
        output_folder.mkdir(exist_ok=True, parents=True)
        for cnt, genome in enumerate(genomes_to_download, 1):
            genome_path = str(motusSearchDB.get_genome_path(genome))
            destpath = str(output_folder) + '/' + str(genome_path).split('/')[-1]
            logging.info(f'Downloading genome ({cnt} / {len(genomes_to_download)}) {genome} to {destpath}')
            urllib.request.urlretrieve(genome_path, destpath)
        logging.info(f'Finished downloading genomes')


def parse_classify():
    parser = argparse.ArgumentParser(usage=f'''Program: motus - a tool for marker gene-based OTU (mOTU) profiling
     Version: {mutils.MOTUS_VERSION}
     Reference: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand 
     taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022). 
     doi: https://doi.org/10.1186/s40168-022-01410-z

     motus classify [options]

         Options:

            -i        Text file with fasta formatted (gzip allowed) genome
                      files which will be associated with existing mOTUs
            -o        Output file. One line per genome with associated mOTU.
            -t        Number of threads (default = 1)

           ''', formatter_class=CapitalisedHelpFormatter, add_help=False)


    parser.add_argument("-i", required=True)
    parser.add_argument("-o", required=True)
    parser.add_argument("-t", default=1, type=int)
    args = parser.parse_args(sys.argv[2:])

    if sys.argv[2:] == []:
        parser.print_usage()
        mutils.shutdown(1)
    mutils.startup()

    genome_files = []
    with open(args.i) as handle:
        for line in handle:
            genome_file = pathlib.Path(line.strip())
            if not genome_file.exists():
                logging.error(f'Genome file: {line.strip()} does not exist. Quitting ...')
                mutils.shutdown(1)
            genome_files.append(genome_file)
    output_file = args.o

    classify(genome_files, output_file, args.t)
    mutils.shutdown(0)



def classify(genome_files: List[pathlib.Path], output_file: pathlib.Path, threads: int = 1):
    """
    Takes a list of genome files and associates them with existing mOTUs.
    A genome can either be:
    - classified with a mOTU (=mOTUXXX)
    - not have enough markergenes to be classified (=notEnoughMGs)
    - have enough markergenes to be classified but can not be represented by an existing mOTU (=Novel)

    Params:
        genome_files: A list with pathlib.Path objects all pointing to an existing genome file
        output_file: The output file where the classification should be recorded
    Returns:
        None
    """

    logging.info(f'Starting mOTUs classify:')
    logging.info(f'\tInput = {len(genome_files)} genomes.')
    logging.info(f'\tOutput will be written to {output_file}')

    logging.info(f'\t Running fetchMGs on genomes.')


    from fetchmgs import fetchmgs
    fetchmgs_tmp_folder = pathlib.Path(str(output_file) + '_classify_tmp')
    if False:
        fetchmgs.extraction_genomes(genome_files, fetchmgs_tmp_folder, 'genome', threads, True)
    logging.info(f'\t Finished running fetchMGs on genomes.')
    logging.info(f'\t Collecting fetchMGs results')


    MOTUS_DB.load_motus_db(mutils.DEFAULT_MOTUS_MGDB_LOCATION, False)


    genome_2_fetchmg_file = {}

    for genome_fetchmgs_file in fetchmgs_tmp_folder.glob('**/*fetchMGs.fna'):
        genome_name = str(genome_fetchmgs_file.name).replace('.fetchMGs.fna', '')
        genome_2_fetchmg_file[genome_name] = genome_fetchmgs_file

    genome_2_markergenes = collections.defaultdict(list)
    motus_marker_genes = set(MOTUS_DB.get_core_motus_mgs())
    for genome_file in genome_files:
        genome_name = str(genome_file.name)
        genome_fetchmgs_file = genome_2_fetchmg_file[genome_name]


        with open(genome_fetchmgs_file) as handle:
            for (header, sequence) in FastaIO.SimpleFastaParser(handle):
                cog = header.rsplit('.', 1)[-1]
                if cog in motus_marker_genes:
                    genome_2_markergenes[genome_name].append((cog, sequence))
    genomes_removed_notenoughmgs = 0
    with open(fetchmgs_tmp_folder.joinpath('motus_classify.fna'), 'w') as seqhandle, open(fetchmgs_tmp_folder.joinpath('motus_classify.tsv'), 'w') as tsvhandle:
        tsvhandle.write('GENOME\tNUM_MGS\tMGS\n')
        for genome_name, cog_2_sequence in genome_2_markergenes.items():
            if len(cog_2_sequence) > 5:
                for cog, sequence in cog_2_sequence:
                    seqhandle.write(f'>{genome_name}.{cog}\n{sequence}\n')
            else:
                genomes_removed_notenoughmgs += 1
            cogs = sorted([x[0] for x in cog_2_sequence])
            cogs_str = ','.join(cogs)
            tsvhandle.write(f'{genome_name}\t{len(cogs)}\t{cogs_str}\n')
    logging.info(f'\t Finished collecting fetchMGs results. Genomes = {len(genome_2_markergenes)}, Genomes with enough MGs = {len(genome_2_markergenes) - genomes_removed_notenoughmgs}')















def parse_downloadDB():
    parser = argparse.ArgumentParser(usage=f'''Program: motus - a tool for marker gene-based OTU (mOTU) profiling
     Version: {mutils.MOTUS_VERSION}
     Reference: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand 
     taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022). 
     doi: https://doi.org/10.1186/s40168-022-01410-z

     motus downloadDB [options]
     
         Options:

            -f        Force download even when database is already present

           ''', formatter_class=CapitalisedHelpFormatter, add_help=False)

    parser.add_argument("-f", action="store_true")
    args = parser.parse_args(sys.argv[2:])

    force_download = False
    if args.f:
        force_download = True

    mutils.startup()

    if mutils.DEFAULT_MOTUS_MGDB_LOCATION_MARKER.exists() and not force_download:
        logging.info('Database already downloaded and -f not set. All good.')
        mutils.shutdown(0)
    if mutils.DEFAULT_MOTUS_MGDB_LOCATION_MARKER.exists() and force_download:
        logging.info('Database already downloaded and -f set. Will delete current database and download again.')
        shutil.rmtree(mutils.DEFAULT_MOTUS_MGDB_LOCATION)

    logging.info('Start downloading mOTUs marker gene database. ~6GB')
    dest_tar_gz_file = mutils.DEFAULT_MOTUS_MGDB_PARENT_LOCATION.joinpath('db_mOTU.tar.gz')
    if dest_tar_gz_file.is_file():
        dest_tar_gz_file.unlink()
    urllib.request.urlretrieve(mutils.MOTUS_MGDB_REMOTE_LOCATION, str(dest_tar_gz_file))
    logging.info('Finished downloading mOTUs marker gene database.')

    logging.info('Start un-taring mOTUs marker gene database.')
    if mutils.DEFAULT_MOTUS_MGDB_LOCATION.exists():
        shutil.rmtree(mutils.DEFAULT_MOTUS_MGDB_LOCATION)

    with tarfile.open(dest_tar_gz_file, 'r') as t:
        t.extractall(mutils.DEFAULT_MOTUS_MGDB_PARENT_LOCATION)
    mutils.DEFAULT_MOTUS_MGDB_LOCATION_MARKER.touch(exist_ok=True)
    logging.info('Finished untaring mOTUs marker gene database.')
    mutils.shutdown(0)


def parse_download():
    parser = argparse.ArgumentParser(usage=f'''Program: motus - a tool for marker gene-based OTU (mOTU) profiling
     Version: {mutils.MOTUS_VERSION}
     Reference: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand 
     taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022). 
     doi: https://doi.org/10.1186/s40168-022-01410-z

     motus download [options]

         Output options:
            -s  FILE  Genome metadata file
            -o  PATH  Genome output folder. Only required when 
                      -l is not set     

         Options:

            -l        Skip genome download. Only create genome report file
            -r        Download only representative genomes
            -w   STR  Keyword: Can be mOTU, genome name or taxonomy 

           ''', formatter_class=CapitalisedHelpFormatter, add_help=False)

    # Output options
    parser.add_argument("-o", required=False)
    parser.add_argument("-s", required=True)

    parser.add_argument("-l", action="store_true")
    parser.add_argument("-r", action="store_true")
    parser.add_argument("-w", type=str, required=True)

    args = parser.parse_args(sys.argv[2:])
    if sys.argv[2:] == []:
        parser.print_usage()
        mutils.shutdown(1)

    output_folder = args.o
    output_file = pathlib.Path(args.s)
    download_genome_metadata_only = False
    download_representative_genomes_only = False
    if args.l:
        download_genome_metadata_only = True
        output_folder = None
    if args.r:
        download_representative_genomes_only = True
    keyword = args.w

    mutils.startup()
    if not download_genome_metadata_only:
        if not output_folder:
            logging.error('Output folder must be set unless -l flag is used. Quitting ...')
            mutils.shutdown(1)
        output_folder = pathlib.Path(output_folder)

    MOTUS_DB.load_motus_db(mutils.DEFAULT_MOTUS_MGDB_LOCATION, load=False)
    motus_search_db = MotusSearchDB(MOTUS_DB.motus_mv_taxonomy_file, MOTUS_DB.genome_metadata_file)
    download_genomes(keyword, motus_search_db, output_folder, output_file, download_representative_genomes_only)
    mutils.shutdown(0)


def parse_taxonomy():
    parser = argparse.ArgumentParser(usage=f'''Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: {mutils.MOTUS_VERSION}
    Reference: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand 
    taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022). 
    doi: https://doi.org/10.1186/s40168-022-01410-z

    motus taxonomy [options]

        Input options:
           -i  FILE  a mOTUs profile, produced by profile, calc_motu or merge

        Output options:
           -o  FILE  output file name       

        Options:

           -t   STR  Taxonomy to use [GTDB]
           -l   STR  Taxonomy level [domain, phylum, class, order, family, genus, species]
           -a        Aggregate values at taxonomic level
           -r        Use taxonomy of representative. If not set, use the consensus taxonomy

          ''', formatter_class=CapitalisedHelpFormatter, add_help=False)

    parser.add_argument("-i", required=True)
    parser.add_argument("-o", required=True)
    parser.add_argument("-a", action="store_true")
    parser.add_argument("-r", action="store_true")
    parser.add_argument("-t", type=str, default='GTDB', choices=['GTDB'])
    parser.add_argument("-l", type=int, default='species', choices=['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'])

    args = parser.parse_args(sys.argv[2:])

    if sys.argv[2:] == []:
        parser.print_usage()
        mutils.shutdown(1)
    input_motu_file = pathlib.Path(args.i)
    output_motu_file = pathlib.Path(args.o)
    aggregate = False
    use_representative_taxonomy = False
    taxonomy_to_use = 'GTDB'
    taxonomic_rank = args.l
    if args.a:
       aggregate = True
    if args.r:
        use_representative_taxonomy = True
    if args.t != 'GTDB':
        logging.error(f'Unknown taxonomy {args.t}. Quitting ...')
        mutils.shutdown(1)

    mutils.startup()
    if not input_motu_file.exists():
        logging.error(f'Input file {input_motu_file} does not exist. Quitting ...')
        mutils.shutdown(1)

    MOTUS_DB.load_motus_db(mutils.DEFAULT_MOTUS_MGDB_LOCATION)
    mf = MotusFile()
    mf.read_mOTUs_file(input_motu_file)
    assign_taxonomy(input_motu_file, output_motu_file, aggregate, use_representative_taxonomy, taxonomy_to_use, taxonomic_rank)
    mutils.shutdown(0)


def parse_calc_motu():
    parser = argparse.ArgumentParser(usage = f'''Program: motus - a tool for marker gene-based OTU (mOTU) profiling
Version: {mutils.MOTUS_VERSION}
Reference: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand 
taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022). 
doi: https://doi.org/10.1186/s40168-022-01410-z
    
motus calc_motu [options]

    Input options:
       -n  STR   sample name [unnamed sample]
       -i  FILE  provide the mgc abundance table (output of motus calc_mgc)
    
    Output options:
       -o  FILE  output file name 
       -c        print result as counts instead of relative abundances
    
    Algorithm options:
       -g   INT   number of marker genes cutoff: 1=higher recall, 6=higher precision, 10=maximum [3]
       -v   INT   verbosity level: 1=error, 2=warning, 3=message, 4+=debugging [3]
       -y  STR    type of read counts [INSERT_SCALED]
                    Values: [INSERT_RAW, INSERT_NORM, INSERT_SCALED, BASE_RAW, BASE_NORM]
      
      ''', formatter_class=CapitalisedHelpFormatter,add_help=False)


    parser.add_argument("-n", type=str, default='unnamed sample')  # sample name
    parser.add_argument("-i", required=True)  # provide the mgc abundance table(output of motus calc_mgc)
    parser.add_argument("-o", required=True)  # output fil name [stdout]
    parser.add_argument("-y", type=str, default='INSERT_SCALED', choices=['INSERT_RAW', 'INSERT_NORM', 'INSERT_SCALED', 'BASE_RAW', 'BASE_NORM'])
    parser.add_argument("-g", type=int, default=3, choices=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10])  # number of marker genes cutoff

    args = parser.parse_args(sys.argv[2:])
    if sys.argv[2:] == []:
        parser.print_usage()
        mutils.shutdown(1)


    mgc_file = pathlib.Path(args.i)
    motu_file = pathlib.Path(args.o)
    mutils.startup()
    samplename = args.n


    MOTUS_DB.load_motus_db(mutils.DEFAULT_MOTUS_MGDB_LOCATION)

    MOTUS_PARAMETERS.set_mgc_file(mgc_file, required_to_exist=True)
    MOTUS_PARAMETERS.set_motu_file(motu_file, required_to_exist=False)
    MOTUS_PARAMETERS.set_sample_name(samplename)
    MOTUS_PARAMETERS.set_threads(1)
    MOTUS_PARAMETERS.set_count_mode(args.y)
    MOTUS_PARAMETERS.set_minimal_number_of_mgcs(args.g)
    calc_motu()
    mutils.shutdown(0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage = f'''Program: motus - a tool for marker gene-based OTU (mOTU) profiling
Version: {mutils.MOTUS_VERSION}
Reference: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand 
taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022). 
doi: https://doi.org/10.1186/s40168-022-01410-z
    
motus <command> [options]
    
    -- Taxonomic profiling
          profile     Perform taxonomic profiling (map_tax + calc_mgc + calc_motu) in a single step

          map_tax     Map reads to the marker gene database
          calc_mgc    Calculate marker gene cluster (MGC) abundance
          calc_motu   Summarize MGC abundances into a mOTU profile
    
    -- Utilities
          download    Download genomes associated with mOTUs
          downloadDB  Download the mOTUs marker gene database
          merge       Merge several taxonomic profiling results into one table
          classify    Add taxonomic information to mOTUs


    Type motus <command> to print the help menu for a specific command
    ''',formatter_class=CapitalisedHelpFormatter,add_help=False)

    parser.add_argument('command', choices=["profile", "map_tax", "calc_mgc", "calc_motu", "download", "merge", "downloadDB", "batch_profile", "classify"])
    args: argparse.Namespace = parser.parse_args(sys.argv[1:2])
    if args.command == 'profile':
        parse_profile()
    if args.command == 'batch_profile':
        parse_batch_profile()
    elif args.command == 'merge':
        parse_merge()
    elif args.command == 'map_tax':
        parse_map_tax()
    elif args.command == 'calc_mgc':
        parse_calc_mgc()
    elif args.command == 'calc_motu':
        parse_calc_motu()
    elif args.command == 'download':
        parse_download()
    elif args.command == 'downloadDB':
        parse_downloadDB()
    # elif args.command == 'taxonomy':
    #     parse_taxonomy()
    #
    elif args.command == 'classify':
        parse_classify()
    # elif args.command == 'prep_long':
    #     logging.error('Command prep_long not implemented yet')
    else:
        parser.print_usage()
        print(f'Unrecognized command {args}')
        mutils.shutdown(1)
    mutils.shutdown(0)
