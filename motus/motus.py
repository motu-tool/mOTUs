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
# Copyright (c) ${2025} ${SunagawaLab}.
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
from motus import mutils, mfind, mentities
from motus.mentities import MOTUS_PARAMETERS
from motus.mentities import MOTUS_DB
import tqdm


__author__ = ('Hans-Joachim Ruscheweyh (hansr@ethz.ch), '
              'Lilith Feer, '
              'Marija Dmitrijeva, '
              'Kang Li, '
              'Florian Ruscheweyh'
              'Daniel Mende, '
              'Georg Zeller, '
              'Shinichi Sunagawa')
__version__ = mutils.MOTUS_VERSION


__date__ = '05 September 2025'
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
            pg_header['CL'] = f'motus.py map_tax -l {minlength}'
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

    logging.info('Sorting BAM file')
    pysam.sort('-n', '-m', '1G', '-@', '1', '-o', str(MOTUS_PARAMETERS.get_alignment_file()),  str(MOTUS_PARAMETERS.get_temporary_alignment_file()))
    logging.info('Finished sorting BAM file')

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
        if False: 
            for mg, count in mg_2_raw_counts.items():
                numerator = float(count) / float(MOTUS_DB.get_length_by_mg(mg))
                norm_count = numerator / denominator
                scaled_count = norm_count * tot_cnt
                scaled_mg_2_counts[mg] = scaled_count
                norm_mg_2_counts[mg] = norm_count
        else:       # norm is actually normalisation by gene length
            for mg, count in mg_2_raw_counts.items():
                numerator = float(count) / float(MOTUS_DB.get_length_by_mg(mg))
                norm_count = numerator
                scaled_count = (norm_count/denominator) * tot_cnt
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



    def check_validity_of_bam_file(self, pg_entries: List[Dict[str, str]]) -> None:
        motus_found = False
        database_tool = None
        minlength = None

        for pg_entry in pg_entries:
            bamid = pg_entry.get('ID', None)
            if bamid == 'mOTUs4':
                motus_found = True
                database_tool = pg_entry.get('VN', None)
                minlength = pg_entry.get('CL', '')
                if '-l' in minlength:
                    minlength = int(minlength.split('-l')[-1].strip())
        if motus_found:
            if MOTUS_DB.get_full_version() != database_tool:
                logging.error('Version of BAM file and databases don\'t match')
                logging.error(f'BAM: {database_tool}')
                logging.error(f'Database/Tool: {MOTUS_DB.get_full_version()}')
                mutils.shutdown(1)
            if minlength > MOTUS_PARAMETERS.get_minimal_alignment_length():
                logging.info(f'Alignments in BAM have a minimum length of {minlength} which is above '
                             f'the provided minimum alignment length of {MOTUS_PARAMETERS.get_minimal_alignment_length()}'
                             f'. Increase minimum alignment length to continue.')
                mutils.shutdown(1)
        else:
            logging.error('BAM file invalid as it was not generated with mOTUs4')
            mutils.shutdown(1)



    def _bam_insert_iterator(self) -> Generator[Tuple[str, Dict[str, List[pysam.AlignedSegment]]], None, None]:
        """Reads through a sorted BAM file and
        finds the best alignment(s) per insert
        """
        alignments = pysam.AlignmentFile(MOTUS_PARAMETERS.get_alignment_file(), 'r')
        pg_entries = alignments.header.get('PG', [])
        self.check_validity_of_bam_file(pg_entries)
        #logging.warning('mOTUs tool/database have changed and bam file is invalid. Lenient mode enabled, will continue but results might be broken ...')

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
        header_line = mutils.create_mgc_header_line()
        #header_line = f'#tool_version={MOTUS_DB.get_tool_version()}\tdatabase_version={MOTUS_DB.get_database_version()}\tmin_alignment_length={MOTUS_PARAMETERS.get_minimal_alignment_length()}'
        handle.write(f'{header_line}\n')
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
#    has_header = False
    with open(mgc_file) as handle:
        header_line = handle.readline().strip()
        mutils.check_validity_of_mgc_header(header_line)
    mgc_2_count = {}
    count_mode = MOTUS_PARAMETERS.get_count_mode()
    with open(mgc_file) as handle:
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
    counts_smf = mentities.SinglemOTUsFile(motu_counts,MOTUS_PARAMETERS.get_minimal_alignment_length(), MOTUS_PARAMETERS.get_min_mgcs(), count_mode, MOTUS_DB.get_database_version(), MOTUS_DB.get_tool_version(), 'counts', MOTUS_PARAMETERS.get_sample_name())
    counts_smf.write_to_file(MOTUS_PARAMETERS.get_motu_file())
    if MOTUS_PARAMETERS._write_relative_abundances:
        relab_smf = counts_smf.get_relative_abundances()
        relab_smf.write_to_file(MOTUS_PARAMETERS.get_motu_file_relab())
    return None


class CapitalisedHelpFormatter(argparse.HelpFormatter):
    def add_usage(self, usage, actions, groups, prefix=None):
        if prefix is None:
            prefix = ''
        return super(CapitalisedHelpFormatter, self).add_usage(usage, actions, groups, prefix)


def parse_map_tax():
    parser = argparse.ArgumentParser(usage = f'''Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: {mutils.MOTUS_VERSION}
    
    {mutils.cite_text()}
        
    Summary:
        The map_tax command takes short read metagenomic sequencing data as input and
        maps reads to the mOTUs marker gene database.


    Usage:
        motus map_tax -f FILE [FILE ...] -r FILE [FILE ...] -s FILE [FILE ...] -o FILE [options]
        motus map_tax -f FILE [FILE ...] -r FILE [FILE ...] -o FILE [options]
        motus map_tax -s FILE [FILE ...] -o FILE [options]


    Input options:
        -f, --forward  FILE [FILE ...]
            Input file(s) for reads in forward orientation, fastQ/A(.gz)-formatted

        -r, --reverse  FILE [FILE ...]
            Input file(s) for reads in reverse orientation, fastQ/A(.gz)-formatted

        -s, --single  FILE [FILE ...]
            Input file(s) for unpaired reads, fastQ/A(.gz)-formatted

    Output options:
        -o, --output-file  FILE
            Output file name [required]

    Algorithm options:
        -l, --alignment-length  INT
            Minimum length of the alignment (bp) (default: 75)

        -t, --threads  INT
            Number of threads (default: 1)

    Database options:
        -db, --database  DIR
            Path to the mOTUs marker gene database directory (default: installation directory)
          ''', formatter_class=CapitalisedHelpFormatter,add_help=False)

    # Input options
    parser.add_argument("-f", "--forward", nargs="+",default=[], dest='f')  # input files(s) for reads in forward orientation, fastq(.gz)-formatted
    parser.add_argument("-r", "--reverse", nargs="+",default=[], dest='r')  # input files(s) for reads in reverse orientation, fastq(.gz)-formatted
    parser.add_argument("-s", "--single", nargs="+", default=[], dest='s')  # input files(s) for unpaired reads, fastq(.gz)-formatted

    # Output options
    parser.add_argument("-o", "--output-file", required=True, dest='o')  # output file name
    #parser.add_argument("-b", action="store_true")  # save the result in BAM format

    # ALgorithm options
    parser.add_argument("-l", "--alignment-length",  type=int, default=75, dest='l')  # min length of the alignment (bp) [75]
    parser.add_argument("-t", "--threads", type=int, default=1, dest='t')  # number of threads

    # Database options (added support for -db/--database)
    parser.add_argument("-db", "--database", type=str, default=None, dest='db')

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
    # Database location
    db_location = pathlib.Path(args.db) if args.db else mutils.DEFAULT_MOTUS_MGDB_LOCATION

    MOTUS_DB.load_motus_db(db_location)



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

    {mutils.cite_text()}

    motus batch_profile [options]

    Input options:
       -f  FILE[ FILE]  input tsv file. First column=samplename, second column=bam file

    Algorithm options:
       -g  INT          number of marker genes cutoff: 1=higher recall, 6=higher precision [3]
       -l  INT          min length of the alignment (bp) [75]
       -y  STR          type of read counts [INSERT_SCALED]
                        Values: [INSERT_RAW, INSERT_NORM, INSERT_SCALED, BASE_RAW, BASE_NORM]

    Database options:
       -db  DIR         path to the mOTUs marker gene database directory (default: installation directory)
    ]''', formatter_class=CapitalisedHelpFormatter, add_help=False)

    # Input options
    parser.add_argument("-f", required=True)  # input file(s) for reads in forward direction

    parser.add_argument("-g", type=int, default=3,
                        choices=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10])  # number of marker genes cutoff
    parser.add_argument("-l", type=int, default=75)  # min length of the alignment (bp) [75]
    parser.add_argument("-t", type=int, default=1)  # number of thread [1]
    parser.add_argument("-y", type=str, default='INSERT_SCALED',
                        choices=['INSERT_RAW', 'INSERT_NORM', 'INSERT_SCALED', 'BASE_RAW', 'BASE_NORM'])
    parser.add_argument("-db", "--database", type=str, default=None, dest='db')

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
    db_location = pathlib.Path(args.db) if args.db else mutils.DEFAULT_MOTUS_MGDB_LOCATION

    MOTUS_DB.load_motus_db(db_location)

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
    
    {mutils.cite_text()}
        
    Summary:
        The profile command in mOTUs is the main function that executes map_tax, calc_mgc,
        and calc_motu in sequence. It takes short read metagenomic sequencing data as input
        and generates a taxonomic profile.


    Usage:
       motus profile -f FILE [FILE ...] -r FILE [FILE ...] -s FILE [FILE ...] -o FILE [options]
       motus profile -f FILE [FILE ...] -r FILE [FILE ...] -o FILE [options]
       motus profile -s FILE [FILE ...] -o FILE [options]


    Input options:
        -f, --forward  FILE [FILE ...]
            Input file(s) for reads in forward orientation, fastQ/A(.gz)-formatted

        -r, --reverse  FILE [FILE ...]
            Input file(s) for reads in reverse orientation, fastQ/A(.gz)-formatted

        -s, --single  FILE [FILE ...]
            Input file(s) for unpaired reads, fastQ/A(.gz)-formatted

        -n, --sample-name  STR
            Sample name (default: 'unnamed sample')

    Output options:
        -o, --output-file  FILE
            Output file name [required]

    Algorithm options:
        -g, --marker-genes  INT
            Required number of marker genes for a mOTU to be called present: 
            1=higher recall, 6=higher precision, 10=maximum (default: 3)

        -l, --alignment-length  INT
            Minimum length of the alignment (bp) (default: 75)

        -t, --threads  INT
            Number of threads (default: 1)

        -y, --counting-mode  STR
            Which scale the abundances are reported in (default: INSERT_SCALED)
            Choices: [INSERT_RAW, INSERT_NORM, INSERT_SCALED, BASE_RAW, BASE_NORM]

    Database options:
        -db, --database  DIR
            Path to the mOTUs marker gene database directory (default: installation directory)
    ''', formatter_class=CapitalisedHelpFormatter,add_help=False)

    # Input options
    parser.add_argument("-f", "--forward", nargs="+", default=[], dest='f')  # input file(s) for reads in forward direction
    parser.add_argument("-r", "--reverse", nargs="+", default=[], dest='r')  # input file(s) for reads in reverse direction
    parser.add_argument("-s", "--single", nargs="+", default=[], dest='s')  # input file(s) for unpaired reads
    parser.add_argument("-n", "--sample-name", type=str, default='unnamed sample', dest='n')  # sample name

    # Output options
    parser.add_argument("-o", "--output-file", required=True, dest='o')
    parser.add_argument("-g", "--marker-genes", type=int, default=3, choices=[1,2,3,4,5,6,7,8,9,10], dest='g')  # number of marker genes cutoff
    parser.add_argument("-l", "--alignment-length", type=int, default=75, dest='l')  # min length of the alignment (bp) [75]
    parser.add_argument("-t", "--threads",  type=int, default=1, dest='t')  # number of thread [1]
    parser.add_argument("-y", "--counting-mode", type=str, default='INSERT_SCALED', choices=['INSERT_RAW', 'INSERT_NORM', 'INSERT_SCALED', 'BASE_RAW', 'BASE_NORM'], dest='y')
    #parser.add_argument("-c", action="store_true", help="Write second output file with relative abundances")

    # Database options
    parser.add_argument("-db", "--database", type=str, default=None, dest='db')

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
    db_location = pathlib.Path(args.db) if args.db else mutils.DEFAULT_MOTUS_MGDB_LOCATION

    MOTUS_DB.load_motus_db(db_location)


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
    MOTUS_PARAMETERS.set_write_relabundances()

    map_tax()
    calc_mgc()
    calc_motu()
    mutils.shutdown(0)


def parse_calc_mgc():
    parser = argparse.ArgumentParser(usage = f'''Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: {mutils.MOTUS_VERSION}
    
    {mutils.cite_text()}
        
    Summary:
        The calc_mgc command takes a file storing the alignments of sequencing reads
        to the mOTUs marker gene database and calculates marker gene cluster abundances.


    Usage:
        motus calc_mgc -i FILE -o FILE [options]


    Input options:
        -i, --input-file  FILE
            Path to BAM file generated after running the motus map_tax command [required]

    Output options:
        -o, --output-file  FILE
            Output file name [required]

    Algorithm options:
        -l, --alignment-length  INT
            Minimum length of the alignment (bp) (default: 75)

    Database options:
        -db, --database  DIR
            Path to the mOTUs marker gene database directory (default: installation directory)
       ''', formatter_class=CapitalisedHelpFormatter,add_help=False)

    parser.add_argument("-i", "--input-file", type=str, required=True, dest='i')  # provide a SAM or BAM input file (or list of files) output of motus map_tax
    parser.add_argument("-o", "--output-file", required=True, dest='o')  # output file name [stdout]
    parser.add_argument("-l", "--alignment-length", type=int, default=75, dest='l')  # min length of the alignment (bp) [75]
    parser.add_argument("-db", "--database", type=str, default=None, dest='db')

    args = parser.parse_args(sys.argv[2:])
    # print usage and exit if no arguments are passed
    if sys.argv[2:] == []:
        parser.print_usage()
        mutils.shutdown(1)

    alignment_file = pathlib.Path(args.i)
    mgc_file = pathlib.Path(args.o)
    inserts_file = pathlib.Path(args.o + '.inserts.gz')


    mutils.startup()
    min_alignment_length = args.l
    db_location = pathlib.Path(args.db) if args.db else mutils.DEFAULT_MOTUS_MGDB_LOCATION

    MOTUS_DB.load_motus_db(db_location)

    MOTUS_PARAMETERS.set_alignment_file(alignment_file, required_to_exist=True)
    MOTUS_PARAMETERS.set_mgc_file(mgc_file,required_to_exist=False)
    MOTUS_PARAMETERS.set_inserts_file(inserts_file, required_to_exist=False)
    MOTUS_PARAMETERS.set_minimal_alignment_length(min_alignment_length)
    MOTUS_PARAMETERS.set_threads(1)
    calc_mgc()
    mutils.shutdown(0)


def merge_profiles(motus_file_paths: List[pathlib.Path], output_motus_file_path: pathlib.Path) -> None:
    merged_motu_file = mentities.MergedmOTUsFile(motus_file_paths)
    merged_motu_file.write_to_file(output_motus_file_path)



def parse_merge():
    parser = argparse.ArgumentParser(usage=f'''Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: {mutils.MOTUS_VERSION}
    
    {mutils.cite_text()}

    Summary:
        The merge command takes multiple profiles produced after running the
        profile command and combines them into a single table.


    Usage:
        motus merge -i FILE [FILE ...] -o FILE


    Input options:
        -i, --input-files  FILE [FILE ...]
            A list of mOTUs profile files or a text file containing the list of profile
            files to be merged, with one line per file [required]

    Output options:
        -o, --output-file  FILE
            Output file name [required]

    Database options:
        -db, --database  DIR
            Path to the mOTUs marker gene database directory (default: installation directory)
          ''', formatter_class=CapitalisedHelpFormatter, add_help=False)


    parser.add_argument("-i", "--input-files", nargs="+", required=True, dest='i')
    parser.add_argument("-o", "--output-file", required=True, dest='o')
    parser.add_argument("-db", "--database", type=str, default=None, dest='db')
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
    db_location = pathlib.Path(args.db) if args.db else mutils.DEFAULT_MOTUS_MGDB_LOCATION

    MOTUS_DB.load_motus_db(db_location)

    merge_profiles(input_mOTUs_files, output_mOTUs_file)

    mutils.shutdown(0)


def download_genomes(genomes_to_download: List[str], genome_locator: mentities.GenomeLocator, output_folder: pathlib.Path, download_representative_genomes_only: bool) -> None:
    """Takes a list of genomes and downloads them to the 
    output folder. Only start downloading if all genome
    names are valid. Otherwise program will terminate
    with an error

    Args:
        genomes_to_download (List[str]): A list of genome names (1-n)
        genome_locator (GenomeLocator): Object which maps genome names to URLs 
        output_folder (pathlib.Path): The output folder. Will be created if it doesnt exist
        download_representative_genomes_only (bool): Only download representative genomes if set to True
    Returns:
        None None: None
    """    
    
    genomes_to_download = set(genomes_to_download)
    logging.info(f'Checking existence of {len(genomes_to_download)} genome names.')
    

    filtered_genomes_to_download = set()
    if download_representative_genomes_only:
        for genome in genomes_to_download:
            if genome_locator.is_represenative_genome(genome):
                filtered_genomes_to_download.add(genome)
        logging.info(f'Filtered for representative genomes. Remaining genomes: {len(filtered_genomes_to_download)}')
    else:
        filtered_genomes_to_download = genomes_to_download


    if len(filtered_genomes_to_download) == 0:
        logging.error('No genomes to download. Quitting')
        mutils.shutdown(1)

    genome_2_url = {}
    for genome in filtered_genomes_to_download:
        genome_2_url[genome] = genome_locator.get_genome_path(genome)
    
    logging.info(f'Downloading {len(genome_2_url)} genomes to {output_folder}')

    output_folder.mkdir(exist_ok=True, parents=True)

    for cnt, (genome, genome_path) in enumerate(genome_2_url.items(), 1):
        destpath = str(output_folder) + '/' + str(genome_path).split('/')[-1]
        logging.info(f'Downloading genome ({cnt} / {len(genome_2_url)}) {genome} to {destpath}')
        urllib.request.urlretrieve(genome_path, destpath)
    logging.info(f'Finished downloading genomes')


def prep_long(input_sequence_file: pathlib.Path, output_sequence_file: pathlib.Path, minlength: int = 50, split_length: int = 300) -> None:
    logging.info('Starting mOTUs - prep_long')
    logging.info(f'Input file: {input_sequence_file}')
    logging.info(f'Output file: {output_sequence_file}')

    total_bases_written = 0
    total_bases_seen = 0
    total_long_reads = 0
    total_short_reads = 0
    with open(output_sequence_file, 'w') as outhandle:
        for header, long_sequence in mutils.yield_reads(input_sequence_file):
            total_long_reads += 1
            total_bases_seen += len(long_sequence)
            substrings = [long_sequence[i:i + split_length - 1] for i in range(0, len(long_sequence), split_length)]
            fill_length = len(str(len(substrings))) + 1
            for cnt, substring in enumerate(substrings, 1):
                if len(substring) >= minlength:
                    header_st = f'{header}_st-{str(cnt).zfill(fill_length)}'
                    outhandle.write(f'>{header_st}\n{substring}\n')
                    total_short_reads += 1
                    total_bases_written += len(substring)
    logging.info(f'{total_long_reads:,} long reads, split into {total_short_reads:,} short reads')
    logging.info(f'Long reads had {total_bases_seen:,} bases. {total_bases_seen - total_bases_written:,} ({round(100.0 * (total_bases_seen - total_bases_written) / (total_bases_seen), 2)}%) bases were removed due to minimum length cutoff.')




def parse_prep_long():
    parser = argparse.ArgumentParser(usage=f'''Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: {mutils.MOTUS_VERSION}

    {mutils.cite_text()}

    Summary:
        The prep_long command takes long-read sequencing data and converts it
        into the appropriate input format to be used by the profile and map_tax commands.


    Usage:
        motus prep_long -i FILE -o FILE [options]


    Input options:
        -i, --input-file  FILE
            Long-read sequencing file to convert, can be in fastQ/A(.gz) format [required]

    Output options:
        -o, --output-file  FILE
            Output file name. This converted file is ready to be used by motus profile [required]

    Algorithm options:
        -sl, --splitting-length  INT
            Target fragment length (in bp) for splitting long reads (default: 300)

        -ml, --minimum-length  INT
            Minimum read length after splitting. Shorter reads are discarded (default: 50)
           ''', formatter_class=CapitalisedHelpFormatter, add_help=False)


    parser.add_argument("-i", "--input-file",  required=True, dest='i')
    parser.add_argument("-o", "--output-file", required=True, dest='o')
    parser.add_argument("-sl", "--splitting-length", default=300, type=int, dest='sl')
    parser.add_argument("-ml", "--minimum-length", default=50, type=int, dest='ml')
    args = parser.parse_args(sys.argv[2:])

    if sys.argv[2:] == []:
        parser.print_usage()
        mutils.shutdown(1)
    mutils.startup()

    input_sequence_file = pathlib.Path(args.i)
    output_sequence_file = pathlib.Path(args.o)
    minlength = args.ml
    split_length = args.sl

    prep_long(input_sequence_file, output_sequence_file, minlength=minlength, split_length=split_length)


def parse_find():
    parser = argparse.ArgumentParser(usage=f'''Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: {mutils.MOTUS_VERSION}

    {mutils.cite_text()}

    Summary:
        The genomes command queries the mOTUs-db based on identifiers, functional,
        or taxonomic annotations and returns a list of genomes matching indicated query.


    Usage:    
        motus genomes -i FILE -o FILE [options]
        motus genomes -i STR [STR ...] -o FILE [options]


    Input options:
        -i, --input-queries  FILE/STR
            Can be either a list of search queries or a text file listing search queries
            with one line per query. Queries can be genome or mOTUs identifiers, PFAM, KEGG, EGGNOG, 
            or GTDB taxonomy names. If the query does not exactly match any database entry,
            alternative queries will be suggested [required]

    Output options:
        -o, --output-file  FILE
            Output file containing a list of genome identifiers matching search queries and their 
            annotations as indicated by the -d parameter. This output file can be used as input
            for the motus download command [required]

        -d, --details  STR [STR ...]
            List of annotations to report. Choose any combination of [KEGG, PFAM, EGGNOG, TAXONOMY],
            for example, -d KEGG PFAM.
           ''', formatter_class=CapitalisedHelpFormatter, add_help=False)


    parser.add_argument("-i", "--input-queries", required=True, nargs="+", dest='i')
    parser.add_argument("-o", "--output-file", required=True, dest='o')
    parser.add_argument("-d", "--details", default=[], nargs="+", dest='d')

    args = parser.parse_args(sys.argv[2:])

    if sys.argv[2:] == []:
        parser.print_usage()
        mutils.shutdown(1)
    mutils.startup()


    output_file = pathlib.Path(args.o)
    search_queries_tmp = args.i
    annotations_to_report_tmp = set(args.d)
    allowed_annotations = ['KEGG', 'EGGNOG', 'PFAM', 'TAXONOMY']
    annotations_to_report = []
    for annotation in annotations_to_report_tmp:
        if len(annotation) == 0:
            continue
        if annotation not in allowed_annotations:
            logging.error(f'Unknown annotation to report: {annotation}. Allowed annotations: {allowed_annotations}')
            mutils.shutdown(1)
        annotations_to_report.append(annotation)

    search_queries = []
    if len(search_queries_tmp) == 1: # can be a file or a query
        search_query = search_queries_tmp[0]
        if pathlib.Path(search_query).exists(): # is a file with search queries
            with open(search_query) as handle:
                for line in handle:
                    search_queries.append(line.strip())
        else:
            search_queries.append(search_query)
    else: # list of genomes
        search_queries = search_queries_tmp

    mfind.find_genomes(search_queries, output_file, annotations_to_report)
    mutils.shutdown(0)


def parse_classify():
    parser = argparse.ArgumentParser(usage=f'''Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: {mutils.MOTUS_VERSION}

    {mutils.cite_text()}

    Summary:
        The classify command takes a list of genome sequence files as input and
        assigns these genomes to existing mOTUs in the database.


    Usage:
        motus classify -i FILE -o FILE [options]


    Input options:
        -i, --input-file  FILE
            Text file listing genome sequence files in fastA(.gz) format to classify.
            One line per genome file [required]

    Output options:
        -o, --output-file  FILE
            Output file name. Each line contains a genome and its associated mOTU [required]

    Algorithm options:
        -t, --threads  INT
            Number of threads (default: 1)

    Database options:
        -db, --database  DIR
            Path to the mOTUs marker gene database directory (default: installation directory)
           ''', formatter_class=CapitalisedHelpFormatter, add_help=False)


    parser.add_argument("-i", "--input-file", required=True, dest='i')
    parser.add_argument("-o", "--output-file", required=True, dest='o')
    parser.add_argument("-t", "--threads", default=1, type=int, dest='t')
    parser.add_argument("-db", "--database", type=str, default=None, dest='db')
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
    db_location = pathlib.Path(args.db) if args.db else mutils.DEFAULT_MOTUS_MGDB_LOCATION

    classify(genome_files, output_file, args.t, db_location)
    mutils.shutdown(0)



def classify(genome_files: List[pathlib.Path], output_file: pathlib.Path, threads: int = 1, db_location: pathlib.Path = None):
    """Takes a list of genome files and associates them with existing mOTUs.
    A genome can either be:
    - classified with a mOTU (=mOTUXXX)
    - not have enough markergenes to be classified (=notEnoughMGs)
    - have enough markergenes to be classified but can not be represented by an existing mOTU (=Novel)

    Args:
        genome_files (List[pathlib.Path]): A list with pathlib.Path objects all pointing to an existing genome file
        output_file (pathlib.Path): Main output file of motus classify
        threads (int, optional): Number of threads. Defaults to 1.
        db_location (pathlib.Path, optional): Path to the mOTUs database directory. Defaults to installation directory.
    """

    if db_location is None:
        db_location = mutils.DEFAULT_MOTUS_MGDB_LOCATION
    genome_files = sorted(genome_files)
    root_tmp_folder =  pathlib.Path(str(output_file) + '_classify_tmp')


    logging.info(f'Starting mOTUs classify:')
    MOTUS_DB.load_motus_db(db_location, True)
    logging.info(f'\tInput = {len(genome_files)} genomes.')
    logging.info(f'\tOutput will be written to {output_file}')
    logging.info(f'\tTemporary files will be written to {root_tmp_folder}')
    logging.info(f'\tRunning fetchMGs on genomes.')

    
    fetchmgs_tmp_folder = root_tmp_folder.joinpath('fetchmgs')
    fetchmgs_marker = root_tmp_folder.joinpath('fetchmgs.done')
    fetchmgs_fna = root_tmp_folder.joinpath('motus_classify.fna')
    fetchmgs_tsv = root_tmp_folder.joinpath('motus_classify.tsv')
    fetchmgs_genomes_done = []
    if fetchmgs_marker.exists():
        with open(fetchmgs_marker) as handle:
            for line in handle:
                fetchmgs_genomes_done.append(line.strip())


    if set(fetchmgs_genomes_done) != set([str(x) for x in genome_files]):
        from fetchmgs import fetchmgs
        fetchmgs.extraction_genomes(genome_files, fetchmgs_tmp_folder, 'genome', threads, True) # TODO rewrite with a Processpool to accelerate for larger number of genomes
        logging.info(f'\tFinished running fetchMGs on genomes.')
        logging.info(f'\tCollecting fetchMGs results')
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
        with open(fetchmgs_fna, 'w') as seqhandle, open(fetchmgs_tsv, 'w') as tsvhandle:
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
        with open(fetchmgs_marker, 'w') as outhandle:
            for genome_file in genome_files:
                outhandle.write(f'{genome_file}\n')
        logging.info(f'\tFinished collecting fetchMGs results. Genomes = {len(genome_2_markergenes)}, Genomes with enough MGs = {len(genome_2_markergenes) - genomes_removed_notenoughmgs}')
    else:
        logging.info('Reusing fetchMGs from previous run')


    logging.info(f'\tMatching genomes against the mOTUs database')
    alignment_m8 = root_tmp_folder.joinpath('motus_classify.align_vs_motu.m8')
    combined_distances_file = root_tmp_folder.joinpath('motus_classify.cd_vs_motu.tsv')
    if True:  # if the marker file with all genomes exists?
        logging.info(f'\tAligning genome marker genes against the mOTUs marker gene database using vsearch')
        vsearch_command = f'vsearch --threads {threads} --usearch_global {str(fetchmgs_fna)} --db {MOTUS_DB.get_bwa_index()} --strand both --id 0.8 --maxaccepts 2000 --maxrejects 2000 --mincols 20 --userout {str(alignment_m8)} --userfields query+target+id+alnlen+mism+ids+ql+tl --mincols 40'
        try:
            subprocess.run(vsearch_command,shell=True, stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Command {vsearch_command} failed") from e
        logging.info(f'\tFinished alignment')

        logging.info(f'\tCombining individual marker gene distances into genome to genome distances')
        query_genome_2_unfilt_hits = collections.defaultdict(list)
        with open(alignment_m8) as handle:
            for line in handle:
                [query, target, percid, alignment_length, mismatches, matches, query_length, target_length] = line.strip().split()
                [genome, cog] = query.rsplit('.', 1)
                if 'unassigned' in target:
                    continue
                query_genome_2_unfilt_hits[genome].append([cog, target, float(percid), int(alignment_length), int(mismatches), int(matches), int(query_length), int(target_length)])
        
        query_genome_2_cogs = collections.defaultdict(dict)
        with open(fetchmgs_fna) as handle:
            for (header, sequence) in FastaIO.SimpleFastaParser(handle):
                [genome, cog] = header.rsplit('.', 1)
                query_genome_2_cogs[genome][cog] = len(sequence)
        motu_2_cogs = MOTUS_DB.get_motu_2_median_mgc_gene_length()
        genome_2_motus = {}
        for query_genome, query_cog_2_length in query_genome_2_cogs.items():
            alignment_unfiltered_hits = query_genome_2_unfilt_hits[query_genome]
            motu_2_cog_2_alignments = collections.defaultdict(lambda : collections.defaultdict(list))
            for [query_cog, motu_mg, percid, alignment_length, mismatches, matches, query_length, target_length] in alignment_unfiltered_hits:
                motu = motu_mg.rsplit('.', 2)[0]
                motu_cog = motu_mg.rsplit('.', 2)[1].split('_')[1]
                if motu_cog != query_cog:
                    continue
                motu_2_cog_2_alignments[motu][motu_cog].append([percid, alignment_length, mismatches, matches, query_length, target_length])
            motu_2_cog_2_bestalignment = {}
            for motu, cog_2_alignments in motu_2_cog_2_alignments.items():
                cog_2_best_aln = {}
                for cog, alignments in cog_2_alignments.items():
                    most_matches = max([x[3] for x in alignments])
                    best_alignment = [x for x in alignments if x[3] == most_matches][0]
                    # remove anything with an alignment length < 80%
                    alignment_length = best_alignment[1]
                    query_length = best_alignment[-2]
                    target_length = best_alignment[-1]
                    shorter_sequence_length = query_length
                    if target_length < query_length:
                        shorter_sequence_length = target_length
                    alignment_coverage = alignment_length * 100.0 / shorter_sequence_length
                    if alignment_coverage >= 80.0:
                        cog_2_best_aln[cog] = [best_alignment[3], shorter_sequence_length, True]
                motu_2_cog_2_bestalignment[motu] = cog_2_best_aln
            # top up alignments by fake alignments
            query_genome_cogs = query_genome_2_cogs[query_genome]
            motu_2_combined_distance = {}
            for motu, cog_2_bestalignment in motu_2_cog_2_bestalignment.items():
                if len(cog_2_bestalignment) < 6:
                    continue
                motu_cogs = motu_2_cogs[motu]
                combined_cogs = collections.Counter(list(motu_cogs.keys()) + list(query_genome_cogs.keys()))
                for cog, count in combined_cogs.items():
                    if count == 1:
                        continue
                    if cog in cog_2_bestalignment:
                        continue
                    motu_cog_length = motu_cogs[cog]
                    genome_cog_length = query_genome_cogs[cog]
                    shorter_sequence_length = motu_cog_length
                    if genome_cog_length< shorter_sequence_length:
                        shorter_sequence_length = genome_cog_length
                    matches = int(shorter_sequence_length * 0.8)
                    cog_2_bestalignment[cog] = [matches, shorter_sequence_length, False]
                tot_matches = 0
                tot_length = 0
                for cog, best_aln in cog_2_bestalignment.items():
                    tot_matches += best_aln[0]
                    tot_length += best_aln[1]
                combined_distance = tot_matches * 100.0 / tot_length
                if combined_distance >= 96.5:
                    motu_2_combined_distance[motu] = combined_distance
            if len(motu_2_combined_distance) == 0:
                motu_2_combined_distance['Unknown'] = -1.0
            genome_2_motus[query_genome] = motu_2_combined_distance
        
        with open(combined_distances_file, 'w') as handle:
            for genome, motus in genome_2_motus.items():
                for (motu, dist) in motus.items():
                    handle.write(f'{genome}\t{motu}\t{round(dist, 2)}\n')
        logging.info(f'\tFinished combining')
    

    genome_2_motus = collections.defaultdict(list)
    with open(combined_distances_file) as handle:
        for line in handle:
            [genome, motu, dist] = line.strip().split('\t')
            dist = float(dist)
            genome_2_motus[genome].append((motu, dist))
    
    genome_2_best_motu = {}
    for genome, motus in genome_2_motus.items():
        best_dist = max([x[1] for x in motus])
        best_motu = sorted([x for x in motus if x[1] == best_dist])[0]
        genome_2_best_motu[genome] = best_motu #(motu, dist)
    with open(fetchmgs_tsv) as handle, open(output_file, 'w') as outhandle:
        outhandle.write('GENOME\tMOTU\tSIMILARITY\tNUM_MGS\n')
        handle.readline()
        for line in handle:
            [genome, num_mgs, mgs] = line.strip().split('\t')
            (motu, dist) = genome_2_best_motu.get(genome, ('<6MGs-no_mOTU', '-1'))
            if motu == 'Unknown':
                motu = 'Novel-no_mOTU'
            tmp = '\t'.join([genome, motu, str(dist), num_mgs])
            outhandle.write(f'{tmp}\n')
            
        
    
    







def parse_downloadDB():
    parser = argparse.ArgumentParser(usage=f'''Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: {mutils.MOTUS_VERSION}

    {mutils.cite_text()}

    Summary:
        The downloadMGDB command downloads the marker gene reference database used
        by the profile and map_tax commands.


    Usage:
        motus downloadMGDB [options]


    Options:
        -f, --force
            Force download even when database is already present

        -db, --database  DIR
            Path to the directory where the database will be downloaded (default: installation directory).
            The database will be placed in a subdirectory named 'db_mOTU' inside this path.
           ''', formatter_class=CapitalisedHelpFormatter, add_help=False)

    parser.add_argument("-f", "--force", action="store_true", dest='f')
    parser.add_argument("-db", "--database", type=str, default=None, dest='db')
    args = parser.parse_args(sys.argv[2:])

    force_download = False
    if args.f:
        force_download = True

    mutils.startup()

    if args.db:
        db_parent_location = pathlib.Path(args.db)
        db_location = db_parent_location.joinpath('db_mOTU')
    else:
        db_parent_location = mutils.DEFAULT_MOTUS_MGDB_PARENT_LOCATION
        db_location = mutils.DEFAULT_MOTUS_MGDB_LOCATION
    db_marker = db_location.joinpath('db_mOTU.downloaded')

    if db_marker.exists() and not force_download:
        logging.info('Database already downloaded and -f not set. All good.')
        mutils.shutdown(0)
    if db_marker.exists() and force_download:
        logging.info('Database already downloaded and -f set. Will delete current database and download again.')
        shutil.rmtree(db_location)

    db_parent_location.mkdir(parents=True, exist_ok=True)
    logging.info('Start downloading mOTUs marker gene database. ~6GB')
    dest_tar_gz_file = db_parent_location.joinpath('db_mOTU.tar.gz')
    if dest_tar_gz_file.is_file():
        dest_tar_gz_file.unlink()


    with urllib.request.urlopen(mutils.MOTUS_MGDB_REMOTE_LOCATION) as response:
        total = int(response.info().get("Content-Length", -1))
        with open(str(dest_tar_gz_file), "wb") as f, tqdm.tqdm(total=total, unit='B', unit_scale=True, desc='Downloading mOTUs marker gene database') as pbar:
            while True:
                chunk = response.read(8192)
                if not chunk:
                    break
                f.write(chunk)
                pbar.update(len(chunk))


    logging.info('Finished downloading mOTUs marker gene database.')
    logging.info('Start un-taring mOTUs marker gene database.')
    if db_location.exists():
        shutil.rmtree(db_location)

    with tarfile.open(dest_tar_gz_file, 'r') as t:
        t.extractall(db_parent_location)

    db_marker.touch(exist_ok=True)
    logging.info('Finished untaring mOTUs marker gene database.')
    mutils.shutdown(0)


def parse_download():
    parser = argparse.ArgumentParser(usage=f'''Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: {mutils.MOTUS_VERSION}
    
    {mutils.cite_text()}

    Summary:
        The download command downloads listed genome files from mOTUs-db.


    Usage:
        motus download -i FILE -o PATH [options]
        motus download -i STR [STR ...] -o PATH [options]


    Input options:
        -i, --input-genomes  FILE/STR
            Can be either a list of genome identifiers separated by spaces or a text file
            listing the identifiers of genomes for download. One line per genome. The output of
            the motus genomes command can be used as input for this command [required]

    Output options:
        -o, --output-folder  PATH
            Path to output folder where the downloaded sequences will be saved [required]

        -r, --representatives
            Download only sequences from representative genomes.

    Database options:
        -db, --database  DIR
            Path to the mOTUs marker gene database directory (default: installation directory)
           ''', formatter_class=CapitalisedHelpFormatter, add_help=False)


    parser.add_argument("-o", "--output-folder", required=True, dest='o')
    parser.add_argument("-i", "--input-genomes", required=True, nargs="+", dest='i')
    parser.add_argument("-r", "--representatives", action="store_true", dest='r')
    parser.add_argument("-db", "--database", type=str, default=None, dest='db')


    args = parser.parse_args(sys.argv[2:])
    if sys.argv[2:] == []:
        parser.print_usage()
        mutils.shutdown(1)
    mutils.startup()


    output_folder = pathlib.Path(args.o)
    input_items = args.i
    download_representative_genomes_only = False
    if args.r:
        download_representative_genomes_only = True

    genomes_to_download = []
    if len(input_items) == 1: # can be a single genome or a file with genomes
        input_item = input_items[0]
        if pathlib.Path(input_item).exists(): # is a file with genome names
            with open(input_item) as handle:
                for line in handle:
                    if line.strip().startswith('GENOME'):
                        continue
                    genomes_to_download.append(line.strip().split('\t')[0])
        else:
            genomes_to_download.append(input_item)
    else: # list of genomes
        genomes_to_download = input_items

    db_location = pathlib.Path(args.db) if args.db else mutils.DEFAULT_MOTUS_MGDB_LOCATION
    MOTUS_DB.load_motus_db(db_location, load=False)
    motus_search_db = mentities.GenomeLocator(MOTUS_DB.genome_metadata_file)
    download_genomes(genomes_to_download, motus_search_db, output_folder, download_representative_genomes_only)
    mutils.shutdown(0)





def parse_calc_motu():
    parser = argparse.ArgumentParser(usage = f'''Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: {mutils.MOTUS_VERSION}
    
    {mutils.cite_text()}
        
    Summary:
        The calc_motu command takes a file containing marker gene cluster
        abundances and generates a taxonomic profile.


    Usage:
        motus calc_motu -i FILE -o FILE [options]


    Input options:
        -i, --input-file  FILE
            MGC abundance table generated by the calc_mgc command [required]

        -n, --sample-name  STR
            Sample name (default: 'unnamed sample')

    Output options:
        -o, --output-file  FILE
            Output file name [required]

    Algorithm options:
        -g, --marker-genes  INT
            Required number of marker genes for a mOTU to be called present: 
            1=higher recall, 6=higher precision, 10=maximum (default: 3)

        -y, --counting-mode  STR
            Which scale the abundances are reported in (default: INSERT_SCALED)
            Choices: [INSERT_RAW, INSERT_NORM, INSERT_SCALED, BASE_RAW, BASE_NORM]

    Database options:
        -db, --database  DIR
            Path to the mOTUs marker gene database directory (default: installation directory)
          ''', formatter_class=CapitalisedHelpFormatter,add_help=False)


    parser.add_argument("-n", "--sample-name", type=str, default='unnamed sample', dest='n')  # sample name
    parser.add_argument("-i", "--input-file",  required=True, dest='i')  # provide the mgc abundance table(output of motus calc_mgc)
    parser.add_argument("-o", "--output-file", required=True, dest='o')  # output fil name [stdout]
    parser.add_argument("-y", "--counting-mode", type=str, default='INSERT_SCALED', choices=['INSERT_RAW', 'INSERT_NORM', 'INSERT_SCALED', 'BASE_RAW', 'BASE_NORM'], dest='y')
    parser.add_argument("-g", "--marker-genes", type=int, default=3, choices=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], dest='g')  # number of marker genes cutoff
    #parser.add_argument("-c", action="store_true",help="Write second output file with relative abundances")
    parser.add_argument("-db", "--database", type=str, default=None, dest='db')

    args = parser.parse_args(sys.argv[2:])
    if sys.argv[2:] == []:
        parser.print_usage()
        mutils.shutdown(1)


    mgc_file = pathlib.Path(args.i)
    motu_file = pathlib.Path(args.o)
    mutils.startup()
    samplename = args.n
    db_location = pathlib.Path(args.db) if args.db else mutils.DEFAULT_MOTUS_MGDB_LOCATION

    MOTUS_DB.load_motus_db(db_location)

    MOTUS_PARAMETERS.set_mgc_file(mgc_file, required_to_exist=True)
    MOTUS_PARAMETERS.set_motu_file(motu_file, required_to_exist=False)
    MOTUS_PARAMETERS.set_sample_name(samplename)
    MOTUS_PARAMETERS.set_threads(1)
    MOTUS_PARAMETERS.set_count_mode(args.y)
    MOTUS_PARAMETERS.set_minimal_number_of_mgcs(args.g)
    MOTUS_PARAMETERS.set_write_relabundances()

    calc_motu()
    mutils.shutdown(0)




if __name__ == "__main__":
    main()

def main():
    parser = argparse.ArgumentParser(usage = f'''Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: {mutils.MOTUS_VERSION}
    
    {mutils.cite_text()}
        
    Usage:
        motus <command> [options]


    Commands:

        -- Taxonomic profiling

            profile       Perform taxonomic profiling (map_tax + calc_mgc + calc_motu) in a single step

            map_tax       Map reads to the marker gene database
            calc_mgc      Calculate marker gene cluster (MGC) abundance
            calc_motu     Summarize MGC abundances into a mOTU profile


        -- Tool utilities

            downloadMGDB  Download the mOTUs marker gene database
            merge         Merge multiple taxonomic profiling results into one table
            classify      Classify user genomes into mOTUs
            prep_long     Prepare long reads to be profiled by mOTUs


        -- Genome accession

            genomes       Search the mOTUs-db by keyword (taxonomic, functional)
            download      Download sequence files from mOTUs-db
    
    
        Type motus <command> to print the help menu for a specific command
        ''',formatter_class=CapitalisedHelpFormatter,add_help=False)

    parser.add_argument('command', choices=["profile", "map_tax", "calc_mgc", "calc_motu", "download", "merge", "downloadMGDB", "batch_profile", "classify", 'prep_long', 'genomes'])
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
    elif args.command == 'downloadMGDB':
        parse_downloadDB()
    elif args.command == 'classify':
        parse_classify()
    elif args.command == 'prep_long':
        parse_prep_long()
    elif args.command == 'genomes':
        parse_find()
    else:
        parser.print_usage()
        print(f'Unrecognized command {args}')
        mutils.shutdown(1)
    mutils.shutdown(0)
