#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import re
from simpletap.motus import MotusTestCase


class TestProfile(MotusTestCase):
    def test_normal_fw_rev(self):
        "motus profile with a forward and reverse file"

        cmd = ("profile",
               "-f", "sample_data/ERR878216_sample_1.fastq.gz",
               "-r", "sample_data/ERR878216_sample_2.fastq.gz")

        stdout, stderr, exit = self.run_motus(cmd)

        self.validate_result(stdout, "sample_data/ERR878216_both.profile")

    def test_normal_fw_rev_counts(self):
        "motus profile with a forward and reverse file to produce raw counts"

        cmd = ("profile",
               "-c",
               "-f", "sample_data/ERR878216_sample_1.fastq.gz",
               "-r", "sample_data/ERR878216_sample_2.fastq.gz")

        stdout, stderr, exit = self.run_motus(cmd)

        self.validate_result(stdout, "sample_data/ERR878216_both_counts.profile")

    def test_single(self):
        "motus profile with a single file"

        cmd = ("profile",
               "-s", "sample_data/ERR878216_sample_1.fastq.gz")

        stdout, stderr, exit = self.run_motus(cmd)

        self.validate_result(stdout, "sample_data/ERR878216_single_1.profile")

    def test_normal_multithread(self):
        "motus profile with multiple threads"

        cmd = ("profile",
               "-f", "sample_data/ERR878216_sample_1.fastq.gz",
               "-r", "sample_data/ERR878216_sample_2.fastq.gz",
               "-t", "8")

        stdout, stderr, exit = self.run_motus(cmd)

        self.validate_result(stdout, "sample_data/ERR878216_both.profile")

    def test_input_bam(self):
        "motus profile by providing a bam file"

        cmd = ("profile",
               "-i", "sample_data/ERR878216_both.bam")

        stdout, stderr, exit = self.run_motus(cmd, allow_error=True)

        self.validate_result(stdout, "sample_data/ERR878216_both.profile")

    def test_input_mgc(self):
        "motus profile by providing a mgc file"

        cmd = ("profile",
               "-m", "sample_data/ERR878216_both.mgc")

        stdout, stderr, exit = self.run_motus(cmd, allow_error=True)

        self.validate_result(stdout, "sample_data/ERR878216_both.profile")

    def test_input_sample_name(self):
        "motus and providing a sample name"

        SAMPLE = "MySuperLongSampleName"

        cmd = ("profile",
               "-n", SAMPLE,
               "-s", "sample_data/ERR878216_sample_1.fastq.gz")

        stdout, stderr, exit = self.run_motus(cmd, allow_error=True)

        self.assertRegexpMatches(stdout,
                                 re.compile("^#consensus_taxonomy\t{0}".format(SAMPLE), re.MULTILINE))

if __name__ == "__main__":
    from simpletap import TAPTestRunner
    unittest.main(testRunner=TAPTestRunner())

# vim: ai sts=4 et sw=4
