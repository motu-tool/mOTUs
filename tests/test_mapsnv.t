#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import unittest
import tempfile
from simpletap.motus import MotusTestCase


class TestMapSNV(MotusTestCase):
    def setUp(self):
        self.tmpbamfh, self.tmpbam = tempfile.mkstemp(".bam")

    def test_normal_fw_rev(self):
        "motus map_snv with a forward and reverse file"

        cmd = ("map_snv",
               "-f", "sample_data/ERR878216_sample_1.fastq.gz",
               "-r", "sample_data/ERR878216_sample_2.fastq.gz",
               "-o", self.tmpbam)

        stdout, stderr, exit = self.run_motus(cmd)

        self.validate_bam(self.tmpbam, "sample_data/ERR878216_both_snv.bam")

    def test_single(self):
        "motus map_snv with a single file"

        cmd = ("map_snv",
               "-s", "sample_data/ERR878216_sample_1.fastq.gz",
               "-o", self.tmpbam)

        stdout, stderr, exit = self.run_motus(cmd)

        self.validate_bam(self.tmpbam, "sample_data/ERR878216_single_1_snv.bam")

    def test_normal_multithread(self):
        "motus map_snv with multiple threads"

        cmd = ("map_snv",
               "-f", "sample_data/ERR878216_sample_1.fastq.gz",
               "-r", "sample_data/ERR878216_sample_2.fastq.gz",
               "-t", "8",
               "-o", self.tmpbam)

        stdout, stderr, exit = self.run_motus(cmd)

        self.validate_bam(self.tmpbam, "sample_data/ERR878216_both_snv.bam")

    def tearDown(self):
        try:
            self.tmpbamfh.close()
        except AttributeError:
            pass

        try:
            os.remove("sample_file")
        except (OSError, ValueError):
            pass

        try:
            os.remove(self.tmpbam)
        except (OSError, ValueError, AttributeError):
            pass


if __name__ == "__main__":
    from simpletap import TAPTestRunner
    unittest.main(testRunner=TAPTestRunner())

# vim: ai sts=4 et sw=4
