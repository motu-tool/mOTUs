#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from simpletap.motus import MotusTestCase


class TestProfileBiom(MotusTestCase):
    def test_normal_fw_rev_biom(self):
        "motus profile with fwd+rev reads for BIOM format"

        cmd = ("profile",
               "-B",
               "-f", "sample_data/ERR878216_sample_1.fastq.gz",
               "-r", "sample_data/ERR878216_sample_2.fastq.gz")

        stdout, stderr, exit = self.run_motus(cmd)

        self.validate_result(stdout, "sample_data/ERR878216_both.biom")

    def test_normal_single_biom(self):
        "motus profile with single reads for BIOM format"

        cmd = ("profile",
               "-B",
               "-s", "sample_data/ERR878216_sample_1.fastq.gz")

        stdout, stderr, exit = self.run_motus(cmd)

        self.validate_result(stdout, "sample_data/ERR878216_single_1.biom")


if __name__ == "__main__":
    from simpletap import TAPTestRunner
    unittest.main(testRunner=TAPTestRunner())

# vim: ai sts=4 et sw=4
