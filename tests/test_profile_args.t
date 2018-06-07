#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from simpletap.motus import MotusTestCase


class TestProfileArgs(MotusTestCase):
    def test_missing_fwd(self):
        "motus profile with a missing forward file"

        cmd = ("profile",
               "-r", "sample_data/ERR878216_sample_1.fastq.gz")

        stdout, stderr, exit = self.run_motus(cmd, allow_error=True)

        self.assertRegexpMatches(stderr, ".*Error: forward reads \(-f\) missing")

    def test_missing_rev(self):
        "motus profile with a missing reverse file"

        cmd = ("profile",
               "-f", "sample_data/ERR878216_sample_1.fastq.gz")

        stdout, stderr, exit = self.run_motus(cmd, allow_error=True)

        self.assertRegexpMatches(stderr, ".*Error: reverse reads \(-r\) missing")


if __name__ == "__main__":
    from simpletap import TAPTestRunner
    unittest.main(testRunner=TAPTestRunner())

# vim: ai sts=4 et sw=4
