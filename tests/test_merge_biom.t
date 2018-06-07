#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import tempfile
import shutil
from simpletap.motus import MotusTestCase


class TestMergeBiom(MotusTestCase):
    def test_normal_multiple_biom(self):
        "motus merge with multiple profiles for BIOM format"

        cmd = ("merge",
               "-B",
               "-i", ",".join([
                   "sample_data/ERR878216_single_1.profile",
                   "sample_data/ERR878216_single_2.profile",
               ]))

        stdout, stderr, exit = self.run_motus(cmd)

        self.validate_result(stdout, "sample_data/ERR878216_single_merge_1_2.biom")

    def test_normal_single_biom(self):
        "motus merge with single profile for BIOM format"

        cmd = ("merge",
               "-B",
               "-i", "sample_data/ERR878216_single_1.profile")

        stdout, stderr, exit = self.run_motus(cmd)

        self.validate_result(stdout, "sample_data/ERR878216_single_1.biom")


class TestSingleMergeBiomDirectory(MotusTestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        shutil.copy("sample_data/ERR878216_single_1.profile", self.tmpdir)

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_normal_single_biom(self):
        "motus merge with single profile in a directory for BIOM format"

        cmd = ("merge",
               "-B",
               "-d", self.tmpdir)

        stdout, stderr, exit = self.run_motus(cmd)

        self.validate_result(stdout, "sample_data/ERR878216_single_1.biom")

class TestMultipleMergeBiomDirectory(MotusTestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        shutil.copy("sample_data/ERR878216_single_1.profile", self.tmpdir)
        shutil.copy("sample_data/ERR878216_single_2.profile", self.tmpdir)

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_normal_multiple_biom(self):
        "motus merge with multiple profiles in a directory for BIOM format"

        cmd = ("merge",
               "-B",
               "-d", self.tmpdir)

        stdout, stderr, exit = self.run_motus(cmd)

        self.validate_result(stdout, "sample_data/ERR878216_single_merge_1_2.biom")


if __name__ == "__main__":
    from simpletap import TAPTestRunner
    unittest.main(testRunner=TAPTestRunner())

# vim: ai sts=4 et sw=4
