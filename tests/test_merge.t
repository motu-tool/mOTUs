#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import tempfile
import shutil
from simpletap.motus import MotusTestCase


class TestMerge(MotusTestCase):
    def test_normal_multiple(self):
        "motus merge with multiple profiles"

        cmd = ("merge",
               "-i", ",".join([
                   "sample_data/ERR878216_single_1.profile",
                   "sample_data/ERR878216_single_2.profile",
               ]))

        stdout, stderr, exit = self.run_motus(cmd)

        self.validate_result(stdout, "sample_data/ERR878216_single_merge_1_2.profile")

    def test_normal_single(self):
        "motus merge with single profile"

        cmd = ("merge",
               "-i", "sample_data/ERR878216_single_1.profile")

        stdout, stderr, exit = self.run_motus(cmd)

        self.validate_result(stdout, "sample_data/ERR878216_single_1.profile")


class TestSingleMergeDirectory(MotusTestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        shutil.copy("sample_data/ERR878216_single_1.profile", self.tmpdir)

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_normal_single(self):
        "motus merge with single profile in a directory"

        cmd = ("merge",
               "-d", self.tmpdir)

        stdout, stderr, exit = self.run_motus(cmd)

        self.validate_result(stdout, "sample_data/ERR878216_single_1.profile")

class TestMultipleMergeDirectory(MotusTestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        shutil.copy("sample_data/ERR878216_single_1.profile", self.tmpdir)
        shutil.copy("sample_data/ERR878216_single_2.profile", self.tmpdir)

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_normal_multiple(self):
        "motus merge with multiple profiles in a directory"

        cmd = ("merge",
               "-d", self.tmpdir)

        stdout, stderr, exit = self.run_motus(cmd)

        self.validate_result(stdout, "sample_data/ERR878216_single_merge_1_2.profile")


if __name__ == "__main__":
    from simpletap import TAPTestRunner
    unittest.main(testRunner=TAPTestRunner())

# vim: ai sts=4 et sw=4
