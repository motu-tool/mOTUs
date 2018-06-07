#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from glob import glob
from subprocess import Popen, PIPE
from simpletap.motus import MotusTestCase


def biom_unavailable():
    "Test if the biom command is available"
    try:
        p = Popen(["biom", "--version"], stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()

        if p.returncode == 0:
            return False
    except Exception:
        # Any failure and we return 'unavailable'
        pass

    return True


class TestBiomValidity(MotusTestCase):
    @unittest.skipIf(biom_unavailable(), "biom command unavailable")
    def test_biom_validity(self):
        "biom validate-table on sample_data/*.biom"

        for biom in glob("sample_data/*.biom"):
            cmd = ("biom", "validate-table", "-i", biom)
            stdout, stderr, exit = self.run_cmd(cmd)

            self.assertRegexpMatches(stdout, "The input file is a valid BIOM-formatted file")


if __name__ == "__main__":
    from simpletap import TAPTestRunner
    unittest.main(testRunner=TAPTestRunner())

# vim: ai sts=4 et sw=4
