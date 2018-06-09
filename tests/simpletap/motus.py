# -*- coding: utf-8 -*-

import os
import sys
from subprocess import Popen, PIPE
from unittest import TestCase

PY3 = sys.version_info.major > 2

ROOT = os.path.dirname(os.path.abspath(__name__))


def as_lines(data1, data2):
    "Breaks a multiline string into a list of lines"
    return data1.split('\n'), data1.split('\n')


def as_sam(file1, file2):
    "Convert BAM to SAM to make them comparable"
    def samview(data):
        cmd = ("samtools", "view", data)
        p = Popen(cmd, stdout=PIPE, stderr=PIPE, cwd=ROOT)
        stdout, stderr = p.communicate()

        return stdout.decode()

    return samview(file1), samview(file2)


def prefilter_biom_format(biom1, biom2):
    "Remove timestamp information from BIOM files in order to make them comparable"
    def remove_date(data):
        keep = []
        for line in data:
            if line.lstrip(" ").startswith('"date": "'):
                continue
            else:
                keep.append(line)

        return keep

    return remove_date(biom1), remove_date(biom2)


def ignore_header(data1, data2):
    "Remove header lines from mOTU profiles to make them comparable"
    def ignore(data):
        keep = []
        for line in data:
            if line.startswith("# "):
                continue
            else:
                keep.append(line)

        return keep

    return ignore(data1), ignore(data2)


class MotusTestCase(TestCase):
    def run_cmd(self, cmd, allow_error=False, cmd_prefix=None):
        "Run specified command and test for successful exit code"

        if cmd_prefix is None:
            cmd = tuple(cmd)
        else:
            cmd = cmd_prefix + tuple(cmd)

        p = Popen(cmd, stdout=PIPE, stderr=PIPE, cwd=ROOT)
        stdout, stderr = p.communicate()

        # On Py3 p.communicate() returns bytes and unittest assertions need unicode
        stdout = stdout.decode()
        stderr = stderr.decode()

        if not allow_error:
            self.assertEquals(
                p.returncode, 0,
                msg="Failed to run '{0}' with exit '{1}'.\nSTDOUT: {2}\nSTDERR: {3}".format(
                    cmd, p.returncode, stdout, stderr))

        return stdout, stderr, p.returncode

    def run_motus(self, cmd, allow_error=False):
        "Call motus with specified arguments"

        # tests should run from the tests/ directory
        cmd_prefix = ("../motus",)

        return self.run_cmd(cmd, allow_error=allow_error, cmd_prefix=cmd_prefix)

    def validate_bam(self, profile_file, expected_file):
        "Validate that the produced bam matches expectation"

        profile, expected = as_sam(profile_file, expected_file)
        self.assertEquals(profile, expected)

    def validate_result(self, profile, expected_file):
        "Validate that the produced profile matches expectation"

        self.assertTrue(os.path.isfile(expected_file),
                        msg="File {0} doesn't exist".format(expected_file))

        with open(expected_file) as fh:
            expected = fh.read()

        expected, profile = as_lines(expected, profile)

        # Remove the header which includes the command used to generate the profile
        expected, profile = ignore_header(expected, profile)

        if expected_file.endswith(".biom"):
            # Remove 'date' timestamps from BIOM format
            expected, profile = prefilter_biom_format(expected, profile)

            # We treat BIOM files as a single blob of text
            self.assertEqual(expected, profile)
        else:
            for line_e, line_p in zip(expected, profile):
                # Comments are to be compared as-is
                if line_e.startswith('#') or line_p.startswith('#'):
                    self.assertEquals(line_e, line_p)
                    continue

                # lines containing values should be compared column by column
                all_e = line_e.split('\t')
                all_p = line_p.split('\t')

                # first compare names only
                self.assertEquals(all_e[0], all_p[0])

                # Then all values but convert them to floats first
                values_e = map(float, all_e[1:])
                values_p = map(float, all_p[1:])

                for e, p in zip(values_e, values_p):
                    # Almost equal means equal at 7 digits of precision
                    self.assertAlmostEqual(e, p)


# vim: ai sts=4 et sw=4
