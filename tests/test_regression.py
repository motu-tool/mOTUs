"""Regression tests for the mOTUs profiling pipeline.

Each test runs a motus command as a subprocess against the toy database and
compares its output file(s) to pre-computed expected files stored in
tests/data/output/. The toy database and all input reads live under
tests/data/; expected output files were generated with run.sh.

To regenerate expected output (e.g. after a tool version bump):
  cd tests/data && bash run.sh

Test categories
---------------
test_profile            – parametrised over 14 `motus profile` invocations
                          covering paired/single-end, all counting modes (-y),
                          different min-MGC thresholds (-g), different
                          min-alignment lengths (-l), and reads with /1,/2
                          suffixes or randomly removed pairs.
test_merge              – merges two compatible default profiles; checks the
                          merged output matches the stored expected file.
test_merge_incompatible_fails
                        – confirms that merging profiles created with different
                          parameters (min_mgcs or count_mode) exits non-zero,
                          preventing silent production of corrupted merged files.
"""

import pathlib
import subprocess

import pytest

TESTS_DIR = pathlib.Path(__file__).parent
DATA_DIR = TESTS_DIR / 'data'
INPUT_DIR = DATA_DIR / 'input'
EXPECTED_DIR = DATA_DIR / 'output'

# ---------------------------------------------------------------------------
# Input file shorthands
# ---------------------------------------------------------------------------
# ERR4507416 – paired-end sample, clean pairs
_S1   = str(INPUT_DIR / 'ERR4507416' / 'ERR4507416_1.motus.fastq.gz')
_S2   = str(INPUT_DIR / 'ERR4507416' / 'ERR4507416_2.motus.fastq.gz')
# ERR4507418 – a second paired-end sample, clean pairs
_T1   = str(INPUT_DIR / 'ERR4507418' / 'ERR4507418_1.motus.fastq.gz')
_T2   = str(INPUT_DIR / 'ERR4507418' / 'ERR4507418_2.motus.fastq.gz')
# ERR4507418-suffix – same reads but read names carry /1 and /2 suffixes
_SUF1 = str(INPUT_DIR / 'ERR4507418-suffix' / 'ERR4507418_1.suffix.motus.fastq.gz')
_SUF2 = str(INPUT_DIR / 'ERR4507418-suffix' / 'ERR4507418_2.suffix.motus.fastq.gz')
# ERR4507418-randomremovedreads – ERR4507418 with some reads randomly deleted,
#   intentionally breaking the paired-end pairing for a subset of reads
_RRR1 = str(INPUT_DIR / 'ERR4507418-randomremovedreads' / 'ERR4507418_1.rrr.motus.fastq.gz')
_RRR2 = str(INPUT_DIR / 'ERR4507418-randomremovedreads' / 'ERR4507418_2.rrr.motus.fastq.gz')

# ---------------------------------------------------------------------------
# Profile test cases
# ---------------------------------------------------------------------------
# Each entry is (case_id, motus_args) where motus_args is the variable portion
# of the command. The common flags (-t 1, -db, -o) are appended in the test.
#
# case_id must match a set of files in EXPECTED_DIR:
#   <case_id>          – counts profile (value_type=counts)
#   <case_id>.relab    – relative-abundance profile (written automatically by
#                        `motus profile` alongside every counts file)
#   <case_id>.mgc      – intermediate marker-gene-cluster counts
#   <case_id>.inserts.gz – per-read insert assignments (checked non-empty only)

PROFILE_CASES = [
    # default paired-end run (INSERT_SCALED counts, min_mgcs=3, min_aln=75)
    ('ERR4507416-default',        ['profile', '-f', _S1,   '-r', _S2,   '-n', 'ERR4507416']),
    # second sample – used as a second input for test_merge
    ('ERR4507418-default',        ['profile', '-f', _T1,   '-r', _T2,   '-n', 'ERR4507418']),
    # reads with /1 and /2 suffixes in their names (--skip-pair-check not needed;
    # mOTUs strips the suffixes before pairing)
    ('ERR4507418-withsuffix',     ['profile', '-f', _SUF1, '-r', _SUF2, '-n', 'ERR4507418']),
    # randomly removed reads processed with --skip-pair-check: orphaned reads are
    # silently assigned as singletons, so the profile differs from ERR4507416-broken
    ('ERR4507416-unbroken',       ['profile', '-f', _RRR1, '-r', _RRR2, '-n', 'ERR4507416', '--skip-pair-check']),
    # -g controls the minimum number of detected MGCs before an mOTU is reported;
    # -g 1 is the most permissive (report if any single MGC detected)
    ('ERR4507416-g1',             ['profile', '-f', _S1,   '-r', _S2,   '-n', 'ERR4507416', '-g', '1']),
    # -g 10 is the most stringent (all 10 COG marker genes must be detected)
    ('ERR4507416-g10',            ['profile', '-f', _S1,   '-r', _S2,   '-n', 'ERR4507416', '-g', '10']),
    # -l sets the minimum alignment length in bp; longer cutoff = fewer alignments
    ('ERR4507416-l150',           ['profile', '-f', _S1,   '-r', _S2,   '-n', 'ERR4507416', '-l', '150']),
    # -l 40 is a relaxed cutoff; more short alignments are accepted
    ('ERR4507416-l40',            ['profile', '-f', _S1,   '-r', _S2,   '-n', 'ERR4507416', '-l', '40']),
    # counting modes: -y selects how raw alignments are aggregated
    # INSERT_RAW  – raw insert count per MGC (integer)
    ('ERR4507416-yINSERT_RAW',    ['profile', '-f', _S1,   '-r', _S2,   '-n', 'ERR4507416', '-y', 'INSERT_RAW']),
    # INSERT_NORM – length-normalised insert count (float, sums to 1 across MGCs)
    ('ERR4507416-yINSERT_NORM',   ['profile', '-f', _S1,   '-r', _S2,   '-n', 'ERR4507416', '-y', 'INSERT_NORM']),
    # BASE_RAW    – raw base count per MGC
    ('ERR4507416-yBASE_RAW',      ['profile', '-f', _S1,   '-r', _S2,   '-n', 'ERR4507416', '-y', 'BASE_RAW']),
    # BASE_NORM   – length-normalised base count
    ('ERR4507416-yBASE_NORM',     ['profile', '-f', _S1,   '-r', _S2,   '-n', 'ERR4507416', '-y', 'BASE_NORM']),
    # single-end mode: only a forward file is provided via -s
    ('ERR4507416-default-single', ['profile', '-s', _S1,                '-n', 'ERR4507416']),
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _cmp(actual: pathlib.Path, expected: pathlib.Path) -> None:
    """Assert that two text files have identical content.

    Provides a descriptive failure message showing both paths so that CI logs
    make it easy to identify which output file diverged.
    """
    assert actual.read_text() == expected.read_text(), (
        f'Output mismatch for {actual.name}\n'
        f'  actual:   {actual}\n'
        f'  expected: {expected}'
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('case_id,motus_args', PROFILE_CASES, ids=[c[0] for c in PROFILE_CASES])
def test_profile(case_id, motus_args, toy_db, tmp_path):
    """Run `motus profile` and compare all output files to expected.

    Four files are checked per case:
      <case_id>           – the counts profile written to the path given by -o
      <case_id>.relab     – the relative-abundance profile written alongside it
                            automatically by every `motus profile` invocation
      <case_id>.mgc       – intermediate MGC counts; validates the alignment and
                            counting stage independently of the mOTU aggregation
      <case_id>.inserts.gz – per-read insert-to-MGC assignment log; only checked
                             for non-empty because the content order can vary with
                             multimapper fractional assignments

    All runs use -t 1 (single BWA thread) to guarantee deterministic alignment
    output across platforms.
    """
    out = tmp_path / case_id
    cmd = ['motus'] + motus_args + ['-t', '1', '-db', str(toy_db), '-o', str(out)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, (
        f'motus exited {result.returncode}\nstderr:\n{result.stderr}'
    )

    exp = EXPECTED_DIR
    _cmp(out,                            exp / case_id)
    _cmp(tmp_path / f'{case_id}.relab', exp / f'{case_id}.relab')
    _cmp(tmp_path / f'{case_id}.mgc',   exp / f'{case_id}.mgc')
    assert (tmp_path / f'{case_id}.inserts.gz').stat().st_size > 0


def test_merge(toy_db, tmp_path):
    """Merge two compatible default profiles and verify the merged output.

    Uses the stored expected output files as merge inputs so that this test
    does not depend on test_profile having run first. The merged file encodes
    both sample columns and the shared header; any change to the merge logic or
    output format will be caught here.
    """
    out = tmp_path / 'merge-default'
    cmd = [
        'motus', 'merge',
        '-i',
        str(EXPECTED_DIR / 'ERR4507418-default'),
        str(EXPECTED_DIR / 'ERR4507416-default'),
        '-o', str(out),
        '-db', str(toy_db),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, (
        f'motus merge failed\nstderr:\n{result.stderr}'
    )
    _cmp(out, EXPECTED_DIR / 'merge-default')


def test_profile_broken_pairs_fails(toy_db, tmp_path):
    """Profiling randomly-removed-reads without --skip-pair-check must fail.

    When paired-end reads are mismatched (some reads missing their pair) and
    --skip-pair-check is not supplied, mOTUs detects the inconsistency during
    BAM processing and exits non-zero. The companion test ERR4507416-unbroken
    covers the same input with --skip-pair-check, where the run succeeds.
    """
    out = tmp_path / 'ERR4507416-broken'
    cmd = [
        'motus', 'profile',
        '-f', _RRR1, '-r', _RRR2,
        '-n', 'ERR4507416', '-t', '1',
        '-db', str(toy_db), '-o', str(out),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode != 0, (
        'Expected non-zero exit when profiling mismatched pairs without --skip-pair-check'
    )


@pytest.mark.parametrize('merge_input_b,reason', [
    ('ERR4507416-g10',          'different min_mgcs (-g 10 vs -g 3)'),
    ('ERR4507416-yINSERT_NORM', 'different count mode (INSERT_NORM vs INSERT_SCALED)'),
], ids=['different-min-mgcs', 'different-count-mode'])
def test_merge_incompatible_fails(merge_input_b, reason, toy_db, tmp_path):
    """Merging profiles created with incompatible parameters must fail.

    mOTUs embeds the key parameters (min_mgcs, count_mode, tool version,
    database version) in each output file's header. The merge command validates
    these headers and must exit non-zero when they do not match, preventing
    silent production of scientifically meaningless combined tables.

    Two incompatibility scenarios are tested:
      different-min-mgcs   – ERR4507418-default (g=3) vs ERR4507416-g10 (g=10)
      different-count-mode – ERR4507418-default (INSERT_SCALED) vs
                             ERR4507416-yINSERT_NORM (INSERT_NORM)
    """
    out = tmp_path / 'merge-out'
    cmd = [
        'motus', 'merge',
        '-i',
        str(EXPECTED_DIR / 'ERR4507418-default'),
        str(EXPECTED_DIR / merge_input_b),
        '-o', str(out),
        '-db', str(toy_db),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode != 0, (
        f'Expected non-zero exit for incompatible merge ({reason}), but got 0'
    )
