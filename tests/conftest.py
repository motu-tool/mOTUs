import pathlib
import subprocess

import pytest

TESTS_DIR = pathlib.Path(__file__).parent
TOY_DB_DIR = TESTS_DIR / 'data' / 'motus4.1-toy-db'


@pytest.fixture(scope='session')
def toy_db():
    """Download the toy database once per test session; reuse on disk if already present.

    The toy database is a lightweight subset of the full mOTUs database used
    exclusively for testing. It is not checked into git (too large), so this
    fixture downloads it on first use via `motus downloadMGDB --toy`.

    The database is cached at tests/data/motus4.1-toy-db/. Subsequent pytest
    runs in the same environment detect the marker file
    (db_mOTU/db_mOTU.downloaded) and skip the download, keeping the test suite
    fast after the first run.

    Returns the Path to the toy-db *parent* folder (the value expected by the
    `-db` flag of all motus commands).
    """
    marker = TOY_DB_DIR / 'db_mOTU' / 'db_mOTU.downloaded'
    if not marker.exists():
        subprocess.run(
            ['motus', 'downloadMGDB', '--toy', '-db', str(TOY_DB_DIR)],
            check=True,
        )
    return TOY_DB_DIR
