"""The gate must prove it found the planted errors, not that it failed.

`make gate` runs scripts/validate_release.py against data/releases/smoke,
a release with four deliberately planted violations. Until 2026-11-07 the
target read ANY non-zero exit as success:

    if python scripts/validate_release.py data/releases/smoke; then
        echo "FAIL: validator accepted planted errors"; exit 1
    else
        echo "OK: validator rejected planted errors"
    fi

An ImportError exits non-zero. A missing smoke directory exits non-zero. A
schema that will not load exits non-zero. All three printed OK. The gate
could not tell a working validator from a broken one, which is the only
thing a negative control exists to establish.

`--expect` makes it a positive claim -- these codes, these counts, nothing
else -- and these tests hold the three properties that claim needs: it
passes on the planted set, it fails when the planted set is not found, and
it fails when the process never got far enough to look.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "validate_release.py"
SMOKE = ROOT / "data" / "releases" / "smoke"

# The census data/releases/smoke is built to produce. Kept here as well as in
# the Makefile on purpose: if a schema change stops one of these firing, the
# Makefile tells whoever runs `make gate` and this tells whoever runs the
# suite, and neither can be quietly edited to match a weaker validator.
PLANTED = "CONSTRAINT=1,MISSING_PROVENANCE=1,MISSING_QUALIFIER=2"


def run(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(SCRIPT), *args],
                          capture_output=True, text=True, cwd=ROOT)


@pytest.mark.skipif(not SMOKE.exists(), reason="data/releases/smoke absent")
def test_the_planted_violations_are_found_exactly():
    r = run(str(SMOKE), "--expect", PLANTED)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "NEGATIVE CONTROL PASS" in r.stdout


@pytest.mark.skipif(not SMOKE.exists(), reason="data/releases/smoke absent")
def test_a_wrong_count_fails_and_names_the_code():
    """The failure has to say WHICH code disagreed. A bare non-zero exit is
    what this whole file exists to stop being treated as information."""
    r = run(str(SMOKE), "--expect", "CONSTRAINT=9,MISSING_PROVENANCE=1,"
                                    "MISSING_QUALIFIER=2")
    assert r.returncode == 1
    assert "NEGATIVE CONTROL FAIL" in r.stdout
    assert "CONSTRAINT" in r.stdout
    assert "expected    9" in r.stdout and "found    1" in r.stdout


@pytest.mark.skipif(not SMOKE.exists(), reason="data/releases/smoke absent")
def test_an_extra_code_fails_even_when_the_planted_ones_are_present():
    """A release that reports the planted four PLUS something new is not the
    control passing. Comparing totals would have let that through."""
    r = run(str(SMOKE), "--expect", PLANTED + ",DANGLING_EDGE=3")
    assert r.returncode == 1
    assert "DANGLING_EDGE" in r.stdout


def test_a_release_with_no_violations_fails_the_control():
    """The inverse of the original bug: a clean release must not satisfy a
    negative control. Exercised on a release the fixture builds itself, so
    this one runs with no data on disk."""
    r = run("data/releases/does-not-exist", "--expect", PLANTED)
    assert r.returncode != 0
    assert "NEGATIVE CONTROL PASS" not in r.stdout


def test_the_makefile_and_this_file_agree_on_what_is_planted():
    """Two copies of the expectation, one check that they match. Editing one
    to silence a failure leaves the other to fail instead."""
    makefile = (ROOT / "Makefile").read_text()
    assert PLANTED in makefile, (
        "tests/test_negative_control.py PLANTED and the --expect in the "
        "Makefile's gate target have drifted apart")
