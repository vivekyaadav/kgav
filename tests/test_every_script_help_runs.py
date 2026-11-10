"""No script may crash while printing its own --help.

WHY. `scripts/hard_negatives.py --help` raised

    TypeError: %o format: an integer is required, not dict

because a help string contained the literal "94.9% of compounds" and
argparse %-formats help text, so `% o` was read as a format spec. The flag
it documents, `--publication-disjoint`, is the one that defines the
headline protocol in docs/evaluation.md's reproduce block, and its help had
never once printed. Found 2026-11-11 by the --help sweep run while
converting fix-list L14, not by a test.

THE RULE IS ABOUT TRACEBACKS, NOT EXIT CODES. Four scripts exit non-zero on
`--help` on purpose and must keep doing so: `compare_results.py` and
`diff_releases.py` print a usage block because they take positional
directories, `diagnose_hop_attrition.py` needs a release, and
`phase1_schema_patch.py` REFUSES to run at all because phase 1 is done and
it would downgrade the schema. Each of those is a deliberate message. An
uncaught exception is not, and that is the difference this test reads --
rather than a list of four names that would go stale the next time a script
is added.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = sorted((ROOT / "scripts").glob("*.py"))


def test_there_are_scripts_to_check():
    """A glob that matches nothing makes every test below vacuous."""
    assert len(SCRIPTS) > 40, f"only {len(SCRIPTS)} scripts found"


@pytest.mark.parametrize("script", SCRIPTS, ids=lambda p: p.name)
def test_help_does_not_raise(script: Path):
    r = subprocess.run([sys.executable, str(script), "--help"],
                       capture_output=True, text=True, timeout=120, cwd=ROOT)
    combined = r.stdout + r.stderr
    assert "Traceback (most recent call last)" not in combined, (
        f"{script.name} --help raised:\n"
        + combined[combined.index("Traceback (most recent call last)"):][:900])
