"""Results live in one place, and `make check-results` audits that place.

FIX-LIST P9. `data/results-corrected` and `data/results-prechange` sat
beside `data/results` with 22 tracked files between them, and nothing read
either: every reference in the repository is a comment describing the
defect they caused, not a use. The Makefile's own note records it --
"AUDIT WHAT THE SCRIPTS WRITE. This pointed at data/results-corrected" --
as do test_provenance.py twice and test_results_match_the_release.py.

WHY DELETING THEM IS NOT HOUSEKEEPING. `check_dir` audits the directory it
is given and the audit is wired to `data/results`, so a sibling directory
of JSON files that look exactly like results is invisible to it. Both were
worse than stale:

    data/results-corrected   release 73ebb09d59aa58b3, 2026-09-17 --
                             a release that is not on disk any more
    data/results-prechange   NO provenance block at all -- it predates
                             the stamp, so no release and no commit can
                             be attached to those numbers

kg_audit_report.txt says data/results is tracked "on purpose: a result that
exists only on the machine that produced it is not evidence for anything".
A result whose release cannot be named is the same claim one step further
on, and these were in the repository rather than on one machine.

They remain in git history, so nothing is lost that was ever attributable.
"""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

#: The one results directory, which is the one `make check-results` audits.
AUDITED = "results"


def test_no_results_directory_sits_outside_the_audit():
    """A new snapshot belongs INSIDE data/results, where the audit will
    stamp and check it, or outside the repository entirely. Beside it is
    the one place where it looks like evidence and is checked by nothing.
    """
    strays = sorted(
        p.name for p in DATA.iterdir()
        if p.is_dir() and p.name.startswith("results") and p.name != AUDITED)
    assert not strays, (
        f"{strays} sit beside data/{AUDITED} and `make check-results` does "
        f"not see them. Move a snapshot out of the repository or into "
        f"data/{AUDITED}, where it is stamped against a release")


def test_the_audited_directory_is_still_there():
    """So the check above cannot pass by everything having been deleted."""
    assert (DATA / AUDITED).is_dir()
    assert any((DATA / AUDITED).glob("*.json")), (
        f"data/{AUDITED} holds no results, which would make the test above "
        f"vacuous")
