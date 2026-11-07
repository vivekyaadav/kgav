"""The committed results must describe the release on disk.

THE GAP THIS CLOSES. data/results holds 28 files, all committed to git and
quoted throughout docs/evaluation.md. The release they were computed from
is NOT committed -- data/releases carries only .gitkeep and the smoke
fixture. So the repository ships figures whose graph no reader has, and
until 2026-11-07 nothing in the test suite checked that those figures still
matched the graph the author did have.

The chain is doc -> results -> release -> code, and each link had a guard
except one:

  doc -> results    five test files re-derive every table cell by cell
                    (test_claim_block, test_evaluation_table,
                     test_sequence_table, test_landscape_table,
                     test_label_overlap_figures) -- 92 tests
  results -> release  provenance.check_dir, run ONLY by `make
                      check-results`, which is not in `make test` and not
                      in CI
  release -> code     provenance.check_layer, via `make rebuild-plan`

Measured on 2026-11-07, immediately after a rebuild: `make check-results`
reported 0 current and 28 stale, and all 92 doc tests PASSED. They pin the
document to the JSON, and the JSON was still internally consistent with the
document; both described a graph that no longer existed. `pytest -q` was
green on a repository whose every published number was stale.

WHY check_dir IS NOT ALREADY EXERCISED HERE. tests/test_provenance.py tests
it thoroughly -- on tmp_path fixtures it builds itself. A unit test asks the
question its author thought to ask; nobody pointed it at the real pair.

SKIPPED WITHOUT A RELEASE, deliberately. CI has no release, so this cannot
run there and must not pretend to: a test that silently passes when it
cannot look is the thing this file exists to stop. It is a local gate, and
`make test` now runs it like any other test.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.provenance import check_dir  # noqa: E402

RESULTS = ROOT / "data" / "results"
RELEASE = ROOT / "data" / "releases" / "v0.1"

needs_both = pytest.mark.skipif(
    not (RELEASE / "edges.jsonl").exists() or not RESULTS.is_dir(),
    reason="needs both data/releases/v0.1 and data/results on disk")


@needs_both
def test_every_committed_result_describes_the_release_on_disk():
    """The assertion `make check-results` makes, inside the suite.

    A stale file describes a graph that no longer exists, and comparing
    against it attributes release drift to whatever change is under test.
    """
    rows = check_dir(RESULTS, RELEASE)
    bad = [(name, status, msg) for name, status, msg in rows
           if status not in ("ok", "refused")]
    assert not bad, (
        f"{len(bad)} of {len(rows)} results files no longer match the "
        f"release:\n" + "\n".join(f"  {n:<44} {s}\n      {m}"
                                  for n, s, m in bad)
        + "\n\nRegenerate with `make reproduce`. Every figure in "
          "docs/evaluation.md is derived from these files, and the five "
          "test files that pin that document compare it against THEM -- "
          "so they stay green while the whole set goes stale.")


@needs_both
def test_the_audit_actually_examined_something():
    """The inverse failure: check_dir globbing an empty or wrong directory
    returns no rows and the test above would pass having checked nothing.
    The Makefile comment for check-results records this exact mistake --
    it was aimed at data/results-corrected, a directory nothing writes."""
    rows = check_dir(RESULTS, RELEASE)
    assert len(rows) >= 20, (
        f"check_dir examined only {len(rows)} file(s) under {RESULTS}; "
        f"an audit aimed at the wrong directory returns a verdict about "
        f"neither")


def test_make_test_runs_this_gate():
    """It is only a gate if it runs. check-results is a separate target
    that is not in `make test` and not in CI, which is why a stale set
    survived in the first place."""
    mk = (ROOT / "Makefile").read_text()
    body = mk.split("\ntest:", 1)[1].split("\n\n", 1)[0]
    assert "pytest" in body, "the test target no longer runs pytest"
