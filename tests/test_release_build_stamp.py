"""The release records what built it, and check_results says when that
is not the code in the tree.

WHY THIS EXISTS. Every results file has carried provenance.code.commit
since provenance.py was written. The release those results describe
carried nothing, so a baseline could say which code computed it while
the graph underneath could not say which code assembled it.

v0.1 was assembled 2026-10-28 from commit 502d81e. Rebuilt the next
day from the same raw inputs it produced a different edges_sha -- the
shipped graph was not reproducible from the tree, every evaluation
rested on it, and no part of the pipeline could report that.
check_results compares results against the release ON DISK and never
against a rebuild, and the release named no commit to compare against.

The stamp does not make a release reproducible. It makes the question
askable.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

import scripts.check_results as CR  # noqa: E402


def _release(tmp_path: Path, code: dict | None, *, manifest=True) -> Path:
    rel = tmp_path / "rel"
    rel.mkdir()
    (rel / "nodes.jsonl").write_text("")
    (rel / "edges.jsonl").write_text("")
    if manifest:
        man: dict = {"release": "v0.1", "schema_version": "0.13.0"}
        if code is not None:
            man["code"] = code
        (rel / "MANIFEST.json").write_text(json.dumps(man))
    return rel


def test_the_assembler_stamps_the_build_commit():
    """Checked at the SOURCE, not on a manifest on disk. A release
    assembled before this change still has no stamp, and reading one
    of those would let the assembler lose the line and pass."""
    src = (ROOT / "scripts" / "assemble_release.py").read_text()
    assert 'manifest["code"] = code_commit(' in src, (
        "assemble_release.py no longer stamps the build commit")
    assert "from kgav.provenance import code_commit" in src


def test_a_matching_commit_says_nothing(tmp_path, capsys):
    """Silence is correct when the release and the tree agree. A line
    printed every run is a line nobody reads."""
    now = CR.code_commit().get("commit")
    CR._report_release_drift(_release(tmp_path, {"commit": now,
                                                 "dirty": False}))
    assert capsys.readouterr().out == ""


def test_a_different_commit_is_reported(tmp_path, capsys):
    CR._report_release_drift(
        _release(tmp_path, {"commit": "0" * 40, "dirty": False}))
    out = capsys.readouterr().out
    assert "000000000000" in out
    assert "is not a failure" in out
    assert "only a rebuild shows that" in out


def test_a_dirty_build_is_called_out(tmp_path, capsys):
    """A commit beside a dirty tree does not identify what built the
    release, and the drift line would otherwise imply it does."""
    CR._report_release_drift(
        _release(tmp_path, {"commit": "0" * 40, "dirty": True}))
    out = capsys.readouterr().out
    assert "DIRTY" in out
    assert "does not\nfully identify what built it" in out


def test_a_clean_build_is_not_called_dirty(tmp_path, capsys):
    """The other direction, so the warning is not simply always
    printed alongside the drift line."""
    CR._report_release_drift(
        _release(tmp_path, {"commit": "0" * 40, "dirty": False}))
    assert "DIRTY" not in capsys.readouterr().out


def test_an_unstamped_release_says_so_rather_than_passing(tmp_path, capsys):
    """Every release built before 2026-10-29 is in this state. Silence
    here would read as agreement, which is the exact failure the stamp
    exists to end."""
    CR._report_release_drift(_release(tmp_path, None))
    out = capsys.readouterr().out
    assert "records no build commit" in out
    assert "Reassemble" in out


def test_a_missing_manifest_is_silent_not_an_error(tmp_path, capsys):
    """The smoke release has no MANIFEST.json and is used throughout
    the suite. Raising here would fail runs that have nothing to do
    with this check."""
    CR._report_release_drift(_release(tmp_path, None, manifest=False))
    assert capsys.readouterr().out == ""


def test_an_unreadable_manifest_does_not_raise(tmp_path, capsys):
    rel = tmp_path / "rel"
    rel.mkdir()
    (rel / "MANIFEST.json").write_text("{not json")
    CR._report_release_drift(rel)
    assert capsys.readouterr().out == ""


def test_the_check_runs_as_part_of_check_results():
    """Pinned at the call site. The function could be correct and
    never invoked, which is how a guard becomes decorative."""
    src = (ROOT / "scripts" / "check_results.py").read_text()
    assert "_report_release_drift(args.release)" in src, (
        "the drift check is defined and not called from main()")
