"""Gate for results provenance.

The failure this closes: metapath_comparison.json was generated the day
before the selectivity layer added 1,568 activity edges, and was still read
as a baseline afterwards. A positive count that had fallen 1,011 -> 889 read
as a rise from 736, because nothing in the file said which graph it described.
"""
import json
import re
from pathlib import Path

from kgav.provenance import (
    PROVENANCE_KEY,
    check_dir,
    check_file,
    release_fingerprint,
    write_results,
)

ROOT = Path(__file__).resolve().parents[1]


def _release(tmp_path, edges="a\n", nodes="n\n"):
    d = tmp_path / "rel"
    d.mkdir(exist_ok=True)
    (d / "nodes.jsonl").write_text(nodes)
    (d / "edges.jsonl").write_text(edges)
    return d


def test_a_results_file_records_the_release_it_came_from(tmp_path):
    rel = _release(tmp_path)
    out = tmp_path / "r" / "x.json"
    write_results(out, {"auc": 0.826}, rel, "0.9.0")
    doc = json.loads(out.read_text())
    assert doc["auc"] == 0.826
    prov = doc[PROVENANCE_KEY]
    assert prov["release"] == release_fingerprint(rel)
    assert prov["schema_version"] == "0.9.0"
    assert "generated_utc" in prov and "code" in prov


def test_a_matching_release_passes(tmp_path):
    rel = _release(tmp_path)
    out = tmp_path / "r" / "x.json"
    write_results(out, {"auc": 0.826}, rel)
    assert check_file(out, rel)[0] == "ok"


def test_a_release_rebuilt_in_place_is_caught(tmp_path):
    """THE ACTUAL BUG. The directory keeps its name; the graph is different.
    Fingerprinting content rather than the name is what catches it."""
    rel = _release(tmp_path)
    out = tmp_path / "r" / "x.json"
    write_results(out, {"auc": 0.826}, rel)
    (rel / "edges.jsonl").write_text("a\nb\n")      # selectivity layer lands
    status, msg = check_file(out, rel)
    assert status == "stale"
    assert "edges" in msg


def test_a_file_with_no_provenance_is_reported_not_assumed_fine(tmp_path):
    rel = _release(tmp_path)
    p = tmp_path / "old.json"
    p.write_text(json.dumps({"auc": 0.826}))
    status, msg = check_file(p, rel)
    assert status == "unstamped"
    assert "cannot tell which release" in msg


def test_check_dir_reports_every_file(tmp_path):
    rel = _release(tmp_path)
    d = tmp_path / "res"
    d.mkdir()
    write_results(d / "good.json", {}, rel)
    (d / "bare.json").write_text("{}")
    rows = {name: status for name, status, _ in check_dir(d, rel)}
    assert rows == {"good.json": "ok", "bare.json": "unstamped"}


def test_nodes_and_edges_are_both_fingerprinted(tmp_path):
    rel = _release(tmp_path)
    out = tmp_path / "r" / "x.json"
    write_results(out, {}, rel)
    (rel / "nodes.jsonl").write_text("n\nn2\n")     # edges untouched
    assert check_file(out, rel)[0] == "stale"


def test_an_unreadable_file_is_stale_not_silently_skipped(tmp_path):
    rel = _release(tmp_path)
    p = tmp_path / "broken.json"
    p.write_text("{not json")
    assert check_file(p, rel)[0] == "stale"


def test_a_missing_release_is_not_reported_as_current(tmp_path):
    """Two Nones compared equal, so a typo'd release path reported every
    results file as current -- a false "current" from the tool that exists to
    catch false "currents"."""
    from kgav.provenance import check_file, write_results

    release = tmp_path / "rel"
    release.mkdir()
    (release / "nodes.jsonl").write_text('{"id":"A","class":"Gene"}')
    (release / "edges.jsonl").write_text("")
    out = tmp_path / "r.json"
    write_results(out, {"n": 1}, release, "0.12.0")
    assert check_file(out, release)[0] == "ok"

    status, msg = check_file(out, tmp_path / "does-not-exist")
    assert status == "stale"
    assert "cannot verify" in msg


# ------------------- the audit must check the directory in actual use
def test_the_results_audit_checks_the_directory_the_scripts_write_to():
    """make check-results pointed at data/results-corrected, a snapshot from
    2026-09-17 that nothing writes to, while every script defaults to
    data/results. It reported 11 stale files in an unused directory and
    nothing about the one in use -- which held eight files with no provenance
    block at all, the exact case provenance.py was written for.

    Compares the Makefile target against the scripts' own defaults, because
    two places that must agree and are edited separately will drift.
    """
    mk = (ROOT / "Makefile").read_text()
    m = re.search(r"check-results:\s*\n\s*python scripts/check_results\.py "
                  r"\$\(or \$\(RESULTS\),([^)]+)\)", mk)
    assert m, "could not find the check-results target"
    audited = m.group(1).strip()

    written = set()
    for script in ("hard_negatives", "run_baselines", "temporal_split",
                   "compare_metapaths"):
        src = (ROOT / "scripts" / f"{script}.py").read_text()
        for d in re.findall(r'default=root / "(data/results[^"]*)"', src):
            if not d.endswith(".json"):      # a file, not an output directory
                written.add(d)

    assert written, "no script declares a results output directory"
    assert written == {audited}, (
        f"make check-results audits {audited!r} but the scripts write to "
        f"{sorted(written)}. An integrity check aimed at the wrong directory "
        f"reports a verdict about neither.")
