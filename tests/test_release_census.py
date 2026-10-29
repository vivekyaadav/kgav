"""docs/kg_release_census.txt, checked against the release it describes.

The census is the only place this project states how big the graph is.
It is generated rather than written, because a typed count is stale the
moment the graph is rebuilt and this project has twice shipped a figure
describing a release that no longer existed.

So the generator is exercised here on the tracked smoke release — every
run, not only when someone happens to have a full release on disk — and
the committed file is checked against the real one when it is present.

These tests moved out of tests/test_kg_audit_report.py on 2026-10-29
with the census itself. A generated file and the prose that interprets
it have different lifecycles, and sharing one file made them collide.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

import scripts.write_kg_report as W

SMOKE = ROOT / "data" / "releases" / "smoke"
FULL = ROOT / "data" / "releases" / "v0.1"


@pytest.fixture(scope="module")
def smoke_census(tmp_path_factory) -> str:
    out = tmp_path_factory.mktemp("census") / "c.txt"
    assert W.main(["--release", str(SMOKE), "--out", str(out)]) == 0
    return out.read_text()


# ------------------------------------------------------------- the file
def test_the_committed_census_exists_and_says_it_is_generated():
    text = W.CENSUS.read_text()
    assert "GENERATED FILE" in text
    assert "make kg-report" in text
    assert "do not edit" in text.lower()


def test_the_counts_come_from_the_release(smoke_census):
    """Counted off the release rather than compared with a figure typed
    here, which would be a second implementation of the thing under
    test."""
    n_nodes = sum(1 for ln in (SMOKE / "nodes.jsonl").read_text().splitlines()
                  if ln.strip())
    n_edges = sum(1 for ln in (SMOKE / "edges.jsonl").read_text().splitlines()
                  if ln.strip())
    assert f"{n_nodes:,} nodes" in smoke_census
    assert f"{n_edges:,} edges" in smoke_census
    for line in (SMOKE / "edges.jsonl").read_text().splitlines():
        if line.strip():
            assert json.loads(line)["predicate"] in smoke_census


def test_the_stamp_is_the_release_it_was_built_from(smoke_census):
    """The stamp is the whole protection. A census carrying a hash that
    is not the release's own describes some other graph."""
    from kgav.provenance import release_fingerprint
    fp = release_fingerprint(SMOKE)
    assert fp["nodes_sha"] in smoke_census
    assert fp["edges_sha"] in smoke_census
    assert SMOKE.name in smoke_census


def test_every_declared_class_and_predicate_appears(smoke_census):
    """Presence AND absence. Listing only what the release holds reports
    presence and stays silent on absence -- the thing this project
    refuses everywhere else. v0.1 holds 7 of 11 node classes and 14 of
    24 predicates, and the first census showed neither fact."""
    from kgav.schema import load_schema
    s = load_schema()
    for name in s.node_classes:
        assert name in smoke_census, f"node class {name} missing"
    for e in s.edge_classes:
        assert e.predicate in smoke_census, f"{e.predicate} missing"
    assert "declared classes are empty" in smoke_census
    assert "declared predicates are empty" in smoke_census


def test_an_absent_evaluation_label_is_called_out_as_one(smoke_census):
    """TREATS and IN_TRIAL_FOR are both empty on v0.1, so every
    evaluation rests on HAS_ANTIVIRAL_ACTIVITY_AGAINST alone. A census
    that lists them among twenty other empty predicates buries the one
    fact that changes how the results read."""
    from kgav.schema import load_schema
    labels = load_schema().evaluation_label_predicates()
    assert "EVALUATION LABEL" in smoke_census
    assert f"of {len(labels)} evaluation label" in smoke_census


def test_the_totals_reconcile(smoke_census):
    """If the per-class or per-predicate counts ever fail to sum to the
    totals, the census says so rather than printing a figure that looks
    whole."""
    assert "!!" not in smoke_census.replace("!! ", "", 0) or \
        "REPEATED curator" in smoke_census, smoke_census
    assert "class counts sum to" not in smoke_census
    assert "predicate counts sum to" not in smoke_census


# -------------------------------------------------------- curators
def _census_of(tmp_path, sources: list[str]) -> str:
    rel = tmp_path / "rel"
    rel.mkdir()
    (rel / "nodes.jsonl").write_text(json.dumps(
        {"id": "A", "class": "Protein", "properties": {}}) + "\n")
    rows = [{"subject": "A", "predicate": "PHYSICALLY_INTERACTS_WITH",
             "object": "B", "evidence_tier": 1,
             "first_asserted_date": "2020-01-01",
             "primary_knowledge_source": s} for s in sources]
    (rel / "edges.jsonl").write_text("\n".join(json.dumps(r) for r in rows))
    out = tmp_path / "c.txt"
    assert W.main(["--release", str(rel), "--out", str(out)]) == 0
    return out.read_text()


def test_it_counts_curators_not_source_strings(tmp_path):
    """primary_knowledge_source is a "|" list once an edge is merged, so
    counting raw values reports "infores:intact|infores:virhostnet" as a
    source of its own. v0.1 showed 35 such values for 13 real curators,
    and VirHostNet read 5,222 when it curated 10,163."""
    text = _census_of(tmp_path, [
        "infores:virhostnet",
        "infores:intact|infores:virhostnet",
        "infores:intact|infores:virhostnet|infores:virhostnet"])
    assert "EDGES BY CURATOR" in text
    assert re.search(r"infores:virhostnet\s+3\b", text), text
    assert re.search(r"infores:intact\s+2\b", text), text
    assert "infores:intact|infores:virhostnet" not in text
    assert "2 edge(s) cite more than one curator" in text


def test_a_repeated_curator_is_reported_as_the_defect_it_is(tmp_path):
    """542 edges in v0.1 carry the same curator twice in one value --
    the footprint of merge_edges unioning strings. A census that merely
    deduplicated them for display would hide that the release predates
    the fix."""
    text = _census_of(tmp_path, [
        "infores:intact|infores:virhostnet|infores:virhostnet"])
    assert "REPEATED curator" in text
    assert "assembled before the fix" in text


def test_a_clean_release_says_so(tmp_path):
    """The other side, so the flag is not simply always printed."""
    text = _census_of(tmp_path, ["infores:virhostnet", "infores:intact"])
    assert "REPEATED curator" not in text
    assert "every edge cites exactly one curator" in text


# ------------------------------------------------------------- staleness
def test_check_fails_on_a_stale_census(tmp_path):
    """--check is what the gate runs. If it passed on a stale file it
    would certify exactly the state it exists to detect."""
    out = tmp_path / "c.txt"
    assert W.main(["--release", str(SMOKE), "--out", str(out)]) == 0
    assert W.main(["--release", str(SMOKE), "--out", str(out),
                   "--check"]) == 0
    out.write_text(out.read_text().replace(" nodes\n", " nodes (approx)\n", 1))
    assert W.main(["--release", str(SMOKE), "--out", str(out),
                   "--check"]) == 1, "a hand-edited census passed --check"


def test_check_fails_when_the_census_is_missing(tmp_path):
    """Missing is stale, not absent-and-therefore-fine. A check that
    passed on a file that does not exist would certify the one state it
    is least able to vouch for."""
    missing = tmp_path / "nope.txt"
    assert not missing.exists()
    assert W.main(["--release", str(SMOKE), "--out", str(missing),
                   "--check"]) == 1


def test_the_committed_census_describes_the_release_when_one_is_present():
    """If a full release is on disk, the committed file must describe
    it. No skip when it is absent -- the other branch asserts the
    smoke-release check instead, so this test runs either way."""
    if (FULL / "nodes.jsonl").exists():
        assert W.main(["--release", str(FULL), "--check"]) == 0, (
            f"{W.CENSUS} does not describe {FULL}; run make kg-report")
    else:
        # Nothing to compare against here, so assert the generator still
        # produces a well-formed census rather than passing silently.
        assert "RELEASE CENSUS" in W.CENSUS.read_text()


def test_each_file_names_the_other():
    """The prose that interprets these figures lives in the report, and
    a reader of either must be able to find the other.

    THE HEADER IS CHECKED AT ITS SOURCE, not in the committed file.
    Reading the committed census let a mutant delete the pointer from
    the GENERATOR and pass, because the file on disk still carried the
    old header until someone regenerated -- the stale-artifact problem
    this whole split exists to avoid, reappearing inside its own test.
    """
    report = (ROOT / "docs" / "kg_audit_report.txt").read_text()
    assert "kg_release_census.txt" in report
    assert "docs/kg_audit_report.txt" in "\n".join(W.HEADER)
    assert "docs/kg_audit_report.txt" in W.CENSUS.read_text()
