"""The release plan's temporal note, re-derived rather than read.

WHY THIS FILE EXISTS. Section D of docs/v1_release_plan.md carried this
sentence from 2026-10-24 to 2026-11-11:

    all 150,432 PHYSICALLY_INTERACTS_WITH edges are dropped from any
    pre-2021 training cut ... M3 and M5 therefore have zero edges in the
    temporal protocol

Four things in it were wrong, and fix-list C1, C2 and V4 had already fixed
or measured every one. The release holds 153,533 such edges, not 150,432.
11,744 of them SURVIVE the cut, not zero. Every survivor is virus-host, so
M3 runs and only M5 cannot -- which is why M5 was missing from every
SARS-CoV-2 temporal table, a fact the note concealed by saying both were
empty. And the two mechanisms it blamed were gone: ingest_vhppi dates each
edge from its own PMID year (C1) and `string_release_date` was renamed
`string_file_mtime` (V4).

A note in a planning document is not load-bearing until someone picks the
work up, which is exactly when a wrong figure does its damage. So the
figures are derived from the release and from the temporal results that
`temporal_split.py` writes for this purpose, and the note has to quote what
they say.
"""
from __future__ import annotations

import json
import re
import sys
from collections import Counter
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

PLAN = ROOT / "docs/v1_release_plan.md"
TEMPORAL = ROOT / "data/results/temporal_2021_include.json"
PREDICATE = "PHYSICALLY_INTERACTS_WITH"


@pytest.fixture(scope="module")
def plan() -> str:
    assert PLAN.exists(), "the release plan is missing"
    return " ".join(PLAN.read_text().split())


def live_claims(plan: str) -> str:
    """The plan with QUOTED material removed.

    WHY THIS IS NEEDED, and the first version of this file did without it
    and failed on its own correction. A document that records what it got
    wrong necessarily contains the wrong figure -- the corrected note
    quotes "all 150,432 PHYSICALLY_INTERACTS_WITH edges are dropped ... M3
    and M5 therefore have zero edges" as the sentence it is replacing. A
    regex over the whole document reads that quotation as an assertion, so
    the absence tests failed on text whose presence is the point.

    scripts/diagnose_doc_figures.py says the same thing about
    docs/evaluation.md at length: "This document records superseded figures
    deliberately -- the value that was wrong, and why -- and those SHOULD
    be absent from the current results."

    So a superseded figure belongs inside quotation marks and a live claim
    does not, and that convention is what these tests enforce. The cost is
    real and is stated: a live claim written inside quotes would be
    invisible here.
    """
    return re.sub(r'"[^"]*"', " ", plan)


@pytest.fixture(scope="module")
def release_ppi() -> Counter:
    """source -> edge count, for every PPI edge in the release."""
    from kgav.pipeline import release_dir
    rd = release_dir(ROOT)
    if not (rd / "edges.jsonl").exists():
        pytest.skip("no assembled release on disk")
    out: Counter = Counter()
    for line in (rd / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        if e["predicate"] == PREDICATE:
            out[e.get("primary_knowledge_source", "?")] += 1
    return out


def test_the_note_quotes_the_releases_own_ppi_total(plan, release_ppi):
    """150,432 was the figure for a release that no longer exists."""
    total = sum(release_ppi.values())
    assert f"{total:,}" in plan, (
        f"the release holds {total:,} {PREDICATE} edges and the plan's "
        f"temporal note does not quote that number")
    assert "150,432" not in live_claims(plan), (
        "the plan asserts the superseded PPI total 150,432 outside a "
        "quotation; recording it as what was wrong is fine, claiming it is "
        "not")


def test_the_note_quotes_the_surviving_count_and_it_is_not_zero(plan):
    """THE HALF THAT MATTERED. "M3 and M5 therefore have zero edges" was
    the claim, and it hid why M5 alone was missing from the tables."""
    if not TEMPORAL.exists():
        pytest.skip("no temporal result on disk")
    d = json.loads(TEMPORAL.read_text())
    surviving = d["training_graph_predicates"].get(PREDICATE, 0)
    assert surviving > 0, (
        "the 2021 training graph holds no PPI edges after all, which would "
        "make the plan's original claim true; re-read this test")
    assert f"{surviving:,}" in plan, (
        f"{surviving:,} {PREDICATE} edges survive the 2021 cut and the "
        f"plan's note does not quote that number")


def test_no_surviving_ppi_edge_is_strings(plan):
    """The reason M3 runs and M5 does not: every survivor is virus-host,
    and STRING -- which is all of the host-host side -- has no per-edge
    date. Derived from the per-source census temporal_split.py writes."""
    if not TEMPORAL.exists():
        pytest.skip("no temporal result on disk")
    d = json.loads(TEMPORAL.read_text())
    per_source = {
        k.split("|", 1)[1]: v
        for k, v in d["training_graph_predicate_sources"].items()
        if k.split("|", 1)[0] == PREDICATE
    }
    assert per_source, "the temporal result records no PPI source census"
    string_side = sum(v for k, v in per_source.items() if "string" in k)
    assert string_side == 0, (
        f"{string_side:,} surviving PPI edges are STRING's, so the plan's "
        f"account of why M5 cannot run is wrong")
    # and the census must add up to the total the same file reports
    total = d["training_graph_predicates"].get(PREDICATE, 0)
    assert sum(per_source.values()) == total, (
        f"the per-source census sums to {sum(per_source.values()):,} and the "
        f"predicate count is {total:,}; one of them is not what it says")


def test_the_dropped_string_count_is_the_releases(plan, release_ppi):
    """141,019 is quoted in three documents; here it is derived."""
    string_edges = sum(v for k, v in release_ppi.items() if "string" in k)
    assert f"{string_edges:,}" in plan, (
        f"{string_edges:,} PPI edges in the release are STRING's and the "
        f"plan does not quote that number")


def test_the_note_no_longer_blames_mechanisms_that_were_fixed(plan):
    """C1 gave every vhppi edge its own PMID year and V4 renamed
    string_release_date. A note still naming those as the cause sends the
    next reader to fix something that is already fixed."""
    assert ("string_release_date" not in live_claims(plan)
            or "renamed" in plan), (
        "the plan names `string_release_date` without saying it was renamed "
        "to string_file_mtime (fix-list V4)")
    m = re.search(r"M3[^.]{0,80}zero edges", live_claims(plan))
    assert not m, (
        "the plan still says M3 has zero edges in the temporal protocol; it "
        "runs, and only M5 cannot")
