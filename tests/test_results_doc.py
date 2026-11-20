"""docs/results.md states what was measured, and must not overstate it.

WHY A SEPARATE DOCUMENT. docs/evaluation.md is the working record: how
each result was arrived at, including the versions that were wrong. The
strongest finding sits four sections deep behind two negatives and a
decomposition, which is right for the people who built it and wrong for
anyone else. docs/results.md leads with the positive results.

THAT IS EXACTLY WHEN A DOCUMENT OVERSTATES, so these tests are about the
qualifications rather than the figures. The figures are generated and
`--check` guards them; what a generator cannot guard is a claim that
quietly loses its caveat, so each one is pinned to the measurement that
forces it:

  - the claimed cell count is the count of cells that clear ABOVE a
    degree baseline that does NOT clear, not the count that clear;
  - the prospective section says that none of its clearing cells beats
    a baseline, because none has a usable comparator;
  - the hit-rate qualification is in the same document, with the cell
    that clears its alpha AND holds no actives at the top of its list;
  - the "not claimed" list names the knowledge graph.

A positive result that needs its caveat moved away from it to look good
is not a positive result.
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "results.md"
SCRIPT = ROOT / "scripts" / "write_results_doc.py"
RESULTS = ROOT / "data" / "results"

needs_results = pytest.mark.skipif(
    not (RESULTS / "structure_scorer_strength_matched.json").exists(),
    reason="no results on disk")


def _text() -> str:
    return DOC.read_text()


def _flat() -> str:
    return " ".join(_text().split())


# ------------------------------------------------------------- currency
@needs_results
@pytest.mark.slow
def test_the_document_is_current():
    """`make results-doc-check` is this. Slow: the generator re-runs
    hit_rate_at_k so the qualification is the measured one."""
    r = subprocess.run([sys.executable, str(SCRIPT), "--check"],
                       capture_output=True, text=True, cwd=ROOT, timeout=1800)
    assert r.returncode == 0, r.stdout + r.stderr


def test_it_says_it_is_generated():
    assert "GENERATED FILE" in _text().splitlines()[0]


def test_no_release_figure_is_typed_into_the_generator():
    """Same rule as every other reporting script here."""
    manifest = ROOT / "data" / "releases" / "v0.1" / "MANIFEST.json"
    if not manifest.exists():
        pytest.skip("no release")
    census = json.loads(manifest.read_text())["nodes_by_class"]
    figures = {v for v in census.values() if v >= 1000}
    figures.add(sum(census.values()))
    src = SCRIPT.read_text().replace(",", "").replace("_", "")
    literals = {int(m.group(0)) for m in re.finditer(r"\b\d{4,}\b", src)}
    typed = sorted(figures & literals)
    assert not typed, f"{SCRIPT.name} types the release figure(s) {typed}"


# ------------------------------------------- the claim is the right one
@needs_results
def test_the_headline_count_is_cells_that_beat_a_non_clearing_baseline():
    """NOT 'cells that clear'. Twelve cells clear something somewhere;
    the claim rests on the ones that clear above a degree baseline that
    does NOT clear, which is a smaller and harder number.

    Re-derived here from the two registered control files.
    """
    import sys as _sys
    _sys.path.insert(0, str(ROOT / "src"))
    from kgav import virus_register as vr

    names = {"NCBITaxon:" + v["taxon"]: v["name"] for v in vr.load()["viruses"]}
    pops = ("primary_chembl", "secondary_ncats", "secondary_ncats_entry")

    def cells(path, channel):
        doc = json.loads((RESULTS / path).read_text())
        out = {}
        for pop in pops:
            for taxon, cell in (doc.get(pop) or {}).items():
                r = cell.get(channel)
                if isinstance(r, dict) and "auc" in r:
                    out[(pop, names.get(taxon, taxon))] = (cell, r)
        return out

    sym = cells("structure_scorer_symmetric.json", "structure_max_k10")
    stg = cells("structure_scorer_strength_matched.json", "structure_max")

    def clears(r):
        return bool(r.get("clears_corrected_alpha")) and r.get(
            "direction") != "below"

    both = [k for k in sym
            if clears(sym[k][1]) and k in stg and clears(stg[k][1])]
    beat = [k for k in both
            if not clears(stg[k][0].get("degree") or {})
            and (stg[k][0].get("degree") or {}).get("auc", 1) < stg[k][1]["auc"]]

    # EVERY OCCURRENCE, NOT THE FIRST. The count appears twice -- in
    # Result 1's summary line and again in "What is claimed" -- so
    # asserting that it APPEARS passes while one of them says something
    # else. The same hole was found in the claim-block tests on
    # 2026-11-20 with a power figure that appeared in two places.
    flat = _flat()
    stated = re.findall(r"(\d+) of (\d+) registered cells", flat)
    assert stated, "the document states no cell count at all"
    assert all(s == (str(len(beat)), str(len(stg))) for s in stated), (
        f"{len(beat)} of {len(stg)} cells clear above a non-clearing "
        f"degree baseline; the document states {set(stated)}")
    assert len(beat) <= len(both), (
        "sanity: the beaten-baseline set cannot exceed the clearing set")


@needs_results
def test_the_prospective_section_says_no_cell_beats_a_baseline():
    """The strongest-sounding result in the document, and the one most
    able to be read as more than it is. Six cells clear; none of them is
    shown to beat a comparator, because none has a usable one."""
    doc = json.loads((RESULTS / "temporal_structure_scorer.json").read_text())
    usable = len(doc.get("comparator_usable") or [])
    n = doc["n_evaluable"]
    flat = _flat()
    assert f"{usable} of {n} cells have a usable degree comparator" in flat
    assert "not shown to beat a baseline" in flat, (
        "the prospective section does not say its clearing cells beat "
        "nothing")
    assert usable == 0 or "none to beat" not in flat, (
        "a comparator became usable and the wording no longer holds")


@needs_results
def test_cells_that_clear_below_chance_are_reported_as_such():
    """Zika and West Nile rank post-cutoff actives reliably LAST. A
    results document that reported them as merely 'not clearing' would
    be hiding a finding."""
    doc = json.loads((RESULTS / "temporal_structure_scorer.json").read_text())
    below = [c["name"] for c in doc["cells"].values()
             if (c.get("primary") or {}).get("clears_corrected_alpha")
             and c["primary"].get("direction") == "below"]
    assert below, "no cell clears below chance; this test asserts nothing"
    text = _text()
    for name in below:
        assert re.search(rf"\|\s*{re.escape(name)}\s*\|.*clears BELOW chance",
                         text), f"{name} is not marked as clearing below"
    assert "clears BELOW chance" in _flat()


# ------------------------------------- the qualification travels with it
@needs_results
def test_the_hit_rate_qualification_is_in_the_same_document():
    """Not deferred to another file. The point of this document is that
    a reader who takes the positive result also takes its limit."""
    flat = _flat()
    assert "An AUC is not a hit rate" in flat
    assert "precision@10" in flat
    # The force of it: a cell that CLEARS and whose top holds nothing.
    assert "and it CLEARS" in flat
    assert "analogue retrieval, not prediction" in flat


@needs_results
def test_the_qualification_names_the_controlled_precision_not_only_the_raw():
    """Reporting only the uncontrolled 50% would invert the finding."""
    flat = _flat()
    assert "no exclusion" in flat and "symmetric k = 10" in flat
    assert "strength-matched" in flat
    # All three rows present means the collapse is visible.
    rows = re.findall(r"\|\s*(no exclusion|symmetric k = 10|"
                      r"strength-matched[^|]*)\s*\|\s*([\d.]+)%", _text())
    assert len(rows) >= 3, f"only {len(rows)} control rows: {rows}"
    raw = float(rows[0][1])
    controlled = [float(v) for _l, v in rows[1:]]
    assert raw > max(controlled), (
        "the uncontrolled precision is not the highest; the paragraph's "
        "argument has changed")


@needs_results
def test_enrichment_is_named_as_the_column_that_matters():
    """A precision of 90% against a 70% base rate is not a result, and a
    document leading with positives has to say so."""
    flat = _flat()
    assert "Enrichment, not precision, is the column that matters" in flat


# ---------------------------------------------- what it refuses to claim
def test_the_not_claimed_list_names_the_knowledge_graph():
    """This is a knowledge-graph project whose knowledge graph does not
    contribute. A results document that omitted that would be the single
    most misleading thing it could do."""
    flat = _flat()
    assert "**Not claimed.**" in flat
    assert "That the knowledge graph contributes to this" in flat
    assert "It does not" in flat


def test_the_not_claimed_list_covers_every_limit_that_has_a_measurement():
    flat = _flat()
    for needle in ("usable hit rate",
                   "beats a baseline prospectively",
                   "Anything clinical",
                   "RESEMBLES a known active"):
        assert needle in flat, f"the not-claimed list is missing {needle!r}"


def test_it_states_the_held_out_predicates():
    """The one methodological fact everything else rests on."""
    flat = _flat()
    assert "held out of every reasoning path" in flat
    for pred in ("HAS_ANTIVIRAL_ACTIVITY_AGAINST", "MEASURED_INACTIVE_AGAINST",
                 "TREATS", "IN_TRIAL_FOR"):
        assert pred in flat, f"{pred} is not named as held out"


@needs_results
def test_every_registration_is_linked():
    """A claim governed by a registration the reader cannot find is a
    claim they have to take on trust."""
    flat = _flat()
    missing = [p.name for p in sorted((ROOT / "docs").glob("prereg_*.md"))
               if p.name not in flat]
    assert not missing, f"registrations not linked: {missing}"
