"""Three of the six prospective cells are analogy, and three are not.

WHAT THIS CLOSES. doc2/registrations/prereg_temporal_structure_scorer.md left H3 open
in as many words -- "a pointer, not a finding" -- and the pointer was
right and understated. The registered AUC recomputed with the k nearest
pre-cutoff actives withheld splits the six clearing cells in two:

    DENV-2, HCoV-229E, SARS-CoV-2   still clear at k = 40
    DENV-1, DENV-3                  stop clearing before k = 10
    DENV-4                          clears at k = 10, not at k = 20

WHY IT MATTERS MORE THAN THE SPLIT ITSELF. The headline AUC is close to
anti-correlated with robustness: DENV-1 has the highest primary of the
six and loses the most, SARS-CoV-2 has the lowest and barely moves. So
"6 of 10 clear" is true and hides two different kinds of result, and a
reader cannot tell which a cell is from its published number.

These tests hold both documents to the diagnostic, and hold the
diagnostic to the registered run it extends: k=0 and k=10 must
reproduce the published primary and secondary columns, or the ladder is
measuring something else and the rest of it means nothing.
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "diagnose_prospective_analogy.py"
PREREG = ROOT / "doc2" / "registrations" / "prereg_temporal_structure_scorer.md"
RESULTS_DOC = ROOT / "doc2" / "08-results.md"
TEMPORAL = ROOT / "data" / "results" / "temporal_structure_scorer.json"
SPLIT = ROOT / "data" / "releases" / "v0.1-train2021"

needs_split = pytest.mark.skipif(
    not (SPLIT / "TEST_LABELS.json").exists(),
    reason="no 2021 temporal split on disk")


@pytest.fixture(scope="module")
def ladder() -> dict:
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "a.json"
        r = subprocess.run([sys.executable, str(SCRIPT), "--out", str(out)],
                           capture_output=True, text=True, cwd=ROOT,
                           timeout=1800)
        assert r.returncode == 0, r.stdout + r.stderr
        return json.loads(out.read_text())


def _amendment() -> str:
    return " ".join(
        PREREG.read_text().split("## Amendment, 2026-11-21 — H3 answered",
                                 1)[1].split())


# ------------------------------------------- the ladder is the real thing
@needs_split
def test_the_rungs_that_are_published_reproduce_the_published_figures(ladder):
    """THE GUARD ON EVERY OTHER NUMBER HERE. k=0 is the registration's
    primary and k=10 its secondary. If those do not match the results
    file, the ladder is scoring something other than the registered
    protocol and the deeper rungs mean nothing."""
    doc = json.loads(TEMPORAL.read_text())
    for c in ladder["cells"]:
        cell = doc["cells"][c["taxon"]]
        assert c["ladder"]["0"]["auc"] == pytest.approx(
            cell["primary"]["auc"], abs=5e-4), f"{c['name']} k=0"
        assert c["ladder"]["10"]["auc"] == pytest.approx(
            cell["symmetric_k10"]["auc"], abs=5e-4), f"{c['name']} k=10"


@needs_split
def test_it_ladders_exactly_the_registered_clearing_cells(ladder):
    doc = json.loads(TEMPORAL.read_text())
    assert {c["taxon"] for c in ladder["cells"]} == set(doc["clears"])


@needs_split
def test_the_split_is_real_and_not_an_artefact_of_one_rung(ladder):
    """Both kinds must exist, or the finding is 'all six behave alike'
    and the amendment says otherwise."""
    deep = [c for c in ladder["cells"] if (c["clears_up_to_k"] or -1) >= 40]
    shallow = [c for c in ladder["cells"] if (c["clears_up_to_k"] or -1) < 10]
    assert deep, "no cell survives to k=40; the amendment claims three do"
    assert shallow, "no cell collapses early; the amendment claims two do"
    # And the deep ones must be deep by a margin, not by rounding.
    for c in deep:
        drop = c["ladder"]["0"]["auc"] - c["ladder"]["40"]["auc"]
        assert drop < 0.2, (
            f"{c['name']} is called robust and drops {drop:.3f}")
    for c in shallow:
        drop = c["ladder"]["0"]["auc"] - c["ladder"]["40"]["auc"]
        assert drop > 0.2, (
            f"{c['name']} is called analogy-driven and drops only {drop:.3f}")


@needs_split
def test_the_headline_auc_does_not_predict_robustness(ladder):
    """The finding, as a property rather than as a sentence: the cell
    with the highest published AUC is not the most robust."""
    by0 = sorted(ladder["cells"], key=lambda c: -c["ladder"]["0"]["auc"])
    by40 = sorted(ladder["cells"], key=lambda c: -c["ladder"]["40"]["auc"])
    assert by0[0]["name"] != by40[0]["name"], (
        "the highest published AUC is now also the most robust; the "
        "amendment's point no longer holds")
    assert by0[0]["name"] == by40[-1]["name"], (
        f"the amendment says the highest primary ({by0[0]['name']}) is "
        f"the LEAST robust; it is now {by40[-1]['name']}")


# --------------------------------------- both documents, against the run
@needs_split
def test_every_ladder_figure_in_the_amendment_is_the_diagnostics(ladder):
    produced = {f"{c['ladder'][str(k)]['auc']:.3f}"
                for c in ladder["cells"] for k in ladder["ks"]}
    produced |= {f"{c['median_tc_positive']:.3f}" for c in ladder["cells"]}
    produced |= {f"{c['median_tc_negative']:.3f}" for c in ladder["cells"]}
    produced |= {f"{c['median_tc_positive'] - c['median_tc_negative']:+.3f}"
                 for c in ladder["cells"]}
    # Drops between any two rungs are derived FROM the ladder, so they
    # are legitimate figures for the amendment to quote -- the check is
    # that a number is traceable to the diagnostic, not that it is a raw
    # rung. SARS-CoV-2's -0.007 is the whole point of the paragraph it
    # appears in.
    produced |= {f"{c['ladder'][str(b)]['auc'] - c['ladder'][str(a)]['auc']:+.3f}"
                 for c in ladder["cells"]
                 for a in ladder["ks"] for b in ladder["ks"]}
    text = PREREG.read_text().split(
        "## Amendment, 2026-11-21 — H3 answered", 1)[1]
    printed = set(re.findall(r"[-+]?\b0\.\d{3}\b", text))
    stray = sorted(printed - produced)
    assert not stray, (
        f"the amendment prints {stray}, which the diagnostic does not "
        f"produce")
    assert len(printed) >= len(ladder["cells"]) * len(ladder["ks"]), (
        "the amendment prints fewer figures than the ladder has rungs")

    # AND THE COUNTS IN PROSE. "3 of 6 survive" is the sentence a reader
    # takes away, and a mutation of it passed every figure check above --
    # the same shape of hole as asserting a figure merely APPEARS.
    flat = _amendment()
    deep = [c for c in ladder["cells"] if (c["clears_up_to_k"] or -1) >= 40]
    shallow = [c for c in ladder["cells"] if (c["clears_up_to_k"] or -1) < 10]
    assert f"{len(deep)} of {len(ladder['cells'])} survive" in flat, (
        f"{len(deep)} of {len(ladder['cells'])} survive to k=40; the "
        f"amendment says otherwise")
    assert f"{len(shallow)} stop clearing before k = 10" in flat, (
        f"{len(shallow)} cells stop clearing before k=10")
    for c in deep:
        assert c["name"] in flat, f"{c['name']} survives and is not named"
    for c in shallow:
        assert c["name"] in flat, f"{c['name']} collapses and is not named"


@needs_split
def test_the_results_document_carries_the_same_split(ladder):
    """A results document that reported '6 of 10 clear' without this
    would be true and misleading, which is the one thing that document
    is for."""
    flat = " ".join(RESULTS_DOC.read_text().split())
    assert "Three of those six are analogy, and three are not" in flat
    for c in ladder["cells"]:
        assert f"{c['ladder']['40']['auc']:.3f}" in flat, (
            f"{c['name']}'s k=40 figure is missing from results.md")
    assert "anti-correlated with robustness" in flat
    assert "does not say which kind of result it is" in flat


def test_both_documents_point_at_the_command():
    for doc in (PREREG, RESULTS_DOC):
        flat = " ".join(doc.read_text().split())
        assert "make prospective-analogy" in flat, doc.name


def test_the_amendment_retracts_nothing():
    flat = _amendment()
    assert "Nothing is retracted" in flat
    assert "H3 is closed" in flat


@needs_split
def test_it_strengthens_the_denv2_choice_for_a_stated_reason(ladder):
    """DENV-2 is the only cell that is most-enriched, powerable AND
    robust. The screen was registered on the first two; this is the
    third, and the amendment must not claim it was known earlier."""
    d2 = next(c for c in ladder["cells"] if c["name"] == "DENV-2")
    assert (d2["clears_up_to_k"] or -1) >= 40
    flat = _amendment()
    assert "strengthens the DENV-2 choice" in flat
    assert "was not known when it was signed" in flat


def test_the_diagnostic_writes_nothing_into_results():
    src = SCRIPT.read_text()
    assert "Not written to\ndata/results" in src or \
           "not written to data/results" in " ".join(src.split()).lower()
    assert not re.search(r'default\s*=\s*.*data.*results', src)
