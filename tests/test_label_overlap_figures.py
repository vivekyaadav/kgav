"""The shared-PMID confound, held to the file it comes from.

WHAT THIS IS FOR. "92.2% of labelled compounds share a publication between
the evidence edge and the label it was supposed to predict" is the single
most load-bearing figure in docs/evaluation.md: it is the reason that
document calls its own M1 result partly self-scoring, and the reason
kg_audit_report.txt splits its two questions the way it does.

Until 2026-11-02 it existed nowhere but the console output of
scripts/diagnose_label_overlap.py, which was read-only. From there it had
been copied into SIX places in evaluation.md and one in the audit report,
and not one of them could be checked -- the same situation as every figure
the gtop release silently moved that day, except that this one underpins
the argument rather than a row of a table.

So the script writes a stamped results file and this holds both documents
to it. Both percentages are derived here from the two integers they are
quotients of: recording a rounded percentage in the file would put a
second rounding rule in a second place, which is how the correction-family
alphas came to be typed beside a correct division.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "evaluation.md"
REPORT = ROOT / "docs" / "kg_audit_report.txt"
RESULTS = ROOT / "data" / "results" / "label_overlap.json"
sys.path.insert(0, str(ROOT / "src"))

#: A RETRACTION IS NOT A CLAIM. evaluation.md keeps a table row comparing
#: the retired 94.9% (1,015/1,070) against the current figure, and a
#: document-wide search for a stale number would fire on the history --
#: forbidding the project from recording what it used to believe. Every
#: check below that forbids a figure is scoped past these.
HISTORY_MARKERS = ("An earlier version", "94.9% shared PMIDs")


def counts() -> dict:
    if not RESULTS.exists():
        pytest.skip(
            "data/results/label_overlap.json is not built in this checkout; "
            "run `python scripts/diagnose_label_overlap.py --out "
            "data/results/label_overlap.json`")
    return json.loads(RESULTS.read_text())


@pytest.fixture(scope="module")
def c() -> dict:
    return counts()


def shared_pct(c: dict) -> str:
    return f"{100 * c['share_a_pub'] / c['both_sides_have_pubs']:.1f}%"


def unwithholdable_pct(c: dict) -> str:
    return f"{100 * c['unwithholdable'] / c['both']:.1f}%"


def live(text: str) -> str:
    """The text with its history paragraphs removed, whitespace-normalised."""
    keep = [p for p in text.split("\n\n")
            if not any(m in p for m in HISTORY_MARKERS)]
    return " ".join("\n\n".join(keep).split())


# ------------------------------------------------- the figure, both halves
def test_the_shared_pmid_claim_is_the_files_own(c):
    """The sentence carries three numbers and all three are in the file:
    the count that share, the denominator, and the percentage."""
    flat = live(DOC.read_text())
    assert (f"{c['both_sides_have_pubs']:,} have a publication on" in flat), (
        f"{c['both_sides_have_pubs']:,} compounds carry publications on both "
        f"sides and the document states a different number")
    assert f"**{c['share_a_pub']:,} -- {shared_pct(c)} -- share a PMID**" in flat, (
        f"{c['share_a_pub']:,} of {c['both_sides_have_pubs']:,} share a PMID "
        f"= {shared_pct(c)}; the document says otherwise")


def test_the_denominator_of_the_claim_is_stated(c):
    """92.2% IS CONDITIONAL and the condition is the hole: it is the share
    among compounds where BOTH sides cite a publication, not among all
    labelled compounds. The script's own comment records quoting it
    without that qualifier as a past mistake, so the document has to
    carry the full population too."""
    flat = live(DOC.read_text())
    assert f"Of {c['both']:,} compounds carrying" in flat, (
        f"the claim's population is {c['both']:,} compounds with both an "
        f"evidence edge and a label; the document states a different number")


def test_the_compounds_outside_the_control_are_counted(c):
    """The ones --publication-disjoint can never withhold, because an edge
    with no PMID is not a candidate for withholding. A control's blind
    spot is part of the measurement."""
    flat = live(DOC.read_text())
    assert (f"**{c['unwithholdable']:,} compounds "
            f"({unwithholdable_pct(c)}) carry no publication") in flat, (
        f"{c['unwithholdable']:,} of {c['both']:,} = "
        f"{unwithholdable_pct(c)} carry no PMID on the evidence side")
    assert f"At {unwithholdable_pct(c)} that does not move row 4" in flat


@pytest.mark.parametrize("where", ["doc", "report"])
def test_every_other_mention_of_the_percentage_agrees(c, where):
    """SIX PLACES IN ONE DOCUMENT. The figure is referred to in passing
    four more times after the paragraph that derives it, plus once in the
    audit report, and a guard on the derivation alone leaves those to
    drift -- which is precisely how section 8 of the audit report came to
    hold four stale pathway figures while evaluation.md's copies were
    checked."""
    text = (DOC if where == "doc" else REPORT).read_text()
    flat = live(text)
    pct = shared_pct(c)
    mentions = re.findall(r"(\d{1,3}\.\d)% (?:of labelled|shared|share)", flat)
    assert mentions, f"no shared-PMID percentage found in {where}"
    wrong = sorted({m for m in mentions if f"{m}%" != pct})
    assert not wrong, (
        f"{where} states {wrong} as the shared-PMID share; the results "
        f"file gives {pct}. Re-derive every mention, not just the first.")


def test_the_summary_table_row_carries_both_figures(c):
    """The one-line summary of the confound, which is what a reader
    skimming takes away."""
    flat = live(DOC.read_text())
    assert (f"{shared_pct(c)} shared PMIDs "
            f"({c['share_a_pub']:,}/{c['both_sides_have_pubs']:,})") in flat
    assert (f"a further {c['unwithholdable']:,} compounds have no") in flat


# ----------------------------------------------- the file itself is honest
def test_the_counts_partition_the_population(c):
    """FOUR DISJOINT CASES, and they must sum. A compound with both an
    evidence edge and a label falls in exactly one: both sides cite a
    paper, only the label does, only the evidence does, or neither. If
    these stop summing, one branch is being counted twice and every
    percentage above is wrong in a way no document check would show."""
    parts = (c["both_sides_have_pubs"] + c["evidence_has_no_pmid"]
             + c["neither_has_a_pmid"] + c["label_has_no_pmid"])
    assert parts == c["both"], (
        f"the four cases sum to {parts:,} and the population is "
        f"{c['both']:,}")


def test_unwithholdable_is_the_two_cases_with_no_evidence_pmid(c):
    """Not "everything outside the 92.2% denominator": a compound whose
    LABEL carries no PMID still has a withholdable evidence edge. Folding
    it in would overstate the control's blind spot."""
    assert c["unwithholdable"] == (c["evidence_has_no_pmid"]
                                  + c["neither_has_a_pmid"])


def test_the_file_names_the_predicates_it_measured(c):
    """The overlap is INHIBITS against the positive label. gtop's 6,765
    new edges are TARGETS, which is why the 2026-11-02 release did not
    move this figure -- but that is a fact about which predicate was
    measured, and a file that does not record it leaves the next reader
    to assume as I did before checking."""
    assert c["evidence_predicate"] == "INHIBITS"
    assert c["label_predicate"] == "HAS_ANTIVIRAL_ACTIVITY_AGAINST"


def test_the_file_is_stamped_against_the_release(c):
    """Same provenance as every other results file, so `make
    check-results` reports it stale when the graph moves. A figure in six
    documents with no release fingerprint is the situation this file
    exists to end."""
    assert "provenance" in c
    prov = c["provenance"]
    assert prov["release"], "no release fingerprint"
    assert set(prov["code"]) == {"commit", "dirty"}


# -------------------------------------------------- the history is scoped
def test_the_retired_figure_is_still_recorded(c):
    """The project records what it used to believe, and the checks above
    must not forbid that. If this stops finding the retraction, `live()`
    is stripping the wrong paragraphs and every check above is reading a
    document with holes in it."""
    text = DOC.read_text()
    assert "94.9% shared PMIDs" in text, (
        "the retired figure is gone from the document, so the history "
        "scoping above is no longer exercised")
    assert "94.9%" not in live(text), (
        "the retraction leaked into the live text, so a stale figure "
        "there would pass the checks above")
