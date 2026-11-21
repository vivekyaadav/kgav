"""doc2/registrations/prereg_temporal_structure_scorer.md holds only operands, derived.

A REGISTRATION IS WORTH WHAT ITS NUMBERS ARE WORTH. This one was signed
before any AUC against a post-2021 label existed, so everything in it is
a count off the split or a constant taken from another registration.
Both kinds are re-derived here: the counts from
`data/releases/v0.1-train2021`, the alpha and z from the same formula
that reproduces the two families already registered.

WHAT THIS FILE MUST NOT DO is check the document against itself. Every
assertion below starts from the split or from arithmetic, and the
document is the thing being tested.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from statistics import NormalDist

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

DOC = ROOT / "doc2" / "registrations" / "prereg_temporal_structure_scorer.md"
SPLIT = ROOT / "data" / "releases" / "v0.1-train2021"
RESULTS = ROOT / "data" / "results" / "temporal_structure_scorer.json"
MIN_CLASS = 10


@pytest.fixture(scope="module")
def doc() -> str:
    return DOC.read_text()


@pytest.fixture(scope="module")
def labels() -> dict:
    f = SPLIT / "TEST_LABELS.json"
    if not f.exists():
        pytest.skip("the 2021 split is not on this checkout")
    return json.loads(f.read_text())


@pytest.fixture(scope="module")
def names() -> dict[str, str]:
    vy = yaml.safe_load((ROOT / "config" / "viruses.yaml").read_text())
    return {f"NCBITaxon:{e['taxon']}": e["name"]
            for lst in ("viruses", "surveyed") for e in (vy.get(lst) or [])}


def _cells(labels) -> dict[str, tuple[int, int, int]]:
    out = {}
    for t in set(labels["test"]) | set(labels["train"]):
        out[t] = (len(labels["train"].get(t, [])),
                  len(labels["test"].get(t, [])),
                  len(labels["test_negatives"].get(t, [])))
    return out


def _evaluable(cells) -> dict[str, tuple[int, int, int]]:
    return {t: c for t, c in cells.items()
            if c[0] >= MIN_CLASS and c[1] >= MIN_CLASS and c[2] >= MIN_CLASS}


# ------------------------------------------------------------- the family
def test_the_family_size_is_the_evaluable_cell_count(labels, doc):
    """TEN, derived. The document writes the number down so a run that
    finds a different one refuses instead of dividing by what it found."""
    n = len(_evaluable(_cells(labels)))
    assert "THE FAMILY IS TEN" in doc or f"family is {n}" in doc.lower()
    assert n == 10, f"the split now yields {n} evaluable cells, not 10"


def test_alpha_and_z_follow_from_the_family(labels, doc):
    """And by the same formula that reproduces the two families already
    registered -- 14 -> 2.9137 and 3 -> 2.3940 -- so the number here is
    not a table lookup nobody can check."""
    n = len(_evaluable(_cells(labels)))
    alpha = 0.05 / n
    z = NormalDist().inv_cdf(1 - alpha / 2)
    assert round(NormalDist().inv_cdf(1 - (0.05 / 14) / 2), 4) == 2.9137
    assert round(NormalDist().inv_cdf(1 - (0.05 / 3) / 2), 4) == 2.3940
    assert f"{alpha:.6f}" in doc, f"alpha {alpha:.6f} is not in the document"
    assert f"{z:.4f}" in doc, f"z {z:.4f} is not in the document"


def test_every_row_of_the_cell_table_matches_the_split(labels, names, doc):
    """All fifteen rows, evaluable or not. A table that lists only the
    cells that made it cannot show what was excluded or why."""
    cells = _cells(labels)
    rows = re.findall(r"^    (\S[^\n]*?)\s+([\d,]+)\s+([\d,]+)\s+([\d,]+)"
                      r"\s+(YES|no)\s*$", doc, re.M)
    assert rows, "no cell table found in the document"
    byname = {names[t]: c for t, c in cells.items() if t in names}
    seen = set()
    for name, tr, tp, tn, mark in rows:
        name = name.strip()
        assert name in byname, f"{name} is not a register virus"
        want = byname[name]
        got = tuple(int(x.replace(",", "")) for x in (tr, tp, tn))
        assert got == want, f"{name}: document {got}, split {want}"
        ok = all(v >= MIN_CLASS for v in want)
        assert (mark == "YES") == ok, f"{name} marked {mark}"
        seen.add(name)
    assert seen == set(byname), (
        f"the table omits {sorted(set(byname) - seen)}")


def test_the_excluded_cells_are_excluded_on_negatives_not_positives(labels,
                                                                    names):
    """The document claims four of the five dropped cells lose their
    negatives rather than their positives, and names MERS-CoV's 6. If a
    rebuild changes that, the paragraph is wrong."""
    cells = _cells(labels)
    dropped = {names[t]: c for t, c in cells.items()
               if t in names and t not in _evaluable(cells)}
    on_neg = [n for n, (_tr, _tp, tn) in dropped.items() if tn < MIN_CLASS]
    assert len(on_neg) >= 4, dropped
    assert dropped["MERS-CoV"][2] == 6, dropped["MERS-CoV"]


def test_the_split_is_not_stale_against_the_release(doc):
    """A registration naming a release it was not built from is a
    registration about nothing."""
    f = SPLIT / "DERIVED_FROM.json"
    if not f.exists():
        pytest.skip("no split on this checkout")
    d = json.loads(f.read_text())
    sha = ((d.get("provenance") or {}).get("release") or {}).get("edges_sha")
    assert sha and sha in doc, (
        f"the split was built from edges_sha {sha}, which the document "
        f"does not name")
    assert d.get("cutoff") == 2021
    assert d.get("undated") == "include"


def test_the_reach_ceiling_the_document_cites_is_the_measured_one(doc):
    """95.2% is the reason this study exists and is quoted as fact."""
    land = ROOT / "data" / "results" / "temporal_2021_include.json"
    if not land.exists():
        pytest.skip("no temporal results on this checkout")
    d = json.loads(land.read_text())
    tot = un = 0
    for a in d.get("audits", []):
        m = re.search(r"ceiling: ([\d,]+)/([\d,]+) test compounds", a)
        if m:
            un += int(m.group(1).replace(",", ""))
            tot += int(m.group(2).replace(",", ""))
    assert tot, "no ceiling audits to derive from"
    assert f"{tot:,}" in doc and f"{un:,}" in doc
    assert f"{un / tot * 100:.1f}%" in doc


# -------------------------------------------------- it has not run yet
def test_the_document_claims_no_auc_while_none_exists(doc):
    """THE RULE THIS PROJECT RUNS ON. A registration may not acquire a
    score before the protocol it registers has produced one. The same
    guard is in tests/test_prereg_clinical_labels.py, and it is what
    makes 'signed before' mean anything."""
    if RESULTS.exists():
        pytest.skip("the protocol has run; the verdict guard takes over")
    body = doc.split("## Predictions", 1)[0]
    stray = [m for m in re.findall(r"\bAUC[^.\n]{0,40}?(0\.\d{3,4})", body)]
    allowed = {"0.4087", "0.7751", "0.6756", "0.48", "0.81"}
    assert not [s for s in stray if s not in allowed], (
        f"this document carries an AUC-shaped figure before the run: "
        f"{stray}")
