"""The strength-matched registration, checked against the runs it derives from.

Every number in docs/prereg_structure_scorer_strength_matched.md is a
function of the two prior runs: the fourteen k_cell values, their
feasibility, the three cells where the new control is weaker, P1's
threshold, the cell counts the predictions are calibrated between, and the
arithmetic showing the plate-matched cell cannot be strength-matched.

THIS IS THE THIRD REGISTRATION IN A CHAIN and the first two each had a
figure wrong on the first draft -- a pooled compound count quoted for a
ChEMBL-only study, and an asymmetry table that showed only the rows
supporting its own hypothesis. Both were caught by derivation, neither by
re-reading.

It computes no score. The protocol is signed and its code does not exist.
"""
from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

PREREG = ROOT / "docs/prereg_structure_scorer_strength_matched.md"
PUB = ROOT / "data/results/structure_scorer.json"
KNN = ROOT / "data/results/structure_scorer_symmetric.json"

pytestmark = pytest.mark.skipif(
    not (PUB.exists() and KNN.exists()),
    reason="both prior structure-scorer runs are needed")


@pytest.fixture(scope="module")
def doc() -> str:
    assert PREREG.exists(), "the registration is missing"
    return " ".join(PREREG.read_text().split())


@pytest.fixture(scope="module")
def runs() -> tuple[dict, dict]:
    return json.loads(PUB.read_text()), json.loads(KNN.read_text())


def cells(pub: dict) -> dict:
    return {v: c for v, c in pub["primary_chembl"].items()
            if c.get("readable")}


def rates(counts: dict) -> tuple[float, float]:
    return (counts["excluded_publication_pos"]
            / max(counts["candidates_pos"], 1),
            counts["excluded_publication_neg"]
            / max(counts["candidates_neg"], 1))


def k_cells(pub: dict) -> dict[str, tuple[int, int, int]]:
    """virus -> (k_cell, references after self-exclusion, what survives)."""
    out = {}
    for v, c in cells(pub).items():
        n = c["structure_max"]["counts"]
        ep, en = rates(n)
        refs = n["references"] - 1
        k = math.ceil(max(ep, en))
        out[v] = (k, refs, refs - k)
    return out


def test_it_says_it_was_signed_before_the_code_and_is_not_blind(doc):
    assert "SIGNED 2026-11-11, BEFORE THE CODE AND BEFORE THE RUN" in doc
    assert "NOT BLIND, AND FOR THE SECOND TIME" in doc, (
        "the third registration in a chain must say it responds to its "
        "predecessor's result")
    src = (ROOT / "scripts/structure_scorer.py").read_text()
    if "k_cell" in src:
        assert "# Verdict" in PREREG.read_text(), (
            "scripts/structure_scorer.py reads a per-cell k, so the run has "
            "happened and this document must carry its verdict")


def test_every_k_cell_is_derived_and_feasible(doc, runs):
    """THE FOURTEEN VALUES, AND THE FEASIBILITY THAT CHOSE THE FORMULA.

    A k that empties a reference set reproduces the first study's P5
    failure, where a cell printed 0.5000 because nothing was ranked.
    """
    import importlib
    D = importlib.import_module("diagnose_structure_scorer_design")
    names = D.virus_names()
    pub, _knn = runs
    table = k_cells(pub)
    for v, (k, refs, left) in table.items():
        assert left > 0, (
            f"{names.get(v, v)} holds {refs} references and k_cell={k} "
            f"empties it")
        assert f"{k:>7} {left:>6,}".strip() in doc or (
            f"{k} " in doc and f"{left:,}" in doc), (
            f"{names.get(v, v)}'s k_cell={k} leaving {left:,} is not in the "
            f"registration's table")
    tightest = min(left for _k, _r, left in table.values())
    assert f"keeps {tightest} references" in doc, (
        f"the tightest cell keeps {tightest} references and the "
        f"registration states otherwise")


def test_the_k_table_is_parsed_and_every_column_matches(doc, runs):
    """PARSED, not looked for anywhere in the document.

    The previous registration's asymmetry-table test passed through two
    worthless drafts -- one requiring only the rows that fit, one satisfied
    by the prose that described the defect -- before it parsed the table.
    This starts there.
    """
    import importlib
    D = importlib.import_module("diagnose_structure_scorer_design")
    names = D.virus_names()
    pub, knn = runs
    row = re.compile(
        r"^\s{4}(.+?)\s{2,}([\d,]+)\s+([\d.]+)\s+([\d.]+)\s+(\d+)\s+"
        r"([\d,]+)\s+([\d.]+)\s*$")
    parsed = {}
    for line in PREREG.read_text().splitlines():
        m = row.match(line)
        if m:
            parsed[m.group(1).strip()] = (
                int(m.group(2).replace(",", "")), float(m.group(3)),
                float(m.group(4)), int(m.group(5)),
                int(m.group(6).replace(",", "")), float(m.group(7)))
    assert parsed, "no k_cell table rows parsed out of the registration"

    table = k_cells(pub)
    expected = {}
    for v, (k, refs, left) in table.items():
        ep, en = rates(pub["primary_chembl"][v]["structure_max"]["counts"])
        expected[names.get(v, v)] = (
            refs, ep, en, k, left,
            knn["primary_chembl"][v]["structure_max"]["auc"])
    assert sorted(parsed) == sorted(expected), (
        f"table rows {sorted(parsed)} against cells {sorted(expected)}")
    for name, want in expected.items():
        got = parsed[name]
        for label, w, g in zip(
                ("refs", "pub/pos", "pub/neg", "k_cell", "left", "at k=10"),
                want, got):
            tol = 0.005 if isinstance(w, float) else 0
            assert abs(w - g) <= tol, (name, label, g, w)


def test_p1s_threshold_is_the_same_derived_one(doc, runs):
    pub, _knn = runs
    was = pub["primary_chembl"]["NCBITaxon:2697049"]["structure_max"]["auc"]
    half = 0.5 + (was - 0.5) / 2
    assert f"{half:.4f}" in doc, (
        f"half the publication-disjoint margin is {half:.4f} and P1 states "
        f"a different threshold")
    k = k_cells(pub)["NCBITaxon:2697049"][0]
    assert f"k_cell = {k} reads at or above" in doc, (
        f"SARS-CoV-2's k_cell is {k} and P1 names a different one")


def test_the_three_weaker_cells_are_named(doc, runs):
    """P3 is calibrated on how many cells get a sterner control, so the
    exceptions have to be named before the run rather than explained
    after it."""
    import importlib
    D = importlib.import_module("diagnose_structure_scorer_design")
    names = D.virus_names()
    pub, _knn = runs
    weaker = {names.get(v, v): k for v, (k, _r, _l) in k_cells(pub).items()
              if k < 10}
    words = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five"}
    assert f"< 10 in {words[len(weaker)]}" in doc, (
        f"{len(weaker)} cells get a weaker control and the registration "
        f"says otherwise")
    for name, k in weaker.items():
        assert f"{name} at {k}" in doc, (
            f"{name} gets k_cell={k}, weaker than the k=10 it is compared "
            f"with, and the registration does not name it")
    sterner = len(k_cells(pub)) - len(weaker)
    assert f">= 10 in {_word(sterner)} cells" in doc, (
        f"{sterner} cells get a sterner control and the registration says "
        f"otherwise")


def _word(n: int) -> str:
    return {0: "zero", 1: "one", 2: "two", 3: "three", 4: "four",
            5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine",
            10: "ten", 11: "eleven", 12: "twelve", 13: "thirteen",
            14: "fourteen"}[n]


def test_p5_is_calibrated_between_the_two_prior_counts(doc, runs):
    """P5 discriminates the two explanations, so its threshold must sit
    between the counts each predicts."""
    pub, knn = runs

    def above(run: dict) -> int:
        return sum(1 for c in run["primary_chembl"].values()
                   if c.get("readable")
                   and c["structure_max"].get("clears_corrected_alpha")
                   and c["structure_max"].get("direction") == "above")
    n_pub, n_knn = above(pub), above(knn)
    m = re.search(r"At least (\d+) of the 14 cells clear", doc)
    assert m, "P5 no longer states a threshold"
    threshold = int(m.group(1))
    assert n_pub < threshold < n_knn, (
        f"P5's threshold {threshold} does not separate the publication "
        f"control's {n_pub} from the symmetric control's {n_knn}, so it "
        f"discriminates nothing")
    assert f"symmetric study's {_word(n_knn)}" in doc
    assert f"control's {_word(n_pub)}" in doc


def test_the_plate_matched_cell_cannot_be_strength_matched(doc, runs):
    """The arithmetic behind the registration's one declared impossibility.
    It is a fact about the PRIOR control -- it removed everything there --
    and if that ever stops being true the paragraph is wrong."""
    pub, _knn = runs
    n = pub["secondary_ncats"]["NCBITaxon:2697049"]["structure_max"]["counts"]
    ep, en = rates(n)
    refs = n["references"] - 1
    k = math.ceil(max(ep, en))
    assert refs - k < 0, (
        f"ceil(max({ep:.2f}, {en:.2f})) = {k} against {refs} references "
        f"no longer empties the cell; strength-matching is possible and the "
        f"registration says it is not")
    assert f"{ep:.2f} references per positive and" in doc
    assert f"{en:.2f} per negative there, out of {refs}" in doc
    assert f"is {k} and leaves {refs - k}" in doc, (
        f"the arithmetic is {k} leaving {refs - k} and the registration "
        f"states something else")
    # and the two k values it is registered to run at must be real
    primary_k = k_cells(pub)["NCBITaxon:2697049"][0]
    assert f"k = 10 and k = {primary_k}" in doc, (
        f"the plate-matched cell should run at the ChEMBL cell's k_cell "
        f"({primary_k}) for comparability and the registration names "
        f"another value")
    assert refs - primary_k > 0, (
        f"k={primary_k} leaves {refs - primary_k} in the plate-matched "
        f"cell, which is not runnable either")


def test_the_confound_bound_it_inherits_is_the_prior_runs(doc, runs):
    """The three cells and three deltas carried over from the symmetric
    verdict, and the 0.22 collection-confound bound."""
    pub, knn = runs
    deltas = []
    for v, c in cells(pub).items():
        ep, _en = rates(c["structure_max"]["counts"])
        if ep < 10:
            deltas.append(knn["primary_chembl"][v]["structure_max"]["auc"]
                          - c["structure_max"]["auc"])
    assert len(deltas) == 3, f"{len(deltas)} cells, not three"
    for d in sorted(deltas):
        assert f"{d:+.4f}" in doc, (
            f"a carried-over delta is {d:+.4f} and the registration quotes "
            f"something else")
    chembl = knn["primary_chembl"]["NCBITaxon:2697049"]["structure_max"]["auc"]
    matched = knn["secondary_ncats"]["NCBITaxon:2697049"]["structure_max"]["auc"]
    assert f"{chembl - matched:.2f} AUC" in doc, (
        f"the collection-confound bound is {chembl - matched:.2f} and the "
        f"registration states otherwise")


def test_mers_is_carried_forward_as_unexplained(doc, runs):
    """The one cell that fell. A registration that stopped mentioning it
    would be the silent-drop this project has relabelled twice."""
    pub, knn = runs
    was = pub["primary_chembl"]["NCBITaxon:1335626"]["structure_max"]["auc"]
    now = knn["primary_chembl"]["NCBITaxon:1335626"]["structure_max"]["auc"]
    assert now < was, "MERS-CoV no longer fell; the paragraph is stale"
    assert f"{was:.4f} to {now:.4f}" in doc
    k, _refs, _left = k_cells(pub)["NCBITaxon:1335626"]
    assert f"`k_cell` is {k}" in doc, (
        f"MERS-CoV's k_cell is {k} and the registration names another")
    n_neg = knn["primary_chembl"]["NCBITaxon:1335626"]["structure_max"]["n_neg"]
    assert f"{n_neg} negatives" in doc
    # AND ITS RANK, DERIVED. The first draft of the registration called this
    # "the smallest negative class in the family"; JEV's 27 is smaller, and
    # JEV moved the opposite way, which is why the claim mattered -- it was
    # the reading that made MERS-CoV look like small-sample noise.
    order = sorted(c["structure_max"]["n_neg"] for c in cells(pub).values())
    rank = order.index(n_neg) + 1
    ordinal = {1: "smallest", 2: "second-smallest", 3: "third-smallest"}
    assert ordinal.get(rank, f"rank {rank}") in doc, (
        f"MERS-CoV has the {ordinal.get(rank, rank)} negative class "
        f"({n_neg}) and the registration describes it differently")


def test_the_family_and_z_are_unchanged_and_derived(doc, runs):
    pub, _knn = runs
    m = re.search(r"alpha = 0\.05/(\d+) = (\d+\.\d+), z = (\d+\.\d+)", doc)
    assert m, "the registration no longer states alpha and z together"
    family, alpha, z = int(m.group(1)), float(m.group(2)), float(m.group(3))
    assert family == len(cells(pub)), (family, len(cells(pub)))
    exact = 0.05 / family
    assert abs(alpha - exact) < 5e-7, (alpha, exact)
    lo, hi = 0.0, 10.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if math.erfc(mid / math.sqrt(2)) > exact:
            lo = mid
        else:
            hi = mid
    assert abs(z - (lo + hi) / 2) < 5e-5, (z, (lo + hi) / 2)
