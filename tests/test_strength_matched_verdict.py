"""The strength-matched verdict, re-derived from the three runs it compares.

Everything in docs/prereg_structure_scorer_strength_matched.md's verdict is
a function of the three results files: the 14-row table, the six prediction
outcomes, the six cells that broke P3 and which of them were named in
advance, the +0.912 ordering correlation, and the plate-matched cell that
stopped clearing.

THE CLAIM THIS FILE MOST NEEDS TO PIN IS A RETRACTION. The previous
verdict's headline -- a channel clearing a corrected threshold on a
plate-matched population -- does not survive k=33, and the paragraph saying
so has to keep matching the data. A retraction that silently became true
again, or silently stopped being warranted, would be worse than the
original claim.
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
SM = ROOT / "data/results/structure_scorer_strength_matched.json"

pytestmark = pytest.mark.skipif(
    not (PUB.exists() and KNN.exists() and SM.exists()),
    reason="all three structure-scorer runs are needed")


@pytest.fixture(scope="module")
def doc() -> str:
    return " ".join(PREREG.read_text().split())


@pytest.fixture(scope="module")
def runs() -> tuple[dict, dict, dict]:
    return (json.loads(PUB.read_text()), json.loads(KNN.read_text()),
            json.loads(SM.read_text()))


def cells(sm: dict) -> dict:
    return {v: c for v, c in sm["primary_chembl"].items()
            if c.get("readable")}


def _spearman(xs, ys) -> float:
    def rank(zs):
        order = sorted(range(len(zs)), key=lambda i: zs[i])
        out = [0] * len(zs)
        for pos, i in enumerate(order):
            out[i] = pos + 1
        return out
    a, b = rank(xs), rank(ys)
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    return (sum((p - ma) * (q - mb) for p, q in zip(a, b))
            / math.sqrt(sum((p - ma) ** 2 for p in a)
                        * sum((q - mb) ** 2 for q in b)))


def test_p2_the_invariant_holds_at_every_cells_own_k(doc, runs):
    """Checked against each cell's k_cell, not against a single k. A
    version comparing everything with 10 would report thirteen of fourteen
    violated and say nothing about the control."""
    _pub, _knn, sm = runs
    for v, c in cells(sm).items():
        k = c["k_cell"]
        n = c["structure_max"]["counts"]
        assert n["excluded_knn_pos"] == k * n["candidates_pos"], (v, "pos")
        assert n["excluded_knn_neg"] == k * n["candidates_neg"], (v, "neg")
        assert n["knn_short"] == 0, (v, "ran short of references")
    m = sm["secondary_ncats"]["NCBITaxon:2697049"]["structure_max"]["counts"]
    k = sm["secondary_ncats"]["NCBITaxon:2697049"]["k_cell"]
    assert m["excluded_knn_pos"] == k * m["candidates_pos"]
    assert "HELD in all 14 cells and the matched one" in doc


def test_p1_rules_out_leniency(doc, runs):
    """The figure, its threshold, and the comparison that makes it mean
    something: k_cell must exceed what the old control removed from both
    classes, or 'sterner on both sides' is false."""
    pub, _knn, sm = runs
    c = cells(sm)["NCBITaxon:2697049"]
    got = c["structure_max"]["auc"]
    was = pub["primary_chembl"]["NCBITaxon:2697049"]["structure_max"]["auc"]
    threshold = 0.5 + (was - 0.5) / 2
    assert got >= threshold, (got, threshold)
    assert f"{got:.4f} [{c['structure_max']['auc_ci_lo_corrected']:.4f}, " \
           f"{c['structure_max']['auc_ci_hi_corrected']:.4f}]" in doc
    n = pub["primary_chembl"]["NCBITaxon:2697049"]["structure_max"]["counts"]
    ep = n["excluded_publication_pos"] / max(n["candidates_pos"], 1)
    en = n["excluded_publication_neg"] / max(n["candidates_neg"], 1)
    k = c["k_cell"]
    assert k >= ep and k >= en, (
        f"k_cell={k} does not exceed both {ep:.2f}/pos and {en:.2f}/neg, so "
        f"the verdict's 'strictly sterner on both sides' is false")
    assert f"{ep:.2f} per positive and {en:.2f} per negative" in doc


def test_p3_failed_and_the_count_and_cells_are_the_runs(doc, runs):
    """THE PREDICTION THAT FAILED. Its count, the six cells that broke it,
    and which were named before the run."""
    import importlib
    D = importlib.import_module("diagnose_structure_scorer_design")
    names = D.virus_names()
    _pub, _knn, sm = runs
    companion = sm["k_companion"]
    lower, higher = 0, []
    for v, c in cells(sm).items():
        a = c["structure_max"]["auc"]
        b = c[f"structure_max_k{companion}"]["auc"]
        if a <= b:
            lower += 1
        else:
            higher.append((names.get(v, v), c["k_cell"], a - b))
    assert f"**FAILED — {lower} of {len(cells(sm))}**" in doc, (
        f"{lower} of {len(cells(sm))} cells fell and the verdict says "
        f"otherwise")
    assert lower < 10, (
        "P3's registered threshold of 10 is now met; the verdict calls it "
        "FAILED and would need rewriting")
    # THE RISER TABLE, PARSED. Substring matching failed here too: the
    # rows right-align the k_cell as "k_cell= 4", so looking for
    # "West Nile k_cell=4" found nothing while the row was correct. That
    # unanchored-substring failure has now happened four times in this
    # chain of registrations, so this parses.
    riser = re.compile(
        r"^\s{4}(.+?)\s{2,}k_cell=\s*(\d+)\s+([+-][\d.]+)\s+"
        r"k [<>] (\d+), (.+?)\s*$")
    parsed = {}
    for line in PREREG.read_text().splitlines():
        m = riser.match(line)
        if m:
            parsed[m.group(1).strip()] = (
                int(m.group(2)), float(m.group(3)), int(m.group(4)),
                m.group(5).strip())
    assert sorted(parsed) == sorted(n for n, _k, _d in higher), (
        f"the riser table lists {sorted(parsed)} and the run gives "
        f"{sorted(n for n, _k, _d in higher)}")
    for name, k, delta in higher:
        p_k, p_delta, p_comp, p_label = parsed[name]
        assert p_k == k, (name, "k_cell", p_k, k)
        assert abs(p_delta - delta) < 5e-5, (name, "delta", p_delta, delta)
        assert p_comp == companion, (name, "companion", p_comp, companion)
        want = "named in advance" if k < companion else "NOT anticipated"
        assert p_label == want, (
            f"{name} has k_cell={k} against a companion of {companion}, so "
            f"its row should read '{want}' and reads '{p_label}'")
    named = [h for h in higher if h[1] < companion]
    unnamed = [h for h in higher if h[1] >= companion]
    words = {1: "One", 2: "Two", 3: "Three", 4: "Four", 5: "Five"}
    assert f"{words[len(named)]} were named before the run" in doc, (
        f"{len(named)} risers were named in advance and the verdict says "
        f"otherwise")
    assert f"other {words[len(unnamed)].lower()} are the failure" in doc, (
        f"{len(unnamed)} risers were unanticipated and the verdict says "
        f"otherwise")


def test_p4_and_p5_are_the_discriminating_counts(doc, runs):
    """P5 is what separates the two explanations, so all three counts --
    this run's, the symmetric study's and the publication control's -- have
    to be the files'."""
    pub, knn, sm = runs

    def counts(run: dict) -> tuple[int, int]:
        above = below = 0
        for c in run["primary_chembl"].values():
            if not c.get("readable"):
                continue
            m = c["structure_max"]
            if not m.get("clears_corrected_alpha"):
                continue
            if m.get("direction") == "above":
                above += 1
            else:
                below += 1
        return above, below

    a_sm, b_sm = counts(sm)
    a_knn, _ = counts(knn)
    a_pub, _ = counts(pub)
    assert b_sm == 0, f"{b_sm} cell(s) clear below chance; P4 failed"
    assert "**HELD — zero**" in doc
    assert f"**HELD — {a_sm}, the symmetric study's count, not the " \
           f"publication control's {a_pub}**" in doc, (
        f"this run clears in {a_sm} cells, the symmetric study in {a_knn} "
        f"and the publication control in {a_pub}; the verdict's P5 row "
        f"states a different combination")
    assert a_sm == a_knn, (
        f"this run clears {a_sm} and the symmetric study {a_knn}; the "
        f"verdict claims they are the same count")
    assert a_sm >= 8, "P5's registered threshold is not met"


def test_p6s_ordering_correlation(doc, runs):
    _pub, knn, sm = runs
    xs, ys = [], []
    for v, c in cells(sm).items():
        xs.append(c["structure_max"]["auc"])
        ys.append(knn["primary_chembl"][v]["structure_max"]["auc"])
    rho = _spearman(xs, ys)
    assert f"**{rho:+.3f}**" in doc, (
        f"the ordering correlation is {rho:+.3f} and the verdict says "
        f"otherwise")
    assert abs(rho) > 0.544


def test_the_retraction_is_warranted_and_stated(doc, runs):
    """THE CLAIM BEING WITHDRAWN. The previous verdict said a channel
    clears a corrected threshold on a plate-matched population. Here it
    does not, and if that ever changes this paragraph is wrong in the
    other direction."""
    _pub, knn, sm = runs
    now = sm["secondary_ncats"]["NCBITaxon:2697049"]["structure_max"]
    before = knn["secondary_ncats"]["NCBITaxon:2697049"]["structure_max"]
    assert before["clears_corrected_alpha"] is True, (
        "the k=10 matched cell no longer cleared, so there is nothing to "
        "retract and the verdict's account of the previous one is wrong")
    assert now["clears_corrected_alpha"] is False, (
        "the matched cell clears at k=33 after all; the verdict retracts a "
        "claim that stands")
    assert (f"{now['auc']:.4f} [{now['auc_ci_lo_corrected']:.4f}, "
            f"{now['auc_ci_hi_corrected']:.4f}]" in doc)
    assert f"At k=10 it read {before['auc']:.4f} and cleared" in doc
    # THE HANDICAP IS READ FROM THE RUN, NOT TYPED. This assertion and the
    # confound-gap one below both carried a literal "k=33", and the L27
    # rebuild moved the matched cell's k_cell to 36 -- so a test written to
    # catch a stale figure was itself the stale figure. That is the "check
    # written with the same literals it is checking" defect, in a file whose
    # docstring is about exactly that.
    k_m = sm["secondary_ncats"]["NCBITaxon:2697049"]["k_cell"]
    assert f"skill {now['auc_skill']:.4f} at k={k_m}" in doc, (
        f"the matched cell's skill is {now['auc_skill']:.4f} and the verdict "
        f"states otherwise")
    assert now["auc"] > 0.5, (
        "the matched cell no longer reads above chance at all; the verdict's "
        "weaker surviving claim is also gone")
    assert "should be read as retracted" in doc


def test_the_confound_gap_at_equal_handicap(doc, runs):
    _pub, knn, sm = runs
    gap_now = (cells(sm)["NCBITaxon:2697049"]["structure_max"]["auc"]
               - sm["secondary_ncats"]["NCBITaxon:2697049"]
               ["structure_max"]["auc"])
    gap_then = (knn["primary_chembl"]["NCBITaxon:2697049"]
                ["structure_max"]["auc"]
                - knn["secondary_ncats"]["NCBITaxon:2697049"]
                ["structure_max"]["auc"])
    k_m = sm["secondary_ncats"]["NCBITaxon:2697049"]["k_cell"]
    assert f"**{gap_now:.4f} at k={k_m}**" in doc, (
        f"the gap at equal handicap is {gap_now:.4f} and the verdict says "
        f"otherwise")
    assert f"against {gap_then:.2f} measured at k=10" in doc


def test_mers_is_carried_across_all_three_runs(doc, runs):
    """The one cell no registration explains, with all three of its
    figures, so it cannot quietly drop out of the third verdict."""
    pub, knn, sm = runs
    v = "NCBITaxon:1335626"
    a = pub["primary_chembl"][v]["structure_max"]["auc"]
    b = knn["primary_chembl"][v]["structure_max"]["auc"]
    c = cells(sm)[v]
    assert (f"{a:.4f} publication, {b:.4f} at k=10, "
            f"{c['structure_max']['auc']:.4f} at k_cell={c['k_cell']}"
            in doc), (
        f"MERS-CoV reads {a:.4f} / {b:.4f} / "
        f"{c['structure_max']['auc']:.4f} and the verdict states otherwise")


def test_the_table_carries_every_cell_and_column(doc, runs):
    """Parsed, with every column compared -- the shape the previous two
    registrations' table tests only reached after several worthless
    drafts."""
    import importlib
    D = importlib.import_module("diagnose_structure_scorer_design")
    names = D.virus_names()
    pub, _knn, sm = runs
    companion = sm["k_companion"]
    row = re.compile(
        r"^\s{4}(.+?)\s{2,}([\d,]+)\s+([\d,]+)\s+(\d+)\s+([\d.]+)[\^v ]\s+"
        r"\[[\d.]+, [\d.]+\]\s+(-?[\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$")
    table = {}
    for line in PREREG.read_text().splitlines():
        m = row.match(line)
        if m:
            table[m.group(1).strip()] = (
                int(m.group(2).replace(",", "")),
                int(m.group(3).replace(",", "")),
                int(m.group(4)), float(m.group(5)), float(m.group(6)),
                float(m.group(7)), float(m.group(8)), float(m.group(9)))
    assert table, "no verdict table rows parsed"
    expected = {}
    for v, c in cells(sm).items():
        m = c["structure_max"]
        expected[names.get(v, v)] = (
            m["n_pos"], m["n_neg"], c["k_cell"], m["auc"], m["auc_skill"],
            c[f"structure_max_k{companion}"]["auc"],
            pub["primary_chembl"][v]["structure_max"]["auc"],
            c["degree"]["auc"])
    assert sorted(table) == sorted(expected), (
        f"table rows {sorted(table)} against cells {sorted(expected)}")
    labels = ("pos", "neg", "k", "auc", "skill", f"k={companion}", "prior",
              "degree")
    for name, want in expected.items():
        for lab, w, g in zip(labels, want, table[name]):
            tol = 5e-5 if isinstance(w, float) else 0
            assert abs(w - g) <= tol, (name, lab, g, w)
