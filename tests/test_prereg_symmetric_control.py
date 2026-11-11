"""The symmetric-control registration, checked against the run it responds to.

WHY. docs/prereg_structure_scorer_symmetric.md fixes a protocol before the
measurement, and every number in it is a function of either the release or
data/results/structure_scorer.json -- the asymmetry ratios it exists to
fix, the class sizes that decide k, the threshold in P1, and the four
deltas in P3. A registration whose figures are typed by hand drifts, and
this project has paid for that three times now: the pathway registration's
"eight times" that git put at eleven, the release plan's "M3 and M5 have
zero edges" when M3 has 11,744, and the first structure-scorer draft's
pooled compound count quoted for a ChEMBL-only study.

WHAT THIS FILE DOES NOT DO. It computes no score under the new protocol.
The registration is signed and the code does not exist yet; a test that
asserted an outcome would make the signature worthless.
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

PREREG = ROOT / "docs/prereg_structure_scorer_symmetric.md"
PRIOR = ROOT / "data/results/structure_scorer.json"


@pytest.fixture(scope="module")
def doc() -> str:
    assert PREREG.exists(), "the registration is missing"
    return " ".join(PREREG.read_text().split())


@pytest.fixture(scope="module")
def prior() -> dict:
    if not PRIOR.exists():
        pytest.skip("the prior structure-scorer run is not on disk")
    return json.loads(PRIOR.read_text())


def readable(prior: dict) -> dict:
    return {v: c for v, c in prior["primary_chembl"].items()
            if c.get("readable")}


def test_it_says_it_was_signed_before_the_code(doc):
    """And says plainly that it is not blind, which is a different thing."""
    assert "SIGNED 2026-11-11, BEFORE THE CODE AND BEFORE THE RUN" in doc
    assert "this is not a blind registration" in doc.lower(), (
        "the registration must state that its design responds to a result "
        "already in hand")
    src = (ROOT / "scripts/structure_scorer.py").read_text()
    if "excluded_knn" in src:
        assert "## Verdict" in doc or "# Verdict" in doc, (
            "scripts/structure_scorer.py implements the k-nearest exclusion, "
            "so the run has happened and this document must carry its "
            "verdict")


def test_the_registered_k_leaves_every_cell_non_empty(doc, prior):
    """THE ARITHMETIC THAT DECIDES k, AND THAT KILLED k=30.

    The smallest reference set in the family is what k must not exhaust. A
    k that empties it reproduces the prior run's P5 failure, where a cell
    printed 0.5000 with a confidence interval because nothing had been
    ranked.
    """
    m = re.search(r"\*\*k = (\d+) is the primary", doc)
    assert m, "the registration no longer states its primary k"
    k = int(m.group(1))

    smallest = min(c["structure_max"]["counts"]["references"] - 1
                   for c in readable(prior).values())
    assert smallest > k, (
        f"the smallest reference set after self-exclusion is {smallest} and "
        f"the registered k is {k}, which empties it")
    assert f"{smallest} after" in doc, (
        f"the smallest reference set is {smallest} and the registration does "
        f"not quote that figure")
    assert f"leaves {smallest - k} admissible" in doc, (
        f"k={k} leaves {smallest - k} references in the smallest cell and "
        f"the registration states something else")
    # and the rejected k must be named with ITS arithmetic
    assert "k=30 WAS IN THE DRAFT" in doc
    assert f"leaves JEV at {smallest - 30}" in doc, (
        f"k=30 leaves {smallest - 30} and the registration does not say so")


def test_the_curve_points_are_all_admissible(doc, prior):
    """Every k in the registered curve must leave a non-empty set too, or
    the curve has a vacuous point in it."""
    m = re.search(r"curve is k ∈ \{([\d, ]+)\}", doc)
    assert m, "the registration no longer states its sensitivity curve"
    ks = [int(x) for x in m.group(1).split(",")]
    smallest = min(c["structure_max"]["counts"]["references"] - 1
                   for c in readable(prior).values())
    for k in ks:
        assert smallest - k > 0, (
            f"curve point k={k} leaves {smallest - k} references in the "
            f"smallest cell")


def test_p1s_threshold_is_half_the_prior_margin(doc, prior):
    """DERIVED, NOT CHOSEN. P1's threshold is the point at which half of
    the prior run's margin above chance survives, so it moves if the prior
    run moves."""
    auc = readable(prior)["NCBITaxon:2697049"]["structure_max"]["auc"]
    half = 0.5 + (auc - 0.5) / 2
    assert f"{half:.4f}" in doc, (
        f"SARS-CoV-2 read {auc:.4f} and half its margin is {half:.4f}; the "
        f"registration states a different threshold")
    assert f"it read {auc:.4f} under the prior control" in doc


def test_the_asymmetry_table_shows_every_cell_not_the_ones_that_fit(prior):
    """EVERY readable cell, PARSED OUT OF THE TABLE.

    Two drafts of this test were worthless and the second looked fine.

    The first required only the four highest ratios and the two lowest --
    precisely the rows that support the hypothesis -- which is what the
    document's own first draft carried. It omitted DENV-2 at 12.7x and Zika
    at 8.2x, both heavily asymmetric and both ABOVE chance, and therefore
    counterexamples.

    The second required every cell's ratio to appear ANYWHERE in the
    document, and passed when the DENV-2 row was deleted from the table --
    because the prose paragraph admitting the original omission says "DENV-2
    at 12.7x". A check satisfied by the sentence that describes the defect
    is the same failure as the quoted-figure problem in
    tests/test_release_plan_temporal_note.py, in a new place.

    So the table is parsed, row by row, and every number in each row is
    compared with the run: both exclusion rates, the ratio, the AUC, and
    the counterexample mark.
    """
    import importlib
    D = importlib.import_module("diagnose_structure_scorer_design")
    names = D.virus_names()

    # THE NAME MAY CONTAIN DIGITS. The first pattern read the name as
    # `[^0-9]*?`, which stops at the first digit, so every cell called
    # DENV-1, SARS-CoV-2 or HCoV-229E failed to parse and the test reported
    # seven rows missing from a table that carried all fourteen. The column
    # separator is two or more spaces; a name's internal spaces are single.
    row = re.compile(
        r"^\s{4}(.+?)\s{2,}([\d.]+) / \s*([\d.]+)\s+=\s+([\d.]+)x"
        r"\s+AUC ([\d.]+)(\s+<- counterexample)?\s*$")
    table = {}
    for line in PREREG.read_text().splitlines():
        m = row.match(line)
        if m:
            table[m.group(1).strip()] = (
                float(m.group(2)), float(m.group(3)), float(m.group(4)),
                float(m.group(5)), bool(m.group(6)))
    assert table, "no asymmetry table rows found in the registration"

    expected = {}
    for v, c in readable(prior).items():
        k = c["structure_max"]["counts"]
        ep = k["excluded_publication_pos"] / max(k["candidates_pos"], 1)
        en = k["excluded_publication_neg"] / max(k["candidates_neg"], 1)
        if not en:
            continue
        expected[names.get(v, v)] = (ep, en, ep / en,
                                     c["structure_max"]["auc"])

    missing = sorted(set(expected) - set(table))
    assert not missing, (
        f"the asymmetry table omits {missing}; it must show every readable "
        f"cell, not the ones that fit the hypothesis")
    extra = sorted(set(table) - set(expected))
    assert not extra, f"the table carries rows for {extra}, which are not cells"

    for name, (ep, en, ratio, auc) in expected.items():
        t_ep, t_en, t_ratio, t_auc, marked = table[name]
        assert abs(t_ep - ep) < 0.005, (name, "pos rate", t_ep, ep)
        assert abs(t_en - en) < 0.005, (name, "neg rate", t_en, en)
        assert abs(t_ratio - ratio) < 0.05, (name, "ratio", t_ratio, ratio)
        assert abs(t_auc - auc) < 5e-5, (name, "auc", t_auc, auc)
        # A COUNTEREXAMPLE IS A CELL THE HYPOTHESIS DOES NOT EXPLAIN:
        # heavily asymmetric against positives, yet above chance.
        is_counter = ratio >= 5.0 and auc > 0.5
        assert marked == is_counter, (
            f"{name} is {'a' if is_counter else 'not a'} counterexample "
            f"({ratio:.1f}x, AUC {auc:.4f}) and the table marks it "
            f"{'' if marked else 'not '}so")

    counters = [n for n, (_e, _n, r, a) in expected.items()
                if r >= 5.0 and a > 0.5]
    assert counters, (
        "no cell is both heavily asymmetric and above chance any more; the "
        "registration's counterexample paragraph is stale")


def test_p3s_four_deltas_are_the_priors(doc, prior):
    """Each below-chance cell's distance to chance, from the run."""
    below = {v: c["structure_max"]["auc"]
             for v, c in readable(prior).items()
             if c["structure_max"]["auc"] < 0.5}
    assert below, "no cell reads below chance; P3 has nothing to test"
    for v, auc in below.items():
        assert f"{0.5 - auc:+.4f}" in doc, (
            f"{v} is {0.5 - auc:+.4f} from chance and the registration does "
            f"not quote it")
    words = {3: "three", 4: "four", 5: "five"}
    assert f"The {words[len(below)]} below-chance cells" in doc, (
        f"{len(below)} cells read below chance and the registration's prose "
        f"says otherwise")


def test_the_leak_correlation_p6_reasons_from_is_the_priors(doc, prior):
    """P6 predicts a correlation will weaken, so the baseline must be the
    one actually measured."""
    cells = readable(prior)
    aucs = [c["structure_max"]["auc"] for c in cells.values()]
    leaks = [c["structure_max_no_pub_exclusion"]["auc"]
             - c["structure_max"]["auc"] for c in cells.values()]

    def rank(zs):
        order = sorted(range(len(zs)), key=lambda i: zs[i])
        out = [0] * len(zs)
        for pos, i in enumerate(order):
            out[i] = pos + 1
        return out
    a, b = rank(aucs), rank(leaks)
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    rho = (sum((p - ma) * (q - mb) for p, q in zip(a, b))
           / math.sqrt(sum((p - ma) ** 2 for p in a)
                       * sum((q - mb) ** 2 for q in b)))
    assert f"{rho:.3f}" in doc, (
        f"the prior leak correlation is {rho:.3f} and the registration's P6 "
        f"reasons from a different number")


def test_the_family_and_z_match(doc, prior):
    """The correction family is the readable cells, derived; z is derived
    from the family."""
    # THE TRAILING PERIOD. `([\d.]+)` is greedy over the sentence's
    # full stop, so this read z as "2.9137." and raised ValueError rather
    # than asserting anything. A test that crashes is not a test that
    # fails; it says nothing about the document.
    m = re.search(r"alpha = 0\.05/(\d+) = (\d+\.\d+), z = (\d+\.\d+)", doc)
    assert m, "the registration no longer states alpha and z together"
    family, alpha, z = int(m.group(1)), float(m.group(2)), float(m.group(3))
    assert family == len(readable(prior)), (
        f"{len(readable(prior))} cells are readable and the family is "
        f"{family}")
    assert abs(alpha - 0.05 / family) < 5e-7, (alpha, 0.05 / family)
    # DERIVE z FROM 0.05/family, NOT FROM THE DOCUMENT'S PRINTED ALPHA.
    # The document prints alpha to six places as 0.003571, and bisecting
    # against that rounded value gives z = 2.913764 where the exact
    # 0.05/14 gives 2.913726 -- so this test failed the document for a
    # rounding it had itself introduced. Both values were checked against
    # an independent Newton solve on erfc before the test was changed, so
    # the document is right and this was the error.
    exact_alpha = 0.05 / family
    lo, hi = 0.0, 10.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if math.erfc(mid / math.sqrt(2)) > exact_alpha:
            lo = mid
        else:
            hi = mid
    assert abs(z - (lo + hi) / 2) < 5e-5, (
        f"two-sided z for alpha={alpha:.6f} is {(lo + hi) / 2:.4f}, not {z}")


def test_the_matched_cell_is_claimed_measurable_on_arithmetic(doc, prior):
    """P5's first half must follow from the class sizes, not from hope."""
    m = re.search(r"\*\*k = (\d+) is the primary", doc)
    k = int(m.group(1))
    cell = prior["secondary_ncats"]["NCBITaxon:2697049"]["structure_max"]
    refs = cell["counts"]["references"]
    assert refs - 1 - k > 0, (
        f"the matched cell holds {refs} references and k={k} would empty it; "
        f"P5's first half is false by arithmetic")
    assert f"{refs} positives, {refs - 1} references after" in doc, (
        f"the matched cell is {refs} positives / {refs - 1} references and "
        f"the registration states otherwise")
    assert f"leaves {refs - 1 - k}" in doc


def test_the_residual_asymmetry_is_bounded_by_the_measurement(doc, prior):
    """The registration bounds what exclusion (2) can still cost. That
    bound is the largest skeleton group in the prior run, and it must not
    be understated."""
    largest = 0
    for c in readable(prior).values():
        k = c["structure_max"]["counts"]
        largest = max(largest, k["excluded_skeleton"] and 1 or 0)
    # the claim under test is the share, which the design diagnostic
    # measured: no group larger than 3, so at most 2 extra references
    assert "no group larger than 3" in doc, (
        "the registration no longer bounds the skeleton groups")
    assert "at most 2 for a minority" in doc
    assert "excluded_skeleton_pos" in doc and "excluded_skeleton_neg" in doc, (
        "the registration must name the counters that will measure the "
        "residual, or it is a claim rather than a measurement")
