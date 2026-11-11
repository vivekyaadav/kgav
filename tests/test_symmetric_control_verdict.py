"""The symmetric-control verdict, re-derived from the run it describes.

Every figure in docs/prereg_structure_scorer_symmetric.md's verdict is a
function of data/results/structure_scorer_symmetric.json and its
comparator data/results/structure_scorer.json: the 14-row curve table, the
six prediction outcomes, the four P3 movements, the two correlations, and
the three cells that bound the strength confound.

THE CORRELATION IN THIS VERDICT CLEARS A THRESHOLD THE PRIOR ONE MISSED, so
it is the figure most likely to be quoted onward and the one that most needs
deriving rather than reading. +0.635 against a ~0.544 critical value at
n=14 is a claim about a specific statistic over a specific set of cells, and
if either moves the claim has to move with it.
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
NEW = ROOT / "data/results/structure_scorer_symmetric.json"
PRIOR = ROOT / "data/results/structure_scorer.json"
CRITICAL_RHO_N14 = 0.544

pytestmark = pytest.mark.skipif(
    not (NEW.exists() and PRIOR.exists()),
    reason="both structure-scorer runs are needed")


@pytest.fixture(scope="module")
def doc() -> str:
    return " ".join(PREREG.read_text().split())


@pytest.fixture(scope="module")
def runs() -> tuple[dict, dict]:
    return json.loads(NEW.read_text()), json.loads(PRIOR.read_text())


def cells(new: dict) -> dict:
    return {v: c for v, c in new["primary_chembl"].items()
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


def test_p2_the_invariant_holds_in_every_cell(doc, runs):
    """THE PREDICTION EVERYTHING ELSE RESTS ON. If the handicap is not k for
    both classes, the control is not symmetric and no AUC in the verdict
    means what it says."""
    new, _prior = runs
    k = new["k_primary"]
    for v, c in cells(new).items():
        n = c["structure_max"]["counts"]
        assert n["excluded_knn_pos"] == k * n["candidates_pos"], (
            v, "positives", n["excluded_knn_pos"], k * n["candidates_pos"])
        assert n["excluded_knn_neg"] == k * n["candidates_neg"], (
            v, "negatives", n["excluded_knn_neg"], k * n["candidates_neg"])
        assert n["knn_short"] == 0, (v, "a candidate ran short of references")
    assert "**HELD in all 14 cells" in doc


def test_p1s_figure_and_its_threshold(doc, runs):
    """P1's threshold came from the prior run; its outcome from this one."""
    new, prior = runs
    got = cells(new)["NCBITaxon:2697049"]["structure_max"]["auc"]
    was = prior["primary_chembl"]["NCBITaxon:2697049"]["structure_max"]["auc"]
    threshold = 0.5 + (was - 0.5) / 2
    assert got >= threshold, (got, threshold)
    assert f"{got:.4f}" in doc and f"{threshold:.4f}" in doc
    # AND THE DIRECTION, which the verdict concedes was unanticipated
    assert got > was, (
        "SARS-CoV-2 no longer rose under the symmetric control; the "
        "verdict's concession about being satisfied from the wrong side is "
        "stale")
    assert "satisfied from the wrong side" in doc


def test_p3s_four_movements_and_the_both_column(doc, runs):
    """Each below-chance cell's move, and the figure that puts it back."""
    import importlib
    D = importlib.import_module("diagnose_structure_scorer_design")
    names = D.virus_names()
    new, prior = runs
    k = new["k_primary"]
    moved = 0
    for v, c in cells(new).items():
        was = prior["primary_chembl"][v]["structure_max"]["auc"]
        if was >= 0.5:
            continue
        now = c["structure_max"]["auc"]
        both = c[f"structure_max_k{k}_both"]["auc"]
        assert f"{was:.4f} -> {now:.4f}" in doc, (
            f"{names.get(v, v)} moved {was:.4f} -> {now:.4f} and the verdict "
            f"does not carry that pair")
        assert f"{now - was:+.4f}" in doc
        # BOUND TO THE VIRUS NAME. `f"{both:.4f}" in doc` passed when
        # "JEV returns to 0.1608" was falsified, because the table's own
        # `both` column still carried 0.1608 -- the second time in this
        # file that an unanchored substring check proved nothing.
        name = names.get(v, v)
        assert (f"{name} returns to {both:.4f}" in doc
                or f"{name} to {both:.4f}" in doc), (
            f"{name}'s both-controls figure is {both:.4f} -- the verdict's "
            f"evidence that the publication exclusion caused the "
            f"anti-ranking -- and the prose attributes a different number "
            f"to it")
        moved += 1
    assert moved == 4, f"{moved} cells were below chance before, not four"


def test_the_asymmetry_correlation_clears_its_critical_value(doc, runs):
    """THE CLAIM MOST LIKELY TO BE QUOTED ONWARD. Prior asymmetry against
    the CHANGE in AUC -- the prior verdict measured it against the LEVEL and
    got -0.490, which did not clear."""
    new, prior = runs
    ratios, deltas = [], []
    for v, c in cells(new).items():
        n = prior["primary_chembl"][v]["structure_max"]["counts"]
        ep = n["excluded_publication_pos"] / max(n["candidates_pos"], 1)
        en = n["excluded_publication_neg"] / max(n["candidates_neg"], 1)
        if not en:
            continue
        ratios.append(math.log(ep / en))
        deltas.append(c["structure_max"]["auc"]
                      - prior["primary_chembl"][v]["structure_max"]["auc"])
    rho = _spearman(ratios, deltas)
    assert f"{rho:+.3f}" in doc, (
        f"the asymmetry-against-change correlation is {rho:+.3f} and the "
        f"verdict states something else")
    assert abs(rho) > CRITICAL_RHO_N14, (
        f"{rho:+.3f} no longer clears {CRITICAL_RHO_N14} at n=14; the "
        f"verdict claims it does and would need rewriting")


def test_p6s_leak_correlation_is_recomputed(doc, runs):
    """P6 predicted decoupling from the prior leak. Both numbers derived."""
    new, prior = runs
    aucs, leaks = [], []
    for v, c in cells(new).items():
        p = prior["primary_chembl"][v]
        aucs.append(c["structure_max"]["auc"])
        leaks.append(p["structure_max_no_pub_exclusion"]["auc"]
                     - p["structure_max"]["auc"])
    rho = _spearman(aucs, leaks)
    assert f"{rho:+.3f}" in doc, (
        f"P6's correlation is {rho:+.3f} and the verdict states otherwise")


def test_the_strength_confound_bound_is_the_runs(doc, runs):
    """The three cells where the NEW control is strictly stronger. This is
    the whole of the evidence that the rise is not just leniency, so it
    cannot be a remembered list."""
    import importlib
    D = importlib.import_module("diagnose_structure_scorer_design")
    names = D.virus_names()
    new, prior = runs
    k = new["k_primary"]
    stronger = []
    for v, c in cells(new).items():
        n = prior["primary_chembl"][v]["structure_max"]["counts"]
        ep = n["excluded_publication_pos"] / max(n["candidates_pos"], 1)
        if ep < k:
            delta = (c["structure_max"]["auc"]
                     - prior["primary_chembl"][v]["structure_max"]["auc"])
            stronger.append((names.get(v, v), ep, delta))
    assert stronger, (
        "no cell has the new control strictly stronger; the verdict's bound "
        "on the strength confound is empty")
    m = re.search(r"In (\w+) cells the publication control removed FEWER",
                  doc)
    assert m, "the verdict no longer states how many such cells there are"
    words = {"one": 1, "two": 2, "three": 3, "four": 4, "five": 5}
    assert words[m.group(1).lower()] == len(stronger), (
        f"{len(stronger)} cells have the new control stronger and the "
        f"verdict says {m.group(1)}")
    for name, ep, delta in stronger:
        assert f"publication {ep:.2f}/pos" in doc, (
            f"{name} had {ep:.2f}/pos removed and the verdict's table of the "
            f"bound does not carry it")
        assert f"{delta:+.4f}" in doc
        assert delta > 0, (
            f"{name}'s AUC fell under a stronger control, which breaks the "
            f"verdict's argument that the direction is not leniency")


def test_the_matched_cell_is_measurable_and_lower(doc, runs):
    """P5, both halves, from the run."""
    new, _prior = runs
    cell = new["secondary_ncats"]["NCBITaxon:2697049"]["structure_max"]
    assert cell.get("measurable") is True, (
        "the matched cell is unmeasurable again; P5's first half fails")
    chembl = cells(new)["NCBITaxon:2697049"]["structure_max"]["auc"]
    assert cell["auc"] < chembl, (cell["auc"], chembl)
    # BOUND TO ITS INTERVAL, not merely present. Asserting the AUC alone
    # passed when the headline "**0.6088 [0.5157, 0.7020]**" was falsified
    # to 0.6188, because the same figure also appears in the sentence
    # comparing it with the ChEMBL cell. "Somewhere in the document" is not
    # a claim about anything.
    assert (f"{cell['auc']:.4f} [{cell['auc_ci_lo_corrected']:.4f}, "
            f"{cell['auc_ci_hi_corrected']:.4f}]" in doc), (
        f"the matched cell is {cell['auc']:.4f} "
        f"[{cell['auc_ci_lo_corrected']:.4f}, "
        f"{cell['auc_ci_hi_corrected']:.4f}] and the verdict's headline "
        f"does not state that")
    assert cell["clears_corrected_alpha"] and cell["direction"] == "above", (
        "the matched cell no longer clears above chance; the verdict's "
        "headline claim depends on it")
    gap = chembl - cell["auc"]
    assert f"{gap:.2f} of AUC" in doc, (
        f"the collection-confound gap is {gap:.2f} and the verdict rounds it "
        f"differently")


def test_the_count_of_clearing_cells_is_the_runs(doc, runs):
    """Ten above, none below -- the verdict's summary line."""
    new, _prior = runs
    above = below = 0
    for c in cells(new).values():
        m = c["structure_max"]
        if not m.get("clears_corrected_alpha"):
            continue
        if m["direction"] == "above":
            above += 1
        else:
            below += 1
    words = {9: "Nine", 10: "Ten", 11: "Eleven", 12: "Twelve"}
    assert f"{words[above]} of fourteen cells clear" in doc, (
        f"{above} cells clear above chance and the verdict says otherwise")
    assert below == 0, (
        f"{below} cell(s) clear below chance; the verdict says none do")
    assert "NO CELL CLEARS IT BELOW" in doc


def test_the_prior_run_was_not_overwritten(runs):
    """THE COMPARATOR MUST SURVIVE. Both verdicts quote the publication
    control's figures, and a target that overwrote them would destroy the
    only thing this study is measured against."""
    _new, prior = runs
    assert prior.get("control") == "publication", prior.get("control")
    assert "structure_max_no_pub_exclusion" in (
        prior["primary_chembl"]["NCBITaxon:2697049"]), (
        "the prior run no longer holds its leak variant")
    # AND THE FIGURE EVERY THRESHOLD CAME FROM MUST STILL BE THERE --
    # CHECKED AGAINST THE REGISTRATION, NOT AGAINST A LITERAL. The first
    # version of this assertion hardcoded 0.7692572068116025, which is not
    # the value in the file (0.7692511756865638) and never was: I typed
    # digits from memory. A check written with a literal it invented is
    # worse than the "same literals it is checking" defect this repository
    # keeps finding, because it fails for a reason that has nothing to do
    # with the thing under test.
    auc = prior["primary_chembl"]["NCBITaxon:2697049"]["structure_max"]["auc"]
    assert f"it read {auc:.4f} under the prior control" in (
        " ".join(PREREG.read_text().split())), (
        f"the prior SARS-CoV-2 figure is {auc:.4f} and the registration "
        f"derived its P1 threshold from a different one")


def test_the_curve_table_carries_every_cell(doc, runs):
    """All three k values for all 14 cells, parsed from the table rather
    than looked for anywhere in the document -- the failure mode the
    registration's own asymmetry-table test went through twice."""
    new, _prior = runs
    # EVERY COLUMN, not just the curve: skill, both, prior and degree are
    # each a claim, and the `both` column is the verdict's evidence for its
    # mechanism.
    row = re.compile(
        r"^\s{4}(.+?)\s{2,}([\d,]+)\s+([\d,]+)\s+([\d.]+)\s+([\d.]+)\s+"
        r"([\d.]+)[\^v ]\s*\[[\d.]+, [\d.]+\]\s+"
        r"(-?[\d.]+\s+[\d.]+\s+[\d.]+\s+[\d.]+)\s*$")
    table = {}
    for line in PREREG.read_text().splitlines():
        m = row.match(line)
        if m:
            table[m.group(1).strip()] = tuple(
                float(m.group(i)) for i in (4, 5, 6)) + tuple(
                float(x) for x in m.group(7).split())
    assert table, "no curve rows parsed out of the verdict table"
    import importlib
    D = importlib.import_module("diagnose_structure_scorer_design")
    names = D.virus_names()
    ks = new["k_curve"]
    for v, c in cells(new).items():
        name = names.get(v, v)
        assert name in table, f"the curve table omits {name}"
        want = tuple(c[f"structure_max_k{k}"]["auc"] for k in ks) + (
            c["structure_max"]["auc_skill"],
            c[f"structure_max_k{new['k_primary']}_both"]["auc"],
            _prior["primary_chembl"][v]["structure_max"]["auc"],
            c["degree"]["auc"],
        )
        got = table[name]
        labels = [f"k={k}" for k in ks] + ["skill", "both", "prior", "degree"]
        assert len(got) == len(want), (name, len(got), len(want))
        for lab, w, g in zip(labels, want, got):
            assert abs(w - g) < 5e-5, (name, lab, g, w)
