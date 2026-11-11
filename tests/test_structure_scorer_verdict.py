"""The structure-scorer verdict, re-derived from the result it describes.

WHY. docs/prereg_structure_scorer.md's verdict states six prediction
outcomes, a 14-row table and two rank correlations. Every one of those is a
function of data/results/structure_scorer.json, and this project has twice
had to go back and fix a figure that was maintained by hand in prose -- the
pathway registration's "eight times" that git put at eleven, and the
release plan's "M3 and M5 have zero edges" when M3 has 11,744.

So the verdict is checked against the run rather than read.

AND THE VERDICT RECORDS A DEFECT THIS FILE PINS. The first run reported the
matched NCATS cell as `0.5000 [0.4103, 0.5897]` -- an AUC with an interval
-- when nothing in that cell had been ranked at all. A vacuous cell that
prints like a null is this repository's recurring defect, and the refusal
that replaced it is asserted here so it cannot quietly come back.
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

PREREG = ROOT / "docs/prereg_structure_scorer.md"
RESULTS = ROOT / "data/results/structure_scorer.json"

pytestmark = pytest.mark.skipif(
    not RESULTS.exists(), reason="the structure scorer has not been run")


@pytest.fixture(scope="module")
def doc() -> str:
    return " ".join(PREREG.read_text().split())


@pytest.fixture(scope="module")
def res() -> dict:
    return json.loads(RESULTS.read_text())


def readable(res: dict) -> dict:
    return {v: c for v, c in res["primary_chembl"].items()
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
    num = sum((p - ma) * (q - mb) for p, q in zip(a, b))
    den = math.sqrt(sum((p - ma) ** 2 for p in a)
                    * sum((q - mb) ** 2 for q in b))
    return num / den


def test_p3_held_the_floor_is_exactly_one_half_everywhere(doc, res):
    """THE PREDICTION THE WHOLE STUDY RESTS ON. If any labelled compound
    failed to fingerprint it takes the imputed floor, coverage drops below
    1 and the floor is no longer 0.5 -- at which point the AUCs are partly
    coverage again, exactly like the metapath channels."""
    for v, c in readable(res).items():
        m = c["structure_max"]
        assert m["coverage_pos"] == 1.0 and m["coverage_neg"] == 1.0, (
            v, m["coverage_pos"], m["coverage_neg"])
        assert m["auc_floor"] == 0.5, (v, m["auc_floor"])
        assert m["auc_ceiling"] == 1.0, (v, m["auc_ceiling"])
    fp = res["fingerprint_stats"]
    assert fp.get("unparseable", 0) == 0 and fp.get("no_smiles", 0) == 0, fp
    assert f"{fp['fingerprinted']:,} of {fp['fingerprinted']:,}" in doc, (
        f"{fp['fingerprinted']:,} compounds fingerprinted and the verdict "
        f"does not quote that figure")
    assert "**HELD exactly, all 14**" in doc


def test_p1_held_and_the_verdict_quotes_the_right_count(doc, res):
    """The leak direction, counted from the run."""
    lower = 0
    for c in readable(res).values():
        primary = c["structure_max"]["auc"]
        no_excl = c["structure_max_no_pub_exclusion"]["auc"]
        if primary < no_excl:
            lower += 1
    total = len(readable(res))
    assert f"**HELD, {lower} of {total}**" in doc, (
        f"the publication-disjoint variant reads lower in {lower} of "
        f"{total} cells and the verdict says otherwise")
    assert lower >= 8, (
        f"P1 registered >= 8 and the run gives {lower}; the verdict's HELD "
        f"is wrong")


def test_the_table_carries_every_cells_auc(doc, res):
    """Fourteen AUCs, each from the file."""
    for v, c in readable(res).items():
        auc = c["structure_max"]["auc"]
        assert f"{auc:.4f}" in doc, (
            f"{v} reads {auc:.4f} and the verdict's table does not carry it")


def test_the_leak_magnitudes_in_the_prose_are_the_runs(doc, res):
    """The four magnitudes the verdict names, and the one null."""
    deltas = {}
    for v, c in readable(res).items():
        deltas[v] = (c["structure_max_no_pub_exclusion"]["auc"]
                     - c["structure_max"]["auc"])
    for v, d in sorted(deltas.items(), key=lambda x: -abs(x[1]))[:4]:
        assert f"{d:+.4f}" in doc, (
            f"{v}'s leak is {d:+.4f} and the verdict does not quote it")


def test_the_two_correlations_are_recomputed(doc, res):
    """Both rank correlations, and the n=14 critical value they are read
    against. The asymmetry one does NOT clear it, and the verdict has to
    keep saying so."""
    cells = readable(res)
    aucs, leaks, ratios = [], [], []
    for c in cells.values():
        m = c["structure_max"]
        k = m["counts"]
        ep = k["excluded_publication_pos"] / max(k["candidates_pos"], 1)
        en = k["excluded_publication_neg"] / max(k["candidates_neg"], 1)
        if not en:
            continue
        aucs.append(m["auc"])
        leaks.append(c["structure_max_no_pub_exclusion"]["auc"] - m["auc"])
        ratios.append(math.log(ep / en))
    rho_leak = _spearman(aucs, leaks)
    rho_asym = _spearman(aucs, ratios)
    assert f"{rho_leak:.3f}" in doc, (
        f"Spearman(AUC, leak) is {rho_leak:.3f} and the verdict says "
        f"something else")
    assert f"{rho_asym:.3f}" in doc, (
        f"Spearman(AUC, log exclusion ratio) is {rho_asym:.3f} and the "
        f"verdict says something else")
    # and the verdict must not promote the weaker one
    assert abs(rho_asym) < 0.544, (
        "the asymmetry correlation now clears the n=14 critical value; the "
        "verdict calls it a hypothesis and would need rewriting")
    assert abs(rho_leak) > 0.544, (
        "the leak correlation no longer clears the critical value the "
        "verdict claims for it")


def test_the_matched_cell_refuses_rather_than_reporting(doc, res):
    """THE DEFECT, PINNED. A cell where nothing was scored must not carry
    an AUC."""
    cell = res["secondary_ncats"]["NCBITaxon:2697049"]["structure_max"]
    assert cell.get("measurable") is False, (
        "the matched NCATS cell now reports as measurable; if that is real "
        "the verdict's P5 section is wrong, and if it is not, the refusal "
        "has been lost")
    assert "auc" not in cell, (
        "the unmeasurable cell carries an AUC again -- 0.5 by arithmetic "
        "over one tie block, which is what the first run printed with a "
        "confidence interval beside it")
    assert cell["counts"]["empty_reference"] == cell["counts"]["candidates"], (
        "not every candidate lost its reference set, so the cell is not "
        "vacuous and should be measured")
    assert "**CANNOT BE EVALUATED**" in doc
    assert f"{cell['counts']['excluded_publication']:,} times" in doc, (
        "the verdict does not quote how many exclusions fired")


def test_p6_failure_count_is_the_runs(doc, res):
    """Degree clearing where the scorer does not -- the count and the
    named cells."""
    import importlib
    sys.path.insert(0, str(ROOT / "scripts"))
    D = importlib.import_module("diagnose_structure_scorer_design")
    names = D.virus_names()
    violations = []
    for v, c in readable(res).items():
        s = c["structure_max"].get("clears_corrected_alpha")
        g = c["degree"].get("clears_corrected_alpha")
        if g and not s:
            violations.append(names.get(v, v))
    m = re.search(r"\*\*FAILED, (\d+) cells\*\*", doc)
    assert m, "the verdict no longer states P6's failure count"
    assert int(m.group(1)) == len(violations), (
        f"degree clears where the scorer does not in {len(violations)} "
        f"cell(s): {sorted(violations)}; the verdict says {m.group(1)}")
    for n in violations:
        assert n in doc, f"{n} is a P6 violation and the verdict omits it"


def test_the_below_chance_cells_are_named_as_unresolved(doc, res):
    """A below-chance AUC under a control with a known asymmetry is not a
    finding about chemistry, and the verdict has to keep saying so."""
    below = [v for v, c in readable(res).items()
             if c["structure_max"]["auc"] < 0.5]
    assert below, "no cell reads below chance; the verdict's caveat is stale"
    assert "are not a finding about chemistry" in doc
    # THE COUNT IN THE PROSE, NOT JUST IN THIS TEST. The first version
    # asserted len(below) == 4 against a literal and never looked at the
    # document, so rewriting "four below-chance cells" as "three" passed.
    # Caught by perturbing the verdict six ways: this was the one that got
    # through.
    words = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five",
             6: "six", 7: "seven"}
    word = words.get(len(below))
    assert word, f"{len(below)} below-chance cells; extend the word table"
    assert f"The {word} below-chance cells" in doc, (
        f"{len(below)} cells read below chance and the verdict's prose does "
        f"not say '{word}'")
    assert f"the {word} where it handicaps positives" in doc, (
        f"the asymmetry paragraph does not say '{word}' either")


def test_the_asymmetry_table_is_whole_and_its_claim_is_true(doc, res):
    """THE SAME CHERRY-PICK, IN A SECOND PLACE, CAUGHT A DAY LATE.

    docs/prereg_structure_scorer_symmetric.md's copy of this table was
    corrected on 2026-11-11 after its test was rewritten three times: the
    first draft required only the rows supporting the hypothesis, the second
    was satisfied by the prose describing the defect, and only the third
    parsed the table. THIS copy kept the cherry-pick -- four rows, an
    elision, two rows -- and the paragraph under it asserted two things
    that the full fourteen rows contradict: that two cells handicap
    negatives more (three do) and that the four cells at 5x or more against
    positives are the four lowest (five are at 5x or more, two of them above
    chance, and the third-lowest cell is at 4.98x).

    So this parses every row, checks every number against the run, checks
    the counterexample marks against the rule the table states, and checks
    the surviving claim about the extremes.
    """
    import importlib
    D = importlib.import_module("diagnose_structure_scorer_design")
    names = D.virus_names()
    row = re.compile(
        r"^\s{4}(.+?)\s{2,}([\d.]+) / ([\d.]+)\s+=\s+([\d.]+)x"
        r"\s+AUC ([\d.]+)(\s+<- counterexample)?\s*$")
    table = {}
    for line in PREREG.read_text().splitlines():
        m = row.match(line)
        if m:
            table[m.group(1).strip()] = (
                float(m.group(2)), float(m.group(3)), float(m.group(4)),
                float(m.group(5)), bool(m.group(6)))
    assert table, "no asymmetry table rows parsed out of the verdict"

    expected = {}
    for v, c in readable(res).items():
        k = c["structure_max"]["counts"]
        ep = k["excluded_publication_pos"] / max(k["candidates_pos"], 1)
        en = k["excluded_publication_neg"] / max(k["candidates_neg"], 1)
        if not en:
            continue
        expected[names.get(v, v)] = (ep, en, ep / en,
                                     c["structure_max"]["auc"])
    assert sorted(table) == sorted(expected), (
        f"the table carries {sorted(table)} against cells "
        f"{sorted(expected)}; every readable cell must be shown, with no "
        f"elision")

    by_auc = sorted(expected, key=lambda n: expected[n][3])
    top2 = set(by_auc[-2:])
    for name, (ep, en, ratio, auc) in expected.items():
        t_ep, t_en, t_ratio, t_auc, marked = table[name]
        assert abs(t_ep - ep) < 0.005, (name, "pos rate", t_ep, ep)
        assert abs(t_en - en) < 0.005, (name, "neg rate", t_en, en)
        assert abs(t_ratio - ratio) < 0.05, (name, "ratio", t_ratio, ratio)
        assert abs(t_auc - auc) < 5e-5, (name, "auc", t_auc, auc)
        # the rule the table's own legend states
        is_counter = ((ratio >= 5.0 and auc > 0.5)
                      or (ratio < 1.0 and name not in top2))
        assert marked == is_counter, (
            f"{name} ({ratio:.2f}x, AUC {auc:.4f}) is "
            f"{'a' if is_counter else 'not a'} counterexample under the "
            f"legend and the table marks it {'' if marked else 'not '}so")

    # THE SURVIVING CLAIM, which is about the extremes only.
    for name in top2:
        assert expected[name][2] < 1.0, (
            f"{name} is one of the two highest cells and is no longer "
            f"handicapped against negatives; the verdict's surviving claim "
            f"about the extremes is wrong")
    for name in by_auc[:2]:
        assert expected[name][2] >= 5.0, (
            f"{name} is one of the two lowest cells and is no longer "
            f"handicapped 5x or more against positives")
    # AND THE CORRECTION PARAGRAPH'S COUNTS ARE DERIVED, not typed here.
    # A literal "three" on both sides would be the "check written with the
    # same literals it is checking" defect, which is how the sentence this
    # test replaces survived in the first place.
    words = {1: "One", 2: "Two", 3: "Three", 4: "Four", 5: "Five",
             6: "Six", 7: "Seven"}
    neg_heavy = [n for n in expected if expected[n][2] < 1.0]
    pos5 = [n for n in expected if expected[n][2] >= 5.0]
    assert f"{words[len(neg_heavy)]} cells handicap negatives more" in doc, (
        f"{len(neg_heavy)} cells handicap negatives more and the verdict's "
        f"correction paragraph says otherwise")
    assert (f"{words[len(pos5)]} cells are handicapped 5x or more against "
            f"positives" in doc), (
        f"{len(pos5)} cells are handicapped 5x or more against positives "
        f"and the verdict's correction paragraph says otherwise")
