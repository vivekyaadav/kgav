"""Every figure in doc2/plan-v2.md, re-derived.

THAT DOCUMENT IS A CRITIQUE OF V1 WRITTEN FROM V1'S OWN ARTIFACTS, which
makes it the one document in the repository with the most to gain from a
figure nobody checks: a plan that overstates the defect it proposes to fix
buys work that is not needed, and a plan that understates it leaves a
published claim standing. It is hand-maintained -- what it records is a
decision and a priority, which no artifact can derive -- so every quantity
in it is held here instead.

The document's own item C8 exists because `09-limits.md` quotes a threshold
figure that is typed in a docstring and derived by nothing, and had drifted.
This file is what stops plan-v2.md from repeating it.

Two of the assertions are ABOUT ABSENCE -- that no paired AUC test exists
(A1) and that the metapath artifacts carry no corrected interval (A3) --
because both are premises of the plan rather than consequences of it. When
v2 lands they will fail, and the failure is the signal that the item is
done and the plan's past tense is wrong.
"""
from __future__ import annotations

import json
import math
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

DOC = ROOT / "doc2" / "plan-v2.md"
RESULTS = ROOT / "data" / "results"
RELEASE = ROOT / "data" / "releases" / "v0.1"

needs_artifacts = pytest.mark.skipif(
    not (RELEASE / "MANIFEST.json").exists(),
    reason="no assembled release")

NAME = {
    "NCBITaxon:11053": "DENV-1", "NCBITaxon:11060": "DENV-2",
    "NCBITaxon:11069": "DENV-3", "NCBITaxon:11070": "DENV-4",
    "NCBITaxon:11072": "JEV", "NCBITaxon:11082": "West Nile",
    "NCBITaxon:11084": "TBEV", "NCBITaxon:11089": "Yellow fever",
    "NCBITaxon:11137": "HCoV-229E", "NCBITaxon:2697049": "SARS-CoV-2",
    "NCBITaxon:64320": "Zika", "NCBITaxon:31631": "HCoV-OC43",
    "NCBITaxon:694009": "SARS-CoV", "NCBITaxon:1335626": "MERS-CoV",
}
LICENSED = ["NCBITaxon:11053", "NCBITaxon:11060", "NCBITaxon:11069",
            "NCBITaxon:11137", "NCBITaxon:11072", "NCBITaxon:2697049",
            "NCBITaxon:11084", "NCBITaxon:11089"]


def doc() -> str:
    return DOC.read_text()


def flat() -> str:
    """The document with runs of whitespace collapsed.

    Every figure below is asserted against this rather than against the
    raw text, because a figure must survive the line wrapping being
    changed -- and reflowing a paragraph is the commonest edit this
    document will get.
    """
    return " ".join(doc().split())


def says(needle: str) -> None:
    assert needle in flat(), f"plan-v2.md no longer says {needle!r}"


def results(name: str) -> dict:
    return json.loads((RESULTS / f"{name}.json").read_text())


def hanley_mcneil_se(auc: float, n_pos: int, n_neg: int) -> float:
    """The SE the project's own auc_interval uses, at a given AUC.

    Imported rather than reimplemented would be better; it is inlined
    because `labels.auc_interval` returns bounds and the plan quotes the
    SE itself. A test asserts the two agree.
    """
    q1 = auc / (2 - auc)
    q2 = 2 * auc * auc / (1 + auc)
    return math.sqrt((auc * (1 - auc) + (n_pos - 1) * (q1 - auc ** 2)
                      + (n_neg - 1) * (q2 - auc ** 2)) / (n_pos * n_neg))


def test_the_se_this_file_computes_is_the_projects_own():
    """If these disagree, every minimum detectable effect below is wrong."""
    from kgav.labels import auc_interval
    lo, hi = auc_interval(0.7, 100, 400)
    se = hanley_mcneil_se(0.7, 100, 400)
    assert abs((hi - lo) / 2 - 1.96 * se) < 1e-9, (
        f"auc_interval's half-width {(hi - lo) / 2} is not 1.96 * "
        f"{se} -- this file's SE is not the project's")


# ---------------------------------------------------------------- A1

def test_no_paired_auc_test_exists():
    """A1's premise. WHEN V2 LANDS THIS FAILS, and that is the point."""
    hit = subprocess.run(
        ["grep", "-rniE", r"delong|paired[_ ]bootstrap|paired_auc",
         "src", "scripts"], cwd=ROOT, capture_output=True, text=True)
    lines = [L for L in hit.stdout.splitlines() if L.strip()]
    assert not lines, (
        "plan-v2.md item A1 says no paired AUC test exists and grep now "
        "finds one:\n" + "\n".join(lines[:10])
        + "\nIf A1 is done, rewrite it in the past tense and delete this.")
    says("returns nothing")


@needs_artifacts
def test_a1_quotes_the_assay_matched_population():
    d = results("hard_negatives")
    cell = d["cross_sectional"]["SARS-CoV-2"]
    comb, deg = cell["COMBINED"], cell["degree"]
    assert d["publication_disjoint"] is True
    assert d["label_source"] == ["infores:ncats-opendata"]
    says(f"COMBINED at {comb['auc']:.4f}")
    says(f"degree at {deg['auc']:.4f}")
    says(f"the same {comb['n_pos']} actives")
    says(f"the same {comb['n_neg']:,} inactives")
    assert comb["n_pos"] == deg["n_pos"] and comb["n_neg"] == deg["n_neg"], (
        "A1's argument is that the two scorers share a population; on this "
        "artifact they do not")


# ---------------------------------------------------------------- A2

@needs_artifacts
def test_a2s_minimum_detectable_effect_and_the_channels_inside_it():
    d = results("hard_negatives")
    cell = d["cross_sectional"]["SARS-CoV-2"]
    n_pos = cell["COMBINED"]["n_pos"]
    n_neg = cell["COMBINED"]["n_neg"]
    se = hanley_mcneil_se(0.5, n_pos, n_neg)
    says(f"SE {se:.4f}")
    says(f"minimum detectable AUC at uncorrected 95% is **{0.5 + 1.96 * se:.4f}**")

    z = results("structure_scorer_strength_matched")["z_corrected"]
    says(f"the bar rises to {0.5 + z * se:.4f}")

    channels = {m: cell[m]["auc"] for m in ("M1", "M2", "M3", "M4", "M5",
                                            "M6", "M7", "M8")}
    floor = 0.5 + 1.96 * se
    above = {m: a for m, a in channels.items() if a >= floor}
    assert not above, (
        f"A2's claim is that ALL EIGHT channels read below the detection "
        f"floor {floor:.4f}; these now read above it: {above}")
    says(f"The highest of them is {max(channels.values()):.4f}")
    says(f"the lowest is {min(channels.values()):.4f}")

    assert 0.5 + z * se > max(cell["COMBINED"]["auc"], cell["degree"]["auc"]), (
        "A2 says neither COMBINED nor degree would clear the corrected bar; "
        "one of them now would")

    comb_floor = cell["COMBINED"]["auc_floor"]
    share = (comb_floor - 0.5) / (cell["COMBINED"]["auc"] - 0.5)
    says(f"{round(share * 100)}% of the")


# ---------------------------------------------------------------- A3

@needs_artifacts
def test_a3s_two_standards_are_still_two():
    """A3's premise, on both sides of the asymmetry it describes."""
    for name in ("hard_negatives", "hard_negatives_row1",
                 "hard_negatives_row2", "hard_negatives_row3",
                 "hard_negatives_row4"):
        d = results(name)
        for virus, cell in (d.get("cross_sectional") or {}).items():
            for scorer, m in cell.items():
                assert "clears_corrected_alpha" not in m, (
                    f"{name} / {virus} / {scorer} now carries a corrected "
                    f"interval -- plan-v2.md item A3 is done; rewrite it")
    sm = results("structure_scorer_strength_matched")["primary_chembl"]
    for virus, cell in sm.items():
        s = cell.get("structure_max")
        if s and s.get("evaluable"):
            assert "clears_corrected_alpha" in s, (
                f"{virus} carries no corrected interval, so A3's "
                f"asymmetry is not the one described")
    says(f"z {results('structure_scorer_strength_matched')['z_corrected']}")


# ---------------------------------------------------------------- B4

@needs_artifacts
def test_b4s_overlap_table_is_the_label_sets():
    from kgav.labels import build_labels
    lab = build_labels(RELEASE)
    rows = []
    for i, a in enumerate(LICENSED):
        for b in LICENSED[i + 1:]:
            A = lab.positives.get(a, set())
            B = lab.positives.get(b, set())
            if not A or not B:
                continue
            inter = len(A & B)
            if not inter:
                continue
            rows.append((inter / len(A | B), inter / min(len(A), len(B)),
                         NAME[a], NAME[b], inter, min(len(A), len(B))))
    rows.sort(reverse=True)
    top = rows[:5]
    assert [f"{r[2]} / {r[3]}" for r in top] == [
        "DENV-1 / DENV-3", "DENV-2 / DENV-3", "DENV-1 / DENV-2",
        "JEV / Yellow fever", "HCoV-229E / SARS-CoV-2"], (
        f"B4's five pairs are no longer the five largest: "
        f"{[f'{r[2]} / {r[3]}' for r in top]}")
    for jac, share, a, b, *_ in top:
        says(f"| {a} / {b} |")
        says(f"{jac:.3f}")
        says(f"{share * 100:.1f}%")
    jev = lab.positives["NCBITaxon:11072"]
    yf = lab.positives["NCBITaxon:11089"]
    says(f"{len(jev & yf)} of JEV's {len(jev)} positives")


# ---------------------------------------------------------------- B5

@needs_artifacts
def test_b5s_margins_and_the_thinnest_licence():
    sm = results("structure_scorer_strength_matched")
    pc = sm["primary_chembl"]
    for v in ("NCBITaxon:11053", "NCBITaxon:11072", "NCBITaxon:11084"):
        c, s = pc[v], pc[v]["structure_max"]
        assert s["clears_corrected_alpha"] is True, (
            f"{NAME[v]} no longer clears, so B5's table is not about "
            f"licensed cells")
        margin = s["auc_ci_lo_corrected"] - 0.5
        says(f"| {NAME[v]} | {c['n_positives']} | {c['n_negatives']} | "
             f"{s['auc']:.4f} | {s['auc_ci_lo_corrected']:.4f} | "
             f"**+{margin:.4f}**" if v == "NCBITaxon:11053" else
             f"| {NAME[v]} | {c['n_positives']} | {c['n_negatives']} | "
             f"{s['auc']:.4f} | {s['auc_ci_lo_corrected']:.4f} | "
             f"+{margin:.4f} |")
    d1 = pc["NCBITaxon:11053"]["structure_max"]
    assert d1["auc_ci_lo_corrected"] - 0.5 < 0.01, (
        "B5 calls DENV-1's margin three thousandths; it has moved above a "
        "hundredth and the sentence is no longer true")

    jev = pc["NCBITaxon:11072"]
    k = jev["k_cell"]
    says(f"{jev['n_positives']} positives, of which k_cell {k}")
    says(f"leaving {jev['n_positives'] - k}")
    se = hanley_mcneil_se(0.5, jev["n_positives"], jev["n_negatives"])
    says(f"alpha is {0.5 + sm['z_corrected'] * se:.4f}")

    from kgav.calibration import MIN_REACHED_POSITIVES
    says(f"`MIN_REACHED_POSITIVES = {MIN_REACHED_POSITIVES}`")


# ---------------------------------------------------------------- B6

@needs_artifacts
def test_b6s_three_control_table():
    runs = {k: results(n)["primary_chembl"] for k, n in (
        ("pub", "structure_scorer"),
        ("sym", "structure_scorer_symmetric"),
        ("str", "structure_scorer_strength_matched"))}
    assert results("structure_scorer")["control"] == "publication"
    assert results("structure_scorer_symmetric")["control"] == "knn"
    assert results("structure_scorer_strength_matched")["control"] == "strength"

    rows = []
    for v in LICENSED:
        a = [runs[k][v]["structure_max"]["auc"] for k in ("pub", "sym", "str")]
        rows.append((max(a) - min(a), NAME[v], a))
    rows.sort(reverse=True)
    assert [r[1] for r in rows] == [
        "JEV", "Yellow fever", "DENV-3", "DENV-2", "TBEV", "DENV-1",
        "HCoV-229E", "SARS-CoV-2"], (
        f"B6's table is no longer ordered by spread: {[r[1] for r in rows]}")
    for spread, name, (p, s, m) in rows:
        says(f"| {name} |")
        for x in (p, s, m):
            says(f"{x:.4f}")
        says(f"{spread:.4f}")

    worst = rows[0]
    assert worst[1] == "JEV", worst
    says(f"a {worst[0]:.4f} swing on a licensed cell")

    prereg = (ROOT / "doc2" / "registrations"
              / "prereg_structure_scorer_symmetric.md").read_text()
    asym = re.search(r"DENV-2\s+[\d.]+\s*/\s*[\d.]+\s*=\s*([\d.]+)x", prereg)
    assert asym, "the symmetric registration no longer states DENV-2's ratio"
    says(f"stripping {asym.group(1)}x more references")


# ---------------------------------------------------------------- B7

@needs_artifacts
def test_b7s_temporal_ladder_figures():
    d = results("temporal_structure_scorer")
    says("`diagnose_prospective_analogy.py`")
    text = flat()
    assert "0.8025 to 0.3574" in text, (
        "B7 quotes DENV-1's fall across the temporal ladder")
    assert "2 of 6" in text or "2 of 6 " in text
    assert d, "temporal_structure_scorer.json is empty"


def test_b7_names_a_script_that_does_not_exist_yet():
    """B7 proposes a script. If it exists, the item is done."""
    assert not (ROOT / "scripts"
                / "diagnose_cross_sectional_analogy.py").exists(), (
        "scripts/diagnose_cross_sectional_analogy.py exists -- plan-v2.md "
        "item B7 is done; rewrite it in the past tense")
    assert (ROOT / "scripts" / "diagnose_prospective_analogy.py").exists(), (
        "B7's argument is that the temporal ladder exists and the "
        "cross-sectional one does not; the temporal one is gone")


# ---------------------------------------------------------------- C8

@needs_artifacts
def test_c8s_threshold_flip_is_derived_here_and_nowhere_else():
    from kgav.labels import build_labels
    at10 = build_labels(RELEASE)
    at1 = build_labels(RELEASE, active_below_nm=1000.0,
                       inactive_above_nm=1000.0)
    pairs = sum(len(v) for v in at10.positives.values())
    kept = sum(len(at10.positives[v] & at1.positives.get(v, set()))
               for v in at10.positives)
    flipped = pairs - kept
    says(f"{pairs:,} compound–virus pairs are positive at 10 µM")
    says(f"{kept:,} of them are still positive at 1 µM")
    says(f"{flipped:,} flip")
    says(f"**{100 * flipped / pairs:.1f}%**")

    A = set().union(*at10.positives.values())
    B = set().union(*at1.positives.values())
    says(f"the figure is {100 * len(A - B) / len(A):.1f}%")

    stale = re.search(r"read \*\"([\d.]+)% of positives at 10 µM", flat())
    assert stale, "C8 no longer quotes the figure 09-limits.md used to carry"
    quoted = float(stale.group(1))
    assert abs(quoted - 100 * flipped / pairs) > 0.5, (
        f"C8's point is that the typed {quoted}% was not the derived "
        f"figure. They now agree, so the drift it reports did not happen")

    # AND THE HALF C8 CALLS FIXED MUST BE FIXED. The stale figure has to
    # be gone from the docstring that held it and from the document that
    # quoted it, or C8's past tense is a claim about work nobody did.
    # THE HISTORY MAY NAME THE NUMBER; THE DOCSTRING MAY NOT ASSERT IT.
    # labels.py records what the figure was and that it drifted, which is
    # the record this project keeps on purpose. What must be gone is the
    # number stated as the quantity -- so the check is on the claim's
    # shape, not on the digits appearing anywhere in the file.
    src = " ".join((ROOT / "src" / "kgav" / "labels.py").read_text().split())
    for shape in (f"{quoted}% of positives", f"{quoted}% of the positive"):
        assert shape not in src, (
            f"src/kgav/labels.py still states {shape!r} as the quantity, "
            f"which C8 says was deleted rather than corrected")
    limits = " ".join((ROOT / "doc2" / "09-limits.md").read_text().split())
    assert f"{quoted}% of positives at 10 µM" not in limits, (
        f"doc2/09-limits.md still prints the stale {quoted}%")
    assert f"{100 * flipped / pairs:.1f}% of positives at 10 µM" in limits, (
        "doc2/09-limits.md does not print the derived flip fraction, so "
        "C8's claim that it is now derived is false")


def test_c8s_premise_that_no_scorer_varies_the_cut():
    """C8's open half: a CENSUS at the lower cut is not a RESULT at it.

    `write_doc2.py` now builds labels at 1 uM to derive the flip
    fraction, which is the half of C8 that is done. The item stays open
    until a scorer runs there, so this excludes the generator by name
    rather than by pattern -- a grep that quietly matched the fix would
    make the premise unfalsifiable.
    """
    hit = subprocess.run(
        ["grep", "-rn", "active_below_nm", "scripts"],
        cwd=ROOT, capture_output=True, text=True)
    varied = [L for L in hit.stdout.splitlines()
              if L.strip() and "DEFAULT" not in L
              and not L.startswith("scripts/write_doc2.py:")]
    assert not varied, (
        "a script other than the doc generator now sets active_below_nm "
        "-- plan-v2.md item C8 may be done:\n" + "\n".join(varied[:10]))
    assert subprocess.run(
        ["grep", "-rn", "active_below_nm", "scripts/write_doc2.py"],
        cwd=ROOT, capture_output=True, text=True).stdout.strip(), (
        "write_doc2.py no longer derives the flip fraction, so C8's "
        "claim that half of it is fixed is false")
    says("is never varied from its default in any script")


# ---------------------------------------------------------------- C9

def test_c9s_premise_that_no_embedding_comparator_exists():
    hit = subprocess.run(
        ["grep", "-rniE", r"transe|rotate|distmult|node2vec|graph neural",
         "src", "scripts"], cwd=ROOT, capture_output=True, text=True)
    lines = [L for L in hit.stdout.splitlines() if L.strip()]
    assert not lines, (
        "plan-v2.md item C9 says no embedding comparator exists:\n"
        + "\n".join(lines[:10]))


# ---------------------------------------------------------------- C10

@needs_artifacts
def test_c10s_two_anti_ranked_cells():
    """THE FIRST VERSION OF THIS TEST TYPED THE TWO AUCS IT HELD, which
    is the defect C8 reports, in the file that exists to prevent it.
    Both the figures and WHICH cells are anti-ranked are derived."""
    d = results("temporal_structure_scorer")
    below = {}
    for cell in d["cells"].values():
        pr = cell.get("primary") or {}
        if pr.get("direction") == "below" and pr.get(
                "clears_corrected_alpha"):
            below[cell["name"]] = pr["auc"]
    assert len(below) == 2, (
        f"C10 is about two cells clearing below chance; the artifact now "
        f"has {len(below)}: {below}")
    for name, auc in sorted(below.items()):
        says(f"{auc:.4f}")
        assert name in flat(), f"C10 does not name {name}"
    first = sorted(below)[0]
    says(f"{first} reads {below[first]:.4f}")


# ------------------------------------------------- the document itself

def test_it_tells_its_reader_it_is_hand_maintained():
    head = " ".join(doc().split()[:60])
    assert "Not generated" in head, (
        "plan-v2.md does not warn its reader that it is hand-maintained")
    assert "tests/test_v2_plan_figures.py" in doc(), (
        "plan-v2.md does not name the test that holds its figures")


def test_it_has_no_template():
    assert not (ROOT / "doc2" / "_templates" / "plan-v2.md.in").exists(), (
        "plan-v2.md has a template and is also on the exemption list")


def test_the_prose_counts_the_items_it_actually_has():
    """IT SAID NINE AND HAD TEN, in four places at once, within an hour
    of being written -- which is this repository's defect 11, `a table
    that is almost complete`, in the document whose whole subject is
    counts that do not match what produced them."""
    word = {9: "nine", 10: "ten", 11: "eleven", 12: "twelve"}
    items = re.findall(r"^### ([A-C]\d+) \u2014 ", doc(), re.M)
    n = len(items)
    assert n in word, f"{n} items; this test needs a word for it"
    says(f"This is what v2 is for: {word[n]} fixes")
    says(f"**None of the {word[n]} items needs new data.**")
    says(f"## The {word[n]} items, by what they cost")
    for wrong in set(word.values()) - {word[n]}:
        assert f"the {wrong} items" not in flat(), (
            f"plan-v2.md has {n} items and somewhere says {wrong!r}")


def test_no_literal_placeholder_survives_into_the_document():
    """doc2 forbids a `{{...}}` anywhere in the folder, because it cannot
    be told apart from a figure that never rendered. This document
    discusses placeholders and so is the one most likely to print one."""
    left = re.findall(r"\{\{[^}]*\}\}", doc())
    assert not left, (
        f"plan-v2.md prints {left}, which tests/test_doc2.py reads as an "
        f"unrendered figure across the whole folder")


def test_every_item_states_a_cost_and_a_remedy():
    """AN ITEM WITH NO COST IS A WISH, and one with no remedy is a
    complaint. `WHAT WOULD CHANGE` is NOT required here: four of the ten
    change what a reader is shown rather than what a number says, and
    the summary table carries that column for all ten -- so requiring
    the heading would buy a heading and not a thought."""
    items = re.findall(r"^### ([A-C]\d+) — (.+)$", doc(), re.M)
    assert len(items) == 10, f"expected ten items, found {items}"
    bodies = re.split(r"^### [A-C]\d+ — ", doc(), flags=re.M)[1:]
    for (code, title), body in zip(items, bodies):
        assert "**COST.**" in body or "**COST**" in body, (
            f"item {code} ({title}) states no cost")
        assert "**WHAT V2 DOES.**" in body, (
            f"item {code} ({title}) says what is wrong and not what to do")


def test_the_total_is_the_sum_of_the_costs_it_quotes():
    """A PLAN'S HEADLINE COST IS A FIGURE LIKE ANY OTHER. It said eight
    to ten over a table summing to nine and a half, which is the same
    defect as the document's item counts and is caught the same way."""
    summary = re.split(r"^## The \w+ items, by what they cost", doc(),
                       flags=re.M)[1].split("## What v2 is not")[0]
    cost = {"hours": (0.1, 0.5), "\u00bd day": (0.5, 0.5),
            "1 day": (1.0, 1.0), "3\u20134 days": (3.0, 4.0)}
    lo = hi = 0.0
    rows = re.findall(r"^\| [A-C]\d+ \| [^|]+ \| ([^|]+) \|", summary, re.M)
    assert len(rows) == len(re.findall(r"^### [A-C]\d+ \u2014 ", doc(), re.M))
    for c in rows:
        c = c.strip()
        assert c in cost, f"the summary table quotes a cost {c!r} this test "\
                          f"has no range for; add it rather than widening it"
        lo += cost[c][0]
        hi += cost[c][1]
    word = {8: "Eight", 9: "Nine", 10: "Ten", 11: "Eleven", 12: "Twelve",
            13: "Thirteen"}
    import math
    stated = re.search(r"^(\w+) to (\w+) working days", summary, re.M)
    assert stated, "the summary states no total"
    want = (word[math.floor(lo)], word[math.ceil(hi)].lower())
    assert (stated.group(1), stated.group(2)) == want, (
        f"the table sums to {lo}-{hi} days, so the total should read "
        f"{want[0]} to {want[1]}, and it reads "
        f"{stated.group(1)} to {stated.group(2)}")


def test_the_summary_table_lists_every_item_once():
    items = [c for c, _ in re.findall(r"^### ([A-C]\d+) — (.+)$", doc(), re.M)]
    summary = re.split(r"^## The \w+ items, by what they cost", doc(),
                       flags=re.M)[1]
    summary = summary.split("## What v2 is not")[0]
    table = re.findall(r"^\| ([A-C]\d+) \| ", summary, re.M)
    assert table == items, (
        f"the summary table lists {table} and the document has {items} -- "
        f"a table that is almost complete is this project's own defect")
