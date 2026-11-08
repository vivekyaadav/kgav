"""MRR is one draw from a tie block, and nothing said so (fix-list E2).

`mrr` is `1/(rank of the first positive)` and `rank()` has no tiebreak, so
when the first positive sits inside a block of equal scores its rank -- and
therefore the whole statistic -- is decided by dict insertion order.
`hits_at_k_range` has covered Hits@k since 2026-10-14; MRR was left on one
draw, in the same sixteen `baselines_*.json` files, for the same reason the
lift was.

MEASURED, AND THE FIX LIST'S FRAMING WAS WRONG. It named "Hits@k/MRR" as one
item. On release v0.1:

    Hits@k     1 of 120 figures undetermined
    MRR       10 of  40 figures undetermined

A TIE AT THE CUT IS NOT AN UNDETERMINED FIGURE: 39 of the 120 Hits@k cells
have a block straddling k and are still determined, because the block holds
no positive. Counting ties instead of indeterminacy overstates that defect
23-fold, and it was first measured here that way.

MRR is the one that was broken, because it reads the FIRST positive's rank:
a tie at the top of the ranking moves the entire statistic, while Hits@k only
moves when a positive sits in the block at the cut. Zika `dwpc` publishes
0.3333 over a true range of [0.2500, 1.0000].
"""
from __future__ import annotations

import glob
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.baselines import hits_at_k_range, mrr, mrr_range  # noqa: E402


def r(*pairs):
    return list(pairs)


# ----------------------------------------------------------- the arithmetic
def test_no_tie_gives_a_point():
    ranked = r(("a", 3.0), ("p", 2.0), ("c", 1.0))
    lo, hi, tied, tied_pos = mrr_range(ranked, {"p"})
    assert (lo, hi) == (0.5, 0.5)
    assert tied == 0 and tied_pos == 0


def test_a_block_of_one_is_not_a_tie():
    """A unique top score reports 0 contested, not 1. The drug is simply
    first, and calling that a tie would flag every determined ranking."""
    lo, hi, tied, _ = mrr_range(r(("p", 3.0), ("a", 2.0)), {"p"})
    assert (lo, hi) == (1.0, 1.0)
    assert tied == 0


def test_the_published_figure_can_be_the_best_case():
    """THE DEFECT, at its sharpest.

    The positive happens to be inserted first in its tie block, so `mrr`
    returns 1.0 -- the best value the block admits -- while the true
    range is [0.5, 1.0]. Nothing in the published file said the figure
    was one side of a coin flip.
    """
    ranked = r(("p", 3.0), ("a", 3.0), ("c", 1.0))
    assert mrr(ranked, {"p"}) == 1.0
    lo, hi, tied, tied_pos = mrr_range(ranked, {"p"})
    assert (lo, hi) == (0.5, 1.0)
    assert tied == 2 and tied_pos == 1


def test_the_published_figure_can_be_the_worst_case():
    """And the other way round, so the fix is not a one-sided correction."""
    ranked = r(("a", 3.0), ("b", 3.0), ("p", 3.0), ("c", 1.0))
    assert mrr(ranked, {"p"}) == pytest.approx(1 / 3)
    lo, hi, _t, _p = mrr_range(ranked, {"p"})
    assert (lo, hi) == (pytest.approx(1 / 3), 1.0)


def test_several_positives_in_the_block_narrow_the_range():
    """With two positives among three tied, the worst case puts the one
    non-positive first, so the first positive is at rank 2 rather than 3."""
    ranked = r(("a", 3.0), ("p", 3.0), ("q", 3.0), ("c", 1.0))
    lo, hi, _t, tied_pos = mrr_range(ranked, {"p", "q"})
    assert (lo, hi) == (0.5, 1.0)
    assert tied_pos == 2


def test_the_block_is_found_even_when_it_starts_before_the_positive():
    """The run of equal scores can begin ahead of the first positive --
    those earlier entries are non-positives with the same score. A
    version that scanned only forwards from the positive reported the
    block as smaller and the range as narrower than it is."""
    ranked = r(("a", 5.0), ("b", 5.0), ("c", 5.0), ("p", 5.0), ("z", 1.0))
    lo, hi, tied, _ = mrr_range(ranked, {"p"})
    assert tied == 4
    assert hi == 1.0, "a positive placed first in the block ranks 1st"
    assert lo == 0.25


def test_no_positives_is_zero_not_an_error():
    assert mrr_range(r(("a", 3.0), ("b", 2.0)), set()) == (0.0, 0.0, 0, 0)


def test_an_empty_ranking_is_zero():
    assert mrr_range([], {"p"}) == (0.0, 0.0, 0, 0)


@pytest.mark.parametrize("ranked,positives", [
    (r(("p", 3.0), ("a", 3.0), ("c", 1.0)), {"p"}),
    (r(("a", 3.0), ("b", 3.0), ("p", 3.0)), {"p"}),
    (r(("a", 3.0), ("p", 2.0), ("c", 1.0)), {"p"}),
    (r(("p", 1.0),), {"p"}),
    (r(("a", 0.0), ("p", 0.0)), {"p"}),
])
def test_the_point_always_lies_inside_the_range(ranked, positives):
    """The property that makes the range trustworthy. A bound computed
    the wrong way round -- and a smaller rank is a LARGER reciprocal, so
    that inversion is easy -- would show up here rather than in a
    results file."""
    lo, hi, _t, _p = mrr_range(ranked, positives)
    assert lo <= hi, "the bounds are inverted"
    assert lo <= mrr(ranked, positives) <= hi


def test_the_hits_range_annotation_matches_what_it_returns():
    """`hits_at_k_range` was annotated tuple[int, int, int, int] and
    returned five values; its own docstring named all five. The
    annotation went stale when `expected` was added."""
    out = hits_at_k_range(r(("p", 3.0), ("a", 2.0)), {"p"}, 1)
    assert len(out) == 5
    import typing
    hints = typing.get_type_hints(hits_at_k_range)
    assert len(typing.get_args(hints["return"])) == 5


# --------------------------------------------------------- the shipped files
def baseline_files() -> list[Path]:
    return sorted(Path(p) for p in
                  glob.glob(str(ROOT / "data/results/baselines_*.json")))


def test_every_published_figure_carries_its_range():
    """The flag E2 asked for, on every cell rather than the ones that
    needed it -- a range present only where it is wide would let a
    reader infer determinacy from its absence, which is the opposite of
    the guarantee."""
    files = baseline_files()
    if not files:
        pytest.skip("the baselines files are not built in this checkout")
    for f in files:
        m = json.loads(f.read_text()).get("metrics") or {}
        assert m, f
        for scorer, cells in m.items():
            for key in ("hits@10_range", "hits@50_range",
                        "hits@100_range", "mrr_range"):
                assert key in cells, f"{f.name} {scorer} has no {key}"
            for k in (10, 50, 100):
                rng = cells[f"hits@{k}_range"]
                assert rng["lo"] <= cells[f"hits@{k}"] <= rng["hi"], (
                    f"{f.name} {scorer} Hits@{k}={cells[f'hits@{k}']} is "
                    f"outside its own range [{rng['lo']}, {rng['hi']}]")
                assert rng["determined"] == (rng["lo"] == rng["hi"])
            mr = cells["mrr_range"]
            assert mr["lo"] <= cells["mrr"] <= mr["hi"], (
                f"{f.name} {scorer} MRR={cells['mrr']} is outside its own "
                f"range [{mr['lo']}, {mr['hi']}]")


def test_mrr_is_the_metric_this_item_was_actually_about():
    """The measured asymmetry, asserted so the docstrings above cannot
    drift from it. If a release ever makes Hits@k the worse of the two,
    the framing in run_baselines.py and in the fix list is wrong and
    should fail here rather than stand."""
    files = baseline_files()
    if not files:
        pytest.skip("the baselines files are not built in this checkout")
    hits_bad = hits_total = mrr_bad = mrr_total = tied_no_pos = 0
    for f in files:
        for cells in (json.loads(f.read_text()).get("metrics") or {}).values():
            for k in (10, 50, 100):
                rng = cells[f"hits@{k}_range"]
                hits_total += 1
                if not rng["determined"]:
                    hits_bad += 1
                elif rng["tied_at_cut"]:
                    tied_no_pos += 1
            mrr_total += 1
            if not cells["mrr_range"]["determined"]:
                mrr_bad += 1
    assert mrr_bad / mrr_total > hits_bad / hits_total, (
        f"MRR is undetermined in {mrr_bad}/{mrr_total} figures and Hits@k "
        f"in {hits_bad}/{hits_total}; the docstrings call MRR the worse of "
        f"the two")
    assert tied_no_pos > hits_bad, (
        f"only {tied_no_pos} Hits@k cells have a tie at the cut with no "
        f"positive in it, against {hits_bad} undetermined -- the "
        f"distinction the docstrings rest on is not present in this data")
