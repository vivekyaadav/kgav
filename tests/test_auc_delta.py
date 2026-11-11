"""The AUC regression gate, proven to fail.

WHY THIS FILE EXISTS. `scripts/compare_results_aucs.py` answers the one
question a rebuild has to answer -- did the scores go down -- and when it was
written it PASSED on the only comparison anyone had run. A gate that has only
ever passed is indistinguishable from a gate that cannot fail, which is the
defect this repository keeps finding; so every refusal it can issue is
exercised here against a planted regression.

IT ALSO PINS THE THREE DISTINCTIONS THE GATE RESTS ON, each of which was
wrong in a first draft and each of which was found by running it rather than
reading it:

  * `auc_skill` is a RATIO over a cell's own span, not a score. On a cell
    whose span is thousandths it moves by most of its range for an AUC move
    of 0.001 -- the largest "fall" in the first run was -0.9286 beside an
    AUC move of +0.0071. Gating on the mixture refuses every rebuild.
  * An ABSENT readability flag is not a pass. `hard_negatives_row1.json`
    carries no `trustworthy`, so the first classifier admitted a cell with
    THREE reached positives, whose AUC moved 0.0438 and fired the gate.
  * A COMPARATOR is not a result. `degree` is the null hypothesis and
    `structure_max_no_pub_exclusion` is the scorer with its control removed;
    both falling means the argument got stronger. The first version of that
    filter read `"no_pub_exclusion" in k` against a TUPLE, which is
    membership and not substring, so half of it silently never matched.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

_spec = importlib.util.spec_from_file_location(
    "compare_results_aucs", ROOT / "scripts" / "compare_results_aucs.py")
C = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(C)


def write(d: Path, name: str, payload: dict) -> None:
    d.mkdir(parents=True, exist_ok=True)
    (d / name).write_text(json.dumps(payload))


def cell(auc, *, lo=None, hi=None, n_pos=100, cov_pos=0.5, skill=None,
         **extra) -> dict:
    """A results cell shaped like the real ones."""
    out = {"auc": auc, "n_pos": n_pos, "n_neg": n_pos,
           "coverage_pos": cov_pos, "coverage_neg": cov_pos}
    if lo is not None:
        out["auc_ci_lo"], out["auc_ci_hi"] = lo, hi
    if skill is not None:
        out["auc_skill"] = skill
    out.update(extra)
    return out


# ------------------------------------------------------- the classification
def test_skill_is_kept_out_of_the_scores(tmp_path):
    write(tmp_path, "r.json", {"x": {"M1": cell(0.60, skill=0.9)}})
    scores, ratios, _reported, _claim = C.aucs(tmp_path)
    assert [k[-1] for k in ratios] == ["auc_skill"]
    assert "auc_skill" not in [k[-1] for k in scores]
    assert ("r.json", "x", "M1", "auc") in scores


def test_a_sample_size_is_not_a_score(tmp_path):
    write(tmp_path, "r.json", {"x": {"M1": cell(0.60)}})
    scores, _r, _rep, _c = C.aucs(tmp_path)
    named = {k[-1] for k in scores}
    assert "auc_n_pos" not in named and "auc_reached_pos" not in named


def test_a_cell_below_the_projects_reach_floor_is_not_reported(tmp_path):
    """The constant is the project's own, imported by the script, so this
    cannot drift from the two protocols that already share it."""
    from kgav.calibration import MIN_REACHED_POSITIVES
    few = (MIN_REACHED_POSITIVES - 1) / 100
    many = (MIN_REACHED_POSITIVES + 1) / 100
    write(tmp_path, "r.json", {"x": {
        "thin": cell(0.60, n_pos=100, cov_pos=few),
        "fat": cell(0.60, n_pos=100, cov_pos=many)}})
    _s, _r, reported, _c = C.aucs(tmp_path)
    assert ("r.json", "x", "fat") in reported
    assert ("r.json", "x", "thin") not in reported, (
        "a cell the project would refuse for reach is being judged")


def test_an_absent_flag_does_not_admit_an_unreadable_cell(tmp_path):
    """`trustworthy: False` and a missing `trustworthy` must not be the same
    thing, which is the bug that fired the gate on three reached positives."""
    write(tmp_path, "r.json", {"x": {
        "flagged": cell(0.60, trustworthy=False),
        "silent": cell(0.60, cov_pos=0.03, n_pos=100)}})
    _s, _r, reported, _c = C.aucs(tmp_path)
    assert ("r.json", "x", "flagged") not in reported
    assert ("r.json", "x", "silent") not in reported


def test_a_claim_is_an_interval_that_excludes_chance(tmp_path):
    write(tmp_path, "r.json", {"x": {
        "clears": cell(0.60, lo=0.55, hi=0.65),
        "spans": cell(0.52, lo=0.48, hi=0.56),
        "below": cell(0.40, lo=0.35, hi=0.45),
        "flagged": cell(0.60, clears_corrected_alpha=True)}})
    _s, _r, _rep, claim = C.aucs(tmp_path)
    assert ("r.json", "x", "clears") in claim
    assert ("r.json", "x", "below") in claim, (
        "an interval excluding chance BELOW 0.5 is a claim too")
    assert ("r.json", "x", "flagged") in claim
    assert ("r.json", "x", "spans") not in claim


# -------------------------------------------------------------- the refusals
def run(old: Path, new: Path, *extra: str) -> int:
    return C.main(["--against-dir", str(old), "--new-dir", str(new),
                   "--fail-below", "-0.01", *extra])


def test_it_refuses_a_claim_carrying_auc_that_fell(tmp_path, capsys):
    before, after = tmp_path / "b", tmp_path / "a"
    write(before, "r.json", {"x": {"M1": cell(0.70, lo=0.65, hi=0.75)}})
    write(after, "r.json", {"x": {"M1": cell(0.65, lo=0.60, hi=0.70)}})
    assert run(before, after) == 1
    out = capsys.readouterr().out
    assert "REFUSED" in out and "r.json/x/M1/auc" in out


def test_it_passes_a_claim_carrying_auc_that_rose(tmp_path, capsys):
    before, after = tmp_path / "b", tmp_path / "a"
    write(before, "r.json", {"x": {"M1": cell(0.65, lo=0.60, hi=0.70)}})
    write(after, "r.json", {"x": {"M1": cell(0.70, lo=0.65, hi=0.75)}})
    assert run(before, after) == 0
    assert "PASS" in capsys.readouterr().out


def test_a_cell_that_stops_carrying_its_claim_is_refused_at_any_size(
        tmp_path, capsys):
    """THE FAILURE A FLOOR ALONE MISSES. A cell can cross from clearing to
    not clearing on a move far smaller than any tolerance -- that is what
    the retraction in the strength-matched verdict was -- so the magnitude
    gate cannot be the only one."""
    before, after = tmp_path / "b", tmp_path / "a"
    write(before, "r.json", {"x": {"M1": cell(0.5100, lo=0.5010, hi=0.5190)}})
    write(after, "r.json", {"x": {"M1": cell(0.5095, lo=0.4995, hi=0.5195)}})
    assert run(before, after) == 1
    out = capsys.readouterr().out
    assert "stopped carrying a claim" in out
    assert "which no floor can excuse" in out


def test_a_falling_comparator_is_not_a_regression(tmp_path, capsys):
    """Both comparators, because the filter that was meant to cover them
    covered only one: `"no_pub_exclusion" in k` tested tuple MEMBERSHIP."""
    for name in ("degree", "structure_max_no_pub_exclusion"):
        before, after = tmp_path / f"b_{name}", tmp_path / f"a_{name}"
        write(before, "r.json", {"x": {name: cell(0.70, lo=0.65, hi=0.75)}})
        write(after, "r.json", {"x": {name: cell(0.60, lo=0.55, hi=0.65)}})
        assert run(before, after) == 0, (
            f"{name} fell 0.10 and the gate called it a regression; it is a "
            f"comparator, and its falling is the argument strengthening")
        assert "COMPARATORS" in capsys.readouterr().out


def test_a_falling_skill_alone_is_not_a_regression(tmp_path, capsys):
    """The ratio moves by half its range while the AUC holds."""
    before, after = tmp_path / "b", tmp_path / "a"
    write(before, "r.json",
          {"x": {"M1": cell(0.7000, lo=0.65, hi=0.75, skill=0.90)}})
    write(after, "r.json",
          {"x": {"M1": cell(0.7001, lo=0.65, hi=0.75, skill=0.10)}})
    assert run(before, after) == 0
    assert "auc_skill" in capsys.readouterr().out


def test_a_vanished_cell_is_named(tmp_path, capsys):
    """A measurement the release no longer supports shows up in no delta,
    which is how four temporal (run, virus) cells went below the
    reachability floor on 2026-11-11 with every figure within 0.003."""
    before, after = tmp_path / "b", tmp_path / "a"
    write(before, "r.json", {"x": {"M1": cell(0.60, lo=0.55, hi=0.65)},
                             "y": {"M1": cell(0.60, lo=0.55, hi=0.65)}})
    write(after, "r.json", {"x": {"M1": cell(0.60, lo=0.55, hi=0.65)}})
    assert run(before, after) == 0
    out = capsys.readouterr().out
    assert "EXIST BEFORE AND NOT AFTER" in out and "r.json / y / M1" in out


# ------------------------------------------------- against the real history
@pytest.mark.slow
def test_the_l27_repair_cost_no_claim_carrying_auc():
    """THE MEASUREMENT THIS SCRIPT WAS WRITTEN FOR, kept as a test so the
    answer cannot rot into a sentence. The identity repair of 2026-11-11
    moved a thousand figures; it must not have weakened a result."""
    import subprocess
    probe = subprocess.run(["git", "cat-file", "-e", "5151de8^{commit}"],
                           cwd=ROOT, capture_output=True)
    if probe.returncode != 0:
        pytest.skip("the pre-repair commit is not in this checkout")
    assert C.main(["--against", "5151de8", "--fail-below", "-0.01"]) == 0
