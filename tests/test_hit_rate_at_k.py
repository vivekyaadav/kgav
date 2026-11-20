"""The hit-rate report, and the one thing it must not be allowed to say.

WHY THIS FILE EXISTS. scripts/repurposing_candidates.py has printed, above
every shortlist it has ever produced, "every row below is a chemical-
analogy suggestion with a measured hit rate" -- while nothing in
data/results recorded a hit rate at all. The banner quoted AUC, which is
the probability a random active outranks a random inactive and is the one
translation that does not give a hit rate.

scripts/hit_rate_at_k.py measures the real quantity. These tests hold it
to three properties:

  1. It types no figure. Same rule as every other reporting script here.
  2. It reports a precision beside the BASE RATE of the cell, always. A
     precision of 90% in a cell where 70% of the screened compounds are
     active is a ranking that added almost nothing, and a report that
     prints the first number without the second is worse than silent.
  3. It is labelled post-hoc, and it does not write into data/results.
     The registration fixed AUC as the statistic before any AUC existed;
     this was computed after reading them, so it is descriptive and must
     not be presentable as a registered result.
"""
from __future__ import annotations

import importlib.util
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "hit_rate_at_k.py"
SPLIT = ROOT / "data" / "releases" / "v0.1-train2021"

needs_split = pytest.mark.skipif(
    not (SPLIT / "TEST_LABELS.json").exists(),
    reason="no 2021 temporal split on disk")


def _module():
    spec = importlib.util.spec_from_file_location(
        "hit_rate_under_test", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------- hygiene
def test_the_script_types_no_figure():
    """Three or more decimal places is a measurement, not a setting."""
    bad = [m.group(0) for m in
           re.finditer(r"\b\d+\.\d{3,}\b", SCRIPT.read_text())]
    assert not bad, f"{SCRIPT.name} types the figure(s) {bad}"


@pytest.mark.slow
@needs_split
def test_running_it_writes_nothing_into_results(tmp_path):
    """data/results is swept by `make check-results`, and everything in it
    reads as a claim-carrying result stamped to a release. A post-hoc
    statistic that landed there would acquire exactly the standing this
    script's own docstring disclaims.

    ASKED BY RUNNING IT, not by reading it. The first version of this test
    searched the source for the string "data/results" and failed on the
    docstring, which names the four files it checked for an existing hit
    rate. A test that cannot tell a prose mention from a write would be
    satisfied by deleting the prose.
    """
    import hashlib

    results = ROOT / "data" / "results"

    def snapshot():
        return {f.name: hashlib.sha256(f.read_bytes()).hexdigest()
                for f in sorted(results.glob("*.json"))}

    before = snapshot()
    assert before, "no results files to protect; this test asserts nothing"
    out = subprocess.run([sys.executable, str(SCRIPT)], capture_output=True,
                         text=True, cwd=ROOT, timeout=1800)
    assert out.returncode == 0, out.stderr
    after = snapshot()

    assert after == before, (
        "running the hit-rate report changed "
        f"{sorted(set(before) ^ set(after)) or [k for k in before if before[k] != after.get(k)]}")


@needs_split
def test_out_is_opt_in_and_goes_where_it_is_told(tmp_path):
    """--out exists, so the table can be kept; it must not default to a
    path, and what it writes must carry the post-hoc label."""
    import json

    dest = tmp_path / "hit_rate.json"
    out = subprocess.run(
        [sys.executable, str(SCRIPT), "--out", str(dest), "--k", "10"],
        capture_output=True, text=True, cwd=ROOT, timeout=1800)
    assert out.returncode == 0, out.stderr
    assert dest.exists(), "--out wrote nothing"
    doc = json.loads(dest.read_text())
    assert doc["post_hoc"] is True
    assert "not_registered" in doc
    assert doc["cells"], "no cells in the written table"


def test_it_says_in_its_own_text_that_it_is_post_hoc():
    """A reader who runs this and pastes the table somewhere must carry
    the caveat with it, so the caveat is in the output and not only in
    the source."""
    src = SCRIPT.read_text()
    assert "post-hoc" in src.lower()
    assert "not_registered" in src or "post_hoc" in src


# ------------------------------------------------- the arithmetic itself
def test_precision_is_hits_over_k_and_enrichment_is_over_the_base_rate():
    """Asserted on a hand-made cell, so the formula is pinned independently
    of any release.

    Ten candidates, four of them active, all distinctly scored: the top
    five hold three actives. precision@5 is 3/5 and the base rate is
    4/10, so enrichment is 0.6 / 0.4.
    """
    mod = _module()
    cand = [f"C{i}" for i in range(10)]
    pos = ["C0", "C1", "C4", "C9"]
    scores = {f"C{i}": 1.0 - i / 100 for i in range(10)}   # C0 best
    r = mod.precision_rows(cand, pos, scores, (5,))
    m = r["per_k"][5]
    assert r["base_rate"] == pytest.approx(0.4)
    assert m["hits_lo"] == m["hits_hi"] == 3, "C0, C1 and C4 are the top five"
    assert m["precision_expected"] == pytest.approx(0.6)
    assert m["enrichment_expected"] == pytest.approx(0.6 / 0.4)


def test_an_unscored_candidate_ranks_last_and_is_not_dropped():
    """score_cell returns no entry for a candidate it could not score.
    Dropping those would shrink the pool and inflate every precision; the
    AUC evaluator imputes them to the floor and so does this."""
    mod = _module()
    cand = ["A", "B", "C", "D"]
    pos = ["D"]
    r = mod.precision_rows(cand, pos, {"A": 0.9, "B": 0.8, "C": 0.7}, (3,))
    assert r["pool"] == 4, "the unscored candidate left the pool"
    assert r["base_rate"] == pytest.approx(0.25)
    assert r["per_k"][3]["hits_lo"] == 0, "D outranked something scored"
    assert mod.UNSCORED < 0, "Tanimoto is non-negative; the floor must be below"


def test_a_tie_straddling_the_cut_is_reported_as_a_range():
    """The project's rule: a point estimate at the centre of a range reads
    as a measurement. Four candidates tied at one score, two positive,
    with two slots -- hits@2 is anywhere in [0, 2]."""
    mod = _module()
    cand = ["A", "B", "C", "D"]
    r = mod.precision_rows(cand, ["A", "B"], dict.fromkeys(cand, 0.5), (2,))
    m = r["per_k"][2]
    assert not m["determined"]
    assert (m["hits_lo"], m["hits_hi"]) == (0, 2)
    assert m["hits_expected"] == pytest.approx(1.0)
    assert m["tied_at_cut"] == 4 and m["tied_at_cut_positives"] == 2


# ------------------------------------------------- against the real split
@pytest.mark.slow
@needs_split
def test_the_report_prints_a_base_rate_beside_every_precision():
    """PROPERTY 2, on the real run. A precision column with no base rate
    beside it is the number that would be quoted on its own."""
    out = subprocess.run([sys.executable, str(SCRIPT)], capture_output=True,
                         text=True, cwd=ROOT, timeout=1800)
    assert out.returncode == 0, out.stderr
    assert "base" in out.stdout
    assert "enrich" in out.stdout
    # The caveat that makes the table readable must be in the output.
    assert "BASE RATE IS DOING MOST OF THE WORK" in out.stdout
    assert "UPPER BOUND" in out.stdout


@pytest.mark.slow
@needs_split
def test_the_plate_ladder_shows_every_exclusion_strength():
    """The measurement that decides the question: how far the precision
    falls when a candidate's nearest known actives are withheld. One rung
    of that ladder, reported alone, is what made the shortlist's banner
    plausible."""
    out = subprocess.run([sys.executable, str(SCRIPT)], capture_output=True,
                         text=True, cwd=ROOT, timeout=1800).stdout
    _, _, plate = out.partition("LIBRARY-LIKE")
    assert plate, "the plate block is gone"
    assert "no exclusion" in plate
    assert "k10" in plate
    assert "k_cell=" in plate
