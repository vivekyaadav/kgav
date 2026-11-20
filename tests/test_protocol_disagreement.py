"""DENV-4 clears prospectively and anti-ranks cross-sectionally.

The pair reads as a contradiction. It is not one: the two protocols use
different REFERENCE SETS, and for a post-cutoff active the
cross-sectional control removes its own chemical series while the
prospective one cannot, because that series is not in the pre-cutoff
reference set at all.

WHY IT NEEDED RECORDING. A reader of
docs/prereg_temporal_structure_scorer.md alone sees DENV-4 at 0.7896,
clearing, and would reasonably take that to license a shortlist. It does
not: the shortlist ranks never-screened compounds against ALL known
actives, which is the cross-sectional setting -- the one where DENV-4 is
below chance and scripts/repurposing_candidates.py refuses it.

These tests hold the amendment's figures to the diagnostic that produced
them, and hold the diagnostic to the two registered results it compares.
The amendment is kept honest in both directions: its numbers must be
re-derivable, and the conclusion it draws must stay true of the data.
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "prereg_temporal_structure_scorer.md"
SCRIPT = ROOT / "scripts" / "diagnose_protocol_disagreement.py"
RESULTS = ROOT / "data" / "results"
SPLIT = ROOT / "data" / "releases" / "v0.1-train2021"
DENV4 = "NCBITaxon:11070"

needs_split = pytest.mark.skipif(
    not (SPLIT / "TEST_LABELS.json").exists(),
    reason="no 2021 temporal split on disk")


def _flat() -> str:
    return " ".join(DOC.read_text().split())


@pytest.fixture(scope="module")
def measured() -> dict:
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "d.json"
        r = subprocess.run(
            [sys.executable, str(SCRIPT), "--out", str(out)],
            capture_output=True, text=True, cwd=ROOT, timeout=1800)
        assert r.returncode == 0, r.stdout + r.stderr
        return json.loads(out.read_text())


# ------------------------------------------- the premise still holds
@needs_split
def test_the_two_protocols_still_disagree(measured):
    """THE GUARD ON EVERYTHING BELOW. If a rebuild moved either result
    the amendment describes a disagreement that no longer exists, and
    re-deriving its figures would be checking a fossil."""
    cross = measured["cross_sectional_by_era"]["all_vs_all"]
    pro = measured["prospective"]["primary"]
    assert cross < 0.5, f"DENV-4 no longer anti-ranks cross-sectionally: {cross}"
    assert measured["prospective"]["clears"], (
        "DENV-4 no longer clears prospectively; the amendment is about a "
        "disagreement that is gone")
    assert pro > cross


@needs_split
def test_the_shortlist_tool_still_refuses_this_cell(measured):
    """The amendment's practical conclusion. If the tool stopped
    refusing DENV-4, 'a prospective pass does not license a shortlist'
    would be advice about a thing that no longer happens."""
    sys.path.insert(0, str(ROOT / "scripts"))
    sys.path.insert(0, str(ROOT / "src"))
    from repurposing_candidates import tier, verdicts

    assert tier(verdicts().get(DENV4)) == "below", (
        "DENV-4 is no longer refused by the shortlist tool")


# ------------------------------------- the amendment's figures, derived
@needs_split
def test_every_auc_in_the_amendment_is_the_diagnostics(measured):
    e, p = measured["cross_sectional_by_era"], measured["prospective"]
    expected = {f"{v:.4f}" for v in
                (e["all_vs_all"], e["post_vs_post"], e["pre_vs_pre"],
                 e["pre_pos_vs_post_neg"], p["primary"], p["symmetric_k10"])}
    flat = _flat()

    # EVERY ONE PRESENT ...
    missing = sorted(v for v in expected if v not in flat)
    assert not missing, f"the amendment omits the measured AUC(s) {missing}"

    # ... AND EVERY FOUR-DECIMAL FIGURE IN THE AMENDMENT MEASURED.
    #
    # Checking only presence passes while a second occurrence says
    # something else: 0.4842 appears twice -- in the era table and in
    # the reference-set comparison -- so perturbing one left the other
    # to satisfy the test. The same hole was found in the claim-block
    # tests and again in the DENV-2 registration's power table, both
    # within two days of this.
    amendment = DOC.read_text().split(
        "## Amendment, 2026-11-21", 1)[1]
    printed = set(re.findall(r"\b0\.\d{4}\b", amendment))
    stray = sorted(printed - expected)
    assert not stray, (
        f"the amendment prints {stray}, which the diagnostic does not "
        f"produce")


@needs_split
def test_the_nesting_claim_is_the_measured_one(measured):
    """'It is not the populations' rests on an exact nesting. A partial
    overlap would make the comparison in the next section invalid."""
    c, n = measured["counts"], measured["nested"]
    assert n["test_positives_in_cross_sectional"] == c["post_cutoff_positives"]
    assert n["references_in_cross_sectional"] == c["pre_cutoff_positives"]
    assert n["train_test_overlap"] == 0, (
        "train and test overlap, so the split's own overlap policy is not "
        "holding and the amendment's nesting claim is wrong")
    assert (c["pre_cutoff_positives"] + c["post_cutoff_positives"]
            == c["cross_sectional_positives"]), (
        "the cross-sectional positives are not exactly the two eras")

    flat = _flat()
    assert f"{n['test_positives_in_cross_sectional']} of" in flat
    assert f"{c['cross_sectional_positives']} positives" in flat


@needs_split
def test_the_era_mate_share_is_the_measured_one(measured):
    """The mechanism. Without it the amendment asserts a cause."""
    s = measured["era_mate_share"]
    flat = _flat()
    assert f"median {s['median']:.0%}" in flat
    assert f"mean\n{s['mean']:.0%}" in DOC.read_text() or \
           f"mean {s['mean']:.0%}" in flat
    assert f"{s['n_compounds']} compounds" in flat
    assert s["median"] > 0.5, (
        "the nearest references are no longer mostly era-mates, so the "
        "mechanism the amendment gives is no longer the mechanism")


@needs_split
def test_the_k_cell_is_the_registered_one(measured):
    """The control's strength is the thing doing the removing, so the
    amendment names it and it must be the run's own."""
    k = json.loads((RESULTS / "structure_scorer_strength_matched.json")
                   .read_text())["k_cells"][DENV4]
    assert measured["k_cell"] == k
    assert f"k_cell {k}" in _flat()


# ------------------------------------------- what it refuses to conclude
def test_the_amendment_forecloses_the_shortlist_reading():
    """The reason it exists. A reader of the registration alone sees
    0.7896 and would take it to license a list."""
    flat = _flat()
    assert "does not license a shortlist" in flat
    assert "rests on cross-era analogy" in flat
    assert "That refusal is correct and stays" in flat


def test_the_amendment_keeps_both_results(measured):
    """Neither is withdrawn. An amendment that quietly retired one would
    be rewriting a signed result rather than explaining it."""
    flat = _flat()
    assert "Both results stand" in flat
    assert "answer different questions" in flat


def test_it_names_the_command_that_reproduces_it():
    flat = _flat()
    assert "diagnose_protocol_disagreement.py" in flat
    assert "make protocol-disagreement" in flat


def test_the_diagnostic_writes_nothing_into_results(measured):
    """Same rule as every post-hoc measurement here: it explains a
    registered result and is not itself one."""
    src = SCRIPT.read_text()
    assert "--out" in src
    assert not re.search(r'default\s*=\s*RESULTS', src), (
        "the diagnostic defaults to writing into data/results")
    assert "Nothing here is written to data/results" in src
