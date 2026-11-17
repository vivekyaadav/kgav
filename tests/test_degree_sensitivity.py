"""docs/evaluation.md's M7 claim states its own weighting sensitivity.

THE CLAIM SECTION RESTS ON ONE CELL. Of thirteen readable cells, exactly
one clears the corrected alpha -- SARS-CoV-2 M7 -- and its DWPC step
weights read `Graph.degree`, which counts the held-out label edges the
channel is scored against. That the claim SURVIVES taking them out is a
fact about the data, not a property anyone could assert, so the document
has to state the alternative figure and this derives it.

The sensitivity is not a protocol: nothing promised that step weights
exclude held-out edges, and baselines.py argues the opposite on purpose.
So the artifact registers nothing, and these tests check the document
reports it rather than that any particular weighting is right.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

DOC = ROOT / "docs" / "evaluation.md"
RESULTS = ROOT / "data" / "results" / "degree_weighting_sensitivity.json"


@pytest.fixture(scope="module")
def run() -> dict:
    if not RESULTS.exists():
        pytest.skip("the sensitivity has not been measured on this checkout")
    return json.loads(RESULTS.read_text())


@pytest.fixture(scope="module")
def doc() -> str:
    return DOC.read_text()


def test_the_claim_cell_is_still_the_only_one_clearing_correction(run):
    """The premise. If a second cell cleared, the document's "only M7"
    sentence would be wrong before any sensitivity mattered."""
    import json as _j
    mp = _j.loads((ROOT / "data" / "results"
                   / "metapath_comparison.json").read_text())["per_virus"]
    readable = [(v, n, m) for v, ch in mp.items() for n, m in ch.items()
                if n not in ("degree", "COMBINED")
                and m.get("hits_determined") and m.get("trustworthy")]
    alpha = 0.05 / len(readable)
    clearing = sorted((v, n) for v, n, m in readable
                      if m.get("p_value") is not None and m["p_value"] < alpha)
    # metapath_comparison.json keys per_virus by LABEL and this artifact
    # by CURIE, so the register does the translation rather than a
    # hand-written pair.
    from kgav import virus_register as vr
    names = {"NCBITaxon:" + v["taxon"]: v["name"]
             for v in vr.load()["viruses"]}
    want = (names[run["claim_cell"]["virus"]], run["claim_cell"]["channel"])
    assert clearing == [want], (
        f"{clearing} clear 0.05/{len(readable)}; the artifact names "
        f"{want} as the claim cell")


def test_both_p_values_for_the_claim_cell_are_in_the_document(run, doc):
    cell = run["claim_cell"]
    for key in ("p_a", "p_b"):
        p = cell[key]
        assert p is not None, f"the artifact carries no {key}"
        assert f"{p:.6f}" in doc, (
            f"the claim cell's {key} = {p:.6f} is not stated in "
            f"docs/evaluation.md")


def test_the_document_states_which_side_of_the_threshold_arm_b_falls_on(
        run, doc):
    """A sensitivity that does not say whether the claim survives is not a
    sensitivity. The multiple of alpha is derived, not typed."""
    import json as _j
    mp = _j.loads((ROOT / "data" / "results"
                   / "metapath_comparison.json").read_text())["per_virus"]
    readable = [m for v, ch in mp.items() for n, m in ch.items()
                if n not in ("degree", "COMBINED")
                and m.get("hits_determined") and m.get("trustworthy")]
    alpha = 0.05 / len(readable)
    p_b = run["claim_cell"]["p_b"]
    assert p_b < alpha, (
        f"the claim no longer clears under arm B ({p_b:.6f} against "
        f"{alpha:.5f}); the document's sensitivity paragraph says it does "
        f"and must be rewritten, not reinterpreted")
    # STATED THE WAY THIS DOCUMENT ALREADY STATES A MARGIN: alpha / p,
    # rounded, as "a factor of N". tests/test_evaluation_table.py holds
    # the published figure to the same convention, so the sensitivity
    # cannot be written in a flattering unit.
    assert f"a factor of {round(alpha / p_b)}" in doc, (
        f"the document does not state arm B's margin of "
        f"{round(alpha / p_b)}x the corrected alpha")


def test_the_number_of_moving_cells_and_the_crossings_are_the_runs(run, doc):
    """FOUR OF SIXTY-THREE, AND NONE CROSSING, is what makes the finding
    small. A document that stated the claim's own move without the
    denominator would read as a general problem."""
    assert f"{run['cells_moved']} of the {run['cells']}" in doc, (
        f"{run['cells_moved']} of {run['cells']} cells move and the "
        f"document states something else")
    assert run["cells_crossing_alpha"] == 0, (
        "a cell now crosses p=0.05 between the arms; that is a result, not "
        "a sensitivity footnote")
    assert "no cell crosses" in doc.lower()


def test_the_hub_set_is_identical_across_the_arms(run, doc):
    """Arm B isolates the weighting only if the hub cut does not move with
    it. It does not, because every held-degree Protein is viral while the
    cut is class-scoped to host proteins."""
    sizes = run["hub_set_size"]
    assert sizes["A"] == sizes["B"] == sizes["C"], sizes
    assert str(sizes["A"]) in doc


def test_the_artifact_registers_nothing_and_says_so(run):
    assert "SENSITIVITY, NOT A PROTOCOL" in run["registers_nothing"]
    assert "registered_by" not in run, (
        "this artifact must not claim a registration; switching the "
        "weighting is a decision that needs its own")


def test_the_document_does_not_present_the_alternative_as_the_figure(doc):
    """The published number stays the published number.

    data/results holds arm A, so the claim section must quote arm A as its
    result and arm B as a sensitivity. Swapping them would restate a
    figure no results file carries.
    """
    section = doc.split("### The one readable cell")[1].split("##")[0]
    a_first = section.index("0.000841") < section.index("0.001788")
    assert a_first, (
        "the alternative weighting is quoted before the published one, so "
        "the section reads as though data/results held arm B")
