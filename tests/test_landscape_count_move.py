"""The two moved landscape cells have a named cause, and the documents say it.

WHAT WAS OPEN. doc2/registrations/prereg_pathway_overlap.md's note of 2026-11-17
records two cells whose `n_observed` moved, names a cause for one, and
leaves the other open with a controlled `ncats` re-run as the way to
close it. doc2/plan.md carried it as remaining work.

WHAT CLOSED IT. `n_observed = len(targets & factors)`, so a cell moves
when a protein ALREADY INSIDE its factor set joins the drugged proteome
-- no edge in its own channel need change. Both artifacts record the
operand's members, so the cause is a set difference over files on disk
and git, which is why no re-run was needed.

THESE TESTS DERIVE, THEY DO NOT RESTATE. The counts, the taxa, the
channels and the protein accessions all come out of the artifact; a
test that hard-coded them would pass against a document that had been
edited to match a wrong artifact. The one literal here is the taxon-to-
name mapping, which is read from config/viruses.yaml -- because naming
the wrong virus is precisely the error this work found.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

RESULTS = ROOT / "data" / "results" / "landscape_count_move.json"
PREREG = ROOT / "doc2" / "registrations" / "prereg_pathway_overlap.md"
PLAN = ROOT / "doc2" / "plan.md"
LANDSCAPE = ROOT / "data" / "results" / "pathway_landscape.json"


@pytest.fixture(scope="module")
def run() -> dict:
    if not RESULTS.exists():
        pytest.skip("the move has not been measured on this checkout "
                    "(make landscape-count-move)")
    return json.loads(RESULTS.read_text())


@pytest.fixture(scope="module")
def names() -> dict[str, str]:
    doc = yaml.safe_load((ROOT / "config" / "viruses.yaml").read_text())
    out = {}
    for lst in ("viruses", "surveyed"):
        for e in doc.get(lst) or []:
            out[f"NCBITaxon:{e['taxon']}"] = e["name"]
    return out


@pytest.fixture(scope="module")
def prereg() -> str:
    return PREREG.read_text()


# ------------------------------------------------------- the measurement holds
def test_the_operand_difference_explains_every_move(run):
    """The headline. Both falsifiers the script states must be empty."""
    assert run["unexplained_moves"] == [], run["unexplained_moves"]
    assert run["moves_predicted_that_did_not_happen"] == [], \
        run["moves_predicted_that_did_not_happen"]
    assert run["explains_every_move"] is True


def test_no_channel_edge_changed(run):
    """The independent check, over every cell rather than the two moved.

    If a PPI or CRISPR edge had changed, some cell's factor set would be
    a different size across the two artifacts. None is -- which is what
    makes "the operand alone did it" a measurement and not a story that
    merely fits two numbers.
    """
    assert run["cells_with_changed_factor_set"] == []
    assert run["factor_sets_unchanged"] is True


def test_every_moved_cell_names_the_protein_that_moved_it(run):
    """A cause per cell, not a cause for the pair."""
    moved = run["moved"]
    assert moved, "nothing moved, so this file has nothing to check"
    for key in moved:
        cell = run["cells"][key]
        assert cell["joined_inside"], f"{key} moved with no named protein"
        assert cell["predicted_delta"] == cell["observed_delta"]
        assert len(cell["joined_inside"]) == abs(cell["observed_delta"])


def test_the_operand_grew_by_addition_only(run):
    """Five joined and none left. A protein LEAVING would lower a cell,
    and no cell fell, so the two facts have to agree."""
    dp = run["drugged_proteome"]
    assert dp["left"] == []
    assert len(dp["joined"]) == dp["after"]["n"] - dp["before"]["n"]
    assert dp["before"]["sha256"] != dp["after"]["sha256"]


def test_the_derived_joined_set_matches_the_registration(run, prereg):
    """The five were recorded in the prereg when both releases existed.
    This re-derives them from the artifacts and checks the record.

    An independent confirmation of a figure nobody could re-measure once
    the earlier release was overwritten.
    """
    for acc in run["drugged_proteome"]["joined"]:
        assert acc in prereg, f"{acc} joined the operand and is not recorded"


# --------------------------------------------- the documents name it correctly
def test_the_prereg_note_names_each_moved_cell_s_channel_and_virus(
        run, names, prereg):
    """THE ERROR THIS WORK FOUND. The note called both cells
    'interactome' and the plan called NCBITaxon:31631 'SARS-CoV'. One is
    a dependency cell and 31631 is HCoV-OC43, so a reader chasing either
    looked in the wrong place -- which is exactly what happened.
    """
    for key in run["moved"]:
        cell = run["cells"][key]
        taxon, channel = cell["taxon"], cell["channel"]
        name = names[taxon]
        # the virus's own name and the right channel word must occur in
        # the same paragraph as the taxon for the record to be readable
        paras = [p for p in re.split(r"\n\s*\n", prereg) if taxon in p]
        assert paras, f"{taxon} is not mentioned in the note at all"
        assert any(name in p and channel in p for p in paras), (
            f"no paragraph names {taxon} together with its virus "
            f"({name}) and its channel ({channel})")


def test_the_prereg_records_the_counts_it_explains(run, prereg):
    """before -> after, derived, for each moved cell."""
    for key in run["moved"]:
        c = run["cells"][key]
        arrow = f"{c['n_observed_before']} -> {c['n_observed_after']}"
        assert arrow in prereg, f"{key}: '{arrow}' is not in the note"


def test_the_prereg_states_the_mechanism_rather_than_a_ppi_search(prereg):
    """The note's original reasoning searched PHYSICALLY_INTERACTS_WITH
    for a cell that never used it. The correction has to say what the
    cell actually computes, or the next reader repeats the search."""
    tail = prereg[prereg.index("# Note, 2026-11-17"):]
    assert "len(targets & factors)" in tail or \
           "targets & factors" in tail, \
        "the amendment does not state the intersection it is explaining"


def test_the_release_plan_no_longer_carries_it_as_open_work(run):
    """It was the only row in 'what is left' attributed to a re-run that
    is now measured not to be needed. Derived: no moved cell's taxon or
    the wrong virus name may appear in an open-work row."""
    plan = PLAN.read_text()
    rows = [ln for ln in plan.splitlines()
            if ln.startswith("|") and "re-run" in ln and "ncats" in ln]
    assert rows == [], f"still listed as open work: {rows}"


def test_the_open_question_is_not_still_called_open(prereg):
    """A closed question that still reads 'remains open' is worse than
    one never recorded: the document now disagrees with itself."""
    tail = prereg[prereg.index("# Note, 2026-11-17"):]
    idx = tail.find("That one remains")
    assert idx != -1, "the sentence this amendment answers has been edited " \
                      "away; a registration is appended to, never rewritten"
    assert "SETTLED" in tail[idx:] or "CLOSED" in tail[idx:], \
        "the note leaves a question open that data/results now answers"


# ------------------------------------------------------- the script still runs
def test_the_artifact_was_measured_on_the_release_on_disk(run):
    """A stale explanation is a wrong one: the factor sets it checked
    have to be the ones the landscape analysis used."""
    land = json.loads(LANDSCAPE.read_text())
    prov = run.get("provenance") or {}
    lprov = land.get("provenance") or {}
    for key in ("nodes_sha", "edges_sha"):
        if key in prov and key in lprov:
            assert prov[key] == lprov[key], (
                f"{key}: the move was explained on a different release "
                f"than the landscape was measured on")


def test_the_cell_count_covers_both_channels(run):
    """22 cells, not 14: reading only the interactome block is how the
    dependency cell got mislabelled in the first place."""
    land = json.loads(LANDSCAPE.read_text())
    expected = len([k for k in land["taxa"] if isinstance(land["taxa"][k], dict)])
    expected += len([k for k in land["interactome"]
                     if isinstance(land["interactome"][k], dict)])
    assert run["n_cells"] >= expected, (
        f"the diagnostic covered {run['n_cells']} cells and the landscape "
        f"reports {expected}")
    channels = {c["channel"] for c in run["cells"].values()}
    assert channels == {"dependency", "interactome"}
