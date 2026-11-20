"""The prospective calibration: loading it, and what the refusal may say.

WHAT WAS WRONG. The agent's candidate-generation refusal said, for every
virus and every question, "this release carries no measured prospective
result" -- while data/results/temporal_structure_scorer.json recorded ten
evaluated cells with six clearing chance.

The cause was a shape mismatch that failed silently.
scripts/kgav_ask.py read the prospective calibration out of
hard_negatives.json's `temporal` section with the cross-sectional loader.
That section is `{}` on this release (hard_negatives.py skipped the
temporal protocol and records that under `temporal_skipped`), so every
calibration came back unusable and the refusal took its no-measurement
branch. Pointing `calibration.load` at the prospective file instead would
not have helped: that file keys `cells` by taxon CURIE with each scorer
beside per-cell counts, and `load` expects {virus: {channel: row}}. It
returns an empty calibration -- no error, no note.

So guards._temporal_sentence had three branches that no caller could
reach, including the one written so that a reader is told about a scorer
that DOES clear chance before being refused.

THE REFUSAL ITSELF DOES NOT MOVE. guards.check returns allowed=False for
every PROSPECT query whatever the calibration says. These tests pin that
too, because a loader that made the agent start answering prospective
questions would be a far worse bug than the text being wrong.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from kgav import calibration as C
from kgav.agent.guards import QueryKind, check

ROOT = Path(__file__).resolve().parents[1]
TEMPORAL = ROOT / "data" / "results" / "temporal_structure_scorer.json"

needs_temporal = pytest.mark.skipif(
    not TEMPORAL.exists(), reason="no prospective results on disk")


def _doc() -> dict:
    return json.loads(TEMPORAL.read_text())


# ================================================== the loader, on real data
@needs_temporal
def test_it_loads_every_cell_the_file_evaluated():
    taxa = C.available_temporal_viruses(TEMPORAL)
    assert len(taxa) == _doc()["n_evaluable"]
    for t in taxa:
        assert C.load_temporal(TEMPORAL, t).usable, t


@needs_temporal
def test_the_discriminating_cells_are_exactly_the_files_own_clears_list():
    """THE PROPERTY THAT MATTERS.

    `discriminates` tests ci_lo > CHANCE, so whichever interval the
    loader puts in ci_lo decides what the agent announces. Handing it the
    95% interval instead of the corrected one would let the agent say
    cells clear chance that the registration's family size says do not --
    the agent claiming more than the study it quotes, out of the same
    file.

    Asserted as set equality against the file's own `clears`, which the
    registered run wrote, so the two cannot drift.
    """
    doc = _doc()
    clears = set(doc["clears"])
    got = {t for t in C.available_temporal_viruses(TEMPORAL)
           if C.load_temporal(TEMPORAL, t).discriminating()}
    assert got == clears, (
        f"the loader says {sorted(got)} clear and the results file says "
        f"{sorted(clears)}")
    assert clears, "no cell clears, so this test asserts nothing"


@needs_temporal
def test_it_uses_the_corrected_interval_and_not_the_95_percent_one():
    """Named directly, so the test above cannot pass for the wrong reason.

    There is at least one cell on this release whose uncorrected interval
    excludes chance and whose corrected one does not; that cell must not
    be reported as clearing.
    """
    doc = _doc()
    # EVERY SCORER, not just `primary`. The first version of this test
    # looked only at the primary column, found no disagreement and
    # skipped -- while DENV-1's symmetric_k10 reads ci_lo 0.5287 and
    # ci_lo_corrected 0.4901, which is precisely the case the test was
    # for. A test that skips is a test that asserts nothing.
    disagreeing = []
    for taxon, cell in doc["cells"].items():
        for scorer, row in cell.items():
            if not isinstance(row, dict) or scorer in C.TEMPORAL_EXCLUDE:
                continue
            lo, lo_c = row.get("auc_ci_lo"), row.get("auc_ci_lo_corrected")
            if lo is None or lo_c is None:
                continue
            if lo > C.CHANCE >= lo_c:
                disagreeing.append((taxon, scorer))
    assert disagreeing, (
        "no scorer on this release clears at 95% and fails the corrected "
        "alpha, so this test cannot distinguish the two intervals; point "
        "it at a run where one does")

    for taxon, scorer in disagreeing:
        cal = C.load_temporal(TEMPORAL, taxon)
        readable = C.TEMPORAL_SCORER_NAMES.get(scorer, scorer)
        assert readable not in cal.discriminating(), (
            f"{taxon}/{scorer} clears at 95% and not at the corrected "
            f"alpha, and the loader reports it as clearing: it read the "
            f"uncorrected interval")


@needs_temporal
def test_the_degree_baseline_is_not_loaded_as_a_scorer():
    """It is the null hypothesis -- what a scorer must beat, not a route
    to rank with. Same exclusion, same reason, as the cross-sectional
    loader's."""
    assert "degree" in C.TEMPORAL_EXCLUDE
    for t in C.available_temporal_viruses(TEMPORAL):
        cal = C.load_temporal(TEMPORAL, t)
        assert not any("degree" in n for n in cal.channels), t


@needs_temporal
def test_a_virus_the_study_never_evaluated_is_refused_not_zeroed():
    """A cell that does not exist is not a cell that scored chance."""
    assert "NCBITaxon:00000" not in _doc()["cells"], "pick a taxon not in the file"
    cal = C.load_temporal(TEMPORAL, "NCBITaxon:00000")
    assert not cal.usable
    assert "no prospective cell" in cal.refused_reason
    assert cal.channels == {}


@needs_temporal
def test_either_key_resolves_and_a_wrong_one_does_not():
    """The file keys by taxon and records a display name; a caller holding
    either must resolve, and a caller holding neither must miss rather
    than receive somebody else's numbers."""
    taxon = sorted(_doc()["cells"])[0]
    name = _doc()["cells"][taxon]["name"]
    by_taxon = C.load_temporal(TEMPORAL, taxon)
    by_name = C.load_temporal(TEMPORAL, name)
    assert by_taxon.usable and by_name.usable
    assert ({c.auc for c in by_taxon.channels.values()}
            == {c.auc for c in by_name.channels.values()})
    assert not C.load_temporal(TEMPORAL, "Marburg virus").usable


@needs_temporal
def test_the_protocol_string_names_the_split_it_came_from():
    """The agent prints this. A split built at another cutoff must not be
    describable with this one's wording."""
    doc = _doc()
    proto = C.temporal_protocol(doc)
    assert str(doc["cutoff"]) in proto
    assert doc["split"] in proto
    assert doc["undated_policy"] in proto


@needs_temporal
def test_the_split_sentence_agrees_with_the_cutoff_on_disk():
    """guards._SPLIT says "published before <year>" as prose, and the
    cutoff is a number in the results file. They are maintained in two
    places, so this is what keeps them from disagreeing: the prose year
    must be the cutoff year plus one, because the cutoff is inclusive.
    """
    import re

    from kgav.agent.guards import _SPLIT

    years = [int(y) for y in re.findall(r"\b((?:19|20)\d{2})\b", _SPLIT)]
    assert years, f"_SPLIT names no year: {_SPLIT!r}"
    assert years[0] == int(_doc()["cutoff"]) + 1, (
        f"_SPLIT says 'before {years[0]}' and the split on disk has "
        f"cutoff {_doc()['cutoff']} (inclusive), so the refusal would "
        f"misdescribe the protocol it quotes")


# ============================================ what the refusal may say
def _cal(auc, lo, hi, *, name="max-Tanimoto to pre-cutoff actives",
         protocol="temporal/v0.1-train2021/cutoff-2021"):
    cal = C.Calibration(virus="Zika", protocol=protocol)
    cal.channels[name] = C.ChannelCalibration(
        name=name, auc=auc, ci_lo=lo, ci_hi=hi, n_pos=567, n_neg=2495,
        coverage_pos=1.0, coverage_neg=1.0)
    return cal


def test_an_interval_below_chance_is_not_reported_as_chance():
    """THE BRANCH ADDED 2026-11-20.

    Zika's strongest prospective scorer reads AUC 0.422 with a corrected
    interval of [0.386, 0.458] -- an upper bound below 0.5, so the
    interval excludes chance on the wrong side. The old wording called
    that "every scorer performs at chance", which is not what was
    measured, and nothing could reach the branch to notice because no
    prospective calibration was ever loaded.
    """
    reason = check(QueryKind.PROSPECT,
                   temporal=_cal(0.422, 0.386, 0.458)).reason
    assert "performs at chance" not in reason
    assert "WORSE" in reason
    assert "wrong side" in reason
    assert "anti-correlated" in reason
    assert "0.422" in reason


def test_inverting_an_anti_correlated_scorer_is_not_offered_as_a_shortcut():
    """The obvious next thought, and it is not supported by anything
    measured: a ranking that anti-correlates on this split is not thereby
    shown to anti-correlate on the next one."""
    reason = check(QueryKind.PROSPECT,
                   temporal=_cal(0.422, 0.386, 0.458)).reason
    assert "inverting it is not a shortcut" in reason


def test_a_wide_interval_containing_chance_says_no_scorer_clears():
    """Not "every scorer performs at chance", which is a claim about every
    scorer. What is measured is that none of their intervals excludes
    chance -- weaker, and true. West Nile's strongest is 0.452
    [0.368, 0.536], an interval containing chance and plenty besides."""
    reason = check(QueryKind.PROSPECT,
                   temporal=_cal(0.452, 0.368, 0.536)).reason
    assert "no scorer clears chance" in reason
    assert "performs at chance" not in reason
    assert "WORSE" not in reason


def test_a_clearing_scorer_is_still_announced():
    """The branch that existed and was unreachable."""
    reason = check(QueryKind.PROSPECT,
                   temporal=_cal(0.802, 0.700, 0.905)).reason
    assert "which clears chance" in reason
    assert "0.802" in reason


@pytest.mark.parametrize("auc,lo,hi", [(0.802, 0.700, 0.905),
                                       (0.452, 0.368, 0.536),
                                       (0.422, 0.386, 0.458)])
def test_the_refusal_is_still_a_refusal_whatever_was_measured(auc, lo, hi):
    """THE SAFETY PROPERTY, pinned across all three branches.

    A loader that made the agent start answering candidate-generation
    questions because one dengue cell cleared chance would be a far worse
    defect than the text being wrong. The coverage argument is what the
    refusal rests on, and it does not depend on any AUC.
    """
    v = check(QueryKind.PROSPECT, temporal=_cal(auc, lo, hi))
    assert v.allowed is False
    assert "novel chemical matter" in v.reason


# ============================================== the wiring, end to end
def test_the_cli_defaults_to_the_prospective_file():
    """The default is the fix. Left at hard_negatives.json's empty
    `temporal` section, every branch above stays unreachable."""
    import importlib.util
    import sys

    spec = importlib.util.spec_from_file_location(
        "kgav_ask_temporal_wiring", ROOT / "scripts" / "kgav_ask.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)

    args = mod.build_parser(ROOT).parse_args(["why"])
    assert args.temporal_results.name == "temporal_structure_scorer.json"
    assert args.temporal_results.parent.name == "results"
    # And it must be a different file from --results, which is the whole
    # point: the prospective study is not a section of the other one.
    assert args.temporal_results != args.results


@needs_temporal
def test_the_prospective_section_of_hard_negatives_is_still_empty():
    """THE PREMISE OF THE FIX, asserted so it cannot rot into a test that
    passes for the wrong reason.

    If hard_negatives.py ever runs its temporal protocol, that section
    becomes populated and this file's claim -- that the old loader found
    nothing -- stops being the explanation. Then the two sources need
    reconciling rather than one silently winning.
    """
    hn = ROOT / "data" / "results" / "hard_negatives.json"
    if not hn.exists():
        pytest.skip("no hard_negatives.json on disk")
    doc = json.loads(hn.read_text())
    assert doc.get("temporal") == {}, (
        "hard_negatives.json now carries a temporal section; kgav_ask.py "
        "loads the prospective calibration from "
        "temporal_structure_scorer.json and would ignore it. Decide which "
        "is authoritative instead of leaving two.")
    assert "temporal_skipped" in doc, (
        "the file used to record WHY its temporal section is empty; that "
        "record is what distinguishes 'skipped' from 'measured nothing'")
