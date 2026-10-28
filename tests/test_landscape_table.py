"""The pathway-landscape table in docs/evaluation.md, checked not sourced.

The sibling of tests/test_evaluation_table.py, and for the same reason.
Every one of the six degree lifts in the cross-sectional table drifted
across releases while the markdown sat unchanged, because a markdown
table has nothing to fail. This table carries twenty-two rows of
permutation p-values; hand-maintaining it across a rebuild is the same
bet, taken twenty-two times.

So the table is re-derived from data/results/pathway_landscape.json and
asserted cell by cell. It hard-codes no figure. The significance
threshold is not restated here either: it is `0.05 / n_tests` read from
the file, so a run that tests a different number of viruses moves the
bar and this test says so.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "evaluation.md"
RESULTS = ROOT / "data" / "results" / "pathway_landscape.json"
REGISTER = ROOT / "config" / "viruses.yaml"

MARKER = "<!-- TABLE: pathway_landscape"
DASH = "—"

# NO skipif ON A MISSING RESULTS FILE. data/results/*.json is tracked, and
# check_dir() globs the directory, so absence is a broken repository rather
# than an optional artifact. A guard that skips itself when its input is
# missing reports green having checked nothing -- the same defect as a
# fixture that never exercises the branch it names.


def names() -> dict[str, str]:
    """CURIE -> display name, from the register rather than a second list.

    The register keys on `taxon` holding a bare NCBI id. A hard-coded map
    here would be a fourth place the virus names live, and the register
    has already outgrown two of the other three.
    """
    reg = yaml.safe_load(REGISTER.read_text())
    return {f"NCBITaxon:{e['taxon']}": e["name"] for e in reg["viruses"]}


def derive() -> dict[str, list[str]]:
    """The table the results file implies, as rows of rendered cells."""
    d = json.loads(RESULTS.read_text())
    alpha = 0.05 / d["n_tests"]
    label = names()
    rows = {}
    for channel in ("taxa", "interactome"):
        for curie, r in d[channel].items():
            hi, lo = r["p_matched"], r["p_matched_less"]
            verdict = ("ENRICHED" if hi < alpha else
                       "DEPLETED" if lo < alpha else DASH)
            rows[f"{channel} | {label.get(curie, curie)}"] = [
                f"{r['n_dependency']:,}", f"{r['n_observed']:,}",
                f"{r['expected']:.2f}", f"{hi:.4f}", f"{lo:.4f}", verdict]
    return rows


def parse() -> dict[str, list[str]]:
    """The table as the document actually prints it."""
    text = DOC.read_text()
    assert MARKER in text, "the table lost its machine-readable marker"
    body = text.split(MARKER, 1)[1]
    lines = [ln.strip() for ln in body.splitlines()]
    start = next(i for i, ln in enumerate(lines) if ln.startswith("| channel |"))
    rows = {}
    for line in lines[start + 2:]:
        if not line.startswith("|"):
            break
        cells = [c.strip() for c in line.strip("|").split("|")]
        assert len(cells) == 8, f"expected 8 columns, got {len(cells)}: {line}"
        rows[f"{cells[0]} | {cells[1]}"] = cells[2:]
    return rows


@pytest.fixture(scope="module")
def tables():
    return parse(), derive()


def test_the_table_covers_exactly_the_rows_that_were_tested(tables):
    """Set equality both ways. A virus that becomes testable cannot appear
    in the results and be left out of the document, and a row dropped from
    the run cannot linger here as a figure nothing re-derives -- which is
    precisely how Zika went missing from the metapath table."""
    doc, got = tables
    assert set(doc) == set(got), (
        f"only in the document: {sorted(set(doc) - set(got))}; "
        f"only in the results: {sorted(set(got) - set(doc))}")


@pytest.mark.parametrize("column", range(6))
def test_every_cell_matches_the_results_file(tables, column):
    doc, got = tables
    labels = ["dependency factors", "observed", "expected",
              "p enriched", "p depleted", "verdict"]
    bad = {k: (doc[k][column], got[k][column])
           for k in sorted(set(doc) & set(got))
           if doc[k][column] != got[k][column]}
    assert not bad, (
        f"{labels[column]} is stale -- document vs re-derived: " +
        "; ".join(f"{k}: {d!r} != {g!r}" for k, (d, g) in bad.items()))


def test_the_stated_threshold_is_the_one_the_verdicts_used():
    """The document prints alpha to five decimals and the row count it
    came from. A run over a different number of viruses changes both, and
    a stale alpha beside live p-values is a table that reads as a result."""
    d = json.loads(RESULTS.read_text())
    text = DOC.read_text()
    alpha = 0.05 / d["n_tests"]
    assert f"n_tests = {d['n_tests']}" in text
    assert f"alpha = {alpha:.5f}" in text
    assert f"`degree_basis = {d['degree_basis']}`" in text
    assert f"{d['permutations']:,} permutations" in text


def test_the_headline_count_is_derived_not_typed():
    """"Zero of twenty-two" is the whole finding. If a row ever clears the
    threshold, this fails before anyone reads the sentence saying none
    did."""
    d = json.loads(RESULTS.read_text())
    alpha = 0.05 / d["n_tests"]
    rows = [r for ch in ("taxa", "interactome") for r in d[ch].values()]
    clearing = [r for r in rows
                if min(r["p_matched"], r["p_matched_less"]) < alpha]
    text = DOC.read_text()
    assert not clearing, (
        "a row now clears the threshold and the document still says none "
        f"does: {clearing}")
    assert f"**0 of {len(rows)}**" in text


def test_the_coverage_claim_counts_the_untested(tables):
    """The table lists what WAS tested, which is the half of the statement
    that flatters. Eight of sixteen viruses have no taxa-channel row at
    all, and a negative result over half a register is a different claim
    from a negative result over a register."""
    d = json.loads(RESULTS.read_text())
    reg = yaml.safe_load(REGISTER.read_text())
    total = len(reg["viruses"])
    text = DOC.read_text()
    assert len(d["taxa"]) + len(d["uncovered"]) == total
    assert len(d["interactome"]) + len(d["interactome_uncovered"]) == total
    assert f"{len(d['uncovered'])} of {total}" in text
    assert f"{len(d['interactome'])} of {total}" in text


def test_the_two_untestable_viruses_are_the_two_the_agent_refuses():
    """A cross-check between two things built separately: the statistics
    cannot reach DENV-4 or Yellow fever, and guards.check refuses a
    ranking for exactly those two. If a later ingest gives either one a
    route, the agent starts answering and this document must stop saying
    it cannot be tested -- so the two are pinned together here."""
    d = json.loads(RESULTS.read_text())
    label = names()
    refused = sorted(label[c] for c in d["interactome_uncovered"])
    assert refused == ["DENV-4", "Yellow fever"], refused
    text = DOC.read_text()
    for name in refused:
        assert name in text


def test_the_naive_contrast_quoted_in_the_prose_is_the_files_own():
    """The paragraph under the table quotes the uncorrected p-value beside
    the degree-matched one for two rows, because that gap IS the artifact
    being reported. Prose, so nothing re-derives it -- except this."""
    d = json.loads(RESULTS.read_text())
    text = DOC.read_text()
    for curie, name in (("NCBITaxon:11053", "DENV-1"),
                        ("NCBITaxon:11060", "DENV-2")):
        r = d["interactome"][curie]
        if "p_naive" not in r:           # older runs did not record it
            continue
        assert f"{r['p_naive']:.4f} naive" in text, name
        assert f"{r['p_matched']:.4f} matched" in text, name


def test_the_table_is_marked_as_machine_checked():
    text = DOC.read_text()
    assert text.count(MARKER) == 1
    assert "tests/test_landscape_table.py" in text
    head = text.split(MARKER, 1)[1].split("-->", 1)[0]
    assert "Do not hand-edit" in head
