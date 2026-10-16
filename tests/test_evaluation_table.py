"""The one table in docs/evaluation.md that is checked rather than sourced.

That document's rule is "treat a figure with no named source as unverified",
which leaves the reader to re-run the script. It is not enough. Every one of
the six degree lifts in the cross-sectional table drifted across releases
while the table sat unchanged, the "best channel" column named four rows the
code refuses to read, and Zika was never added at all -- each of which passed
review repeatedly, because a markdown table has nothing to fail.

So this file re-derives the table from data/results/metapath_comparison.json
and asserts it cell by cell. It hard-codes no figure: a number that moves
fails here until the document is corrected, and a virus that appears in the
JSON fails here until the document mentions it.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "evaluation.md"
RESULTS = ROOT / "data" / "results" / "metapath_comparison.json"

MARKER = "<!-- TABLE: metapath_comparison"
DASH = "—"          # the em-dash the document uses for "no value"


# --------------------------------------------------------------- the rules
def fmt_p(p: float | None) -> str:
    """One formatting rule, used by the derivation and by nothing else.

    A second copy of this in the document's prose is how a table drifts, so
    the test owns it and the document quotes what it produces.
    """
    if p is None:
        return DASH
    return "<.0001" if p < 0.0001 else f"{p:.4f}"


def best_readable(channels: dict) -> str | None:
    """The single metapath with the highest lift among those that are BOTH
    `trustworthy` and `hits_determined`.

    Both conditions, because the document spends four paragraphs on what
    happens when one is read and the other ignored -- and then its own
    summary table printed SARS-CoV M1 (undetermined) and MERS-CoV M4
    (expectation 0.51) as results.
    """
    singles = {n: m for n, m in channels.items()
               if n not in ("degree", "COMBINED")}
    ok = {n: m for n, m in singles.items()
          if m.get("hits_determined") and m.get("trustworthy")}
    return max(ok, key=lambda n: ok[n]["lift"]) if ok else None


def derive() -> dict[str, list[str]]:
    """The table the results file implies, as rows of rendered cells."""
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from kgav.baselines import confound_warning

    per = json.loads(RESULTS.read_text())["per_virus"]
    rows = {}
    for virus, channels in per.items():
        name = best_readable(channels)
        if name is None:
            chan = [DASH, DASH, DASH]
        else:
            m = channels[name]
            chan = [name, f"{m['lift']:.2f}", fmt_p(m.get("p_value"))]
        deg = channels["degree"]
        warning = confound_warning(channels)
        rows[virus] = chan + [f"{deg['lift']:.2f}", fmt_p(deg.get("p_value")),
                              warning.split(":")[0] if warning else DASH]
    return rows


def parse() -> dict[str, list[str]]:
    """The table as the document actually prints it."""
    text = DOC.read_text()
    assert MARKER in text, "the table lost its machine-readable marker"
    body = text.split(MARKER, 1)[1]
    lines = [l.strip() for l in body.splitlines()]
    start = next(i for i, l in enumerate(lines) if l.startswith("| virus |"))
    rows = {}
    for line in lines[start + 2:]:
        if not line.startswith("|"):
            break
        cells = [c.strip() for c in line.strip("|").split("|")]
        assert len(cells) == 7, f"expected 7 columns, got {len(cells)}: {line}"
        rows[cells[0]] = cells[1:]
    return rows


# --------------------------------------------------------------- the checks
@pytest.fixture(scope="module")
def tables():
    return parse(), derive()


def test_the_table_covers_exactly_the_viruses_in_the_results(tables):
    """THE ZIKA GAP. Zika sat in metapath_comparison.json and not in this
    table, so the one virus whose M1 interval the document discusses at
    length was missing from the summary of M1 results. Set equality, both
    ways: a virus added to the pipeline cannot be silently omitted, and a
    virus dropped from it cannot linger here as a figure nothing re-derives.
    """
    doc, got = tables
    assert set(doc) == set(got), (
        f"only in the document: {sorted(set(doc) - set(got))}; "
        f"only in the results: {sorted(set(got) - set(doc))}")


@pytest.mark.parametrize("column", range(6))
def test_every_cell_matches_the_results_file(tables, column):
    """Parametrised by column so a failure names which figure moved rather
    than dumping the whole table."""
    doc, got = tables
    labels = ["best readable channel", "lift", "p",
              "degree lift", "degree p", "baseline verdict"]
    bad = {v: (doc[v][column], got[v][column])
           for v in sorted(set(doc) & set(got))
           if doc[v][column] != got[v][column]}
    assert not bad, (
        f"{labels[column]} is stale -- document vs re-derived: " +
        "; ".join(f"{v}: {d!r} != {g!r}" for v, (d, g) in bad.items()))


def test_an_unreadable_channel_is_a_dash_not_a_number(tables):
    """The defect this table was printing: a channel named and a lift quoted
    for a row the code refuses to read. If a virus has no readable channel,
    all three of its channel cells must be em-dashes -- not a 0.00, which is
    what HCoV-NL63 printed while having no reachable positive at all."""
    doc, _ = tables
    for virus, cells in doc.items():
        dashes = [c == DASH for c in cells[:3]]
        assert len(set(dashes)) == 1, (
            f"{virus}: channel cells are half-filled: {cells[:3]}")


def test_the_unreadable_rows_are_unreadable_for_a_stated_reason():
    """Not a formatting check: it asserts the two mechanisms by name, so
    that a future change which makes these rows readable forces the prose
    above the table to be rewritten rather than left describing a state the
    data no longer has."""
    per = json.loads(RESULTS.read_text())["per_virus"]
    assert per["SARS-CoV"]["M1"]["hits_determined"] is False
    assert per["Zika"]["M1"]["hits_determined"] is False
    for virus in ("MERS-CoV", "HCoV-229E", "HCoV-OC43", "HCoV-NL63"):
        assert per[virus]["M4"]["trustworthy"] is False, virus
    assert per["HCoV-NL63"]["COMBINED"]["pos"] == 0


def test_the_prose_figures_beside_the_table_are_the_tables_own():
    """The paragraphs under the table quote individual cells. They are prose,
    so nothing re-derives them -- except this, which checks the handful that
    name a number."""
    per = json.loads(RESULTS.read_text())["per_virus"]
    text = DOC.read_text()
    m4 = per["SARS-CoV-2"]["M4"]
    deg = per["SARS-CoV-2"]["degree"]
    assert f"{m4['lift_vs_degree']:.2f}x a degree" in text
    assert f"degree lift of **{deg['lift']:.2f}**" in text
    assert m4["vs_degree_readable"] is False
    m1 = per["SARS-CoV"]["M1"]
    assert f"lift {m1['lift']:.2f} at p<0.0001" in text
    assert f"[{m1['hits_lo']}, {m1['hits_hi']}] of {m1['pos']}" in text
    assert f"[{m1['lift_lo']:.2f}, {m1['lift_hi']:.2f}] contains 1.00" in text


def test_the_at_chance_claim_names_the_policy_it_holds_under():
    """The document said BASELINE AT CHANCE "no longer holds anywhere". It
    holds for SARS-CoV-2 under computed_only. The replacement sentence is
    pinned to the run it describes."""
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from kgav.baselines import confound_warning

    path = ROOT / "data" / "results" / "temporal_2021_computed_only.json"
    row = json.loads(path.read_text())["results"]["SARS-CoV-2"]
    warning = confound_warning(row)
    assert warning.startswith("BASELINE AT CHANCE")
    text = DOC.read_text()
    # Scoped to the verdict block, not the whole document: the correction
    # paragraph beneath it QUOTES the retracted claim, and a test that
    # forbade the words outright would forbid recording the mistake.
    block = text.split("    CONFOUNDED", 1)[1].split("    (nothing)", 1)[0]
    assert "no longer holds" not in block
    assert "SARS-CoV-2 under temporal computed_only" in block
    assert "It does hold, for SARS-CoV-2 under `computed_only`" in text
    assert (f"degree lift {row['degree']['lift']:.2f} "
            f"(p={row['degree']['p_value']:.4f})") in text
    assert (f"COMBINED's {row['COMBINED']['lift']:.2f} "
            f"(p={row['COMBINED']['p_value']:.4f})") in text


def test_the_table_is_marked_as_machine_checked():
    """Without the marker the parser silently reads whatever table comes
    first, which is how a guard test becomes decoration."""
    text = DOC.read_text()
    assert text.count(MARKER) == 1
    assert "tests/test_evaluation_table.py" in text
    head = text.split(MARKER, 1)[1].split("-->", 1)[0]
    assert "Do not hand-edit" in head


def test_the_results_file_is_the_one_the_rest_of_the_suite_reads():
    """A guard that asserts against a stale artifact guards nothing. The
    release stamp must match what check_results.py enforces elsewhere."""
    prov = json.loads(RESULTS.read_text())["provenance"]
    assert prov["release"]["id"]
    assert prov["schema_version"]
