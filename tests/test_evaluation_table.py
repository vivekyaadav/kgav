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


WORD = {0: "No", 1: "One", 2: "Two", 3: "Three", 4: "Four", 5: "Five",
        6: "Six", 7: "Seven", 8: "Eight", 9: "Nine", 10: "Ten",
        11: "Eleven", 12: "Twelve", 13: "Thirteen", 14: "Fourteen",
        15: "Fifteen", 16: "Sixteen"}


def test_the_unreadable_count_is_derived_not_typed(tables):
    """"Four of seven" survived two releases as prose while the table
    under it changed. The headline is the count, so the count is
    derived."""
    doc, _ = tables
    n = sum(1 for cells in doc.values() if cells[0] == DASH)
    text = " ".join(DOC.read_text().split())
    assert (f"{WORD[n]} of {WORD[len(doc)].lower()} viruses have no "
            f"readable channel") in text, (
        f"{n} of {len(doc)} rows name no channel and the document does "
        f"not say so in those words")


def test_every_virus_with_no_reachable_positive_is_named():
    """A virus whose COMBINED pool holds no positive cannot be scored at
    all, which is a different statement from scoring badly. Naming them
    individually is the only way a reader can tell which rows of the
    table are weak evidence and which are no evidence."""
    per = json.loads(RESULTS.read_text())["per_virus"]
    text = DOC.read_text()
    unreachable = sorted(v for v, ch in per.items()
                         if ch.get("COMBINED", {}).get("pos") == 0)
    assert unreachable, "no virus has an empty positive set; prose stale"
    claim = (f"{WORD[len(unreachable)]} of {WORD[len(per)].lower()} viruses "
             f"have no positive reachable")
    # SCOPED TO THE PARAGRAPH MAKING THE CLAIM, not the whole document.
    # Every one of these viruses also has a row in the table above, so a
    # document-wide `v in text` was satisfied by the row and a name
    # deleted from this list passed unnoticed.
    paras = [" ".join(p.split()) for p in text.split("\n\n")]
    hit = [p for p in paras if claim in p]
    assert hit, f"the document does not state: {claim}"
    assert len(hit) == 1, "the claim is made in more than one paragraph"
    # SET EQUALITY ON THE ENUMERATION ITSELF, not membership in the
    # paragraph. The paragraph goes on to discuss several of these
    # viruses by name, so a name deleted from the list was still found
    # a sentence later and the mutant passed. The list runs from the
    # marker to the first full stop.
    marker = "(`COMBINED.pos == 0`):"
    assert marker in hit[0], "the list lost its machine-readable marker"
    listed_text = hit[0].split(marker, 1)[1].split(".", 1)[0]
    listed = {w.strip() for w in listed_text
              .replace(" and ", ", ").split(",") if w.strip()}
    assert listed == set(unreachable), (
        f"the enumeration is {sorted(listed)} and the results file says "
        f"{sorted(unreachable)}")
    # The two shapes are not alike and the document must not merge them.
    empty_pool = sorted(v for v, ch in per.items()
                        if ch.get("COMBINED", {}).get("pool") == 0)
    for v in empty_pool:
        assert f"{v} and" in hit[0] or f"and {v}" in hit[0] \
            or f"{v} reaches" in hit[0], (
            f"{v} reaches no compound at all, which the paragraph must "
            f"distinguish from reaching compounds with no label")


def test_every_significant_readable_cell_is_named():
    """The document names each significant readable cell with its own
    figures, and states how many there are.

    IT USED TO HARD-CODE ONE. `assert len(sig) == 1` held from the day
    it was written until 2026-10-30, when raising the path cap
    completed M4's traversal and DENV-2's M4 came in at p=0.019 --
    a second cell, and the assertion failed on the arrival of a
    RESULT rather than on a defect. The document had to change; so did
    the test, which is a rewrite that deriving would have avoided.

    at_chance_runs() three functions above carries the same lesson in
    its own docstring, learned there and not here: "Deriving it means
    the next move is caught without a rewrite, which is the property
    the document asks of its own figures."

    So the count is derived too. A third cell now fails this by
    pointing at a sentence that states the wrong number, not by
    tripping a literal.
    """
    per = json.loads(RESULTS.read_text())["per_virus"]
    readable = [(v, n, m) for v, ch in per.items() for n, m in ch.items()
                if n not in ("degree", "COMBINED")
                and m.get("hits_determined") and m.get("trustworthy")]
    sig = [(v, n, m) for v, n, m in readable
           if m.get("p_value") is not None and m["p_value"] < 0.05]
    flat = " ".join(DOC.read_text().split())
    assert f"{WORD[len(readable)].lower()} cells are both" in flat, (
        f"{len(readable)} readable cells; the document states a different "
        f"number")
    assert sig, (
        "no readable cell is significant any more; the document builds a "
        "section on at least one and must be corrected rather than this "
        "test relaxed")
    assert f"{WORD[len(sig)]} of them are nominally significant" in flat, (
        f"{len(sig)} readable cells are significant; the document states "
        f"a different number")
    for v, n, m in sorted(sig, key=lambda x: x[2]["p_value"]):
        assert (f"{v} {n}, lift {m['lift']:.2f}, "
                f"p={m['p_value']:.6f}") in flat, (
            f"{v} {n} is significant at p={m['p_value']:.6f} and the "
            f"document does not name it with its figures")


def test_the_correction_families_are_arithmetic():
    """The alpha table is four divisions. A typed alpha beside a live
    p-value is a verdict nothing re-derives."""
    flat = " ".join(DOC.read_text().split())
    for n in (13, 91, 120, 150):
        assert f"| {n} | {0.05 / n:.5f} |" in flat, (
            f"the alpha row for n={n} is not 0.05/{n} = {0.05 / n:.5f}")


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
    # SARS-CoV M1 IS NOT DETERMINED, SO NO POINT IS QUOTED FOR IT.
    # This asserted f"lift {m1['lift']:.2f} at p<0.0001", which held
    # until 2026-10-30: writing the release in a canonical order
    # changed nothing about the graph and moved the figure from 2.08
    # to 0.83, because rank() sorts on score with no tiebreak and
    # Python's sort is stable, so which of 23 tied compounds fall
    # inside Hits@100 follows the order edges.jsonl was written in.
    #
    # The test was demanding that the document publish one draw from
    # [0.21, 2.39] and then say, two clauses later, that the interval
    # contains 1.00. It now requires the RANGE and forbids the point.
    m1 = per["SARS-CoV"]["M1"]
    assert m1["hits_determined"] is False, (
        "SARS-CoV M1 is determined now; this paragraph is written for "
        "an undetermined cell and needs rewriting, not relaxing")
    para = text.split("Note SARS-CoV especially.", 1)[1].split("\n\n##", 1)[0]
    assert "**NO POINT ESTIMATE" in para, (
        "the correction paragraph is gone, so the scoping below is "
        "silently checking the whole paragraph")
    # EVERY ASSERTION SCOPED TO THE CLAIM, not the paragraph. The
    # correction beneath it also names [0.21, 2.39], so a check
    # against the paragraph passed with the interval deleted from the
    # claim -- satisfied by an unrelated occurrence a few lines down,
    # which is the defect this suite has found most often and which
    # appeared here in the test written to fix it.
    claim = para.split("**NO POINT ESTIMATE", 1)[0]
    assert f"[{m1['hits_lo']}, {m1['hits_hi']}] of {m1['pos']}" in claim
    assert f"[{m1['lift_lo']:.2f}, {m1['lift_hi']:.2f}]" in claim
    assert "contains 1.00" in claim
    # The correction QUOTES the retracted "lift 2.08" on purpose, so
    # forbidding the point outright would forbid recording the
    # mistake -- and if a later reorder lands the figure back on 2.08,
    # an unscoped check would trip on the history rather than on the
    # claim. Same reason at_chance_runs scopes its block.
    assert f"lift {m1['lift']:.2f}" not in claim, (
        f"the paragraph quotes lift {m1['lift']:.2f} as a point estimate "
        f"for a cell whose hits are undetermined across "
        f"[{m1['hits_lo']}, {m1['hits_hi']}]")


def at_chance_runs() -> dict[str, dict]:
    """Every (run, virus) whose verdict is BASELINE AT CHANCE, derived.

    The first version of this test hard-coded the run: SARS-CoV-2 under
    temporal computed_only. On the 16-virus release that run fell to NOT
    DETERMINED and the condition moved to the cross-sectional table, so
    the test failed for the right reason and then had to be rewritten
    anyway. Deriving it means the next move is caught without a rewrite,
    which is the property the document asks of its own figures.
    """
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from kgav.baselines import confound_warning

    runs = {"cross-sectional": ("metapath_comparison.json", "per_virus")}
    for policy in ("computed_only", "include", "exclude"):
        runs[f"temporal {policy}"] = (
            f"temporal_2021_{policy}.json", "results")
    out = {}
    for label, (name, key) in runs.items():
        path = ROOT / "data" / "results" / name
        if not path.exists():
            continue
        for virus, row in json.loads(path.read_text())[key].items():
            warning = confound_warning(row)
            if warning and warning.startswith("BASELINE AT CHANCE"):
                out[f"{label} / {virus}"] = row
    return out


def test_the_at_chance_claim_names_the_run_it_holds_under():
    """The document said BASELINE AT CHANCE "no longer holds anywhere".
    It does hold. The replacement sentence is pinned to whichever run
    holds it now, and to the channel carrying it -- naming the verdict
    without the channel is how the table came to print M4 beside a
    verdict about M7."""
    holding = at_chance_runs()
    assert holding, (
        "no run reads BASELINE AT CHANCE any more; the document still "
        "says one does and must be corrected rather than this test "
        "relaxed")
    text = DOC.read_text()
    # Whitespace-normalised, because the verdict block is a wrapped
    # indented paragraph: "M7's 1.34\n   (p=0.0004)" is the figure being
    # present, and a raw substring search calls it missing.
    flat = " ".join(text.split())
    # Scoped to the verdict block: the correction paragraph beneath it
    # QUOTES the retracted claim, and a test forbidding the words
    # outright would forbid recording the mistake.
    block = text.split("    CONFOUNDED", 1)[1].split("    (nothing)", 1)[0]
    assert "no longer holds" not in block
    for label, row in holding.items():
        virus = label.split("/")[1].strip()
        assert virus in block, f"{label} holds AT CHANCE and is not named"
        deg = row["degree"]
        assert (f"degree lift {deg['lift']:.2f} "
                f"(p={deg['p_value']:.4f})") in flat, label
        # The channel that carries it, named with its own figures: the
        # verdict text itself quotes one, so the document must agree.
        from kgav.baselines import confound_warning
        carried = confound_warning(row).split("so ", 1)[1].split("'s", 1)[0]
        m = row[carried]
        assert (f"{carried}'s {m['lift']:.2f} "
                f"(p={m['p_value']:.4f})") in flat, (
            f"{label} is carried by {carried} and the document does not "
            f"name it with its figures")


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
