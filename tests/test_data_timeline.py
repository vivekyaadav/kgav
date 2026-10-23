"""The two timelines, kept apart.

A file downloaded yesterday can carry evidence from 1998. Reporting only
the download date says this is a current graph; reporting only the edge
dates says nothing about whether a newer release exists. The script reports
both and labels which is which.
"""
from __future__ import annotations

import io
import json
import sys
from collections import Counter
from contextlib import redirect_stdout
from importlib import import_module
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

dt = import_module("data_timeline")


def _edge(src, date):
    return {"subject": "a", "predicate": "TARGETS", "object": "b",
            "primary_knowledge_source": src, "first_asserted_date": date}


def _release(tmp: Path, edges: list[dict]) -> Path:
    d = tmp / "data" / "releases" / "v0.1"
    d.mkdir(parents=True)
    (d / "edges.jsonl").write_text(
        "\n".join(json.dumps(e) for e in edges))
    (d / "MANIFEST.json").write_text("{}")
    return d


def _run(*argv) -> str:
    buf = io.StringIO()
    with redirect_stdout(buf):
        dt.main(list(argv))
    return buf.getvalue()


# --------------------------------------------------------------- acquired
def test_the_input_list_comes_from_the_pipeline_not_a_second_copy():
    """A layer added to the rebuild driver must appear here without a
    second edit, and cannot be quietly missing from one of them."""
    from kgav.pipeline import LAYERS
    rows = dt.acquired(ROOT)
    assert {r["layer"] for r in rows} == {x.name for x in LAYERS}


def test_an_absent_input_is_listed_rather_than_omitted(tmp_path):
    """Dropping it would make a missing source invisible in the one table
    a reader checks for what the graph is built from."""
    rows = dt.acquired(tmp_path)
    assert rows and all(r["state"] == "absent" for r in rows)
    assert {r["layer"] for r in rows}


def test_a_date_with_no_stated_version_is_labelled_mtime_only(tmp_path):
    """mtime is when the file arrived here, or when something touched it.
    Presenting it as a release date would be a claim the file cannot
    support."""
    d = tmp_path / "data/raw/vhppi"
    d.mkdir(parents=True)
    (d / "vhn_1.tab27").write_bytes(b"x")
    row = next(r for r in dt.acquired(tmp_path) if r["layer"] == "vhppi")
    assert row["provenance"] == "mtime only"
    assert row["version"] == "-"


def test_a_stated_release_is_read_from_the_path(tmp_path):
    """ChEMBL numbers its releases and puts the number in the directory, so
    that is a fact rather than an mtime."""
    d = tmp_path / "data/raw/chembl/chembl_37/chembl_37_sqlite"
    d.mkdir(parents=True)
    (d / "chembl_37.db").write_bytes(b"x")
    row = next(r for r in dt.acquired(tmp_path) if r["layer"] == "chembl")
    assert row["version"] == "chembl_37"


def test_a_sidecar_is_preferred_over_the_filename(tmp_path):
    """filter_intact and fetch_raw record what a file was derived from."""
    d = tmp_path / "data/raw/gtop"
    d.mkdir(parents=True)
    f = d / "interactions.csv"
    f.write_bytes(b"x")
    (d / "interactions.csv.source.json").write_text(json.dumps(
        {"derived_from": "intact.zip", "written": "2026-10-23T00:00:00Z"}))
    row = next(r for r in dt.acquired(tmp_path) if r["layer"] == "gtop")
    assert row["provenance"] == "sidecar"
    assert row["version"] == "intact.zip"


# --------------------------------------------------------------- covered
def test_edge_years_come_from_the_assembled_release(tmp_path):
    """Not from any layer's own report: assembly merges one fact from two
    sources and keeps the EARLIEST date, so a per-layer count would
    double-count and date it late."""
    rel = _release(tmp_path, [_edge("infores:a", "2019-05-05"),
                              _edge("infores:a", "2021-01-01"),
                              _edge("infores:b", "2019-01-01")])
    years, by_source, _dates, _c = dt.covered(rel)
    assert years == {"2019": 2, "2021": 1}
    assert set(by_source) == {"infores:a", "infores:b"}
    assert by_source["infores:a"] == {"2019": 1, "2021": 1}


def test_a_merged_edge_is_counted_under_both_sources(tmp_path):
    """assemble writes "infores:string|infores:virhostnet" for one fact two
    sources assert. Counting that as one unnamed source would hide both."""
    rel = _release(tmp_path, [_edge("infores:string|infores:virhostnet",
                                    "2020-02-02")])
    _years, by_source, _dates, _c = dt.covered(rel)
    assert set(by_source) == {"infores:string", "infores:virhostnet"}


def test_an_undated_edge_is_counted_as_undated_not_dropped(tmp_path):
    """19,220 censored values and a temporal split at 2021 both turn on
    what has no date. Silently excluding them would make the covered range
    look cleaner than the evidence is."""
    rel = _release(tmp_path, [_edge("infores:a", None),
                              _edge("infores:a", ""),
                              _edge("infores:a", "2020-01-01")])
    years, _by, _dates, _c = dt.covered(rel)
    assert years["undated"] == 2


def test_the_undated_share_is_reported_with_the_split_policy(tmp_path):
    rel = _release(tmp_path, [_edge("infores:a", None),
                              _edge("infores:a", "2020-01-01")])
    out = _run("--root", str(tmp_path), "--release", str(rel), "--evidence")
    assert "carry NO usable date" in out
    assert "--undated policy" in out


def test_a_malformed_date_is_undated_rather_than_a_year(tmp_path):
    rel = _release(tmp_path, [_edge("infores:a", "n/a"),
                              _edge("infores:a", "20x1-01-01")])
    years, _by, _dates, _c = dt.covered(rel)
    assert years == {"undated": 2}


# --------------------------------------------------------------- the pair
def test_the_two_timelines_are_labelled_differently(tmp_path):
    """The whole point. One is when a file arrived, the other is when the
    claim entered the literature."""
    rel = _release(tmp_path, [_edge("infores:a", "2019-05-05")])
    out = _run("--root", str(tmp_path), "--release", str(rel))
    assert "ACQUIRED" in out and "COVERED" in out
    assert "not when the source published" in out
    assert "not the date it was downloaded" in out


def test_evidence_without_a_release_says_so_rather_than_printing_nothing(
        tmp_path):
    out = _run("--root", str(tmp_path), "--release", str(tmp_path / "none"),
               "--acquired")
    assert "ACQUIRED" in out


# ------------------- the placeholder, and the blanket release stamp
def test_the_1970_placeholder_is_undated_not_a_year(tmp_path):
    """IT READ AS 6,287 EDGES OF 1970 EVIDENCE. Every ingest writes
    1970-01-01 when it has no date -- src/kgav/temporal.py names it
    PLACEHOLDER_DATE and edge_year returns None for it -- so reading the
    year off the string invented an oldest year and reported `undated` as
    0 for every source."""
    from kgav.temporal import PLACEHOLDER_DATE
    rel = _release(tmp_path, [_edge("infores:a", PLACEHOLDER_DATE),
                              _edge("infores:a", "2020-01-01")])
    years, by_source, _d, _c = dt.covered(rel)
    assert "1970" not in years
    assert years["undated"] == 1
    assert by_source["infores:a"]["undated"] == 1


def test_the_undated_definition_is_temporals_not_a_second_one():
    """assemble._earliest special-cases the same string, and its comment is
    about exactly this duplication. A second answer here would drift."""
    from kgav.temporal import edge_year
    src = (ROOT / "scripts" / "data_timeline.py").read_text()
    assert "from kgav.temporal import" in src
    assert "edge_year(e)" in src
    assert edge_year({"first_asserted_date": "1970-01-01"}) is None


def test_a_single_edge_source_is_not_called_a_one_date_source(tmp_path):
    """One edge trivially has one date. Calling that a stamp reported that
    chembl stamps its edges, on a fixture holding a single chembl row."""
    rel = _release(tmp_path, [_edge("infores:chembl", "2001-01-01")])
    out = _run("--root", str(tmp_path), "--release", str(rel), "--evidence")
    assert "ONE DATE" not in out


def test_a_source_stamping_one_date_on_every_edge_is_marked(tmp_path):
    """ingest_host writes string_release_date to all 141,019 STRING edges
    and ingest_vhppi defaults all of its own to 2024-01-01. That is the
    FILE's release date, not when a fact was established, and drawn in the
    histogram it towers over every real year as a discovery burst."""
    rel = _release(tmp_path, [_edge("infores:string", "2023-08-28")] * 5
                   + [_edge("infores:chembl", "2001-01-01"),
                      _edge("infores:chembl", "2004-01-01")])
    out = _run("--root", str(tmp_path), "--release", str(rel), "--evidence")
    assert "ONE DATE" in out
    assert "all dated 2023-08-28" in out
    assert "the file's release" in out


def test_a_stamped_year_is_kept_out_of_the_histogram(tmp_path):
    rel = _release(tmp_path, [_edge("infores:string", "2023-08-28")] * 5
                   + [_edge("infores:chembl", "2001-01-01"),
                      _edge("infores:chembl", "2004-01-01")])
    out = _run("--root", str(tmp_path), "--release", str(rel), "--evidence")
    # just the bars: the by-source table below legitimately shows 2023
    # for the ONE DATE row, which is the point of that column.
    hist = out[out.index("per-edge dates only"):out.index("by source")]
    assert "2001" in hist
    assert "2023" not in hist, hist


def test_a_source_with_real_spread_is_not_called_a_stamp(tmp_path):
    rel = _release(tmp_path, [_edge("infores:chembl", "2001-01-01"),
                              _edge("infores:chembl", "2004-01-01")])
    out = _run("--root", str(tmp_path), "--release", str(rel), "--evidence")
    assert "per-edge" in out
    assert "ONE DATE" not in out


def test_a_fully_undated_source_is_called_none_not_a_stamp(tmp_path):
    """kgav-computed's 3,384 similarity edges have no date at all. Calling
    that a stamp would imply a release date exists for them."""
    from kgav.temporal import PLACEHOLDER_DATE
    rel = _release(tmp_path,
                   [_edge("infores:kgav-computed", PLACEHOLDER_DATE)] * 3)
    out = _run("--root", str(tmp_path), "--release", str(rel), "--evidence")
    assert "none" in out
    assert "ONE DATE" not in out


def test_a_one_date_source_stays_so_when_it_also_has_placeholders(tmp_path):
    """The placeholder must not count as a second distinct date. STRING
    happens to have none, but a source carrying a release stamp on its
    dated rows and 1970-01-01 on the rest would otherwise look like it
    dates its edges -- and be drawn in the histogram."""
    from kgav.temporal import PLACEHOLDER_DATE
    rel = _release(tmp_path,
                   [_edge("infores:string", "2023-08-28")] * 3
                   + [_edge("infores:string", PLACEHOLDER_DATE)]
                   + [_edge("infores:chembl", "2001-01-01"),
                      _edge("infores:chembl", "2004-01-01")])
    out = _run("--root", str(tmp_path), "--release", str(rel), "--evidence")
    assert "ONE DATE" in out
    hist = out[out.index("per-edge dates only"):out.index("by source")]
    assert "2023" not in hist, hist


def test_a_year_a_one_date_source_lands_in_is_not_itself_deleted(tmp_path):
    """THE FIX FOR THE HISTOGRAM DELETED THE MOST IMPORTANT BAR IN IT.

    Excluding every year a one-date source touched dropped 2020, 2023 and
    2024 outright -- ORCS dates its screens per edge and runs to 2024,
    ChEMBL to 2025, and 2020 is the year most SARS-CoV-2 evidence was
    published. What must be excluded is the stamped source's EDGES.
    """
    rel = _release(tmp_path,
                   [_edge("infores:string", "2023-08-28")] * 4
                   + [_edge("infores:biogrid-orcs", "2023-06-01"),
                      _edge("infores:biogrid-orcs", "2024-06-01")])
    out = _run("--root", str(tmp_path), "--release", str(rel), "--evidence")
    hist = out[out.index("per-edge dates only"):out.index("by source")]
    assert "2023" in hist, "ORCS's own 2023 edge was deleted with STRING's"
    assert "2024" in hist
    # and STRING's four are still out of the bars
    assert "        4  " not in hist


def test_the_histogram_says_how_many_edges_it_covers(tmp_path):
    """A bar chart that silently omits part of the graph is the shape this
    whole section keeps getting wrong. The count makes it auditable."""
    rel = _release(tmp_path,
                   [_edge("infores:string", "2023-08-28")] * 4
                   + [_edge("infores:chembl", "2001-01-01"),
                      _edge("infores:chembl", "2004-01-01")])
    out = _run("--root", str(tmp_path), "--release", str(rel), "--evidence")
    assert "per-edge dates only (2 edges)" in out


def test_an_edge_merged_from_a_stamped_and_a_dated_source_is_kept(tmp_path):
    """assemble writes "infores:string|infores:virhostnet" for one fact two
    sources assert. An edge whose sources include one that dates its edges
    has a real date and belongs in the bars."""
    rel = _release(tmp_path,
                   [_edge("infores:string|infores:chembl", "2001-01-01")]
                   + [_edge("infores:string", "2023-08-28")] * 3)
    _y, _b, _d, combos = dt.covered(rel)
    assert combos[("2001", ("infores:chembl", "infores:string"))] == 1
    out = _run("--root", str(tmp_path), "--release", str(rel), "--evidence")
    hist = out[out.index("per-edge dates only"):out.index("by source")]
    assert "2001" in hist


def test_a_merged_edge_is_not_double_counted_in_the_histogram(tmp_path):
    """by_source counts it once per source; the histogram must count the
    EDGE once, or a bar exceeds the number of edges in that year."""
    rel = _release(tmp_path,
                   [_edge("infores:chembl|infores:uniprot", "2001-01-01")])
    out = _run("--root", str(tmp_path), "--release", str(rel), "--evidence")
    assert "per-edge dates only (1 edges)" in out


# ------------------------------------------- the rule, tested directly
def test_an_edge_with_one_per_edge_source_is_in_the_bars():
    """ANY, not ALL. assemble keeps the EARLIEST date, so when a one-date
    source and a per-edge source assert the same fact, the per-edge source
    is what set the date."""
    combos = Counter({("2001", ("infores:chembl", "infores:string")): 7})
    assert dt.histogram(combos, {"infores:string": "2023-08-28"}) == {"2001": 7}


def test_an_edge_all_of_whose_sources_stamp_is_excluded():
    combos = Counter({("2023", ("infores:string", "infores:virhostnet")): 5})
    stamped = {"infores:string": "2023-08-28",
               "infores:virhostnet": "2024-01-01"}
    assert dt.histogram(combos, stamped) == {}


def test_years_accumulate_across_source_combinations():
    combos = Counter({("2001", ("infores:chembl",)): 3,
                      ("2001", ("infores:uniprot",)): 4,
                      ("2001", ("infores:string",)): 9})
    assert dt.histogram(combos, {"infores:string": "x"}) == {"2001": 7}
