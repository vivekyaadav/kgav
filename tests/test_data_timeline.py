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
    years, by_source = dt.covered(rel)
    assert years == {"2019": 2, "2021": 1}
    assert set(by_source) == {"infores:a", "infores:b"}
    assert by_source["infores:a"] == {"2019": 1, "2021": 1}


def test_a_merged_edge_is_counted_under_both_sources(tmp_path):
    """assemble writes "infores:string|infores:virhostnet" for one fact two
    sources assert. Counting that as one unnamed source would hide both."""
    rel = _release(tmp_path, [_edge("infores:string|infores:virhostnet",
                                    "2020-02-02")])
    _years, by_source = dt.covered(rel)
    assert set(by_source) == {"infores:string", "infores:virhostnet"}


def test_an_undated_edge_is_counted_as_undated_not_dropped(tmp_path):
    """19,220 censored values and a temporal split at 2021 both turn on
    what has no date. Silently excluding them would make the covered range
    look cleaner than the evidence is."""
    rel = _release(tmp_path, [_edge("infores:a", None),
                              _edge("infores:a", ""),
                              _edge("infores:a", "2020-01-01")])
    years, _ = dt.covered(rel)
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
    years, _ = dt.covered(rel)
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
