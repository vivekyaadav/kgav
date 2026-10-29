"""scripts/compare_results.py, against planted differences.

The tool exists to answer "did reassembling the graph move any
figure". It was written on 2026-10-29 after the first attempt at that
comparison was run before anything had been regenerated, so it
compared a backup with its own unchanged source and reported no
differences. A check that passes because nothing happened looks
exactly like a check that passes because nothing changed, and the
whole value of this tool is in telling those apart -- so the
byte-identical case is tested first and hardest.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import scripts.compare_results as C

BASE = {
    "virus": "NCBITaxon:2697049",
    "n_positives": 121,
    "metrics": {"auc": 0.5161, "pool": 333},
    "top": {"degree": [{"drug": "A", "score": 2.0},
                       {"drug": "B", "score": 2.0}]},
    "provenance": {"generated_utc": "2026-10-28T00:00:00Z",
                   "code": {"commit": "aaaa", "dirty": False},
                   "release": {"id": "v0.1", "nodes_sha": "n1",
                               "edges_sha": "e1"}},
}


def _pair(tmp_path, mutate=None, stamp=True):
    b, a = tmp_path / "before", tmp_path / "after"
    b.mkdir()
    a.mkdir()
    doc = json.loads(json.dumps(BASE))
    (b / "r.json").write_text(json.dumps(doc, indent=1))
    after = json.loads(json.dumps(BASE))
    if stamp:
        after["provenance"]["release"]["edges_sha"] = "e2"
    if mutate:
        mutate(after)
    (a / "r.json").write_text(json.dumps(after, indent=1))
    return b, a


def test_byte_identical_is_a_problem_not_agreement(tmp_path, capsys):
    """The failure that motivated the tool. Nothing was regenerated,
    so the stamp did not move either -- and silence here would be read
    as 'the rebuild changed nothing'."""
    b, a = _pair(tmp_path, stamp=False)
    assert C.main([str(b), str(a)]) == 1
    out = capsys.readouterr().out
    assert "BYTE-IDENTICAL" in out
    assert "not regenerated" in out


def test_a_moved_stamp_alone_is_clean(tmp_path, capsys):
    b, a = _pair(tmp_path)
    assert C.main([str(b), str(a)]) == 0
    out = capsys.readouterr().out
    assert "1 file(s) differ in the stamp only" in out
    assert "NUMERIC CHANGE" not in out


@pytest.mark.parametrize("field,mutate", [
    ("a top-level int", lambda d: d.update(n_positives=122)),
    ("a nested float", lambda d: d["metrics"].update(auc=0.5161000001)),
    ("a pool count", lambda d: d["metrics"].update(pool=332)),
    ("a name in a list", lambda d: d["top"]["degree"][0].update(drug="Z")),
    ("a list length", lambda d: d["top"]["degree"].pop()),
])
def test_any_real_change_is_caught(tmp_path, capsys, field, mutate):
    b, a = _pair(tmp_path, mutate)
    assert C.main([str(b), str(a)]) == 1, field
    assert "NUMERIC CHANGE" in capsys.readouterr().out


def test_every_provenance_field_is_ignored(tmp_path, capsys):
    """All five volatile fields at once. Ignoring only some of them
    would report a clean rebuild as changed, which is the false alarm
    that looks identical to the real one."""
    def mutate(d):
        d["provenance"] = {"generated_utc": "2030-01-01T00:00:00Z",
                           "code": {"commit": "zzzz", "dirty": True},
                           "release": {"id": "v9", "nodes_sha": "n9",
                                       "edges_sha": "e9"}}
    b, a = _pair(tmp_path, mutate, stamp=False)
    assert C.main([str(b), str(a)]) == 0
    assert "NUMERIC CHANGE" not in capsys.readouterr().out


def test_a_missing_file_is_a_problem(tmp_path, capsys):
    b, a = _pair(tmp_path)
    (a / "r.json").unlink()
    assert C.main([str(b), str(a)]) == 1
    assert "MISSING AFTER" in capsys.readouterr().out


def test_comparing_nothing_is_not_a_pass(tmp_path, capsys):
    """Two empty directories agree on everything, and that must not
    read as success."""
    b, a = tmp_path / "b", tmp_path / "a"
    b.mkdir()
    a.mkdir()
    assert C.main([str(b), str(a)]) == 1
    assert "NOTHING WAS COMPARED" in capsys.readouterr().out


def test_fields_names_the_value_and_groups_it(tmp_path, capsys):
    """--fields is the diagnosis, not decoration: it separates one
    field moving everywhere from every file moving its own way."""
    b, a = _pair(tmp_path, lambda d: d["metrics"].update(pool=332))
    C.main([str(b), str(a), "--fields"])
    out = capsys.readouterr().out
    assert "metrics.pool" in out
    assert "before  333" in out
    assert "after   332" in out
    assert "FIELDS BY HOW MANY VALUES MOVED" in out


def test_a_list_index_collapses_in_the_histogram(tmp_path, capsys):
    """top.degree[0].drug and top.degree[1].drug are the same FIELD.
    Counting them separately would bury the signal the histogram
    exists to show."""
    def mutate(d):
        d["top"]["degree"][0]["drug"] = "Y"
        d["top"]["degree"][1]["drug"] = "Z"
    b, a = _pair(tmp_path, mutate)
    C.main([str(b), str(a), "--fields"])
    out = capsys.readouterr().out
    assert "2  top.degree[].drug" in out
    assert "top.degree[0].drug" not in out.split(
        "FIELDS BY HOW MANY VALUES MOVED")[1]
