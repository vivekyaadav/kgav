"""What scripts/kgav_ask.py actually hands its dependencies.

THE BUG THIS FILE EXISTS FOR. labels.measured_reach takes the results
DIRECTORY and reads metapath_comparison.json from it. --results names a
results FILE -- hard_negatives.json, the calibration run. kgav_ask passed
--results straight through, so measured_reach looked for

    data/results/hard_negatives.json/metapath_comparison.json

which is not a path that can exist. It returned None, the agent printed
"reach is unmeasured", and every question ever asked of the CLI carried
that caveat -- including for the five viruses whose reach IS measured and
IS zero, which is the distinction measured_reach was written to make.

tests/test_labels.py covers measured_reach thoroughly and all of it passed,
because every one of those tests passes a directory: the author of the test
knew what the function wanted. Nothing asserted what the CLI passes. A unit
test asks the question its author thought to ask, so this file asks the
other one -- about the defaults the parser declares, not about a copy of
them kept here.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

from kgav.labels import measured_reach

ROOT = Path(__file__).resolve().parents[1]


def _kgav_ask():
    """The script as a module. It is not importable as a package member."""
    spec = importlib.util.spec_from_file_location(
        "kgav_ask_under_test", ROOT / "scripts" / "kgav_ask.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _defaults(root: Path):
    """The parsed namespace for a bare invocation, with root substituted."""
    return _kgav_ask().build_parser(root).parse_args(["why"])


def test_results_defaults_to_a_file_not_a_directory(tmp_path):
    """The premise of the bug, asserted so the test below cannot go green
    by --results quietly becoming a directory."""
    args = _defaults(tmp_path)
    assert args.results.name.endswith(".json")
    assert args.results.parent.name == "results"


def test_the_reach_lookup_is_given_the_directory(tmp_path):
    """The fix, end to end: build a release and a results directory, point
    the parser's own defaults at them, and require the reach to come back.

    Before the fix this returned None no matter what was on disk."""
    rel = tmp_path / "data" / "releases" / "v0.1"
    rel.mkdir(parents=True)
    (rel / "nodes.jsonl").write_text('{"id": "NCBITaxon:1"}\n')
    (rel / "edges.jsonl").write_text(
        '{"subject": "a", "predicate": "P", "object": "b"}\n')

    from kgav.provenance import write_results
    res = tmp_path / "data" / "results"
    res.mkdir(parents=True)
    write_results(res / "metapath_comparison.json",
                  {"per_virus": {"DENV-1": {"COMBINED": {"pos": 0}},
                                 "DENV-2": {"COMBINED": {"pos": 117}}}},
                  rel)

    args = _defaults(tmp_path)
    assert args.release == rel, "the release default moved; fix this test"
    # EXACTLY WHAT main() DOES. Spelled the same way on purpose: if the
    # call site stops deriving the directory, this stops finding the file.
    assert measured_reach(args.results.parent, args.release) == {
        "DENV-1": 0, "DENV-2": 117}


def test_handing_measured_reach_a_file_is_an_error_not_an_absence(tmp_path):
    """None means "cannot be known". A caller that passed the wrong kind of
    path has not established that, and must not be able to say it has."""
    rel = tmp_path / "rel"
    rel.mkdir()
    (rel / "nodes.jsonl").write_text('{"id": "NCBITaxon:1"}\n')
    (rel / "edges.jsonl").write_text(
        '{"subject": "a", "predicate": "P", "object": "b"}\n')
    f = tmp_path / "hard_negatives.json"
    f.write_text("{}")
    with pytest.raises(ValueError, match="DIRECTORY"):
        measured_reach(f, rel)


def test_an_absent_directory_is_still_an_absence(tmp_path):
    """The guard above must not turn a legitimately missing results
    directory into a crash -- that is the state on a fresh clone."""
    rel = tmp_path / "rel"
    rel.mkdir()
    (rel / "nodes.jsonl").write_text('{"id": "NCBITaxon:1"}\n')
    (rel / "edges.jsonl").write_text(
        '{"subject": "a", "predicate": "P", "object": "b"}\n')
    assert measured_reach(tmp_path / "no-results-here", rel) is None
