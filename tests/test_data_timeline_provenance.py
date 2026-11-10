"""The provenance column must describe the inputs, not the sidecars.

WHY. `make data-timeline` printed `host  10 files  328,848,975  mtime only`
for a layer whose five files had all carried sha256 sidecars since
2026-10-25. The host layer declares the glob `data/raw/host/*`, so each
`<file>.source.json` this project writes was counted as an input; the
newest input was therefore a sidecar, written the day provenance was
recorded, and a sidecar has no sidecar of its own.

The report whose job is to say whether provenance exists was pushed from
"sidecar" to "mtime only" BY RECORDING PROVENANCE. docs/v1_release_plan.md
item A2 quotes the correct figures -- five files, 328,846,356 bytes -- and
was written before the sidecars existed; it is the report that moved.

Two defects, tested separately: a sidecar is not an input, and the verdict
is over all inputs rather than whichever happens to be newest.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import data_timeline as DT  # noqa: E402
from kgav.pipeline import LAYERS  # noqa: E402


class _Layer:
    def __init__(self, name, files=(), globs=()):
        self.name, self.files, self.globs = name, files, globs
    any_globs = ()


def _raw(tmp_path, names_and_sidecars):
    d = tmp_path / "data" / "raw" / "x"
    d.mkdir(parents=True)
    for name, with_sidecar in names_and_sidecars:
        (d / name).write_text("payload")
        if with_sidecar:
            (d / f"{name}.source.json").write_text(json.dumps(
                {"sha256": "0" * 64, "url": None, "provenance": "reconstructed"}))
    return tmp_path


def test_a_sidecar_is_not_counted_as_an_input(tmp_path):
    """THE BUG, in its smallest form: two data files, two sidecars, and a
    glob that sees all four. Before the fix inputs_for returned 4."""
    root = _raw(tmp_path, [("a.gz", True), ("b.gz", True)])
    layer = _Layer("x", globs=("data/raw/x/*",))
    got = [p.name for p in DT.inputs_for(layer, root)]
    assert got == ["a.gz", "b.gz"], got


def test_the_verdict_reads_sidecar_when_every_input_has_one(tmp_path):
    """And it must not depend on which file is newest, which is how the
    host row read `mtime only`: the sidecars were the newest paths."""
    root = _raw(tmp_path, [("a.gz", True), ("b.gz", True)])
    DT.LAYERS = [_Layer("x", globs=("data/raw/x/*",))]
    try:
        row, = DT.acquired(root)
    finally:
        DT.LAYERS = LAYERS
    assert row["provenance"] == "sidecar"
    assert row["path"] == "2 files"


def test_partial_coverage_is_reported_as_partial(tmp_path):
    """Neither "sidecar" nor "mtime only" is true of a layer where some
    inputs are stamped; both readings mislead in opposite directions."""
    root = _raw(tmp_path, [("a.gz", True), ("b.gz", False), ("c.gz", False)])
    DT.LAYERS = [_Layer("x", globs=("data/raw/x/*",))]
    try:
        row, = DT.acquired(root)
    finally:
        DT.LAYERS = LAYERS
    assert row["provenance"] == "sidecar 1/3", row["provenance"]


def test_the_host_layer_on_disk_now_reads_sidecar():
    """The real instance. Skipped rather than failed if data/raw is absent,
    because data/raw is gitignored and hand-populated."""
    import pytest
    if not (ROOT / "data/raw/host").exists():
        pytest.skip("data/raw/host not populated in this checkout")
    rows = {r["layer"]: r for r in DT.acquired(ROOT)}
    host = rows["host"]
    assert host["provenance"] == "sidecar", host
    assert host["path"] == "5 files", host
    assert host["bytes"] == 328_846_356, host["bytes"]
