"""run_baselines evaluates every virus from one DWPC computation.

dwpc_scores returns {metapath: {drug: {virus: score}}} -- every virus
at once -- and this script used ONE virus's slice. The Reproducing
block invoked it sixteen times, so the same traversal over 25,383
drugs ran sixteen times and fifteen-sixteenths of each result was
discarded. At ~91 seconds a virus that is roughly 24 of the block's
60 minutes.

The virus list also came from sixteen taxon ids hand-written into the
document. That list was eight until 2026-10-23, when the remaining
flaviviruses were promoted, so it has already gone stale once.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import scripts.run_baselines as RB  # noqa: E402

SMOKE = ROOT / "data" / "releases" / "smoke"


def test_the_register_is_the_source_of_the_virus_list():
    """Read from config/viruses.yaml, not typed anywhere."""
    got = RB.register_viruses(ROOT)
    doc = yaml.safe_load((ROOT / "config" / "viruses.yaml").read_text())
    assert got == [f"NCBITaxon:{v['taxon']}" for v in doc["viruses"]]
    assert len(got) == len(doc["viruses"])
    assert len(set(got)) == len(got), "the register has a duplicate taxon"


def test_a_promoted_virus_needs_no_document_edit(tmp_path):
    """The point of reading the register. A seventeenth virus must
    appear without anyone touching docs/evaluation.md."""
    doc = yaml.safe_load((ROOT / "config" / "viruses.yaml").read_text())
    doc["viruses"].append({"taxon": "999999", "name": "Test virus"})
    fake = tmp_path / "config"
    fake.mkdir()
    (fake / "viruses.yaml").write_text(yaml.safe_dump(doc))
    assert "NCBITaxon:999999" in RB.register_viruses(tmp_path)


def test_the_default_is_still_one_virus():
    """Nine other scripts and every ad-hoc invocation rely on it."""
    assert RB.DEFAULT_VIRUS == "NCBITaxon:2697049"


def _run(out: Path, *extra: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "run_baselines.py"),
         "--release", str(SMOKE), "--out", str(out),
         "--db", "/nonexistent", *extra],
        capture_output=True, text=True, timeout=900, cwd=ROOT)


def _payload(p: Path) -> dict:
    d = json.loads(p.read_text())
    d.pop("provenance", None)
    return d


@pytest.mark.slow
def test_all_viruses_gives_the_same_file_as_one_at_a_time(tmp_path):
    """THE WHOLE CLAIM. Batching must change nothing but the time it
    takes -- so the file a single invocation writes and the file the
    batch writes for the same virus are compared field by field, with
    only the provenance stamp dropped.

    Run on the tracked smoke release, so it needs no built graph.
    """
    one, many = tmp_path / "one", tmp_path / "many"
    one.mkdir()
    many.mkdir()
    v = RB.DEFAULT_VIRUS
    r1 = _run(one, "--virus", v)
    assert r1.returncode == 0, r1.stdout[-2000:] + r1.stderr[-2000:]
    r2 = _run(many, "--all-viruses")
    assert r2.returncode == 0, r2.stdout[-2000:] + r2.stderr[-2000:]

    name = f"baselines_{v.replace(':', '_')}.json"
    assert _payload(one / name) == _payload(many / name)
    assert len(list(many.glob("baselines_*.json"))) == \
        len(RB.register_viruses(ROOT))


@pytest.mark.slow
def test_repeating_the_flag_evaluates_each_one(tmp_path):
    out = tmp_path / "two"
    out.mkdir()
    a, b = RB.register_viruses(ROOT)[:2]
    r = _run(out, "--virus", a, "--virus", b)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    assert {p.name for p in out.glob("baselines_*.json")} == {
        f"baselines_{a.replace(':', '_')}.json",
        f"baselines_{b.replace(':', '_')}.json"}


@pytest.mark.slow
def test_the_two_selectors_cannot_be_combined(tmp_path):
    """--all-viruses with --virus is a request whose meaning is not
    obvious -- all of them, or that one? It errors rather than
    picking."""
    out = tmp_path / "x"
    out.mkdir()
    r = _run(out, "--all-viruses", "--virus", RB.DEFAULT_VIRUS)
    assert r.returncode != 0
    assert "mutually exclusive" in r.stderr
