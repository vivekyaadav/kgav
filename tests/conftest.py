"""Shared fixtures, so the assembled release is parsed ONCE per session.

WHY. Four test files each loaded data/releases/v0.1 from scratch --
GraphTools(REAL) twice in test_agent_tools, a full nodes+edges json parse in
test_assemble, and Graph.load(REAL) in test_baselines -- all function-scoped
and uncached. That is 102 MB of edges.jsonl read and 371,117 lines parsed four
or more times per run, and it was the whole difference between a suite that
takes 8 seconds without a release on disk and 85 seconds with one.

Nothing about the code prevented this. The reads were simply repeated.

THESE FIXTURES ARE READ-ONLY BY CONTRACT. Session scope means one object is
shared by every test that asks for it, so a test that mutates a graph or a
node list changes what later tests see -- and the failure would surface as an
unrelated test failing depending on collection order. Query them; do not
modify them. A test that needs to mutate should build its own small release
with tmp_path, which is what every other test here already does.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
REAL_RELEASE = ROOT / "data" / "releases" / "v0.1"


def _have_release() -> bool:
    return (REAL_RELEASE / "nodes.jsonl").exists()


needs_release = pytest.mark.skipif(
    not _have_release(),
    reason="assembled release not built -- run scripts/assemble_release.py")


@pytest.fixture(scope="session")
def real_release() -> Path:
    if not _have_release():
        pytest.skip("assembled release not built")
    return REAL_RELEASE


@pytest.fixture(scope="session")
def real_jsonl(real_release) -> tuple[list[dict], list[dict]]:
    """(nodes, edges) parsed once. READ-ONLY -- see the module docstring."""
    def load(name):
        return [json.loads(line) for line
                in (real_release / name).read_text().splitlines() if line.strip()]
    return load("nodes.jsonl"), load("edges.jsonl")


@pytest.fixture(scope="session")
def real_graph(real_release):
    """baselines.Graph over the real release, built once. READ-ONLY."""
    from kgav.baselines import Graph
    from kgav.schema import load_schema
    s = load_schema()
    return Graph.load(real_release, skip_predicates=s.held_out_predicates(),
                      symmetry=s.symmetry())


@pytest.fixture(scope="session")
def real_tools(real_release):
    """agent.GraphTools over the real release, built once. READ-ONLY."""
    from kgav.agent.tools import GraphTools
    return GraphTools(real_release)
