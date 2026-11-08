"""The release MANIFEST is in version control (fix-list P12).

THE PATTERN COULD NEVER HAVE WORKED. `.gitignore` held

    data/releases/*
    !data/releases/**/MANIFEST.json

and git does not descend into an ignored DIRECTORY, so the negation had
nothing to apply to: `data/releases/*` ignores `data/releases/v0.1` itself.
Zero MANIFEST.json files were tracked, for however long the two lines had
sat there looking like they cooperated.

WHAT THAT COSTS. MANIFEST.json is the only artifact recording a release's
LAYER hashes, its per-class node counts, its per-predicate edge counts, its
merge census and the code that built it. docs/kg_release_census.txt is
committed and holds the top-level nodes_sha and edges_sha and nothing per
layer -- so a question like "which layer changed between these two
releases" had no answer in the repository at all, while every commit
message in this session quotes layer hashes to answer exactly that.

EVERY ASSERTION HERE GOES THROUGH `git check-ignore`, not through reading
the pattern. Reading the pattern is what produced the broken one: both
lines are individually correct and the pair is inert. A test that parses
.gitignore would have agreed with the author.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RELEASE = ROOT / "data" / "releases" / "v0.1"


def _is_repo() -> bool:
    return (ROOT / ".git").exists()


def _ignored(path: Path) -> bool:
    """What git itself says, which is the only authority on this."""
    r = subprocess.run(["git", "check-ignore", "-q", str(path)],
                       cwd=ROOT, capture_output=True)
    return r.returncode == 0


def _tracked(path: Path) -> bool:
    r = subprocess.run(["git", "ls-files", "--error-unmatch", str(path)],
                       cwd=ROOT, capture_output=True)
    return r.returncode == 0


pytestmark = pytest.mark.skipif(not _is_repo(), reason="not a git checkout")


def test_the_manifest_is_not_ignored():
    manifest = RELEASE / "MANIFEST.json"
    if not manifest.exists():
        pytest.skip("release v0.1 is not built in this checkout")
    assert not _ignored(manifest), (
        "data/releases/v0.1/MANIFEST.json is ignored again. A negation "
        "cannot reach inside an ignored directory, so the parent pattern "
        "must match one level deeper than the directory itself")


def test_the_manifest_is_actually_tracked():
    """Not ignored is not the same as committed. The point of P12 is that
    the file is IN the repository, so a reader without the build can still
    see which layers produced a release."""
    manifest = RELEASE / "MANIFEST.json"
    if not manifest.exists():
        pytest.skip("release v0.1 is not built in this checkout")
    assert _tracked(manifest), (
        "the manifest is no longer ignored but was never committed; "
        "`git add data/releases/v0.1/MANIFEST.json`")


@pytest.mark.parametrize("name", ["nodes.jsonl", "edges.jsonl"])
def test_the_big_files_are_still_ignored(name):
    """The other half, and the one a careless widening breaks: these are
    tens of megabytes and must stay out. A rule that tracked the manifest
    by tracking everything would pass the tests above."""
    path = RELEASE / name
    if not path.exists():
        pytest.skip("release v0.1 is not built in this checkout")
    assert _ignored(path), f"{name} is no longer ignored"


def test_a_layer_releases_contents_are_still_ignored():
    """Layer releases carry BUILT_BY.json rather than a MANIFEST, and
    nothing in them is meant to be tracked."""
    layer = ROOT / "data" / "releases" / "v0.1-folds"
    if not (layer / "edges.jsonl").exists():
        pytest.skip("the folds layer is not built in this checkout")
    assert _ignored(layer / "edges.jsonl")


def test_the_manifest_names_the_layers_it_is_tracked_for():
    """A tracked manifest that recorded nothing per layer would satisfy
    this file while leaving the question it exists for unanswerable."""
    import json
    manifest = RELEASE / "MANIFEST.json"
    if not manifest.exists():
        pytest.skip("release v0.1 is not built in this checkout")
    doc = json.loads(manifest.read_text())
    layers = doc.get("layers")
    assert layers, "the manifest records no layers"
    assert len(layers) > 1
    for name, layer in layers.items():
        assert layer.get("nodes_sha") and layer.get("edges_sha"), (
            f"layer {name!r} carries no content hash, so tracking the "
            f"manifest answers nothing about which layer moved")
