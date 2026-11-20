"""docs/figures/meta.dot describes the release it names.

WHY THIS FILE EXISTS. meta.dot was hand-written on 2026-11-03 and nothing
checked it. By 2026-11-20 it carried thirteen stale figures and omitted
two predicates:

    nodes 78,789 -> 78,649, unique pairs 390,684 -> 393,927,
    SmallMolecule 25,383 -> 25,243, HAS_ANTIVIRAL 24,496 -> 27,928,
    TARGETS 17,982 -> 17,635, INHIBITS 3,943 -> 3,874,
    CHEMICALLY_SIMILAR_TO 3,384 -> 3,035, SARS-CoV-2 degree
    17,493 -> 19,498, and six more.

    absent entirely: IN_TRIAL_FOR (568) and TREATS (9)

The omissions are the part a stale-number check would never have caught,
and they are the worse half: that is the clinical axis, added 2026-11-17,
and one of the four findings v1 publishes. A figure that is redrawn by
hand acquires a predicate only if somebody remembers to draw it.

So the generator derives the TOPOLOGY as well as the figures, and these
tests hold it to that: every block in the release must be an arrow, and
every arrow must be a block.

The rendered PNG/SVG are deliberately NOT compared. Graphviz versions
differ in their output, so a byte comparison would fail on a machine with
a different graphviz and tell nobody anything. `make metagraph-figure`
regenerates both from the .dot in one command.
"""
from __future__ import annotations

import collections
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOT = ROOT / "docs" / "figures" / "meta.dot"
RELEASE = ROOT / "data" / "releases" / "v0.1"
SCRIPT = ROOT / "scripts" / "write_metagraph_figure.py"

needs_release = pytest.mark.skipif(
    not (RELEASE / "edges.jsonl").exists(),
    reason="no assembled release on disk")


@pytest.fixture(scope="module")
def blocks() -> dict:
    """{(subject class, predicate, object class): distinct (s, o) pairs}."""
    cls = {}
    for line in (RELEASE / "nodes.jsonl").read_text().splitlines():
        if line.strip():
            n = json.loads(line)
            cls[n["id"]] = n.get("class")
    out = collections.defaultdict(set)
    for line in (RELEASE / "edges.jsonl").read_text().splitlines():
        if line.strip():
            e = json.loads(line)
            out[(cls.get(e["subject"]), e["predicate"],
                 cls.get(e["object"]))].add((e["subject"], e["object"]))
    return out


@pytest.fixture(scope="module")
def arrows() -> dict:
    """{predicate: nnz as the figure prints it}."""
    text = DOT.read_text()
    found = {}
    for m in re.finditer(r'\[label="([A-Z_]+)[^"]*?· ([\d,]+) ·', text):
        found[m.group(1)] = int(m.group(2).replace(",", ""))
    return found


# ------------------------------------------------------------- currency
@needs_release
def test_the_figure_is_current():
    """`make metagraph-figure-check` is this, in the suite so a rebuild
    that moves the graph fails here rather than in a reviewer's eye."""
    r = subprocess.run([sys.executable, str(SCRIPT), "--check"],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr


@needs_release
def test_check_mode_fails_on_a_perturbed_file():
    """Proof the gate can fail, which is the only way to know it is one."""
    original = DOT.read_text()
    try:
        DOT.write_text(original.replace("25,243", "25,999"))
        r = subprocess.run([sys.executable, str(SCRIPT), "--check"],
                           cwd=ROOT, capture_output=True, text=True)
        assert r.returncode == 1, r.stdout + r.stderr
        assert "STALE" in r.stdout, r.stdout + r.stderr
    finally:
        DOT.write_text(original)
    # And restored byte-for-byte, or the next test reads a mutant.
    assert DOT.read_text() == original


def test_it_says_it_is_generated():
    head = DOT.read_text().splitlines()[0]
    assert "GENERATED FILE" in head and "Do not edit" in head


# ------------------------------------------------------------- topology
@needs_release
def test_every_block_in_the_release_is_an_arrow(blocks, arrows):
    """THE CHECK THAT WOULD HAVE CAUGHT THE REAL DEFECT.

    IN_TRIAL_FOR and TREATS were in the graph and not in the figure for
    three days. Comparing the numbers that WERE drawn could never find
    that; only asking the release what it holds can.
    """
    missing = sorted({p for (_s, p, _o) in blocks} - set(arrows))
    assert not missing, (
        f"the release holds {missing} and the figure does not draw them")


@needs_release
def test_every_arrow_is_a_block_in_the_release(blocks, arrows):
    """The other direction: a predicate that disappears must stop being
    drawn, or the figure asserts an edge type the graph no longer has."""
    extra = sorted(set(arrows) - {p for (_s, p, _o) in blocks})
    assert not extra, (
        f"the figure draws {extra} and the release holds no such edges")


@needs_release
def test_every_nnz_is_distinct_pairs_not_file_lines(blocks, arrows):
    """nnz is a matrix entry, so it is DISTINCT (subject, object).

    edges.jsonl holds more lines than pairs -- two sources can assert the
    same edge -- and the figure calls each arrow a sparse block matrix,
    so counting lines would be the wrong quantity under a label that
    says otherwise.
    """
    by_pred = collections.defaultdict(set)
    for (_s, pred, _o), pairs in blocks.items():
        by_pred[pred] |= pairs
    wrong = {p: (arrows[p], len(v)) for p, v in by_pred.items()
             if p in arrows and arrows[p] != len(v)}
    assert not wrong, f"figure vs release, (drawn, actual): {wrong}"


@needs_release
def test_the_held_out_predicates_are_drawn_as_held_out(blocks):
    """The red, dashed arrows are the label axis. A held-out predicate
    drawn like an ordinary one invites a reader to count it as evidence,
    which is the distinction the rest of this project is built on."""
    sys.path.insert(0, str(ROOT / "src"))
    from kgav.schema import load_schema

    text = DOT.read_text()
    held = set(load_schema().held_out_predicates())
    present = {p for (_s, p, _o) in blocks} & held
    assert present, "no held-out predicate is in the release"
    for pred in sorted(present):
        line = next((ln for ln in text.splitlines()
                     if f'label="{pred} ' in ln), None)
        assert line, f"{pred} is not drawn"
        assert "(HELD OUT)" in line, f"{pred} is drawn without the marking"
        assert "#c00000" in line, f"{pred} is not drawn in the held-out colour"


@needs_release
def test_the_figure_names_the_release_it_was_built_from():
    """A correct count against the wrong release is still wrong, and this
    figure goes into a paper."""
    sys.path.insert(0, str(ROOT / "src"))
    from kgav.provenance import release_fingerprint

    fp = release_fingerprint(RELEASE)
    text = DOT.read_text()
    assert fp["nodes_sha"] in text
    assert fp["edges_sha"] in text


@needs_release
def test_the_node_counts_are_the_manifests(blocks):
    """Read from MANIFEST.json rather than recounted, so this disagrees
    with the generator if either drifts from the release."""
    census = json.loads((RELEASE / "MANIFEST.json").read_text())["nodes_by_class"]
    text = DOT.read_text()
    for cls, n in census.items():
        assert f"{cls}\\n{n:,}" in text or f"{n:,}" in text, (
            f"{cls} has {n:,} nodes and the figure does not print it")
    total = sum(census.values())
    assert f"{total:,} × {total:,}" in text, (
        f"the adjacency is {total:,} square and the figure says otherwise")


@needs_release
def test_no_release_figure_is_typed_into_the_generator():
    """Same rule as every other reporting script here. The thirteen stale
    numbers this replaced were all literals in a hand-written file.

    ASKED AGAINST THE RELEASE, NOT AS A DIGIT PATTERN. The first version
    of this forbade any three-digit run outside the module docstring and
    flagged four things, none of them a figure: the headline taxon's
    CURIE (NCBITaxon:2697049) and the numbers in `annotation`'s own
    docstring, which explains the point using "17,635 TARGETS edges
    could be 17,635 drugs or 400". A rule that reports an identifier and
    an explanation as defects gets relaxed until it reports nothing.

    So: no figure this release actually reports may appear as a literal
    in the generator. Four digits or more, with the smaller classes
    named in the failure text -- a stylesheet-free Python file has few
    innocent four-digit numbers, and every census figure that matters
    here clears the floor.
    """
    census = json.loads(
        (RELEASE / "MANIFEST.json").read_text())["nodes_by_class"]
    figures = {v: k for k, v in census.items() if v >= 1000}
    figures[sum(census.values())] = "the node total"
    assert figures, "no census figure is large enough to protect"

    # COMMAS STRIPPED FIRST. The generator writes every figure with
    # {:,}, so a hand-typed one would most likely be typed that way too,
    # and "25,243" does not match a plain run-of-digits pattern. Without
    # this the check could not catch the shape of the defect it is for.
    #
    # The module docstring is removed because it is the RECORD of the
    # drift -- it quotes all thirteen stale figures beside their
    # replacements on purpose -- and comments go with it for the same
    # reason.
    import ast

    src = SCRIPT.read_text()
    tree = ast.parse(src)
    doc = ast.get_docstring(tree, clean=False) or ""
    body = src.replace(doc, "", 1)
    body = re.sub(r"#.*", "", body)
    body = body.replace(",", "").replace("_", "")

    literals = {int(m.group(0)) for m in re.finditer(r"\b\d{4,}\b", body)}
    typed = sorted(f"{figures[n]}={n:,}" for n in figures.keys() & literals)
    assert not typed, (
        f"{SCRIPT.name} contains {typed}, which this release reports. "
        f"Every figure it writes must be counted from the release.")
