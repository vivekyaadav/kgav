"""The architecture figures state facts about this repository.

WHY THIS FILE EXISTS. docs/figures/arch.dot and arch_a4.dot were drawn on
2026-11-03 and nothing checked them. By 2026-11-20 they carried twelve
wrong items, and `docs/agent_guide.md` embeds them under "How the pieces
fit together" — presented as the current architecture, red boxes and all:

    nodes                78,789  ->  78,649
    edges               392,907  ->  396,223
    results files            28  ->  36
    stale results dirs   names results-corrected/ and results-prechange/,
                         neither of which exists
    train2021 variants        3  ->  4
    web/                 "Manus/shadcn template, static Express server,
                         hand-typed demo graph, NOT connected to anything"
                         -- a directory deleted on 2026-11-20
    defect: 0 PPI in train2021   -> it holds 11,744
    defect: CI ruff 281 findings -> ruff is pinned and lint is green
    defect: reach caveat never loads -> closed 2026-11-07
    inactive share          69%  ->  77.8%
    unreachable labelled    71%  ->  72.9%
    degree vs channels  0.62/0.58 -> 0.570/0.540

THESE ARE NOT GENERATED, AND SHOULD NOT BE. An architecture diagram is
editorial: the grouping, the ordering and the choice of what to draw are
judgements, and a generator would flatten them. What CAN be checked is
every claim in them that is a fact about the repository, which is what
this file does.

A DEFECT DRAWN IN RED THAT IS CLOSED IS THE WORST CASE. It is the one a
reader acts on, and it is the one nobody re-reads. arch_a4.dot marks
closed entries grey rather than deleting them, because that figure is a
dated audit snapshot and deleting them would redraw history instead of
recording it.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "docs" / "figures"
ARCH = FIGURES / "arch.dot"
ARCH_A4 = FIGURES / "arch_a4.dot"
RELEASE = ROOT / "data" / "releases" / "v0.1"

needs_release = pytest.mark.skipif(
    not (RELEASE / "edges.jsonl").exists(),
    reason="no assembled release on disk")

BOTH = pytest.mark.parametrize("dot", [ARCH, ARCH_A4],
                               ids=["arch.dot", "arch_a4.dot"])


def _digits(text: str) -> set[int]:
    """Every comma-grouped figure in the file, as integers."""
    return {int(m.group(0).replace(",", ""))
            for m in re.finditer(r"\b\d{1,3}(?:,\d{3})+\b", text)}


@BOTH
def test_the_figure_exists(dot):
    assert dot.exists(), f"{dot.name} is gone"


@BOTH
@needs_release
def test_the_release_counts_are_the_releases(dot):
    """The two headline figures, which both files print in their title or
    their release cylinder. They were stale for seventeen days."""
    manifest = json.loads((RELEASE / "MANIFEST.json").read_text())
    nodes = sum(manifest["nodes_by_class"].values())
    edges = sum(1 for line in (RELEASE / "edges.jsonl").read_text().splitlines()
                if line.strip())
    text = dot.read_text()
    assert f"{nodes:,} nodes" in text, (
        f"{dot.name} does not say {nodes:,} nodes")
    assert f"{edges:,} edges" in text, (
        f"{dot.name} does not say {edges:,} edges")
    # And the superseded pair must not survive anywhere in the file.
    for stale in (78789, 392907):
        assert stale not in _digits(text), (
            f"{dot.name} still carries the superseded figure {stale:,}")


@BOTH
def test_the_results_census_is_the_directory(dot):
    """'28 files' outlived eight more being written."""
    n = len(list((ROOT / "data" / "results").glob("*.json")))
    text = dot.read_text()
    m = re.search(r"(\d+) (?:JSON files|files,)", text)
    assert m, f"{dot.name} no longer states a results count"
    assert int(m.group(1)) == n, (
        f"{dot.name} says {m.group(1)} results files; there are {n}")


@BOTH
def test_no_box_names_a_path_that_does_not_exist(dot):
    """A figure that draws a deleted component is worse than one with a
    stale number: the number is checkable by eye and the box is not.

    Only paths the figure spells with a directory separator are checked,
    so prose is not mistaken for a claim about the tree.
    """
    text = dot.read_text()
    claimed = set()
    for m in re.finditer(r"\b((?:src|scripts|docs|data|config|tests|web)/"
                         r"[A-Za-z0-9_./*-]*)", text):
        path = m.group(1).rstrip("/.,")
        # A GLOB OR A PLACEHOLDER IS NOT A CLAIM. The figure writes
        # "data/releases/v0.1-<layer>/" for the per-layer releases and
        # "docs/*.md" for the docs; neither names a file that must
        # exist. The first version of this flagged both and would have
        # been "fixed" by deleting the prose.
        if "*" in path or text[m.end():m.end() + 1] == "<":
            continue
        claimed.add(path)
    assert claimed, f"{dot.name} names no repository path at all"

    def present(rel: str) -> bool:
        # A MODULE MAY BE NAMED WITHOUT ITS EXTENSION. The evaluation
        # cluster is labelled "src/kgav/baselines, labels, temporal,
        # calibration, landscape", so the first of them arrives here
        # bare. That is a real reference to a real file.
        return (ROOT / rel).exists() or (ROOT / (rel + ".py")).exists()

    missing = sorted(p for p in claimed if not present(p))
    assert not missing, (
        f"{dot.name} draws {missing}, which this repository does not have")


def test_the_deleted_frontend_is_not_described_as_present():
    """THE BOX THAT WOULD HAVE SHIPPED IN THE PAPER.

    arch.dot described web/ as a 'Manus/shadcn template, static Express
    server, hand-typed demo graph, NOT connected to anything' — true
    until 2026-11-20 and false from the commit that deleted it. What is
    there now is one HTML file served by the agent's own server, and the
    figure said the opposite of the thing that matters: whether the
    frontend is wired to anything.
    """
    text = ARCH.read_text()
    for gone in ("Manus", "shadcn", "Express", "NOT connected to anything"):
        assert gone not in text, (
            f"arch.dot still describes the deleted frontend ({gone!r})")
    assert "chat.html" in text, "arch.dot does not draw the frontend that exists"
    assert "agent_server.py" in text, "it does not say what serves it"


@needs_release
def test_no_closed_defect_is_drawn_as_open():
    """Three of arch.dot's seven red boxes were fixed or false.

    Each is pinned to the measurement that closed it, so a regression
    puts the box back rather than leaving a quiet contradiction between
    the figure and the repository.
    """
    import collections

    text = ARCH.read_text()

    # d1: "train2021 has 0 PPI edges"
    ppi = collections.Counter()
    split = ROOT / "data" / "releases" / "v0.1-train2021" / "edges.jsonl"
    if split.exists():
        for line in split.read_text().splitlines():
            if line.strip():
                ppi[json.loads(line)["predicate"]] += 1
        assert ppi["PHYSICALLY_INTERACTS_WITH"] > 0, (
            "the 2021 training graph really has no PPI edges; put the box back")
        assert "0 PPI edges" not in text, (
            f"arch.dot says the 2021 training graph has 0 PPI edges and it "
            f"has {ppi['PHYSICALLY_INTERACTS_WITH']:,}")

    # d6: "CI runs full ruff -> 281 findings -> red"
    assert "281 findings" not in text, (
        "arch.dot still shows the unpinned-ruff defect; pyproject pins "
        "ruff and declares its rule set")
    assert 'ruff==' in (ROOT / "pyproject.toml").read_text(), (
        "ruff is no longer pinned, so that defect is open again")

    # d7: "reach caveat never loads (file passed where dir expected)"
    assert "reach caveat never loads" not in text, (
        "arch.dot still shows the reach-caveat defect, closed 2026-11-07")
    assert (ROOT / "tests" / "test_agent_cli_wiring.py").exists(), (
        "the test that guards the reach wiring is gone, so that defect is "
        "open again")


def test_the_a4_key_marks_closed_entries_rather_than_deleting_them():
    """It is a dated audit snapshot. Deleting a closed entry would redraw
    history; greying it records what was found and what became of it."""
    text = ARCH_A4.read_text()
    assert "CLOSED" in text, "the A4 key marks nothing as closed"
    assert text.count("CLOSED") >= 4, (
        "four of the eight audit entries are closed and the key should say so")
    assert "grey = closed since the audit" in text, (
        "the title does not explain the grey entries")


@BOTH
def test_the_figure_says_when_its_numbers_are_from(dot):
    """A correct count against an unnamed date is a count nobody can
    check later. Both files name the release, and the A4 names the day
    its figures were refreshed."""
    text = dot.read_text()
    assert "v0.1" in text
    assert re.search(r"2026-11-\d\d", text), (
        f"{dot.name} carries no date for its figures")


@needs_release
def test_the_open_defect_boxes_carry_the_measured_figures():
    """THE GAP A MUTATION FOUND.

    Reverting 72.9% to the superseded 71% passed the whole file: the
    checks above cover the headline counts, the paths and the closed
    boxes, and nothing looked at the figures inside the boxes that are
    still open. Those are the numbers a reader of the diagram acts on.

    Derived here from the release and the results file, not restated.
    """
    text = ARCH.read_text()

    # The share of labelled compounds with no edge but their label.
    labelled, other = set(), set()
    held = {"HAS_ANTIVIRAL_ACTIVITY_AGAINST", "MEASURED_INACTIVE_AGAINST"}
    for line in (RELEASE / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        if e["predicate"] == "HAS_ANTIVIRAL_ACTIVITY_AGAINST":
            labelled.add(e["subject"])
        if e["predicate"] not in held:
            other.add(e["subject"])
            other.add(e["object"])
    share = len(labelled - other) / len(labelled)
    assert f"{share:.1%} of labelled compounds" in text, (
        f"{share:.1%} of labelled compounds have no other edge and the "
        f"figure says otherwise")

    # The degree baseline against the strongest reasoning channel.
    hn = ROOT / "data" / "results" / "hard_negatives.json"
    if not hn.exists():
        pytest.skip("hard_negatives.json absent")
    cells = json.loads(hn.read_text())["cross_sectional"]["SARS-CoV-2"]
    degree = cells["degree"]["auc"]
    best = max(c["auc"] for k, c in cells.items()
               if k not in ("degree", "COMBINED")
               and isinstance(c, dict) and c.get("evaluable"))
    assert f"AUC {degree:.3f} vs ≤{best:.3f}" in text, (
        f"degree is {degree:.3f} and the strongest channel {best:.3f}; "
        f"the figure prints something else")
    assert degree > best, (
        "the strongest channel now beats degree, so that box is no longer "
        "the defect it describes")


def test_the_inactive_share_agrees_with_the_scope_statement():
    """One figure, one home. arch.dot's '77.8% measured INACTIVE' is the
    same quantity docs/scope_and_limits.md states and derives; if they
    disagree, one of them is stale and a reader cannot tell which."""
    scope = (ROOT / "docs" / "scope_and_limits.md").read_text()
    m = re.search(r"([\d.]+)%\)\*\*, [\d,]+ active", scope)
    assert m, "scope_and_limits.md no longer states the inactive share"
    assert f"{m.group(1)}% of labelled activity edges" in ARCH.read_text(), (
        f"scope_and_limits.md says {m.group(1)}% inactive and arch.dot "
        f"disagrees")


def test_the_a4_pdf_is_one_page_that_fits_a4():
    """docs/agent_guide.md promises "one page, fit to A4", and nothing
    checked it.

    graphviz adds about 0.2 inches of margin on top of -Gsize, so
    rendering with -Gsize=8.27,11.69 -- A4 exactly -- produces a page
    11.89in tall, which does not fit A4. A re-render during the
    2026-11-20 correction did exactly that. `make arch-figures` carries
    the flags that land inside the page.
    """
    import re as _re

    pdf = FIGURES / "kgav_architecture_a4.pdf"
    assert pdf.exists(), "the A4 PDF is gone"
    data = pdf.read_bytes()
    pages = data.count(b"/Type /Page") - data.count(b"/Type /Pages")
    assert pages == 1, f"the A4 PDF has {pages} pages; the guide says one"

    box = _re.findall(rb"/MediaBox\s*\[([^\]]+)\]", data)
    assert box, "the A4 PDF declares no MediaBox"
    n = [float(x) for x in box[0].split()]
    w, h = (n[2] - n[0]) / 72, (n[3] - n[1]) / 72
    assert w <= 8.27 and h <= 11.69, (
        f"the A4 PDF is {w:.2f} x {h:.2f} in and A4 is 8.27 x 11.69; it "
        f"would not print on one sheet without scaling")


def test_the_a4_png_is_still_print_resolution():
    """A careless re-render dropped it from 2480x3508 to 930x1726."""
    import struct

    png = FIGURES / "kgav_architecture_a4.png"
    assert png.exists()
    w, h = struct.unpack(">II", png.read_bytes()[16:24])
    assert h >= 3000, (
        f"the A4 PNG is {w}x{h}; at 300 dpi an A4 page is about 2480x3508, "
        f"so this has been re-rendered at screen resolution")
