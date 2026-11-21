"""The figures in doc2/plan.md, re-derived (fix-list D4).

THAT DOCUMENT ASSERTED A FINDING THE RESULTS DO NOT SUPPORT. It read
"drugged human proteins overlap CRISPR dependency factors at 0.86-0.94 times
the rate a degree-matched random protein set does, in every one of the eight
taxa ... with matched p from 0.73 to 1.00". On the release of 2026-11-09 the
ratios run 0.00 to 1.36, TWO taxa are above 1.0 -- HCoV-NL63 at 1.36 and
SARS-CoV-2 at 1.07 are enriched, not depleted -- and matched p starts at
0.071. Its node and edge counts were four rebuilds stale (78,782 / 386,192
against 78,789 / 392,624).

It is the document the paper is written from, and it was the last results
table in this repository with nothing re-deriving it. Every other one has a
test file; the reason this one did not is that it lives in a planning
document rather than in `evaluation.md`, which is not a reason.

The conclusion survives and weakens: SARS-CoV-2's naive p of 3.96e-04 still
disappears under degree matching, so the control works. What cannot be said
is that the overlap is DEPLETED -- it is indistinguishable from the matched
background in both directions.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

DOC = ROOT / "doc2" / "plan.md"
RESULTS = ROOT / "data" / "results" / "pathway_landscape.json"
RELEASE = ROOT / "data" / "releases" / "v0.1"
MARKER = "<!-- TABLE: landscape_by_taxon"
HEADER = ["virus", "dependency factors", "observed", "matched expectation",
          "ratio", "p naive", "p matched"]


def fmt_p(p: float) -> str:
    """The document's rule for a p-value, owned here and nowhere else."""
    return f"{p:.4g}" if p >= 1e-4 else f"{p:.2e}"


def virus_labels() -> dict[str, str]:
    """taxon id -> the label config/viruses.yaml gives it."""
    import yaml
    cfg = yaml.safe_load((ROOT / "config" / "viruses.yaml").read_text())
    out: dict[str, str] = {}

    def harvest(o):
        if isinstance(o, dict):
            taxon, name = o.get("taxon"), o.get("label") or o.get("name")
            if taxon and name:
                out[f"NCBITaxon:{str(taxon).split(':')[-1]}"] = name
            for v in o.values():
                harvest(v)
        elif isinstance(o, list):
            for v in o:
                harvest(v)

    harvest(cfg["viruses"])
    return out


def derive() -> dict[str, list[str]]:
    """virus label -> the rendered cells the results file implies."""
    taxa = json.loads(RESULTS.read_text())["taxa"]
    labels = virus_labels()
    rows = {}
    for taxon, v in taxa.items():
        if not isinstance(v, dict) or "expected" not in v:
            continue
        # THE RATIO IS DERIVED THE WAY landscape.TaxonOverlap DERIVES IT.
        # A stored ratio is a third number that can disagree with the two
        # it comes from, which is why that class refuses to store one.
        ratio = (v["n_observed"] / v["expected"] if v["expected"]
                 else (float("inf") if v["n_observed"] else 0.0))
        rows[labels.get(taxon, taxon)] = [
            f"{v['n_dependency']:,}",
            str(v["n_observed"]),
            f"{v['expected']:.1f}",
            f"{ratio:.2f}",
            fmt_p(v["p_naive"]),
            f"{v['p_matched']:.4g}",
        ]
    return rows


def parse() -> dict[str, list[str]]:
    text = DOC.read_text()
    assert MARKER in text, "the landscape table lost its machine-readable marker"
    lines = [ln.strip() for ln in text.split(MARKER, 1)[1].splitlines()]
    start = next(i for i, ln in enumerate(lines) if ln.startswith("| virus |"))
    header = [c.strip() for c in lines[start].strip("|").split("|")]
    assert header == HEADER, header
    rows = {}
    for line in lines[start + 2:]:
        if not line.startswith("|"):
            break
        cells = [c.strip().strip("*") for c in line.strip("|").split("|")]
        assert len(cells) == 7, f"expected 7 columns: {line}"
        assert cells[0] not in rows, f"{cells[0]} appears twice"
        rows[cells[0]] = cells[1:]
    return rows


@pytest.fixture(scope="module")
def tables():
    if not RESULTS.exists():
        pytest.skip("pathway_landscape.json is not built in this checkout")
    return parse(), derive()


def test_the_table_covers_exactly_the_taxa_in_the_results(tables):
    doc, got = tables
    assert set(doc) == set(got), (
        f"only in the document: {sorted(set(doc) - set(got))}; "
        f"only in the results: {sorted(set(got) - set(doc))}")


@pytest.mark.parametrize("column", range(6))
def test_every_cell_matches_the_results_file(tables, column):
    doc, got = tables
    labels = ["dependency factors", "observed", "matched expectation",
              "ratio", "p naive", "p matched"]
    bad = {v: (doc[v][column], got[v][column])
           for v in sorted(set(doc) & set(got))
           if doc[v][column] != got[v][column]}
    assert not bad, (
        f"{labels[column]} is stale -- document vs re-derived: " +
        "; ".join(f"{v}: {d!r} != {g!r}" for v, (d, g) in bad.items()))


def test_the_depletion_claim_is_not_reinstated(tables):
    """The sentence this item is about.

    While any taxon sits above 1.0 the document may not say the overlap
    is depleted everywhere, in those words or any others that assert a
    one-sided band. Checked against the derived ratios rather than
    against a date, so reinstating the claim fails on data.
    """
    _doc, got = tables
    ratios = {v: float(cells[3]) for v, cells in got.items()}
    enriched = {v: r for v, r in ratios.items() if r > 1.0}
    flat = " ".join(DOC.read_text().split())
    # A WITHDRAWN CLAIM MAY BE QUOTED; IT MAY NOT BE ASSERTED. The
    # correction paragraph quotes the old sentence verbatim in order to
    # retract it, so a flat substring search calls the retraction the
    # defect -- which it did on the first run of this file, exactly as it
    # did in test_temporal_table.py for "degree is at chance". Quoted
    # spans are stripped, and the retraction is required to survive that.
    unquoted = re.sub(r'"[^"]*"', " ", flat)
    if enriched:
        assert "0.86–0.94 times the rate" not in unquoted, (
            f"{sorted(enriched)} have a ratio above 1.0 and the document "
            f"asserts a 0.86-0.94 band in every taxon outside a quotation")
        assert "0.86–0.94 times the rate" in flat, (
            "the retraction of that claim has been deleted rather than "
            "the claim corrected; the document should record what it "
            "used to assert")
        assert "ENRICHED, not depleted" in flat, (
            f"{sorted(enriched)} are enriched and the document does not "
            f"say so")
    lo, hi = min(ratios.values()), max(ratios.values())
    assert f"{lo:.2f} to {hi:.2f}" in flat, (
        f"the ratios run {lo:.2f} to {hi:.2f} and the document states a "
        f"different range")


def test_the_control_still_works_and_the_document_says_which_way(tables):
    """The surviving half of the finding, asserted so it is not lost in
    the correction: SARS-CoV-2's naive p is significant and its matched p
    is not. If a release ever makes the matched p significant too, the
    degree-matched control has stopped controlling and the paragraph has
    to be rewritten rather than kept."""
    taxa = json.loads(RESULTS.read_text())["taxa"]
    sc2 = taxa["NCBITaxon:2697049"]
    assert sc2["p_naive"] < 0.05, "the naive p is no longer significant"
    assert sc2["p_matched"] > 0.05, (
        "the matched p is now significant too; degree matching no longer "
        "removes this signal and the release plan says it does")
    flat = " ".join(DOC.read_text().split())
    assert fmt_p(sc2["p_naive"]) in flat
    assert f"{sc2['p_matched']:.4g}" in flat


def test_no_superseded_landscape_p_is_still_asserted_in_prose(tables):
    """THE HALF THE TEST ABOVE COULD NOT FAIL ON, found 2026-11-17.

    It asserts the CURRENT figures are present, and they were -- in the
    machine-checked table. Two sentences underneath went on asserting the
    PREVIOUS run's: "matched p runs from **0.05964**" and "SARS-CoV-2's
    naive p of 3.96e-04", each twice, after the L27 rebuild moved them to
    0.05909 and 0.0003728. A presence check passes while a stale claim
    sits beside the live one, which is this repository's own shape: a
    check that cannot fail in the way that matters.

    SCOPED TO THE SECTION THAT DISCUSSES THE TABLE, because a bare decimal
    elsewhere in the document is not a p-value -- the first draft of this
    test flagged "+0.0050", which is an AUC delta in the scorer note.

    A superseded figure belongs inside quotation marks -- or after the
    word "not", which is the same statement in another form -- recorded as
    what was claimed; that is the convention tests/test_release_plan_temporal_
    note.py enforces, and quoted material is stripped here for the same
    reason. The cost is the same and is stated there: a live claim written
    inside quotes would be invisible to this test.
    """
    import re

    text = DOC.read_text()
    start = text.index("### The finding the paper is about")
    section = text[start:text.index("## A. Blocking", start)]
    # Quoted material is what was CLAIMED, and so is the "is X, not Y"
    # correction form -- "the naive p is 0.0003728, not 1.28e-04" marks
    # 1.28e-04 as superseded exactly as a quotation would. Both are
    # removed before the scan. The cost is stated in the sibling test and
    # is the same here: a live claim written inside quotes, or immediately
    # after the word "not", would be invisible to this check.
    scrubbed = re.sub(r'"[^"]*"', " ", section)
    scrubbed = re.sub(r"\bnot\s+\*{0,2}[\d.]+(?:e-\d+)?\*{0,2}", " ", scrubbed)
    unquoted = " ".join(scrubbed.split())

    taxa = json.loads(RESULTS.read_text())["taxa"]
    live = {fmt_p(r["p_naive"]) for r in taxa.values()}
    live |= {f"{r['p_matched']:.4g}" for r in taxa.values()}
    live |= {"1", "0.73", "1.00"}          # the quoted band, and p = 1 exactly

    found = re.findall(r"\b\d\.\d{2,3}e-0\d\b|\b0\.0\d{3,5}\b", unquoted)
    stale = sorted({f for f in found if f not in live})
    assert not stale, (
        f"these p-values are asserted in the landscape section and match "
        f"no taxon in {RESULTS.name}: {stale}. A superseded figure goes "
        f"inside quotes, recorded as what was claimed")


def test_the_node_and_edge_counts_match_the_release():
    """Four rebuilds stale, in the first line of the document."""
    manifest = RELEASE / "MANIFEST.json"
    if not manifest.exists():
        pytest.skip("release v0.1 is not built in this checkout")
    m = json.loads(manifest.read_text())
    text = DOC.read_text()
    want = f"{m['nodes']:,} nodes, {m['edges']:,} edges"
    assert want in text, (
        f"the release holds {want} and the document's opening line states "
        f"something else: "
        f"{re.search(r'Release v0.1:.*', text).group(0) if re.search(r'Release v0.1:.*', text) else '(no such line)'}")


def test_the_scorer_row_counts_the_cells_the_runs_hold():
    """THE ROW SAID "10 of 15" AND THE RUNS HELD 13 OF 17.

    This document's re-derived census table exists because its own prose
    went stale; the row recording the scorer's state then went stale in
    the same way, in the same table, because it was typed. The entry-assay
    ingest added two plate cells and the degree-comparator fix moved which
    cells clear, and nothing could fail on either.

    So the count is derived. Both halves of it: how many cells clear their
    own registered alpha, and how many of those do it above a degree
    baseline that does not clear -- which is the claim
    doc2/08-results.md's claim block makes, and the only one of the two
    that is evidence about chemistry.
    """
    import json

    sym = ROOT / "data" / "results" / "structure_scorer_symmetric.json"
    if not sym.exists():
        pytest.skip("the symmetric scorer run is absent")
    doc = json.loads(sym.read_text())
    total = clearing = strong = 0
    for pop in ("primary_chembl", "secondary_ncats", "secondary_ncats_entry"):
        for cell in (doc.get(pop) or {}).values():
            if not isinstance(cell, dict):
                continue
            s, g = cell.get("structure_max_k10"), cell.get("degree")
            if not (isinstance(s, dict) and "auc" in s
                    and isinstance(g, dict) and "auc" in g):
                continue
            total += 1
            if not s.get("clears_corrected_alpha"):
                continue
            clearing += 1
            if s["auc"] > g["auc"] and not g.get("clears_corrected_alpha"):
                strong += 1

    text = " ".join(DOC.read_text().split())
    assert f"**{clearing} of {total}** cells clearing their alpha" in text, (
        f"{clearing} of {total} cells clear and the plan says otherwise")
    assert f"{strong} of them above a non-clearing degree baseline" in text, (
        f"{strong} cells clear above a non-clearing baseline and the plan "
        f"says otherwise")
