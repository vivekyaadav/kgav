"""The scope table in doc2/09-limits.md, re-derived cell by cell.

Following tests/test_evaluation_table.py, whose docstring records why:
"a markdown table has nothing to fail". Six degree lifts drifted across
releases while the table sat unchanged, and Zika was never added at all.

A scope statement is worse to get wrong than a result table. It is the
part of the paper a reader uses to decide what the result applies to, and
every figure in it moves whenever the graph is rebuilt. So this derives
the whole table from the release and the register and hard-codes nothing:
a count that moves fails here until the document is corrected, and a
virus in the register fails here until the document lists it.
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "doc2" / "09-limits.md"
RELEASE = ROOT / "data" / "releases" / "v0.1"
CONFIG = ROOT / "config" / "viruses.yaml"

START = "<!-- TABLE: scope_by_virus"
#: NO END MARKER SINCE THE TABLE MOVED INTO doc2. The generated
#: document has one table under the marker and the parser stops at the
#: first line that is not a row, which the hand-written predecessor
#: could not rely on.
END = None

needs_release = pytest.mark.skipif(
    not (RELEASE / "edges.jsonl").exists(),
    reason="no assembled release on disk")


def _rows_in_doc() -> dict[str, dict]:
    """The table as the document states it, keyed by taxon."""
    src = DOC.read_text()
    assert START in src, "the scope table lost its marker"
    out = {}
    for line in src.split(START, 1)[1].splitlines():
        line = line.strip()
        if not line.startswith("|"):
            if out:
                break
            continue
        if line.startswith("|:") or line.startswith("|---"):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if len(cells) != 8 or cells[1] == "taxon":
            continue
        out[cells[1]] = {
            "name": cells[0], "family": cells[2],
            "dependency": int(cells[3].replace(",", "")),
            "restriction": int(cells[4].replace(",", "")),
            "antiviral": int(cells[5].replace(",", "")),
            "chain_source": cells[6],
            "in_overlap": cells[7].replace("*", ""),
        }
    assert out, "found the table marker but parsed no rows"
    return out


def _derive() -> dict[str, dict]:
    """The same table, computed from the release and the register."""
    cfg = yaml.safe_load(CONFIG.read_text())
    dep: dict[str, int] = defaultdict(int)
    res: dict[str, int] = defaultdict(int)
    av: dict[str, int] = defaultdict(int)
    with open(RELEASE / "edges.jsonl") as fh:
        for line in fh:
            if not line.strip():
                continue
            e = json.loads(line)
            obj = e.get("object", "")
            if not obj.startswith("NCBITaxon:"):
                continue
            t = obj.split(":", 1)[1]
            if e["predicate"] == "HOST_FACTOR_FOR":
                d = (e.get("qualifiers") or {}).get("direction")
                (dep if d == "dependency" else res)[t] += 1
            elif e["predicate"] == "HAS_ANTIVIRAL_ACTIVITY_AGAINST":
                av[t] += 1
    out = {}
    for v in cfg["viruses"]:
        t = v["taxon"]
        cs = v.get("chain_source") or {}
        out[t] = {
            "name": v["name"],
            "family": v.get("family") or "Coronaviridae",
            "dependency": dep[t], "restriction": res[t], "antiviral": av[t],
            "chain_source": ("species" if cs.get("taxon") in (None, t)
                             else f"strain {cs['taxon']}"),
            "in_overlap": "yes" if dep[t] else "no",
        }
    return out


# ------------------------------- the derivation the writer reads (P10)
MARKER = START          #: how the cell-by-cell parser finds the table
KEEP_COLUMNS = 2        #: name and taxon lead each row


def derive() -> dict[str, list[str]]:
    """The table the release implies: "name | taxon" -> the six cells.

    A SECOND DERIVATION, independent of the generator's. This table's
    header has said
    "generated" since it was written and NOTHING GENERATED IT -- the rows
    were maintained by hand against `test_each_column_matches_the_release`,
    which is the state fix-list P10 is about. Two leading columns are the
    key because the name is what a reader looks up and the taxon is what
    the register keys on; both are the document's and neither is rewritten.
    """
    return {f"{r['name']} | {t}": [
                r["family"], f"{r['dependency']:,}", f"{r['restriction']:,}",
                f"{r['antiviral']:,}", r["chain_source"],
                "yes" if r["in_overlap"] == "yes" else "**no**"]
            for t, r in _derive().items()}


# ------------------------------------------------------------- the table
@needs_release
def test_every_built_virus_is_in_the_table():
    """A virus the build contains and the scope statement omits is the
    worst failure available to this document."""
    missing = sorted(set(_derive()) - set(_rows_in_doc()))
    assert not missing, f"built but not described: {missing}"


@needs_release
def test_the_table_names_no_virus_the_build_lacks():
    extra = sorted(set(_rows_in_doc()) - set(_derive()))
    assert not extra, f"described but not built: {extra}"


@needs_release
@pytest.mark.parametrize(
    "field", ["name", "family", "dependency", "restriction", "antiviral",
              "chain_source", "in_overlap"])
def test_each_column_matches_the_release(field):
    doc, real = _rows_in_doc(), _derive()
    wrong = [f"{real[t]['name']} ({t}): doc {doc[t][field]!r}, "
             f"release {real[t][field]!r}"
             for t in sorted(set(doc) & set(real))
             if doc[t][field] != real[t][field]]
    assert not wrong, "\n  " + "\n  ".join(wrong)


@needs_release
def test_the_stated_totals_match_the_table():
    """Prose beside a table drifts from it. These three numbers appear in
    the sentence under the table and are the ones a reader quotes."""
    real = _derive()
    dep = sum(r["dependency"] for r in real.values())
    res = sum(r["restriction"] for r in real.values())
    av = sum(r["antiviral"] for r in real.values())
    src = DOC.read_text()
    for n in (dep, res, av):
        assert f"{n:,}" in src, f"{n:,} is not stated in the document"


@needs_release
def test_the_count_of_viruses_without_dependency_factors_is_stated():
    real = _derive()
    without = sum(1 for r in real.values() if not r["dependency"])
    # DERIVED, not spelled out. An earlier version asserted the literal
    # words "Eight of the sixteen", which is the hardcoded figure this
    # file exists to refuse: promote one virus and the sentence is wrong
    # while the assertion still reads as checking it.
    assert f"{without} of the {len(real)} viruses have no dependency " \
        f"factor" in DOC.read_text(), (
            f"{without} of {len(real)} viruses have no dependency factor; "
            f"the document must state that, in numerals a test can check")


# ----------------------------------------- claims that are not in the table
def test_the_document_names_the_reason_each_absent_virus_is_absent():
    """Three different reasons, and collapsing them would misdescribe the
    field: no screen in this release, no screen anywhere, or no screen
    reaching that coronavirus."""
    src = DOC.read_text()
    assert "No genome-wide screen in the literature at all" in src
    assert "JEV and TBEV" in src
    assert "No screen in this ORCS release" in src


def test_the_underpowered_rows_are_marked_as_such():
    """DENV-2 at 0.26 and West Nile at 0.00 invite a reading as depletion.
    An expectation of 3.1 observed as 0 is one Poisson draw."""
    src = DOC.read_text()
    assert "must not be read as depletion" in src
    assert "Poisson" in src


@needs_release
def test_the_direction_limit_is_stated_with_its_share():
    """DERIVED, not pinned. This asserted the literal "84.5%" and
    "12,155 of 14,381" -- the pre-GtoPdb split -- which stayed true as
    a sentence about history while the live figure moved to a different
    pair. The live pair is counted here and the historical one is left
    to the document to record as history.
    """
    import json as _json
    total = undirected = 0
    with open(RELEASE / "edges.jsonl") as fh:
        for line in fh:
            if '"TARGETS"' not in line:
                continue
            e = _json.loads(line)
            if e["predicate"] != "TARGETS":
                continue
            total += 1
            if (e.get("qualifiers") or {}).get("direction") in (
                    None, "", "unknown"):
                undirected += 1
    src = DOC.read_text()
    assert f"{undirected:,} of {total:,}" in src, (
        f"{undirected:,} of {total:,} TARGETS edges state no direction and "
        f"the document states a different pair")
    assert f"{100 * undirected / total:.1f}%" in src, (
        f"that is {100 * undirected / total:.1f}% and the document says "
        f"otherwise")


def test_the_evaluation_label_gap_is_stated():
    """One label means no independent check on any result."""
    src = DOC.read_text()
    assert "TREATS" in src and "IN_TRIAL_FOR" in src
    assert "rests on" in src


def test_the_stated_cause_of_the_label_gap_is_the_current_one():
    """A PRESENCE CHECK CANNOT FAIL ON A WRONG EXPLANATION.

    The test above asserts both predicate names appear and that the
    document says every evaluation "rests on" one label. Both stayed
    true while the sentence beneath them went stale: it attributed the
    empty labels to 13 of 16 viruses lacking a MONDO id, which
    doc2/records/data_audit_2026-11-07.md recorded as already false. The names
    were present, the claim was wrong, and nothing could fail.

    So this derives both figures from the register and the release.
    """
    import json
    from kgav import virus_register as vr

    src = DOC.read_text()
    viruses = vr.load()["viruses"]
    missing = [v for v in viruses if not v.get("mondo")]
    distinct = {v["mondo"] for v in viruses if v.get("mondo")}
    nodes = RELEASE / "nodes.jsonl"
    diseases = sum(1 for ln in nodes.open() if '"Disease"' in ln
                   and json.loads(ln)["class"] == "Disease")

    assert not missing, (
        f"{len(missing)} register entries carry mondo: null again; this "
        "section's cause has to be restated, not restored")
    assert len(distinct) == diseases, (
        f"{len(distinct)} distinct MONDO ids in the register but {diseases} "
        "Disease nodes in the release")
    assert f"{_WORDS[diseases]} `Disease` nodes" in src or \
        f"{diseases} `Disease` nodes" in src, (
        f"the document does not state {diseases} Disease nodes")
    # The superseded sentence is KEPT, in quotation marks, as the record
    # of what this section said -- the same convention the landscape
    # section uses. So the scan is over the prose with quoted material
    # removed, which is what distinguishes a stale claim from a quoted one.
    import re as _re
    unquoted = _re.sub(r'"[^"]*"', "", src, flags=_re.S)
    assert "13 of 16" not in unquoted, (
        "the superseded cause is asserted here outside quotation marks")
    assert "13 of 16" in src, (
        "the superseded sentence was deleted rather than quoted; this "
        "section's correction is then unreadable and this test is vacuous")
    assert "drug_indication" in src, (
        "the document states no current cause for the empty labels")


#: Spelled out, because the prose does.
_WORDS = {3: "Three", 10: "Ten", 16: "Sixteen"}
