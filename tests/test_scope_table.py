"""The scope table in docs/scope_and_limits.md, re-derived cell by cell.

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
DOC = ROOT / "docs" / "scope_and_limits.md"
RELEASE = ROOT / "data" / "releases" / "v0.1"
CONFIG = ROOT / "config" / "viruses.yaml"

START = "<!-- TABLE: scope_by_virus"
END = "<!-- END TABLE: scope_by_virus"

needs_release = pytest.mark.skipif(
    not (RELEASE / "edges.jsonl").exists(),
    reason="no assembled release on disk")


def _rows_in_doc() -> dict[str, dict]:
    """The table as the document states it, keyed by taxon."""
    src = DOC.read_text()
    block = src[src.index(START):src.index(END)]
    out = {}
    for line in block.splitlines():
        if not line.startswith("|") or line.startswith("|---"):
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
    assert out, "found the table markers but parsed no rows"
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


def test_the_direction_limit_is_stated_with_its_share():
    src = DOC.read_text()
    assert "84.5%" in src and "12,155 of 14,381" in src


def test_the_evaluation_label_gap_is_stated():
    """One label means no independent check on any result."""
    src = DOC.read_text()
    assert "TREATS" in src and "IN_TRIAL_FOR" in src
    assert "rests on" in src
