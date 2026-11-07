"""A year is written as a year, and a placeholder is counted.

Fix-list C10. Four ingests had only a publication year and wrote
`f"{year}-01-01"`. Measured on the release before this change: 55,183
edges carried a date ending -01-01, and 2,757 of them -- UniProt entry
dates -- GENUINELY fell on 1 January. A reader could not tell a real
1 January from a fabricated one, and neither could any check.

    infores:chembl                              33,848
    infores:biogrid-orcs                         8,263
    infores:intact|infores:virhostnet            4,656
    infores:chembl|infores:guidetopharmacology   3,023
    infores:uniprot                              2,757   <- real dates
    infores:intact                               1,745
    ... nine more sources

The shorter form changes no verdict. temporal.edge_year,
labels.build_labels and agent/tools all read `date[:4]`, and
assemble._earliest compares ISO strings with min(), where "2024" sorts
before "2024-06-01" because it is a prefix -- the correct answer for
"the earliest point in that year".

The 14,888 placeholder edges are a separate count and are reported per
source by the drivers, because "undated" spread across seven sources is
seven different gaps with seven different remedies.
"""
from __future__ import annotations

import datetime as _dt
from pathlib import Path

import pytest

from kgav.labels import build_labels
from kgav.schema import _is_iso_date
from kgav.temporal import PLACEHOLDER_DATE, edge_year

ROOT = Path(__file__).resolve().parents[1]


# ---------------------------------------------------------- the validator
@pytest.mark.parametrize("value,ok", [
    ("2024", True),            # a date at YEAR precision
    ("2024-01-01", True),
    ("1970-01-01", True),
    ("0001", True),
    ("9999", True),
    ("0000", False),           # not a year
    ("20244", False),
    ("abcd", False),
    ("24", False),
    ("", False),
    (None, False),
])
def test_a_bare_year_is_a_valid_date(value, ok):
    assert _is_iso_date(value) is ok


def test_a_real_datetime_is_still_accepted():
    assert _is_iso_date("2024-03-15") is True
    assert _is_iso_date(str(_dt.date(2024, 3, 15))) is True


# ------------------------------------------------- nothing downstream moves
@pytest.mark.parametrize("date,year", [
    ("2024", 2024),
    ("2024-01-01", 2024),
    ("2024-06-30", 2024),
    (PLACEHOLDER_DATE, None),
])
def test_edge_year_reads_both_forms_identically(date, year):
    assert edge_year({"first_asserted_date": date}) == year


def test_the_earliest_merge_orders_mixed_precision_correctly():
    """assemble._earliest compares ISO strings. "2024" is a PREFIX of
    "2024-06-01", so it sorts first -- which is the right answer, because
    a year-precision date means "somewhere in 2024" and the earliest
    point in 2024 precedes June."""
    from kgav.assemble import _earliest
    assert _earliest("2024", "2024-06-01") == "2024"
    assert _earliest("2023-05-01", "2024") == "2023-05-01"
    assert _earliest("2024", "2025") == "2024"
    # and the placeholder still loses to any real date
    assert _earliest(PLACEHOLDER_DATE, "2024") == "2024"
    assert _earliest("2024", PLACEHOLDER_DATE) == "2024"


def test_a_year_precision_label_lands_in_the_same_temporal_window(tmp_path):
    """The split is what this must not move."""
    rel = tmp_path / "rel"
    rel.mkdir()
    (rel / "nodes.jsonl").write_text("")
    def edge(date):
        return {"subject": "INCHIKEY:A", "predicate":
                "HAS_ANTIVIRAL_ACTIVITY_AGAINST", "object": "NCBITaxon:1",
                "qualifiers": {"ec50_nm": 50.0, "relation": "=",
                               "assay_type": "cell_based_antiviral"},
                "primary_knowledge_source": "s", "evidence_tier": 1,
                "first_asserted_date": date}
    import json
    for date in ("2020", "2020-01-01"):
        (rel / "edges.jsonl").write_text(json.dumps(edge(date)))
        got = build_labels(rel, year_to=2021)
        assert got.positives.get("NCBITaxon:1") == {"INCHIKEY:A"}, date
        assert build_labels(rel, year_from=2022).stats["out_of_window"] == 1


# ------------------------------------------- no ingest fabricates 1 January
@pytest.mark.parametrize("module,needle", [
    ("src/kgav/chembl.py", 'str(int(r["year"]))'),
    ("src/kgav/host_targets.py", 'str(int(r["year"]))'),
    ("src/kgav/vhppi.py", 'str(int(r["pub_year"]))'),
    ("src/kgav/orcs.py", "return m.group(1) if m else PLACEHOLDER_DATE"),
])
def test_the_ingest_writes_a_year_not_a_first_of_january(module, needle):
    src = (ROOT / module).read_text()
    assert needle in src, f"{module} no longer writes a bare year"


def test_no_ingest_builds_a_date_by_appending_01_01():
    """The pattern itself, across every module in the package. A new layer
    reintroducing it is the regression this guards.

    PARSED, NOT GREPPED. A substring search flagged two DOCSTRINGS that
    quote the defect verbatim -- including this fix's own explanation of
    it -- so the guard would have forbidden describing the bug. An f-string
    is an ast.JoinedStr; a docstring mentioning one is a plain constant,
    and the parser tells them apart where a regex cannot.
    """
    import ast
    offenders = []
    for path in sorted((ROOT / "src" / "kgav").glob("*.py")):
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if not isinstance(node, ast.JoinedStr):
                continue
            text = "".join(v.value for v in node.values
                           if isinstance(v, ast.Constant)
                           and isinstance(v.value, str))
            if "-01-01" in text:
                offenders.append(f"{path.name}:{node.lineno}")
    assert not offenders, (
        f"these f-strings build a date ending -01-01, which cannot be "
        f"told from a real 1 January: {offenders}")
