"""The README, checked against the thing it describes.

It is the front door: the first file anyone outside this project reads,
and the one nothing was testing. On 2026-10-29 it described a
seven-virus Coronaviridae graph with 10 node classes, 21 edge classes
and 6 metapaths, pointed at a schema path that does not exist, carried a
four-row build table from days 3-6, said the four seasonal coronaviruses
had no MONDO term when all sixteen viruses carry one, and credited
DrugCentral as a source the graph has never ingested.

Every one of those was checkable against the repository. So they are
checked here, by the same rule the audit report follows: structural
claims are re-derived, and release-dependent counts are not restated at
all -- the README points at `make audit` instead.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

README = ROOT / "README.md"
REGISTER = ROOT / "config" / "viruses.yaml"


@pytest.fixture(scope="module")
def text() -> str:
    return README.read_text()


@pytest.fixture(scope="module")
def flat(text) -> str:
    return " ".join(text.split())


@pytest.fixture(scope="module")
def schema():
    from kgav.schema import load_schema
    return load_schema()


def test_the_schema_path_it_names_exists(flat):
    """It named schema/kg_schema.yaml. The file is kg_schema.yaml at the
    repository root, so the one instruction a new reader would follow
    led nowhere."""
    paths = re.findall(r"`([a-zA-Z0-9_./]+\.(?:yaml|md|txt))`", flat)
    assert paths, "the README names no files"
    for rel in paths:
        assert (ROOT / rel).exists(), f"README names {rel}, which does not exist"


def test_the_schema_counts_match_the_schema(flat, schema):
    assert f"{len(schema.node_classes)} node classes" in flat
    assert f"{len(schema.edge_classes)} edge classes" in flat
    assert f"{len(schema.metapaths)} metapaths" in flat


def test_every_virus_in_the_register_is_listed(flat):
    """It listed seven coronaviruses. The register holds sixteen viruses
    across two families, and nine of them were absent from the front
    page of their own repository."""
    from kgav.virus_register import family_of
    doc = yaml.safe_load(REGISTER.read_text())
    # SCOPED TO THE LIST ITSELF. Several of these viruses are named
    # again further down -- "DENV-4 and Yellow fever have no curated
    # virus-host interaction" -- so a document-wide membership check
    # passed with Yellow fever deleted from the list that is supposed
    # to enumerate them. Seventh time this exact shape has slipped
    # through in this project.
    listing = flat.split("Layout:", 1)[0]
    for v in doc["viruses"]:
        assert v["name"] in listing, (
            f"{v['name']} is in the register and not in the README's "
            f"virus list")
    fams: dict[str, int] = {}
    for v in doc["viruses"]:
        fam = family_of(v, doc) or "?"
        fams[fam] = fams.get(fam, 0) + 1
    assert f"{len(doc['viruses'])} viruses" in flat
    for fam, n in fams.items():
        assert f"{fam} ({n})" in flat, (
            f"the register holds {n} {fam} and the README says otherwise")


def test_it_does_not_restate_a_release_count(text):
    """The same rule the audit report follows, and for the same reason:
    the build table here quoted node and edge counts from day 4 of the
    project. `make audit` is the live figure.

    SCOPED TO THE DESCRIPTIVE HALF. "What the build actually taught us"
    is a record of past states and legitimately quotes the figures those
    bugs produced -- "12,756 edges claimed no measured value existed" is
    about a defect that was fixed, not a claim about the release. A rule
    that forbade those would forbid recording what happened.
    """
    head = text.split("## What the build actually taught us", 1)[0]
    bad = re.findall(r"\b\d{1,3},\d{3}\s+(?:nodes|edges)\b", head)
    assert not bad, f"the README restates a release count: {bad}"


def test_it_points_at_the_live_figures(flat):
    for cmd in ("make audit", "make kg-report", "make check-results"):
        assert cmd in flat, f"the README does not mention {cmd}"
    assert "docs/kg_audit_report.txt" in flat
    assert "docs/evaluation.md" in flat


def test_it_does_not_credit_a_source_the_graph_never_ingested(flat):
    """It said the graph "is built on DrugCentral and ChEMBL".
    DrugCentral appears in this repository only in two comments
    suggesting it as a future source. A licence note naming a database
    that was never ingested is a provenance claim that is simply
    untrue."""
    srcs = set()
    for path in list((ROOT / "src" / "kgav").rglob("*.py")) + \
            list((ROOT / "scripts").glob("*.py")):
        srcs |= set(re.findall(r"infores:([a-z0-9][a-z0-9._-]*)",
                               path.read_text()))
    assert "drugcentral" not in srcs, (
        "DrugCentral is now a real source; the README may credit it")
    assert "DrugCentral" not in flat, (
        "the README credits DrugCentral and nothing ingests it")


def test_the_mondo_claim_matches_the_register(flat):
    """It said the four seasonal coronaviruses have no MONDO term and
    are left unmapped. All sixteen carry one."""
    doc = yaml.safe_load(REGISTER.read_text())
    without = [v["name"] for v in doc["viruses"] if not v.get("mondo")]
    if without:
        assert "no MONDO disease term" in flat, (
            f"{without} carry no MONDO term and the README does not say so")
    else:
        assert "have no MONDO disease term" not in flat, (
            "every virus carries a MONDO term and the README says some "
            "do not")


def test_the_empty_label_predicates_are_named_if_they_are_empty():
    """The README states that two of three label predicates are empty.
    That is a fact about the release, so it is checked against the
    census the release generates rather than left as prose -- and if a
    later ingest populates them, this fails until the claim is removed.
    """
    import scripts.write_kg_report as W
    block = W.CENSUS.read_text()
    flat_readme = " ".join(README.read_text().split())
    if "NO ASSEMBLED RELEASE" in block:
        # NOT A SKIP. A test that reports green having found nothing to
        # look at is the defect this suite keeps finding in itself, so
        # the placeholder is asserted instead of waved through.
        assert "make kg-report" in block
        return
    m = re.search(r"(\d+) of (\d+) evaluation label", block)
    assert m, "the census no longer reports absent label predicates"
    n_absent = int(m.group(1))
    if n_absent:
        assert f"{'Two' if n_absent == 2 else n_absent} of three evaluation " \
               f"label predicates carry no edges" in flat_readme or \
               f"{n_absent} of three" in flat_readme, (
            f"{n_absent} label predicates are empty and the README does "
            f"not say so")
    else:
        assert "carry no edges" not in flat_readme
