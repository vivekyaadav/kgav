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
    # ANCHORED ON THE HEADING, NOT ON A WORD OF PROSE. This read
    # `flat.split("Layout:")[0]`, which scoped the check to everything
    # above a paragraph that happened to sit below the list. Moving that
    # paragraph -- which a restructure did -- silently widened the scope
    # to the whole file, which is the exact failure the scoping exists
    # to prevent.
    assert "## The register" in flat, (
        "the README has no `## The register` section; the virus list has "
        "nowhere to be scoped to")
    listing = flat.split("## The register", 1)[1].split("## Layout", 1)[0]
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

    # WIDENED. This forbade only `N,NNN nodes|edges`, so six AUCs, two
    # confidence intervals and four other release counts -- 568
    # IN_TRIAL_FOR, 9 TREATS, 6,765 gtop edges, 12 assembled layers --
    # sat in the descriptive half of a file whose own policy paragraph
    # says a count here goes stale silently. Three forms are forbidden
    # now, each one a shape only a measurement takes:
    #
    #   0.dddd                an AUC or a p-value
    #   [lo, hi]              a confidence interval
    #   N,NNN                 any comma-grouped count; every release
    #                         figure in this graph is at least a thousand
    #
    # SMALL BARE INTEGERS ARE STILL ALLOWED, because the structural
    # counts are not release-dependent and are checked against their
    # source by the tests above: 11 node classes, 24 edge classes and 8
    # metapaths against the frozen schema, 16 viruses and the two family
    # sizes against the register.
    forbidden = {
        "a four-decimal figure (an AUC or a p-value)":
            re.findall(r"\b\d\.\d{4}\b", head),
        "a confidence interval":
            re.findall(r"\[\s*\d\.\d+\s*,\s*\d\.\d+\s*\]", head),
        "a comma-grouped count":
            re.findall(r"\b\d{1,3},\d{3}\b", head),
    }
    bad = {why: hits for why, hits in forbidden.items() if hits}
    assert not bad, (
        "the descriptive half of the README restates a measured figure, "
        "which its own policy paragraph forbids: "
        + "; ".join(f"{why}: {hits}" for why, hits in bad.items())
        + ". Put the figure in doc2/, where it is derived, and point at it.")


def test_every_path_it_links_or_quotes_exists(text):
    """IT POINTED AT A SCHEMA PATH THAT DOES NOT EXIST, which is what
    `test_the_schema_path_it_names_exists` was written for -- one path,
    checked, while the file now carries eleven markdown links and a
    dozen backticked paths. A front door whose links 404 is worse than
    one with no links."""
    flat_text = " ".join(text.split())
    missing = []
    for target in re.findall(r"\]\(([^)]+)\)", text):
        if target.startswith(("http://", "https://", "#")):
            continue
        if not (ROOT / target).exists():
            missing.append(target)
    assert not missing, f"the README links to paths that do not exist: {missing}"

    # AND THE BACKTICKED ONES. A path in backticks is a path a reader
    # will type, and nothing was checking those at all. Only the ones
    # that look like paths into this repository are tested, and a glob
    # is resolved by its parent directory.
    #
    # ONE PATH IS EXEMPT AND IT HAS TO EARN IT: `docs/` is named here
    # precisely because it does NOT exist any more, and the exemption is
    # void unless the file says so. An allowlist that outlives its
    # reason is defect shape 3 in doc2/14-defect-ledger.md.
    retired = {"docs/": "retired on 2026-11-21 and named here as retired"}
    for name, why in retired.items():
        assert f"`{name}` was retired" in flat_text, (
            f"{name} is exempt from this check because it is {why}, and "
            f"the README no longer says it was retired")
        assert not (ROOT / name).exists(), (
            f"{name} exists again; delete its exemption here")

    for quoted in sorted(set(re.findall(r"`([A-Za-z0-9_./*-]+)`", text))):
        if "/" not in quoted or quoted in retired:
            continue
        probe = Path(quoted)
        if "*" in quoted:
            probe = probe.parent
        if not (ROOT / probe).exists():
            missing.append(quoted)
    assert not missing, (
        f"the README quotes paths that do not exist: {missing}")


def test_the_no_counts_rule_can_actually_fail(text):
    """A CHECK THAT CANNOT FAIL IS WORSE THAN NO CHECK -- defect shape 1
    in doc2/14-defect-ledger.md, and the predecessor of this rule was an
    instance of it: it matched `N,NNN nodes|edges` and the README's
    measured figures were AUCs, so it passed over all of them for weeks.

    Each forbidden shape is planted in a copy of the descriptive half and
    the rule must reject it. Nothing on disk is touched -- an earlier
    test in this repository proved a check could fail by editing a
    tracked file in place, which under `pytest -n auto` another worker
    reads mid-edit.
    """
    head = text.split("## What the build actually taught us", 1)[0]
    for planted in ("the scorer reads 0.7803 on that cell",
                    "0.5666 [0.5044, 0.6289]",
                    "the release holds 78,649 nodes"):
        mutated = head + "\n" + planted
        hits = (re.findall(r"\b\d\.\d{4}\b", mutated)
                + re.findall(r"\[\s*\d\.\d+\s*,\s*\d\.\d+\s*\]",
                             mutated)
                + re.findall(r"\b\d{1,3},\d{3}\b", mutated))
        assert hits, (
            f"the no-counts rule does not reject {planted!r}, so it would "
            f"not have caught the figures it was written to remove")


def test_the_no_counts_rule_is_stated_where_a_reader_will_find_it(flat):
    """A RULE ENFORCED BY A TEST AND STATED NOWHERE is a rule the next
    person to edit the file breaks without being told. The policy and
    the check have to travel together, so the test asserts the prose."""
    assert "A count in a README is pinned to one release and goes stale" \
        in flat, "the README no longer states the rule this file enforces"
    assert "tests/test_readme.py" in flat, (
        "the README does not name the test that enforces its no-counts "
        "rule, so a reader cannot find out why a figure was removed")


def test_the_dated_record_declares_that_its_figures_are_frozen(text):
    """The no-counts rule stops at one heading, and the exemption has to
    say why in the document rather than only in this file."""
    head, _, rest = text.partition("## What the build actually taught us")
    assert rest, "the dated-record section is gone; the rule above has no end"
    opening = " ".join(rest.split()[:80])
    assert "DATED RECORD" in opening.upper(), (
        "the section exempted from the no-counts rule does not tell its "
        "reader that its figures are as-measured-then")
    assert "doc2/records/" in rest, (
        "the exemption claims to follow doc2/records/'s convention and "
        "does not name it")


def test_it_points_at_the_live_figures(flat):
    for cmd in ("make audit", "make kg-report", "make check-results"):
        assert cmd in flat, f"the README does not mention {cmd}"
    assert "doc2/kg_audit_report.txt" in flat
    assert "doc2/08-results.md" in flat


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
    # THE CENSUS NOW PRINTS BOTH CASES, and this test is why it has to.
    # It read only "N of M evaluation label predicate(s) carry no edges",
    # which the writer emitted only while one WAS empty. The indications
    # layer filled the last two on 2026-11-17, the line vanished, and
    # this assertion fired on a missing match -- correct, and uninformative
    # about which way the world had moved. write_kg_report.py prints the
    # positive form too now, so both branches below read a statement
    # rather than an absence.
    absent = re.search(r"(\d+) of (\d+) evaluation label", block)
    filled = re.search(r"all (\d+) evaluation label predicate\(s\) carry "
                       r"edges", block)
    assert absent or filled, (
        "the census reports neither an absent nor a filled evaluation "
        "label; the one fact that changes how every result reads is now "
        "unstated")
    if absent:
        n_absent = int(absent.group(1))
        assert f"{'Two' if n_absent == 2 else n_absent} of three evaluation " \
               f"label predicates carry no edges" in flat_readme or \
               f"{n_absent} of three" in flat_readme, (
            f"{n_absent} label predicates are empty and the README does "
            f"not say so")
        return
    # Filled. The README must say so -- and the superseded sentence is
    # kept in quotation marks as the record, the same convention
    # doc2/09-limits.md uses, so the scan that forbids the old
    # claim runs over the prose with quoted material removed.
    assert "All three evaluation label predicates carry edges" in flat_readme, (
        "every label predicate carries edges and the README still reads as "
        "though two were empty")
    unquoted = re.sub(r'"[^"]*"', "", flat_readme)
    assert "carry no edges" not in unquoted, (
        "the superseded claim is asserted outside quotation marks")
