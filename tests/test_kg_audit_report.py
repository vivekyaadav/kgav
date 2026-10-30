"""docs/kg_audit_report.txt, checked against the code it describes.

A narrative report is the worst kind of figure holder: it reads as
authoritative, it is the file someone outside the project is handed
first, and nothing in a .txt fails when the thing it describes changes.
This project's own count is six degree lifts that drifted across
releases while the table sat unchanged, a "best channel" column naming
four rows the code refuses to read, and two reports quoting a figure
from a previous graph.

So every STRUCTURAL claim in that report is re-derived here: the node
classes, the predicates and their held-out flags, the metapath set and
each path's shape, the pipeline layers, and the provenance identifiers.
A predicate renamed, a metapath added, a layer removed or a source
dropped fails here until the report is corrected.

WHAT IS DELIBERATELY NOT CHECKED: release-dependent counts. The report
does not restate them -- it points at `make audit` -- and a test
asserting a node count would pin this file to one release, which is the
defect it exists to prevent. The few run-stamped figures it does quote
are marked with their run in the text.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

REPORT = ROOT / "docs" / "kg_audit_report.txt"
REGISTER = ROOT / "config" / "viruses.yaml"


@pytest.fixture(scope="module")
def text() -> str:
    assert REPORT.exists(), f"{REPORT} is missing"
    return REPORT.read_text()


@pytest.fixture(scope="module")
def flat(text) -> str:
    """Whitespace-normalised. The report is hard-wrapped at 78 columns,
    so every multi-word phrase in it spans a line break somewhere and a
    raw substring search reports present text as absent."""
    return " ".join(text.split())



def section(blob: str, title: str, until: str | None = None) -> str:
    """The text of one numbered section, found by TITLE not by number.

    Tests used to split on "5. METAPATHS". Inserting two sections
    renumbered everything and five tests started raising IndexError on
    a document that was perfectly correct. The number is the one part
    of a heading guaranteed to move.
    """
    m = re.search(rf"\d+\.\s+{re.escape(title)}", blob)
    assert m, f"no section titled {title!r}"
    rest = blob[m.end():]
    if until:
        m2 = re.search(rf"\d+\.\s+{re.escape(until)}", rest)
        assert m2, f"no section titled {until!r} after {title!r}"
        return rest[:m2.start()]
    return rest


@pytest.fixture(scope="module")
def schema():
    from kgav.schema import load_schema
    return load_schema()


# ------------------------------------------------------------------ schema
def test_the_schema_version_is_the_one_described(text, schema):
    assert f"Schema version  {schema.version}" in text


def test_every_node_class_is_listed(flat, schema):
    for name in schema.node_classes:
        assert name in flat, f"node class {name} is not in the report"


def test_the_node_class_count_is_right(flat, schema):
    assert f"{len(schema.node_classes)} node classes" in flat


def test_every_predicate_is_listed(flat, schema):
    for e in schema.edge_classes:
        assert e.predicate in flat, f"{e.predicate} is not in the report"


def test_the_predicate_count_is_right(flat, schema):
    assert f"{len(schema.edge_classes)} predicates" in flat


def test_the_held_out_set_is_exactly_the_schema_s(flat, schema):
    """The four held-out predicates are the whole basis for calling the
    pathway work label-free. If a fifth is added, or one stops being
    held out, that paragraph is wrong.

    SET EQUALITY ON THE ENUMERATION, not membership in the document.
    Every held-out predicate is also named elsewhere -- in the chemistry
    list, in the ChEMBL section -- so a document-wide `pred in flat`
    passed with MEASURED_INACTIVE_AGAINST deleted from the held-out
    list, which is exactly the claim the test existed to protect.
    """
    held = set(schema.held_out_predicates())
    marker = "are HELD OUT of every metapath and are never evidence:"
    assert marker in flat, "the held-out enumeration lost its lead-in"
    block = flat.split(marker, 1)[1].split("Three of those", 1)[0]
    listed = {e.predicate for e in schema.edge_classes
              if e.predicate in block}
    assert listed == held, (
        f"the enumeration lists {sorted(listed)} and the schema holds "
        f"{sorted(held)}")
    assert f"{len(held)} predicates are HELD OUT" in flat or \
        "Four predicates are HELD OUT" in flat, \
        "the report does not state how many predicates are held out"


def test_each_predicate_s_endpoints_match_the_schema(flat, schema):
    """The report draws every predicate as subject -> object. A reversed
    arrow is the kind of error that reads fine and misleads completely:
    HOST_FACTOR_FOR runs Protein -> OrganismTaxon, and the opposite
    would describe a different graph."""
    by_pred = {e.predicate: e for e in schema.edge_classes}
    for pred in ("INHIBITS", "TARGETS", "HOST_FACTOR_FOR",
                 "PARTICIPATES_IN", "PHYSICALLY_INTERACTS_WITH",
                 "ENCODED_BY", "BELONGS_TO", "CAUSES",
                 "HAS_ANTIVIRAL_ACTIVITY_AGAINST", "FOLD_SIMILAR_TO"):
        e = by_pred[pred]
        assert re.search(rf"{pred}\s+{e.subject}\s+->\s+{e.object}", flat), (
            f"the report does not draw {pred} as "
            f"{e.subject} -> {e.object}")


def test_every_symmetric_predicate_is_marked(flat, schema):
    for e in schema.edge_classes:
        if not e.symmetric or e.predicate == "SAME_AS":
            continue
        # ANY mention, not the FIRST. This read
        # flat.split(e.predicate, 1)[1][:120] -- the text after the
        # first occurrence -- so adding a paragraph that merely LISTS
        # the predicate earlier in the document broke it, which is what
        # happened on 2026-10-30 when the Biolink mapping note landed
        # above the predicate table. The claim is that the report marks
        # the predicate symmetric somewhere, not that the first time it
        # appears is that place.
        mentions = flat.split(e.predicate)[1:]
        assert any("symmetric" in m[:120] for m in mentions), (
            f"{e.predicate} is symmetric and the report does not say so "
            f"within 120 characters of any of its {len(mentions)} "
            f"mention(s)")


# --------------------------------------------------------------- metapaths
def test_every_metapath_is_described(flat, schema):
    for name in schema.metapaths:
        assert re.search(rf"\b{name}\b", flat), f"{name} is not described"


def test_the_metapath_count_is_right(flat, schema):
    assert f"{len(schema.metapaths)} reasoning routes" in flat


def test_each_metapath_names_its_own_predicates(flat, schema):
    """Not the shape drawn in ASCII -- the PREDICATES it walks. A
    metapath described with a predicate it does not traverse is a
    reader being told the wrong mechanism."""
    body = section(flat, "METAPATHS").split("HUB EXCLUSION", 1)[0]
    for name, spec in schema.metapaths.items():
        seg = body.split(name, 1)
        assert len(seg) == 2, f"{name} missing from the metapath section"
        # Up to the next metapath heading, or the end.
        rest = seg[1]
        nxt = [rest.index(o) for o in schema.metapaths
               if o != name and o in rest]
        block = rest[:min(nxt)] if nxt else rest
        for step in spec["path"]:
            if step.isupper() or step.startswith("<"):
                pred = step.lstrip("<")
                assert pred in block, (
                    f"{name} walks {pred} and its description does not "
                    f"name it")


def test_the_dependency_constraint_is_stated_for_every_path_that_has_it(
        flat, schema):
    """direction=dependency is the difference between a therapeutic
    target and one it would be harmful to inhibit. Every metapath
    constrained on it must say so."""
    body = section(flat, "METAPATHS").split("HUB EXCLUSION", 1)[0]
    for name, spec in schema.metapaths.items():
        want = (spec.get("constraints") or {}).get(
            "HOST_FACTOR_FOR.direction")
        if want != "dependency":
            continue
        rest = body.split(name, 1)[1]
        nxt = [rest.index(o) for o in schema.metapaths
               if o != name and o in rest]
        block = rest[:min(nxt)] if nxt else rest
        assert "direction=dependency" in block, (
            f"{name} constrains HOST_FACTOR_FOR.direction=dependency and "
            f"its description does not state it")


def test_hub_exclusion_is_stated_for_the_paths_that_use_it(flat, schema):
    excluding = sorted(n for n, s in schema.metapaths.items()
                       if s.get("exclude_hubs"))
    assert excluding, "no metapath excludes hubs; the report says some do"
    for name in excluding:
        if name == "M8":          # described under its own gate instead
            continue
        assert name in flat.split("HUB EXCLUSION", 1)[1][:200] or \
            "hubs excluded" in flat.split(name, 1)[1][:400], (
            f"{name} excludes hubs and the report does not say so")


# ----------------------------------------------------------------- layers
def test_every_pipeline_layer_is_listed(flat):
    from kgav.pipeline import LAYERS
    for layer in LAYERS:
        assert layer.name in flat, f"layer {layer.name} is not in the report"


def test_the_layer_count_is_right(flat):
    from kgav.pipeline import LAYERS
    assert f"{len(LAYERS)} LAYERS" in flat


def test_each_layer_is_marked_required_or_optional_as_declared(text):
    """Two layers are required and the rest are not. A report calling an
    optional layer required implies a build fails without it.

    READS THE RAW TABLE LINES, not the flattened text. Flattened, the
    first occurrence of "similarity" is inside "fold similarity
    [InterPro]" several lines earlier, so a substring search landed on
    prose and reported the layer table as wrong.
    """
    from kgav.pipeline import LAYERS
    body = section(text, "BUILD PIPELINE").split("WHY hosttargets", 1)[0]
    rows = {}
    for line in body.splitlines():
        m = re.match(r"\s{4}(\w+)\s+(required|optional)\b", line)
        if m:
            rows[m.group(1)] = m.group(2)
    assert rows, "no layer rows parsed out of the pipeline table"
    for layer in LAYERS:
        want = "required" if layer.required else "optional"
        assert rows.get(layer.name) == want, (
            f"{layer.name} is {want} in the pipeline and the table says "
            f"{rows.get(layer.name)!r}")


# ----------------------------------------------------------------- sources
def test_every_provenance_identifier_used_in_the_code_is_described(flat):
    """An infores id that edges carry and the report never mentions is a
    source nobody outside the project knows was used."""
    srcs = set()
    for path in list((ROOT / "src" / "kgav").rglob("*.py")) + \
            list((ROOT / "scripts").glob("*.py")):
        srcs |= set(re.findall(r"infores:[a-z0-9][a-z0-9._-]*",
                               path.read_text()))
    srcs = {s.rstrip(".") for s in srcs if s != "infores:"}
    missing = sorted(s for s in srcs if s not in flat)
    assert not missing, (
        f"these sources appear on edges and not in the report: {missing}")


def test_the_register_sources_are_described(flat):
    doc = yaml.safe_load(REGISTER.read_text())
    for key in ("proteome_source", "taxonomy_source", "disease_source"):
        src = doc["provenance"][key]
        assert src in flat, f"register {key}={src} is not in the report"


def test_a_source_surveyed_and_rejected_says_why(flat):
    """"No data" and "cannot be read" are different claims, and this
    project got that exact pair wrong about BioGRID before measuring
    it."""
    block = flat.split("SURVEYED AND DELIBERATELY NOT USED", 1)[1]
    assert "Entrez Gene" in block, (
        "BioGRID is excluded because it keys on Entrez Gene ids, and the "
        "report must say so rather than implying it holds no data")
    assert "redundant with IntAct" in block


# ---------------------------------------------------------------- register
def test_the_virus_counts_match_the_register(flat):
    from kgav.virus_register import family_of
    doc = yaml.safe_load(REGISTER.read_text())
    n = len(doc["viruses"])
    fams: dict[str, int] = {}
    for v in doc["viruses"]:
        fams[family_of(v, doc) or "?"] = fams.get(family_of(v, doc) or "?",
                                                  0) + 1
    assert f"{n} viruses" in flat
    for fam, count in fams.items():
        assert f"{count} {fam}" in flat, (
            f"the register holds {count} {fam} and the report says "
            f"otherwise")


def test_every_virus_in_the_register_is_placed_in_a_tier(flat):
    """Section 7 sorts all sixteen into four tiers. A virus added to the
    register and not to that list is a virus the report silently omits,
    which is how Zika went missing from the metapath table."""
    doc = yaml.safe_load(REGISTER.read_text())
    block = section(flat, "PER-VIRUS STATE", until="KNOWN GAPS")
    for v in doc["viruses"]:
        assert v["name"] in block, (
            f"{v['name']} is in the register and is in no tier in "
            f"section 7")


def _virus_pathways_body(tools_src: str) -> str:
    """The method's own source, sliced out so a claim about it is not
    satisfied by an identifier appearing somewhere else in the file."""
    m = re.search(r"\n    def virus_pathways\(.*?(?=\n    def |\nclass |\Z)",
                  tools_src, re.S)
    return m.group(0) if m else ""


def test_the_virus_only_pathway_claim_matches_the_orchestrator(flat):
    """The report said "for any virus, the graph can list the Reactome
    pathways its dependency factors sit in ... the agent exposes it as
    a PATHWAY query". True of the graph, false of the agent at the
    time, so the claim was rewritten as a gap and this test pinned it
    to the branch.

    IT PINNED THE WRONG BRANCH. The predicate was

        re.search(r"if intent is guards.QueryKind.PATHWAY and compounds:")

    -- the COMPOUND branch, which adding a virus-only form does not
    touch. virus_pathways shipped, the elif beside it appeared, that
    line stayed exactly as it was, and the test went on demanding that
    the report call the capability a gap. Report and guard agreed with
    each other and both described a version of the agent that no
    longer existed; the suite passed green through all of it.

    The lesson is not about this branch. A guard keyed on a line the
    change does not have to touch is a guard that cannot see the
    change. So the predicate now keys on the three things that must
    ALL exist for the capability to be real -- the elif, the call, the
    method -- and deleting any one of them flips this test back to
    demanding the gap text.
    """
    src = (ROOT / "src" / "kgav" / "agent" / "orchestrator.py").read_text()
    tools_src = (ROOT / "src" / "kgav" / "agent" / "tools.py").read_text()
    # CONJUNCTION, deliberately. Three arms that must all hold, not a
    # disjunction where one covers for another.
    has_virus_only = bool(
        re.search(r"elif intent is guards\.QueryKind\.PATHWAY and virus:",
                  src)
        and re.search(r"self\.tools\.virus_pathways\(", src)
        and re.search(r"\n    def virus_pathways\(", tools_src))
    # TWO SCOPED ASSERTIONS, NOT ONE DISJUNCTION. Written as
    # `A in flat or B in flat`, deleting the statement from section 8
    # passed because the KNOWN GAPS entry satisfied the other arm --
    # each arm covering for the other is the same shape as a membership
    # check satisfied by an unrelated mention.
    body = section(flat, "HOW THE GRAPH SUPPORTS PATHWAY ANALYSIS",
                   until="PER-VIRUS STATE")
    gaps = section(flat, "KNOWN GAPS")
    if not has_virus_only:
        assert "no virus-only form of this query" in body.lower(), (
            "the PATHWAY branch requires a compound and section 8 does "
            "not say so where it describes the capability")
        assert "falls through to a node profile" in body
        assert "has no form of" in gaps, (
            "the gap is described in section 8 and missing from KNOWN GAPS")
        return

    assert "no virus-only form of this query" not in body.lower(), (
        "a virus-only PATHWAY form now exists and section 8 still "
        "calls it a gap")
    assert "has no form of" not in gaps, (
        "virus_pathways exists and KNOWN GAPS still lists its absence")
    assert "virus_pathways" in body, (
        "section 8 no longer calls it a gap and does not name the tool "
        "that closed it either")


def test_section_8_does_not_credit_virus_pathways_with_a_p_value(flat):
    """The method reports a ratio against a size-matched expectation
    and NO tail probability, because the pathways it ranks are chosen
    by the same dependency set it then scores them on -- selection and
    outcome are one event. That omission is the most load-bearing
    thing about it and the easiest for a document to quietly upgrade
    into "significant", so it is checked against the method's own
    source rather than against the surrounding prose."""
    tools_src = (ROOT / "src" / "kgav" / "agent" / "tools.py").read_text()
    body_src = _virus_pathways_body(tools_src)
    assert body_src, "virus_pathways not found; this test is vacuous"
    emits_p = bool(re.search(r"p_value|pvalue|\bp_adj\b", body_src))
    body = section(flat, "HOW THE GRAPH SUPPORTS PATHWAY ANALYSIS",
                   until="PER-VIRUS STATE")
    span = body[body.find("virus_pathways"):] if "virus_pathways" in body \
        else ""
    if not emits_p:
        assert span, "section 8 does not describe virus_pathways"
        assert "reports no p-value" in span.lower(), (
            "virus_pathways emits no p-value and section 8 does not "
            "say so where it describes it")
        assert "selection and outcome are one event" in span.lower()
    else:
        assert "reports no p-value" not in span.lower(), (
            "virus_pathways now emits a p-value and section 8 still "
            "says it does not")


def test_section_8_carries_the_names_the_virus_warning(flat):
    """Reactome curates virus-specific pathways from that virus's own
    literature -- the literature the screens come from -- so such a
    pathway scores high because of where its members were found. The
    code marks them; a report that described the ranking without the
    mark would be describing a cleaner instrument than exists."""
    tools_src = (ROOT / "src" / "kgav" / "agent" / "tools.py").read_text()
    body_src = _virus_pathways_body(tools_src)
    assert body_src, "virus_pathways not found; this test is vacuous"
    if "names_the_virus" not in body_src:
        return  # covered by tests/test_virus_pathways.py
    body = section(flat, "HOW THE GRAPH SUPPORTS PATHWAY ANALYSIS",
                   until="PER-VIRUS STATE")
    assert "NAMES THIS VIRUS" in body, (
        "the tool flags virus-named pathways and section 8 does not "
        "mention the flag")
    assert "marked rather than dropped" in body.lower(), (
        "section 8 mentions the flag without saying the row is kept")


def test_the_report_does_not_claim_the_agent_walks_m4(flat):
    """tools.pathways is a census, not a metapath traversal, and its own
    docstring says _find_paths walks M1 and M2 only. A report implying a
    pathway-mediated ranking would be describing a capability that does
    not exist."""
    src = (ROOT / "src" / "kgav" / "agent" / "tools.py").read_text()
    assert "_find_paths walks M1 and M2 only" in src, (
        "the agent's traversal set changed; re-check what the report "
        "claims about M4")
    assert "The agent never walks M4" in flat


# ------------------------------- figures, now that the census is elsewhere
def test_the_report_restates_no_release_total(text):
    """With the census in its own file, the report carries no release
    count at all -- so the rule is simple again: any total here is
    typed, and a typed total is stale the moment the graph is rebuilt.
    That is the defect that produced two reports describing a release
    that no longer existed.

    The historical lessons are exempt by scope, not by wording: section
    5 and after may quote the figures a past defect produced, because
    those are a record of what happened rather than a claim about now.
    """
    head = text.split("5. DATA SOURCES", 1)[0]
    bad = re.findall(r"\b\d{1,3},\d{3}\s+(?:nodes|edges)\b", head)
    assert not bad, f"the report restates a release total: {bad}"


def undescribed_curators() -> list[str]:
    """infores ids the census shows carrying edges that the report
    never mentions. Factored out so it can be tested on a synthetic
    census as well as on the committed one."""
    import scripts.write_kg_report as W
    census = W.CENSUS.read_text()
    return _undescribed(census, REPORT.read_text())


def _undescribed(census: str, prose: str) -> list[str]:
    ids = {i.rstrip(".") for i in
           re.findall(r"infores:[a-z0-9][a-z0-9._-]*", census)}
    return sorted(i for i in ids if i not in prose)


def test_the_curator_check_catches_one_the_report_never_names():
    """Exercised on a synthetic census, so it runs whether or not a
    release has been assembled here. A test that returns early when its
    input is missing is a test that reports green having checked
    nothing."""
    prose = "infores:uniprot is described in the prose."
    census = "  infores:uniprot 10\n  infores:nosuchdb 1\n"
    assert _undescribed(census, prose) == ["infores:nosuchdb"]
    assert _undescribed(census, prose + " infores:nosuchdb too.") == []


def test_every_curator_in_the_census_is_described_in_the_report():
    """THE GAP THIS CLOSES. The sibling test reads infores ids out of
    the SOURCE CODE, and several -- infores:dip, infores:mint,
    infores:hpidb -- are slugged at runtime from a curator's own name in
    MITAB column 12, so they are never literals anywhere in the tree.
    v0.1 carries edges from all three while the report said those
    endpoints returned nothing usable. The census is the only place
    their names appear, so it is what the report must agree with."""
    missing = undescribed_curators()
    assert not missing, (
        f"these curators carry edges in the release and are not "
        f"described anywhere in the report: {missing}")


def test_section_9_describes_the_backstop_that_exists(flat):
    """Section 9 claims the agent refuses candidate-seeking questions
    that are on no keyword list, and that the prospective refusal's
    closing offer narrows for a routeless virus. Both are claims about
    code, so both are checked against it.

    KEYED ON WHAT THE CHANGE WOULD HAVE TO TOUCH. The lesson from the
    virus-only pathway guard, one patch earlier, is that a predicate
    reading a line the change does not have to modify cannot see the
    change. Removing the backstop means removing _seeks_candidates or
    its call; narrowing the offer means the routes branch inside the
    PROSPECT arm. Those are the lines pinned here.
    """
    src = (ROOT / "src" / "kgav" / "agent" / "guards.py").read_text()
    has_backstop = bool(
        re.search(r"\ndef _seeks_candidates\(", src)
        and re.search(r"if _seeks_candidates\(t\):", src))
    has_route_aware_offer = bool(
        re.search(r"if routes is not None and not routes\.routes:\s*\n\s+"
                  r"offer = \(", src))
    body = section(flat, "PER-VIRUS STATE", until="KNOWN GAPS")

    if has_backstop:
        assert "classify_intent" in body, (
            "a structural backstop exists and section 9 does not say how "
            "a question reaches the refusal at all")
        assert "failing closed" in body.lower()
        # THE CONTRAST, NOT THE PHRASE. `"which compounds should" in
        # body` passed with the claim deleted, because the report
        # quotes that query verbatim four lines above -- a membership
        # check satisfied by an unrelated occurrence of the same
        # string, which is the defect this suite has now found ten
        # times and which I reintroduced while writing a comment about
        # it. Caught by the mutant, not by reading.
        assert ('"which drug should" was on it and "which compounds '
                'should" was not') in body, (
            "the pair of phrasings IS the evidence that a keyword list "
            "leaks; without the contrast the claim is an assertion")
    else:
        assert "failing closed" not in body.lower(), (
            "the backstop is gone and section 9 still describes it")

    if has_route_aware_offer:
        assert "route-aware" in body.lower(), (
            "the prospective refusal narrows its offer and section 9 "
            "does not say so")
    else:
        assert "route-aware" not in body.lower()


def test_section_9_does_not_claim_a_model_classifies_intent(flat):
    """classify_intent is keyword-plus-structure and deliberately not a
    model. A report implying otherwise would describe a refusal that
    can vary between runs, which is the property this layer exists to
    not have."""
    src = (ROOT / "src" / "kgav" / "agent" / "guards.py").read_text()
    assert "Deliberately keyword-based rather than model-based" in src, (
        "classify_intent's basis changed; re-check what section 9 claims")
    body = section(flat, "PER-VIRUS STATE", until="KNOWN GAPS")
    assert "keyword-based on purpose" in body


def test_section_6_reproducibility_claims_match_the_code(flat):
    """Section 6 states the release is reproducible. That rests on
    three mechanisms, each in a different file, and the claim must not
    outlive any of them.

    KEYED ON WHAT REMOVING THE BEHAVIOUR WOULD HAVE TO TOUCH -- the
    sort in Assembly.write, the tie-break in similarity_edges, and the
    manifest's code stamp. Not on a line a future change can leave
    alone, which is how the virus-only pathway guard went blind.
    """
    asm = (ROOT / "src" / "kgav" / "assemble.py").read_text()
    sim = (ROOT / "src" / "kgav" / "similarity.py").read_text()
    rel = (ROOT / "scripts" / "assemble_release.py").read_text()

    sorted_write = bool(
        re.search(r"nodes = sorted\(self\.nodes\.values\(\)", asm)
        and re.search(r"edges = sorted\(self\.edges\.values\(\)", asm)
        and asm.count("json.dumps(n, sort_keys=True)")
        and asm.count("json.dumps(e, sort_keys=True)"))
    tiebreak = bool(
        re.search(r"hits\.sort\(key=lambda kv: \(-kv\[1\], kv\[0\]\)\)", sim))
    stamped = 'manifest["code"] = code_commit(' in rel

    body = section(flat, "BUILD PIPELINE", until="METAPATHS")
    claims = "the graph itself is reproducible" in body.lower()

    if sorted_write and tiebreak and stamped:
        assert claims, (
            "all three mechanisms are in place and section 6 does not "
            "say the release is reproducible")
        # The mechanisms, named, so the claim can be checked by a reader
        assert "no canonical form" in body
        assert "tie-break" in body
        assert "MANIFEST.json carries the commit" in body
    else:
        assert not claims, (
            "section 6 claims reproducibility and at least one of the "
            "three mechanisms behind it is gone: "
            f"sorted_write={sorted_write} tiebreak={tiebreak} "
            f"stamped={stamped}")


def test_section_6_does_not_resolve_what_it_could_not_resolve(flat):
    """The 542-vs-568 discrepancy concerns a release that no longer
    exists. A later edit that quietly turns it into an explanation
    would be claiming evidence nobody has."""
    body = section(flat, "BUILD PIPELINE", until="METAPATHS")
    assert "542" in body and "568" in body
    assert "unexplained historical discrepancy" in body
    assert "overwritten before anyone thought to keep" in body


def test_section_1_counts_the_channels_against_the_results(flat):
    """Section 1 states two numbers that are derivable, and nothing
    derived them.

    It said "exactly one readable cell in the whole table is
    significant". On 2026-10-30 raising the path cap completed M4's
    traversal and DENV-2's M4 arrived at p=0.019, making two -- and no
    test anywhere noticed, because the figure was typed. The same
    claim in docs/evaluation.md HAS been guarded since it was written;
    this copy of it was not.

    Both numbers are read off metapath_comparison.json, so the next
    move is caught wherever it happens rather than wherever someone
    remembered to look.
    """
    import json

    res = ROOT / "data" / "results" / "metapath_comparison.json"
    assert res.exists(), (
        "metapath_comparison.json is missing, so section 1's counts "
        "cannot be checked -- and an unchecked count is how this one "
        "went stale")
    per = json.loads(res.read_text())["per_virus"]
    readable = [(v, n, m) for v, ch in per.items() for n, m in ch.items()
                if n not in ("degree", "COMBINED")
                and m.get("hits_determined") and m.get("trustworthy")]
    sig = [m for _, _, m in readable
           if m.get("p_value") is not None and m["p_value"] < 0.05]
    with_channel = {v for v, _, _ in readable}
    without = len(per) - len(with_channel)

    assert f"{without} of {len(per)} scored viruses have no" in flat, (
        f"{without} of {len(per)} viruses have no readable channel; "
        f"section 1 states a different number")
    assert f"{len(sig)} readable cells in the whole" in flat, (
        f"{len(sig)} readable cells are significant uncorrected; "
        f"section 1 states a different number")


def test_biolink_conformance_is_claimed_only_if_something_checks_it(flat):
    """"Biolink-conformant" sat in line 3 from the first version of
    this report until 2026-10-30. The word `biolink` appeared NOWHERE
    else in the repository -- no linkml, no biolink-model dependency,
    no validation -- and measured against the upstream model it was
    8 of 11 node classes and 8 of 24 predicates by name.

    So the claim may be made only when a mapping exists to check it
    against. Until then the report says ALIGNED and shows its working,
    which is what every other claim in it does.
    """
    mapping_exists = bool(
        list(ROOT.glob("schema/*biolink*"))
        or "biolink_mapping" in (ROOT / "schema" / "kg_schema.yaml").read_text())
    # THE TITLE, not the document. The section-3 note says "BIOLINK-
    # ALIGNED, NOT BIOLINK-CONFORMANT" and explains why -- a check
    # over the whole text fires on that negation, which is a
    # membership check satisfied by an unrelated occurrence, written
    # here by the same hand that keeps finding it elsewhere.
    title = flat.split("Repository", 1)[0]
    conformant = "biolink-conformant" in title.lower()
    if not mapping_exists:
        assert not conformant, (
            "the report claims Biolink CONFORMANCE and nothing in the "
            "repository maps a single class or predicate to the model; "
            "say ALIGNED, or add schema mappings and a check")
        assert "biolink-aligned" in title.lower(), (
            "the report should say what its relationship to Biolink is, "
            "not go silent on it")


def test_the_biolink_note_keeps_the_two_divergences_that_are_deliberate(flat):
    """MEASURED_INACTIVE_AGAINST and HOST_FACTOR_FOR.direction are not
    gaps to close. Biolink cannot express a measurement that found no
    activity, and `contributes to` drops the dependency/restriction
    sign -- the one that separates a therapy from a harm. A later
    mapping effort that quietly folded them in would lose the
    hard-negative protocol and the direction, so the reasons are
    recorded where someone doing that work would read them."""
    block = section(flat, "SCHEMA SPECIFICATION", until="RELEASE CENSUS")
    assert "MEASURED_INACTIVE_AGAINST" in block
    assert "no way to say that" in block
    assert "HOST_FACTOR_FOR" in block
    assert "dependency versus restriction" in block
