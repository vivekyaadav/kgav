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


def test_every_make_target_the_report_cites_exists(text):
    """The report attributes figures to commands. On 2026-11-02 section 9
    opened "From `make audit`" -- and there is no `audit` target in the
    Makefile, and never has been. The figures under it were correct; the
    provenance line was a command nobody could run, which is the one
    thing a reader would use to check them.

    This is the cheapest guard in the file and it covers every section:
    a target renamed in the Makefile breaks every document that cites it
    by the old name, silently.
    """
    # RETRACTIONS ARE NOT CITATIONS. This report records what it used to
    # say -- three paragraphs open "An earlier version" -- and one of
    # them quotes the `make audit` that this test was written to catch.
    # Forbidding the words outright would forbid recording the mistake,
    # which is the same scoping every other guard in this suite applies
    # to its own claim block.
    live = "\n\n".join(p for p in text.split("\n\n")
                       if "An earlier version" not in p)
    assert "An earlier version" in text, (
        "no retraction paragraph in the report, so the exclusion above "
        "is not being exercised and this test is checking the whole text "
        "under a rule that looks narrower than it is")
    cited = set(re.findall(r"make ([a-z][a-z0-9-]*)", live))
    targets = set(re.findall(r"^([A-Za-z][A-Za-z0-9_-]*):",
                             (ROOT / "Makefile").read_text(), re.M))
    # "make the" and the like: prose, not a target. Only a citation that
    # LOOKS like a target name is checked, which is why this reads the
    # Makefile rather than guessing.
    prose = {"the", "a", "an", "it", "sure", "this", "that"}
    missing = sorted(c for c in cited - targets - prose)
    assert not missing, (
        f"the report cites make target(s) that the Makefile does not "
        f"define: {missing}")


def tier_membership() -> dict[str, set[str]]:
    """The four tiers of section 9, derived.

    Section 9 named them as output of `make audit`, which does not
    exist. They are in fact a function of two tracked files, so they are
    derived from those: COMBINED's pool and pos in
    metapath_comparison.json, plus the register for the virus that is
    scored at all.
    """
    import json

    per = json.loads((ROOT / "data" / "results"
                      / "metapath_comparison.json").read_text())["per_virus"]
    register = {v["name"] for v in yaml.safe_load(REGISTER.read_text())
                ["viruses"]}
    combined = {v: ch.get("COMBINED", {}) for v, ch in per.items()}
    return {
        # A route that reaches at least one held-out positive.
        "EVALUABLE": {v for v, c in combined.items() if c.get("pos", 0) > 0},
        # A route that reaches compounds, none of them labelled.
        "ROUTES REACH NO LABEL": {v for v, c in combined.items()
                                  if c.get("pos") == 0 and c.get("pool", 0) > 0},
        # Labels exist and no route terminates at the virus.
        "LABELS, NO FEATURE ROUTE": {v for v, c in combined.items()
                                     if c.get("pool") == 0},
        # Scored nowhere, because there is nothing to score against.
        "FEATURES, NOTHING TO SCORE": register - set(per),
    }


@pytest.mark.parametrize("tier", ["EVALUABLE", "ROUTES REACH NO LABEL",
                                  "LABELS, NO FEATURE ROUTE",
                                  "FEATURES, NOTHING TO SCORE"])
def test_each_tier_holds_the_viruses_the_results_file_puts_in_it(text, tier):
    """SET EQUALITY ON THE TIER'S OWN LIST, not membership in the
    section. The existing check that every register virus appears
    somewhere in section 9 is satisfied by any of the four tiers, so a
    virus in the wrong tier passed it -- and the tiers are not degrees
    of weak evidence, they are different kinds of absence. Putting a
    virus one tier up turns "nothing can mark this ranking right or
    wrong" into "this was evaluated".
    """
    want = tier_membership()[tier]
    body = section(" ".join(text.split()), "PER-VIRUS STATE",
                   until="KNOWN GAPS")
    m = re.search(rf"{re.escape(tier)} \((\d+)\)", body)
    assert m, f"section 9 has no {tier!r} tier with a count"
    assert int(m[1]) == len(want), (
        f"the {tier} tier is headed ({m[1]}) and the results file puts "
        f"{len(want)} viruses in it: {sorted(want)}")
    # The names run from the tier heading to the NEXT TIER HEADING, found
    # by the four known tier names rather than by a shape. A generic
    # "capitals then (N)" boundary matched inside the tier's own virus
    # list -- "JEV, TBEV LABELS, NO FEATURE ROUTE (2)" is all capitals
    # and commas -- so the span stopped two viruses early and the test
    # failed on a correct document. The boundary is the thing being
    # parsed, so it is named, not guessed.
    rest = body[m.end():]
    bounds = [b.start() for name in tier_membership()
              for b in [re.search(rf"{re.escape(name)} \(\d+\)", rest)] if b]
    bounds += [b.start() for b in [re.search(r"THE AGENT ENFORCES", rest)] if b]
    span = rest[:min(bounds)] if bounds else rest
    listed = {v for v in {x["name"] for x in
                          yaml.safe_load(REGISTER.read_text())["viruses"]}
              if re.search(rf"(?<![\w-]){re.escape(v)}(?![\w-])", span)}
    assert listed == want, (
        f"the {tier} tier lists {sorted(listed)} and the results file "
        f"gives {sorted(want)}")


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


def test_section_8_quotes_the_landscape_run_it_names(flat):
    """SECTION 8 KEEPS A SECOND COPY OF THE PATHWAY RESULT, and on
    2026-11-02 every figure in it was stale: the run date, the
    uncorrected count, and both DENV naive/matched pairs.

    docs/evaluation.md holds the same four figures and has been guarded
    against the results file since the landscape table was written.
    THIS copy was not -- the same split that let section 1's channel
    counts drift while the identical claim in evaluation.md was
    checked. A figure duplicated into a second document needs a second
    guard, or the document that is not checked becomes the one people
    quote.
    """
    import json

    res = ROOT / "data" / "results" / "pathway_landscape.json"
    assert res.exists(), (
        "pathway_landscape.json is missing, so section 8's result "
        "paragraph cannot be checked -- and an unchecked figure here is "
        "how all four of them went stale")
    d = json.loads(res.read_text())
    body = section(flat, "HOW THE GRAPH SUPPORTS PATHWAY ANALYSIS",
                   until="PER-VIRUS STATE")

    assert f"run {d['measured_on']}" in body, (
        f"the run is dated {d['measured_on']} in the results file and "
        f"section 8 names a different date")
    assert f"degree_basis={d['degree_basis']}" in body
    assert f"{d['permutations']:,} permutations" in body
    assert f"seed {d['seed']}" in body
    assert f"{d['n_tests']} tests" in body

    rows = [r for ch in ("taxa", "interactome") for r in d[ch].values()]
    alpha = 0.05 / d["n_tests"]
    clearing = [r for r in rows
                if min(r["p_matched"], r["p_matched_less"]) < alpha]
    assert not clearing, (
        "a row now clears Bonferroni and section 8 still says zero do")
    assert f"{len(rows)} rows tested" in body

    # THE UNCORRECTED COUNT IS THE POINT OF THE PARAGRAPH. It said four
    # and the file gives six.
    WORDS = {0: "zero", 1: "one", 2: "two", 3: "three", 4: "four",
             5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine",
             10: "ten", 11: "eleven", 12: "twelve"}
    uncorrected = [r for r in rows
                   if min(r["p_matched"], r["p_matched_less"]) < 0.05]
    assert (f"would have reported {WORDS[len(uncorrected)]}") in body, (
        f"{len(uncorrected)} rows clear 0.05 uncorrected and section 8 "
        f"states a different number")

    # Both DENV pairs, scoped to this section.
    for curie, name in (("NCBITaxon:11053", "DENV-1"),
                        ("NCBITaxon:11060", "DENV-2")):
        r = d["interactome"][curie]
        assert f"{r['p_naive']:.4f}" in body, (
            f"{name}'s naive p is {r['p_naive']:.4f} and section 8 does "
            f"not carry it")
        assert f"{r['p_matched']:.4f}" in body, (
            f"{name}'s degree-matched p is {r['p_matched']:.4f} and "
            f"section 8 does not carry it")


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
    # The noun agrees with the derived count, for the reason the same rule
    # was added to test_evaluation_table.py on 2026-11-02: at one cell this
    # demanded "1 readable cells in the whole table", and a test that can
    # only be satisfied by bad grammar is a test the document will be
    # edited around. The rule is here; the report quotes what it produces.
    noun = "cell" if len(sig) == 1 else "cells"
    assert f"{len(sig)} readable {noun} in the whole" in flat, (
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


# ------------------------- the prose layer, re-derived from the tree
# THE AUDIT'S CLOSING DIAGNOSIS, verbatim: "the prose layer moves slower
# than the code, and nothing tests the prose that describes METHOD; the
# existing tests pin structure well and pin no method sentence." These
# are the sentences that were wrong on 2026-11-02, each pinned to the
# thing it describes so correcting it once is enough.
SCHEMA_FILE = ROOT / "schema" / "kg_schema.yaml"
MAKEFILE = ROOT / "Makefile"


def test_the_layout_block_names_the_schemas_real_directory(text):
    """It said `config/` held the schema. The schema is in schema/, and
    schema/ was absent from the block entirely."""
    assert SCHEMA_FILE.exists()
    block = text.split("DIRECTORY LAYOUT")[1][:1800]
    assert "schema/" in block, block[:400]
    assert "viruses.yaml (the register), schema," not in block, \
        "config/ is still described as holding the schema"


def test_the_module_count_is_the_number_of_modules(text):
    """Re-derived. It said 28 where there are 27 plus __init__.py, and
    enumerated "one per source layer" -- which is true of scripts/ and
    not of this directory, because there is no spine.py."""
    mods = sorted(p.name for p in (ROOT / "src" / "kgav").glob("*.py")
                  if p.name != "__init__.py")
    # Whitespace-normalised: the block is wrapped prose, so a phrase
    # that spans a line break is invisible to a raw search. The first
    # version of this test looked for "no spine.py" in the raw slice and
    # failed against text that contained it across two lines.
    block = " ".join(text.split("DIRECTORY LAYOUT")[1][:2200].split())
    assert f"{len(mods)} modules" in block, (
        f"{len(mods)} modules on disk; the report says otherwise")
    assert not (ROOT / "src" / "kgav" / "spine.py").exists()
    assert "no spine.py" in block


def test_the_entry_point_count_is_the_number_of_entry_points(text):
    """61 .py plus tier_spotcheck.sh. It said 57."""
    n = len([p for p in (ROOT / "scripts").iterdir()
             if p.suffix in (".py", ".sh")])
    block = text.split("DIRECTORY LAYOUT")[1][:1800]
    assert f"{n} entry points" in block, (
        f"{n} entry points in scripts/; the report says otherwise")


def test_no_predicate_qualifier_is_named_that_the_schema_lacks(text):
    """`pmid` was listed as a qualifier and appears NOWHERE in
    kg_schema.yaml. Publications ride on the edge's `publications` list."""
    schema_text = SCHEMA_FILE.read_text()
    assert "pmid" not in schema_text
    assert "detection_method, pmid" not in text, \
        "pmid is still listed as a qualifier"


def test_hub_exclusion_is_claimed_for_exactly_the_gated_metapaths(text):
    """Re-derived from `exclude_hubs` in the schema. The sentence named
    M4 and M5; M8 carries the flag too and was missing."""
    doc = yaml.safe_load(SCHEMA_FILE.read_text())
    gated = sorted(name for name, mp in (doc.get("metapaths") or {}).items()
                   if isinstance(mp, dict) and mp.get("exclude_hubs"))
    assert gated, "no metapath declares exclude_hubs; this test is inert"
    line = next(ln for ln in text.splitlines() if ln.startswith("HUB EXCLUSION"))
    for name in gated:
        assert name in line, f"{name} gates hubs and is not in: {line}"


def test_the_manifest_claim_matches_what_git_tracks(text):
    """THIS TEST USED TO ASSERT ONE OF THE TWO STATES, AND THE STATE MOVED.

    It was `assert not tracked` plus `"except MANIFEST.json" not in text`,
    written when the report claimed a manifest was tracked and none was --
    `data/releases/*` ignored the directory and a negation cannot reach
    inside an ignored one, so `!data/releases/**/MANIFEST.json` was inert.
    P12 fixed the depth on 2026-11-09, the report's original claim became
    TRUE, and this test failed for the sentence being right.

    A document check belongs between the document and the world, not
    against a state that was true on the day it was written. So it now
    compares the two in whichever direction they point, and fails only
    when they disagree.
    """
    import subprocess
    out = subprocess.run(["git", "ls-files", "data/releases"], cwd=ROOT,
                         capture_output=True, text=True).stdout
    tracked = [ln for ln in out.splitlines() if "MANIFEST.json" in ln]
    claims_tracked = "except the smoke release and each release's" in text
    denies_tracked = "NO MANIFEST.json is" in text
    assert claims_tracked != denies_tracked, (
        "the report says both or neither about MANIFEST.json; this check "
        "cannot compare a claim it cannot read")
    if tracked:
        assert claims_tracked, (
            f"{tracked} is tracked and the report still says no manifest is")
    else:
        assert denies_tracked, (
            "no MANIFEST.json is tracked and the report claims one is -- "
            "which is the state this test was originally written for")


def test_the_provenance_requirement_is_credited_where_it_lives(text):
    """emit.py is 42 lines and holds no validator. Enforcement is
    schema.Schema.validate_batch, and a layer is checked only if its
    ingest calls it -- four of twelve do not, host among them."""
    emit = (ROOT / "src" / "kgav" / "emit.py").read_text()
    assert "validate_batch" not in emit, \
        "emit.py validates now; the report paragraph needs revisiting"
    skipping = sorted(
        p.name for p in (ROOT / "scripts").glob("ingest_*.py")
        if "validate_batch" not in p.read_text())
    assert skipping, "every ingest validates now; correct the report"
    for name in skipping:
        assert name[: -len(".py")] in text, (
            f"{name} skips validate_batch and the report does not say so")
    # AND THE CLAIM ITSELF, not just its surroundings. The first version
    # of this test asserted only that the four skipping ingests were
    # NAMED somewhere, so restoring "Schema validation is enforced in
    # emit.py" survived it -- a test checking something adjacent to the
    # sentence rather than the sentence.
    para = " ".join(text.split("ONE WRITER")[1][:900].split())
    assert "validate_batch" in para, (
        "the paragraph does not say where the rules actually live")
    for claim in ("enforced in emit.py", "enforced there"):
        assert claim not in para, (
            f"the report still credits emit.py with enforcement: {claim!r}")


def test_no_make_target_is_denied_that_the_makefile_defines(text):
    """D3 and the shape of it. "There is no `audit` target in the
    Makefile and there never has been" stood in the present tense after
    Makefile:141 defined one -- inside a RETRACTION paragraph, which is
    the one part of this report nothing checks, because the exclusion
    that lets a retraction quote the old wrong text also exempts its own
    claims. Checked on the WHOLE text for this class only."""
    targets = set(re.findall(r"^([A-Za-z][A-Za-z0-9_-]*):",
                             MAKEFILE.read_text(), re.M))
    assert "audit" in targets, "the audit target went away"
    flat = " ".join(text.split())
    for t in sorted(targets):
        for claim in (f"no `{t}` target", f"no {t} target"):
            assert claim not in flat, (
                f"the report denies the `{t}` target, which the Makefile "
                f"defines")


def test_a_built_layer_is_not_described_as_not_built(text):
    """gtop. The census -- generated from the release and stamped with
    its content hashes -- lists it among the assembled layers with 6,765
    edges, and three documents said it was not built."""
    census = (ROOT / "docs" / "kg_release_census.txt")
    if not census.exists():
        pytest.skip("census absent")
    body = census.read_text()
    assembled = body.split("layer(s) assembled")[1].split("\n\n")[0]
    built = {ln.strip() for ln in assembled.splitlines() if ln.strip()
             and not ln.strip().startswith("hub")}
    assert "gtop" in built, "gtop is no longer assembled; recheck the prose"
    # THE LAYER TABLE, which is where the claim actually sits. Restoring
    # "gtop  optional  second TARGETS source  NOT BUILT" survived the
    # first version of this test, because that version looked only for
    # the phrase "needs the" in prose.
    for line in text.splitlines():
        if "NOT BUILT" not in line:
            continue
        layer = line.split()[0]
        assert layer not in built, (
            f"the layer table marks {layer} NOT BUILT and the census "
            f"lists it among the assembled layers: {line.strip()}")
    for doc in ("docs/kg_audit_report.txt", "README.md",
                "docs/v1_release_plan.md"):
        d = (ROOT / doc).read_text()
        for layer in sorted(built):
            # The claim shape, not the word: a retraction may say a layer
            # "WAS in this sentence and is built".
            assert f"{layer} needs the" not in d, (
                f"{doc} says the built layer {layer} still needs something")
            assert f"`{layer}` needs a" not in d, (
                f"{doc} says the built layer {layer} still needs something")


# ------------------------- the TARGETS direction figure, in SEVEN places
# D2. The audit named one file. The pre-gtop figure "12,155 of 14,381
# TARGETS edges state no direction (84.5%)" was quoted in seven places --
# gtop.py, assemble.py, guards.py, ingest_gtop.py twice,
# v1_release_plan.md and scope_and_limits.md -- every one of them
# labelling it as v0.1's. It is not: it describes the release BEFORE the
# gtop layer landed. Correct numbers attributed to the wrong release,
# which is harder to spot than a wrong number.
#
# On the assembled v0.1: 10,395 of 17,982 unstated (57.8%), directed
# 2,226 -> 7,587. Both halves of that moved, so quoting either one
# without its release is ambiguous.
CENSUS = ROOT / "docs" / "kg_release_census.txt"
RELEASE = ROOT / "data" / "releases" / "v0.1"

#: Files that quote a TARGETS total. Kept as a list so a new quote in a
#: file nobody thought of still has to pick a side.
QUOTING_FILES = (
    "src/kgav/gtop.py",
    "src/kgav/assemble.py",
    "src/kgav/agent/guards.py",
    "scripts/ingest_gtop.py",
    "docs/v1_release_plan.md",
    "docs/scope_and_limits.md",
)


def _census_targets() -> int:
    for line in CENSUS.read_text().splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[0] == "TARGETS":
            return int(parts[1].replace(",", ""))
    raise AssertionError("no TARGETS row in the census")


def test_the_census_still_holds_a_targets_total():
    assert _census_targets() > 0


@pytest.mark.parametrize("path", QUOTING_FILES)
def test_no_file_quotes_the_pre_gtop_total_as_this_release(path):
    """The denominator pinned against the committed census, which is
    stamped against the release and so needs no release on disk.

    A file MAY still quote 14,381 -- the pre-gtop state is the reason the
    layer exists and deleting it would lose the argument -- but not
    beside a claim that it is this release. The check is the pairing, not
    the number."""
    total = _census_targets()
    text = (ROOT / path).read_text()
    flat = " ".join(text.split())
    for stale in ("14,381", "12,155", "2,226"):
        if stale not in flat:
            continue
        # Wherever a pre-gtop number appears, the sentence around it must
        # say so, and the current total must appear somewhere in the file.
        assert f"{total:,}" in flat, (
            f"{path} quotes the pre-gtop figure {stale} and never quotes "
            f"this release's {total:,}")
        for claim in (f"{stale} on v0.1", f"{stale} in v0.1",
                      f"v0.1's {stale}"):
            assert claim not in flat, (
                f"{path} attributes the pre-gtop figure to v0.1: {claim!r}")


def test_the_direction_split_matches_the_release_when_one_is_present():
    """RE-DERIVED where it can be. The quoted numerator cannot come from
    the census -- it holds no qualifier breakdown -- so on a machine with
    the release this recounts it, and elsewhere the test above pins the
    denominator. Adding the split to the census would be better and needs
    a regeneration step, so it is a follow-up rather than this patch."""
    if not (RELEASE / "edges.jsonl").exists():
        pytest.skip("no full release on disk")
    import json
    total = unknown = 0
    with (RELEASE / "edges.jsonl").open() as fh:
        for line in fh:
            e = json.loads(line)
            if e.get("predicate") != "TARGETS":
                continue
            total += 1
            if ((e.get("qualifiers") or {}).get("direction") or "unknown") \
                    == "unknown":
                unknown += 1
    assert total == _census_targets(), (
        f"the release holds {total:,} TARGETS and the census says "
        f"{_census_targets():,}; run make kg-report")
    scope = (ROOT / "docs" / "scope_and_limits.md").read_text()
    assert f"{unknown:,} of {total:,}" in " ".join(scope.split()), (
        f"scope_and_limits.md does not state {unknown:,} of {total:,}")
    pct = f"{unknown / total:.1%}"
    assert pct in scope, f"the quoted percentage is not {pct}"
