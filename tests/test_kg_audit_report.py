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
        after = flat.split(e.predicate, 1)[1][:120]
        assert "symmetric" in after, (
            f"{e.predicate} is symmetric and the report does not say so")


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


def test_the_tier_totals_add_up_to_the_register(flat):
    """(8) + (5) + (2) + (1) must be the register's size, or one virus is
    counted twice or not at all."""
    doc = yaml.safe_load(REGISTER.read_text())
    block = section(flat, "PER-VIRUS STATE", until="KNOWN GAPS")
    totals = [int(n) for n in re.findall(r"\((\d+)\) —", block)]
    assert len(totals) == 4, f"expected four tiers, found {totals}"
    assert sum(totals) == len(doc["viruses"]), (
        f"tiers sum to {sum(totals)} and the register holds "
        f"{len(doc['viruses'])}")


# ------------------------------------------------------------------ claims
def test_a_release_total_appears_only_inside_the_generated_block(text):
    """Sharpened, not dropped.

    This forbade release totals ANYWHERE, which is why the report
    carried none -- and that was the wrong rule. An audit report is a
    stamped snapshot of one release; omitting its size makes it useless
    to a reader and protects nothing. The defect it was aimed at is a
    total that is TYPED: unstamped, unchecked, and stale the moment the
    graph is rebuilt. That is two reports in this project describing a
    graph that no longer existed.

    So a total inside the generated census is the point of the census.
    A total outside it is the defect, and still fails here.
    """
    import scripts.write_kg_report as W
    outside = text.split(W.BEGIN, 1)[0] + text.split(W.END, 1)[1]
    bad = re.findall(r"\b\d{1,3},\d{3}\s+(?:nodes|edges)\b", outside)
    assert not bad, (
        f"a release total is typed outside the generated census: {bad}. "
        f"Totals belong in the census, where they carry the release hash "
        f"they were computed from")


def test_a_generated_census_carries_the_totals(text):
    """The other direction. A census block that somehow holds no totals
    is a section that looks generated and says nothing, which reads
    worse than the omission it replaced."""
    import scripts.write_kg_report as W
    block = text.split(W.BEGIN, 1)[1].split(W.END, 1)[0]
    if "not yet generated" in block:
        return                       # covered by the placeholder test
    assert re.search(r"\b[\d,]+\s+nodes\b", block), block[:300]
    assert re.search(r"\b[\d,]+\s+edges\b", block), block[:300]
    assert "NODES BY CLASS" in block and "EDGES BY PREDICATE" in block


def test_the_report_points_at_the_live_audit(flat):
    assert "make audit" in flat
    assert "make check-results" in flat


def test_the_retracted_finding_is_recorded_as_retracted(flat):
    """The whole-graph degree matching was circular and its uniform
    depletion was an artifact. A report that simply omitted it would
    read as though the error never happened."""
    assert "circular" in flat
    assert "retracted" in flat


# ------------------------------------------------- the generated census
def test_the_census_block_has_its_markers(text):
    """Without them write_kg_report.py cannot splice and the section
    silently stops being regenerable."""
    import scripts.write_kg_report as W
    assert text.count(W.BEGIN) == 1
    assert text.count(W.END) == 1
    assert text.index(W.BEGIN) < text.index(W.END)


def test_the_census_is_generated_from_a_release_not_typed(tmp_path):
    """Runs the generator against the tracked smoke release and checks
    the block it produces, so the generator is exercised on every run
    rather than only when someone has a full release on disk.

    Counting the smoke release here by hand would be a second
    implementation of the thing under test, so the expected numbers are
    read off the release itself.
    """
    import json
    import shutil

    import scripts.write_kg_report as W

    rel = ROOT / "data" / "releases" / "smoke"
    assert (rel / "nodes.jsonl").exists(), "the smoke release is tracked"
    report = tmp_path / "r.txt"
    shutil.copy(REPORT, report)
    rc = W.main(["--release", str(rel), "--report", str(report)])
    assert rc == 0
    block = report.read_text().split(W.BEGIN, 1)[1].split(W.END, 1)[0]

    n_nodes = sum(1 for ln in (rel / "nodes.jsonl").read_text().splitlines()
                  if ln.strip())
    n_edges = sum(1 for ln in (rel / "edges.jsonl").read_text().splitlines()
                  if ln.strip())
    assert f"{n_nodes:,} nodes" in block
    assert f"{n_edges:,} edges" in block
    for line in (rel / "edges.jsonl").read_text().splitlines():
        if line.strip():
            assert json.loads(line)["predicate"] in block
    assert "nodes_sha" in block and "edges_sha" in block


def test_the_census_stamp_is_the_release_it_was_built_from(tmp_path):
    """The stamp is the whole protection. A block carrying a hash that
    is not the release's own describes some other graph, which is the
    defect `make check-results` exists to catch everywhere else."""
    import shutil

    import scripts.write_kg_report as W
    from kgav.provenance import release_fingerprint

    rel = ROOT / "data" / "releases" / "smoke"
    report = tmp_path / "r.txt"
    shutil.copy(REPORT, report)
    W.main(["--release", str(rel), "--report", str(report)])
    block = report.read_text().split(W.BEGIN, 1)[1].split(W.END, 1)[0]
    fp = release_fingerprint(rel)
    assert fp["nodes_sha"] in block
    assert fp["edges_sha"] in block


def test_check_mode_fails_on_a_stale_block(tmp_path):
    """--check is what a gate would run. If it passed on a stale block
    it would certify exactly the state it exists to detect."""
    import shutil

    import scripts.write_kg_report as W

    rel = ROOT / "data" / "releases" / "smoke"
    report = tmp_path / "r.txt"
    shutil.copy(REPORT, report)
    assert W.main(["--release", str(rel), "--report", str(report),
                   "--check"]) == 1, "a placeholder block passed --check"
    W.main(["--release", str(rel), "--report", str(report)])
    assert W.main(["--release", str(rel), "--report", str(report),
                   "--check"]) == 0
    # And a hand-edit of a generated block is caught.
    t = report.read_text().replace(" nodes\n", " nodes (approx)\n", 1)
    report.write_text(t)
    assert W.main(["--release", str(rel), "--report", str(report),
                   "--check"]) == 1, "a hand-edited block passed --check"


def test_the_committed_census_matches_the_release_when_one_is_present():
    """If a full release is on disk, the block in the repository must
    describe it. No skip when it is absent -- the assertion below is
    about the placeholder, so this test runs either way."""
    import scripts.write_kg_report as W

    rel = ROOT / "data" / "releases" / "v0.1"
    if (rel / "nodes.jsonl").exists():
        assert W.main(["--release", str(rel), "--check"]) == 0, (
            "docs/kg_audit_report.txt's census does not describe "
            "data/releases/v0.1; run make kg-report")
    else:
        block = REPORT.read_text().split(W.BEGIN, 1)[1].split(W.END, 1)[0]
        assert "not yet generated" in block, (
            "no release is on disk and the census is not marked as "
            "ungenerated, so it describes a graph nobody can check")


# -------------------------------------------------------- architecture
def test_the_module_counts_in_the_architecture_section_are_right(flat):
    """The layout block states how many modules and scripts there are.
    Those move, and a count in prose is exactly what this project has
    repeatedly found stale."""
    mods = sorted((ROOT / "src" / "kgav").glob("*.py"))
    scripts = sorted((ROOT / "scripts").glob("*.py"))
    assert f"{len(mods)} modules" in flat, (
        f"src/kgav holds {len(mods)} modules")
    assert f"{len(scripts)} entry points" in flat, (
        f"scripts/ holds {len(scripts)} entry points")


def test_every_agent_module_is_named_in_the_architecture_section(flat):
    """The serve-side pipeline is described stage by stage. A module
    added to the agent and not to that description is a stage nobody
    outside the project knows runs."""
    agent = sorted(p.stem for p in (ROOT / "src" / "kgav" / "agent")
                   .glob("*.py") if p.stem != "__init__")
    block = section(flat, "ARCHITECTURE", until="SCHEMA SPECIFICATION")
    for mod in agent:
        assert f"{mod}.py" in block, (
            f"src/kgav/agent/{mod}.py is not described in the "
            f"architecture section")


def test_the_three_stages_are_named_in_the_dataflow_itself(flat):
    """SCOPED TO THE DIAGRAM, not the section.

    The fifth time in this project a membership check was satisfied by
    an unrelated occurrence of the same string. "INGEST -> ASSEMBLE ->
    SERVE" appears in the sentence above the diagram, so renaming the
    ASSEMBLE box inside it passed. The diagram is what a reader follows;
    the sentence is a summary of it.
    """
    block = section(flat, "ARCHITECTURE", until="SCHEMA SPECIFICATION")
    diagram = block.split("no stage reaches past the one before it.", 1)[1] \
        .split("WHY LAYERS ARE INDEPENDENT", 1)[0]
    for stage in ("INGEST", "ASSEMBLE", "SERVE", "SCORE"):
        assert stage in diagram, (
            f"{stage} is not a labelled box in the data-flow diagram")
    assert "data/raw" in diagram and "data/releases" in diagram


def test_the_directory_layout_describes_each_tree(flat):
    """Also scoped: data/results/ appears in the diagram above, so a
    layout entry could be renamed and still be "present" in the
    section."""
    block = section(flat, "ARCHITECTURE", until="SCHEMA SPECIFICATION")
    layout = block.split("DIRECTORY LAYOUT", 1)[1]
    for tree in ("config/", "src/kgav/", "src/kgav/agent/", "scripts/",
                 "tests/", "data/raw/", "data/releases/", "data/results/",
                 "docs/"):
        assert tree in layout, f"{tree} is not in the directory layout"
    # data/results is tracked on purpose and the layout must say why.
    assert "TRACKED" in layout


def test_every_declared_class_and_predicate_appears_in_a_census(tmp_path):
    """Presence AND absence. Listing only what the release holds reports
    presence and stays silent on absence -- the thing this project
    refuses everywhere else ("a source this ingest cannot read is not a
    source with no data"). v0.1 holds 7 of 11 node classes and 14 of 24
    predicates, and the first census showed neither fact."""
    import shutil

    import scripts.write_kg_report as W
    rel = ROOT / "data" / "releases" / "smoke"
    report = tmp_path / "r.txt"
    shutil.copy(REPORT, report)
    W.main(["--release", str(rel), "--report", str(report)])
    block = report.read_text().split(W.BEGIN, 1)[1].split(W.END, 1)[0]
    s = __import__("kgav.schema", fromlist=["load_schema"]).load_schema()
    for name in s.node_classes:
        assert name in block, f"node class {name} missing from the census"
    for e in s.edge_classes:
        assert e.predicate in block, f"{e.predicate} missing from the census"
    assert "declared classes are empty" in block
    assert "declared predicates are empty" in block


def test_an_absent_evaluation_label_is_called_out_as_one(tmp_path):
    """TREATS and IN_TRIAL_FOR are both empty on v0.1, so every
    evaluation in this project rests on HAS_ANTIVIRAL_ACTIVITY_AGAINST
    alone. A census that lists them among twenty other empty predicates
    buries the one fact that changes how the results read."""
    import shutil

    import scripts.write_kg_report as W
    rel = ROOT / "data" / "releases" / "smoke"
    report = tmp_path / "r.txt"
    shutil.copy(REPORT, report)
    W.main(["--release", str(rel), "--report", str(report)])
    block = report.read_text().split(W.BEGIN, 1)[1].split(W.END, 1)[0]
    s = __import__("kgav.schema", fromlist=["load_schema"]).load_schema()
    labels = s.evaluation_label_predicates()
    assert "EVALUATION LABEL" in block, (
        "an absent label predicate is not marked as a label")
    assert f"of {len(labels)} evaluation label" in block


def test_the_census_totals_reconcile(tmp_path):
    """The per-class and per-predicate counts must sum to the totals. If
    they ever do not, the census is quietly describing a subset and the
    block says so rather than printing a figure that looks whole."""
    import shutil

    import scripts.write_kg_report as W
    rel = ROOT / "data" / "releases" / "smoke"
    report = tmp_path / "r.txt"
    shutil.copy(REPORT, report)
    W.main(["--release", str(rel), "--report", str(report)])
    block = report.read_text().split(W.BEGIN, 1)[1].split(W.END, 1)[0]
    assert "!!" not in block, f"the census failed to reconcile:\n{block}"


def undescribed_sources(text: str) -> list[str]:
    """infores ids the census shows carrying edges that the prose never
    mentions. Factored out so it can be tested on a synthetic block as
    well as on the committed one."""
    import re as _re

    import scripts.write_kg_report as W
    block = text.split(W.BEGIN, 1)[1].split(W.END, 1)[0]
    prose = text.split(W.BEGIN, 1)[0] + text.split(W.END, 1)[1]
    ids = {i.rstrip(".") for i in
           _re.findall(r"infores:[a-z0-9][a-z0-9._-]*", block)}
    return sorted(i for i in ids if i not in prose)


def test_the_source_check_catches_a_source_the_report_never_names():
    """Exercised on a synthetic census, so it runs whether or not a
    release has been assembled here. A test that returns early when its
    input is a placeholder is a test that reports green having checked
    nothing, which is the defect this suite keeps finding in itself."""
    import scripts.write_kg_report as W
    fake = (f"infores:uniprot is described in the prose.\n{W.BEGIN}\n"
            f"      infores:uniprot    10\n"
            f"      infores:nosuchdb    1\n{W.END}\ntail\n")
    assert undescribed_sources(fake) == ["infores:nosuchdb"]
    ok = fake.replace("tail", "infores:nosuchdb is described here.")
    assert undescribed_sources(ok) == []


def test_every_source_in_the_census_is_described_in_the_report(text):
    """THE GAP THIS CLOSES. The sibling test reads infores ids out of
    the SOURCE CODE, and several -- infores:dip, infores:mint,
    infores:hpidb -- are derived at runtime from a curator's own name in
    MITAB column 12, so they are never literals anywhere in the tree.
    v0.1 carries edges from all three while the report said those
    endpoints returned nothing usable. The census is the only place
    their names appear, so it is what the report must agree with.
    """
    import scripts.write_kg_report as W
    block = text.split(W.BEGIN, 1)[1].split(W.END, 1)[0]
    if "not yet generated" in block:
        # Not silent: the placeholder is itself asserted, so this test
        # can never pass by finding nothing to look at.
        assert "run `make kg-report`" in block
        return
    missing = undescribed_sources(text)
    assert not missing, (
        f"these sources carry edges in the release and are not described "
        f"anywhere in the report: {missing}")


def test_the_census_counts_curators_not_source_strings(tmp_path):
    """primary_knowledge_source is a "|" list once an edge is merged, so
    counting raw values reports "infores:intact|infores:virhostnet" as a
    source of its own. v0.1 showed 35 such values for 13 real curators,
    and VirHostNet read 5,222 when it curated 10,163."""
    import json
    import shutil

    import scripts.write_kg_report as W
    rel = tmp_path / "rel"
    rel.mkdir()
    (rel / "nodes.jsonl").write_text(json.dumps(
        {"id": "A", "class": "Protein", "properties": {}}) + "\n")
    e = {"subject": "A", "predicate": "PHYSICALLY_INTERACTS_WITH",
         "object": "B", "evidence_tier": 1,
         "first_asserted_date": "2020-01-01"}
    rows = [dict(e, primary_knowledge_source=s) for s in (
        "infores:virhostnet",
        "infores:intact|infores:virhostnet",
        "infores:intact|infores:virhostnet|infores:virhostnet")]
    (rel / "edges.jsonl").write_text(
        "\n".join(json.dumps(r) for r in rows))
    report = tmp_path / "r.txt"
    shutil.copy(REPORT, report)
    W.main(["--release", str(rel), "--report", str(report)])
    block = report.read_text().split(W.BEGIN, 1)[1].split(W.END, 1)[0]
    assert "EDGES BY CURATOR" in block
    # virhostnet on all three, intact on two -- not three "sources".
    assert re.search(r"infores:virhostnet\s+3\b", block), block
    assert re.search(r"infores:intact\s+2\b", block), block
    assert "infores:intact|infores:virhostnet" not in block
    assert "2 edge(s) cite more than one curator" in block


def test_a_repeated_curator_is_reported_as_the_defect_it_is(tmp_path):
    """542 edges in v0.1 carry the same curator twice in one value --
    the footprint of merge_edges unioning strings. A census that merely
    deduplicated them for display would hide that the release predates
    the fix."""
    import json
    import shutil

    import scripts.write_kg_report as W
    rel = tmp_path / "rel"
    rel.mkdir()
    (rel / "nodes.jsonl").write_text(json.dumps(
        {"id": "A", "class": "Protein", "properties": {}}) + "\n")
    (rel / "edges.jsonl").write_text(json.dumps(
        {"subject": "A", "predicate": "PHYSICALLY_INTERACTS_WITH",
         "object": "B", "evidence_tier": 1,
         "first_asserted_date": "2020-01-01",
         "primary_knowledge_source":
             "infores:intact|infores:virhostnet|infores:virhostnet"}))
    report = tmp_path / "r.txt"
    shutil.copy(REPORT, report)
    W.main(["--release", str(rel), "--report", str(report)])
    block = report.read_text().split(W.BEGIN, 1)[1].split(W.END, 1)[0]
    assert "REPEATED curator" in block
    assert "assembled before the fix" in block
