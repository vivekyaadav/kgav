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
    body = flat.split("5. METAPATHS", 1)[1].split("HUB EXCLUSION", 1)[0]
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
    body = flat.split("5. METAPATHS", 1)[1].split("HUB EXCLUSION", 1)[0]
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
    body = text.split("4. BUILD PIPELINE", 1)[1].split("WHY hosttargets",
                                                       1)[0]
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
    block = flat.split("7. PER-VIRUS STATE", 1)[1].split("8. KNOWN GAPS",
                                                         1)[0]
    for v in doc["viruses"]:
        assert v["name"] in block, (
            f"{v['name']} is in the register and is in no tier in "
            f"section 7")


def test_the_tier_totals_add_up_to_the_register(flat):
    """(8) + (5) + (2) + (1) must be the register's size, or one virus is
    counted twice or not at all."""
    doc = yaml.safe_load(REGISTER.read_text())
    block = flat.split("7. PER-VIRUS STATE", 1)[1].split("8. KNOWN GAPS",
                                                         1)[0]
    totals = [int(n) for n in re.findall(r"\((\d+)\) —", block)]
    assert len(totals) == 4, f"expected four tiers, found {totals}"
    assert sum(totals) == len(doc["viruses"]), (
        f"tiers sum to {sum(totals)} and the register holds "
        f"{len(doc['viruses'])}")


# ------------------------------------------------------------------ claims
def test_the_report_does_not_restate_a_release_node_or_edge_total(text):
    """The one rule it sets for itself. A total quoted here is pinned to
    one release and goes stale silently; `make audit` is the live
    figure. This is the defect that produced two reports describing a
    graph that no longer existed."""
    bad = re.findall(r"\b\d{1,3},\d{3}\s+(?:nodes|edges)\b", text)
    assert not bad, f"the report restates a release total: {bad}"


def test_the_report_points_at_the_live_audit(flat):
    assert "make audit" in flat
    assert "make check-results" in flat


def test_the_retracted_finding_is_recorded_as_retracted(flat):
    """The whole-graph degree matching was circular and its uniform
    depletion was an artifact. A report that simply omitted it would
    read as though the error never happened."""
    assert "circular" in flat
    assert "retracted" in flat
