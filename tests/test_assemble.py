"""Day 8 gate: cross-layer problems are caught, not merged away."""
import json
from pathlib import Path

import pytest

from kgav.assemble import (
    LAYER_PRECEDENCE,
    assemble,
    connectivity_report,
    merge_edges,
    orphan_nodes,
)
from kgav.schema import load_schema

ROOT = Path(__file__).resolve().parents[1]
ASSEMBLED = ROOT / "data" / "releases" / "v0.1"
needs_release = pytest.mark.skipif(not (ASSEMBLED / "nodes.jsonl").exists(),
                                   reason="assembled release not built")

NSP5 = "UniProtKB:PRO_1"
ACE2 = "UniProtKB:Q9BYF1"
DRUG = "INCHIKEY:AAAAAAAAAAAAAA-BBBBBBBBBB-C"


def _prot(pid, viral, **extra):
    p = {"taxon_id": "NCBITaxon:2697049" if viral else "NCBITaxon:9606",
         "is_viral": viral, "sequence_hash": "h", "reviewed": True}
    p.update(extra)
    return {"id": pid, "class": "Protein", "properties": p}


def _edge(s, p, o, src, date, tier=1, q=None, pmids=None):
    e = {"subject": s, "predicate": p, "object": o, "qualifiers": q or {},
         "primary_knowledge_source": src, "evidence_tier": tier,
         "first_asserted_date": date}
    if pmids:
        e["publications"] = pmids
    return e


@pytest.fixture
def layers(tmp_path):
    def write(name, nodes, edges):
        d = tmp_path / f"v0.1-{name}"
        d.mkdir(parents=True)
        (d / "nodes.jsonl").write_text("\n".join(json.dumps(x) for x in nodes))
        (d / "edges.jsonl").write_text("\n".join(json.dumps(x) for x in edges))
        return d

    virus = {"id": "NCBITaxon:2697049", "class": "OrganismTaxon",
             "properties": {"label": "SARS-CoV-2", "family": "Coronaviridae",
                            "baltimore_class": "IV", "is_enveloped": True}}
    out = {}
    out["spine"] = write("spine",
                         [virus, _prot(NSP5, True, protein_family="nsp5",
                                       mature_peptide="3CLpro")], [])
    out["host"] = write("host",
                        [_prot(ACE2, False, gene_symbol="ACE2"),
                         _prot(NSP5, True, sequence_hash="DIFFERENT")],
                        [_edge(ACE2, "PHYSICALLY_INTERACTS_WITH", NSP5,
                               "infores:string", "2020-01-01",
                               q={"detection_method": "string:x"})])
    out["vhppi"] = write("vhppi", [],
                         [_edge(ACE2, "PHYSICALLY_INTERACTS_WITH", NSP5,
                                "infores:virhostnet", "2019-05-05",
                                q={"detection_method": "MI:1314"}, pmids=["PMID:1"])])
    out["chembl"] = write("chembl",
                          [{"id": DRUG, "class": "SmallMolecule",
                            "properties": {"smiles": "CC", "inchikey_skel": "A" * 14,
                                           "is_approved": True, "salt_collapsed": False,
                                           "stereo_collapsed": False}}],
                          [_edge(DRUG, "INHIBITS", NSP5, "infores:chembl", "2021-03-04",
                                 q={"assay_type": "biochemical", "ic50_nm": 25.0,
                                    "unquantified": True})])
    return out


@pytest.fixture
def built(layers):
    return assemble(layers)


# ------------------------------------------------------------ node conflicts
def test_property_conflict_is_reported_not_resolved_silently(built):
    """Two layers giving the same node different sequence_hash values is a
    normalisation bug. Picking a winner quietly hides it."""
    assert built.stats["property_conflict:sequence_hash"] == 1
    assert any("sequence_hash" in c for c in built.conflicts)


def test_earlier_layer_wins_on_conflict(built):
    """The spine is hand-curated; chembl knows least about what it did not
    create. Precedence must not depend on dict ordering."""
    assert built.nodes[NSP5]["properties"]["sequence_hash"] == "h"


def test_gaps_are_filled_from_later_layers(built):
    """Non-conflicting properties merge: the host layer adds what the spine
    omitted."""
    assert built.nodes[ACE2]["properties"]["gene_symbol"] == "ACE2"


def test_layer_order_is_deterministic(layers):
    """Passing the dict in a different order must not change the result."""
    a = assemble(layers)
    b = assemble({k: layers[k] for k in reversed(list(layers))})
    assert a.nodes[NSP5]["properties"]["sequence_hash"] == \
        b.nodes[NSP5]["properties"]["sequence_hash"]
    assert set(a.edges) == set(b.edges)


def test_precedence_covers_every_layer():
    assert set(LAYER_PRECEDENCE) == {"spine", "host", "vhppi", "orcs", "chembl",
                                     "hosttargets"}


# -------------------------------------------------------------- edge merging
def test_same_fact_from_two_sources_becomes_one_edge(built):
    """Stacking them would double the fact's weight in every path count."""
    ppi = [e for k, e in built.edges.items() if k[1] == "PHYSICALLY_INTERACTS_WITH"]
    assert len(ppi) == 1
    assert built.stats["edge_merged"] == 1


def test_merged_edge_keeps_both_sources_and_all_publications(built):
    e = built.edges[(ACE2, "PHYSICALLY_INTERACTS_WITH", NSP5)]
    assert e["primary_knowledge_source"] == "infores:string|infores:virhostnet"
    assert e["publications"] == ["PMID:1"]


def test_merged_edge_takes_the_earliest_date(built):
    """A temporal split must use when a fact was FIRST established, not when
    the second source repeated it. Getting this wrong leaks the future."""
    e = built.edges[(ACE2, "PHYSICALLY_INTERACTS_WITH", NSP5)]
    assert e["first_asserted_date"] == "2019-05-05"


def test_merged_edge_takes_the_strongest_tier():
    a = _edge("A", "P", "B", "s1", "2020-01-01", tier=4)
    b = _edge("A", "P", "B", "s2", "2020-01-01", tier=1)
    assert merge_edges(a, b)["evidence_tier"] == 1


def test_merge_handles_a_missing_date():
    a = _edge("A", "P", "B", "s1", "2020-01-01")
    b = dict(_edge("A", "P", "B", "s2", "2019-01-01"), first_asserted_date=None)
    assert merge_edges(a, b)["first_asserted_date"] == "2020-01-01"


# -------------------------------------------------------------- graph checks
def test_no_dangling_references_after_assembly(built):
    """Per-layer, edges-only releases dangle by design. After assembly a
    dangling edge is a join that failed."""
    assert built.dangling() == []


def test_dangling_is_detected_when_present(layers, tmp_path):
    d = tmp_path / "v0.1-broken"
    d.mkdir()
    (d / "nodes.jsonl").write_text("")
    (d / "edges.jsonl").write_text(json.dumps(
        _edge("UniProtKB:NOPE", "INHIBITS", NSP5, "s", "2020-01-01",
              q={"assay_type": "biochemical"})))
    a = assemble({**layers, "broken": d})
    assert len(a.dangling()) == 1


def test_orphan_nodes_are_counted_by_class(built):
    orphans = orphan_nodes(built)
    assert orphans["OrganismTaxon"] == 1
    assert orphans["SmallMolecule"] == 0


def test_connectivity_report_shows_endpoint_classes(built):
    rep = connectivity_report(built)
    assert rep["INHIBITS"]["SmallMolecule->Protein"] == 1


def test_assembled_graph_validates(built):
    violations = load_schema().validate_batch(built.nodes.values(), built.edges.values())
    assert violations == [], [str(v) for v in violations[:6]]


def test_manifest_records_degree_and_merges(built):
    m = built.manifest({}, "0.3.1")
    assert m["nodes"] == len(built.nodes)
    assert set(m["degree"]) == {"median", "p99", "p99_9", "max"}
    assert "property_conflict:sequence_hash" in m["merges"]


# ----------------------------------------------------------- real assembly
@needs_release
def test_real_assembled_release_validates():
    nodes = [json.loads(x) for x in (ASSEMBLED / "nodes.jsonl").read_text().splitlines() if x.strip()]
    edges = [json.loads(x) for x in (ASSEMBLED / "edges.jsonl").read_text().splitlines() if x.strip()]
    assert len(nodes) > 50_000 and len(edges) > 300_000
    ids = {n["id"] for n in nodes}
    dangling = [e for e in edges if e["subject"] not in ids or e["object"] not in ids]
    assert not dangling, f"{len(dangling)} dangling edges in the assembled release"


# ------------------------------------------------- identity-defining qualifiers
IDENTITY = {"HOST_FACTOR_FOR": ["direction"]}


def _hf(protein, virus, direction, source, date):
    return _edge(protein, "HOST_FACTOR_FOR", virus, source, date,
                 q={"direction": direction, "screen_type": "CRISPRko",
                    "cell_line": "CVCL_0574"})


def test_opposite_directions_stay_separate_facts(tmp_path):
    """A dependency edge and a restriction edge for the same (protein, virus)
    are different facts. Merging on the triple alone let whichever source was
    read first decide the direction -- 320 ORCS pairs are contested, ACE2
    among them at 34 dependency to 5 restriction."""
    d = tmp_path / "v0.1-orcs"
    d.mkdir()
    (d / "nodes.jsonl").write_text("")
    (d / "edges.jsonl").write_text("\n".join(json.dumps(x) for x in [
        _hf(ACE2, "NCBITaxon:2697049", "dependency", "infores:biogrid-orcs", "2020-01-01"),
        _hf(ACE2, "NCBITaxon:2697049", "restriction", "infores:biogrid-orcs", "2021-01-01"),
    ]))
    a = assemble({"orcs": d}, IDENTITY)
    hf = [k for k in a.edges if k[1] == "HOST_FACTOR_FOR"]
    assert len(hf) == 2, "opposite directions collapsed into one fact"
    assert {k[3] for k in hf} == {"dependency", "restriction"}


def test_same_direction_from_many_screens_records_support_count(tmp_path):
    """Fifteen screens agreeing is stronger evidence than one. Without a count
    DWPC weights a lone noisy hit exactly like ACE2."""
    d = tmp_path / "v0.1-orcs"
    d.mkdir()
    (d / "nodes.jsonl").write_text("")
    (d / "edges.jsonl").write_text("\n".join(json.dumps(
        _hf(ACE2, "NCBITaxon:2697049", "dependency", "infores:biogrid-orcs",
            f"20{20+i}-01-01")) for i in range(4)))
    a = assemble({"orcs": d}, IDENTITY)
    e = a.edges[(ACE2, "HOST_FACTOR_FOR", "NCBITaxon:2697049", "dependency")]
    assert e["support_count"] == 4
    assert e["first_asserted_date"] == "2020-01-01"


def test_identity_qualifiers_do_not_affect_other_predicates(built):
    """INHIBITS declares none, so its key stays the plain triple."""
    assert (DRUG, "INHIBITS", NSP5) in built.edges


def test_support_count_defaults_to_one(built):
    e = built.edges[(DRUG, "INHIBITS", NSP5)]
    assert e.get("support_count", 1) == 1
