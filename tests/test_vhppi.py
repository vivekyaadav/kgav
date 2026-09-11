"""Day 5 gate: the bridge layer joins, and alias resolution is what makes it."""
import json
from pathlib import Path

import pytest
import yaml

from kgav.aliases import AliasIndex
from kgav.emit import Emit
from kgav.schema import load_schema
from kgav.vhppi import ingest_mitab, parse_identifier, parse_row

ROOT = Path(__file__).resolve().parents[1]
VHPPI = ROOT / "data" / "raw" / "vhppi"
SPINE = ROOT / "data" / "releases" / "v0.1-spine"
needs_data = pytest.mark.skipif(not (VHPPI.exists() and SPINE.exists()),
                                reason="vhppi sources or spine release missing")

NSP5 = "UniProtKB:PRO_0000449623"
NSP5_PP1A = "UniProtKB:PRO_0000449639"      # pp1a id for the same protein
ORF7A = "UniProtKB:P0DTC7"
ORF7A_CHAIN = "UniProtKB:PRO_0000449654"
ACE2, TMPRSS2 = "UniProtKB:Q9BYF1", "UniProtKB:O15393"


def _row(a, b, method, pmid, ta, tb, itype, score):
    f = ["-"] * 15
    f[0], f[1] = a, b
    f[6] = f'psi-mi:"{method}"(x)'
    f[8] = f"pubmed:{pmid}"
    f[9], f[10] = f"taxid:{ta}", f"taxid:{tb}"
    f[11] = f'psi-mi:"{itype}"(y)'
    f[14] = f"virhostnet-miscore:{score}"
    return "\t".join(f)


@pytest.fixture
def world(tmp_path):
    def prot(pid, viral, **extra):
        p = {"taxon_id": "NCBITaxon:2697049" if viral else "NCBITaxon:9606",
             "is_viral": viral, "sequence_hash": pid[-4:], "reviewed": True}
        p.update(extra)
        return {"id": pid, "class": "Protein", "properties": p}

    spine_nodes = [
        {"id": "NCBITaxon:2697049", "class": "OrganismTaxon",
         "properties": {"label": "SARS-CoV-2", "family": "Coronaviridae",
                        "baltimore_class": "IV", "is_enveloped": True}},
        prot(NSP5, True, mature_peptide="3C-like proteinase", protein_family="nsp5"),
        prot(NSP5_PP1A, True), prot(ORF7A, True), prot(ORF7A_CHAIN, True),
    ]
    same_as = lambda s, o, rule: {
        "subject": s, "predicate": "SAME_AS", "object": o,
        "qualifiers": {"merge_rule": rule},
        "primary_knowledge_source": "infores:uniprot",
        "evidence_tier": 1, "first_asserted_date": "2020-04-22"}
    spine_edges = [same_as(NSP5_PP1A, NSP5, "shared_chain_between_polyproteins"),
                   same_as(ORF7A_CHAIN, ORF7A, "single_chain_equals_parent_protein")]
    host_nodes = [prot(ACE2, False), prot(TMPRSS2, False)]

    for name, nodes, edges in [("spine", spine_nodes, spine_edges),
                               ("host", host_nodes, [])]:
        d = tmp_path / name
        d.mkdir()
        (d / "nodes.jsonl").write_text("\n".join(json.dumps(x) for x in nodes))
        (d / "edges.jsonl").write_text("\n".join(json.dumps(x) for x in edges))

    rows = [
        _row(f"uniprotkb:P0DTD1-{NSP5_PP1A.split(':')[1]}", f"uniprotkb:{ACE2.split(':')[1]}",
             "MI:1314", "32838362", 2697049, 9606, "MI:0915", "0.51"),
        _row(f"uniprotkb:P0DTC7-{ORF7A_CHAIN.split(':')[1]}", f"uniprotkb:{TMPRSS2.split(':')[1]}",
             "MI:0018", "33110215", 2697049, 9606, "MI:0407", "0.72"),
        _row(f"uniprotkb:P0DTD1-{NSP5_PP1A.split(':')[1]}", f"uniprotkb:{ACE2.split(':')[1]}",
             "MI:1314", "32838362", 2697049, 9606, "MI:0915", "0.51"),      # duplicate
        _row(f"uniprotkb:{ACE2.split(':')[1]}", f"uniprotkb:{TMPRSS2.split(':')[1]}",
             "MI:0400", "1234", 9606, 9606, "MI:0915", "0.3"),              # host-host
        _row(f"uniprotkb:{NSP5.split(':')[1]}", "uniprotkb:P0DTC7",
             "MI:0400", "1235", 2697049, 2697049, "MI:0915", "0.3"),        # viral-viral
        _row("uniprotkb:NOPE9", f"uniprotkb:{ACE2.split(':')[1]}",
             "MI:0400", "1236", 2697049, 9606, "MI:0915", "0.3"),           # unresolved
        _row(f"uniprotkb:P0DTD1-{NSP5_PP1A.split(':')[1]}", f"uniprotkb:{ACE2.split(':')[1]}",
             "MI:0400", "1237", 99999, 9606, "MI:0915", "0.4"),             # unmapped taxon
    ]
    (tmp_path / "vhn_test.tab27").write_text("\n".join(rows))
    (tmp_path / "viruses.yaml").write_text(
        yaml.safe_dump({"viruses": [{"taxon": "2697049", "name": "SARS-CoV-2"}]}))
    return tmp_path


@pytest.fixture
def built(world):
    index = AliasIndex.from_releases(world / "spine", world / "host")
    em = Emit()
    stats = ingest_mitab(em, [world / "vhn_test.tab27"], index,
                         {"2697049": "2697049"}, "infores:virhostnet", "2024-01-01")
    return em, stats, index


# ------------------------------------------------------------------- parsing
def test_chain_reference_resolves_to_the_chain_not_the_parent():
    """The source is saying WHICH mature peptide it measured. Collapsing to the
    polyprotein throws that away and loses all target specificity in M3."""
    assert parse_identifier("uniprotkb:P0DTD1-PRO_0000449624") == "UniProtKB:PRO_0000449624"
    assert parse_identifier("uniprotkb:P0DTD1") == "UniProtKB:P0DTD1"
    assert parse_identifier("-") is None


def test_row_parsing_extracts_psi_mi_codes():
    r = parse_row(_row("uniprotkb:A", "uniprotkb:B", "MI:0018", "999",
                       2697049, 9606, "MI:0407", "0.9"))
    assert r["detection_method"] == "MI:0018"
    assert r["interaction_type"] == "MI:0407"
    assert r["pmid"] == "999" and r["score"] == "0.9"


# ------------------------------------------------------------------- aliases
def test_alias_index_follows_same_as(world):
    ix = AliasIndex.from_releases(world / "spine", world / "host")
    assert ix.resolve(NSP5_PP1A) == NSP5
    assert ix.resolve(ORF7A_CHAIN) == ORF7A
    assert ix.resolve(NSP5) == NSP5
    assert ix.resolve("UniProtKB:NOPE") is None


def test_alias_cycle_returns_none(tmp_path):
    """A cycle is a data bug. Looping forever would hide it."""
    d = tmp_path / "r"
    d.mkdir()
    nodes = [{"id": f"X:{i}", "class": "Protein",
              "properties": {"taxon_id": "NCBITaxon:9606", "is_viral": False,
                             "sequence_hash": "z", "reviewed": True}} for i in (1, 2)]
    edges = [{"subject": "X:1", "predicate": "SAME_AS", "object": "X:2",
              "qualifiers": {"merge_rule": "r"}, "primary_knowledge_source": "s",
              "evidence_tier": 1, "first_asserted_date": "2020-01-01"},
             {"subject": "X:2", "predicate": "SAME_AS", "object": "X:1",
              "qualifiers": {"merge_rule": "r"}, "primary_knowledge_source": "s",
              "evidence_tier": 1, "first_asserted_date": "2020-01-01"}]
    (d / "nodes.jsonl").write_text("\n".join(json.dumps(x) for x in nodes))
    (d / "edges.jsonl").write_text("\n".join(json.dumps(x) for x in edges))
    assert AliasIndex.from_releases(d).resolve("X:1") is None


# -------------------------------------------------------------------- bridge
def test_aliased_endpoints_join(built):
    """The whole point: without alias resolution both of these rows drop."""
    em, stats, _ = built
    assert stats["edges"] == 2
    objects = {e["object"] for e in em.edges}
    assert objects == {NSP5, ORF7A}


def test_edges_are_host_to_viral(built):
    em, _, index = built
    for e in em.edges:
        assert index.is_viral(e["subject"]) is False
        assert index.is_viral(e["object"]) is True


def test_same_species_pairs_are_excluded(built):
    """host-host belongs to the STRING pass; viral-viral is out of scope."""
    _, stats, _ = built
    assert stats["host_host"] == 1 and stats["viral_viral"] == 1


def test_duplicates_unresolved_and_unmapped_taxa_are_counted(built):
    _, stats, _ = built
    assert stats["duplicate"] == 1
    assert stats["unresolved_id"] == 1
    assert stats["unmapped_viral_taxon"] == 1


def test_assay_distinction_is_preserved_not_averaged(built):
    """Proximity and binary evidence must stay distinguishable at query time."""
    em, stats, _ = built
    assert stats["direct"] == 1 and stats["co_complex"] == 1
    # interaction_type is deliberately absent: VirHostNet reports MI:0915 for
    # every row regardless of assay, so storing it would fabricate a
    # distinction the source does not make.
    assert all("interaction_type" not in e["qualifiers"] for e in em.edges)
    methods = {e["qualifiers"]["detection_method"] for e in em.edges}
    assert methods == {"MI:1314", "MI:0018"}


def test_all_bridge_edges_are_tier_one(built):
    """Tier means how the fact was established, not how good the assay was.
    BioID is curated experimental; demoting it would corrupt the Day 12
    tier ablation."""
    em, _, _ = built
    assert {e["evidence_tier"] for e in em.edges} == {1}


def test_publications_and_score_are_carried(built):
    em, _, _ = built
    for e in em.edges:
        assert e["publications"] and e["publications"][0].startswith("PMID:")
        assert isinstance(e["source_score"], float)


def test_output_validates_against_schema(built):
    em, _, index = built
    violations = load_schema().validate_batch(index.nodes.values(), em.edges)
    assert violations == [], [str(v) for v in violations[:8]]


# ----------------------------------------------------------------- real data
@needs_data
def test_real_bridge_resolves_nearly_everything(tmp_path):
    index = AliasIndex.from_releases(SPINE, ROOT / "data" / "releases" / "v0.1-host")
    cfg = yaml.safe_load((ROOT / "config" / "viruses.yaml").read_text())
    tmap = {}
    for v in cfg["viruses"]:
        tmap[str(v["taxon"])] = str(v["taxon"])
        for iso in v.get("isolate_taxa") or []:
            tmap[str(iso)] = str(v["taxon"])

    em = Emit()
    stats = ingest_mitab(em, sorted(VHPPI.glob("vhn_*.tab27")), index, tmap,
                         "infores:virhostnet", "2024-01-01")
    assert stats["edges"] > 3_000, f"only {stats['edges']} bridge edges"
    # Unresolved ids are the silent-failure mode this layer exists to avoid.
    assert stats["unresolved_id"] / max(stats["rows"], 1) < 0.25
