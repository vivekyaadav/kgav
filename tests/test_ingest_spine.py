"""Day 3 gate: every viral protein resolves to a virus, no orphans, no phantom
duplicate mature peptides from the pp1a/pp1ab overlap.
"""
import gzip
import json
import sys
from pathlib import Path

import pytest
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from ingest_spine import build

from kgav.schema import load_schema

ROOT = Path(__file__).resolve().parents[1]
PROTEOMES = ROOT / "data" / "raw" / "proteomes"
CONFIG = ROOT / "config" / "viruses.yaml"
needs_data = pytest.mark.skipif(not PROTEOMES.exists(), reason="proteomes not downloaded")


def _chain(pid, desc, start):
    return {"type": "Chain", "featureId": pid, "description": desc,
            "location": {"start": {"value": start}, "end": {"value": start + 100}}}


@pytest.fixture
def fake_proteome(tmp_path):
    """pp1ab (16 chains) + pp1a (12 overlapping) + spike (1 chain)."""
    pp1ab = [_chain(f"PRO_00004496{10+i}", f"Non-structural protein {i+1} (nsp{i+1})", 1 + i * 180)
             for i in range(16)]
    pp1a = [_chain(f"PRO_00004497{10+i}", f"Non-structural protein {i+1} (nsp{i+1})", 1 + i * 180)
            for i in range(12)]
    entries = [
        {"primaryAccession": "P0DTD1", "entryType": "UniProtKB reviewed (Swiss-Prot)",
         "proteinDescription": {"recommendedName": {"fullName": {"value": "Replicase polyprotein 1ab"}}},
         "genes": [{"geneName": {"value": "rep"}}], "sequence": {"value": "MESLVPGFNEK"},
         "entryAudit": {"firstPublicDate": "2020-04-22"}, "features": pp1ab},
        {"primaryAccession": "P0DTC1", "entryType": "UniProtKB reviewed (Swiss-Prot)",
         "proteinDescription": {"recommendedName": {"fullName": {"value": "Replicase polyprotein 1a"}}},
         "genes": [{"geneName": {"value": "rep"}}], "sequence": {"value": "MESLVPGFNE"},
         "entryAudit": {"firstPublicDate": "2020-04-22"}, "features": pp1a},
        {"primaryAccession": "P0DTC2", "entryType": "UniProtKB reviewed (Swiss-Prot)",
         "proteinDescription": {"recommendedName": {"fullName": {"value": "Spike glycoprotein"}}},
         "genes": [{"geneName": {"value": "S"}}], "sequence": {"value": "MFVFLVLLPLV"},
         "entryAudit": {"firstPublicDate": "2020-04-22"},
         "features": [_chain("PRO_0000449647", "Spike glycoprotein", 1)]},
    ]
    d = tmp_path / "proteomes"
    d.mkdir()
    (d / "UP000464024.json.gz").write_bytes(gzip.compress(json.dumps({"results": entries}).encode()))
    return d


@pytest.fixture
def one_virus_config(tmp_path):
    cfg = yaml.safe_load(CONFIG.read_text())
    cfg["viruses"] = [v for v in cfg["viruses"] if v["taxon"] == "2697049"]
    p = tmp_path / "viruses.yaml"
    p.write_text(yaml.safe_dump(cfg))
    return p


@pytest.fixture
def spine(fake_proteome, one_virus_config, tmp_path):
    return build(one_virus_config, fake_proteome, tmp_path / "out")


def test_shared_chains_are_not_duplicated(spine):
    """pp1a and pp1ab share nsp1-nsp11. Sixteen mature peptides, not 28."""
    mature = [n for n in spine.nodes.values()
              if n["class"] == "Protein" and n["properties"].get("mature_peptide")]
    assert len(mature) == 16, f"got {len(mature)} mature peptides"
    assert spine.notes["chain_shared_between_polyproteins"] == 12


def test_single_chain_entries_become_aliases_not_targets(spine):
    """Spike is one chain: the chain IS the protein, so it must not become a
    second target. But sources cite the chain id (VirHostNet references ORF7a
    and ORF9b that way), so the identifier has to resolve. It exists as an
    alias -- no mature_peptide, no protein_family, SAME_AS to the parent."""
    cid = "UniProtKB:PRO_0000449647"
    assert cid in spine.nodes, "chain id must resolve, or source edges silently drop"
    props = spine.nodes[cid]["properties"]
    assert not props.get("mature_peptide"), "alias must not count as a target"
    assert not props.get("protein_family")
    aliases = {(e["subject"], e["object"]) for e in spine.edges if e["predicate"] == "SAME_AS"}
    assert (cid, "UniProtKB:P0DTC2") in aliases


def test_nsp_family_is_extracted(spine):
    fams = {n["properties"].get("protein_family") for n in spine.nodes.values()}
    assert "nsp5" in fams, "Mpro (nsp5) not labelled — direct-acting targets unreachable"
    assert "nsp12" in fams, "RdRp (nsp12) not labelled"


def test_every_viral_protein_reaches_its_virus(spine):
    """No orphans: each protein has a path to the taxon via its gene."""
    genes = {e["object"] for e in spine.edges if e["predicate"] == "ENCODED_BY"}
    to_taxon = {e["subject"] for e in spine.edges if e["predicate"] == "BELONGS_TO"}
    assert genes and genes <= to_taxon, "gene(s) with no BELONGS_TO edge"


def test_all_proteins_flagged_viral(spine):
    for n in spine.nodes.values():
        if n["class"] == "Protein":
            assert n["properties"]["is_viral"] is True


def test_output_validates_against_schema(spine):
    schema = load_schema()
    violations = schema.validate_batch(spine.nodes.values(), spine.edges)
    assert violations == [], [str(v) for v in violations[:10]]


def test_no_dangling_edges(spine):
    ids = set(spine.nodes)
    for e in spine.edges:
        assert e["subject"] in ids and e["object"] in ids, e


def test_every_edge_has_provenance_and_date(spine):
    for e in spine.edges:
        assert e["primary_knowledge_source"]
        assert e["first_asserted_date"]
        assert e["evidence_tier"] in (1, 2, 3, 4)


def test_build_is_deterministic(fake_proteome, one_virus_config, tmp_path):
    a = build(one_virus_config, fake_proteome, tmp_path / "a")
    b = build(one_virus_config, fake_proteome, tmp_path / "b")
    assert list(a.nodes) == list(b.nodes)
    assert a.edges == b.edges


# ------------------------------------------------------------- real data
@needs_data
def test_real_spine_validates(tmp_path):
    spine = build(CONFIG, PROTEOMES, tmp_path / "real")
    schema = load_schema()
    violations = schema.validate_batch(spine.nodes.values(), spine.edges)
    assert violations == [], [str(v) for v in violations[:10]]


@needs_data
def test_real_spine_has_all_seven_viruses(tmp_path):
    spine = build(CONFIG, PROTEOMES, tmp_path / "real")
    taxa = [n for n in spine.nodes.values() if n["class"] == "OrganismTaxon"]
    assert len(taxa) == 7, f"expected 7 coronaviruses, got {len(taxa)}"


@needs_data
def test_real_spine_exposes_mpro_and_rdrp(tmp_path):
    """The two targets every coronavirus antiviral programme aims at."""
    spine = build(CONFIG, PROTEOMES, tmp_path / "real")
    fams: dict[str, int] = {}
    for n in spine.nodes.values():
        f = n["properties"].get("protein_family")
        if f:
            fams[f] = fams.get(f, 0) + 1
    assert fams.get("nsp5", 0) >= 1, "no Mpro found in any proteome"
    assert fams.get("nsp12", 0) >= 1, "no RdRp found in any proteome"


@needs_data
def test_every_virus_exposes_the_core_nsps(tmp_path):
    """All seven coronaviruses must label the same nsps.

    Alphacoronaviruses name mature peptides by function ("3C-like proteinase")
    where betacoronaviruses use the number. If the alias table regresses, this
    produces no error anywhere -- just a graph where cross-viral metapaths
    cannot transfer, and a leave-one-virus-out evaluation that underperforms
    for reasons nothing reports.
    """
    spine = build(CONFIG, PROTEOMES, tmp_path / "nsp")
    core = {"nsp3", "nsp5", "nsp12", "nsp13", "nsp14", "nsp15", "nsp16"}
    by_taxon: dict[str, set[str]] = {}
    labels: dict[str, str] = {}
    for n in spine.nodes.values():
        if n["class"] == "OrganismTaxon":
            labels[n["id"]] = n["properties"]["label"]
        elif n["class"] == "Protein" and n["properties"].get("protein_family"):
            by_taxon.setdefault(n["properties"]["taxon_id"], set()).add(
                n["properties"]["protein_family"])

    assert len(labels) == 7
    for tid, name in labels.items():
        missing = core - by_taxon.get(tid, set())
        assert not missing, f"{name} missing {sorted(missing)}"
