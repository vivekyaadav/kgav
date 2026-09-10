"""Day 4 gate: STRING identifiers map correctly, edges are deduped, and the
threshold policy is the one recorded in the module docstring.
"""
import gzip
import json
from collections import Counter
from pathlib import Path

import pytest

from kgav.emit import Emit
from kgav.host import build_ensp_map, ingest_proteome, ingest_string, load_human_proteome
from kgav.schema import load_schema

ROOT = Path(__file__).resolve().parents[1]
HOST = ROOT / "data" / "raw" / "host"
needs_data = pytest.mark.skipif(not HOST.exists(), reason="host sources not downloaded")

PROV = {"proteome_source": "infores:uniprot", "pathway_source": "infores:uniprot",
        "ppi_source": "infores:string", "string_release_date": "2023-08-28"}

HDR = ("protein1 protein2 neighborhood fusion cooccurence coexpression "
       "experimental database textmining combined_score")


@pytest.fixture
def fake_host(tmp_path):
    entries = []
    for acc, sym, gid in [("Q9BYF1", "ACE2", "59272"), ("O15393", "TMPRSS2", "7113"),
                          ("P08183", "ABCB1", "5243")]:
        entries.append({
            "primaryAccession": acc, "entryType": "UniProtKB reviewed (Swiss-Prot)",
            "proteinDescription": {"recommendedName": {"fullName": {"value": f"{sym} protein"}}},
            "genes": [{"geneName": {"value": sym}}], "sequence": {"value": "MKT" + acc},
            "entryAudit": {"firstPublicDate": "2000-05-01"},
            "uniProtKBCrossReferences": [
                {"database": "GeneID", "id": gid},
                {"database": "Reactome", "id": "R-HSA-9679191",
                 "properties": [{"key": "PathwayName", "value": "Potential therapeutics for SARS"}]},
                {"database": "GO", "id": "GO:0006508",
                 "properties": [{"key": "GoTerm", "value": "P:proteolysis"}]},
                {"database": "GO", "id": "GO:0005886",
                 "properties": [{"key": "GoTerm", "value": "C:plasma membrane"}]},
            ]})
    (tmp_path / "human_proteome.json.gz").write_bytes(
        gzip.compress(json.dumps({"results": entries}).encode()))

    aliases = ["9606.ENSP1\tQ9BYF1\tUniProt_AC", "9606.ENSP2\tO15393\tUniProt_AC",
               "9606.ENSP3\tP08183\tUniProt_AC", "9606.ENSP3\tQ9BYF1\tUniProt_AC",
               "9606.ENSP9\tZZZZZZ\tUniProt_AC", "9606.ENSP1\tACE2_HUMAN\tUniProt_ID"]
    (tmp_path / "9606.protein.aliases.v12.0.txt.gz").write_bytes(
        gzip.compress("\n".join(aliases).encode()))

    rows = [HDR,
            "9606.ENSP1 9606.ENSP2 0 0 0 100 900 0 200 950",
            "9606.ENSP2 9606.ENSP1 0 0 0 100 900 0 200 950",
            "9606.ENSP1 9606.ENSP3 0 0 0 0 100 0 800 850",
            "9606.ENSP1 9606.ENSP9 0 0 0 0 900 0 0 900",
            "9606.ENSP1 9606.ENSP1 0 0 0 0 900 0 0 900"]
    (tmp_path / "9606.protein.links.detailed.v12.0.txt.gz").write_bytes(
        gzip.compress("\n".join(rows).encode()))
    return tmp_path


@pytest.fixture
def built(fake_host):
    em = Emit()
    pstats = ingest_proteome(em, load_human_proteome(fake_host / "human_proteome.json.gz"), PROV)
    accs = {n["id"].split(":", 1)[1] for n in em.nodes.values() if n["class"] == "Protein"}
    mapping = build_ensp_map(fake_host / "9606.protein.aliases.v12.0.txt.gz", accs, em.notes)
    sstats = ingest_string(em, fake_host / "9606.protein.links.detailed.v12.0.txt.gz",
                           mapping, PROV)
    return em, pstats, sstats, mapping


def test_only_biological_process_go_terms_become_pathways(built):
    """C: (component) and F: (function) are not pathways."""
    em, _pstats, _, _ = built
    labels = {n["properties"]["label"] for n in em.nodes.values() if n["class"] == "Pathway"}
    assert "proteolysis" in labels
    assert "plasma membrane" not in labels


def test_ambiguous_ensp_is_dropped_not_guessed(built):
    """ENSP3 maps to two reviewed accessions; a wrong edge beats no edge only
    if you never have to defend it."""
    em, _, _, mapping = built
    assert "9606.ENSP3" not in mapping
    assert em.notes["ensp_ambiguous_after_filter"] == 1


def test_reciprocal_duplicates_are_collapsed(built):
    """STRING lists A-B and B-A. Without dedup every degree doubles."""
    _, _, sstats, _ = built
    assert sstats["reciprocal_duplicate"] == 1
    assert sstats["edges"] == 1


def test_weak_experimental_evidence_is_filtered(built):
    """combined=850 but experimental=100: passes the conventional cut, fails ours."""
    em, _, _, _ = built
    ppi = [e for e in em.edges if e["predicate"] == "PHYSICALLY_INTERACTS_WITH"]
    endpoints = {frozenset((e["subject"], e["object"])) for e in ppi}
    assert frozenset(("UniProtKB:Q9BYF1", "UniProtKB:P08183")) not in endpoints


def test_unmapped_and_self_loops_are_reported(built):
    _, _, sstats, _ = built
    assert sstats["unmapped_endpoint"] == 1
    assert sstats["self_loop"] == 1


def test_edges_are_canonically_ordered(built):
    em, _, _, _ = built
    for e in em.edges:
        if e["predicate"] == "PHYSICALLY_INTERACTS_WITH":
            assert e["subject"] < e["object"]


def test_channel_scores_are_preserved(built):
    """The combined>=700 variant must be reconstructible without re-ingest."""
    em, _, _, _ = built
    ppi = next(e for e in em.edges if e["predicate"] == "PHYSICALLY_INTERACTS_WITH")
    assert ppi["qualifiers"]["string_score"] == 950.0
    assert "experimental=900" in ppi["qualifiers"]["channel"]


def test_host_proteins_are_not_flagged_viral(built):
    em, _, _, _ = built
    for n in em.nodes.values():
        if n["class"] == "Protein":
            assert n["properties"]["is_viral"] is False


def test_output_validates_against_schema(built):
    em, _, _, _ = built
    violations = load_schema().validate_batch(em.nodes.values(), em.edges)
    assert violations == [], [str(v) for v in violations[:10]]


# ----------------------------------------------------------------- real data
@needs_data
def test_real_host_layer_scale(tmp_path):
    """Sanity bands. A tenfold miss here means the mapping silently broke."""
    em = Emit()
    ingest_proteome(em, load_human_proteome(HOST / "human_proteome.json.gz"), PROV)
    accs = {n["id"].split(":", 1)[1] for n in em.nodes.values() if n["class"] == "Protein"}
    assert 18_000 < len(accs) < 23_000, f"{len(accs)} reviewed human proteins"

    mapping = build_ensp_map(HOST / "9606.protein.aliases.v12.0.txt.gz", accs, em.notes)
    assert len(mapping) > 15_000, f"only {len(mapping)} ENSP mapped — check UniProt_AC parsing"

    sstats = ingest_string(em, HOST / "9606.protein.links.detailed.v12.0.txt.gz", mapping, PROV)
    assert 100_000 < sstats["edges"] < 400_000, f"{sstats['edges']} PPI edges"

    counts = Counter(n["class"] for n in em.nodes.values())
    assert counts["Pathway"] > 10_000, f"{counts['Pathway']} pathways — GO branch may be silently skipping"
