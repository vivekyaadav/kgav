"""Gate for curated host-factor screens.

The layer exists so a published screen can be added by transcription
rather than by new code, and the thing that would ruin it is a row that
becomes a WRONG edge instead of an error. src/kgav/orcs.py's header
records why: a sign error here produces a confidently wrong
recommendation rather than a missing one.
"""
from __future__ import annotations

import sys
from importlib import import_module
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from kgav.curated import (
    CURATED_SOURCE,
    CURATED_TIER,
    CuratedRowError,
    ingest_curated,
    read_table,
    validate_row,
)
from kgav.emit import Emit

HEAD = "ncbi_gene_id\ttaxon\tdirection\tpmid\tcell_line\tscreen_type\teffect_size"


def _row(**kw):
    base = {"ncbi_gene_id": "7157", "taxon": "11060",
            "direction": "dependency", "pmid": "27545878",
            "cell_line": "Huh7.5.1", "screen_type": "CRISPR KO",
            "effect_size": "2.4"}
    base.update(kw)
    return base


# ------------------------------------------------------------- validation
def test_a_good_row_normalises():
    got = validate_row(_row(), 2)
    assert got["gene"] == "NCBIGene:7157"
    assert got["taxon"] == "11060"
    assert got["direction"] == "dependency"
    assert got["pmid"] == "PMID:27545878"
    assert got["effect_size"] == 2.4


def test_a_missing_direction_is_refused_not_defaulted():
    """THE RULE THE LAYER EXISTS FOR. ORCS infers a sign because screen
    conventions disagree; a curator read the paper and knows. Defaulting
    would file half a susceptibility screen backwards in silence."""
    with pytest.raises(CuratedRowError, match="direction is empty"):
        validate_row(_row(direction=""), 2)
    with pytest.raises(CuratedRowError, match="never inferred"):
        validate_row(_row(direction="unknown"), 2)
    with pytest.raises(CuratedRowError, match="never inferred"):
        validate_row(_row(direction="proviral"), 7)


def test_a_gene_symbol_is_refused_because_it_will_not_map():
    """build_gene_to_proteins keys on NCBIGene CURIEs. A symbol silently
    maps to nothing, so the row would vanish into gene_not_in_graph and
    look like missing biology."""
    with pytest.raises(CuratedRowError, match="will not map"):
        validate_row(_row(ncbi_gene_id="TP53"), 3)


@pytest.mark.parametrize("field,value,msg", [
    ("taxon", "NCBITaxon:11060", "bare NCBI taxon id"),
    ("pmid", "Marceau et al", "not numeric"),
    ("effect_size", "strong", "not a number"),
    ("pmid", "", "pmid is empty"),
])
def test_every_malformed_field_names_its_row(field, value, msg):
    with pytest.raises(CuratedRowError, match=msg):
        validate_row(_row(**{field: value}), 9)


def test_a_pmid_prefix_is_accepted_and_normalised():
    assert validate_row(_row(pmid="PMID:27545878"), 2)["pmid"] == "PMID:27545878"


def test_a_table_missing_a_required_column_is_refused_whole(tmp_path):
    p = tmp_path / "bad.tsv"
    p.write_text("ncbi_gene_id\ttaxon\tpmid\n7157\t11060\t27545878\n")
    with pytest.raises(CuratedRowError, match="direction"):
        read_table(p)


def test_a_well_formed_table_reads(tmp_path):
    p = tmp_path / "ok.tsv"
    p.write_text(HEAD + "\n7157\t11060\tdependency\t27545878\tHuh7.5.1\tKO\t2.4\n")
    rows = read_table(p)
    assert len(rows) == 1
    assert validate_row(rows[0], 2)["gene"] == "NCBIGene:7157"


# ---------------------------------------------------------------- emission
def test_edges_carry_their_own_source_and_a_weaker_tier():
    """Hand-transcribed from a PDF by one person with no second reader is
    weaker provenance than a database export, and anything downstream must
    be able to tell them apart and discount them."""
    em = Emit()
    ingest_curated(em, [validate_row(_row(), 2)],
                   {"NCBIGene:7157": ["UniProtKB:P04637"]}, known_taxa={"11060"})
    assert len(em.edges) == 1
    e = em.edges[0]
    assert e["primary_knowledge_source"] == CURATED_SOURCE != "infores:biogrid-orcs"
    assert e["evidence_tier"] == CURATED_TIER == 2
    assert e["predicate"] == "HOST_FACTOR_FOR"
    assert e["object"] == "NCBITaxon:11060"
    assert e["qualifiers"]["direction"] == "dependency"
    assert e["publications"] == ["PMID:27545878"]


def test_a_row_for_an_unbuilt_taxon_is_skipped_not_dangled(tmp_path):
    """An edge pointing at a taxon the spine never built validates as a
    dangling reference much later, and reads as a graph fault rather than
    a curation one. DENV-1, DENV-3, DENV-4, JEV and TBEV sit in
    `surveyed:` precisely like this."""
    em = Emit()
    s = ingest_curated(em, [validate_row(_row(taxon="11053"), 2)],
                       {"NCBIGene:7157": ["UniProtKB:P04637"]},
                       known_taxa={"11060"})
    assert em.edges == []
    assert s["taxon_not_built"] == 1


def test_a_gene_absent_from_the_graph_is_counted(tmp_path):
    em = Emit()
    s = ingest_curated(em, [validate_row(_row(), 2)], {}, known_taxa={"11060"})
    assert em.edges == []
    assert s["gene_not_in_graph"] == 1


def test_one_gene_with_several_proteins_gets_the_edge_on_each():
    """The screen knocked out the gene, so the claim applies to every
    product -- the same rule build_gene_to_proteins states."""
    em = Emit()
    ingest_curated(em, [validate_row(_row(), 2)],
                   {"NCBIGene:7157": ["UniProtKB:A", "UniProtKB:B"]},
                   known_taxa={"11060"})
    assert {e["subject"] for e in em.edges} == {"UniProtKB:A", "UniProtKB:B"}


def test_the_same_hit_from_one_paper_is_not_double_counted():
    em = Emit()
    rows = [validate_row(_row(), 2), validate_row(_row(), 3)]
    s = ingest_curated(em, rows, {"NCBIGene:7157": ["UniProtKB:P04637"]},
                       known_taxa={"11060"})
    assert len(em.edges) == 1
    assert s["duplicate"] == 1


def test_two_papers_reporting_the_same_hit_both_land():
    """Independent replication is evidence, not duplication, so the key
    includes the PMID."""
    em = Emit()
    rows = [validate_row(_row(), 2), validate_row(_row(pmid="30578351"), 3)]
    ingest_curated(em, rows, {"NCBIGene:7157": ["UniProtKB:P04637"]},
                   known_taxa={"11060"})
    assert len(em.edges) == 2


def test_opposite_directions_from_one_paper_both_land():
    """A screen can report a gene as restricting in one arm and depending
    in another; collapsing them would silently drop one claim."""
    em = Emit()
    rows = [validate_row(_row(), 2),
            validate_row(_row(direction="restriction"), 3)]
    ingest_curated(em, rows, {"NCBIGene:7157": ["UniProtKB:P04637"]},
                   known_taxa={"11060"})
    assert len(em.edges) == 2


# ------------------------------------------------------------------ driver
def test_the_driver_refuses_a_bad_table_and_writes_nothing(tmp_path, capsys):
    D = import_module("ingest_curated")
    tables = tmp_path / "t"
    tables.mkdir()
    (tables / "bad.tsv").write_text(
        HEAD + "\n7157\t11060\t\t27545878\tHuh7.5.1\tKO\t2.4\n")
    out = tmp_path / "out"
    assert D.main(["--tables", str(tables), "--out", str(out)]) == 1
    assert not out.exists()
    printed = capsys.readouterr().out
    assert "direction is empty" in printed
    assert "REFUSED ROWS ABOVE" in printed


def test_an_empty_directory_is_not_an_error(tmp_path, capsys):
    D = import_module("ingest_curated")
    assert D.main(["--tables", str(tmp_path / "none")]) == 1
    assert "not an error" in capsys.readouterr().out
