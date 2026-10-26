"""Gate for the Guide to PHARMACOLOGY layer.

Built against a published column list rather than a downloaded file, so
the failure that matters is a release whose columns have shifted being
read anyway -- a field in the wrong place becoming a direction is worse
than no layer at all.
"""
from __future__ import annotations

import sys
from importlib import import_module
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from kgav.emit import Emit
from kgav.gtop import (
    ACTION_TO_DIRECTION,
    GTOP_SOURCE,
    INTERACTION_COLUMNS,
    GtopFormatError,
    _rows,
    direction_of,
    ingest_interactions,
    ligand_inchikeys,
)

IHEAD = ("ligand_id,target_uniprot,target_species,action,pubmed_id,"
         "original_affinity_median_nm")


def _inter(**kw):
    base = {"ligand_id": "1", "target_uniprot": "P04637",
            "target_species": "Human", "action": "Inhibitor",
            "pubmed_id": "12345678", "original_affinity_median_nm": "40"}
    base.update(kw)
    return base


KEYS = {"1": "AAAAAAAAAAAAAA-BBBBBBBBBB-C"}
PROTS = {"UniProtKB:P04637"}
CMPDS = {"INCHIKEY:AAAAAAAAAAAAAA-BBBBBBBBBB-C"}


# ------------------------------------------------------------- the columns
def test_a_file_missing_a_used_column_is_refused_by_name(tmp_path):
    """Real files drift from docs. A shifted field read as a direction is
    worse than no layer, so the reader names what it looked for."""
    p = tmp_path / "interactions.csv"
    p.write_text("#GtP v2026.1\nligand_id,target_uniprot\n1,P04637\n")
    with pytest.raises(GtopFormatError, match="missing column"):
        _rows(p, INTERACTION_COLUMNS)


def test_the_leading_version_comment_is_skipped(tmp_path):
    p = tmp_path / "interactions.csv"
    p.write_text("#GtP v2026.1 2026-01-01\n" + IHEAD +
                 "\n1,P04637,Human,Inhibitor,12345678,40\n")
    rows = _rows(p, INTERACTION_COLUMNS)
    assert len(rows) == 1 and rows[0]["target_uniprot"] == "P04637"


def test_a_ligand_file_with_no_inchikey_column_is_refused(tmp_path):
    p = tmp_path / "ligands.csv"
    p.write_text("ligand_id,name\n1,aspirin\n")
    with pytest.raises(GtopFormatError, match="no InChIKey column"):
        ligand_inchikeys(p)


def test_inchikeys_load_under_either_spelling(tmp_path):
    for col in ("inchikey", "InChIKey"):
        p = tmp_path / f"l_{col}.csv"
        p.write_text(f"ligand_id,{col}\n1,KEY-1\n2,\n")
        got = ligand_inchikeys(p)
        assert got == {"1": "KEY-1"}, col


# ----------------------------------------------------------- the direction
@pytest.mark.parametrize("action,want", [
    ("Inhibitor", "inhibitor"), ("Antagonist", "inhibitor"),
    ("Channel blocker", "inhibitor"), ("Agonist", "agonist"),
    ("Partial agonist", "agonist"), ("Activator", "agonist"),
    ("Allosteric modulator", "modulator"),
])
def test_known_actions_map_into_the_schemas_enum(action, want):
    assert direction_of(action) == want


@pytest.mark.parametrize("action", ["Binding", "Unknown", "", "None",
                                    "Competitive", None])
def test_a_binder_is_not_recorded_as_an_inhibitor(action):
    """THE POINT OF THE LAYER IS AN HONEST DIRECTION. Recording a binder
    as an inhibitor would manufacture the mechanistic claim this exists to
    supply properly."""
    assert direction_of(action) == "unknown"


def test_the_direction_vocabulary_is_the_schemas(tmp_path):
    """One vocabulary, so a direction means the same thing whichever
    source wrote it. host_targets.py populates the same enum from ChEMBL."""
    import yaml
    enum = yaml.safe_load(
        (ROOT / "schema" / "kg_schema.yaml").read_text()
    )["enums"]["target_direction"]
    assert set(ACTION_TO_DIRECTION.values()) <= set(enum)
    assert "unknown" in enum


# ------------------------------------------------------------- the emission
def test_a_human_interaction_becomes_a_targets_edge():
    em = Emit()
    s = ingest_interactions(em, [_inter()], KEYS, PROTS, CMPDS)
    assert s["edges"] == 1
    e = em.edges[0]
    assert e["subject"] == "INCHIKEY:AAAAAAAAAAAAAA-BBBBBBBBBB-C"
    assert e["predicate"] == "TARGETS"
    assert e["object"] == "UniProtKB:P04637"
    assert e["qualifiers"]["direction"] == "inhibitor"
    assert e["qualifiers"]["affinity_nm"] == 40.0
    assert e["publications"] == ["PMID:12345678"]
    assert e["primary_knowledge_source"] == GTOP_SOURCE


def test_non_human_rows_are_dropped():
    """GtP curates rat, mouse and others alongside human. A rat target is
    not a human one and the host layer holds only human proteins."""
    em = Emit()
    s = ingest_interactions(em, [_inter(target_species="Rat")], KEYS,
                            PROTS, CMPDS)
    assert em.edges == [] and s["not_human"] == 1


def test_a_complex_accession_list_is_dropped_not_split():
    """A pipe-separated accession means a protein COMPLEX. Taking the
    first would assert the ligand targets one subunit, which the source
    does not say."""
    em = Emit()
    s = ingest_interactions(em, [_inter(target_uniprot="P04637|Q00987")],
                            KEYS, PROTS, CMPDS)
    assert em.edges == [] and s["multi_accession"] == 1


def test_a_ligand_without_an_inchikey_is_skipped_and_counted():
    """Without one it cannot be matched to a compound already in the
    graph, and minting on a name would duplicate a drug already there."""
    em = Emit()
    s = ingest_interactions(em, [_inter(ligand_id="99")], KEYS, PROTS, CMPDS)
    assert em.edges == [] and s["ligand_without_inchikey"] == 1


def test_a_compound_not_in_the_graph_is_skipped_by_default():
    em = Emit()
    s = ingest_interactions(em, [_inter()], KEYS, PROTS, known_compounds=set())
    assert em.edges == [] and s["compound_not_in_graph"] == 1


def test_minting_is_possible_but_off_by_default():
    """Passing known_compounds=None adds ligands the graph lacks. It is
    opt-in because a node keyed on an InChIKey the chemistry layer never
    normalised sits beside its own duplicate."""
    em = Emit()
    s = ingest_interactions(em, [_inter()], KEYS, PROTS, known_compounds=None)
    assert s["edges"] == 1


def test_a_protein_not_in_the_graph_is_skipped():
    em = Emit()
    s = ingest_interactions(em, [_inter()], KEYS, set(), CMPDS)
    assert em.edges == [] and s["protein_not_in_graph"] == 1


def test_the_same_pair_twice_is_one_edge():
    em = Emit()
    s = ingest_interactions(em, [_inter(), _inter(action="Antagonist")],
                            KEYS, PROTS, CMPDS)
    assert s["edges"] == 1 and s["duplicate"] == 1


def test_an_unparseable_affinity_is_counted_not_crashed():
    em = Emit()
    s = ingest_interactions(em, [_inter(original_affinity_median_nm="n/a")],
                            KEYS, PROTS, CMPDS)
    assert s["edges"] == 1
    assert s["unparseable_affinity"] == 1
    assert "affinity_nm" not in em.edges[0]["qualifiers"]


# ---------------------------------------------------------------- placement
def test_gtop_outranks_hosttargets_because_only_it_knows_direction():
    from kgav.assemble import LAYER_PRECEDENCE
    assert LAYER_PRECEDENCE.index("gtop") < LAYER_PRECEDENCE.index("hosttargets")


def test_the_driver_refuses_without_the_files(tmp_path, capsys):
    D = import_module("ingest_gtop")
    assert D.main(["--raw", str(tmp_path / "none")]) == 1
    assert "guidetopharmacology.org/download" in capsys.readouterr().out


def test_the_driver_writes_nothing_when_no_edge_survives(tmp_path, capsys):
    """An empty layer directory assembles silently and looks like
    coverage."""
    D = import_module("ingest_gtop")
    raw = tmp_path / "raw"
    raw.mkdir()
    (raw / "interactions.csv").write_text(
        "#v\n" + IHEAD + "\n1,P04637,Rat,Inhibitor,1,40\n")
    (raw / "ligands.csv").write_text("ligand_id,inchikey\n1,KEY-1\n")
    rel = tmp_path / "rel"
    rel.mkdir()
    import json as _json
    (rel / "nodes.jsonl").write_text(_json.dumps(
        {"id": "UniProtKB:P04637", "class": "Protein",
         "properties": {"is_viral": False}}))
    out = tmp_path / "out"
    assert D.main(["--raw", str(raw), "--release", str(rel),
                   "--out", str(out)]) == 1
    assert not out.exists()
    assert "would assemble silently" in capsys.readouterr().out


# --------------------------------------------- the whole table, not a sample
# Sign words first: a term carrying a direction is that direction even when
# it also names a mechanism, which is how host_targets reads "negative
# allosteric modulator" as an inhibitor. Mechanism words only decide when
# nothing signed the term.
def _expected_direction(action: str) -> str:
    a = action.upper()
    if "ANTAGONIST" in a or "INHIBIT" in a or "BLOCKER" in a or "NEGATIVE" in a:
        return "inhibitor"
    if "AGONIST" in a or "ACTIVATOR" in a or "POSITIVE" in a:
        return "agonist"
    if "MODULATOR" in a or "ENHANCER" in a or a == "MIXED":
        return "modulator"
    return "unknown"


@pytest.mark.parametrize("action", sorted(ACTION_TO_DIRECTION))
def test_every_action_in_the_table_means_what_it_says(action):
    """Derived from the table, so a key cannot be added without a test.

    The hand-listed cases above covered "Allosteric modulator" and missed
    the bare "MODULATOR" key, plus MIXED, ALLOSTERIC ENHANCER, NEGATIVE,
    POSITIVE, BIASED AGONIST, FULL AGONIST, GATING INHIBITOR,
    IRREVERSIBLE INHIBITION and ANTISENSE INHIBITION. Mapping MODULATOR to
    "inhibitor" left all thirty tests passing -- and an allosteric
    modulator recorded as an inhibitor is exactly the manufactured
    mechanism this module's docstring says it exists to avoid, now carried
    into an agent that reasons about direction.
    """
    assert ACTION_TO_DIRECTION[action] == _expected_direction(action)


def test_the_consistency_rule_is_not_vacuous():
    """It must disagree with a wrong table, or it is asserting nothing."""
    assert _expected_direction("MODULATOR") == "modulator"
    assert _expected_direction("INHIBITOR") != "modulator"
    assert _expected_direction("NEGATIVE ALLOSTERIC MODULATOR") == "inhibitor"


def test_gtop_and_host_targets_agree_where_their_vocabularies_meet():
    """One direction vocabulary across sources, or a direction means
    different things depending on who wrote it."""
    from kgav.host_targets import ACTION_TO_DIRECTION as HT
    shared = set(ACTION_TO_DIRECTION) & set(HT)
    assert shared, "the two tables must overlap or this proves nothing"
    for k in sorted(shared):
        assert ACTION_TO_DIRECTION[k] == HT[k], k


def test_a_non_numeric_pubmed_id_is_dropped():
    """PMID:abc is not a citation. The isdigit filter had no test, so
    removing it changed nothing."""
    em = Emit()
    rows = [_inter(pubmed_id="12345678|not-a-pmid|  |87654321")]
    ingest_interactions(em, rows, {"1": "KEYAAA"})
    assert em.edges[0]["publications"] == ["PMID:12345678", "PMID:87654321"]


def test_an_edge_with_no_usable_pmid_carries_none_rather_than_junk():
    em = Emit()
    ingest_interactions(em, [_inter(pubmed_id="n/a")], {"1": "KEYAAA"})
    assert not (em.edges[0].get("publications") or [])
