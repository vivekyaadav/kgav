"""Gate for the Guide to PHARMACOLOGY layer.

A field in the wrong place becoming a direction is worse than no layer at
all, so the reader refuses a file whose columns are not what it expects.

THAT GUARD WAS THE ONLY THING THAT WORKED. Until 2026-11-01 both the
reader and every fixture here used snake_case column names taken from the
published column list -- `ligand_id`, `target_uniprot`, `action` -- and
GtoPdb serves `"Ligand ID"`, `"Target UniProt ID"`, `"Action"`. The tests
matched the reader, the reader matched the documentation, and nothing
matched the data. 286 lines of them passed; the first real download was
refused on the first column check.

So the fixtures here are the download's own header, verbatim. A test that
invents its subject cannot find out that the subject is wrong.
"""
from __future__ import annotations

import csv
import sys
from importlib import import_module
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from kgav.emit import Emit
from kgav.gtop import (
    ACTION_TO_DIRECTION,
    C_TYPE,
    GTOP_SOURCE,
    INTERACTION_COLUMNS,
    SKELETON_LEN,
    TYPE_TO_DIRECTION,
    CompoundIndex,
    GtopFormatError,
    _rows,
    direction_of,
    file_version,
    ingest_interactions,
    ligand_chembl_ids,
    ligand_inchikeys,
)

#: interactions.csv's header, VERBATIM from the GtoPdb 2026.3.1 download,
#: and its quoted version line.
#:
#: NOT A HAND-WRITTEN SUBSET, which is what this file had before: snake_case
#: names invented from the published column list. Every fixture matched the
#: reader, the reader matched the documentation, and the documentation
#: matched nothing the site serves -- so 286 lines of tests passed while the
#: reader could not open a real file. A fixture taken from the data is the
#: only thing that could have caught it.
#:
#: Kept in full rather than trimmed to the six columns used, because the
#: pair that matters is in the tail: `Type` and `Ligand Type` both exist and
#: mean different things, and a loose column match would take either.
GTOP_VERSION_LINE = '"# GtoPdb Version: 2026.3.1 - published: 2026-10-12"'

IHEAD = (
    '"Target","Target ID","Target Subunit IDs","Target Gene Symbol",'
    '"Target UniProt ID","Target Ensembl Gene ID","Target Ligand",'
    '"Target Ligand ID","Target Ligand Subunit IDs",'
    '"Target Ligand Gene Symbol","Target Ligand UniProt ID",'
    '"Target Ligand Ensembl Gene ID","Target Ligand PubChem SID",'
    '"Target Species","Ligand ID","Ligand","Ligand Type",'
    '"Ligand Subunit IDs","Ligand Gene Symbol","Ligand Species",'
    '"Ligand PubChem SID","Approved","Type","Action","Action comment",'
    '"Selectivity","Endogenous","Primary Target","concentration Range",'
    '"Affinity Units","Affinity High","Affinity Median","Affinity Low",'
    '"Original Affinity Units","Original Affinity Low nm",'
    '"Original Affinity Median nm","Original Affinity High nm",'
    '"Original Affinity Relation","Assay Description","Receptor Site",'
    '"Ligand Context","PubMed ID","Webpage URLs","Patent Numbers"'
)

LHEAD = (
    '"Ligand ID","Name","Species","Type","Approved","Withdrawn","Labelled",'
    '"Radioactive","PubChem SID","PubChem CID","UniProt ID","Ensembl ID",'
    '"ChEMBL ID","Ligand Subunit IDs","Ligand Subunit Name",'
    '"Ligand Subunit UniProt IDs","Ligand Subunit Ensembl IDs",'
    '"IUPAC name","INN","Synonyms","SMILES","InChIKey","InChI",'
    '"GtoImmuPdb","GtoMPdb","Antibacterial"'
)


def _row_text(head: str, **vals) -> str:
    """One CSV row under `head`, with named columns filled and the rest
    blank -- so a fixture names the columns it means instead of counting
    commas, and stays correct when GtoP adds a column."""
    cols = next(csv.reader([head]))
    return ",".join(f'"{vals.get(c, "")}"' for c in cols)


def _inter(**kw):
    """One interaction row, keyed by the real column names.

    UNKNOWN KEYS RAISE. When the columns were renamed, call sites still
    passed `target_species="Rat"` -- a key nothing reads, so the row stayed
    Human and the test asserted a drop that had not happened. A fixture
    that silently ignores an argument is a test that passes for the wrong
    reason, so the names are checked against the download's own header.
    """
    base = {"Ligand ID": "1", "Target UniProt ID": "P04637",
            "Target Species": "Human", "Type": "Inhibitor",
            "Action": "Inhibition", "PubMed ID": "12345678",
            "Original Affinity Median nm": "40"}
    bad = sorted(set(kw) - set(next(csv.reader([IHEAD]))))
    if bad:
        raise KeyError(
            f"{bad} are not columns of interactions.csv. Pass the real "
            f"names, e.g. _inter(**{{'Target Species': 'Rat'}}).")
    base.update(kw)
    return base


KEYS = {"1": "AAAAAAAAAAAAAA-BBBBBBBBBB-C"}
PROTS = {"UniProtKB:P04637"}
CMPDS = {"INCHIKEY:AAAAAAAAAAAAAA-BBBBBBBBBB-C"}
IDX = CompoundIndex(frozenset(CMPDS))


# ------------------------------------------------------------- the columns
def test_the_real_download_header_is_accepted(tmp_path):
    """THE TEST THAT WAS MISSING. Everything else here could pass against
    an invented header; this one cannot."""
    p = tmp_path / "interactions.csv"
    p.write_text(GTOP_VERSION_LINE + "\n" + IHEAD + "\n"
                 + _row_text(IHEAD, **_inter()) + "\n")
    rows = _rows(p, INTERACTION_COLUMNS)
    assert len(rows) == 1
    assert rows[0]["Target UniProt ID"] == "P04637"
    assert rows[0]["Type"] == "Inhibitor"


def test_the_quoted_version_comment_is_skipped_and_kept(tmp_path):
    """GtoPdb QUOTES its comment line, so it does not start with '#'. The
    old filter was `not l.startswith("#")` and the old fixture wrote an
    unquoted '#GtP v2026.1' -- the test passed and the reader could not
    open a real export. The version is returned too: DATA_LICENSES.md
    wants the release in the MANIFEST, and this is the release."""
    p = tmp_path / "interactions.csv"
    p.write_text(GTOP_VERSION_LINE + "\n" + IHEAD + "\n"
                 + _row_text(IHEAD, **_inter()) + "\n")
    assert file_version(p) == "GtoPdb Version: 2026.3.1 - published: 2026-10-12"
    assert len(_rows(p, INTERACTION_COLUMNS)) == 1


def test_a_file_with_no_comment_line_keeps_its_first_row(tmp_path):
    """Only the FIRST line is a comment. Filtering every line that looks
    like one would silently drop data if an export ever stopped quoting."""
    p = tmp_path / "interactions.csv"
    p.write_text(IHEAD + "\n" + _row_text(IHEAD, **_inter()) + "\n")
    assert file_version(p) is None
    assert len(_rows(p, INTERACTION_COLUMNS)) == 1


def test_a_file_missing_a_used_column_is_refused_by_name(tmp_path):
    """Real files drift from docs. A shifted field read as a direction is
    worse than no layer, so the reader names what it looked for."""
    p = tmp_path / "interactions.csv"
    p.write_text(GTOP_VERSION_LINE + '\n"Ligand ID","Target UniProt ID"\n'
                 '"1","P04637"\n')
    with pytest.raises(GtopFormatError, match="missing column"):
        _rows(p, INTERACTION_COLUMNS)


def test_the_snake_case_names_the_reader_used_to_expect_are_refused(tmp_path):
    """The format that was assumed for weeks. If this ever passes, the
    column constants have drifted back to the documentation."""
    p = tmp_path / "interactions.csv"
    p.write_text("ligand_id,target_uniprot,target_species,action,pubmed_id\n"
                 "1,P04637,Human,Inhibition,12345678\n")
    with pytest.raises(GtopFormatError, match="missing column"):
        _rows(p, INTERACTION_COLUMNS)


def test_the_interaction_type_is_read_and_not_the_ligand_type():
    """interactions.csv has BOTH `Type` and `Ligand Type`. `Type` is the
    interaction; `Ligand Type` is the molecule class, which says nothing
    about direction. Reading the wrong one is a whole-layer error that
    looks like a working layer."""
    assert C_TYPE == "Type"
    assert "Ligand Type" in IHEAD and '"Type"' in IHEAD
    em = Emit()
    ingest_interactions(em, [_inter(**{"Ligand Type": "Antibody"})],
                        KEYS, PROTS, IDX)
    assert em.edges[0]["qualifiers"]["direction"] == "inhibitor"


def test_a_ligand_file_with_no_inchikey_column_is_refused(tmp_path):
    p = tmp_path / "ligands.csv"
    p.write_text('"Ligand ID","Name"\n"1","aspirin"\n')
    with pytest.raises(GtopFormatError, match="no InChIKey column"):
        ligand_inchikeys(p)


def test_inchikeys_load_under_either_spelling(tmp_path):
    for col in ("inchikey", "InChIKey"):
        p = tmp_path / f"l_{col}.csv"
        p.write_text(f'"Ligand ID","{col}"\n"1","KEY-1"\n"2",""\n')
        assert ligand_inchikeys(p) == {"1": "KEY-1"}, col


def test_the_real_ligand_header_yields_its_inchikeys(tmp_path):
    p = tmp_path / "ligands.csv"
    p.write_text(GTOP_VERSION_LINE + "\n" + LHEAD + "\n"
                 + _row_text(LHEAD, **{"Ligand ID": "1",
                                       "InChIKey": "KEY-1"}) + "\n")
    assert ligand_inchikeys(p) == {"1": "KEY-1"}


# ----------------------------------------------------------- the direction
@pytest.mark.parametrize("type_value,want", [
    ("Inhibitor", "inhibitor"), ("Antagonist", "inhibitor"),
    ("Channel blocker", "inhibitor"), ("Gating inhibitor", "inhibitor"),
    ("Agonist", "agonist"), ("Activator", "agonist"),
    ("Allosteric modulator", "modulator"),
])
def test_every_type_value_maps_into_the_schemas_enum(type_value, want):
    assert direction_of(type_value) == want


@pytest.mark.parametrize("type_value", ["None", "Antibody", "Fusion protein",
                                        "Subunit-specific", "", None])
def test_a_type_that_states_no_direction_is_unknown(type_value):
    """1,600 rows of 24,968 in 2026.3.1. An antibody may antagonise or may
    only bind, and recording a guess would manufacture the mechanistic
    claim this layer exists to supply properly."""
    assert direction_of(type_value) == "unknown"


def test_the_action_column_alone_decides_nothing():
    """`Action` is read to refuse, never to assign. Were it assigning, its
    largest value -- `Inhibition`, 8,159 rows -- would have to map, and
    the bug this replaces was exactly that it did not."""
    assert direction_of("", "Inhibition") == "unknown"
    assert direction_of("None", "Inhibition") == "unknown"


@pytest.mark.parametrize("type_value,action,want", [
    ("Antagonist", "Partial agonist", "unknown"),   # 77 rows in 2026.3.1
    ("Agonist", "Inverse agonist", "unknown"),      # 59 rows
    ("Inhibitor", "Inhibition", "inhibitor"),       # agreement
    ("Agonist", "Full agonist", "agonist"),
    ("Inhibitor", "Binding", "inhibitor"),          # Action uninformative
])
def test_a_row_that_contradicts_itself_asserts_no_direction(
        type_value, action, want):
    """Taking `Type` where the same row's `Action` says the opposite would
    assert the opposite of what the source states beside it."""
    assert direction_of(type_value, action) == want


@pytest.mark.parametrize("action", ["Positive", "Negative", "Inhibition",
                                    "Inverse agonist", "Potentiation", ""])
def test_an_allosteric_modulator_stays_a_modulator(action):
    """DECIDED 2026-11-01, and the cost is real. For allosteric rows the
    two columns are complementary -- `Type` is the site, `Action` is the
    direction, Positive 364 rows and Negative 355 -- so this is not a
    contradiction. A positive allosteric modulator has no intrinsic
    activity, so filing it as an agonist would overstate GtoP's claim; the
    price is that M2, M4 and M5 cannot tell a PAM from a NAM."""
    assert direction_of("Allosteric modulator", action) == "modulator"


def test_the_direction_vocabulary_is_the_schemas():
    """One vocabulary, so a direction means the same thing whichever
    source wrote it. host_targets.py populates the same enum from ChEMBL."""
    import yaml
    enum = yaml.safe_load(
        (ROOT / "schema" / "kg_schema.yaml").read_text()
    )["enums"]["target_direction"]
    assert set(TYPE_TO_DIRECTION.values()) <= set(enum)
    assert set(ACTION_TO_DIRECTION.values()) <= set(enum)
    assert "unknown" in enum


# ------------------------------------------------------------- the emission
def test_a_human_interaction_becomes_a_targets_edge():
    em = Emit()
    s = ingest_interactions(em, [_inter()], KEYS, PROTS, IDX)
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
    s = ingest_interactions(em, [_inter(**{"Target Species": "Rat"})], KEYS,
                            PROTS, IDX)
    assert em.edges == [] and s["not_human"] == 1


def test_a_complex_accession_list_is_dropped_not_split():
    """A pipe-separated accession means a protein COMPLEX. Taking the
    first would assert the ligand targets one subunit, which the source
    does not say."""
    em = Emit()
    s = ingest_interactions(em, [_inter(**{"Target UniProt ID": "P04637|Q00987"})],
                            KEYS, PROTS, IDX)
    assert em.edges == [] and s["multi_accession"] == 1


def test_a_ligand_without_an_inchikey_is_skipped_and_counted():
    """Without one it cannot be matched to a compound already in the
    graph, and minting on a name would duplicate a drug already there."""
    em = Emit()
    s = ingest_interactions(em, [_inter(**{"Ligand ID": "99"})], KEYS, PROTS, IDX)
    assert em.edges == [] and s["ligand_without_inchikey"] == 1


def test_a_compound_not_in_the_graph_is_skipped_by_default():
    em = Emit()
    s = ingest_interactions(em, [_inter()], KEYS, PROTS, compounds=CompoundIndex(frozenset()))
    assert em.edges == [] and s["compound_not_in_graph"] == 1


def test_minting_is_possible_but_off_by_default():
    """Passing compounds=None adds ligands the graph lacks. It is
    opt-in because a node keyed on an InChIKey the chemistry layer never
    normalised sits beside its own duplicate."""
    em = Emit()
    s = ingest_interactions(em, [_inter()], KEYS, PROTS, compounds=None)
    assert s["edges"] == 1


def test_a_protein_not_in_the_graph_is_skipped():
    em = Emit()
    s = ingest_interactions(em, [_inter()], KEYS, set(), IDX)
    assert em.edges == [] and s["protein_not_in_graph"] == 1


def test_the_same_pair_twice_is_one_edge():
    em = Emit()
    s = ingest_interactions(em, [_inter(), _inter(Type="Antagonist")],
                            KEYS, PROTS, IDX)
    assert s["edges"] == 1 and s["duplicate"] == 1


def test_an_unparseable_affinity_is_counted_not_crashed():
    em = Emit()
    s = ingest_interactions(em, [_inter(**{"Original Affinity Median nm": "n/a"})],
                            KEYS, PROTS, IDX)
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
        GTOP_VERSION_LINE + "\n" + IHEAD + "\n"
        + _row_text(IHEAD, **_inter(**{"Target Species": "Rat"})) + "\n")
    (raw / "ligands.csv").write_text(
        GTOP_VERSION_LINE + "\n" + LHEAD + "\n"
        + _row_text(LHEAD, **{"Ligand ID": "1", "InChIKey": "KEY-1"}) + "\n")
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
    # INVERSE AGONIST FIRST, and this is not cosmetic ordering: it contains
    # "AGONIST", so the agonist rule would claim it and an inverse agonist
    # reduces constitutive activity -- the opposite. The 2026.3.1 vocabulary
    # also brought ACTIVATION (not ACTIVATOR) and POTENTIATION, neither of
    # which any rule here matched, so both read "unknown" while the table
    # said "agonist".
    if "INVERSE AGONIST" in a:
        return "inhibitor"
    if "ANTAGONIST" in a or "INHIBIT" in a or "BLOCKER" in a or "NEGATIVE" in a:
        return "inhibitor"
    if ("AGONIST" in a or "ACTIVAT" in a or "POSITIVE" in a
            or "POTENTIAT" in a):
        return "agonist"
    if "MODULATOR" in a or "ENHANCER" in a or a == "MIXED":
        return "modulator"
    return "unknown"


@pytest.mark.parametrize("term", sorted(TYPE_TO_DIRECTION))
def test_every_type_in_the_table_means_what_it_says(term):
    """The same independent re-derivation, over the column that decides."""
    assert TYPE_TO_DIRECTION[term] == _expected_direction(term)


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
    rows = [_inter(**{"PubMed ID": "12345678|not-a-pmid|  |87654321"})]
    ingest_interactions(em, rows, {"1": "KEYAAA"})
    assert em.edges[0]["publications"] == ["PMID:12345678", "PMID:87654321"]


def test_an_edge_with_no_usable_pmid_carries_none_rather_than_junk():
    em = Emit()
    ingest_interactions(em, [_inter(**{"PubMed ID": "n/a"})], {"1": "KEYAAA"})
    assert not (em.edges[0].get("publications") or [])


# ------------------------------------------------- the compound join routes
#
# 9,827 human interactions in 2026.3.1 have an InChIKey the graph does not
# hold exactly. 661 of them ARE compounds the graph holds, under a
# different form of the same key, because the chemistry layer salt- and
# stereo-collapses; 24 more match by ChEMBL id. The remaining 9,142 are
# genuinely absent, and minting those is a different decision with a
# different cost -- unlabelled molecules enlarge every `pool` in
# lift = hits * pool / (k * positives) and add to no numerator.
SKEL = "AAAAAAAAAAAAAA"
NODE = "INCHIKEY:AAAAAAAAAAAAAA-BBBBBBBBBB-C"


def _node(nid, skel=None, chembl=None, cls="SmallMolecule"):
    props = {}
    if skel:
        props["inchikey_skel"] = skel
    if chembl:
        props["chembl_id"] = chembl
    return {"id": nid, "class": cls, "properties": props}


def test_an_exact_key_resolves_as_exact():
    idx = CompoundIndex.from_nodes([_node(NODE, SKEL, "CHEMBL1")])
    assert idx.resolve("AAAAAAAAAAAAAA-BBBBBBBBBB-C") == (NODE, "exact")


def test_a_different_salt_or_stereo_form_resolves_by_skeleton():
    """The 661. Same molecule, different key -- the chemistry layer already
    collapsed it, so the graph holds it under the collapsed form."""
    idx = CompoundIndex.from_nodes([_node(NODE, SKEL)])
    assert idx.resolve("AAAAAAAAAAAAAA-ZZZZZZZZZZ-N") == (NODE, "skeleton")


def test_a_ligand_with_no_key_in_the_graph_resolves_by_chembl_id():
    idx = CompoundIndex.from_nodes([_node(NODE, SKEL, "CHEMBL1")])
    got = idx.resolve("QQQQQQQQQQQQQQ-WWWWWWWWWW-N", "CHEMBL1")
    assert got == (NODE, "chembl_id")


def test_the_exact_key_is_preferred_over_the_looser_ones():
    """Order matters: a compound present exactly must never be resolved to
    a different node that merely shares its skeleton."""
    idx = CompoundIndex.from_nodes([_node(NODE, SKEL),
                                    _node("INCHIKEY:OTHER", "OTHERSKELETON")])
    assert idx.resolve("AAAAAAAAAAAAAA-BBBBBBBBBB-C")[1] == "exact"


def test_a_skeleton_shared_by_two_compounds_resolves_to_neither():
    """A SKELETON THAT NAMES TWO COMPOUNDS NAMES NEITHER. The graph holds
    24,051 skeletons over 25,383 compounds, so some are shared, and
    resolving to whichever node was read first would attach an interaction
    to an arbitrary one of two molecules."""
    idx = CompoundIndex.from_nodes([_node("INCHIKEY:ONE", SKEL),
                                    _node("INCHIKEY:TWO", SKEL)])
    assert idx.by_skeleton == {}
    assert idx.resolve("AAAAAAAAAAAAAA-ZZZZZZZZZZ-N") == (None, "absent")


def test_a_chembl_id_shared_by_two_compounds_resolves_to_neither():
    idx = CompoundIndex.from_nodes([_node("INCHIKEY:ONE", None, "CHEMBL1"),
                                    _node("INCHIKEY:TWO", None, "CHEMBL1")])
    assert idx.by_chembl == {}


def test_a_genuinely_absent_compound_stays_absent():
    idx = CompoundIndex.from_nodes([_node(NODE, SKEL, "CHEMBL1")])
    assert idx.resolve("ZZZZZZZZZZZZZZ-YYYYYYYYYY-N", "CHEMBL999") == (
        None, "absent")


def test_an_empty_inchikey_resolves_to_nothing():
    """A blank key's first 14 characters are a blank skeleton, which would
    match any node whose skeleton was also blank."""
    idx = CompoundIndex.from_nodes([_node(NODE, SKEL)])
    assert idx.resolve("") == (None, "absent")
    assert idx.resolve("  ") == (None, "absent")


def test_only_small_molecules_enter_the_index():
    idx = CompoundIndex.from_nodes(
        [_node("UniProtKB:P1", None, None, cls="Protein"), _node(NODE, SKEL)])
    assert idx.exact == frozenset({NODE})


def test_the_skeleton_length_is_the_inchikeys_first_block():
    assert SKELETON_LEN == 14
    first_block = "AAAAAAAAAAAAAA-BBBBBBBBBB-C".partition("-")[0]
    assert len(first_block) == SKELETON_LEN


@pytest.mark.parametrize("key,route", [
    ("AAAAAAAAAAAAAA-ZZZZZZZZZZ-N", "compound_by_skeleton"),
    ("QQQQQQQQQQQQQQ-WWWWWWWWWW-N", "compound_by_chembl_id"),
])
def test_each_looser_route_is_counted_separately(key, route):
    """A join on a looser key has to be auditable. If either count ever
    dwarfs the exact matches, the looser key is matching things it should
    not."""
    idx = CompoundIndex.from_nodes([_node(NODE, SKEL, "CHEMBL1")])
    em = Emit()
    s = ingest_interactions(em, [_inter()], {"1": key}, PROTS, idx,
                            ligand_chembl={"1": "CHEMBL1"})
    assert s["edges"] == 1 and s[route] == 1
    assert em.edges[0]["subject"] == NODE


def test_an_exact_match_is_not_counted_as_a_looser_route():
    idx = CompoundIndex.from_nodes([_node(NODE, SKEL, "CHEMBL1")])
    em = Emit()
    s = ingest_interactions(em, [_inter()], KEYS, PROTS, idx)
    assert s["edges"] == 1
    assert not s["compound_by_skeleton"] and not s["compound_by_chembl_id"]


def test_the_chembl_route_needs_the_ligand_map_to_be_passed():
    """ligand_chembl is optional, so a driver that forgets it must lose the
    route rather than crash -- and the edge count then says so."""
    idx = CompoundIndex.from_nodes([_node(NODE, SKEL, "CHEMBL1")])
    em = Emit()
    s = ingest_interactions(em, [_inter()], {"1": "QQQQQQQQQQQQQQ-WW-N"},
                            PROTS, idx)
    assert s["edges"] == 0 and s["compound_not_in_graph"] == 1


def test_the_chembl_id_column_is_read_from_the_real_ligand_header(tmp_path):
    p = tmp_path / "ligands.csv"
    p.write_text(GTOP_VERSION_LINE + "\n" + LHEAD + "\n"
                 + _row_text(LHEAD, **{"Ligand ID": "1",
                                       "ChEMBL ID": "CHEMBL25"}) + "\n")
    assert ligand_chembl_ids(p) == {"1": "CHEMBL25"}
