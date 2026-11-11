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
    AFFINITY_UNIT_TO_QUALIFIER,
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
            "Original Affinity Median nm": "40",
            "Original Affinity Units": "pKi"}
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
    assert e["qualifiers"]["ki_nm"] == 40.0
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
    assert "ki_nm" not in em.edges[0]["qualifiers"]


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


# ------------------------------------------- one pair stated more than once
#
# GtoP states a (compound, target) pair repeatedly -- 987 times on 2026.3.1,
# different papers for one interaction -- and the reader kept whichever row
# the CSV listed first. On 61 of those pairs the rows DISAGREE about
# direction, so 24 edges read `agonist` or `inhibitor` on no better basis
# than file order. Reconciled through assemble.resolve_qualifier, the same
# function that reconciles two layers.
def _two(d1, d2, **kw):
    """Two rows for one pair, differing only in Type (and anything in kw).

    ACTION IS CLEARED. _inter defaults it to "Inhibition", which
    contradicts any Type that is not inhibitor-like -- so direction_of
    returned `unknown` for both rows and these tests passed or failed on
    the single-row contradiction rule instead of the pair rule they are
    about. The interaction itself is tested by
    test_a_row_that_contradicts_itself_asserts_no_direction.
    """
    return [_inter(Type=d1, Action="", **kw),
            _inter(Type=d2, Action="", **kw)]


def test_two_rows_that_disagree_outright_assert_a_contested_direction():
    """THE 24. Both orders, because the point is that order stops
    mattering.

    The value is `contested`, not `unknown`: the first says the sources
    disagree and the second says nobody stated one, and a consumer that
    cannot tell them apart reports a gap in the data where there is a
    contradiction in it."""
    for rows in (_two("Agonist", "Inhibitor"), _two("Inhibitor", "Agonist")):
        em = Emit()
        s = ingest_interactions(em, rows, KEYS, PROTS, IDX)
        assert s["edges"] == 1
        assert em.edges[0]["qualifiers"]["direction"] == "contested"
        assert s["qualifier_irreconcilable:direction"] == 1


def test_a_row_with_no_direction_does_not_erase_one_that_has_it():
    """THE 32. `Antibody` maps to unknown; a real direction on the other row
    must survive it whichever came first."""
    for rows in (_two("Agonist", "Antibody"), _two("Antibody", "Agonist")):
        em = Emit()
        ingest_interactions(em, rows, KEYS, PROTS, IDX)
        assert em.edges[0]["qualifiers"]["direction"] == "agonist"


def test_a_graded_direction_does_not_replace_a_specific_one():
    for rows in (_two("Inhibitor", "Allosteric modulator"),
                 _two("Allosteric modulator", "Inhibitor")):
        em = Emit()
        ingest_interactions(em, rows, KEYS, PROTS, IDX)
        assert em.edges[0]["qualifiers"]["direction"] == "inhibitor"


def test_the_direction_is_the_same_whatever_order_the_rows_arrive():
    """The whole reason for reconciling rather than first-wins: the export
    lists rows in whatever order it likes, and a release built twice must
    not state two different mechanisms."""
    import itertools
    seen = set()
    for perm in itertools.permutations(["Agonist", "Inhibitor", "Antibody"]):
        em = Emit()
        ingest_interactions(em, [_inter(Type=t, Action="") for t in perm],
                            KEYS, PROTS, IDX)
        seen.add(em.edges[0]["qualifiers"]["direction"])
    assert len(seen) == 1, seen


def test_publications_from_every_row_survive():
    """987 rows' PMIDs were being discarded with the rows. The evidence for
    one interaction is every paper that reported it."""
    em = Emit()
    rows = [_inter(**{"PubMed ID": "111"}), _inter(**{"PubMed ID": "222"})]
    s = ingest_interactions(em, rows, KEYS, PROTS, IDX)
    assert s["edges"] == 1 and s["duplicate"] == 1
    assert em.edges[0]["publications"] == ["PMID:111", "PMID:222"]


def test_the_affinity_is_the_most_potent_reported_not_the_first_listed():
    """A stated rule rather than row order. One qualifier cannot hold
    several papers' measurements; the publications say which papers."""
    for order in (["90", "40"], ["40", "90"]):
        em = Emit()
        ingest_interactions(em, [_inter(**{"Original Affinity Median nm": v})
                                 for v in order], KEYS, PROTS, IDX)
        assert em.edges[0]["qualifiers"]["ki_nm"] == 40.0


def test_an_affinity_on_only_one_row_still_reaches_the_edge():
    em = Emit()
    ingest_interactions(em, [_inter(**{"Original Affinity Median nm": ""}),
                             _inter(**{"Original Affinity Median nm": "40"})],
                        KEYS, PROTS, IDX)
    assert em.edges[0]["qualifiers"]["ki_nm"] == 40.0


def test_an_edge_with_no_affinity_carries_no_affinity_qualifier():
    em = Emit()
    ingest_interactions(em, [_inter(**{"Original Affinity Median nm": ""})],
                        KEYS, PROTS, IDX)
    assert "ki_nm" not in em.edges[0]["qualifiers"]


def test_a_changed_direction_is_counted_so_the_driver_can_report_it():
    em = Emit()
    s = ingest_interactions(em, _two("Agonist", "Inhibitor"), KEYS, PROTS, IDX)
    assert s["duplicate_changed_direction"] == 1
    em = Emit()
    s = ingest_interactions(em, _two("Agonist", "Agonist"), KEYS, PROTS, IDX)
    assert not s["duplicate_changed_direction"]


def test_the_direction_counters_count_edges_and_not_rows():
    """They are emitted once per pair now. Counting per row would double
    every figure the driver prints."""
    em = Emit()
    s = ingest_interactions(em, _two("Inhibitor", "Inhibitor"),
                            KEYS, PROTS, IDX)
    assert s["edges"] == 1 and s["direction_inhibitor"] == 1


def test_distinct_pairs_are_still_separate_edges():
    """The accumulator is keyed on the pair; a second protein must not be
    folded into the first."""
    em = Emit()
    s = ingest_interactions(
        em, [_inter(), _inter(**{"Target UniProt ID": "Q00987"})],
        KEYS, PROTS | {"UniProtKB:Q00987"}, IDX)
    assert s["edges"] == 2 and not s["duplicate"]


def test_the_reconciliation_is_the_assemblers_own_function():
    """One implementation. A second copy of this rule in a layer is how the
    graph would come to answer the same question two ways."""
    import inspect

    import kgav.gtop as G
    from kgav.assemble import resolve_qualifier as canonical
    assert G.resolve_qualifier is canonical
    assert "resolve_qualifier(" in inspect.getsource(G.ingest_interactions)


# ------------------------------------------------- the edges must be legal
def test_the_emitted_edges_validate_against_the_schema():
    """THE TEST THAT WAS NOT HERE, and its absence cost a rebuild.

    102 tests in this file built edges with Emit() and never ran one
    through the validator, so `affinity_nm` -- a qualifier gtop.py invented
    and the schema has never declared on TARGETS -- could only surface in
    `make rebuild`: 4,270 violations, 8 seconds into a phase-4 assembly,
    9 minutes in. tests/test_host_targets.py has validated its own output
    all along; this layer simply never borrowed the idea.

    TARGETS already declares kd_nm, ic50_nm, ki_nm and ec50_nm. Which
    measurement a number is forms part of the claim, which is why there is
    no untyped `affinity` field to put it in.
    """
    from kgav.schema import load_schema
    em = Emit()
    ingest_interactions(em, [_inter()], KEYS, PROTS, IDX)
    assert em.edges, "nothing was emitted, so nothing was validated"
    # COMPLETE stub nodes. Incomplete ones raise MISSING_PROPERTY, which
    # buries the edge-level violation this test is for -- the first run of
    # it reported nine node defects and not the qualifier.
    nodes = [{"id": NODE, "class": "SmallMolecule",
              "properties": {"smiles": "CCO", "inchikey_skel": SKEL,
                             "is_approved": False, "salt_collapsed": False,
                             "stereo_collapsed": False, "tautomer_collapsed": True}},
             {"id": "UniProtKB:P04637", "class": "Protein",
              "properties": {"taxon_id": "NCBITaxon:9606", "is_viral": False,
                             "sequence_hash": "ab", "reviewed": True}}]
    violations = load_schema().validate_batch(
        nodes + list(em.nodes.values()), em.edges)
    assert violations == [], [str(v) for v in violations[:6]]


def test_two_measurement_types_are_kept_apart():
    """A Ki and an IC50 are not comparable, so the minimum is taken PER
    TYPE. Taking the smaller of the two would invent a measurement no
    paper reported."""
    em = Emit()
    ingest_interactions(em, [
        _inter(**{"Original Affinity Median nm": "40",
                  "Original Affinity Units": "pKi"}),
        _inter(**{"Original Affinity Median nm": "90",
                  "Original Affinity Units": "pIC50"})],
        KEYS, PROTS, IDX)
    q = em.edges[0]["qualifiers"]
    assert q["ki_nm"] == 40.0 and q["ic50_nm"] == 90.0


@pytest.mark.parametrize("unit", ["pA2", "pKB", "-", "", "Hill coefficient"])
def test_an_unrecognised_unit_emits_no_potency_rather_than_a_guess(unit):
    """FAIL TO SILENCE, NOT TO A GUESS. AFFINITY_UNIT_TO_QUALIFIER was
    written from GtoP's documented vocabulary and not a measured one --
    which is the mistake that produced this module's entire history. An
    unmapped unit drops its number and says so, rather than landing it
    under a qualifier it does not mean."""
    em = Emit()
    s = ingest_interactions(em, [_inter(**{"Original Affinity Units": unit})],
                            KEYS, PROTS, IDX)
    assert s["edges"] == 1 and s["affinity_unit_unmapped"] == 1
    assert em.edges[0]["qualifiers"] == {"direction": "inhibitor"}


def test_every_mapped_unit_names_a_qualifier_the_schema_declares():
    """The table is only safe because the schema agrees with it. A
    qualifier here that TARGETS does not declare is precisely the defect
    this replaced, and it would again surface only in a rebuild."""
    import yaml
    spec = yaml.safe_load((ROOT / "schema" / "kg_schema.yaml").read_text())
    targets = next(e for e in spec["edge_classes"]
                   if e["predicate"] == "TARGETS")
    declared = set(targets["qualifiers"])
    assert set(AFFINITY_UNIT_TO_QUALIFIER.values()) <= declared, (
        set(AFFINITY_UNIT_TO_QUALIFIER.values()) - declared)


def test_a_typed_potency_still_validates_against_the_schema():
    """The validation test above runs with _inter's default pKi, so this
    pins that the TYPED path is the one being validated -- not a row that
    happens to emit no potency at all."""
    from kgav.schema import load_schema
    em = Emit()
    ingest_interactions(em, [_inter()], KEYS, PROTS, IDX)
    assert "ki_nm" in em.edges[0]["qualifiers"], em.edges[0]["qualifiers"]
    nodes = [{"id": NODE, "class": "SmallMolecule",
              "properties": {"smiles": "CCO", "inchikey_skel": SKEL,
                             "is_approved": False, "salt_collapsed": False,
                             "stereo_collapsed": False, "tautomer_collapsed": True}},
             {"id": "UniProtKB:P04637", "class": "Protein",
              "properties": {"taxon_id": "NCBITaxon:9606", "is_viral": False,
                             "sequence_hash": "ab", "reviewed": True}}]
    assert load_schema().validate_batch(nodes, em.edges) == []


# ------------------------------------------- the unmapped affinity unit
# Fix-list E7. The edge STILL SHIPS when the affinity unit is not in
# AFFINITY_UNIT_TO_QUALIFIER -- with a direction and no number -- so it is
# indistinguishable from an edge where GtoPdb reported no affinity at all,
# and those are different facts. On release 2026.3.1, 50 of 24,968 rows
# carry a value whose unit is unmapped: 43 KB (the antagonist dissociation
# constant from a functional assay, a real measure this schema has no
# qualifier for), 6 "-", 1 blank.
def test_an_unmapped_affinity_unit_is_named_not_merely_counted():
    em = Emit()
    row = _inter(**{"Original Affinity Median nm": "42",
                    "Original Affinity Units": "KB"})
    stats = ingest_interactions(em, [row], KEYS, PROTS, IDX)
    assert stats["affinity_unit_unmapped"] == 1
    assert stats["affinity_unit:KB"] == 1
    # The edge ships, and carries no potency -- which is why the counter
    # has to exist for anyone to know.
    assert stats["edges"] == 1
    q = em.edges[0]["qualifiers"]
    assert not any(k.endswith("_nm") for k in q)


def test_no_affinity_reported_is_distinguished_from_unreadable():
    em = Emit()
    row = _inter(**{"Original Affinity Median nm": "",
                    "Original Affinity Units": ""})
    stats = ingest_interactions(em, [row], KEYS, PROTS, IDX)
    assert stats["no_affinity_reported"] == 1
    assert stats["affinity_unit_unmapped"] == 0


# ------------------------- P3: why an accession did not join, not just that
#
# One counter stood where the protein lookup failed -- `protein_not_in_graph`
# -- and it could not tell "this protein is not in our proteome" from
# "UniProt renamed it and we never noticed". normalize_bio.UniProtResolver
# returns exactly that distinction and was reached by nothing but its own
# test, while data/raw/refs/sec_ac.txt (1,159,114 accessions) and
# delac_sp.txt (1,912) were downloaded on every build and read by no code.
#
# THE AUDIT'S OWN REMEDY WOULD HAVE RECOVERED NOTHING. It said "use the
# alias index in chembl, host_targets, gtop"; four of the five ingests
# already use AliasIndex, and all 19,313 human single-accession rows in
# GtoPdb 2026.3.1 join directly -- 0 of the graph's 20,888 UniProtKB nodes
# carry a secondary or deleted accession. The behaviour below is therefore
# unchanged on this release, verified byte-identical, and what is added is
# the ability to SAY SO next time.


def _resolver(secondary=None, deleted=()):
    """A UniProtResolver with its lists injected rather than downloaded."""
    from kgav.normalize_bio import UniProtResolver
    r = UniProtResolver()
    r.secondary = {k: list(v) for k, v in (secondary or {}).items()}
    r.deleted = set(deleted)
    return r


def test_a_secondary_accession_joins_through_its_primary():
    """The row GtoP cites by an old name. The protein IS in the graph; it
    was being dropped for a reason that is not about the protein."""
    em = Emit()
    s = ingest_interactions(
        em, [_inter(**{"Target UniProt ID": "Q00001"})], KEYS, PROTS, IDX,
        accessions=_resolver(secondary={"Q00001": ["P04637"]}))
    assert s["edges"] == 1, dict(s)
    assert s["recovered_via_secondary_accession"] == 1
    assert s["not_in_graph_accession_secondary"] == 1
    assert s["protein_not_in_graph"] == 0
    assert em.edges[0]["object"] == "UniProtKB:P04637", (
        "the edge has to point at the canonical node, not the old id")


def test_a_deleted_accession_is_named_as_deleted_and_dropped():
    em = Emit()
    s = ingest_interactions(
        em, [_inter(**{"Target UniProt ID": "Q00002"})], KEYS, PROTS, IDX,
        accessions=_resolver(deleted={"Q00002"}))
    assert s["edges"] == 0
    assert s["not_in_graph_accession_deleted"] == 1
    assert s["protein_not_in_graph"] == 1, (
        "a withdrawn entry is still an unjoined row and must stay in the "
        "accounting")
    assert s["recovered_via_secondary_accession"] == 0


def test_a_demerged_accession_is_not_resolved_silently():
    """UniProtResolver returns no curie for an accession that maps to
    several primaries, because picking one depends on which isoform the
    source meant. The row is dropped and the reason is named."""
    em = Emit()
    s = ingest_interactions(
        em, [_inter(**{"Target UniProt ID": "Q00003"})], KEYS, PROTS, IDX,
        accessions=_resolver(secondary={"Q00003": ["P04637", "P04638"]}))
    assert s["edges"] == 0
    assert s["not_in_graph_accession_demerged"] == 1
    assert s["recovered_via_secondary_accession"] == 0


def test_an_absent_protein_is_still_just_absent():
    """Not every miss is an alias. A primary accession the proteome does
    not hold must read as a missing protein, which is what the single
    counter always meant and now means only this."""
    em = Emit()
    s = ingest_interactions(
        em, [_inter(**{"Target UniProt ID": "P99999"})], KEYS, PROTS, IDX,
        accessions=_resolver())
    assert s["edges"] == 0
    assert s["protein_not_in_graph"] == 1
    assert s["not_in_graph_accession_primary"] == 1
    assert s["recovered_via_secondary_accession"] == 0


def test_without_the_lists_the_behaviour_is_exactly_the_old_one():
    """The lists are a download. An absent file must degrade to the old
    counter rather than change what joins -- which is also what makes the
    layer byte-identical on the release of 2026-11-09."""
    em = Emit()
    s = ingest_interactions(
        em, [_inter(**{"Target UniProt ID": "Q00001"})], KEYS, PROTS, IDX)
    assert s["edges"] == 0
    assert s["protein_not_in_graph"] == 1
    assert not [k for k in s if k.startswith("not_in_graph_accession_")], (
        "a classification was reported with no accession lists loaded")


def test_the_resolver_is_no_longer_reached_only_by_its_own_test():
    """The defect in one line. `UniProtResolver` is imported by gtop now;
    before 2026-11-09 the only module that named it was tests/."""
    import kgav.gtop as gtop_module
    assert gtop_module.UniProtResolver is not None


# --------------------------------------- the spread behind the point (P14)
def test_a_single_row_states_a_range_of_one():
    """EMITTED ON EVERY EDGE WITH A POTENCY, including the one-row case.
    Emitting the four only for replicates would make their absence mean two
    things -- one measurement, or a layer that predates this -- which is the
    ambiguity host_targets' own comment says this project keeps paying for.
    """
    em = Emit()
    ingest_interactions(em, [_inter()], KEYS, PROTS, IDX)
    q = em.edges[0]["qualifiers"]
    assert q["potency_spread_field"] == "ki_nm"
    assert (q["potency_n"], q["potency_min_nm"], q["potency_max_nm"]) == (
        1, 40.0, 40.0)


def test_duplicate_rows_become_the_range_and_leave_the_point_alone():
    """The point is the most potent value by the rule stated above this
    layer's accumulator; P14 adds what that rule was choosing BETWEEN. The
    point does not move -- that would be a second change hiding inside this
    one -- it is simply no longer the only thing on the edge."""
    from collections import Counter
    em = Emit()
    s = ingest_interactions(em, [_inter(**{"Original Affinity Median nm": v,
                                           "PubMed ID": p})
                                 for v, p in (("40", "1"), ("4000", "2"),
                                              ("400", "3"))],
                            KEYS, PROTS, IDX)
    assert s["edges"] == 1
    q = em.edges[0]["qualifiers"]
    assert q["ki_nm"] == 40.0, "the point moved"
    assert (q["potency_n"], q["potency_min_nm"], q["potency_max_nm"]) == (
        3, 40.0, 4000.0)
    assert s["edges_collapsing_several_rows"] == 1
    assert s["edges_whose_rows_disagree"] == 1
    assert isinstance(s, Counter)


def test_rows_that_agree_are_counted_apart_from_rows_that_do_not():
    """Three repeats of one number are not a disagreement, and a counter
    that cannot tell them apart says nothing about the data."""
    em = Emit()
    s = ingest_interactions(em, [_inter(**{"PubMed ID": p}) for p in "123"],
                            KEYS, PROTS, IDX)
    assert s["edges_collapsing_several_rows"] == 1
    assert not s["edges_whose_rows_disagree"]
    q = em.edges[0]["qualifiers"]
    assert q["potency_n"] == 3 and q["potency_min_nm"] == q["potency_max_nm"]


def test_the_named_field_is_the_one_with_rows_behind_it():
    """ONE RANGE SLOT, TWO MEASUREMENTS. A field with a single row has no
    spread to report, so naming it would spend the edge's only range slot
    on the measurement that needed it least."""
    em = Emit()
    rows = [_inter(**{"Original Affinity Units": "IC50",
                      "Original Affinity Median nm": v, "PubMed ID": p})
            for v, p in (("10", "1"), ("1000", "2"))]
    rows.append(_inter(**{"Original Affinity Units": "pKi",
                          "Original Affinity Median nm": "40",
                          "PubMed ID": "9"}))
    em2 = Emit()
    s = ingest_interactions(em, rows, KEYS, PROTS, IDX)
    # and the same answer with the rows the other way round
    s2 = ingest_interactions(em2, list(reversed(rows)), KEYS, PROTS, IDX)
    for emit, stats in ((em, s), (em2, s2)):
        q = emit.edges[0]["qualifiers"]
        assert q["ic50_nm"] == 10.0 and q["ki_nm"] == 40.0, (
            "both measurements belong on the edge; only the RANGE is one")
        assert q["potency_spread_field"] == "ic50_nm"
        assert (q["potency_min_nm"], q["potency_max_nm"]) == (10.0, 1000.0)
        assert stats["edge_measures_several_fields_one_spread_kept"] == 1
        # the ki has one row, so no spread was lost with it
        assert not stats["spread_not_recorded_for_another_measured_field"]


def test_a_spread_that_could_not_be_recorded_is_counted():
    em = Emit()
    rows = [_inter(**{"Original Affinity Units": u,
                      "Original Affinity Median nm": v, "PubMed ID": p})
            for u, v, p in (("IC50", "10", "1"), ("IC50", "30", "2"),
                            ("IC50", "90", "3"),
                            ("pKi", "40", "4"), ("pKi", "4000", "5"))]
    s = ingest_interactions(em, rows, KEYS, PROTS, IDX)
    q = em.edges[0]["qualifiers"]
    assert q["potency_spread_field"] == "ic50_nm"
    assert s["spread_not_recorded_for_another_measured_field"] == 1, (
        "the Ki spans 40-4000 nM and the edge cannot say so; an edge that "
        "loses a measured spread has to be countable")


def test_a_tie_is_broken_by_the_declared_order_not_by_dict_order():
    """Two fields with one row each cannot be separated by row count, and
    `max` over a dict would then return whichever key was inserted first --
    which is file order, the defect this layer's accumulator exists for."""
    from kgav.assemble import POTENCY_FIELD_ORDER
    named = []
    for order in (("pKd", "IC50"), ("IC50", "pKd")):
        em = Emit()
        ingest_interactions(
            em, [_inter(**{"Original Affinity Units": u,
                           "Original Affinity Median nm": "7",
                           "PubMed ID": str(i)})
                 for i, u in enumerate(order)], KEYS, PROTS, IDX)
        named.append(em.edges[0]["qualifiers"]["potency_spread_field"])
    assert named == ["ic50_nm", "ic50_nm"], (
        f"the named field followed row order: {named}")
    assert POTENCY_FIELD_ORDER.index("ic50_nm") < POTENCY_FIELD_ORDER.index(
        "kd_nm"), "the tie-break above no longer matches the declared order"


def test_every_unit_this_layer_maps_is_a_field_the_range_can_name():
    """_attach_spread indexes POTENCY_FIELD_ORDER to break ties, so a unit
    mapped to a field outside it raises ValueError nine minutes into a
    rebuild -- the same shape of failure the unit table's own comment
    records from `affinity_nm`."""
    from kgav.assemble import POTENCY_FIELD_ORDER
    assert set(AFFINITY_UNIT_TO_QUALIFIER.values()) <= set(POTENCY_FIELD_ORDER)


def test_an_affinity_stated_only_as_a_range_is_named_not_folded_away():
    """940 rows on 2026.3.1 state a low and/or a high with no median. The
    layer carries no number for them -- a point out of a bare interval
    would be invented -- but `no_affinity_reported` read as if GtoPdb had
    said nothing about them."""
    em = Emit()
    s = ingest_interactions(em, [_inter(**{
        "Original Affinity Median nm": "",
        "Original Affinity Low nm": "10",
        "Original Affinity High nm": "100"})], KEYS, PROTS, IDX)
    assert s["edges"] == 1
    assert s["no_affinity_reported"] == 1
    assert s["affinity_range_without_median"] == 1
    q = em.edges[0]["qualifiers"]
    assert "ki_nm" not in q, "a point was invented from an interval"
    assert "potency_min_nm" not in q, (
        "a range with no point is what assemble drops: potency_spread_field "
        "has to name a qualifier that is on the edge")
    # a row with neither is still only the first counter
    em2 = Emit()
    s2 = ingest_interactions(em2, [_inter(**{
        "Original Affinity Median nm": ""})], KEYS, PROTS, IDX)
    assert s2["no_affinity_reported"] == 1
    assert not s2["affinity_range_without_median"]


def test_the_export_never_states_a_median_and_a_range_on_one_row():
    """THE MEASUREMENT _attach_spread RESTS ON, re-derived from the
    download rather than quoted from it. If GtoPdb ever ships both on one
    row, the range should come from the source's own interval instead of
    from counting duplicates -- and `potency_n` would stop meaning the same
    thing it means in host_targets, which is the reason this is a test and
    not a comment."""
    import re
    from kgav.gtop import _attach_spread
    path = ROOT / "data" / "raw" / "gtop" / "interactions.csv"
    if not path.exists():
        pytest.skip("the GtoP download is not in this checkout")
    both = medians = ranges_only = 0
    for r in _rows(path, INTERACTION_COLUMNS):
        med = (r.get("Original Affinity Median nm") or "").strip()
        span = ((r.get("Original Affinity Low nm") or "").strip()
                or (r.get("Original Affinity High nm") or "").strip())
        unit = (r.get("Original Affinity Units") or "").strip().upper()
        if med and unit in AFFINITY_UNIT_TO_QUALIFIER:
            medians += 1
            both += bool(span)
        elif not med and span:
            ranges_only += 1
    assert both == 0, (
        f"{both} rows state a median AND an interval; _attach_spread says "
        f"none do, and builds the range by counting duplicates instead")
    doc = " ".join(_attach_spread.__doc__.split())
    m = re.search(r"of the ([\d,]+) rows that state a usable median, "
                  r"\*\*(\d+)\*\* state a low or a high\. The ([\d,]+) rows",
                  doc)
    assert m, "the docstring no longer states the census in the form read here"
    if file_version(path) == "2026.3.1":
        got = (medians, both, ranges_only)
        want = tuple(int(g.replace(",", "")) for g in m.groups())
        assert got == want, f"docstring says {want}, the download gives {got}"
