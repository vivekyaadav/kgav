"""The declared graph against the built one.

The 2026-10-24 audit found four declared node classes with no nodes and
ten predicates with no edges, including BOTH the schema marks
evaluation_label. None of it was visible in assemble_release's output,
which lists what IS there: a class with no nodes prints no line, and a
reader counting lines sees a complete graph.
"""
from __future__ import annotations

import io
import json
import sys
from contextlib import redirect_stdout
from importlib import import_module
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

da = import_module("data_audit")
SCHEMA = yaml.safe_load((ROOT / "schema" / "kg_schema.yaml").read_text())


def _rel(tmp: Path, nodes: list[dict], edges: list[dict]) -> Path:
    d = tmp / "rel"
    d.mkdir(parents=True, exist_ok=True)
    (d / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in nodes))
    (d / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in edges))
    return d


def _run(rel: Path, *extra) -> str:
    buf = io.StringIO()
    with redirect_stdout(buf):
        da.main(["--release", str(rel), *extra])
    return buf.getvalue()


MOL = {"id": "INCHIKEY:C0", "class": "SmallMolecule"}
PRO = {"id": "UniProtKB:P0", "class": "Protein"}
TAX = {"id": "NCBITaxon:2697049", "class": "OrganismTaxon"}


# ------------------------------------------------ declared vs populated
def test_a_declared_class_with_no_nodes_is_reported(tmp_path):
    """assemble_release prints one line per class it HAS. A class with no
    nodes prints nothing, so the absence is invisible."""
    out = _run(_rel(tmp_path, [MOL], []))
    assert "VirusStrain" in out and "EMPTY" in out


def test_a_declared_predicate_with_no_edges_is_reported(tmp_path):
    out = _run(_rel(tmp_path, [MOL], []))
    assert "CONFERS_RESISTANCE_TO" in out


def test_nothing_is_hardcoded_from_the_schema(tmp_path):
    """A predicate added to the schema and never ingested must appear on
    the next run, not when someone notices."""
    out = _run(_rel(tmp_path, [MOL], []))
    for p in (e["predicate"] for e in SCHEMA["edge_classes"]):
        assert p in out, p
    for c in SCHEMA["node_classes"]:
        assert c in out, c


def test_something_in_the_release_but_not_the_schema_is_flagged(tmp_path):
    """The audit must look both ways, or an ingest writing an undeclared
    predicate passes unmentioned."""
    rel = _rel(tmp_path, [dict(MOL, **{"class": "Sandwich"})],
               [{"subject": "a", "predicate": "ENJOYS", "object": "b"}])
    out = _run(rel)
    assert "NOT IN THE SCHEMA" in out
    assert "Sandwich" in out and "ENJOYS" in out


# ------------------------------------------- an empty class, two reasons
def test_a_class_stored_as_a_qualifier_says_so(tmp_path):
    """CellLine has no nodes and needs none: cell lines live in a
    `cell_line` qualifier. Reporting that as a data gap would send someone
    looking for a source that is not wanted."""
    rel = _rel(tmp_path, [PRO, TAX],
               [{"subject": "UniProtKB:P0", "predicate": "HOST_FACTOR_FOR",
                 "object": "NCBITaxon:2697049",
                 "qualifiers": {"direction": "dependency",
                                "cell_line": "CVCL_0063"}}])
    out = _run(rel)
    assert "stored as a `cell_line` qualifier" in out


def test_a_predicate_empty_because_its_class_is_empty_says_which(tmp_path):
    """HAS_VARIANT cannot have edges while Variant has no nodes. Listing
    it beside an unrelated gap would imply two problems."""
    out = _run(_rel(tmp_path, [MOL], []))
    line = next(x for x in out.splitlines() if x.strip().startswith("HAS_VARIANT"))
    assert "Variant is empty" in line


# ------------------------------------------------- the evaluation labels
def test_the_evaluation_labels_are_called_out_separately(tmp_path):
    """THE FINDING THAT PROMPTED THIS. Both predicates the schema marks
    evaluation_label were empty, so the project's only label was
    HAS_ANTIVIRAL_ACTIVITY_AGAINST and nothing said so."""
    out = _run(_rel(tmp_path, [MOL], []))
    assert "EVALUATION LABELS" in out
    labels = [e["predicate"] for e in SCHEMA["edge_classes"]
              if e.get("evaluation_label")]
    assert len(labels) >= 2
    section = out[out.index("EVALUATION LABELS"):out.index("METAPATHS")]
    for p in labels:
        assert p in section


def test_one_live_label_is_distinguished_from_none(tmp_path):
    none = _run(_rel(tmp_path, [MOL], []))
    assert "NO label carries data" in none
    one = _run(_rel(tmp_path, [MOL, TAX],
                    [{"subject": "INCHIKEY:C0",
                      "predicate": "HAS_ANTIVIRAL_ACTIVITY_AGAINST",
                      "object": "NCBITaxon:2697049", "qualifiers": {}}]))
    assert "Only infores" not in one
    assert "Only HAS_ANTIVIRAL_ACTIVITY_AGAINST carries data" in one


# ------------------------------------------------------------- channels
def test_a_metapath_needing_an_empty_predicate_is_dead(tmp_path):
    out = _run(_rel(tmp_path, [MOL], []))
    assert "DEAD" in out
    assert "needs" in out


def test_a_metapath_whose_predicates_all_have_edges_walks(tmp_path):
    rel = _rel(tmp_path, [MOL, PRO, TAX],
               [{"subject": "INCHIKEY:C0", "predicate": "TARGETS",
                 "object": "UniProtKB:P0",
                 "qualifiers": {"direction": "inhibitor"}},
                {"subject": "UniProtKB:P0", "predicate": "HOST_FACTOR_FOR",
                 "object": "NCBITaxon:2697049",
                 "qualifiers": {"direction": "dependency"}}])
    out = _run(rel)
    assert "M2    walks" in out


# ------------------------------------- a required enum that says nothing
def test_a_required_enum_that_is_mostly_unknown_is_flagged(tmp_path):
    """12,155 of 14,381 TARGETS edges state direction `unknown`. That
    passes validation -- the qualifier is present -- and carries no
    information, so no schema check can see it."""
    edges = [{"subject": f"INCHIKEY:C{i}", "predicate": "TARGETS",
              "object": "UniProtKB:P0",
              "qualifiers": {"direction": "unknown"}} for i in range(9)]
    edges.append({"subject": "INCHIKEY:CX", "predicate": "TARGETS",
                  "object": "UniProtKB:P0",
                  "qualifiers": {"direction": "inhibitor"}})
    out = _run(_rel(tmp_path, [MOL, PRO], edges))
    assert "TARGETS.direction" in out
    assert "90.0% unknown" in out
    assert "(required)" in out


def test_a_well_populated_enum_is_not_flagged(tmp_path):
    edges = [{"subject": f"INCHIKEY:C{i}", "predicate": "TARGETS",
              "object": "UniProtKB:P0",
              "qualifiers": {"direction": "inhibitor"}} for i in range(9)]
    out = _run(_rel(tmp_path, [MOL, PRO], edges))
    assert "TARGETS.direction" not in out.split("ENUM QUALIFIERS")[1].split(
        "PER-TAXON")[0]


# ------------------------------------------------------ taxon coverage
def test_a_taxon_with_no_edges_appears_with_zeros(tmp_path):
    """THE MOST USEFUL ROW IN THE TABLE, and one built from edges alone
    cannot contain it. JEV and TBEV are built and have no host factor,
    no PPI and no CRISPR screen; that must be visible."""
    rel = _rel(tmp_path,
               [TAX, {"id": "NCBITaxon:11072", "class": "OrganismTaxon"}],
               [{"subject": "NCBITaxon:2697049", "predicate": "CAUSES",
                 "object": "MONDO:1", "qualifiers": {}}])
    out = _run(rel)
    section = out[out.index("PER-TAXON"):]
    assert "NCBITaxon:11072" in section, "a taxon with no edge vanished"
    assert "NCBITaxon:2697049" in section


def test_the_gap_count_is_not_a_judgement_about_size(tmp_path):
    out = _run(_rel(tmp_path, [MOL], []))
    assert "not a judgement that a" in out


def test_a_missing_release_is_an_error_not_an_empty_report(tmp_path):
    assert da.main(["--release", str(tmp_path / "none")]) == 2


# ------------------------------------------------------------- the notes
def test_every_note_names_a_real_class_or_predicate():
    """A note keyed on something the schema no longer has would never
    print, and would read as recorded knowledge that is silently lost."""
    known = set(SCHEMA["node_classes"]) | {
        e["predicate"] for e in SCHEMA["edge_classes"]}
    assert not set(da.NOTES) - known, set(da.NOTES) - known


def test_notes_are_dated_because_they_are_judgements_not_derivations():
    for key, (when, text) in da.NOTES.items():
        assert len(when) == 10 and when[4] == "-", (key, when)
        assert text.strip(), key


@pytest.mark.parametrize("key", ["TREATS", "IN_TRIAL_FOR", "CAUSES"])
def test_the_clinical_remedies_are_recorded(key):
    """Worked out once on 2026-10-24: ChEMBL's drug_indication table is
    already on disk and fills both evaluation labels, blocked behind the
    Disease axis. A remedy worked out and lost is a remedy worked out
    twice."""
    assert key in da.NOTES
    assert "drug_indication" in da.NOTES[key][1] or "mondo" in da.NOTES[key][1]


def test_a_note_does_not_contradict_the_line_above_it(tmp_path):
    """VirusStrain printed "stored as a `virus_strain` qualifier" and then
    "Genuinely absent" directly beneath. Both were true for different
    reasons -- the strain is named in free text and cannot be joined on --
    but as printed they disagreed, which is the label-against-evidence
    defect this repository keeps finding, here in the tool built to find
    it."""
    rel = _rel(tmp_path, [MOL, TAX],
               [{"subject": "INCHIKEY:C0",
                 "predicate": "HAS_ANTIVIRAL_ACTIVITY_AGAINST",
                 "object": "NCBITaxon:2697049",
                 "qualifiers": {"virus_strain": "USA-WA1/2020",
                                "assay_type": "cell_based_antiviral"}}])
    out = _run(rel)
    block = out[out.index("VirusStrain"):out.index("Variant")]
    assert "stored as a `virus_strain` qualifier" in block
    assert "Genuinely absent" not in block
    assert "cannot be joined on" in block


def test_a_truly_vestigial_class_still_says_so(tmp_path):
    """The fix must not have removed the distinction. CellLine IS
    vestigial; VirusStrain is not."""
    rel = _rel(tmp_path, [PRO, TAX],
               [{"subject": "UniProtKB:P0", "predicate": "HOST_FACTOR_FOR",
                 "object": "NCBITaxon:2697049",
                 "qualifiers": {"direction": "dependency",
                                "cell_line": "CVCL_0063"}}])
    out = _run(rel)
    block = out[out.index("CellLine"):out.index("Publication")]
    assert "none expected" in block
