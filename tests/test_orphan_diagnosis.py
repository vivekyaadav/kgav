"""Telling the three orphan causes apart, on built fixtures.

v0.1 ships one viral orphan -- a SARS-CoV-2 protein in the graph that no
edge touches. "Protein 634" hid it for every release until the host/viral
split; naming it and saying WHICH failure it is needs the per-layer files,
because the whole question is whether something was present before assembly
and absent after.
"""
from __future__ import annotations

import json
import sys
from importlib import import_module
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
V = "NCBITaxon:2697049"


def _diag():
    sys.path.insert(0, str(ROOT / "scripts"))
    return import_module("diagnose_orphans")


def _prot(nid, **props):
    return {"id": nid, "class": "Protein",
            "properties": {"is_viral": True, "taxon_id": V, **props}}


# ----------------------------------------------------------- untouched()
def test_a_node_in_object_position_is_not_an_orphan():
    """The easy wrong version checks subjects only. Every ENCODED_BY points
    AT a Gene, and every INHIBITS points at a Protein, so a subjects-only
    sweep would declare most of the graph orphaned."""
    d = _diag()
    nodes = [_prot("P:sub"), {"id": "G:obj", "class": "Gene", "properties": {}}]
    edges = [{"subject": "P:sub", "predicate": "ENCODED_BY", "object": "G:obj"}]
    assert d.untouched(nodes, edges) == []


def test_a_node_named_by_no_edge_is_an_orphan():
    d = _diag()
    nodes = [_prot("P:lone"), _prot("P:joined")]
    edges = [{"subject": "P:joined", "predicate": "X", "object": "P:joined"}]
    assert [n["id"] for n in d.untouched(nodes, edges)] == ["P:lone"]


# ------------------------------------------------------------- classify()
def test_a_layer_that_emitted_an_edge_means_assembly_dropped_it():
    """The layer file is the evidence: the edge existed and the assembled
    graph has none, so precedence, merging or an edge-key collision ate it."""
    from collections import Counter
    d = _diag()
    assert d.classify(_prot("P:x"), Counter({"vhppi": 2}), []) == \
        "DROPPED AT ASSEMBLY"


def test_a_connected_sibling_means_the_identifier_was_discarded():
    """This repo's README calls a discarded identifier "the dominant bug
    class" and counts four. The evidence is not missing -- it is filed
    under a second id for the same protein."""
    from collections import Counter
    d = _diag()
    sibs = [("sequence_hash", "P:twin", True)]
    assert d.classify(_prot("P:x"), Counter(), sibs) == "DISCARDED IDENTIFIER"


def test_a_sibling_that_is_also_orphaned_is_not_a_discarded_identifier():
    """Two orphaned copies are two orphans, not evidence filed elsewhere.
    Reporting them as a discarded identifier would send someone looking for
    edges that do not exist."""
    from collections import Counter
    d = _diag()
    sibs = [("sequence_hash", "P:twin", False)]
    assert d.classify(_prot("P:x"), Counter(), sibs) == "NEVER REFERENCED"


def test_a_dropped_edge_outranks_a_sibling():
    """If both hold, the layer edge is the concrete actionable one: a named
    layer emitted a named edge for THIS id and the graph lost it."""
    from collections import Counter
    d = _diag()
    sibs = [("sequence_hash", "P:twin", True)]
    assert d.classify(_prot("P:x"), Counter({"chembl": 1}), sibs) == \
        "DROPPED AT ASSEMBLY"


# ------------------------------------------------------------- siblings()
def test_a_sibling_must_share_a_class_as_well_as_a_hash():
    """A Gene and a Protein can carry the same label without being the same
    thing, and calling that a discarded identifier would be a false lead."""
    d = _diag()
    me = _prot("P:x", sequence_hash="h")
    other = {"id": "G:x", "class": "Gene",
             "properties": {"sequence_hash": "h"}}
    assert d.siblings(me, [me, other], {"G:x"}) == []


def test_an_empty_identity_property_matches_nobody():
    """Most nodes have no sequence_hash at all. Matching empty against empty
    would make every one of them everyone else's twin."""
    d = _diag()
    me = _prot("P:x", sequence_hash="", mature_peptide="")
    others = [_prot(f"P:{i}", sequence_hash="", mature_peptide="")
              for i in range(3)]
    assert d.siblings(me, [me, *others], set()) == []


def test_a_sibling_is_found_on_either_identity_property():
    d = _diag()
    me = _prot("P:x", mature_peptide="nsp16")
    twin = _prot("P:y", mature_peptide="nsp16")
    found = d.siblings(me, [me, twin], {"P:y"})
    assert found == [("mature_peptide", "P:y", True)]


# ------------------------------------------------- end to end on files
def test_the_three_causes_are_told_apart_in_one_run(tmp_path, capsys,
                                                   monkeypatch):
    """Built from files, not from mocks, because the point of the script is
    to read a release and the per-layer directories beside it."""
    d = _diag()
    rel = tmp_path / "v0.1"
    rel.mkdir()
    (tmp_path / "v0.1-vhppi").mkdir()

    nodes = [
        _prot("P:drop", gene_symbol="orf9b", sequence_hash="h1"),
        _prot("P:twinA", sequence_hash="h2"),
        _prot("P:twinB", sequence_hash="h2"),
        _prot("P:lone", gene_symbol="", sequence_hash="h3"),
        _prot("P:ok", gene_symbol="S", sequence_hash="h4"),
        {"id": "H:1", "class": "Protein", "properties": {"is_viral": False}},
    ]
    edges = [
        {"subject": "P:twinB", "predicate": "PHYSICALLY_INTERACTS_WITH",
         "object": "P:ok"},
    ]
    for p, rows in ((rel / "nodes.jsonl", nodes), (rel / "edges.jsonl", edges),
                    (tmp_path / "v0.1-vhppi" / "edges.jsonl",
                     [{"subject": "P:drop", "predicate": "X",
                       "object": "P:ok"}])):
        p.write_text("".join(json.dumps(r) + "\n" for r in rows))

    # monkeypatch, not a bare assignment: tests share a process and a
    # leaked sys.argv would reach whatever runs next in this worker.
    monkeypatch.setattr(sys, "argv",
                        ["diagnose_orphans", "--release", str(rel),
                         "--releases", str(tmp_path)])
    assert d.main() == 0
    out = capsys.readouterr().out
    assert "P:drop" in out and "DROPPED AT ASSEMBLY" in out
    assert "P:twinA" in out and "DISCARDED IDENTIFIER" in out
    assert "P:lone" in out and "NEVER REFERENCED" in out
    assert "P:ok" not in out      # connected, so never listed as an orphan
    # the host orphan is counted but not itemised without --all
    assert "1 host/other orphans not listed" in out


def test_a_release_with_no_viral_orphan_says_so_rather_than_printing_nothing(
        tmp_path, capsys, monkeypatch):
    """An empty section reads as "not checked". The benign case has to be
    stated."""
    d = _diag()
    rel = tmp_path / "v0.1"
    rel.mkdir()
    (rel / "nodes.jsonl").write_text(json.dumps(
        {"id": "H:1", "class": "Protein",
         "properties": {"is_viral": False}}) + "\n")
    (rel / "edges.jsonl").write_text("")
    monkeypatch.setattr(sys, "argv",
                        ["diagnose_orphans", "--release", str(rel),
                         "--releases", str(tmp_path)])
    assert d.main() == 0
    assert "No viral orphan" in capsys.readouterr().out


# ------------------------- the correction: a label is not identity
def test_the_same_label_in_a_different_virus_is_not_the_same_protein():
    """THE BUG THIS FILE SHIPPED WITH.

    Run against v0.1 the script called SARS-CoV-2's orphaned pp1a a
    DISCARDED IDENTIFIER on six matches -- every one a different
    coronavirus whose pp1a chain carries the same description, "Replicase
    polyprotein 1a". ingest_spine dedups chain labels in `seen_chains`,
    which is local to ingest_virus, so one label per virus is expected and
    proves nothing across viruses. Six different proteins, reported as one
    protein filed twice.

    My first test for siblings() used _prot() for both nodes, so both
    carried the same taxon and the cross-virus case was never exercised:
    a fixture correct about itself and silent about the world, which is
    the failure this project keeps finding.
    """
    d = _diag()
    me = {"id": "UniProtKB:PRO_0000449634", "class": "Protein",
          "properties": {"is_viral": True, "taxon_id": "NCBITaxon:2697049",
                         "mature_peptide": "Replicase polyprotein 1a"}}
    others = [{"id": f"UniProtKB:PRO_000033{i}", "class": "Protein",
               "properties": {"is_viral": True, "taxon_id": f"NCBITaxon:{t}",
                              "mature_peptide": "Replicase polyprotein 1a"}}
              for i, t in enumerate(("694009", "1335626", "11137", "277944",
                                     "31631", "290028"))]
    connected = {o["id"] for o in others}
    assert d.siblings(me, [me, *others], connected) == []

    from collections import Counter
    assert d.classify(me, Counter(), []) == "NEVER REFERENCED"


def test_the_same_label_in_the_SAME_virus_is_still_reported():
    """The taxon clause must not switch the check off. Within one virus the
    spine's own dedup should already have collapsed a repeated label, so a
    same-taxon match is a finding about that dedup."""
    d = _diag()
    me = _prot("P:a", mature_peptide="nsp16")
    twin = _prot("P:b", mature_peptide="nsp16")      # same taxon via _prot
    assert d.siblings(me, [me, twin], {"P:b"}) == \
        [("mature_peptide", "P:b", True)]


def test_a_chain_sequence_hash_is_built_from_the_id_not_the_sequence():
    """Anchored to ingest_spine's source, because the script's advice
    depends on it: a chain's sequence_hash is _sha(accession:chain_id), so
    it is injective and CANNOT match another chain's. Reading its absence
    as "not a duplicate" would be wrong, and the diagnostic says so only
    while this stays true.
    """
    src = (ROOT / "scripts" / "ingest_spine.py").read_text()
    assert 'sequence_hash=_sha(f"{acc}:{cid}")' in src   # chain nodes
    assert "sequence_hash=_sha(seq)" in src              # entry-level nodes


def test_the_chain_hash_caveat_is_printed_for_a_PRO_node(tmp_path, capsys,
                                                         monkeypatch):
    """A reader who sees NEVER REFERENCED needs to know the strongest
    identity check was unavailable, not that it passed."""
    d = _diag()
    rel = tmp_path / "v0.1"
    rel.mkdir()
    (rel / "nodes.jsonl").write_text(json.dumps(
        _prot("UniProtKB:PRO_0000449634",
              mature_peptide="Replicase polyprotein 1a")) + "\n")
    (rel / "edges.jsonl").write_text("")
    monkeypatch.setattr(sys, "argv",
                        ["diagnose_orphans", "--release", str(rel),
                         "--releases", str(tmp_path)])
    assert d.main() == 0
    out = capsys.readouterr().out
    assert "NEVER REFERENCED" in out
    assert "_sha(accession:chain_id)" in out
