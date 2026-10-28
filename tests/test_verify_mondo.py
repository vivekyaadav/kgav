"""Gate for the MONDO verification.

ingest_spine emits a CAUSES edge straight from the register's declared
term, so an unchecked identifier becomes an unchecked edge. Three terms
were curated by hand in 2026-09 and nothing re-read them, and two of
the three recorded a LABEL MATCH as their basis -- the one method
normalize_bio.MondoResolver bans.
"""
from __future__ import annotations

import json
import sys
from importlib import import_module
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

V = import_module("verify_mondo")
from kgav.normalize_bio import MondoResolver


def _ontology(tmp_path: Path, terms) -> Path:
    """A tiny MONDO json. terms: [(curie, label, [xrefs], deprecated)]"""
    nodes = []
    for curie, label, xrefs, dep in terms:
        meta = {"xrefs": [{"val": x} for x in xrefs]}
        if dep:
            meta["deprecated"] = True
        nodes.append({"id": "http://purl.obolibrary.org/obo/"
                            + curie.replace(":", "_"),
                      "lbl": label, "meta": meta})
    p = tmp_path / "mondo.json"
    p.write_text(json.dumps({"graphs": [{"nodes": nodes}]}))
    return p


def _res(tmp_path, terms):
    return MondoResolver(_ontology(tmp_path, terms))


LIVE = [("MONDO:0005502", "dengue disease", ["DOID:12205"], False)]


def test_a_term_backed_by_its_xref_passes(tmp_path):
    r = _res(tmp_path, LIVE)
    status, _d = V.check({"name": "DENV-1", "mondo": "MONDO:0005502",
                          "mondo_source": "DOID:12205 xref, verified"}, r)
    assert status == ""


def test_a_deprecated_term_fails(tmp_path):
    """MONDO deprecates terms between releases. A CAUSES edge to one
    points at a node the ontology has withdrawn."""
    r = _res(tmp_path, [("MONDO:0005502", "dengue disease", ["DOID:12205"],
                         True)])
    status, detail = V.check({"name": "X", "mondo": "MONDO:0005502",
                              "mondo_source": "DOID:12205 xref"}, r)
    assert status == "DEPRECATED" and "deprecated" in detail


def test_a_term_not_in_the_ontology_fails(tmp_path):
    """A typo in a CURIE produces a Disease node that joins to nothing."""
    r = _res(tmp_path, LIVE)
    status, _d = V.check({"name": "X", "mondo": "MONDO:9999999",
                          "mondo_source": "DOID:12205 xref"}, r)
    assert status == "UNKNOWN"


def test_a_label_match_is_not_an_acceptable_basis(tmp_path):
    """The basis two of the three curated terms recorded. Both were
    right; the reason given for believing them was the banned one."""
    r = _res(tmp_path, LIVE)
    status, detail = V.check(
        {"name": "X", "mondo": "MONDO:0005502",
         "mondo_source": "dengue disease; sole exact match, verified"}, r)
    assert status == "NO XREF BASIS"
    assert "no cross-reference" in detail


def test_citing_the_mondo_id_itself_is_not_a_basis(tmp_path):
    """Circular: the identifier cannot justify itself."""
    r = _res(tmp_path, LIVE)
    status, _d = V.check({"name": "X", "mondo": "MONDO:0005502",
                          "mondo_source": "MONDO:0005502, verified"}, r)
    assert status == "NO XREF BASIS"


def test_an_xref_pointing_somewhere_else_fails(tmp_path):
    """The declared term and its stated basis must agree, or one of
    them is a copy-paste from the entry above."""
    r = _res(tmp_path, LIVE + [("MONDO:0002282", "West Nile fever",
                                ["DOID:2366"], False)])
    status, _d = V.check({"name": "X", "mondo": "MONDO:0002282",
                          "mondo_source": "DOID:12205 xref"}, r)
    assert status == "XREF DISAGREES"


def test_an_xref_claimed_by_two_terms_fails(tmp_path):
    """MondoResolver calls it demerged and refuses to pick."""
    r = _res(tmp_path, [("MONDO:1", "a", ["DOID:1"], False),
                        ("MONDO:2", "b", ["DOID:1"], False)])
    status, _d = V.check({"name": "X", "mondo": "MONDO:1",
                          "mondo_source": "DOID:1 xref"}, r)
    assert status == "DEMERGED"


def test_a_missing_term_is_reported_but_does_not_fail(tmp_path, capsys):
    """An absent term is honest where a guessed one is not."""
    r = _res(tmp_path, LIVE)
    status, _d = V.check({"name": "X", "mondo": None}, r)
    assert status == "NO TERM"


def test_the_run_exits_nonzero_when_a_declared_term_fails(tmp_path, capsys):
    cfg = tmp_path / "v.yaml"
    cfg.write_text(yaml.safe_dump({"viruses": [
        {"name": "Bad", "taxon": "1", "mondo": "MONDO:9999999",
         "mondo_source": "DOID:12205 xref"}]}))
    rc = V.main(["--config", str(cfg),
                 "--mondo", str(_ontology(tmp_path, LIVE))])
    assert rc == 1
    assert "FAILED" in capsys.readouterr().out


def test_a_missing_ontology_says_how_to_get_it(tmp_path, capsys):
    cfg = tmp_path / "v.yaml"
    cfg.write_text(yaml.safe_dump({"viruses": []}))
    rc = V.main(["--config", str(cfg), "--mondo", str(tmp_path / "nope.json")])
    out = capsys.readouterr().out
    assert rc == 2 and "monarch-initiative" in out


# -------------------------------------------------- the live register
def test_every_register_term_states_an_xref_basis():
    """Checked without the ontology, which is a 103 MB download: the
    BASIS is in the register and can be audited from the repository
    alone. A term with no xref cited is one nothing can re-verify."""
    doc = yaml.safe_load((ROOT / "config/viruses.yaml").read_text())
    missing = [v["name"] for v in doc["viruses"]
               if v.get("mondo") and not V.stated_xrefs(v.get("mondo_source"))]
    assert not missing, missing


def test_the_shared_terms_say_why_they_are_shared():
    """Four coronaviruses share `common cold` and four dengue serotypes
    share `dengue disease`, because MONDO has no finer term. A reader
    seeing one disease on four viruses must find the reason in the
    register, not infer a merge."""
    doc = yaml.safe_load((ROOT / "config/viruses.yaml").read_text())
    by_term: dict[str, list[dict]] = {}
    for v in doc["viruses"]:
        if v.get("mondo"):
            by_term.setdefault(v["mondo"], []).append(v)
    shared = {t: vs for t, vs in by_term.items() if len(vs) > 1}
    assert shared, "the fixture of interest must exist"
    # AN EXACT PHRASE, not a keyword soup. The first version required
    # "no" and "term" somewhere in the note, and passed on a note
    # saying MONDO *has* a specific term -- because "not a merge"
    # contains "no" and the sentence contains "term". A commitment the
    # register makes should be checkable as a commitment.
    REASON = "MONDO has no"
    for term, vs in shared.items():
        for v in vs:
            note = v["mondo_source"]
            assert REASON in note, (
                f"{v['name']} shares {term} without stating that MONDO "
                f"offers nothing finer: {note}")


def test_an_unknown_mondo_curie_does_not_resolve_as_primary(tmp_path):
    """MondoResolver checked the OBSOLETE set and not the live one, so
    any well-formed MONDO CURIE came back "primary" -- including one
    the ontology has never heard of. A typo in a curated identifier
    resolved cleanly and became a Disease node joining to nothing.
    """
    r = _res(tmp_path, LIVE)
    assert r.resolve("MONDO:0005502").status == "primary"
    bad = r.resolve("MONDO:9999999")
    assert bad.status == "unmapped" and bad.curie is None


def test_an_empty_ontology_does_not_reject_everything(tmp_path):
    """The existence check keys off a loaded ontology. With none
    loaded, MondoResolver is a pass-through and must stay one rather
    than rejecting every identifier."""
    r = MondoResolver(None)
    assert r.resolve("MONDO:0005502").status == "primary"
