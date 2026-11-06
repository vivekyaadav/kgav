"""Disease nodes carry the disease's name, not the last virus's.

Every Disease node in v0.1 was labelled `f"{virus} disease"`, built from
the virus being ingested. Emit.node UPDATES properties on a repeated id,
and two of the ten configured MONDO ids are shared by four viruses each,
so the surviving label was whichever virus config/viruses.yaml listed
last:

    MONDO:0005709   was "HCoV-229E/NL63/OC43/HKU1 disease", last wins
                    -> shipped as "HCoV-HKU1 disease"
                    -> MONDO calls it "common cold"
    MONDO:0005502   shipped as "DENV-4 disease"
                    -> MONDO calls it "dengue disease"

Six of sixteen viruses pointed at a disease node named after a different
virus. The other four labels were invented where a real one existed:
"SARS-CoV-2 disease" for what the ontology calls COVID-19.

causative_taxon had the same problem and worse consequences: one string,
overwritten per virus, so the common cold node asserted a single
causative agent while carrying four CAUSES edges. The node contradicted
its own edges.

data/raw/mondo/mondo.json had been on disk since 2026-10-28 and the
2026-11-07 audit found that no ingest opened it.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from kgav.mondo import load_labels  # noqa: E402

REAL_MONDO = ROOT / "data" / "raw" / "mondo" / "mondo.json"


def _tiny_mondo(tmp_path: Path, terms: dict[str, str]) -> Path:
    p = tmp_path / "mondo.json"
    p.write_text(json.dumps({"graphs": [{"nodes": [
        {"id": f"http://purl.obolibrary.org/obo/{k.replace(':', '_')}",
         "lbl": v, "type": "CLASS"} for k, v in terms.items()]}]}))
    return p


# ------------------------------------------------------------- the reader
def test_it_reads_the_label_for_an_id_that_was_asked_for(tmp_path):
    p = _tiny_mondo(tmp_path, {"MONDO:0005709": "common cold",
                               "MONDO:0005502": "dengue disease"})
    assert load_labels(p, {"MONDO:0005709"}) == {"MONDO:0005709": "common cold"}


def test_an_id_the_release_does_not_define_is_absent_not_guessed(tmp_path):
    p = _tiny_mondo(tmp_path, {"MONDO:0005709": "common cold"})
    assert load_labels(p, {"MONDO:9999999"}) == {}


def test_a_term_with_no_label_is_not_returned_as_empty(tmp_path):
    p = tmp_path / "m.json"
    p.write_text(json.dumps({"graphs": [{"nodes": [
        {"id": "http://purl.obolibrary.org/obo/MONDO_0005709"}]}]}))
    assert load_labels(p, {"MONDO:0005709"}) == {}


def test_non_mondo_terms_are_ignored(tmp_path):
    """MONDO imports from other ontologies; their IRIs are in the same file."""
    p = tmp_path / "m.json"
    p.write_text(json.dumps({"graphs": [{"nodes": [
        {"id": "http://purl.obolibrary.org/obo/NCBITaxon_9606", "lbl": "human"},
        {"id": "http://purl.obolibrary.org/obo/MONDO_0005709",
         "lbl": "common cold"}]}]}))
    assert load_labels(p, {"MONDO:0005709", "NCBITaxon:9606"}) == {
        "MONDO:0005709": "common cold"}


# ---------------------------------------------------- what the spine emits
VIRUSES = [
    {"taxon": "11137", "name": "HCoV-229E", "mondo": "MONDO:0005709"},
    {"taxon": "277944", "name": "HCoV-NL63", "mondo": "MONDO:0005709"},
    {"taxon": "31631", "name": "HCoV-OC43", "mondo": "MONDO:0005709"},
    {"taxon": "290028", "name": "HCoV-HKU1", "mondo": "MONDO:0005709"},
    {"taxon": "2697049", "name": "SARS-CoV-2", "mondo": "MONDO:0100096"},
    {"taxon": "11060", "name": "DENV-2"},          # no mondo configured
]


def _facts(tmp_path):
    from ingest_spine import disease_facts
    p = _tiny_mondo(tmp_path, {"MONDO:0005709": "common cold",
                               "MONDO:0100096": "COVID-19"})
    return disease_facts(VIRUSES, p)


def test_the_label_is_the_disease_s_not_the_last_virus_s(tmp_path):
    f = _facts(tmp_path)
    assert f["MONDO:0005709"]["label"] == "common cold"
    assert f["MONDO:0100096"]["label"] == "COVID-19"


def test_a_disease_four_viruses_cause_asserts_no_single_agent(tmp_path):
    """causative_taxon holds ONE string. The common cold node carried four
    CAUSES edges and a causative_taxon naming whichever virus was ingested
    last, so the node contradicted its own edges."""
    f = _facts(tmp_path)
    assert f["MONDO:0005709"]["causative_taxon"] is None


def test_a_disease_one_virus_causes_keeps_its_agent(tmp_path):
    f = _facts(tmp_path)
    assert f["MONDO:0100096"]["causative_taxon"] == "NCBITaxon:2697049"


def test_a_virus_with_no_configured_mondo_contributes_nothing(tmp_path):
    assert len(_facts(tmp_path)) == 2


def test_an_unresolvable_id_is_left_out_so_the_caller_sees_the_gap(tmp_path):
    from ingest_spine import disease_facts
    p = _tiny_mondo(tmp_path, {"MONDO:0100096": "COVID-19"})
    f = disease_facts(VIRUSES, p)
    assert set(f) == {"MONDO:0100096"}


def test_without_the_ontology_the_node_is_skipped_not_invented(tmp_path):
    """The fallback IS the bug. A build missing ten nodes is visible; one
    carrying ten wrong labels is not."""
    import gzip

    from ingest_spine import ingest_virus
    from kgav.emit import Emit
    # ingest_virus returns early when the proteome file is absent, BEFORE it
    # reaches the disease block, so an empty one has to exist for this test
    # to be about diseases at all.
    with gzip.open(tmp_path / "UP000000000.json.gz", "wt") as fh:
        json.dump({"results": []}, fh)
    em = Emit()
    stats = ingest_virus(
        em, {"taxon": "290028", "name": "HCoV-HKU1", "family": "Coronaviridae",
             "baltimore_class": "IV", "is_enveloped": True,
             "proteome_id": "UP000000000",   # absent, so no protein is read
             "mondo": "MONDO:0005709"},
        {"family": "Coronaviridae"},
        {"proteome_source": "infores:uniprot", "disease_source": "infores:mondo",
         "evidence_tier": 1},
        tmp_path, diseases=None)
    assert stats["diseases"] == 0
    assert stats["disease_unlabelled"] == 1
    assert not [n for n in em.nodes.values() if n["class"] == "Disease"]
    assert any("no_mondo_label" in k for k in em.notes)
    assert not any("HKU1 disease" in json.dumps(n) for n in em.nodes.values())


# ------------------------------------------------------- against the real file
@pytest.mark.skipif(not REAL_MONDO.exists(), reason="no MONDO release on disk")
def test_every_configured_id_resolves_in_the_real_release():
    """A pinned id MONDO has since deprecated would silently cost a Disease
    node and the CAUSES edge with it."""
    import yaml
    from ingest_spine import disease_facts
    cfg = yaml.safe_load((ROOT / "config" / "viruses.yaml").read_text())
    want = {v["mondo"] for v in cfg["viruses"] if v.get("mondo")}
    got = disease_facts(cfg["viruses"], REAL_MONDO)
    assert set(got) == want, f"unresolved: {sorted(want - set(got))}"


@pytest.mark.skipif(not REAL_MONDO.exists(), reason="no MONDO release on disk")
def test_no_label_is_a_virus_name_with_disease_appended():
    """The exact shape of the bug, asserted against the real ontology."""
    import yaml
    from ingest_spine import disease_facts
    cfg = yaml.safe_load((ROOT / "config" / "viruses.yaml").read_text())
    names = {v["name"] for v in cfg["viruses"]}
    for mondo, f in disease_facts(cfg["viruses"], REAL_MONDO).items():
        assert f["label"] not in {f"{n} disease" for n in names}, \
            f"{mondo} is still labelled from a virus name"
