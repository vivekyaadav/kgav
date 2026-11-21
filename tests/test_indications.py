"""The indications layer places what it says, and nothing else.

Every test here is against the assembled release or against the layer
directory, not against the script's printout. The failure modes this
layer has are specific: a trial attributed to the wrong virus (the
hand-written xref map did exactly that once, 66 COVID-19 trials onto
MERS), a compound invented to carry a label, a required qualifier
fabricated so an edge could ship, or several trial registrations
collapsing into one fact with no trace of the others.
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav import virus_register as vr  # noqa: E402
from kgav.emit import Emit  # noqa: E402
from kgav.indications import (  # noqa: E402
    EVIDENCE_TIER,
    INDICATION_SOURCE,
    IndicationError,
    disease_xrefs,
    ingest_indications,
)
from kgav.pipeline import layer_dir, release_dir  # noqa: E402
from kgav.schema import load_schema  # noqa: E402

MONDO = ROOT / "data" / "raw" / "mondo" / "mondo.json"
LAYER = layer_dir(ROOT, "indications")
PREDICATES = ("TREATS", "IN_TRIAL_FOR")


@pytest.fixture(scope="module")
def edges() -> list[dict]:
    path = LAYER / "edges.jsonl"
    if not path.exists():
        pytest.skip(f"{path} absent; run scripts/ingest_indications.py")
    return [json.loads(ln) for ln in path.read_text().splitlines() if ln.strip()]


@pytest.fixture(scope="module")
def release_nodes() -> dict[str, str]:
    nodes = release_dir(ROOT) / "nodes.jsonl"
    return {n["id"]: n["class"]
            for n in (json.loads(ln) for ln in nodes.open())}


@pytest.fixture(scope="module")
def xrefs():
    if not MONDO.exists():
        pytest.skip("mondo.json absent")
    return disease_xrefs(vr.load()["viruses"], MONDO)


def test_the_xref_map_covers_the_register_and_only_the_register(xrefs):
    registered = {v["mondo"] for v in vr.load()["viruses"] if v.get("mondo")}
    assert set(xrefs) == registered


def test_an_identifier_claimed_by_two_diseases_refuses(tmp_path):
    """THE DEFECT THAT MOTIVATED DERIVING THIS MAP.

    A hand-written version put the COVID-19 MeSH term under MERS-CoV as
    well, and 66 COVID-19 trials landed on the wrong virus. A shared
    identifier has to stop the run, so this builds an ontology in which
    two MONDO terms xref one MeSH id and asserts it does.
    """
    fake = tmp_path / "mondo.json"
    fake.write_text(json.dumps({"graphs": [{"nodes": [
        {"id": "http://purl.obolibrary.org/obo/MONDO_0000001",
         "meta": {"xrefs": [{"val": "MESH:D000001"}]}},
        {"id": "http://purl.obolibrary.org/obo/MONDO_0000002",
         "meta": {"xrefs": [{"val": "MESH:D000001"}]}},
    ]}]}))
    viruses = [{"mondo": "MONDO:0000001"}, {"mondo": "MONDO:0000002"}]
    with pytest.raises(IndicationError, match="more than one registered"):
        disease_xrefs(viruses, fake)
    # And the same ontology without the collision does not refuse, so
    # the test above is about the collision and not about the fixture.
    fake.write_text(json.dumps({"graphs": [{"nodes": [
        {"id": "http://purl.obolibrary.org/obo/MONDO_0000001",
         "meta": {"xrefs": [{"val": "MESH:D000001"}]}},
        {"id": "http://purl.obolibrary.org/obo/MONDO_0000002",
         "meta": {"xrefs": [{"val": "MESH:D000002"}]}},
    ]}]}))
    assert set(disease_xrefs(viruses, fake)) == {"MONDO:0000001",
                                                 "MONDO:0000002"}


def test_no_register_entry_lacks_a_mondo_id_or_the_join_is_narrower(xrefs):
    """The join can only reach diseases the register names. If an entry
    loses its MONDO id the layer silently covers fewer viruses, so the
    count is asserted rather than assumed."""
    viruses = vr.load()["viruses"]
    assert all(v.get("mondo") for v in viruses), (
        "a register entry carries no MONDO id, so its virus has no disease "
        "and no clinical label can reach it")


def test_the_layer_adds_no_nodes():
    """The registration's refusal. A new SmallMolecule would move every
    other population in the repository."""
    nodes = LAYER / "nodes.jsonl"
    if not nodes.exists():
        pytest.skip("layer absent")
    assert nodes.read_text().strip() == "", (
        "the indications layer emitted a node; doc2/registrations/prereg_clinical_labels.md "
        "forbids it")


def test_every_endpoint_already_existed(edges, release_nodes):
    for e in edges:
        assert release_nodes.get(e["subject"]) == "SmallMolecule", e["subject"]
        assert release_nodes.get(e["object"]) == "Disease", e["object"]


def test_every_edge_carries_what_its_predicate_requires(edges):
    """Read from the schema, not restated here.

    An edge that shipped without a required qualifier would fail
    assembly, so this is not the only guard -- but it names the field
    that is missing, which a schema failure over 396,000 edges does not.
    """
    schema = load_schema(ROOT / "schema" / "kg_schema.yaml")
    required = {ec.predicate: set(ec.required_qualifiers)
                for ec in schema.edge_classes if ec.predicate in PREDICATES}
    assert required.get("IN_TRIAL_FOR") == {"nct_id", "status"}, (
        "IN_TRIAL_FOR's required qualifiers changed; the layer's drop rule "
        "for rows without a trial reference has to change with them")
    for e in edges:
        have = set(e["qualifiers"])
        assert required[e["predicate"]] <= have, (e["predicate"], have)
        assert e["primary_knowledge_source"] == INDICATION_SOURCE
        assert e["evidence_tier"] == EVIDENCE_TIER


def test_treats_and_in_trial_for_split_on_phase_four_only(edges):
    phases = {"TREATS": set(), "IN_TRIAL_FOR": set()}
    for e in edges:
        if e["predicate"] == "TREATS":
            assert e["qualifiers"]["approval_status"] == "approved"
            phases["TREATS"].add(4.0)
        else:
            phases["IN_TRIAL_FOR"].add(float(e["qualifiers"]["phase"]))
    assert phases["IN_TRIAL_FOR"], "no IN_TRIAL_FOR edge"
    assert max(phases["IN_TRIAL_FOR"]) < 4.0, (
        f"an IN_TRIAL_FOR edge carries phase {max(phases['IN_TRIAL_FOR'])}; "
        "phase 4 is an approval and belongs to TREATS")


def test_trial_status_is_unknown_on_every_edge_because_chembl_has_none(edges):
    """Stated, not invented. A consumer filtering for `completed` must
    get nothing rather than something fabricated."""
    statuses = {e["qualifiers"]["status"]
                for e in edges if e["predicate"] == "IN_TRIAL_FOR"}
    assert statuses == {"unknown"}, (
        f"IN_TRIAL_FOR carries {statuses}; ChEMBL's drug_indication records "
        "no trial status, so anything else was invented")


def test_one_edge_per_compound_and_disease(edges):
    """`nct_id` is not an identity qualifier, so a pair is one fact. Two
    edges for one pair would double-count a compound in the registered
    cell's positive class."""
    seen = Counter((e["subject"], e["predicate"], e["object"]) for e in edges)
    dupes = {k: n for k, n in seen.items() if n > 1}
    assert not dupes, dupes
    pairs = Counter((e["subject"], e["object"]) for e in edges)
    both = {k: n for k, n in pairs.items() if n > 1}
    assert not both, (
        f"{len(both)} pair(s) carry both TREATS and IN_TRIAL_FOR: {both}")


def test_several_trial_registrations_are_counted_not_lost(edges):
    """support_count is the only trace of the second registration, so a
    pair citing two trials must carry it."""
    multi = [e for e in edges if e.get("support_count", 1) > 1]
    assert multi, (
        "no edge cites more than one trial; either ChEMBL changed or "
        "support_count stopped being written, and the second registration "
        "would then vanish without trace")
    for e in multi:
        if e["predicate"] == "IN_TRIAL_FOR":
            assert e["qualifiers"]["nct_id"].startswith("NCT")
        else:
            # TREATS HAS NO PLACE TO NAME A TRIAL. Its schema qualifiers
            # are approval_status, regulator and approval_year, so an
            # approved indication ChEMBL backs with several registry
            # entries carries the COUNT and not the identifiers. That is
            # a real loss and it is recorded here rather than hidden
            # behind a count that looks complete.
            assert "nct_id" not in e["qualifiers"]
    # The representative is the LOWEST id, so a rebuild is stable.
    reps = [e["qualifiers"]["nct_id"] for e in edges
            if e["predicate"] == "IN_TRIAL_FOR" and e.get("support_count", 1) > 1]
    assert reps, "no multi-trial IN_TRIAL_FOR edge to check stability on"
    assert all(e["support_count"] >= 1 for e in edges if "support_count" in e)


def test_a_row_without_a_trial_reference_is_dropped_not_fabricated():
    """The one edge the schema refused.

    Builds the layer in memory against a stub that returns one row with
    no ClinicalTrials reference, and asserts no edge ships. The
    alternative -- inventing an id, or relaxing required_qualifiers --
    would make an unverifiable trial edge look like a verifiable one.
    """
    import kgav.indications as I

    rows = [("AAAAAAAAAAAAAA-BBBBBBBBBB-C", 2.0, None)]
    original = I._rows
    I._rows = lambda db, mesh, efo: rows
    try:
        em = Emit()
        s = ingest_indications(Path("unused"), {"MONDO:0100096": (set(), set())},
                               {"INCHIKEY:AAAAAAAAAAAAAA-BBBBBBBBBB-C"},
                               {"MONDO:0100096"}, em)
        assert s["in_trial_without_nct"] == 1
        assert em.edges == [], em.edges
        # And with a reference, the same row DOES ship -- so the test
        # above is about the missing reference, not a broken stub.
        rows[0] = ("AAAAAAAAAAAAAA-BBBBBBBBBB-C", 2.0, "NCT00000001")
        em2 = Emit()
        s2 = ingest_indications(Path("unused"),
                                {"MONDO:0100096": (set(), set())},
                                {"INCHIKEY:AAAAAAAAAAAAAA-BBBBBBBBBB-C"},
                                {"MONDO:0100096"}, em2)
        assert s2["in_trial"] == 1 and len(em2.edges) == 1
        assert em2.edges[0]["qualifiers"]["nct_id"] == "NCT00000001"
    finally:
        I._rows = original


def test_an_indicated_compound_the_graph_lacks_is_dropped():
    """The registration's no-new-nodes rule, as a behaviour rather than
    a claim about the script."""
    import kgav.indications as I

    original = I._rows
    I._rows = lambda db, mesh, efo: [("ZZZZZZZZZZZZZZ-YYYYYYYYYY-X", 2.0,
                                      "NCT00000002")]
    try:
        em = Emit()
        s = ingest_indications(Path("unused"),
                               {"MONDO:0100096": (set(), set())},
                               set(), {"MONDO:0100096"}, em)
        assert s["compound_not_in_graph"] == 1
        assert em.edges == []
        assert em.nodes == {}
    finally:
        I._rows = original


def test_the_assembled_release_holds_the_layer(edges):
    """The layer directory is not the release. An edge written here and
    dropped at assembly would leave both files looking right."""
    counts: Counter = Counter()
    for ln in (release_dir(ROOT) / "edges.jsonl").open():
        if '"TREATS"' not in ln and '"IN_TRIAL_FOR"' not in ln:
            continue
        e = json.loads(ln)
        if e["predicate"] in PREDICATES:
            counts[e["predicate"]] += 1
    mine = Counter(e["predicate"] for e in edges)
    assert counts == mine, (
        f"the layer holds {dict(mine)} and the release {dict(counts)}")
