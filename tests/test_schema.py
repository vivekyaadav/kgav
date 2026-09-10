"""Day 1 acceptance gate.

The gate: the validator loads the YAML and rejects a deliberately malformed
triple. Everything else guards invariants the rest of the pipeline assumes.
"""
import pytest

from kgav.schema import load_schema


@pytest.fixture(scope="module")
def schema():
    return load_schema()


def _mol(nid="INCHIKEY:RWWYLEGWBNMMLJ-YSOARWBDSA-N", **over):
    props = {"smiles": "CCO", "inchikey_skel": "RWWYLEGWBNMM", "is_approved": True,
             "salt_collapsed": True, "stereo_collapsed": False}
    props.update(over)
    return {"id": nid, "class": "SmallMolecule", "properties": props}


def _protein(nid, is_viral, taxon="NCBITaxon:2697049", **over):
    props = {"taxon_id": taxon, "is_viral": is_viral, "sequence_hash": "ab" * 16,
             "reviewed": True}
    props.update(over)
    return {"id": nid, "class": "Protein", "properties": props}


def _virus(nid="NCBITaxon:2697049"):
    return {"id": nid, "class": "OrganismTaxon",
            "properties": {"label": "SARS-CoV-2", "family": "Coronaviridae",
                           "baltimore_class": "IV", "is_enveloped": True}}


def _edge(subject, predicate, obj, qualifiers=None, **over):
    e = {"subject": subject, "predicate": predicate, "object": obj,
         "qualifiers": qualifiers or {},
         "primary_knowledge_source": "infores:chembl",
         "evidence_tier": 1,
         "first_asserted_date": "2021-03-04"}
    e.update(over)
    return e


VIRAL_P = "UniProtKB:P0DTD1"
HOST_P = "UniProtKB:O15393"


@pytest.fixture
def world():
    return {"nodes": [_mol(), _protein(VIRAL_P, True),
                      _protein(HOST_P, False, taxon="NCBITaxon:9606"), _virus()]}


def _index(nodes):
    return {n["id"]: n for n in nodes}


def test_schema_loads(schema):
    assert schema.version
    assert schema.node_classes and schema.edge_classes
    assert "first_asserted_date" in schema.required_provenance


def test_schema_is_internally_consistent(schema):
    classes = set(schema.node_classes)
    preds = {ec.predicate for ec in schema.edge_classes}
    for name, mp in schema.metapaths.items():
        path = mp["path"]
        assert len(path) % 2 == 1, f"{name}: path must alternate node/edge/node"
        for i, step in enumerate(path):
            bucket = classes if i % 2 == 0 else preds
            assert step in bucket, f"{name}: {step!r} undeclared"


def test_enums_referenced_by_qualifiers_exist(schema):
    for ec in schema.edge_classes:
        for qname, spec in ec.qualifiers.items():
            if "enum" in spec:
                assert spec["enum"] in schema.enums, f"{ec.predicate}.{qname}"


def test_held_out_predicates_are_declared(schema):
    assert {"TREATS", "IN_TRIAL_FOR", "HAS_ANTIVIRAL_ACTIVITY_AGAINST"} <= schema.held_out_predicates()


def test_valid_direct_acting_edge(schema, world):
    e = _edge(_mol()["id"], "INHIBITS", VIRAL_P,
              {"ec50_nm": 120.0, "cc50_nm": 24000.0, "assay_type": "biochemical",
               "lifecycle_stage": "genome_replication"})
    assert schema.validate_edge(e, _index(world["nodes"])) == []


def test_valid_host_directed_edge(schema, world):
    e = _edge(HOST_P, "HOST_FACTOR_FOR", "NCBITaxon:2697049",
              {"direction": "dependency", "screen_type": "CRISPRko",
               "cell_line": "CVCL_0574", "effect_size": -2.4})
    assert schema.validate_edge(e, _index(world["nodes"])) == []


def test_valid_node(schema):
    assert schema.validate_node(_mol()) == []


def test_gate_rejects_malformed_triple(schema, world):
    """A drug TARGETS edge pointed at a VIRAL protein must be rejected.

    TARGETS is the host-directed predicate; INHIBITS is direct-acting.
    Topologically the triple looks fine. Only the schema catches it.
    """
    bad = _edge(_mol()["id"], "TARGETS", VIRAL_P, {"direction": "inhibitor"})
    v = schema.validate_edge(bad, _index(world["nodes"]))
    assert v, "validator accepted a triple that violates the object constraint"
    assert any(x.code == "CONSTRAINT" for x in v), [str(x) for x in v]


def test_gate_rejects_unknown_predicate(schema, world):
    bad = _edge(_mol()["id"], "CURES_EVERYTHING", "NCBITaxon:2697049")
    assert any(x.code == "UNKNOWN_PREDICATE" for x in schema.validate_edge(bad, _index(world["nodes"])))


def test_gate_rejects_undeclared_class_combination(schema, world):
    bad = _edge("NCBITaxon:2697049", "CAUSES", HOST_P, {})
    v = schema.validate_edge(bad, _index(world["nodes"]))
    assert any(x.code == "UNDECLARED_TRIPLE" for x in v), [str(x) for x in v]


def test_edge_without_date_is_rejected(schema, world):
    bad = _edge(HOST_P, "HOST_FACTOR_FOR", "NCBITaxon:2697049",
                {"direction": "dependency", "screen_type": "CRISPRko", "cell_line": "CVCL_0574"},
                first_asserted_date=None)
    v = schema.validate_edge(bad, _index(world["nodes"]))
    assert any(x.code == "MISSING_PROVENANCE" and "first_asserted_date" in x.message for x in v)


def test_edge_without_source_is_rejected(schema, world):
    bad = _edge(_mol()["id"], "INHIBITS", VIRAL_P, {"assay_type": "biochemical"},
                primary_knowledge_source=None)
    assert any(x.code == "MISSING_PROVENANCE" for x in schema.validate_edge(bad, _index(world["nodes"])))


def test_bad_date_format_is_rejected(schema, world):
    bad = _edge(_mol()["id"], "INHIBITS", VIRAL_P, {"assay_type": "biochemical"},
                first_asserted_date="March 2021")
    assert any(x.code == "TYPE" for x in schema.validate_edge(bad, _index(world["nodes"])))


def test_host_factor_without_direction_is_rejected(schema, world):
    """The single most dangerous omission in the whole schema."""
    bad = _edge(HOST_P, "HOST_FACTOR_FOR", "NCBITaxon:2697049",
                {"screen_type": "CRISPRko", "cell_line": "CVCL_0574"})
    v = schema.validate_edge(bad, _index(world["nodes"]))
    assert any(x.code == "MISSING_QUALIFIER" and "direction" in x.message for x in v)


def test_bad_enum_value_is_rejected(schema, world):
    bad = _edge(HOST_P, "HOST_FACTOR_FOR", "NCBITaxon:2697049",
                {"direction": "helpful", "screen_type": "CRISPRko", "cell_line": "CVCL_0574"})
    assert any(x.code == "ENUM" for x in schema.validate_edge(bad, _index(world["nodes"])))


def test_undeclared_qualifier_is_rejected(schema, world):
    bad = _edge(_mol()["id"], "INHIBITS", VIRAL_P, {"assay_type": "biochemical", "vibes": "good"})
    assert any(x.code == "UNDECLARED_QUALIFIER" for x in schema.validate_edge(bad, _index(world["nodes"])))


def test_similarity_below_threshold_is_rejected(schema):
    nodes = [_mol(), _mol("INCHIKEY:AAAAAAAAAAAAAA-BBBBBBBBBB-C")]
    bad = _edge(nodes[0]["id"], "CHEMICALLY_SIMILAR_TO", nodes[1]["id"], {"tanimoto_ecfp4": 0.31})
    assert any(x.code == "RANGE" for x in schema.validate_edge(bad, _index(nodes)))


def test_model_output_cannot_carry_evidence_tier(schema, world):
    bad = _edge(_mol()["id"], "PREDICTED_INDICATION", "NCBITaxon:2697049",
                {"score": 0.91, "rank": 3, "model_version": "dwpc-0.1"})
    assert any(x.code == "MODEL_OUTPUT_IN_EVIDENCE" for x in schema.validate_edge(bad, _index(world["nodes"])))


def test_wrong_prefix_for_class_is_rejected(schema):
    assert any(x.code == "PREFIX" for x in schema.validate_node(_mol(nid="CHEMBL:CHEMBL521")))


def test_missing_required_property_is_rejected(schema):
    bad = _mol()
    del bad["properties"]["is_approved"]
    assert any(x.code == "MISSING_PROPERTY" for x in schema.validate_node(bad))


def test_dangling_edge_is_rejected(schema, world):
    bad = _edge(_mol()["id"], "INHIBITS", "UniProtKB:NOPE", {"assay_type": "biochemical"})
    assert any(x.code == "DANGLING" for x in schema.validate_edge(bad, _index(world["nodes"])))


def test_batch_validation_catches_duplicate_nodes(schema):
    assert any(x.code == "DUPLICATE_NODE" for x in schema.validate_batch([_mol(), _mol()], []))
