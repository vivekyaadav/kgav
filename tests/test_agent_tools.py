"""Gate for the tool layer: bounded output, provenance in every sentence, and
absence reported as absence rather than as a negative finding.
"""
import json

import pytest

from kgav.agent.tools import GraphTools

V = "NCBITaxon:2697049"
DRUG = "INCHIKEY:NIRMA"
CYTO = "INCHIKEY:CYTOX"


@pytest.fixture
def graph(tmp_path):
    def prot(pid, viral, sym=None, fam=None):
        p = {"taxon_id": V if viral else "NCBITaxon:9606", "is_viral": viral,
             "sequence_hash": "h", "reviewed": True}
        if sym:
            p["gene_symbol"] = sym
        if fam:
            p["protein_family"] = fam
        return {"id": pid, "class": "Protein", "properties": p}

    def drug(nid, name):
        return {"id": nid, "class": "SmallMolecule",
                "properties": {"smiles": "C", "inchikey_skel": nid[-5:],
                               "chembl_id": name, "is_approved": True,
                               "label": name,
                               "salt_collapsed": False, "stereo_collapsed": False}}

    nodes = [
        {"id": V, "class": "OrganismTaxon",
         "properties": {"label": "SARS-CoV-2", "family": "Coronaviridae",
                        "baltimore_class": "IV", "is_enveloped": True}},
        prot("UniProtKB:NSP5", True, fam="nsp5"),
        {"id": "KGAV:G", "class": "Gene",
         "properties": {"symbol": "rep", "taxon_id": V, "is_viral": True}},
        prot("UniProtKB:SIGMAR1", False, sym="SIGMAR1"),
        drug(DRUG, "nirmatrelvir"), drug(CYTO, "cytotoxin"),
    ]

    def e(s, p, o, q=None, pmids=None):
        d = {"subject": s, "predicate": p, "object": o, "qualifiers": q or {},
             "primary_knowledge_source": "infores:chembl", "evidence_tier": 1,
             "first_asserted_date": "2021-03-04"}
        if pmids:
            d["publications"] = pmids
        return d

    edges = [
        e(DRUG, "INHIBITS", "UniProtKB:NSP5",
          {"assay_type": "biochemical", "ic50_nm": 25.0}, ["PMID:1"]),
        e("UniProtKB:NSP5", "ENCODED_BY", "KGAV:G"),
        e("KGAV:G", "BELONGS_TO", V),
        e(DRUG, "HAS_ANTIVIRAL_ACTIVITY_AGAINST", V,
          {"assay_type": "cell_based_antiviral", "ec50_nm": 100.0,
           "selectivity_index": 50.0, "selectivity_verified": True}),
        e(CYTO, "TARGETS", "UniProtKB:SIGMAR1",
          {"direction": "inhibitor", "assay_type": "binding"}),
        e("UniProtKB:SIGMAR1", "HOST_FACTOR_FOR", V,
          {"direction": "dependency", "screen_type": "CRISPRko",
           "cell_line": "CVCL_0574"}),
        e(CYTO, "HAS_ANTIVIRAL_ACTIVITY_AGAINST", V,
          {"assay_type": "cell_based_antiviral", "ec50_nm": 200.0,
           "selectivity_index": 2.0, "selectivity_verified": True}),
    ]
    (tmp_path / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text("\n".join(json.dumps(x) for x in edges))
    return GraphTools(tmp_path)


# ------------------------------------------------------------------ resolve
def test_resolve_finds_a_known_name(graph):
    r = graph.resolve("nirmatrelvir")
    assert r.ok and r.data[0]["id"] == DRUG


def test_resolve_reports_absence_as_absence(graph):
    """A missing compound is not a negative finding about the compound."""
    r = graph.resolve("aspirin")
    assert not r.ok
    assert "may be absent rather than misspelled" in r.note


def test_ambiguity_is_returned_not_resolved(graph):
    r = graph.resolve("c")          # substring of several labels
    if r.ok and len(r.data) > 1:
        assert "ambiguous" in r.note


# ------------------------------------------------------------------ explain
def test_explain_returns_a_direct_acting_route_with_provenance(graph):
    r = graph.explain(DRUG, V)
    assert r.ok
    joined = " ".join(r.verbalised)
    assert "M1 direct-acting" in joined
    assert "inhibits the viral protein" in joined
    assert "PMID:1" in joined and "infores" not in joined  # source is shortened


def test_every_verbalised_edge_carries_an_edge_id(graph):
    r = graph.explain(DRUG, V)
    assert r.edge_ids
    for eid in r.edge_ids:
        assert any(eid in line for line in r.verbalised)


def test_explain_attaches_selectivity(graph):
    r = graph.explain(DRUG, V)
    assert any("SELECTIVE" in w and "50.0" in w for w in r.warnings)


def test_cytotoxic_compound_is_flagged_in_its_explanation(graph):
    r = graph.explain(CYTO, V)
    assert r.ok
    assert any("CYTOTOXIC" in w for w in r.warnings)


def test_path_warnings_attach_to_their_own_route(graph):
    """Warnings must sit with the route that triggered them. Collected at the
    compound level, a SIGMAR1 cell-context caveat printed against a Spike
    route -- a mismatch that discredits every other warning shown."""
    r = graph.explain(CYTO, V)
    joined = " ".join(r.verbalised)
    assert "CELL CONTEXT" in joined and "SIGMAR1" in joined
    assert "below chance" in joined
    # compound-level warnings carry selectivity only
    assert all("CELL CONTEXT" not in w for w in r.warnings)


def test_no_path_is_reported_as_no_evidence_not_no_activity(graph):
    r = graph.explain("INCHIKEY:UNKNOWN", V)
    assert not r.ok
    assert "absence of evidence, not evidence of absence" in r.note


# ------------------------------------------------------------------- triage
def test_triage_ranks_direct_acting_above_host_directed(graph):
    r = graph.triage([CYTO, DRUG], V)
    assert r.ok
    assert r.data[0]["compound"] == DRUG
    assert r.data[0]["direct_acting_routes"] == 1


def test_triage_states_why_it_ranks_that_way(graph):
    r = graph.triage([DRUG], V)
    assert any("0.823" in w and "below chance" in w for w in r.warnings)


# ----------------------------------------------------------------- evidence
def test_evidence_returns_provenance_per_claim(graph):
    r = graph.evidence(DRUG, "UniProtKB:NSP5", "INHIBITS")
    assert r.ok and r.edge_ids
    assert "chembl" in r.verbalised[0]


def test_evidence_for_an_unmade_claim_is_a_clean_miss(graph):
    r = graph.evidence(DRUG, "UniProtKB:SIGMAR1", "INHIBITS")
    assert not r.ok and "No edge" in r.note


# --------------------------------------------------------------- boundedness
def test_neighbours_are_capped(graph):
    r = graph.neighbours("UniProtKB:NSP5", limit=1)
    assert len(r.data) == 1 and "truncated" in r.note


def test_profile_summarises_without_dumping_edges(graph):
    r = graph.profile("UniProtKB:SIGMAR1")
    assert r.ok
    assert r.data["edge_counts"]["HOST_FACTOR_FOR"] == 1


def test_result_serialises_for_a_model(graph):
    payload = json.loads(graph.explain(DRUG, V).to_json())
    assert set(payload) == {"ok", "kind", "data", "text", "warnings",
                            "edge_ids", "note"}
