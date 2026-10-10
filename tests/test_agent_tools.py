"""Gate for the tool layer: bounded output, provenance in every sentence, and
absence reported as absence rather than as a negative finding.
"""
import json
from pathlib import Path

import pytest

from kgav.assemble import fact_id
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


@pytest.fixture
def graph_dir(graph, tmp_path):
    """The directory the `graph` fixture wrote, for tests that rebuild it."""
    return tmp_path


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
    # The host-directed caveat must be attached. WHAT it says comes from the
    # calibration -- asserting the literal "below chance" is how this test
    # passed while the figure inside it went stale.
    assert "HOST-DIRECTED ROUTE" in joined


def test_path_warnings_are_also_required_warnings(graph):
    """Printing a caveat is not enforcing it.

    This assertion used to read `all("CELL CONTEXT" not in w for w in
    r.warnings)` -- it locked in the defect. path_warnings went only into
    `verbalised`, so cell_context and host_directed_route never reached
    ToolResult.warnings, never reached the brief, and verify()'s
    warning-survival check never required them. A model that dropped the
    SIGMAR1 caveat still passed verification.

    Sitting with its route and being required are both true of the same
    caveat: the test above checks the first, this one the second.
    """
    r = graph.explain(CYTO, V)
    assert any("CELL CONTEXT" in w and "SIGMAR1" in w for w in r.warnings)
    assert any("HOST-DIRECTED ROUTE" in w for w in r.warnings)


def test_direct_acting_route_requires_no_cell_context_warning(graph):
    """The converse, so the fix above cannot become 'warn on everything'."""
    r = graph.explain(DRUG, V)
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
    # NOT "0.826": that figure came from a protocol which filtered positives
    # and left negatives whole, and the same metapath reads 0.499 under the
    # publication-disjoint control. The caveat states what was measured, or
    # that nothing was.
    assert any("direct-acting" in w for w in r.warnings)
    assert not any("0.826" in w for w in r.warnings)


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


REAL = __import__("pathlib").Path(__file__).resolve().parents[1] / "data/releases/v0.1"


@pytest.mark.skipif(not (REAL / "nodes.jsonl").exists(), reason="release not built")
def test_triage_ranks_the_named_controls_correctly():
    """The ranking has been wrong three ways, and each time the aggregate
    output looked reasonable while the controls were inverted:

      v1 counting routes        chloroquine above nirmatrelvir
      v2 potency before SI      chloroquine above remdesivir
      v3 a 1 uM potency gate    remdesivir demoted below chloroquine for
                                sitting the wrong side of an arbitrary line

    The rule may use only quantities this project measured: direct-acting
    evidence, then selectivity, then potency as a tiebreak.
    """
    g = GraphTools(REAL)
    ids = {n: g.resolve(n).data[0]["id"]
           for n in ("nirmatrelvir", "remdesivir", "chloroquine",
                     "hydroxychloroquine")}
    r = g.triage(list(ids.values()), "NCBITaxon:2697049")
    order = [row["label"] for row in r.data]

    assert order.index("nirmatrelvir") < order.index("chloroquine")
    assert order.index("remdesivir") < order.index("chloroquine")
    assert order.index("chloroquine") < order.index("hydroxychloroquine")


@pytest.mark.skipif(not (REAL / "nodes.jsonl").exists(), reason="release not built")
def test_triage_reports_compounds_it_could_not_rank():
    """Silently returning three of four requested compounds is a failure the
    user cannot detect."""
    g = GraphTools(REAL)
    real_id = g.resolve("nirmatrelvir").data[0]["id"]
    r = g.triage([real_id, "INCHIKEY:DOESNOTEXIST"], "NCBITaxon:2697049")
    assert len(r.data) == 1
    assert any("not in this graph" in w for w in r.warnings)


# -------------------------------------------- citations that hold still
def _edges_of(release: Path) -> list[dict]:
    return [json.loads(l) for l in
            (release / "edges.jsonl").read_text().splitlines() if l.strip()]


def _reordered(src: Path, dst: Path, edges: list[dict]) -> GraphTools:
    """Same nodes, same facts, a different line order in edges.jsonl."""
    dst.mkdir(parents=True, exist_ok=True)
    (dst / "nodes.jsonl").write_text((src / "nodes.jsonl").read_text())
    (dst / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in edges))
    return GraphTools(dst)


def test_without_an_edge_id_a_citation_is_a_line_number(tmp_path, graph_dir):
    """THE BUG THIS FIXES. tools.py falls back to f"E{i}", the edge's position
    in edges.jsonl. Reorder the file -- which reassembly does whenever a layer
    changes -- and the same fact answers to a different id, so a citation in a
    saved answer silently comes to point at something else.
    """
    raw = _edges_of(graph_dir)
    for e in raw:
        e.pop("edge_id", None)
    a = _reordered(graph_dir, tmp_path / "a", raw)
    b = _reordered(graph_dir, tmp_path / "b", list(reversed(raw)))
    ia, ib = a.explain(DRUG, V).edge_ids, b.explain(DRUG, V).edge_ids
    assert ia and ib
    assert ia != ib, "line-index ids agreed by luck; the test proves nothing"


def test_with_an_edge_id_the_citation_survives_a_reorder(tmp_path, graph_dir):
    """Same facts, same ids, whatever order the file is in."""
    raw = _edges_of(graph_dir)
    for e in raw:
        e["edge_id"] = fact_id((e["subject"], e["predicate"], e["object"]))
    c = _reordered(graph_dir, tmp_path / "c", raw)
    d = _reordered(graph_dir, tmp_path / "d", list(reversed(raw)))
    ic, idd = c.explain(DRUG, V).edge_ids, d.explain(DRUG, V).edge_ids
    assert ic and set(ic) == set(idd)
    assert all(x.startswith("E") for x in ic)


def test_the_agents_ids_pass_its_own_verifier(graph):
    """An id the verifier cannot parse makes every answer fail as an invented
    citation, so the two have to agree on the shape."""
    from kgav.agent.verify import EDGE_ID_RE
    for eid in graph.explain(DRUG, V).edge_ids:
        assert EDGE_ID_RE.findall(f"claim [{eid}].") == [eid], eid
