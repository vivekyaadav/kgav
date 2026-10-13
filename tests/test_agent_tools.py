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
        # chembl.py sets target_resolution on EVERY INHIBITS edge, both
        # branches, and the schema requires it -- so M1 constrains on it and
        # qualifier_ok reads a missing qualifier as a failed match. A fixture
        # omitting it is not a simpler fixture; it is an edge the
        # direct-acting channel refuses to walk, in the scorer and now in the
        # agent too.
        if p == "INHIBITS":
            q = dict(q or {})
            q.setdefault("target_resolution", "chain_from_assay_description")
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


def test_triage_ranks_the_named_controls_correctly(real_tools):
    """The ranking has been wrong three ways, and each time the aggregate
    output looked reasonable while the controls were inverted:

      v1 counting routes        chloroquine above nirmatrelvir
      v2 potency before SI      chloroquine above remdesivir
      v3 a 1 uM potency gate    remdesivir demoted below chloroquine for
                                sitting the wrong side of an arbitrary line

    The rule may use only quantities this project measured: direct-acting
    evidence, then selectivity, then potency as a tiebreak.
    """
    g = real_tools
    ids = {n: g.resolve(n).data[0]["id"]
           for n in ("nirmatrelvir", "remdesivir", "chloroquine",
                     "hydroxychloroquine")}
    r = g.triage(list(ids.values()), "NCBITaxon:2697049")
    order = [row["label"] for row in r.data]

    assert order.index("nirmatrelvir") < order.index("chloroquine")
    assert order.index("remdesivir") < order.index("chloroquine")
    assert order.index("chloroquine") < order.index("hydroxychloroquine")


def test_triage_reports_compounds_it_could_not_rank(real_tools):
    """Silently returning three of four requested compounds is a failure the
    user cannot detect."""
    g = real_tools
    real_id = g.resolve("nirmatrelvir").data[0]["id"]
    r = g.triage([real_id, "INCHIKEY:DOESNOTEXIST"], "NCBITaxon:2697049")
    assert len(r.data) == 1
    assert any("not in this graph" in w for w in r.warnings)


# ------------------------------- the agent and the scorer must walk one graph
@pytest.fixture
def resolution_release(tmp_path):
    """One resolved edge and one unresolved, both reaching the virus.

    The real asymmetry on v0.1: nirmatrelvir's activity resolved to the
    mature nsp5 chain, pomotrelvir's stayed on the parent polyprotein as
    parent_unresolved. M1 requires the former.
    """
    V = "NCBITaxon:2697049"
    nodes = [
        {"id": V, "class": "OrganismTaxon",
         "properties": {"label": "SARS-CoV-2", "family": "Coronaviridae",
                        "baltimore_class": "IV", "is_enveloped": True}},
        {"id": "UniProtKB:CHAIN", "class": "Protein",
         "properties": {"taxon_id": V, "is_viral": True,
                        "protein_family": "nsp5", "sequence_hash": "h",
                        "reviewed": True}},
        {"id": "UniProtKB:PARENT", "class": "Protein",
         "properties": {"taxon_id": V, "is_viral": True,
                        "protein_family": "pp1ab", "sequence_hash": "h2",
                        "reviewed": True}},
        {"id": "KGAV:GENE_rep", "class": "Gene",
         "properties": {"symbol": "rep", "taxon_id": V, "is_viral": True}},
    ]
    for d in ("INCHIKEY:RESOLVED", "INCHIKEY:VAGUE"):
        nodes.append({"id": d, "class": "SmallMolecule",
                      "properties": {"smiles": "C", "inchikey_skel": d[-6:],
                                     "chembl_id": "CH" + d[-3:],
                                     "salt_collapsed": False,
                                     "stereo_collapsed": False}})

    def e(s_, p_, o_, q=None):
        return {"subject": s_, "predicate": p_, "object": o_,
                "qualifiers": dict(q or {}),
                "primary_knowledge_source": "infores:chembl",
                "evidence_tier": 1, "first_asserted_date": "2021-03-04"}

    edges = [
        e("INCHIKEY:RESOLVED", "INHIBITS", "UniProtKB:CHAIN",
          {"assay_type": "biochemical", "ic50_nm": 25.0,
           "target_resolution": "chain_from_assay_description"}),
        e("INCHIKEY:VAGUE", "INHIBITS", "UniProtKB:PARENT",
          {"assay_type": "biochemical", "ic50_nm": 25.0,
           "target_resolution": "parent_unresolved"}),
        e("UniProtKB:CHAIN", "ENCODED_BY", "KGAV:GENE_rep"),
        e("UniProtKB:PARENT", "ENCODED_BY", "KGAV:GENE_rep"),
        e("KGAV:GENE_rep", "BELONGS_TO", V),
    ]
    (tmp_path / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text(
        "\n".join(json.dumps(x) for x in edges))
    return tmp_path


def test_the_agent_refuses_the_route_the_scorer_refuses(resolution_release):
    """THE DIVERGENCE THIS CLOSES. _find_paths enumerated M1 as its own loop
    and applied no constraints, so when M1 gained
    INHIBITS.target_resolution=chain_from_assay_description the scorer
    stopped walking parent_unresolved edges and the agent carried on
    explaining routes through them -- able to narrate a mechanism the ranking
    refuses. Same shape as the agent holding one calibration for every virus:
    a fact declared in one place and consumed in one of two.
    """
    from kgav.baselines import Graph, walk_metapath
    from kgav.schema import load_schema
    schema = load_schema()
    m1 = schema.metapaths["M1"]
    g = Graph.load(resolution_release,
                   skip_predicates=schema.held_out_predicates(),
                   symmetry=schema.symmetry())
    t = GraphTools(resolution_release)
    V = "NCBITaxon:2697049"
    for drug, expected in (("INCHIKEY:RESOLVED", True),
                           ("INCHIKEY:VAGUE", False)):
        scorer = bool(walk_metapath(g, drug, m1["path"],
                                    m1.get("constraints"), set()))
        agent = any(r.startswith("M1")
                    for r, _es in t._find_paths(drug, V, 10))
        assert scorer is expected, f"scorer disagrees for {drug}"
        assert agent is expected, f"agent disagrees for {drug}"
        assert scorer == agent, f"agent and scorer diverge on {drug}"


def test_m2s_constraint_comes_from_the_schema_not_a_hand_copy(real_tools):
    """It was written out here as `q.get("direction") == "dependency"`. A
    hand-copied constraint is a second copy: it happened to match, and it is
    the mechanism by which M1's did not."""
    t = real_tools
    assert t._constraints["M2"] == {"HOST_FACTOR_FOR.direction": "dependency"}
    assert t._constraints["M1"] == {
        "INHIBITS.target_resolution": "chain_from_assay_description"}


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


# --------------------------------------------------------------------------
# one calibration per virus
# --------------------------------------------------------------------------
ZIKA = "NCBITaxon:64320"


class _Ch:
    """A channel that DISCRIMINATES, so a leak takes the most confident branch."""
    def __init__(self, auc=0.675, lo=0.660, hi=0.690):
        self.auc, self.ci_lo, self.ci_hi = auc, lo, hi
        self.n_pos, self.n_neg = 2020, 10956
        self.coverage_pos, self.coverage_neg = 0.368, 0.018

    @property
    def n_reached_pos(self):
        return round(self.coverage_pos * self.n_pos)

    evaluable = True
    discriminates = True


class _Cal:
    def __init__(self, virus, protocol="single-screen"):
        self.virus, self.protocol = virus, protocol
        self.channels = {"M1": _Ch(), "M4": _Ch()}
        self.refused_reason = ""
        self.usable = True


@pytest.fixture
def two_virus_graph(tmp_path):
    """The `graph` fixture plus a Zika taxon, so a virus exists with no
    calibration of its own."""
    def drug(nid, name):
        return {"id": nid, "class": "SmallMolecule",
                "properties": {"smiles": "C", "inchikey_skel": nid[-5:],
                               "chembl_id": name, "is_approved": True,
                               "label": name, "salt_collapsed": False,
                               "stereo_collapsed": False}}
    nodes = [
        {"id": V, "class": "OrganismTaxon",
         "properties": {"label": "SARS-CoV-2", "family": "Coronaviridae",
                        "baltimore_class": "IV", "is_enveloped": True}},
        {"id": ZIKA, "class": "OrganismTaxon",
         "properties": {"label": "Zika", "family": "Flaviviridae",
                        "baltimore_class": "IV", "is_enveloped": True}},
        {"id": "UniProtKB:NS3", "class": "Protein",
         "properties": {"taxon_id": ZIKA, "is_viral": True, "reviewed": True,
                        "sequence_hash": "h", "protein_family": "NS3"}},
        {"id": "KGAV:GZ", "class": "Gene",
         "properties": {"symbol": "POLY", "taxon_id": ZIKA, "is_viral": True}},
        drug(DRUG, "nirmatrelvir"),
    ]

    def e(s, p, o, q=None):
        # chembl.py sets target_resolution on EVERY INHIBITS edge, both
        # branches, and the schema requires it -- so M1 constrains on it and
        # qualifier_ok reads a missing qualifier as a failed match. A fixture
        # omitting it is not a simpler fixture; it is an edge the
        # direct-acting channel refuses to walk, in the scorer and now in the
        # agent too.
        if p == "INHIBITS":
            q = dict(q or {})
            q.setdefault("target_resolution", "chain_from_assay_description")
        return {"subject": s, "predicate": p, "object": o, "qualifiers": q or {},
                "primary_knowledge_source": "infores:chembl",
                "evidence_tier": 1, "first_asserted_date": "2021-03-04"}

    edges = [
        e(DRUG, "INHIBITS", "UniProtKB:NS3",
          {"assay_type": "biochemical", "ic50_nm": 25.0}),
        e("UniProtKB:NS3", "ENCODED_BY", "KGAV:GZ"),
        e("KGAV:GZ", "BELONGS_TO", ZIKA),
        e(DRUG, "HAS_ANTIVIRAL_ACTIVITY_AGAINST", ZIKA,
          {"assay_type": "cell_based_antiviral", "ec50_nm": 100.0}),
    ]
    (tmp_path / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text("\n".join(json.dumps(x) for x in edges))
    return tmp_path


def test_a_single_calibration_is_not_reused_for_another_virus(two_virus_graph):
    """THE BUG. GraphTools held one calibration and explain()/triage() both
    took a virus argument without using it, so an instance built with
    SARS-CoV-2's numbers answered Zika questions with SARS-CoV-2's AUC."""
    t = GraphTools(two_virus_graph, calibration=_Cal("SARS-CoV-2"))
    assert t.calibration_for(V) is not None, "its own virus still resolves"
    assert t.calibration_for(ZIKA) is None, \
        "another virus must resolve to None, never to these numbers"


def test_the_curie_and_the_label_name_the_same_virus(two_virus_graph):
    """calibration.load stores whatever key the CALLER passed -- kgav_ask
    passes "SARS-CoV-2", a caller holding the release passes the CURIE -- so
    neither form can be assumed. The release's own labels bridge them."""
    by_label = GraphTools(two_virus_graph, calibration=_Cal("SARS-CoV-2"))
    assert by_label.calibration_for(V) is not None
    by_curie = GraphTools(two_virus_graph, calibration=_Cal(V))
    assert by_curie.calibration_for("SARS-CoV-2") is not None
    assert by_curie.calibration_for(ZIKA) is None


def test_a_mapping_resolves_each_virus_to_its_own(two_virus_graph):
    t = GraphTools(two_virus_graph, calibration={
        "SARS-CoV-2": _Cal("SARS-CoV-2"), "Zika": _Cal("Zika")})
    assert t.calibration_for(V).virus == "SARS-CoV-2"
    assert t.calibration_for(ZIKA).virus == "Zika"


def test_a_mis_keyed_mapping_raises(two_virus_graph):
    """Keying one virus's calibration under another reintroduces the bug by a
    different route, so it is refused rather than tolerated."""
    with pytest.raises(ValueError, match="measured on"):
        GraphTools(two_virus_graph, calibration={"Zika": _Cal("SARS-CoV-2")})


def test_triage_on_an_uncalibrated_virus_quotes_no_auc(two_virus_graph):
    """END TO END, and this is the test that would have caught it. The channel
    above DISCRIMINATES, so a leak renders the module's most confident
    sentence -- "separates measured actives from measured inactives at AUC
    0.675" -- about a virus whose real M1 is 0.511, at chance."""
    t = GraphTools(two_virus_graph, calibration=_Cal("SARS-CoV-2"))
    out = t.triage([DRUG], ZIKA)
    text = (out.note or "") + " " + " ".join(out.warnings or [])
    assert "0.675" not in text, f"SARS-CoV-2's AUC leaked into a Zika answer:\n{text}"
    assert "separates measured actives" not in text
    assert "no measured performance" in text or "not been evaluated" in text, \
        f"expected the uncalibrated branch, got:\n{text}"


def test_triage_on_the_calibrated_virus_still_quotes_its_auc(two_virus_graph):
    """The fix must not silence the virus that WAS measured."""
    t = GraphTools(two_virus_graph, calibration=_Cal("SARS-CoV-2"))
    out = t.triage([DRUG], V)
    text = (out.note or "") + " " + " ".join(out.warnings or [])
    assert "0.675" in text, f"the measured virus lost its figure:\n{text}"


# --------------------------- a route contradicted by its own source says so
@pytest.fixture
def contradicted_release(tmp_path):
    """One drug, one nsp5 chain, and BOTH verdicts on the same pair.

    INHIBITS and MEASURED_INACTIVE_AGAINST come from the same ChEMBL rows
    with opposite verdicts, and assembly keys edges on subject, predicate
    and object -- so both survive as separate edges on one pair.
    """
    V = "NCBITaxon:2697049"
    nodes = [
        {"id": V, "class": "OrganismTaxon",
         "properties": {"label": "SARS-CoV-2", "family": "Coronaviridae",
                        "baltimore_class": "IV", "is_enveloped": True}},
        {"id": "UniProtKB:NSP5", "class": "Protein",
         "properties": {"taxon_id": V, "is_viral": True,
                        "protein_family": "nsp5", "sequence_hash": "h",
                        "reviewed": True}},
        {"id": "KGAV:GENE_rep", "class": "Gene",
         "properties": {"symbol": "rep", "taxon_id": V, "is_viral": True}},
        {"id": "INCHIKEY:TWOFACED", "class": "SmallMolecule",
         "properties": {"smiles": "C", "inchikey_skel": "TWOFAC",
                        "chembl_id": "CHEMBL9", "salt_collapsed": False,
                        "stereo_collapsed": False}},
    ]

    def e(s_, p_, o_, q=None):
        q = dict(q or {})
        if p_ == "INHIBITS":
            q.setdefault("target_resolution", "chain_from_assay_description")
        return {"subject": s_, "predicate": p_, "object": o_, "qualifiers": q,
                "primary_knowledge_source": "infores:chembl",
                "evidence_tier": 1, "first_asserted_date": "2021-03-04"}

    edges = [
        e("INCHIKEY:TWOFACED", "INHIBITS", "UniProtKB:NSP5",
          {"assay_type": "biochemical", "ic50_nm": 25.0}),
        e("INCHIKEY:TWOFACED", "MEASURED_INACTIVE_AGAINST", "UniProtKB:NSP5",
          {"assay_type": "biochemical", "ic50_nm": 50000.0, "relation": ">"}),
        e("UniProtKB:NSP5", "ENCODED_BY", "KGAV:GENE_rep"),
        e("KGAV:GENE_rep", "BELONGS_TO", V),
    ]
    (tmp_path / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text(
        "\n".join(json.dumps(x) for x in edges))
    return tmp_path


def test_a_contradicted_route_carries_the_caveat(contradicted_release):
    """The agent loads EVERY edge, held-out predicates included, so the
    contradicting measurement was already in hand. It reported "X inhibits
    nsp5" from a graph that also records "X was measured not to inhibit
    nsp5", and said nothing."""
    t = GraphTools(contradicted_release)
    r = t.explain("INCHIKEY:TWOFACED", "NCBITaxon:2697049")
    assert r.ok
    assert any("CONTRADICTED" in w for w in r.warnings)
    said = next(w for w in r.warnings if "CONTRADICTED" in w)
    assert "MEASURED INACTIVE" in said
    assert "50000" in said or "50,000" in said or "> 50000" in said


def test_the_caveat_lands_in_the_route_block_that_triggered_it(
        contradicted_release):
    """Per route, like the cell-context caveat: a contradiction naming one
    protein printed against a route through another discredits every warning
    on the page."""
    t = GraphTools(contradicted_release)
    r = t.explain("INCHIKEY:TWOFACED", "NCBITaxon:2697049")
    text = "\n".join(r.verbalised)
    assert "CONTRADICTED" in text


def test_an_uncontradicted_route_carries_no_such_caveat(two_virus_graph):
    """So the caveat marks a contradiction rather than appearing on every
    direct-acting route."""
    t = GraphTools(two_virus_graph)
    for virus in ("NCBITaxon:2697049", "NCBITaxon:64320"):
        r = t.explain(DRUG, virus)
        if r.ok:
            assert not any("CONTRADICTED" in w for w in r.warnings)


def test_the_contradiction_is_enforced_like_every_other_caveat(
        contradicted_release):
    """It is in CAVEAT_LABELS, so verify.dropped_warnings requires its label
    in the answer. A caveat the model may drop for brevity is not a caveat."""
    from kgav.agent.verify import caveat_labels, dropped_warnings
    t = GraphTools(contradicted_release)
    r = t.explain("INCHIKEY:TWOFACED", "NCBITaxon:2697049")
    assert "CONTRADICTED" in caveat_labels(r.warnings)
    assert any(getattr(w, "key", "") == "contradicted"
               for w in dropped_warnings("an answer mentioning nothing",
                                         r.warnings))
