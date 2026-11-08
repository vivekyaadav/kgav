"""Gate for the tool layer: bounded output, provenance in every sentence, and
absence reported as absence rather than as a negative finding.
"""
import json
from pathlib import Path

import pytest

from kgav.assemble import fact_id
from kgav.agent.tools import GraphTools
from kgav.schema import load_schema

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

    def discriminating(self):
        """Matches Calibration's interface.

        This double had no such method: the only caller sat behind a
        `virus != SUPPORTED_VIRUS` gate and SARS-CoV-2 is the virus
        these tests calibrate, so the branch was unreachable and the
        gap invisible.
        """
        return [n for n, c in sorted(self.channels.items())
                if getattr(c, "discriminates", False)]


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
        # SARS-CoV-2 NEEDS A ROUTE. The fixture gave it a calibration and
        # no inhibitor, host factor or virus-host interaction, so once the
        # tools began acting on the no-route verdict (2026-11-02, audit
        # A6) triage refused it and the test read as a regression. A virus
        # the agent will answer about has at least one route; a virus with
        # none is the Yellow fever tier and is refused by design.
        {"id": "UniProtKB:NSP5", "class": "Protein",
         "properties": {"taxon_id": V, "is_viral": True, "reviewed": True,
                        "sequence_hash": "h", "protein_family": "NSP5"}},
        {"id": "UniProtKB:SIGMAR1", "class": "Protein",
         "properties": {"taxon_id": "NCBITaxon:9606", "is_viral": False,
                        "reviewed": True, "sequence_hash": "h"}},
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
        e("UniProtKB:NSP5", "PHYSICALLY_INTERACTS_WITH",
          "UniProtKB:SIGMAR1"),
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


# ------------------------------------------- the answer sheet, marked
# A11. VirusRoutes' docstring says of the label count: "never to answer
# from: quoting a label as evidence is the leak every other guard here
# exists to prevent." Three tools did exactly that, silently --
# evidence(), profile() and NEIGHBOURS(), which the external audit did
# not name. evidence() was worst: for a label edge the only warning it
# attached was evidence_note's "EVIDENCE: Curated experimental result",
# so the answer sheet arrived carrying the tier-1 reassurance.
#
# The edges are MARKED, not suppressed. For a virus with no route the
# agent's own closing offer is to report what the graph stores and
# retrieve evidence for a named claim, and for Yellow fever the graph
# stores little else but these labels -- suppressing them would empty
# the one offer left.
HELD_OUT = sorted(load_schema().held_out_predicates())


def test_the_schema_declares_more_than_one_held_out_predicate():
    """The audit named HAS_ANTIVIRAL_ACTIVITY_AGAINST. The schema holds
    out four, and a fix covering one would have left three leaking."""
    assert len(HELD_OUT) >= 4, HELD_OUT
    assert "HAS_ANTIVIRAL_ACTIVITY_AGAINST" in HELD_OUT
    assert "MEASURED_INACTIVE_AGAINST" in HELD_OUT


@pytest.mark.parametrize("predicate", HELD_OUT)
def test_every_held_out_predicate_is_marked_when_rendered(tmp_path, predicate):
    """Parametrised over the SCHEMA, not over a list kept here. A
    predicate newly declared held_out is covered by this test the moment
    it is declared."""
    nodes = [{"id": "CHEMBL:X", "class": "SmallMolecule",
              "properties": {"label": "testdrug"}},
             {"id": "NCBITaxon:1", "class": "OrganismTaxon",
              "properties": {"label": "Testvirus", "family": "Flaviviridae",
                             "baltimore_class": "IV", "is_enveloped": True}},
             {"id": "MONDO:1", "class": "Disease",
              "properties": {"label": "testdisease"}}]
    obj = "MONDO:1" if predicate in ("TREATS", "IN_TRIAL_FOR") \
        else "NCBITaxon:1"
    edges = [{"subject": "CHEMBL:X", "predicate": predicate, "object": obj,
              "qualifiers": {}, "primary_knowledge_source": "infores:chembl",
              "evidence_tier": 1, "first_asserted_date": "2020-01-01"}]
    (tmp_path / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text(
        "\n".join(json.dumps(x) for x in edges))
    t = GraphTools(tmp_path)
    r = t.evidence("CHEMBL:X")
    assert r.ok
    assert all(ln.startswith("[HELD-OUT LABEL]") for ln in r.verbalised), \
        r.verbalised
    assert "ANSWER SHEET" in [w.label for w in r.warnings]


def test_a_feature_edge_is_not_marked(graph):
    """The marker must distinguish, not decorate everything."""
    r = graph.evidence(DRUG, "UniProtKB:NSP5")
    assert r.ok
    assert not any("HELD-OUT LABEL" in ln for ln in r.verbalised), r.verbalised
    assert "ANSWER SHEET" not in [getattr(w, "label", "") for w in r.warnings]


def test_a_mixed_result_marks_only_the_label_lines(graph):
    r = graph.evidence(DRUG)
    assert r.ok
    marked = [ln for ln in r.verbalised if ln.startswith("[HELD-OUT LABEL]")]
    plain = [ln for ln in r.verbalised if not ln.startswith("[HELD-OUT LABEL]")]
    assert marked and plain, r.verbalised
    assert all("shows activity against" in ln for ln in marked), marked
    assert all("inhibits the viral protein" in ln for ln in plain), plain


def test_the_answer_sheet_caveat_comes_before_the_tier_reassurance(graph):
    """ORDER MATTERS HERE. The tier note says "Curated experimental
    result", which raises confidence; it was the ONLY warning a label
    edge carried."""
    labels = [getattr(w, "label", "") for w in graph.evidence(DRUG).warnings]
    assert "ANSWER SHEET" in labels and "EVIDENCE" in labels, labels
    assert labels.index("ANSWER SHEET") < labels.index("EVIDENCE"), labels


def test_neighbours_marks_them_too(graph):
    """The tool the audit did not name. It had the same hole and NO
    warning of any kind -- found by probing every public tool against
    one fixture instead of the two the audit listed."""
    r = graph.neighbours(DRUG)
    assert r.ok
    assert any(ln.startswith("[HELD-OUT LABEL]") for ln in r.verbalised)
    assert "ANSWER SHEET" in [getattr(w, "label", "") for w in r.warnings]


def test_profile_lists_labels_apart_from_features(graph):
    """One flat list of predicate counts is how a Yellow fever profile
    read "1,030 antiviral measurements" -- which looks like 1,030
    supporting facts and is 1,030 rows of the answer sheet."""
    r = graph.profile(V)
    assert r.ok
    lines = r.verbalised
    heading = next(i for i, ln in enumerate(lines)
                   if "held-out evaluation labels" in ln)
    # EXACTLY ONCE, AND BELOW THE HEADING. The first version of this
    # test asserted only that the heading and the data split existed,
    # and a mutant that left the labels in the FEATURE rows as well --
    # listing each label predicate twice, once as a feature -- survived
    # it. The verbalised block is what the model is shown, so a label
    # counted among the features is the leak back in full.
    for pred in HELD_OUT:
        hits = [i for i, ln in enumerate(lines) if ln.strip().startswith(pred)]
        assert len(hits) <= 1, f"{pred} listed {len(hits)} times: {lines}"
        for i in hits:
            assert i > heading, f"{pred} listed among the features: {lines}"
    assert any(ln.strip().startswith("HAS_ANTIVIRAL_ACTIVITY_AGAINST")
               for ln in lines[heading + 1:]), lines
    assert "HAS_ANTIVIRAL_ACTIVITY_AGAINST" in r.data["held_out_label_counts"]
    assert "HAS_ANTIVIRAL_ACTIVITY_AGAINST" not in r.data[
        "feature_edge_counts"]
    said = next(w for w in r.warnings
                if getattr(w, "label", "") == "ANSWER SHEET")
    # The COUNT and the predicate are facts and are pinned. The rest of
    # the prose is deliberately not: guards.py's own contract is that
    # "the prose below it can be rewritten freely without changing what
    # is checked", and the label is what verification enforces.
    assert "HAS_ANTIVIRAL_ACTIVITY_AGAINST" in said
    assert str(r.data["held_out_label_counts"][
        "HAS_ANTIVIRAL_ACTIVITY_AGAINST"]) in said


def test_profile_keeps_edge_counts_as_the_full_census(graph):
    """The split is ADDITIONAL. edge_counts is read elsewhere as the
    node's degree, and redefining a key in place is how a consumer
    starts under-counting without failing."""
    r = graph.profile(V)
    total = sum(r.data["edge_counts"].values())
    assert total == sum(r.data["feature_edge_counts"].values()) \
        + sum(r.data["held_out_label_counts"].values())
    assert total > 0


def test_a_node_with_no_label_gets_no_heading_and_no_caveat(graph):
    """Said only when true. A heading reading "held-out evaluation
    labels:" above nothing teaches a reader to skip it."""
    r = graph.profile("UniProtKB:NSP5")
    assert r.ok
    assert "held-out evaluation labels" not in " ".join(r.verbalised)
    assert not r.data["held_out_label_counts"]
    assert "ANSWER SHEET" not in [getattr(w, "label", "") for w in r.warnings]


def test_copying_a_marked_fact_does_not_satisfy_the_caveat(graph):
    """THE REASON THE TWO MARKERS DIFFER. verify.dropped_warnings
    enforces a caveat by searching the answer for its LABEL. The
    per-line marker is "[HELD-OUT LABEL]" and the caveat's label is
    "ANSWER SHEET": had they matched, a model that merely copied a
    marked fact line would satisfy the caveat without ever stating what
    it means, and the enforcement could never bite."""
    from kgav.agent import verify
    r = graph.evidence(DRUG)
    copied = next(ln for ln in r.verbalised
                  if ln.startswith("[HELD-OUT LABEL]"))
    dropped = [getattr(w, "label", "") for w in
               verify.dropped_warnings(copied, r.warnings)]
    assert "ANSWER SHEET" in dropped, \
        "copying a marked fact satisfied the caveat that explains it"
    kept = verify.dropped_warnings(
        copied + " ANSWER SHEET: these are held-out evaluation labels.",
        r.warnings)
    assert "ANSWER SHEET" not in [getattr(w, "label", "") for w in kept]


# ------------------------------------- [NAMES THIS VIRUS], by longest match
# A12. Measured against release v0.1 before the fix: 6 Reactome pathways
# naming SARS-CoV-2 were flagged as naming SARS-CoV, and 6 dengue
# pathways were flagged for nothing. A circularity warning on the WRONG
# virus is worse than none -- it tells a reader the overlap for SARS-CoV
# was decided in advance when nothing about SARS-CoV was involved.
@pytest.fixture
def namer(tmp_path):
    """A release holding four taxa and nothing else, for label matching."""
    def taxon(nid, label, family):
        return {"id": nid, "class": "OrganismTaxon",
                "properties": {"label": label, "family": family,
                               "baltimore_class": "IV", "is_enveloped": True}}
    nodes = [taxon("NCBITaxon:2697049", "SARS-CoV-2", "Coronaviridae"),
             taxon("NCBITaxon:694009", "SARS-CoV", "Coronaviridae"),
             taxon("NCBITaxon:11053", "DENV-1", "Flaviviridae"),
             taxon("NCBITaxon:11060", "DENV-2", "Flaviviridae")]
    (tmp_path / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text(
        '{"subject": "a", "predicate": "P", "object": "b"}')
    return GraphTools(tmp_path)


SARS2, SARS1 = "NCBITaxon:2697049", "NCBITaxon:694009"
DEN1, DEN2 = "NCBITaxon:11053", "NCBITaxon:11060"

#: Every label here is a REAL Reactome pathway label from release v0.1,
#: except the two negatives and the two-virus case. The six SARS-CoV-2
#: rows are exactly the six the probe found mis-flagged.
NAMING_CASES = [
    ("Replication of the SARS-CoV-2 genome", {SARS2: True, SARS1: False}),
    ("SARS-CoV-2 modulates autophagy", {SARS2: True, SARS1: False}),
    ("SARS-CoV-2 modulates host translation machinery",
     {SARS2: True, SARS1: False}),
    ("SARS-CoV-2 targets PDZ proteins in cell-cell junction",
     {SARS2: True, SARS1: False}),
    ("SARS-CoV-2 targets host intracellular signalling and regulation",
     {SARS2: True, SARS1: False}),
    ("SARS-CoV-2 activates/modulates innate and adaptive immune responses",
     {SARS2: True, SARS1: False}),
    ("Dengue Virus-Host Interactions", {DEN1: True, DEN2: True, SARS2: False}),
    ("Dengue Virus Attachment and Entry", {DEN1: True, DEN2: True}),
    ("Assembly and Release of Dengue Virus Virions",
     {DEN1: True, DEN2: True}),
    ("Uncoating of the HIV Virion", {SARS2: False, DEN1: False}),
    ("Glycolysis", {SARS2: False, DEN1: False}),
]


@pytest.mark.parametrize("label,want", NAMING_CASES,
                         ids=[c[0][:34] for c in NAMING_CASES])
def test_a_pathway_label_names_only_the_virus_it_names(namer, label, want):
    got = {v: namer._names_virus(label, v) for v in want}
    assert got == want, f"{label!r}: {got} != {want}"


def test_a_label_naming_two_viruses_names_both(namer):
    """SHADOWING IS BY SPAN, NOT GLOBAL. The eight-character "sars-cov"
    at offset 0 is not shadowed by the ten-character "sars-cov-2" at
    offset 13. Suppressing a name whenever a longer one matched ANYWHERE
    in the label would drop the first virus -- which is the opposite
    error from the one being fixed, and just as wrong."""
    label = "SARS-CoV and SARS-CoV-2 entry"
    assert namer._names_virus(label, SARS1)
    assert namer._names_virus(label, SARS2)


def test_a_name_is_found_at_every_occurrence_not_only_the_first(namer):
    """THE DISCRIMINATING ORDER, and the one my first set of cases could
    not produce. In "SARS-CoV and SARS-CoV-2 entry" the eight-character
    name happens to come first and is unshadowed there, so a matcher
    recording only a name's FIRST occurrence passes that case. Reverse
    the order and the first "sars-cov" sits at offset 0 inside
    "sars-cov-2" where it IS shadowed, while the real mention at offset
    16 is never looked at -- so the label names SARS-CoV plainly and the
    flag is dropped. A surviving mutant found this; no test did."""
    label = "SARS-CoV-2 and SARS-CoV entry"
    assert namer._names_virus(label, SARS2)
    assert namer._names_virus(label, SARS1), \
        "only the first occurrence of the name was searched"


def test_the_serotype_specific_name_wins_over_the_shared_one(namer):
    """The two mechanisms compose. "dengue virus" names all four
    serotypes and "dengue virus 1" names one, so a label carrying the
    longer name must flag DENV-1 alone."""
    assert namer._names_virus("Dengue Virus 1 Attachment", DEN1)
    assert not namer._names_virus("Dengue Virus 1 Attachment", DEN2)


def test_the_threshold_is_the_registers_and_not_a_second_copy():
    """It was a literal `len(n) >= 4` in tools.py. One owner."""
    import inspect
    from kgav.agent import tools as T
    src = inspect.getsource(T.GraphTools._label_name_map)
    assert "MIN_SEARCH_NAME" in src
    assert "len(lab) >= 4" not in src


def test_a_taxon_absent_from_the_release_is_not_matched(namer):
    """Register names are only matched for taxa this graph holds, the
    same rule resolution already applies: flagging a pathway for a virus
    with no node behind it reports a confound that cannot exist here."""
    assert "NCBITaxon:11089" not in namer.nodes       # Yellow fever, unbuilt
    assert not namer._names_virus("Yellow fever virus entry",
                                  "NCBITaxon:11089")


# --------------------------- a contested direction is not a stated one
# S9. assemble.resolve_qualifier now writes "contested" where two curated
# sources assert opposite directions, instead of "unknown". The two mean
# different things -- nobody said, against the sources disagree -- and a
# consumer that cannot tell them apart reports a gap in the data where
# there is a contradiction in it.
def _targets_graph(tmp_path, direction):
    nodes = [{"id": "CHEMBL:D", "class": "SmallMolecule",
              "properties": {"label": "testdrug"}},
             {"id": "UniProtKB:SIGMAR1", "class": "Protein",
              "properties": {"taxon_id": "NCBITaxon:9606", "is_viral": False,
                             "sequence_hash": "h", "reviewed": True,
                             "gene_symbol": "SIGMAR1"}}]
    q = {} if direction is None else {"direction": direction}
    edges = [{"subject": "CHEMBL:D", "predicate": "TARGETS",
              "object": "UniProtKB:SIGMAR1", "qualifiers": q,
              "primary_knowledge_source": "infores:chembl",
              "evidence_tier": 1, "first_asserted_date": "2021-01-01"}]
    (tmp_path / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text(
        "\n".join(json.dumps(e) for e in edges))
    return GraphTools(tmp_path)


def test_a_contested_direction_is_not_rendered_as_a_mechanism(tmp_path):
    """It read `as an {direction}` for anything that was not None or
    "unknown", which produced "engages the host protein SIGMAR1 as an
    contested" -- a sentence asserting a mechanism named contested."""
    t = _targets_graph(tmp_path, "contested")
    line = t.neighbours("CHEMBL:D").verbalised[0]
    assert "as an contested" not in line, line
    assert "sources disagree on whether it inhibits or activates" in line


@pytest.mark.parametrize("direction,shows", [
    ("inhibitor", "as an inhibitor"),
    ("agonist", "as an agonist"),
    ("unknown", "mechanism not stated"),
    (None, "mechanism not stated"),
])
def test_the_other_directions_render_as_before(tmp_path, direction, shows):
    t = _targets_graph(tmp_path, direction)
    assert shows in t.neighbours("CHEMBL:D").verbalised[0]


@pytest.mark.parametrize("direction,counts", [
    ("inhibitor", True), ("agonist", True), ("modulator", True),
    ("unknown", False), ("contested", False), (None, False),
])
def test_only_a_stated_direction_counts_as_directed(tmp_path, direction,
                                                    counts):
    """The count behind TARGET DIRECTION answers "does this edge say what
    the drug DOES to the target". Two sources asserting opposite actions
    do not, and "unknown" was excluded only by being falsy -- luck rather
    than a rule, which is why "contested" slipped through it."""
    import inspect

    from kgav.agent import tools as T
    src = inspect.getsource(T.GraphTools.pathways)
    assert 'not in (None, "", "unknown", "contested")' in src
    t = _targets_graph(tmp_path, direction)
    edges = [e for e in t.edges if e["predicate"] == "TARGETS"]
    n = sum(1 for e in edges
            if (e.get("qualifiers") or {}).get("direction")
            not in (None, "", "unknown", "contested"))
    assert bool(n) is counts


# ------------------------- A3/A4: what triage ranks on, and what it calls it
#
# A4. `_best_direct_potency` took the minimum over EVERY INHIBITS edge the
# compound has, whatever virus the protein belonged to, across ic50_nm,
# ki_nm and kd_nm together -- and the caller printed the result as "best
# direct IC50" beside a named virus. Measured on v0.1 of 2026-11-09 over
# 3,823 (compound, virus) pairs: 133 reported a potency measured against a
# DIFFERENT virus, 4 reported one for a virus whose own edges hold no
# potency at all, and 439 reported a Ki or a Kd under the name IC50.
#
# A3. triage orders compounds by the selectivity index, which is a
# qualifier of the selectivity layer's edges -- every one of them a
# HAS_ANTIVIRAL_ACTIVITY_AGAINST edge, the evaluation's own labels.
# `_held_out_warnings` covers a tool that RETURNS such an edge; nothing
# covered a tool that merely SORTS on one.

V2 = "NCBITaxon:11137"          # HCoV-229E, a second virus in the register


def _two_virus_graph(tmp_path, *, with_selectivity: bool,
                     ic50_on_v2: float | None = 2560.0,
                     ic50_on_v1: float | None = None):
    """One compound inhibiting a protein of each of two viruses.

    The numbers are the shape of the worst real case:
    FZNWILPXXNKETO-YADHBBJMSA-N was reported at 18 nM against HCoV-229E
    from a Kd measured on a SARS-CoV-2 protein, where its own scoped
    functional IC50 is 2,560 nM.
    """
    release = tmp_path / "rel"
    release.mkdir()

    def prot(pid, taxon, fam):
        return {"id": pid, "class": "Protein",
                "properties": {"taxon_id": taxon, "is_viral": True,
                               "sequence_hash": "h", "reviewed": True,
                               "protein_family": fam}}

    def gene(gid, taxon):
        return {"id": gid, "class": "Gene",
                "properties": {"symbol": "rep", "taxon_id": taxon,
                               "is_viral": True}}

    def taxon(tid, label, fam):
        return {"id": tid, "class": "OrganismTaxon",
                "properties": {"label": label, "family": fam,
                               "baltimore_class": "IV", "is_enveloped": True}}
    nodes = [
        taxon(V, "SARS-CoV-2", "Coronaviridae"),
        taxon(V2, "HCoV-229E", "Coronaviridae"),
        prot("UniProtKB:P1", V, "nsp5"), gene("KGAV:G1", V),
        prot("UniProtKB:P2", V2, "nsp5"), gene("KGAV:G2", V2),
        {"id": DRUG, "class": "SmallMolecule",
         "properties": {"smiles": "C", "inchikey_skel": "NIRMA",
                        "chembl_id": "thecompound", "is_approved": True,
                        "label": "thecompound", "salt_collapsed": False,
                        "stereo_collapsed": False}},
    ]

    def e(s, p, o, q=None):
        if p == "INHIBITS":
            q = dict(q or {})
            q.setdefault("target_resolution", "chain_from_assay_description")
        return {"subject": s, "predicate": p, "object": o,
                "qualifiers": q or {},
                "primary_knowledge_source": "infores:chembl",
                "evidence_tier": 1, "first_asserted_date": "2021-03-04"}
    edges = [
        # an 18 nM Kd against virus 1 -- the flattering number
        e(DRUG, "INHIBITS", "UniProtKB:P1",
          {"assay_type": "biochemical", "kd_nm": 18.0}),
        e("UniProtKB:P1", "ENCODED_BY", "KGAV:G1"),
        e("KGAV:G1", "BELONGS_TO", V),
        e("UniProtKB:P2", "ENCODED_BY", "KGAV:G2"),
        e("KGAV:G2", "BELONGS_TO", V2),
    ]
    # virus 2's own edge: present with an IC50, present without a potency
    # at all, or absent entirely -- the three cases the measurements found.
    if ic50_on_v2 is not None:
        edges.insert(1, e(DRUG, "INHIBITS", "UniProtKB:P2",
                          {"assay_type": "biochemical",
                           "ic50_nm": ic50_on_v2}))
    else:
        edges.insert(1, e(DRUG, "INHIBITS", "UniProtKB:P2",
                          {"assay_type": "biochemical"}))
    # a second readout on virus 1's OWN protein, to test the preference
    # order where both quantities describe the same virus.
    if ic50_on_v1 is not None:
        edges.append(e(DRUG, "INHIBITS", "UniProtKB:P1",
                       {"assay_type": "biochemical", "ic50_nm": ic50_on_v1}))
    (release / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in nodes))
    (release / "edges.jsonl").write_text(
        "\n".join(json.dumps(x) for x in edges))
    if with_selectivity:
        layer = tmp_path / "rel-selectivity"
        layer.mkdir()
        layer.joinpath("edges.jsonl").write_text(json.dumps(
            e(DRUG, "HAS_ANTIVIRAL_ACTIVITY_AGAINST", V2,
              {"assay_type": "cell_based_antiviral", "ec50_nm": 100.0,
               "selectivity_index": 50.0, "selectivity_verified": True})))
    return GraphTools(release)


def test_the_potency_is_the_one_measured_against_this_virus(tmp_path):
    """The compound has an 18 nM Kd on a SARS-CoV-2 protein and a 2,560 nM
    IC50 on an HCoV-229E protein. Asked about HCoV-229E it must report
    2,560, and asked about SARS-CoV-2 it must report 18 -- not 18 for
    both, which is what a minimum over every edge gives."""
    g = _two_virus_graph(tmp_path, with_selectivity=False)
    assert g._best_direct_potency(DRUG, V2) == (2560.0, "ic50_nm")
    assert g._best_direct_potency(DRUG, V) == (18.0, "kd_nm")


def test_a_binding_constant_is_not_printed_as_an_ic50(tmp_path):
    """"best direct IC50 18 nM" for a Kd is a different measurement under a
    borrowed name, and 18 against 2,560 is the size of the difference."""
    g = _two_virus_graph(tmp_path, with_selectivity=False)
    line = g.triage([DRUG], V2).verbalised[0]
    assert "best direct IC50 2560 nM" in line, line
    other = g.triage([DRUG], V).verbalised[0]
    assert "best direct Kd 18 nM" in other, other
    assert "IC50" not in other, "a Kd is still being called an IC50"


def test_an_ic50_is_preferred_over_a_lower_binding_constant(tmp_path):
    """A PREFERENCE ORDER, NOT A MINIMUM -- labels.POTENCY_FIELDS' own
    convention, where `_potency` reads the first field present rather than
    the smallest number. Both readouts describe virus 1 here: the 18 nM Kd
    and a 9,000 nM IC50. The IC50 wins although it is five hundred times
    larger, because they are not the same quantity and the smaller one is
    not the better measurement."""
    g = _two_virus_graph(tmp_path, with_selectivity=False, ic50_on_v1=9000.0)
    assert g._best_direct_potency(DRUG, V) == (9000.0, "ic50_nm")
    assert "best direct IC50 9000 nM" in g.triage([DRUG], V).verbalised[0]


def test_the_potency_has_no_value_for_a_virus_with_no_measurement(tmp_path):
    """4 real (compound, virus) pairs had no potency on any edge of their
    own and were given one from another virus. Absence must read as
    absence."""
    g = _two_virus_graph(tmp_path, with_selectivity=False, ic50_on_v2=None)
    # The route to virus 2 exists -- so triage runs rather than refusing --
    # and its edge carries no potency field at all.
    assert g._best_direct_potency(DRUG, V2) is None
    assert "no direct-acting data" in g.triage([DRUG], V2).verbalised[0]


def test_triage_says_the_order_came_from_the_answer_sheet(tmp_path):
    """A3. The selectivity index that decides the order is a qualifier of a
    held-out label edge, and triage printed no caveat about it."""
    from kgav.agent.verify import caveat_labels
    g = _two_virus_graph(tmp_path, with_selectivity=True)
    assert g.si_predicate_counts == {"HAS_ANTIVIRAL_ACTIVITY_AGAINST": 1}
    r = g.triage([DRUG], V2)
    assert "ANSWER SHEET" in caveat_labels(r.warnings), (
        f"triage ranked on a held-out label's qualifier without saying so: "
        f"{caveat_labels(r.warnings)}")
    text = " ".join(str(w) for w in r.warnings)
    assert "HAS_ANTIVIRAL_ACTIVITY_AGAINST" in text
    assert "answer sheet" in text


def test_no_answer_sheet_caveat_when_the_order_used_no_selectivity(tmp_path):
    """Not unconditional decoration: with no selectivity layer there is no
    label-derived quantity in the sort, and the caveat must be absent. A
    caveat attached to everything is a caveat a reader stops reading."""
    from kgav.agent.verify import caveat_labels
    g = _two_virus_graph(tmp_path, with_selectivity=False)
    r = g.triage([DRUG], V2)
    assert "ANSWER SHEET" not in caveat_labels(r.warnings)


def test_the_caveat_follows_the_schema_and_not_a_hardcoded_predicate(tmp_path):
    """Derived, so a release whose selectivity index sits on a predicate
    the schema does NOT hold out produces no claim about labels."""
    from kgav.agent.verify import caveat_labels
    g = _two_virus_graph(tmp_path, with_selectivity=True)
    assert set(g.si_predicate_counts) <= set(g._held_out), (
        "this fixture's selectivity index is not on a held-out predicate, "
        "so the test below would pass vacuously")
    g._held_out = frozenset()
    assert "ANSWER SHEET" not in caveat_labels(g.triage([DRUG], V2).warnings)


def test_the_release_contains_the_cases_a4_is_about():
    """The fixtures above are miniatures; this is the measurement.

    Every figure is parsed out of `_best_direct_potency`'s own docstring
    and recomputed from the release, for the reason test_landscape.py
    states: a check written with the same literals it is checking cannot
    survive the event it exists for. A rebuild that moves these counts
    fails here and names the docstring.
    """
    import re
    from collections import defaultdict
    from kgav.agent.tools import GraphTools, POTENCY_PREFERENCE

    release = Path(__file__).resolve().parents[1] / "data/releases/v0.1"
    if not (release / "edges.jsonl").exists():
        pytest.skip("release v0.1 is not built in this checkout")

    inhibits = defaultdict(list)
    encoded, belongs = defaultdict(list), defaultdict(list)
    for line in (release / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        edge = json.loads(line)
        pred = edge["predicate"]
        if pred == "INHIBITS":
            inhibits[edge["subject"]].append(
                (edge["object"], edge.get("qualifiers") or {}))
        elif pred == "ENCODED_BY":
            encoded[edge["subject"]].append(edge["object"])
        elif pred == "BELONGS_TO":
            belongs[edge["subject"]].append(edge["object"])
    virus_of = {p: {t for g in gs for t in belongs.get(g, [])}
                for p, gs in encoded.items()}

    def old_best(rows):
        """What the code did: one min over every field and every virus."""
        vals = [(q[k], k) for _p, q in rows for k in POTENCY_PREFERENCE
                if q.get(k) is not None]
        return min(vals) if vals else None

    pairs = cross_virus = no_own_potency = not_an_ic50 = changed = 0
    genuinely_not_an_ic50 = 0
    for _compound, rows in inhibits.items():
        reported = old_best(rows)
        if reported is None:
            continue
        for virus in {v for p, _q in rows for v in virus_of.get(p, ())}:
            scoped = [(p, q) for p, q in rows if virus in virus_of.get(p, ())]
            pairs += 1
            best = None
            for name in POTENCY_PREFERENCE:
                vals = [q[name] for _p, q in scoped if q.get(name) is not None]
                if vals:
                    best = (min(vals), name)
                    break
            if best is None:
                no_own_potency += 1
                continue
            # TWO DISTINCT QUANTITIES, and conflating them is how this test
            # first disagreed with the docstring by 125. "Measured against a
            # different virus" is the OLD rule applied to the scoped rows
            # against the old rule applied to all of them -- the scope alone.
            # "Changed" is the new value against the old one, which also
            # moves wherever the preference order overrides a smaller
            # binding constant.
            scoped_old = old_best(scoped)
            if scoped_old[0] != reported[0]:
                cross_virus += 1
            if best[0] != reported[0]:
                changed += 1
            # Again two quantities. Under the OLD rule a binding constant
            # won on value and was printed as an IC50, which is the defect's
            # size; under the new one only the pairs with no IC50 anywhere
            # report a Ki or Kd, and those are named.
            if scoped_old[1] != "ic50_nm":
                not_an_ic50 += 1
            if best[1] != "ic50_nm":
                genuinely_not_an_ic50 += 1

    doc = " ".join(GraphTools._best_direct_potency.__doc__.split())
    m = re.search(r"over ([\d,]+) \(compound, virus\) pairs: \*\*([\d,]+) "
                  r"report a potency measured against a different virus\*\*, "
                  r"and ([\d,]+) report one", doc)
    assert m, "the docstring no longer states the pair census as read here"
    stated = tuple(int(g.replace(",", "")) for g in m.groups())
    assert (pairs, cross_virus, no_own_potency) == stated, (
        f"the docstring states {stated} (pairs, cross-virus, no own potency) "
        f"and the release gives {(pairs, cross_virus, no_own_potency)}")
    m2 = re.search(r"\*\*([\d,]+) of the ([\d,]+) pairs have a best value "
                   r"that is not an IC50", doc)
    assert m2, "the docstring no longer states the non-IC50 count"
    assert (not_an_ic50, pairs) == tuple(
        int(g.replace(",", "")) for g in m2.groups()), (
        f"the docstring states {m2.groups()} and the release gives "
        f"{(not_an_ic50, pairs)}")
    m3 = re.search(r"([\d,]+) of the pairs change the number they report", doc)
    assert m3, "the docstring no longer states how many figures move"
    assert changed == int(m3.group(1).replace(",", "")), (
        f"the docstring says {m3.group(1)} pairs change and the release "
        f"gives {changed:,}")
    m4 = re.search(r"([\d,]+) pairs have no IC50 at all", doc)
    assert m4, "the docstring no longer states how many are genuinely not IC50s"
    assert genuinely_not_an_ic50 == int(m4.group(1).replace(",", "")), (
        f"the docstring says {m4.group(1)} pairs have no IC50 and the "
        f"release gives {genuinely_not_an_ic50:,}")
    assert cross_virus and not_an_ic50, (
        "neither defect is reachable on this release any more; the "
        "measurement in that docstring is no longer about anything")
