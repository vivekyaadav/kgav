"""Gate for the baselines: held-out edges are never traversed, DWPC
down-weights promiscuous intermediates, proximity is degree-corrected.
"""
import json
from pathlib import Path

import pytest

from kgav.baselines import (
    Graph,
    combine,
    degree_ranking,
    dwpc_scores,
    hits_at_k,
    label_publication_index,
    mrr,
    network_proximity,
    rank,
    walk_metapath,
)
from kgav.schema import load_schema

V = "NCBITaxon:2697049"
D1 = "INCHIKEY:D1"
D2 = "INCHIKEY:D2"
D3 = "INCHIKEY:D3"
D4 = "INCHIKEY:D4"


def _prot(pid, viral, sym=None, fam=None):
    p = {"taxon_id": V if viral else "NCBITaxon:9606", "is_viral": viral,
         "sequence_hash": "h", "reviewed": True}
    if sym:
        p["gene_symbol"] = sym
    if fam:
        p["protein_family"] = fam
    return {"id": pid, "class": "Protein", "properties": p}


def _drug(nid):
    return {"id": nid, "class": "SmallMolecule",
            "properties": {"smiles": "C", "inchikey_skel": nid[-2:], "is_approved": True,
                           "chembl_id": "CHEMBL_" + nid[-2:],
                           "salt_collapsed": False, "stereo_collapsed": False}}


def _e(s, p, o, q=None, date="2020-01-01"):
    return {"subject": s, "predicate": p, "object": o, "qualifiers": q or {},
            "primary_knowledge_source": "s", "evidence_tier": 1,
            "first_asserted_date": date}


@pytest.fixture
def release(tmp_path):
    nodes = [
        {"id": V, "class": "OrganismTaxon",
         "properties": {"label": "SARS-CoV-2", "family": "Coronaviridae",
                        "baltimore_class": "IV", "is_enveloped": True}},
        _prot("UniProtKB:NSP5", True, fam="nsp5"),
        {"id": "KGAV:GENE_rep", "class": "Gene",
         "properties": {"symbol": "rep", "taxon_id": V, "is_viral": True}},
        _prot("UniProtKB:ACE2", False, sym="ACE2"),
        _prot("UniProtKB:HUB", False, sym="HUB"),
        _prot("UniProtKB:DEP2", False, sym="DEP2"),
        _drug(D1), _drug(D2), _drug(D3),
    ]
    edges = [
        _e(D1, "INHIBITS", "UniProtKB:NSP5",
           {"assay_type": "biochemical", "ic50_nm": 25.0}),
        _e("UniProtKB:NSP5", "ENCODED_BY", "KGAV:GENE_rep"),
        _e("KGAV:GENE_rep", "BELONGS_TO", V),
        _e(D2, "TARGETS", "UniProtKB:ACE2",
           {"direction": "inhibitor", "assay_type": "binding", "ic50_nm": 50.0}),
        _e("UniProtKB:ACE2", "HOST_FACTOR_FOR", V,
           {"direction": "dependency", "screen_type": "CRISPRko", "cell_line": "x"}),
        _e(D3, "TARGETS", "UniProtKB:HUB",
           {"direction": "unknown", "assay_type": "binding", "ic50_nm": 10.0}),
        _e("UniProtKB:HUB", "PHYSICALLY_INTERACTS_WITH", "UniProtKB:DEP2",
           {"detection_method": "m"}),
        _e("UniProtKB:DEP2", "HOST_FACTOR_FOR", V,
           {"direction": "dependency", "screen_type": "CRISPRko", "cell_line": "x"}),
        # ground truth, held out of traversal
        _e(D1, "HAS_ANTIVIRAL_ACTIVITY_AGAINST", V,
           {"assay_type": "cell_based_antiviral", "ec50_nm": 100.0,
            "relation": "="}),
    ]
    (tmp_path / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in edges))
    return tmp_path


@pytest.fixture
def schema():
    return load_schema()


@pytest.fixture
def g(release, schema):
    return Graph.load(release, skip_predicates=schema.held_out_predicates(),
                      symmetry=schema.symmetry())


# ------------------------------------------------------- held-out predicates
def test_held_out_edges_are_not_traversable(g):
    """HAS_ANTIVIRAL_ACTIVITY_AGAINST is the label. If it were traversable the
    baselines would be a lookup, not a prediction."""
    assert ("HAS_ANTIVIRAL_ACTIVITY_AGAINST", D1) not in g.out


def test_held_out_edges_still_count_toward_degree(g):
    """Degree is a property of the graph. Excluding label edges from it would
    make the degree baseline easier to beat than it really is."""
    assert g.degree[D1] >= 2


# ------------------------------------------------------------------ metapaths
def test_m1_direct_acting_path_is_found(g, schema):
    hits = walk_metapath(g, D1, schema.metapaths["M1"]["path"], None, set())
    assert V in hits and hits[V] > 0


def test_m2_requires_dependency_direction(g, schema):
    """A restriction factor must not satisfy M2: inhibiting one would HELP the
    virus."""
    mp = schema.metapaths["M2"]
    assert V in walk_metapath(g, D2, mp["path"], mp.get("constraints"), set())
    flipped = {"HOST_FACTOR_FOR.direction": "restriction"}
    assert walk_metapath(g, D2, mp["path"], flipped, set()) == {}


def test_m5_two_hop_host_path_is_found(g, schema):
    mp = schema.metapaths["M5"]
    assert V in walk_metapath(g, D3, mp["path"], mp.get("constraints"), set())


# ------------------------------------------------ reverse-oriented PPI (H1)
@pytest.fixture
def reverse_ppi(tmp_path):
    """HUB2 holds a PPI in EACH direction: an outgoing one to a protein that
    is not a host factor, and an incoming one from a dependency factor.

    That shape is the H1 bug. STRING stores each interaction once, ordered by
    accession, so 9,016 host proteins in v0.1 carry edges both ways; a walker
    that follows outgoing edges and consults incoming ones only when there
    are none sees the decoy, stops, and never finds the route.
    """
    nodes = [
        {"id": V, "class": "OrganismTaxon",
         "properties": {"label": "SARS-CoV-2", "family": "Coronaviridae",
                        "baltimore_class": "IV", "is_enveloped": True}},
        _prot("UniProtKB:HUB2", False, sym="HUB2"),
        _prot("UniProtKB:DECOY", False, sym="DECOY"),
        _prot("UniProtKB:DEP3", False, sym="DEP3"),
        _drug(D4),
    ]
    edges = [
        _e(D4, "TARGETS", "UniProtKB:HUB2",
           {"direction": "unknown", "assay_type": "binding", "ic50_nm": 20.0}),
        # forward, and a dead end
        _e("UniProtKB:HUB2", "PHYSICALLY_INTERACTS_WITH", "UniProtKB:DECOY",
           {"detection_method": "m"}),
        # stored DEP3 -> HUB2: against the direction M5 walks it
        _e("UniProtKB:DEP3", "PHYSICALLY_INTERACTS_WITH", "UniProtKB:HUB2",
           {"detection_method": "m"}),
        _e("UniProtKB:DEP3", "HOST_FACTOR_FOR", V,
           {"direction": "dependency", "screen_type": "CRISPRko", "cell_line": "x"}),
    ]
    (tmp_path / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in edges))
    return tmp_path


@pytest.fixture
def g_rev(reverse_ppi, schema):
    return Graph.load(reverse_ppi, skip_predicates=schema.held_out_predicates(),
                      symmetry=schema.symmetry())


def test_m5_finds_a_ppi_stored_against_the_walk_direction(g_rev, schema):
    """THE H1 REGRESSION. HUB2 holds one outgoing interaction (to a protein
    that is not a host factor) and one incoming (from a dependency factor).
    Inferring direction from `does this node have outgoing edges` sees only
    the decoy and loses the route entirely."""
    mp = schema.metapaths["M5"]
    hits = walk_metapath(g_rev, D4, mp["path"], mp.get("constraints"), set())
    assert V in hits, "inbound half of a symmetric predicate was not traversed"


def test_symmetric_predicates_are_walked_both_ways(g_rev):
    """Both adjacency maps, because which side a pair was stored on is an
    ingest artifact -- STRING canonicalises each interaction to one row."""
    both = g_rev.neighbours("PHYSICALLY_INTERACTS_WITH", "UniProtKB:HUB2")
    assert {n for n, _e in both} == {"UniProtKB:DECOY", "UniProtKB:DEP3"}


def test_directed_predicates_are_not_walked_backwards(g_rev):
    """TARGETS is directed: a protein must not reach back to its drug."""
    assert g_rev.neighbours("TARGETS", "UniProtKB:HUB2") == []
    assert [n for n, _e in g_rev.neighbours("TARGETS", D4)] == ["UniProtKB:HUB2"]


def test_a_reverse_hop_reads_the_inverse_map(g_rev):
    """M4's second PARTICIPATES_IN and M8's second MEMBER_OF_CLASS are
    reverse hops over a directed predicate; they used to work only because
    the walker guessed."""
    fwd = g_rev.neighbours("HOST_FACTOR_FOR", "UniProtKB:DEP3")
    rev = g_rev.neighbours("HOST_FACTOR_FOR", V, reverse=True)
    assert [n for n, _e in fwd] == [V]
    assert "UniProtKB:DEP3" in {n for n, _e in rev}


def test_traversing_an_undeclared_predicate_raises(reverse_ppi, schema):
    """Silence is what H1 was: a predicate with no declared symmetry must
    stop the walk, not fall back to inferring direction."""
    bare = Graph.load(reverse_ppi, symmetry={})
    with pytest.raises(ValueError, match="no declared symmetry"):
        bare.neighbours("PHYSICALLY_INTERACTS_WITH", "UniProtKB:HUB2")


def test_blocked_intermediates_are_skipped_but_endpoints_are_not(g, schema):
    """Hub exclusion applies to path connectors. The endpoint is the answer."""
    mp = schema.metapaths["M5"]
    assert walk_metapath(g, D3, mp["path"], mp.get("constraints"),
                         {"UniProtKB:HUB"}) == {}
    # blocking the endpoint must NOT remove the path
    assert V in walk_metapath(g, D3, mp["path"], mp.get("constraints"), {V})


def test_paths_do_not_revisit_nodes(g, schema):
    """A path looping through its own intermediate is not evidence."""
    mp = schema.metapaths["M5"]
    hits = walk_metapath(g, D3, mp["path"], mp.get("constraints"), set())
    assert all(v > 0 for v in hits.values())


# ----------------------------------------------------------------- weighting
def test_dwpc_downweights_high_degree_intermediates(g, schema):
    """Without this, path counts measure how well-studied the intermediates
    are rather than how specific the route is."""
    low = walk_metapath(g, D1, schema.metapaths["M1"]["path"], None, set(), w=0.0)
    high = walk_metapath(g, D1, schema.metapaths["M1"]["path"], None, set(), w=0.8)
    assert high[V] < low[V]
    assert low[V] == pytest.approx(1.0)


def test_combine_takes_the_best_route_not_the_sum(g, schema):
    """Summing across metapaths rewards APPEARING IN MANY POOLS rather than
    scoring well in any. Under sum-of-percentiles a promiscuous kinase
    inhibitor in all five pools reached ~4.0 while nirmatrelvir, with one
    clean M1 route, was capped at 1.0 and ranked 1,428."""
    d = dwpc_scores(g, schema.metapaths, schema.retrieval_limits)
    total = combine(d, V)
    assert set(total) == {D1, D2, D3}
    # max-percentile: every score is a percentile, so nothing exceeds 1.0
    assert all(0.0 <= v <= 1.0 for v in total.values())


def test_single_occupant_pool_scores_one_not_zero(g, schema):
    """Each fixture metapath has one drug. Percentile is undefined for a
    single-element pool, but the drug IS the top of it."""
    d = dwpc_scores(g, schema.metapaths, schema.retrieval_limits)
    total = combine(d, V)
    assert total[D1] == 1.0


# ----------------------------------------------------------------- baselines
def test_degree_ranking_ignores_the_virus(g):
    scores = degree_ranking(g)
    assert set(scores) == {D1, D2, D3}
    assert scores[D1] == float(g.degree[D1])


def test_network_proximity_returns_z_scores(g, schema):
    prox = network_proximity(g, V, schema.retrieval_limits, n_random=20)
    assert prox, "no drugs scored"
    assert all(isinstance(v, float) for v in prox.values())


def test_proximity_is_empty_without_dependency_factors(g, schema):
    assert network_proximity(g, "NCBITaxon:999", schema.retrieval_limits) == {}


# ------------------------------------------------------------------- metrics
def test_hits_and_mrr():
    ranked = [("a", 3.0), ("b", 2.0), ("c", 1.0)]
    assert hits_at_k(ranked, {"b"}, 1) == 0
    assert hits_at_k(ranked, {"b"}, 2) == 1
    assert mrr(ranked, {"b"}) == pytest.approx(0.5)
    assert mrr(ranked, {"zzz"}) == 0.0


def test_rank_is_descending(g):
    r = rank(degree_ranking(g))
    assert [s for _d, s in r] == sorted((s for _d, s in r), reverse=True)


# --------------------------------------------------------------- real release
REAL = Path(__file__).resolve().parents[1] / "data" / "releases" / "v0.1"


def test_real_graph_has_traversable_host_paths(real_graph):
    g = real_graph
    targets = sum(len(v) for (p, _n), v in g.out.items() if p == "TARGETS")
    assert targets > 5_000, f"only {targets} TARGETS edges -- M2-M6 need these"


# --------------------------------------------------------------------------
# A declared metapath that reaches nothing must be a REPORTED VALUE, not a
# missing key. M8 has produced zero paths since schema 0.8.0 declared it.
# --------------------------------------------------------------------------
def test_metapath_reach_reports_a_metapath_that_reaches_nothing(tmp_path):
    """The check that would have caught M8 on the day it was added.

    compare_metapaths printed "no paths" and then omitted the metapath from
    its results file; hard_negatives skipped it silently in both. So a
    metapath contributing nothing looked exactly like one that was never
    declared.
    """
    from kgav.baselines import Graph, dwpc_scores, metapath_reach

    nodes = [
        {"id": "NCBITaxon:1", "class": "OrganismTaxon",
         "properties": {"label": "V", "family": "Coronaviridae",
                        "baltimore_class": "IV", "is_enveloped": True}},
        {"id": "UniProtKB:P1", "class": "Protein",
         "properties": {"taxon_id": "NCBITaxon:1", "is_viral": True,
                        "sequence_hash": "h", "reviewed": True,
                        "protein_family": "nsp5"}},
        {"id": "KGAV:G", "class": "Gene",
         "properties": {"symbol": "rep", "taxon_id": "NCBITaxon:1",
                        "is_viral": True}},
        {"id": "INCHIKEY:D", "class": "SmallMolecule",
         "properties": {"smiles": "C", "inchikey_skel": "D", "is_approved": True,
                        "salt_collapsed": False, "stereo_collapsed": False}},
    ]

    def e(s, p, o, q=None):
        return {"subject": s, "predicate": p, "object": o, "qualifiers": q or {},
                "primary_knowledge_source": "s", "evidence_tier": 1,
                "first_asserted_date": "2020-01-01"}

    # M1 is walkable; LIVE_ONLY has a predicate with no edges in this graph.
    edges = [e("INCHIKEY:D", "INHIBITS", "UniProtKB:P1", {"assay_type": "biochemical"}),
             e("UniProtKB:P1", "ENCODED_BY", "KGAV:G"),
             e("KGAV:G", "BELONGS_TO", "NCBITaxon:1")]
    rel = tmp_path / "rel"
    rel.mkdir()
    (rel / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in nodes))
    (rel / "edges.jsonl").write_text("\n".join(json.dumps(x) for x in edges))

    metapaths = {
        "M1": {"path": ["SmallMolecule", "INHIBITS", "Protein", "ENCODED_BY",
                        "Gene", "BELONGS_TO", "OrganismTaxon"]},
        "MDEAD": {"path": ["SmallMolecule", "INHIBITS", "Protein",
                           "MEMBER_OF_CLASS", "TargetClass", "FOLD_SIMILAR_TO",
                           "TargetClass", "<MEMBER_OF_CLASS", "Protein",
                           "ENCODED_BY", "Gene", "BELONGS_TO", "OrganismTaxon"]},
    }
    symmetry = {"INHIBITS": False, "ENCODED_BY": False, "BELONGS_TO": False,
                "MEMBER_OF_CLASS": False, "FOLD_SIMILAR_TO": True}
    g = Graph.load(rel, symmetry=symmetry)
    dwpc = dwpc_scores(g, metapaths, {})
    reach = metapath_reach(g, metapaths, dwpc)

    # Every declared metapath appears, including the dead one.
    assert set(reach) == {"M1", "MDEAD"}
    assert reach["M1"]["drugs_reached"] == 1
    assert reach["MDEAD"]["drugs_reached"] == 0
    # And it says WHY, so nobody has to reconstruct it.
    assert set(reach["MDEAD"]["missing_predicates"]) == {
        "MEMBER_OF_CLASS", "FOLD_SIMILAR_TO"}
    assert reach["M1"]["missing_predicates"] == []


def test_metapath_reach_distinguishes_missing_edges_from_broken_traversal():
    """A self-loop is the M8 case: the predicate EXISTS, so the break is in the
    edges, not in an absent layer. The diagnosis must not claim otherwise."""
    from kgav.baselines import Graph, metapath_reach

    g = Graph(symmetry={"FOLD_SIMILAR_TO": True})
    g.out[("FOLD_SIMILAR_TO", "INTERPRO:X")] = [("INTERPRO:X", {})]
    mp = {"M": {"path": ["TargetClass", "FOLD_SIMILAR_TO", "TargetClass"]}}
    reach = metapath_reach(g, mp, {"M": {}})
    assert reach["M"]["drugs_reached"] == 0
    assert reach["M"]["missing_predicates"] == []


# --------------------------------------------------------------------------
# Lift is uninformative when the pool is no larger than k.
# --------------------------------------------------------------------------
def test_lift_on_a_pool_no_larger_than_k_is_exactly_pool_over_k():
    """Not an approximation -- an identity, and it was reported as a result.

    ranked[:k] is the whole pool, so hits == every positive in it and
    lift = p / (k*p/pool) = pool/k, independent of the ordering. On v0.1:
    MERS M1 pool 24 -> 0.24, HCoV-229E M1 pool 17 -> 0.17, MERS M7 pool 14 ->
    0.14, HCoV-229E M7 pool 7 -> 0.07. None was flagged, because the guard
    tested only expectation >= 1 and those expectations run to 57.
    """
    from kgav.baselines import lift_row

    for pool_size, n_pos in [(24, 2), (17, 7), (14, 8), (7, 1), (100, 19)]:
        scores = {f"D{i}": float(pool_size - i) for i in range(pool_size)}
        pos = {f"D{i}" for i in range(n_pos)}
        m = lift_row(scores, pos, k=100)
        assert m["lift"] == pytest.approx(pool_size / 100, abs=1e-9), \
            f"pool {pool_size}: lift should be the pool-size identity"
        assert m["lift_is_pool_artifact"]
        assert not m["trustworthy"], \
            f"pool {pool_size} <= k must never be reported as trustworthy"


def test_a_pool_larger_than_k_is_judged_on_expectation():
    """The guard must not become 'distrust everything'."""
    from kgav.baselines import lift_row

    scores = {f"D{i}": float(500 - i) for i in range(500)}
    strong = lift_row(scores, {f"D{i}" for i in range(50)}, k=100)
    assert strong["trustworthy"] and not strong["lift_is_pool_artifact"]
    assert strong["lift"] > 1.0          # ranking was genuinely exercised

    thin = lift_row(scores, {"D0", "D499"}, k=100)
    assert not thin["trustworthy"]       # expectation 0.4, still noise
    assert not thin["lift_is_pool_artifact"]


# ------------------------------------------- publication-disjoint evaluation
@pytest.fixture
def coincident(tmp_path):
    """Two drugs with identical topology; only the provenance differs.

    D1's INHIBITS edge cites the same PMID as its own antiviral label -- one
    paper reporting a target IC50 and a cell-based EC50 together, which is how
    94.9% of v0.1's labelled compounds look. D2's cite different papers.

    A publication-disjoint run must withhold D1's edge and keep D2's. If it
    withholds both, or neither, the filter is not reading provenance.
    """
    nodes = [_drug(D1), _drug(D2), _prot("UniProtKB:VIRAL", True),
             {"id": "NCBIGene:V", "class": "Gene",
              "properties": {"symbol": "orf1ab", "taxon_id": V, "is_viral": True}},
             {"id": V, "class": "OrganismTaxon",
              "properties": {"label": "SARS-CoV-2", "family": "Coronaviridae",
                             "baltimore_class": "IV", "is_enveloped": True}}]

    def pub(e, pmids):
        e["publications"] = pmids
        return e

    edges = [
        # D1: edge and label from ONE paper -> leakage
        pub(_e(D1, "INHIBITS", "UniProtKB:VIRAL"), ["PMID:111"]),
        pub(_e(D1, "HAS_ANTIVIRAL_ACTIVITY_AGAINST", V), ["PMID:111"]),
        # D2: edge and label from DIFFERENT papers -> sound evidence
        pub(_e(D2, "INHIBITS", "UniProtKB:VIRAL"), ["PMID:222"]),
        pub(_e(D2, "HAS_ANTIVIRAL_ACTIVITY_AGAINST", V), ["PMID:333"]),
        _e("UniProtKB:VIRAL", "ENCODED_BY", "NCBIGene:V"),
        _e("NCBIGene:V", "BELONGS_TO", V),
    ]
    (tmp_path / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in edges))
    return tmp_path


def test_label_publication_index_collects_only_label_edges(coincident, schema):
    idx = label_publication_index(coincident, schema.held_out_predicates())
    assert idx == {D1: {"PMID:111"}, D2: {"PMID:333"}}, idx
    # the INHIBITS PMIDs must NOT appear: the index describes labels, and
    # folding evidence into it would make every edge its own exclusion.
    assert "PMID:222" not in idx[D2]


def test_same_paper_evidence_is_withheld_and_other_papers_are_kept(coincident, schema):
    held = schema.held_out_predicates()
    idx = label_publication_index(coincident, held)
    g = Graph.load(coincident, skip_predicates=held, symmetry=schema.symmetry(),
                   drop_pubs_shared_with=idx)
    assert g.neighbours("INHIBITS", D1) == []          # same PMID as its label
    assert len(g.neighbours("INHIBITS", D2)) == 1      # different papers
    assert g.stats["dropped_publication_coincident"] == 1
    assert g.stats["dropped_INHIBITS"] == 1


def test_without_the_flag_both_drugs_are_traversable(coincident, schema):
    """The guard: this is the behaviour every published number was computed
    with, so it must still be reachable."""
    held = schema.held_out_predicates()
    g = Graph.load(coincident, skip_predicates=held, symmetry=schema.symmetry())
    assert len(g.neighbours("INHIBITS", D1)) == 1
    assert len(g.neighbours("INHIBITS", D2)) == 1
    assert g.stats["dropped_publication_coincident"] == 0


def test_withheld_edge_still_counts_toward_degree(coincident, schema):
    """Degree is a property of the graph. Dropping the edge from degree would
    make the degree baseline easier to beat, which is the opposite of a
    control -- same reasoning as held-out predicates."""
    held = schema.held_out_predicates()
    idx = label_publication_index(coincident, held)
    g = Graph.load(coincident, skip_predicates=held, symmetry=schema.symmetry(),
                   drop_pubs_shared_with=idx)
    assert g.degree[D1] == g.degree[D2] == 2


def test_m1_loses_the_leaky_drug_but_not_the_sound_one(coincident, schema):
    """The whole point, at the metapath level: D1 scores only because its
    evidence and its answer came from one paper."""
    held = schema.held_out_predicates()
    idx = label_publication_index(coincident, held)
    m1 = schema.metapaths["M1"]
    g_leak = Graph.load(coincident, skip_predicates=held, symmetry=schema.symmetry())
    g_clean = Graph.load(coincident, skip_predicates=held, symmetry=schema.symmetry(),
                         drop_pubs_shared_with=idx)
    for d in (D1, D2):
        assert walk_metapath(g_leak, d, m1["path"], m1.get("constraints"), set())
    assert not walk_metapath(g_clean, D1, m1["path"], m1.get("constraints"), set())
    assert walk_metapath(g_clean, D2, m1["path"], m1.get("constraints"), set())


# ------------------- lift needs a null, not just a readability heuristic
def test_the_lift_pvalue_matches_exact_integer_arithmetic():
    """Computed in log space because pool 4,316 choose 100 is a 240-digit
    integer. Checked against the exact rational value on the rows the 2021
    temporal split actually produced, so a precision regression shows up as a
    wrong p and not as a silent drift.
    """
    from math import comb

    from kgav.baselines import lift_pvalue

    def exact(N, K, k, h):
        return sum(comb(K, i) * comb(N - K, k - i)
                   for i in range(h, min(K, k) + 1)) / comb(N, k)

    for N, K, k, h in ((296, 8, 100, 6), (186, 4, 100, 3), (3246, 50, 100, 4),
                       (4316, 231, 100, 5), (2682, 5, 100, 3), (203, 2, 100, 2)):
        got, want = lift_pvalue(N, K, k, h), exact(N, K, k, h)
        assert abs(got - want) < 1e-9, f"N={N} K={K} h={h}: {got} vs {want}"


def test_the_null_is_hypergeometric_not_poisson():
    """Taking the top k from a pool is sampling WITHOUT replacement, so the
    null has the finite-population correction and a TIGHTER tail than Poisson.

    The direction matters and is counter-intuitive: Poisson is CONSERVATIVE
    here, so using it would discard real results rather than invent them. My
    first draft of this test asserted the opposite and failed, which is the
    reason the direction is pinned by a test at all.
    """
    from math import exp, factorial

    from kgav.baselines import lift_pvalue

    N, K, k, h = 120, 20, 100, 20          # every positive inside the top k
    lam = k * K / N
    poisson = 1.0 - sum(exp(-lam) * lam ** i / factorial(i) for i in range(h))
    exact = lift_pvalue(N, K, k, h)

    assert exact < 0.05, "all 20 positives in the top 100 of 120 is unlikely"
    assert poisson > 0.2, "and Poisson calls the same row unremarkable"
    assert poisson > 10 * exact, f"Poisson is {poisson / exact:.0f}x too high"


def test_an_empty_question_gets_no_pvalue():
    """None rather than 1.0: 'every draw is identical' is not 'consistent with
    chance', and a reader given 1.0 would think the question was asked.
    """
    from kgav.baselines import lift_pvalue
    assert lift_pvalue(100, 5, 100, 5) is None, "k covers the whole pool"
    assert lift_pvalue(91, 29, 100, 29) is None, "k exceeds the pool"
    assert lift_pvalue(500, 0, 100, 0) is None, "no positives to find"


def test_the_heuristic_and_the_pvalue_are_reported_separately():
    """MERS-CoV M4 is the case: pool 2,682, 5 positives, 3 hits, expectation
    0.19 -- flagged untrustworthy by `exp >= 1`, and the least likely row in
    the whole run at p=0.0005. Folding the p-value into `trustworthy` would
    have hidden it; dropping the heuristic would have let lift 16.09 be read
    as an effect size. Both are reported.
    """
    from kgav.baselines import lift_pvalue, lift_row

    scores = {f"D{i}": 1.0 / (i + 1) for i in range(2682)}
    pos = {"D0", "D1", "D2", "D2000", "D2500"}        # 3 of 5 inside top 100
    m = lift_row(scores, pos, k=100)

    assert m["pos"] == 5 and m["hits"] == 3
    assert m["expected"] < 1.0
    assert m["trustworthy"] is False, "the readability heuristic still fires"
    assert m["p_value"] < 0.001, "while the exact null says it is surprising"
    assert m["p_value"] == lift_pvalue(2682, 5, 100, 3)


# ------------- a p-value is against random, not against the confound
def _row(lift, p):
    return {"lift": lift, "p_value": p}


def test_the_confound_warning_fires_when_degree_is_itself_significant():
    """SARS-CoV cross-sectional on v0.1: M1 lift 2.10 p<0.0001, and degree
    lift 9.51 p<0.0001 on the SAME labels. In every cross-sectional virus
    where a channel cleared 0.05, degree cleared it too -- 4 of 4, no
    exception. A reader shown only M1's p would call that a discovery.
    """
    from kgav.baselines import confound_warning

    per = {"M1": _row(2.10, 0.00001), "M3": _row(1.90, 0.1424),
           "degree": _row(9.51, 0.00001)}
    warn = confound_warning(per)
    assert warn is not None
    assert "CONFOUNDED" in warn
    assert "no channel exceeds it" in warn, \
        "degree outranks M1 4.5x here, which is the whole point"


def test_the_baseline_read_states_the_good_case_out_loud():
    """SARS-CoV-2 temporal: M1 lift 2.22 p=0.020 with degree at 0.93 p=0.63.
    Degree at chance is the condition every real result here depends on, and
    it holds in exactly one place in the project. Printing it only as an
    ABSENCE of warning left the reader to notice it, so it is now stated.
    """
    from kgav.baselines import confound_warning

    per = {"M1": _row(2.22, 0.0196), "M6": _row(1.40, 0.3693),
           "degree": _row(0.93, 0.6276)}
    warn = confound_warning(per)
    assert warn is not None and "BASELINE AT CHANCE" in warn
    assert "M1's lift 2.22" in warn and "not explained by study volume" in warn
    assert "CONFOUNDED" not in warn


def test_a_baseline_at_chance_with_nothing_significant_says_nothing():
    """SARS-CoV-2 cross-sectional: degree 0.45 p=0.99, best channel M4 1.99
    p=0.14. Neither confounded nor a result, so there is nothing to report.
    """
    from kgav.baselines import confound_warning

    per = {"M1": _row(0.86, 0.8141), "M4": _row(1.99, 0.1401),
           "degree": _row(0.45, 0.9902)}
    assert confound_warning(per) is None


def test_a_ratio_against_a_chance_baseline_is_marked_unreadable():
    """M4 reads 4.44x a degree lift of 0.45 on SARS-CoV-2 cross-sectional,
    which sounds decisive and is lift 1.99 at p=0.14 -- the ratio divides by
    noise. My first threshold for this was `base >= 1.0`, which wrongly
    flagged the clean temporal row where degree sits at 0.93: that is degree
    at chance, the condition a real result needs, not a defect in it. The
    test is whether the DENOMINATOR is a real effect.
    """
    from kgav.baselines import degree_relative

    noisy = {"M4": _row(1.99, 0.1401), "degree": _row(0.45, 0.9902)}
    degree_relative(noisy)
    assert noisy["M4"]["vs_degree_readable"] is False
    assert abs(noisy["M4"]["lift_vs_degree"] - 1.99 / 0.45) < 1e-12

    real = {"M1": _row(2.10, 0.0001), "degree": _row(9.51, 0.0001)}
    degree_relative(real)
    assert real["M1"]["vs_degree_readable"] is True


def test_the_warning_names_a_channel_that_does_beat_degree():
    from kgav.baselines import confound_warning

    per = {"M4": _row(11.86, 0.00001), "degree": _row(10.87, 0.00001)}
    warn = confound_warning(per)
    assert "M4 exceeds it at 11.86" in warn, warn


def test_degree_relative_annotates_every_row_against_the_baseline():
    from kgav.baselines import degree_relative

    per = {"M1": _row(2.10, 0.0001), "degree": _row(9.51, 0.0001)}
    degree_relative(per)
    assert abs(per["M1"]["lift_vs_degree"] - 2.10 / 9.51) < 1e-12
    assert per["M1"]["beats_degree"] is False
    assert per["degree"]["lift_vs_degree"] is None, \
        "the baseline is not compared to itself"


def test_a_zero_baseline_reports_the_win_instead_of_hiding_it():
    """HCoV-229E temporal computed_only: degree lift 0.00 and M2 lift 4.65.
    The ratio is undefined, and the first version of this sent it down the
    no-baseline branch, so the one row where a channel beat degree outright
    printed as "-" -- indistinguishable from the baseline row itself.

    Carried by beats_degree and NOT by float("inf"), which json.dumps writes
    as `Infinity` and would make every results file invalid JSON.
    """
    import json
    import math

    from kgav.baselines import degree_relative

    per = {"M2": _row(4.65, 0.2151), "M4": _row(0.0, 1.0),
           "degree": _row(0.0, 1.0)}
    degree_relative(per)

    assert per["M2"]["beats_degree"] is True
    assert per["M2"]["lift_vs_degree"] is None
    assert per["M4"]["beats_degree"] is False
    assert per["degree"]["beats_degree"] is None, "the baseline row is distinct"

    for m in per.values():
        v = m["lift_vs_degree"]
        assert v is None or math.isfinite(v)
    json.loads(json.dumps(per))          # must stay valid JSON


def test_no_baseline_row_at_all_raises_nothing():
    """HCoV-HKU1 has no measured activity, so a table can arrive with no
    degree row."""
    from kgav.baselines import confound_warning, degree_relative

    per = {"M1": _row(2.0, 0.01)}
    degree_relative(per)
    assert per["M1"]["lift_vs_degree"] is None
    assert per["M1"]["beats_degree"] is None
    assert confound_warning(per) is None

    # A baseline with no p-value cannot be judged at chance OR confounded,
    # so neither verdict is given. This also pins a crash: the at-chance
    # branch formats dp, and lift_pvalue returns None for a row with no
    # positives or where k covers the pool, so this raised TypeError.
    unstamped = {"M1": _row(2.0, 0.01), "degree": {"lift": 1.0,
                                                   "p_value": None}}
    degree_relative(unstamped)
    assert unstamped["M1"]["vs_degree_readable"] is False
    assert confound_warning(unstamped) is None
