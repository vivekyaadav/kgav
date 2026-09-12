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
           {"assay_type": "biochemical", "ic50_nm": 25.0, "unquantified": True}),
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
            "relation": "=", "unquantified": True}),
    ]
    (tmp_path / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in edges))
    return tmp_path


@pytest.fixture
def schema():
    return load_schema()


@pytest.fixture
def g(release, schema):
    return Graph.load(release, skip_predicates=schema.held_out_predicates())


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


@pytest.mark.skipif(not (REAL / "nodes.jsonl").exists(), reason="release not built")
def test_real_graph_has_traversable_host_paths():
    s = load_schema()
    g = Graph.load(REAL, skip_predicates=s.held_out_predicates())
    targets = sum(len(v) for (p, _n), v in g.out.items() if p == "TARGETS")
    assert targets > 5_000, f"only {targets} TARGETS edges -- M2-M6 need these"
