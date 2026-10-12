"""Gate for the baselines: held-out edges are never traversed, DWPC
down-weights promiscuous intermediates, proximity is degree-corrected.
"""
import json
import math
from pathlib import Path

import pytest

from kgav.baselines import (
    Graph,
    channel_percentiles,
    combine,
    confound_warning,
    degree_ranking,
    dwpc_scores,
    hits_at_k,
    hits_at_k_range,
    hits_cell,
    label_publication_index,
    lift_row,
    mrr,
    network_proximity,
    potency_key,
    potency_tiebreak,
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
    """A fixture edge, mirroring one production invariant deliberately.

    chembl.py sets target_resolution on EVERY INHIBITS edge it emits -- both
    branches, resolved and unresolved -- and the schema now requires it,
    because M1/M7/M8 constrain on it and a missing qualifier reads as a
    failed match. A fixture omitting it is not a simpler fixture; it is an
    edge the direct-acting channel silently refuses to walk. Defaulted here
    so that stays true as fixtures are added, and overridable so a test can
    say parent_unresolved and mean it.
    """
    quals = dict(q or {})
    if p == "INHIBITS":
        quals.setdefault("target_resolution", "chain_from_assay_description")
    return {"subject": s, "predicate": p, "object": o, "qualifiers": quals,
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


# ------------------------------------------------------------- tied percentiles
def _one_channel(values: dict[str, float]) -> dict:
    """A single metapath pool holding exactly these raw DWPC values."""
    return {"M1": {drug: {V: v} for drug, v in values.items()}}


def test_identical_dwpc_gets_an_identical_percentile():
    """`pct = i / n` read the POSITIONAL index, so equal scores were ordered
    by whatever order sorted() left them in -- for equal keys, insertion
    order. Measured on v0.1: nirmatrelvir and chloroquine both scored
    0.000242 on M1, both by the single path through nsp5, and were handed
    0.0180 and 0.0020. The reported rank 6,432 was an artifact of iteration
    order, not a measurement.
    """
    tied = combine(_one_channel({"a": 1.0, "nirmatrelvir": 2.0,
                                 "chloroquine": 2.0, "d": 3.0}), V)
    assert tied["nirmatrelvir"] == tied["chloroquine"]


def test_reinserting_the_same_pool_in_another_order_changes_nothing():
    """The property the positional index broke: the output must depend on the
    scores, not on the order the drugs were added to the dict."""
    vals = {"a": 1.0, "b": 2.0, "c": 2.0, "d": 2.0, "e": 3.0}
    forward = combine(_one_channel(vals), V)
    backward = combine(_one_channel(dict(reversed(vals.items()))), V)
    assert forward == backward


def test_a_tied_group_gets_the_mean_of_the_positions_it_spans():
    """Midrank, so the tie neither gains nor loses against the drugs around
    it: three drugs tied across positions 1, 2 and 3 of 5 all score 2/4."""
    out = combine(_one_channel({"lo": 1.0, "t1": 2.0, "t2": 2.0, "t3": 2.0,
                                "hi": 3.0}), V)
    assert out["lo"] == pytest.approx(0.0)
    assert out["hi"] == pytest.approx(1.0)
    for d in ("t1", "t2", "t3"):
        assert out[d] == pytest.approx(0.5)


def test_untied_pools_are_unchanged():
    """The fix must not move a pool that had no ties, or every published
    COMBINED figure would shift for no reason."""
    out = combine(_one_channel({"a": 1.0, "b": 2.0, "c": 3.0}), V)
    assert out["a"] == pytest.approx(0.0)
    assert out["b"] == pytest.approx(0.5)
    assert out["c"] == pytest.approx(1.0)


def test_a_wholly_degenerate_pool_hands_nobody_a_perfect_score():
    """If every drug in a channel scores the same, the channel separates
    nothing. The positional index awarded 1.0 to whichever drug sorted last
    and 0.0 to whichever sorted first -- a top-of-pool score, drawn from a
    pool with no top. All of them now sit at the midpoint."""
    out = combine(_one_channel({c: 7.0 for c in "abcdef"}), V)
    assert len(out) == 6
    assert all(v == pytest.approx(0.5) for v in out.values())
    assert 1.0 not in out.values() and 0.0 not in out.values()


def test_combine_cannot_diverge_from_channel_percentiles():
    """THE DIVERGENCE THIS EXISTS TO PREVENT. The percentile arithmetic once
    lived inside combine while scripts/diagnose_controls.py carried a second
    copy so it could print the percentile beside the combined score. combine
    moved to midranks; the copy did not. The diagnostic then reported a
    percentile of 0.0182 next to a combined score of 0.4041 for a drug whose
    only route is a single channel -- two numbers that cannot both be right,
    printed by the script written to make the scoring checkable.

    On one channel with no weighting, the combined score IS the channel
    percentile. Asserting equality here fails the moment either side grows
    its own arithmetic again.
    """
    pool = {"a": 1.0, "b": 2.0, "c": 2.0, "d": 2.0, "e": 9.0, "f": 0.5}
    assert combine(_one_channel(pool), V) == channel_percentiles(pool)


def test_channel_percentiles_handles_the_degenerate_shapes():
    """Empty, single-occupant, and all-identical, since each one used to be a
    special case written separately in each copy."""
    assert channel_percentiles({}) == {}
    assert channel_percentiles({"only": 4.0}) == {"only": 1.0}
    flat = channel_percentiles({c: 7.0 for c in "abcd"})
    assert all(v == pytest.approx(0.5) for v in flat.values())


# ------------------------------------ the baseline comparison needs determinacy
def _cw_row(lift, p, det=True, lo=0, hi=0, **kw):
    return {"lift": lift, "p_value": p, "trustworthy": True,
            "hits": hi, "hits_lo": lo, "hits_hi": hi,
            "hits_determined": det, "tied_at_cut": 0,
            "tied_at_cut_positives": 0, **kw}


def test_an_undetermined_baseline_cannot_license_a_channel():
    """THE REAL SARS-CoV-2 train2021 SHAPE. confound_warning printed
    "BASELINE AT CHANCE: degree lift 0.97 (p=0.5909) is indistinguishable
    from random here, so M1's lift 2.49 (p=0.0033) is not explained by study
    volume" -- where M1's Hits@100 was 0-8 of 8 positives and degree's was
    3-5. A sentence comparing two intervals as though they were points, and
    it licensed the only non-null in the project.
    """
    per = {"M1": _cw_row(2.49, 0.0033, det=False, lo=0, hi=8),
           "degree": _cw_row(0.97, 0.5909, det=False, lo=3, hi=5,
                          vs_degree_readable=True)}
    msg = confound_warning(per)
    assert msg is not None
    assert "BASELINE NOT DETERMINED" in msg
    assert "AT CHANCE" not in msg


def test_a_determined_baseline_still_licenses_a_determined_channel():
    """So the guard is a condition and not a blanket refusal -- the good case
    must still be sayable."""
    per = {"M1": _cw_row(2.49, 0.0033),
           "degree": _cw_row(0.97, 0.5909, vs_degree_readable=True)}
    msg = confound_warning(per)
    assert msg is not None and "BASELINE AT CHANCE" in msg


def test_an_undetermined_channel_is_not_promoted():
    """Third time this function has read one field and ignored the one beside
    it: it read p and ignored `trustworthy`, promoting HCoV-OC43 M4 on a
    single positive; then read both and ignored `hits_determined`."""
    per = {"M1": _cw_row(2.49, 0.0033, det=False, lo=0, hi=8),
           "degree": _cw_row(0.97, 0.5909, vs_degree_readable=True)}
    assert confound_warning(per) is None


def test_an_undetermined_exceedance_is_flagged_in_the_confounded_branch():
    """When degree is itself significant, a channel said to exceed it must
    not be reported as a clean exceedance if its own hits are an interval."""
    per = {"M4": _cw_row(16.09, 0.0005, det=False, lo=0, hi=5),
           "degree": _cw_row(4.84, 0.0231, vs_degree_readable=True)}
    msg = confound_warning(per)
    assert "CONFOUNDED" in msg and "flagged (!)" in msg


# ------------------------------------------------- potency as a tie-break
def _blocky():
    """A pool shaped like the real M1: one big tie block with groups above
    and below it. 13 entries, so the block holds positions 1..10 of n=12 and
    its span is 1/12 .. 10/12."""
    return {"lo": 0.1, **{f"t{i}": 0.5 for i in range(10)},
            "hi1": 0.9, "hi2": 1.5}


def test_nothing_outside_a_tie_block_moves():
    """The design's central promise, per docs/prereg_m1_potency.md Amendment
    1: block boundaries come from the untouched raw DWPC sort, so every
    cross-target ordering -- remdesivir above the nsp5 group included -- is
    preserved exactly, and COMBINED does not change."""
    pool = _blocky()
    before = channel_percentiles(pool)
    after = channel_percentiles(pool, {f"t{i}": 5.0 + i for i in range(6)})
    for d in ("lo", "hi1", "hi2"):
        assert before[d] == after[d]


def test_the_block_keeps_its_own_span():
    """The keyed drugs spread across i..j and no further, so the block
    occupies the same percentile range it did. On v0.1 that range is roughly
    0.002 to 0.806 -- which is why Amendment 1 states that nirmatrelvir
    cannot exceed ~0.81 by construction."""
    pool = _blocky()
    n = len(pool) - 1
    after = channel_percentiles(pool, {f"t{i}": 5.0 + i for i in range(6)})
    blk = [after[f"t{i}"] for i in range(10)]
    assert min(blk) == pytest.approx(1 / n)
    assert max(blk) == pytest.approx(10 / n)


def test_an_unmeasured_drug_keeps_the_score_it_already_has():
    """A strict refinement: no imputation decides anything, and "unmeasured"
    is never read as "weak" -- the error chembl.py refuses when it declines
    to coerce a censored bound into a number."""
    pool = _blocky()
    mid = channel_percentiles(pool)["t0"]
    after = channel_percentiles(pool, {f"t{i}": 5.0 + i for i in range(6)})
    for i in range(6, 10):
        assert after[f"t{i}"] == mid


def test_the_block_mean_is_unchanged_so_it_neither_gains_nor_loses():
    """Spreading a block must redistribute it, not shift it. If the mean
    moved, tie-breaking would be smuggling in a score change."""
    pool = _blocky()
    before = channel_percentiles(pool)
    after = channel_percentiles(pool, {f"t{i}": 5.0 + i for i in range(6)})
    b = sum(before[f"t{i}"] for i in range(10)) / 10
    a = sum(after[f"t{i}"] for i in range(10)) / 10
    assert a == pytest.approx(b)


def test_keyed_drugs_come_out_monotone_in_potency():
    pool = _blocky()
    after = channel_percentiles(pool, {f"t{i}": 5.0 + i for i in range(6)})
    seq = [after[f"t{i}"] for i in range(6)]
    assert seq == sorted(seq)


def test_a_single_keyed_drug_cannot_rank_a_block():
    """One value is not an ordering. Everyone keeps the midrank rather than
    one drug being promoted for being the only one measured."""
    pool = _blocky()
    before = channel_percentiles(pool)
    after = channel_percentiles(pool, {"t3": 9.0})
    assert after == before


def test_ties_in_the_key_itself_take_their_own_midrank():
    """The same rule one level up. Equal potency is equal evidence."""
    pool = _blocky()
    after = channel_percentiles(pool, {f"t{i}": 5.0 for i in range(6)})
    assert len({round(after[f"t{i}"], 9) for i in range(6)}) == 1


def test_potency_key_reads_the_declared_precedence_not_the_best_value():
    """nirmatrelvir on v0.1 carries ic50_nm=75, ki_nm=3 and kd_nm=190. Taking
    the most flattering would read 3; taking the worst would read 190. The
    FIELD is fixed in advance -- ic50, then ki, then kd -- because
    selectivity.py records both of those as "the same error: selecting a
    study by its answer"."""
    e = _e("INCHIKEY:N", "INHIBITS", "UniProtKB:V",
           {"assay_type": "biochemical", "ic50_nm": 75.0, "ki_nm": 3.0,
            "kd_nm": 190.0})
    px, field = potency_key(e)
    assert field == "ic50_nm"
    assert px == pytest.approx(9 - math.log10(75.0))


def test_a_censored_bound_is_not_a_measurement():
    """chembl.py: "IC50 > 10000 nM is a NEGATIVE result ... coercing it to a
    potency number inverts its meaning." A bound cannot break a tie."""
    for rel in (">", "<", ">=", "<=", "~"):
        e = _e("INCHIKEY:N", "INHIBITS", "UniProtKB:V",
               {"assay_type": "biochemical", "ic50_nm": 10000.0,
                "relation": rel})
        assert potency_key(e) is None


def test_an_edge_with_no_potency_field_yields_no_key():
    assert potency_key(_e("INCHIKEY:N", "INHIBITS", "UniProtKB:V",
                          {"assay_type": "biochemical"})) is None


def test_a_zero_potency_is_skipped_rather_than_logged():
    """The schema permits ic50_nm >= 0 and log10(0) does not exist. Skipped,
    not crashed and not clamped to a number that would outrank everything."""
    assert potency_key(_e("INCHIKEY:N", "INHIBITS", "UniProtKB:V",
                          {"assay_type": "biochemical",
                           "ic50_nm": 0.0})) is None


def test_potency_tiebreak_keeps_the_worst_of_several_measurements(g, schema):
    """selectivity.py's rule for the same situation: "among several exact
    measurements the WORST is kept, since a compound shown cytotoxic in a
    real measurement is cytotoxic." The direction that cannot flatter."""
    keys, fields = potency_tiebreak(
        g, V, schema.metapaths["M1"]["path"],
        schema.metapaths["M1"].get("constraints"))
    # the fixture's D1 inhibits NSP5 with ic50_nm 25.0
    assert keys[D1] == pytest.approx(9 - math.log10(25.0))
    assert fields["ic50_nm"] == 1


def test_potency_tiebreak_ignores_edges_the_channel_cannot_walk(g, schema):
    """Keyed on the evidence the channel uses, not on every INHIBITS edge in
    the graph: a protein with no onward route to this virus is not part of
    M1 whatever its qualifiers say."""
    keys, _f = potency_tiebreak(
        g, "NCBITaxon:99999", schema.metapaths["M1"]["path"],
        schema.metapaths["M1"].get("constraints"))
    assert keys == {}


def test_combine_applies_a_tiebreak_per_channel_not_globally(schema):
    """A tiebreak is a claim about one kind of evidence. Potency against a
    viral protein orders M1; it says nothing about the host channels, and
    nothing about M7, which borrows a NEIGHBOUR's target."""
    d = {"M1": {f"t{i}": {V: 0.5} for i in range(4)},
         "M2": {f"t{i}": {V: 0.5} for i in range(4)}}
    keys = {f"t{i}": 5.0 + i for i in range(4)}
    out = combine(d, V, tiebreak={"M1": keys})
    # M1 ordered -> t3 strictly above t0; M2 still a flat tie at 0.5, and
    # combine takes the MAX, so t0 keeps M2's midrank
    assert out["t3"] > out["t0"]
    assert out["t0"] == pytest.approx(0.5)


def test_omitting_the_tiebreak_reproduces_every_published_number(g, schema):
    """The default path must be untouched, or the change would silently move
    figures it was never argued to move."""
    d = dwpc_scores(g, schema.metapaths, schema.retrieval_limits)
    assert combine(d, V) == combine(d, V, tiebreak=None)
    assert combine(d, V) == combine(d, V, tiebreak={})


# ------------------------------------------- ties the midrank does NOT reach
def test_hits_at_k_is_decided_by_dict_order_when_a_tie_straddles_k():
    """A SEPARATE DEFECT FROM THE PERCENTILE ONE, one level up and still open.

    rank() is `sorted(scores.items(), key=lambda kv: -kv[1])` with no
    tiebreak, and Python's sort is stable, so equal scores come out in
    insertion order. hits_at_k reads ranked[:k]. Where a tie block straddles
    k, which members land inside the cut is decided by insertion order rather
    than by any score -- so Hits@k, MRR and lift are reproducible for one
    dict and not properties of the graph.

    This test documents the defect rather than fixing it: the fix is to
    report Hits@k as the range it could take, and that changes every
    published metric, so it waits on the measurement.
    """
    scores_a = {"hit": 5.0, "miss1": 5.0, "miss2": 5.0, "low": 1.0}
    scores_b = {"miss1": 5.0, "miss2": 5.0, "hit": 5.0, "low": 1.0}
    pos = {"hit"}
    assert hits_at_k(rank(scores_a), pos, 1) == 1
    assert hits_at_k(rank(scores_b), pos, 1) == 0
    assert mrr(rank(scores_a), pos) == pytest.approx(1.0)
    assert mrr(rank(scores_b), pos) == pytest.approx(1 / 3)


def test_hits_at_k_range_brackets_every_tie_ordering():
    """Checked against brute force over every permutation of the tie block,
    because an interval that is merely plausible is worse than no interval:
    it would license the figure it was meant to qualify."""
    import itertools
    scores = {"top": 9.0, "t1": 5.0, "t2": 5.0, "t3": 5.0, "t4": 5.0,
              "low": 1.0}
    pos = {"t4"}
    lo, hi, n_tied, n_pos = hits_at_k_range(rank(scores), pos, 4)
    seen = set()
    for perm in itertools.permutations(["t1", "t2", "t3", "t4"]):
        d = {"top": 9.0}
        for t in perm:
            d[t] = 5.0
        d["low"] = 1.0
        seen.add(hits_at_k(rank(d), pos, 4))
    assert (lo, hi) == (min(seen), max(seen))
    assert (n_tied, n_pos) == (4, 1)


def test_a_cut_in_no_tie_is_determined():
    """So the interval is a measurement and not a blanket disclaimer. On v0.1
    COMBINED is determined at k=10, 50 and 100 while M1 is not at any of
    them -- the reason this is reported per row."""
    scores = {"a": 5.0, "b": 4.0, "c": 3.0, "d": 2.0}
    lo, hi, n_tied, _p = hits_at_k_range(rank(scores), {"a"}, 2)
    assert lo == hi == 1
    assert n_tied == 0, "a block of one is not a tie"


def test_a_pool_no_larger_than_k_has_nothing_to_contest():
    """ranked[:k] is then the whole pool, so every positive in it is counted
    whatever the order -- the same condition lift_row already flags as
    lift_is_pool_artifact."""
    scores = {"a": 1.0, "b": 1.0, "c": 1.0}
    assert hits_at_k_range(rank(scores), {"a", "b"}, 10) == (2, 2, 0, 0)


def test_lift_row_carries_the_range_beside_the_figure():
    """Reported the way lift_is_pool_artifact is: alongside, not folded in.
    The figure is still the one computed; the range says whether it was
    determined, and the lift inherits the same spread."""
    scores = {"top": 9.0, "t1": 5.0, "t2": 5.0, "t3": 5.0, "t4": 5.0,
              "low": 1.0}
    m = lift_row(scores, {"t1", "t4"}, k=2)
    assert m["hits_determined"] is False
    assert m["hits_lo"] <= m["hits"] <= m["hits_hi"]
    assert m["lift_lo"] <= m["lift"] <= m["lift_hi"]
    assert m["tied_at_cut"] == 4 and m["tied_at_cut_positives"] == 2


def test_hits_cell_prints_a_range_only_when_undetermined():
    """One formatter, imported by both evaluation scripts. This project has
    twice had one calculation in two files and watched them drift."""
    assert hits_cell({"hits": 6, "hits_determined": True}, 5).strip() == "6"
    undet = {"hits": 7, "hits_determined": False, "hits_lo": 2, "hits_hi": 12}
    assert hits_cell(undet, 5).strip() == "2-12"


def test_ties_are_merged_only_when_exactly_equal():
    """Exact equality, not a tolerance: identical topology yields
    bit-identical products, while a tolerance would merge genuinely
    different scores and add a fitted parameter to a no-training baseline."""
    out = combine(_one_channel({"a": 1.0, "b": 2.0, "c": 2.0 + 1e-12,
                                "d": 3.0}), V)
    assert out["b"] != out["c"]


# ----------------------------------------------------------------- baselines
def test_degree_ranking_ignores_the_virus(g):
    scores = degree_ranking(g)
    assert set(scores) == {D1, D2, D3}
    assert scores[D1] == float(g.degree[D1])


def test_degree_ranking_covers_drugs_only_and_g_degree_covers_proteins(g):
    """degree_ranking is a DRUG ranking. A protein looked up in it returns
    nothing, so `degree_ranking(g).get(protein, 0)` reports 0 for a node that
    has a perfectly good degree -- a silent zero, not an error.

    The distinction is load-bearing for anything explaining DWPC, whose
    denominator is the degree of the INTERMEDIATE nodes a path traverses:
    those are proteins, so they must be read off g.degree. scripts/
    diagnose_controls.py printed a 0-degree column for every target until
    this was caught, which would have argued the opposite of the truth --
    that nsp5 is unstudied.
    """
    prot = "UniProtKB:NSP5"
    assert g.nodes[prot]["class"] == "Protein"
    assert g.degree[prot] > 0, "the fixture protein should carry edges"
    assert prot not in degree_ranking(g)
    assert degree_ranking(g).get(prot, 0) == 0  # the trap, pinned


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
    edges = [e("INCHIKEY:D", "INHIBITS", "UniProtKB:P1",
               {"assay_type": "biochemical", "target_resolution": "chain_from_assay_description"}),
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


# ------------------------------------------------- unresolved target exclusion
@pytest.fixture
def resolution(tmp_path):
    """Two drugs, same biology, filed at two levels of detail.

    RESOLVED inhibits the mature chain, as nirmatrelvir's ChEMBL activity
    does. VAGUE inhibits the parent polyprotein because its assay
    description named no protein, as pomotrelvir's does. The parent is given
    a LOWER degree than the chain, which is the real situation on v0.1
    (38 against 3,117) and the reason DWPC preferred the vaguer edge.
    """
    chain, parent = "UniProtKB:PRO_CHAIN", "UniProtKB:PARENT"
    nodes = [_drug("INCHIKEY:RESOLVED"), _drug("INCHIKEY:VAGUE"),
             _prot(chain, True, fam="nsp5"), _prot(parent, True, fam="pp1ab"),
             {"id": "KGAV:GENE_rep", "class": "Gene",
              "properties": {"symbol": "rep", "taxon_id": V, "is_viral": True}},
             {"id": V, "class": "OrganismTaxon",
              "properties": {"label": "SARS-CoV-2", "family": "Coronaviridae",
                             "baltimore_class": "IV", "is_enveloped": True}}]
    edges = [
        _e("INCHIKEY:RESOLVED", "INHIBITS", chain, {"assay_type": "biochemical"}),
        _e("INCHIKEY:VAGUE", "INHIBITS", parent,
           {"assay_type": "biochemical",
            "target_resolution": "parent_unresolved"}),
        _e(chain, "ENCODED_BY", "KGAV:GENE_rep"),
        _e(parent, "ENCODED_BY", "KGAV:GENE_rep"),
        _e("KGAV:GENE_rep", "BELONGS_TO", V),
    ]
    # Padding so the CHAIN outranks the PARENT in degree, as it does on v0.1:
    # the parent's degree is low precisely BECAUSE resolution succeeded
    # elsewhere, so it measures how few activities failed to resolve.
    for i in range(12):
        d = f"INCHIKEY:PAD{i:02d}"
        nodes.append(_drug(d))
        edges.append(_e(d, "INHIBITS", chain, {"assay_type": "biochemical"}))
    rel = tmp_path / "rel"
    rel.mkdir()
    (rel / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in nodes))
    (rel / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in edges))
    return rel


def test_an_unresolved_target_does_not_reach_the_direct_acting_channel(
        resolution, schema):
    """ChEMBL has no chain-level coronavirus targets, so chembl.py recovers
    the nsp from the assay description and otherwise files the edge on the
    parent accession as parent_unresolved -- its docstring's stated purpose
    being that such an edge is "never mistaken for a measured one".
    Traversal was where it was mistaken for one.
    """
    g = Graph.load(resolution, skip_predicates=schema.held_out_predicates(),
                   symmetry=schema.symmetry())
    m1 = schema.metapaths["M1"]
    assert g.degree["UniProtKB:PARENT"] < g.degree["UniProtKB:PRO_CHAIN"]
    assert walk_metapath(g, "INCHIKEY:RESOLVED", m1["path"],
                         m1.get("constraints"), set())
    assert not walk_metapath(g, "INCHIKEY:VAGUE", m1["path"],
                             m1.get("constraints"), set())


def test_without_the_constraint_the_vaguer_edge_outscores_the_specific_one(
        resolution, schema):
    """The defect, kept measurable. Dropping the constraint restores the old
    behaviour, and the drug whose target was never pinned down scores HIGHER
    than the one whose target is known -- purely because the parent node
    carries fewer edges."""
    g = Graph.load(resolution, skip_predicates=schema.held_out_predicates(),
                   symmetry=schema.symmetry())
    path = schema.metapaths["M1"]["path"]
    vague = walk_metapath(g, "INCHIKEY:VAGUE", path, None, set())[V]
    resolved = walk_metapath(g, "INCHIKEY:RESOLVED", path, None, set())[V]
    assert vague > resolved


def test_the_gap_is_exactly_the_degree_ratio_raised_to_w(resolution, schema):
    """The mechanism, pinned as arithmetic rather than described in prose.

    Two drugs, one path each, through different single targets to the same
    gene and taxon: every factor cancels but the target's own degree, so the
    DWPC ratio is (deg_a / deg_b) ** w and nothing else. Measured on v0.1
    this predicted all three one-channel controls to within 0.04% --
    pomotrelvir, remdesivir and nirmatrelvir, ranks 17, 396 and 6,432, in
    inverse order of how well studied their targets are.
    """
    w = 0.4
    g = Graph.load(resolution, skip_predicates=schema.held_out_predicates(),
                   symmetry=schema.symmetry())
    path = schema.metapaths["M1"]["path"]
    vague = walk_metapath(g, "INCHIKEY:VAGUE", path, None, set(), w=w)[V]
    resolved = walk_metapath(g, "INCHIKEY:RESOLVED", path, None, set(), w=w)[V]
    ratio = (g.degree["UniProtKB:PRO_CHAIN"]
             / g.degree["UniProtKB:PARENT"]) ** w
    assert vague / resolved == pytest.approx(ratio, rel=1e-9)


def test_every_direct_acting_metapath_carries_the_constraint(schema):
    """M7 is M1 with a similarity hop in front and M8 is M1 across a fold
    boundary. A constraint on M1 alone would leave two unguarded doors into
    the same channel."""
    for name in ("M1", "M7", "M8"):
        c = schema.metapaths[name].get("constraints") or {}
        assert c.get("INHIBITS.target_resolution") == \
            "chain_from_assay_description", f"{name} is unguarded"


def test_an_inhibits_edge_without_target_resolution_is_refused(schema):
    """Why the qualifier is REQUIRED rather than merely declared: _qualifier_ok
    reads a missing qualifier as a failed match, so without this the omission
    would drop the edge from the direct-acting channel in silence."""
    e = _e("INCHIKEY:X", "INHIBITS", "UniProtKB:V", {"assay_type": "biochemical"})
    e["qualifiers"].pop("target_resolution")
    nodes = {"INCHIKEY:X": _drug("INCHIKEY:X"),
             "UniProtKB:V": _prot("UniProtKB:V", True)}
    codes = [v.code for v in schema.validate_edge(e, nodes)]
    assert "MISSING_QUALIFIER" in codes


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
def _row(lift, p, trustworthy=True):
    return {"lift": lift, "p_value": p, "trustworthy": trustworthy}


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


def test_an_untrustworthy_row_is_never_promoted_by_the_baseline_read():
    """HCoV-OC43 temporal computed_only: M4 pool 2,545, ONE positive, lift
    25.45, p=0.0393 -- and 100/2545 = 0.0393 exactly, so that p IS the
    sampling fraction for a single compound. The row already carried (!) for
    expectation < 1, and the first version of this branch read p, ignored the
    field beside it, and announced the result anyway.
    """
    from kgav.baselines import confound_warning

    per = {"M4": _row(25.45, 0.0393, trustworthy=False),
           "degree": _row(1.32, 0.5358, trustworthy=False)}
    assert confound_warning(per) is None

    # the same table with a trustworthy channel does report
    per["M4"]["trustworthy"] = True
    assert "BASELINE AT CHANCE" in confound_warning(per)


def test_a_row_with_no_trustworthy_verdict_is_not_promoted():
    """`is True`, not `is not False`: a row that never got a verdict has not
    earned a promotion either."""
    from kgav.baselines import confound_warning

    per = {"M4": {"lift": 25.45, "p_value": 0.0393},
           "degree": _row(1.32, 0.5358)}
    assert confound_warning(per) is None


def test_the_confounded_tail_marks_an_untrustworthy_challenger():
    """MERS-CoV temporal: M4 exceeds degree at 16.09 on five positives with
    expectation 0.19. Naming it without the flag read as a channel beating the
    baseline."""
    from kgav.baselines import confound_warning

    per = {"M4": _row(16.09, 0.0005, trustworthy=False),
           "degree": _row(4.62, 0.0260)}
    warn = confound_warning(per)
    assert "M4 exceeds it at 16.09" in warn
    assert "flagged (!)" in warn, warn
