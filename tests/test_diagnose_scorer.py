"""Gate for the scorer comparison.

The measurement it reports decides whether the default scoring should change,
so the two things that must not be wrong are the statistic and the silence:
a rank correlation computed without tie correction would be meaningless on
the channel that puts 80.8% of its pool in one block, and a run that reaches
nothing must say so rather than printing a header and two paragraphs.
"""
from __future__ import annotations

import json
import sys
from importlib import import_module
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

D = import_module("diagnose_scorer")


# ------------------------------------------------------------- the statistic
def test_ranks_average_within_a_tie_block():
    """1..n in sort order would separate equal scores by a whole rank. On M1
    that is most of the pool."""
    assert D.ranks([10.0, 20.0, 30.0]) == [1.0, 2.0, 3.0]
    assert D.ranks([5.0, 5.0, 9.0]) == [1.5, 1.5, 3.0]
    assert D.ranks([7.0, 7.0, 7.0]) == [2.0, 2.0, 2.0]
    assert D.ranks([2.0, 1.0, 2.0, 1.0]) == [3.5, 1.5, 3.5, 1.5]


@pytest.mark.parametrize("xs,ys,want", [
    # Closed forms, not values copied from a run. The tie case was also
    # checked against scipy.stats.spearmanr over 400 random tie-heavy
    # vectors with zero mismatches; scipy is NOT a dependency (pyproject:
    # pyyaml, pandas, pyarrow), which is why these are pinned rather than
    # imported.
    ([1, 2, 3], [1, 2, 3], 1.0),
    ([1, 2, 3], [3, 2, 1], -1.0),
    ([1, 2, 2, 3], [1, 1, 2, 2], 0.7071067811865476),
])
def test_spearman_matches_the_closed_form(xs, ys, want):
    assert D.spearman(xs, ys) == pytest.approx(want, abs=1e-12)


def test_a_constant_vector_has_no_correlation_rather_than_zero():
    """0.00 would read as "measured, and unrelated". A channel whose every
    drug scores the same is the degenerate case this project actually has."""
    assert D.spearman([1.0, 1.0, 1.0], [1.0, 2.0, 3.0]) is None
    assert D.spearman([1.0], [2.0]) is None


def test_spearman_refuses_unpaired_vectors():
    with pytest.raises(ValueError):
        D.spearman([1.0, 2.0], [1.0])


def test_degeneracy_reports_distinct_largest_and_share_in_any_tie():
    assert D.degeneracy([1.0, 1.0, 1.0, 2.0]) == (2, 3, 75.0)
    assert D.degeneracy([1.0, 2.0, 3.0]) == (3, 1, 0.0)
    assert D.degeneracy([]) == (0, 0, 0.0)


# ------------------------------------------------------------------ end to end
def _release(tmp_path: Path, extra_targets: int = 0) -> Path:
    """A release that completes M1: SmallMolecule -INHIBITS-> Protein
    -ENCODED_BY-> Gene -BELONGS_TO-> OrganismTaxon.

    `extra_targets` gives compound B that many additional routes, so the two
    compounds differ in degree and nothing else.
    """
    d = tmp_path / "rel"
    d.mkdir(parents=True, exist_ok=True)
    nodes = [{"id": "NCBITaxon:2697049", "class": "OrganismTaxon",
              "properties": {}}]
    edges = []

    def protein(i):
        nodes.append({"id": f"UniProtKB:P{i}", "class": "Protein",
                      "properties": {"taxon_id": "NCBITaxon:2697049",
                                     "is_viral": True}})
        nodes.append({"id": f"GENE:{i}", "class": "Gene", "properties": {}})
        edges.append((f"UniProtKB:P{i}", "ENCODED_BY", f"GENE:{i}", {}))
        edges.append((f"GENE:{i}", "BELONGS_TO", "NCBITaxon:2697049", {}))

    for i in range(1 + extra_targets):
        protein(i)
    for c in ("INCHIKEY:AAA", "INCHIKEY:BBB"):
        nodes.append({"id": c, "class": "SmallMolecule", "properties": {}})
    q = {"target_resolution": "chain_from_assay_description"}
    edges.append(("INCHIKEY:AAA", "INHIBITS", "UniProtKB:P0", q))
    for i in range(1 + extra_targets):
        edges.append(("INCHIKEY:BBB", "INHIBITS", f"UniProtKB:P{i}", q))

    (d / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in nodes))
    (d / "edges.jsonl").write_text("\n".join(json.dumps(
        {"subject": s, "predicate": p, "object": o, "qualifiers": qu,
         "primary_knowledge_source": "infores:chembl", "evidence_tier": 1,
         "first_asserted_date": "2021-01-01"}) for s, p, o, qu in edges))
    return d


def test_it_reports_a_populated_table(tmp_path, capsys):
    """The empty cases below are easy to pass and prove nothing on their own.
    This asserts the table fills in when a metapath actually completes."""
    rel = _release(tmp_path, extra_targets=4)
    assert D.main(["--release", str(rel)]) == 0
    out = capsys.readouterr().out
    assert "NO VIRUS REACHED" not in out
    assert "=== NCBITaxon:2697049" in out
    m1 = next(l for l in out.splitlines() if l.strip().startswith("M1"))
    assert "no drug reached" not in m1
    assert "2" in m1.split()[1]          # both compounds reached


def test_the_legacy_scoring_tracks_degree_more_closely_than_pdp(tmp_path):
    """THE MEASUREMENT. Two compounds identical but for edge count: under the
    legacy form the gap is the degree gap, under pdp it is damped. Asserted
    on the scores rather than the printed table so a formatting change cannot
    silently turn this green."""
    from kgav.baselines import (SCORING_LEGACY, SCORING_PDP, Graph,
                                dwpc_scores)
    from kgav.schema import load_schema
    schema = load_schema()
    g = Graph.load(_release(tmp_path, extra_targets=9),
                   skip_predicates=schema.held_out_predicates(),
                   symmetry=schema.symmetry())
    got = {}
    for scoring in (SCORING_LEGACY, SCORING_PDP):
        per = dwpc_scores(g, schema.metapaths, schema.retrieval_limits,
                          scoring=scoring)["M1"]
        got[scoring] = {d: h["NCBITaxon:2697049"] for d, h in per.items()}
    for scoring in (SCORING_LEGACY, SCORING_PDP):
        assert set(got[scoring]) == {"INCHIKEY:AAA", "INCHIKEY:BBB"}
    ratio = {s: got[s]["INCHIKEY:BBB"] / got[s]["INCHIKEY:AAA"] for s in got}
    assert ratio[SCORING_LEGACY] > ratio[SCORING_PDP], (
        f"pdp must damp the compound's own degree: {ratio}")


def test_a_run_that_reaches_nothing_says_so(tmp_path, capsys):
    """It printed a header, two trailing paragraphs and not one table, which
    reads as "nothing to flag" rather than "nothing was measured"."""
    rel = tmp_path / "bare"
    rel.mkdir()
    (rel / "nodes.jsonl").write_text(json.dumps(
        {"id": "INCHIKEY:AAA", "class": "SmallMolecule", "properties": {}}))
    (rel / "edges.jsonl").write_text("")
    assert D.main(["--release", str(rel), "--all-viruses"]) == 1
    out = capsys.readouterr().out
    assert "NO VIRUS REACHED" in out
    assert "0 completed walks" in out


def test_an_unreached_virus_is_distinguished_from_a_null(tmp_path, capsys):
    """An empty table for a virus nothing reaches is not a negative result,
    and the two must not look alike."""
    rel = _release(tmp_path)
    D.main(["--release", str(rel), "--virus", "NCBITaxon:64320"])
    out = capsys.readouterr().out
    assert "no metapath reaches NCBITaxon:64320" in out
    assert "empty by construction, not by a null result" in out


def test_a_missing_release_exits_two(tmp_path, capsys):
    assert D.main(["--release", str(tmp_path / "nope")]) == 2
    assert "no release at" in capsys.readouterr().out


# ------------------------------------------------------------- what it is not
def test_it_writes_nothing_and_changes_no_default():
    src = (ROOT / "scripts" / "diagnose_scorer.py").read_text()
    assert "write_results" not in src
    assert "Read-only. Writes nothing" in src
    assert "is a separate registered decision" in src.replace("\n", " ")


def test_the_default_scoring_is_still_the_published_results_one():
    """If someone flips the default, every number in data/results silently
    describes a scorer that produced none of them."""
    import inspect

    from kgav.baselines import SCORING_LEGACY, dwpc_scores, walk_metapath
    for fn in (walk_metapath, dwpc_scores):
        assert inspect.signature(fn).parameters["scoring"].default == \
            SCORING_LEGACY, fn.__name__
