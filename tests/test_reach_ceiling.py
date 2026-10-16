"""Gate for the reach-ceiling diagnostic.

It answers one question -- is a channel's ceiling set by missing data or by
a number someone chose -- and the way it goes wrong is by being read as a
recommendation to change the number. The other way is by quietly applying
limits to a channel that does not declare them, which would invent a cost
that is not there.
"""
from __future__ import annotations

import json
import sys
from importlib import import_module
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

R = import_module("diagnose_reach_ceiling")


def _release(tmp_path: Path, pathway_members: int = 3,
             hub_degree: int = 2) -> Path:
    """A release completing M4: SmallMolecule -TARGETS-> Protein
    -PARTICIPATES_IN-> Pathway -<PARTICIPATES_IN- Protein
    -HOST_FACTOR_FOR-> OrganismTaxon.
    """
    d = tmp_path / "rel"
    d.mkdir(parents=True, exist_ok=True)
    nodes = [{"id": "NCBITaxon:2697049", "class": "OrganismTaxon",
              "properties": {}},
             {"id": "REACT:BIG", "class": "Pathway", "properties": {}},
             {"id": "UniProtKB:HOST", "class": "Protein",
              "properties": {"taxon_id": "NCBITaxon:9606",
                             "is_viral": False}}]
    edges = [("UniProtKB:HOST", "PARTICIPATES_IN", "REACT:BIG", {}),
             ("UniProtKB:HOST", "HOST_FACTOR_FOR", "NCBITaxon:2697049",
              {"direction": "dependency"})]
    for i in range(pathway_members):
        t = f"UniProtKB:T{i}"
        nodes.append({"id": t, "class": "Protein",
                      "properties": {"taxon_id": "NCBITaxon:9606",
                                     "is_viral": False}})
        edges.append((t, "PARTICIPATES_IN", "REACT:BIG", {}))
        for j in range(hub_degree):
            c = f"INCHIKEY:C{i}_{j}"
            nodes.append({"id": c, "class": "SmallMolecule",
                          "properties": {}})
            edges.append((c, "TARGETS", t, {}))
    (d / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in nodes))
    (d / "edges.jsonl").write_text("\n".join(json.dumps(
        {"subject": s, "predicate": p, "object": o, "qualifiers": q,
         "primary_knowledge_source": "infores:string", "evidence_tier": 1,
         "first_asserted_date": "2021-01-01"}) for s, p, o, q in edges))
    return d


def _graph(rel: Path):
    from kgav.baselines import Graph
    from kgav.schema import load_schema
    s = load_schema()
    return Graph.load(rel, skip_predicates=s.held_out_predicates(),
                      symmetry=s.symmetry()), s


def test_survivors_counts_nodes_not_paths(tmp_path):
    """The ceiling is how many COMPOUNDS can arrive. A hop that multiplies
    paths without widening the frontier does not raise it, and counting
    paths would make a thin layer look wide."""
    g, s = _graph(_release(tmp_path, pathway_members=3, hub_degree=2))
    rows = R.survivors(g, s.metapaths["M4"]["path"],
                       s.metapaths["M4"].get("constraints"), set(),
                       list(g.drugs()))
    labels = [r[0] for r in rows]
    assert labels[0] == "SmallMolecule"
    assert labels[-1] == "compounds arriving"
    assert rows[0][1] == 6               # 3 targets x 2 compounds
    assert rows[-1][1] == 6              # all of them arrive
    # The pathway hop narrows to ONE pathway even though six compounds feed
    # it: that is the node-count behaviour being asserted.
    assert any(n == 1 for _label, n in rows[1:-1])


def test_the_constraint_on_the_last_hop_is_honoured(tmp_path):
    """M4 requires HOST_FACTOR_FOR.direction == dependency. Dropping the
    qualifier check would inflate every reach figure here."""
    rel = _release(tmp_path)
    edges = [json.loads(l) for l in
             (rel / "edges.jsonl").read_text().splitlines()]
    for e in edges:
        if e["predicate"] == "HOST_FACTOR_FOR":
            e["qualifiers"] = {"direction": "restriction"}
    (rel / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in edges))
    g, s = _graph(rel)
    rows = R.survivors(g, s.metapaths["M4"]["path"],
                       s.metapaths["M4"].get("constraints"), set(),
                       list(g.drugs()))
    assert rows[-1][1] == 0


def test_blocking_an_intermediate_lowers_reach_and_blocking_the_end_does_not(
        tmp_path):
    """Exclusion applies to INTERMEDIATE nodes only -- the endpoint is the
    answer, not a shortcut. walk_metapath says so; survivors must agree or
    the two disagree about what the channel does."""
    g, s = _graph(_release(tmp_path))
    mp = s.metapaths["M4"]
    free = R.survivors(g, mp["path"], mp.get("constraints"), set(),
                       list(g.drugs()))[-1][1]
    blocked_mid = R.survivors(g, mp["path"], mp.get("constraints"),
                              {"REACT:BIG"}, list(g.drugs()))[-1][1]
    blocked_end = R.survivors(g, mp["path"], mp.get("constraints"),
                              {"NCBITaxon:2697049"}, list(g.drugs()))[-1][1]
    assert free > 0
    assert blocked_mid == 0
    assert blocked_end == free


def test_reach_under_agrees_with_the_real_walker(tmp_path):
    """reach_under calls walk_metapath rather than reimplementing it, so
    this pins that the two stay the same function."""
    from kgav.baselines import walk_metapath
    g, s = _graph(_release(tmp_path))
    mp = s.metapaths["M4"]
    got = R.reach_under(g, mp, set(), "NCBITaxon:2697049")
    by_hand = sum(1 for d in g.drugs()
                  if "NCBITaxon:2697049" in walk_metapath(
                      g, d, mp["path"], mp.get("constraints"), set()))
    assert got == by_hand == 6


def test_a_channel_without_exclude_hubs_reports_that_the_limits_do_not_apply(
        tmp_path, capsys):
    """M2 does not declare exclude_hubs, so its counterfactuals are all the
    same number. Printing four identical figures without saying why would
    read as "the limits cost nothing here", which is a different claim."""
    R.main(["--release", str(_release(tmp_path)), "--channel", "M2"])
    out = capsys.readouterr().out
    assert "LIMITS DO NOT APPLY" in out
    assert "that is the finding" in out


def test_m4_declares_exclude_hubs_so_the_question_is_live():
    """If this ever flips, the script is measuring a cost that cannot
    exist and the near-miss cell needs a different diagnosis."""
    from kgav.schema import load_schema
    s = load_schema()
    assert s.metapaths["M4"].get("exclude_hubs") is True
    assert s.metapaths["M5"].get("exclude_hubs") is True
    assert not s.metapaths["M2"].get("exclude_hubs")


def test_the_measured_gap_matches_the_label_headroom_run():
    """The gap this script explains is pinned, not remembered. If
    label_headroom's figures move, this constant is stale and the verdict
    line it feeds would compare against a number nobody produced."""
    assert R.MEASURED_GAP[("NCBITaxon:2697049", "M4")] == (2982, 3869)
    assert R.MEASURED_GAP[("NCBITaxon:2697049", "M5")] == (2192, 6504)


def test_it_frames_the_counterfactuals_as_cost_not_advice(tmp_path, capsys):
    """Raising either limit changes every number in data/results. A
    diagnostic that reads as permission to do that is worse than none."""
    R.main(["--release", str(_release(tmp_path))])
    out = capsys.readouterr().out
    assert "COST, not a recommendation" in out
    assert "needs a registration" in out
    src = (ROOT / "scripts" / "diagnose_reach_ceiling.py").read_text()
    assert "Read-only. Writes nothing, changes no limit." in src


def test_an_unknown_channel_is_refused_with_the_list(tmp_path, capsys):
    assert R.main(["--release", str(_release(tmp_path)),
                   "--channel", "M99"]) == 2
    assert "no metapath M99" in capsys.readouterr().out


def test_a_missing_release_exits_two(tmp_path, capsys):
    assert R.main(["--release", str(tmp_path / "nope")]) == 2
    assert "no release at" in capsys.readouterr().out


def test_the_gap_verdict_is_withheld_on_another_release(tmp_path, capsys):
    """The pinned figures are v0.1's. Comparing a reach computed on another
    graph against them would be a verdict with one number from each --
    exactly the stale-before-column defect provenance.py exists for."""
    rel = _release(tmp_path)          # directory is named "rel", not v0.1
    R.main(["--release", str(rel)])
    out = capsys.readouterr().out
    assert "no gap verdict is given" in out
    assert "measured on v0.1" in out
    assert "does NOT close the gap" not in out


def test_the_gap_verdict_fires_on_the_release_it_was_measured_on(tmp_path,
                                                                 capsys):
    rel = _release(tmp_path)
    target = rel.parent / R.MEASURED_ON
    rel.rename(target)
    R.main(["--release", str(target)])
    out = capsys.readouterr().out
    assert "no gap verdict is given" not in out
    assert "does NOT close the gap" in out
    assert "it is the data" in out
