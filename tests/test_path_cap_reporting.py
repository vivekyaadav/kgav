"""A truncated walk says so.

walk_metapath stops expanding the frontier once it exceeds max_paths
(20,000) and returns what it has. Which paths it saw by then depends
on the order g.neighbours yields them, which is the order
edges.jsonl was written in. Until 2026-10-29 that happened silently:
a truncated score was indistinguishable from a complete one in the
return value and in every results file.

MEASURED ON v0.1. 114 of M4's 2,983 reaching drugs truncate; no other
metapath truncates at all. Shuffling edges.jsonl and recomputing
moved M4's Hits@100 for SARS-CoV-2 from 2 to 1 and its lift from 0.99
to 0.50, while M5 and M7 did not change by a bit.

hits_determined cannot see this. It detects EXACT ties at the cut,
and two differently-truncated compounds are not tied -- they were
computed from different sets of paths. So this is a separate flag,
not a widening of that one.

These tests do NOT assert that a truncated channel is excluded from
anything. Whether truncation should make a channel untrustworthy is a
decision about what the evaluation claims; reporting it is not.
"""
from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.baselines import (  # noqa: E402
    Graph,
    dwpc_scores,
    walk_metapath,
    walk_metapath_detail,
)
from kgav.schema import load_schema  # noqa: E402

PATH = ["SmallMolecule", "TARGETS", "Protein"]
# TWO HOPS, because the cap's `break` leaves the loop over the
# frontier -- with a single-node frontier it skips nothing and the
# score is complete. A one-hop fixture reports the cap as reached and
# loses no paths, which is exactly the distinction the flag has to
# make.
PATH2 = ["SmallMolecule", "TARGETS", "Protein",
         "PARTICIPATES_IN", "Pathway"]


def _graph2(tmp_path: Path, n_prot: int, n_path: int) -> Graph:
    """n proteins, each in n_path pathways, so the frontier after hop
    one is n_prot and expanding it yields n_prot * n_path."""
    import json
    rel = tmp_path / "rel2"
    rel.mkdir()
    nodes = [{"id": "INCHIKEY:D", "class": "SmallMolecule",
              "properties": {"smiles": "C", "inchikey_skel": "DD",
                             "is_approved": True, "salt_collapsed": False,
                             "stereo_collapsed": False}}]
    edges = []
    for i in range(n_prot):
        pid = f"UniProtKB:P{i:04d}"
        nodes.append({"id": pid, "class": "Protein",
                      "properties": {"taxon_id": "NCBITaxon:9606",
                                     "is_viral": False, "sequence_hash": "h",
                                     "reviewed": True}})
        edges.append({"subject": "INCHIKEY:D", "predicate": "TARGETS",
                      "object": pid, "qualifiers": {"direction": "unknown"},
                      "evidence_tier": 1, "first_asserted_date": "2020-01-01",
                      "primary_knowledge_source": "infores:chembl"})
        for j in range(n_path):
            wid = f"REACT:R{i:03d}{j:03d}"
            nodes.append({"id": wid, "class": "Pathway",
                          "properties": {"label": wid, "source": "reactome",
                                         "n_members": 5}})
            edges.append({"subject": pid, "predicate": "PARTICIPATES_IN",
                          "object": wid, "qualifiers": {},
                          "evidence_tier": 1,
                          "first_asserted_date": "2020-01-01",
                          "primary_knowledge_source": "infores:reactome"})
    (rel / "nodes.jsonl").write_text("\n".join(json.dumps(x) for x in nodes))
    (rel / "edges.jsonl").write_text("\n".join(json.dumps(x) for x in edges))
    return Graph.load(rel, symmetry=load_schema().symmetry())


def _graph(tmp_path: Path, n_targets: int) -> Graph:
    """One drug wired to n host proteins, so the frontier after the
    first hop is exactly n and the cap can be placed either side."""
    import json
    rel = tmp_path / "rel"
    rel.mkdir()
    nodes = [{"id": "INCHIKEY:D", "class": "SmallMolecule",
              "properties": {"smiles": "C", "inchikey_skel": "DD",
                             "is_approved": True, "salt_collapsed": False,
                             "stereo_collapsed": False}}]
    edges = []
    for i in range(n_targets):
        pid = f"UniProtKB:P{i:04d}"
        nodes.append({"id": pid, "class": "Protein",
                      "properties": {"taxon_id": "NCBITaxon:9606",
                                     "is_viral": False, "sequence_hash": "h",
                                     "reviewed": True}})
        edges.append({"subject": "INCHIKEY:D", "predicate": "TARGETS",
                      "object": pid, "qualifiers": {"direction": "unknown"},
                      "evidence_tier": 1, "first_asserted_date": "2020-01-01",
                      "primary_knowledge_source": "infores:chembl"})
    (rel / "nodes.jsonl").write_text("\n".join(json.dumps(x) for x in nodes))
    (rel / "edges.jsonl").write_text("\n".join(json.dumps(x) for x in edges))
    return Graph.load(rel, symmetry=load_schema().symmetry())


def test_a_complete_walk_reports_not_truncated(tmp_path):
    g = _graph(tmp_path, 10)
    scores, cut = walk_metapath_detail(g, "INCHIKEY:D", PATH, None, set(),
                                       max_paths=1000)
    assert len(scores) == 10
    assert cut is False


def test_a_capped_walk_reports_truncated(tmp_path):
    """The cap is placed BELOW the frontier size on purpose. A fixture
    too small to reach the cap would pass whatever the flag did --
    which is how this went unreported for as long as it did."""
    g = _graph2(tmp_path, 20, 10)
    full, uncut = walk_metapath_detail(g, "INCHIKEY:D", PATH2, None, set(),
                                       max_paths=10_000)
    assert uncut is False and len(full) == 200
    scores, cut = walk_metapath_detail(g, "INCHIKEY:D", PATH2, None, set(),
                                       max_paths=15)
    assert cut is True, "the cap fired and the walk did not say so"
    assert len(scores) < len(full), (
        "nothing was actually dropped; this fixture cannot detect "
        "truncation")


def test_reaching_the_cap_on_the_last_frontier_entry_is_not_truncation(
        tmp_path):
    """The break leaves the loop over the frontier. Hitting the cap
    while expanding its LAST entry drops nothing, and a warning
    printed on a complete figure is the same as no warning.

    One protein, fifty pathways: the cap fires, the frontier had one
    entry, every pathway is still returned."""
    g = _graph2(tmp_path, 1, 50)
    scores, cut = walk_metapath_detail(g, "INCHIKEY:D", PATH2, None, set(),
                                       max_paths=5)
    assert len(scores) == 50, "paths were dropped; the premise is wrong"
    assert cut is False, (
        "the cap was reached but nothing was lost, and it reported "
        "truncation anyway")


def test_the_plain_wrapper_still_returns_only_scores(tmp_path):
    """Nine callers read walk_metapath positionally. Returning the
    tuple from it would break every one of them."""
    g = _graph(tmp_path, 10)
    out = walk_metapath(g, "INCHIKEY:D", PATH, None, set())
    assert isinstance(out, dict)
    assert len(out) == 10


def test_dwpc_scores_counts_truncated_drugs(tmp_path):
    g = _graph2(tmp_path, 20, 10)
    mp = {"MX": {"path": PATH2}}
    seen: Counter = Counter()
    dwpc_scores(g, mp, {}, truncated=seen)
    assert seen == Counter(), "nothing should truncate at the default cap"

    import kgav.baselines as B
    original = B.walk_metapath_detail

    def tiny(*a, **k):
        k["max_paths"] = 15
        return original(*a, **k)

    B.walk_metapath_detail = tiny
    try:
        seen2: Counter = Counter()
        dwpc_scores(g, mp, {}, truncated=seen2)
    finally:
        B.walk_metapath_detail = original
    assert seen2["MX"] == 1, seen2


def test_the_counter_is_optional(tmp_path):
    """Omitted, the signature behaves exactly as the nine existing
    callers expect."""
    g = _graph(tmp_path, 10)
    out = dwpc_scores(g, {"MX": {"path": PATH}}, {})
    assert "MX" in out and "INCHIKEY:D" in out["MX"]


@pytest.mark.parametrize("needle", [
    'truncated=truncated',
    'per[name]["truncated_walks"]',
    'per[name]["order_dependent"]',
    '~capped',
    'PATH CAP REACHED',
])
def test_compare_metapaths_records_and_prints_it(needle):
    """Pinned at the call site. A flag the publishing script never
    asks for is a flag nobody sees, which is the state this replaces."""
    src = (ROOT / "scripts" / "compare_metapaths.py").read_text()
    assert needle in src, needle


# ----------------------------------------------- the cap is a measured number
def test_the_cap_clears_the_largest_frontier_this_graph_produces():
    """20,000 had no recorded reason, and on v0.1 M4's largest
    frontier is 230,176 -- 11.5x over, truncating 114 of its 2,983
    reaching drugs. M3 (16,726) and M5 (15,281) cleared it by 16% and
    24%, which is luck rather than design: a layer touching those
    paths would have started truncating them silently.

    Pinned so the constant cannot drift back under the measurement
    without this failing. The figures are in the constant's comment,
    where a reader meets them."""
    from kgav.baselines import MAX_PATHS
    measured_peak = 230_176          # M4, v0.1, 2026-10-30
    assert MAX_PATHS > measured_peak, (
        f"MAX_PATHS {MAX_PATHS:,} no longer clears the largest frontier "
        f"measured on v0.1 ({measured_peak:,}); M4 would truncate again")
    assert MAX_PATHS >= 4 * measured_peak, (
        "the headroom is for layers not yet built -- gtop and curated "
        "both add edges on M4's path")


def test_the_cap_is_still_a_cap():
    """Not removed. A pathological graph should stop rather than
    exhaust memory, and the truncation flag exists so that a release
    which reaches the cap says so instead of returning one draw."""
    import inspect

    from kgav.baselines import MAX_PATHS, walk_metapath_detail
    assert MAX_PATHS < 10 ** 9, "the cap was effectively removed"
    src = inspect.getsource(walk_metapath_detail)
    assert "if len(nxt) > max_paths:" in src
    assert "truncated = True" in src


def test_the_measured_peaks_are_recorded_where_the_constant_is():
    """A number chosen from a measurement is only defensible while the
    measurement is written down beside it."""
    src = (ROOT / "src" / "kgav" / "baselines.py").read_text()
    head = src.split("MAX_PATHS = ")[0][-2000:]
    assert "230,176" in head, "M4's measured peak is not recorded"
    assert "16,726" in head and "15,281" in head, (
        "M3 and M5's near-misses are not recorded, and they are the "
        "reason the old cap was luck rather than design")
    assert "77 seconds" in head, "the cost of a complete traversal is "\
        "not recorded, and it is why the cap was raised"


def test_both_entry_points_actually_default_to_the_measured_cap():
    """PINNING THE CONSTANT IS NOT PINNING ITS USE. Two mutants
    survived the tests above by reverting either signature to the
    literal 20000 while MAX_PATHS stayed at 1,000,000 -- the constant
    correct, documented, measured, and not read by anything.

    That is the defect this suite has found in its own guards more
    than once: keyed on a line the change does not have to touch.
    """
    import inspect

    from kgav.baselines import MAX_PATHS, walk_metapath, walk_metapath_detail
    for fn in (walk_metapath, walk_metapath_detail):
        got = inspect.signature(fn).parameters["max_paths"].default
        assert got == MAX_PATHS, (
            f"{fn.__name__} defaults to {got!r}, not MAX_PATHS "
            f"({MAX_PATHS:,})")


def test_a_frontier_over_the_old_cap_is_not_truncated_by_default(tmp_path):
    """Behavioural, not structural. 30,000 clears the old cap of
    20,000 and sits far under the new one, so a default call must
    walk it whole. This is the case that was silently truncated
    before 2026-10-30."""
    g = _graph2(tmp_path, 300, 100)
    scores, cut = walk_metapath_detail(g, "INCHIKEY:D", PATH2, None, set())
    assert cut is False, "the default cap truncated a 30,000 frontier"
    assert len(scores) == 30_000, len(scores)
