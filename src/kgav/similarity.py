"""Chemical similarity edges: the cold-start route.

WHY. 953 of 1,846 SARS-CoV-2 positives have no metapath at all, and 947 of
those have NO target edges whatsoever -- they carry a measured activity and
nothing else. No amount of better scoring reaches them: they are not in any
pool. CHEMICALLY_SIMILAR_TO is the only route that can, by letting a compound
borrow the target annotations of a structural neighbour. All 953 have a SMILES,
so the fingerprints are computable.

M6  Drug -CHEMICALLY_SIMILAR_TO-> Drug -TARGETS-> Protein -HOST_FACTOR_FOR-> Virus

SCOPE IS BIPARTITE ON PURPOSE. All-pairs over 12,821 compounds is 82M
comparisons, and most of the resulting edges would connect compounds that
already have paths -- adding degree without adding reach. Only
(path-less compound -> annotated compound) pairs are computed, which is both
cheaper and avoids flooding the graph with edges that serve no metapath.

WHAT THIS CANNOT FIX. Similarity inherits a neighbour's annotations wholesale,
including its artifacts. A structural analogue of a promiscuous kinase
inhibitor inherits the promiscuity, and an analogue of chloroquine inherits the
SIGMAR1 route that made chloroquine look promising in Vero E6. So M6 raises
COVERAGE, which is not the same as accuracy, and the check that matters after
adding it is whether the named controls still behave.
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path

from kgav.normalize_chem import skeleton_of_id

DEFAULT_THRESHOLD = 0.7
FP_RADIUS = 2
FP_BITS = 2048


def fingerprints(smiles_by_id: dict[str, str]) -> tuple[dict[str, object], Counter]:
    """ECFP4 (Morgan radius 2, 2048 bits) per node id.

    Unparseable structures are counted, not silently skipped: a compound that
    cannot be fingerprinted is unreachable by any method and that number
    belongs in the report.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFingerprintGenerator

    RDLogger.DisableLog("rdApp.*")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_BITS)

    out: dict[str, object] = {}
    stats: Counter = Counter()
    for nid, smi in smiles_by_id.items():
        if not smi:
            stats["no_smiles"] += 1
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            stats["unparseable"] += 1
            continue
        out[nid] = gen.GetFingerprint(mol)
        stats["fingerprinted"] += 1
    return out, stats


def similarity_edges(query_fps: dict[str, object], ref_fps: dict[str, object],
                     threshold: float = DEFAULT_THRESHOLD,
                     top_k: int | None = 10,
                     stats: Counter | None = None
                     ) -> list[tuple[str, str, float]]:
    """(query, reference, tanimoto) above threshold.

    ONE COMPOUND IS NOT ITS OWN NEIGHBOUR. This excluded a pair only when
    the two node ids were equal, and the graph uses the FULL InChIKey as a
    node id so that a salt and its parent, or two stereoisomers, remain
    separate nodes. 434 of the 3,384 CHEMICALLY_SIMILAR_TO edges in v0.1
    join two keys of ONE molecule -- same 14-character skeleton, different
    suffix -- almost all at tanimoto 1.0, because they are the same
    structure.

    That is not a near-miss. normalize_chem.same_as_rule states that a
    shared skeleton means one compound, and gtop.py uses those same 14
    characters as PROOF that two keys are one compound when joining
    ligands to the graph. Only this function read id inequality as proof
    they were two. A metapath stepping M6 (Drug -CHEM_SIMILAR-> Drug
    -TARGETS-> ...) across such an edge has not moved to another compound;
    it has walked in a circle and scored the result as independent
    evidence.

    `stats`, if given, is filled with the counters this makes -- notably
    same_skeleton_excluded, which was 434 and silent.

    top_k caps how many neighbours a query keeps. Without it a compound in a
    dense scaffold series acquires hundreds of near-identical neighbours, and
    DWPC would read that as overwhelming evidence when it is one scaffold
    counted many times.
    """
    from rdkit import DataStructs

    ref_ids = list(ref_fps)
    ref_list = [ref_fps[r] for r in ref_ids]
    if not ref_list:
        return []

    ref_skel = [skeleton_of_id(r) for r in ref_ids]

    edges: list[tuple[str, str, float]] = []
    for qid, qfp in query_fps.items():
        qskel = skeleton_of_id(qid)
        sims = DataStructs.BulkTanimotoSimilarity(qfp, ref_list)
        hits = []
        for i, sim in enumerate(sims):
            if sim < threshold or ref_ids[i] == qid:
                continue
            # BOTH must have a skeleton for a match to mean anything. A
            # minted KGAV: id yields None, and None == None must not read
            # as "same compound".
            if qskel is not None and ref_skel[i] == qskel:
                if stats is not None:
                    stats["same_skeleton_excluded"] += 1
                continue
            hits.append((ref_ids[i], sim))
        # TIE-BREAK ON THE ID, AND IT IS LOAD-BEARING. This sorted on
        # similarity alone. Python's sort is stable, so tied
        # neighbours kept the order of ref_ids -- which is
        # list(ref_fps), the dict's insertion order, built by
        # ingest_similarity.py from a SET. Set iteration order for
        # strings varies per process, so when a tie straddled position
        # top_k, which neighbour survived the cut depended on the
        # hash seed.
        #
        # Measured: two builds of v0.1 from identical raw inputs
        # differed by exactly 8 CHEMICALLY_SIMILAR_TO edges out of
        # 389,306, with the total unchanged -- eight swapped for
        # eight. That made edges_sha unreproducible, and nothing could
        # say whether a rebuilt release was a different graph.
        #
        # baselines.py documents this same trap for rank() and
        # deliberately does NOT break the tie there, reporting the
        # range over tie orderings instead, because inventing an order
        # would invent a ranking. Here it is the opposite: a top_k cut
        # MUST choose, keeping every tie would change the edge count,
        # and so the only question is whether the choice is declared
        # or accidental. The id is arbitrary and says nothing about
        # which neighbour is better -- it is chosen because it is a
        # property of the content, so the same graph yields the same
        # edges anywhere.
        hits.sort(key=lambda kv: (-kv[1], kv[0]))
        if top_k:
            # THE CUT IS COUNTED (fix-list E7, and E8 asks for exactly
            # this). A neighbour above the threshold that falls outside
            # top_k is a real similarity the graph does not carry, and
            # discarding it silently means the cap's cost is unknown --
            # the docstring argues for the cap on principle and had no
            # number behind it.
            if stats is not None:
                over = len(hits) - top_k
                if over > 0:
                    stats["above_threshold_beyond_top_k"] += over
                    stats["queries_truncated_by_top_k"] += 1
            hits = hits[:top_k]
        edges.extend((qid, rid, float(s)) for rid, s in hits)
    return edges


def similarity_distribution(query_fps: dict[str, object],
                            ref_fps: dict[str, object],
                            bins: tuple[float, ...] = (0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
                            ) -> Counter:
    """How many query compounds have at least one neighbour at each threshold.

    Run before choosing a cut: 0.7 is conventional but arbitrary, and the right
    value depends on how much reach it buys.
    """
    from rdkit import DataStructs

    ref_list = list(ref_fps.values())
    if not ref_list:
        return Counter()
    out: Counter = Counter()
    for qfp in query_fps.values():
        best = max(DataStructs.BulkTanimotoSimilarity(qfp, ref_list), default=0.0)
        for b in bins:
            if best >= b:
                out[b] += 1
    return out


def build(em, smiles_by_id: dict[str, str], query_ids: set[str], ref_ids: set[str],
          source: str, threshold: float = DEFAULT_THRESHOLD,
          top_k: int | None = 10, date: str = "1970-01-01") -> Counter:
    """Emit CHEMICALLY_SIMILAR_TO edges from query compounds to annotated ones.

    evidence_tier 3: computed, not asserted by any source. The date is 1970
    rather than today's, deliberately -- a computed edge has no assertion date,
    and stamping it with the build date would make it look newer than every
    real fact and corrupt any temporal split.

    SORTED HERE, WHERE THE SETS ARE ACTUALLY READ (fix-list P15). Both
    arguments are SETS, and string hashing is randomised per process, so
    iterating them carried the hash seed into the emission order: the layer
    wrote the same 3,154 edges in a different order on every build, and
    `v0.1-similarity`'s edges_sha moved every time -- fdb80f57, 262cfc5e,
    6cd90b09 over three rebuilds, which made the layer stamp unable to tell
    a real change from the noise it exists to filter.

    ingest_similarity.py sorts too, and its comment said that sort was kept
    "because the emission order should not depend on the hash seed either".
    It could not do that: those sorted dicts feed the survey printout, and
    this function was handed the raw sets and re-fingerprinted them. A
    protection one file away from the code it protects is the defect this
    module's own history is made of, so the sort is now where the iteration
    happens.

    The EDGE SET never depended on it -- similarity_edges breaks the top-k
    tie on the id, which is what closed the 8-edge discrepancy -- so this
    changes bytes and ordering, not content, and the test asserts exactly
    that.
    """
    q_fps, q_stats = fingerprints({i: smiles_by_id.get(i, "")
                                   for i in sorted(query_ids)})
    r_fps, r_stats = fingerprints({i: smiles_by_id.get(i, "")
                                   for i in sorted(ref_ids)})

    stats: Counter = Counter()
    stats["query_fingerprinted"] = q_stats["fingerprinted"]
    stats["ref_fingerprinted"] = r_stats["fingerprinted"]
    stats["unparseable"] = q_stats["unparseable"] + r_stats["unparseable"]
    stats["no_smiles"] = q_stats["no_smiles"] + r_stats["no_smiles"]

    for qid, rid, sim in similarity_edges(q_fps, r_fps, threshold, top_k,
                                          stats=stats):
        em.edge(qid, "CHEMICALLY_SIMILAR_TO", rid,
                source=source, date=date, tier=3,
                quals={"tanimoto_ecfp4": round(sim, 4)})
        stats["edges"] += 1
        stats["queries_connected"] = len(
            {e["subject"] for e in em.edges if e["predicate"] == "CHEMICALLY_SIMILAR_TO"})
    return stats


def load_compounds(release: Path) -> dict[str, str]:
    """node id -> SMILES for every SmallMolecule in a release."""
    import json
    out: dict[str, str] = {}
    for line in (Path(release) / "nodes.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        n = json.loads(line)
        if n["class"] == "SmallMolecule":
            out[n["id"]] = n["properties"].get("smiles") or ""
    return out
