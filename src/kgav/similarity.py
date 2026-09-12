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
                     top_k: int | None = 10) -> list[tuple[str, str, float]]:
    """(query, reference, tanimoto) above threshold.

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

    edges: list[tuple[str, str, float]] = []
    for qid, qfp in query_fps.items():
        sims = DataStructs.BulkTanimotoSimilarity(qfp, ref_list)
        hits = [(ref_ids[i], s) for i, s in enumerate(sims)
                if s >= threshold and ref_ids[i] != qid]
        hits.sort(key=lambda kv: -kv[1])
        if top_k:
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
    """
    q_fps, q_stats = fingerprints({i: smiles_by_id.get(i, "") for i in query_ids})
    r_fps, r_stats = fingerprints({i: smiles_by_id.get(i, "") for i in ref_ids})

    stats: Counter = Counter()
    stats["query_fingerprinted"] = q_stats["fingerprinted"]
    stats["ref_fingerprinted"] = r_stats["fingerprinted"]
    stats["unparseable"] = q_stats["unparseable"] + r_stats["unparseable"]
    stats["no_smiles"] = q_stats["no_smiles"] + r_stats["no_smiles"]

    for qid, rid, sim in similarity_edges(q_fps, r_fps, threshold, top_k):
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
