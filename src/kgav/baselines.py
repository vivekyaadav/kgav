"""Three non-learned baselines. Build these before any model.

    1. degree ranking      the null hypothesis
    2. DWPC over M1-M6     degree-weighted path count, Hetionet-style
    3. network proximity   drug targets to host factors in the interactome

A learned scorer that does not clearly beat all three is not a result. Most
published knowledge-graph repurposing work never runs this comparison, which
is why so much of it is unreliable.

HELD-OUT PREDICATES ARE NOT TRAVERSED. HAS_ANTIVIRAL_ACTIVITY_AGAINST,
TREATS and IN_TRIAL_FOR are evaluation labels. They are also the only edges
connecting most compounds directly to a virus, so excluding them is what makes
this a prediction rather than a lookup: the question is whether mechanistic
paths recover the cell-based antiviral activity we deliberately withheld.

DWPC WEIGHTING. Each path contributes the product of its nodes' degrees raised
to -w (w=0.4, as in Hetionet). A route through a promiscuous protein counts for
less than one through a specific target. Without this, path counts simply
measure how well-studied the intermediates are.
"""
from __future__ import annotations

import json
import math
import random
from collections import Counter, defaultdict
from pathlib import Path

DEFAULT_W = 0.4


class Graph:
    """Adjacency built for metapath walking.

    out[(predicate, node)] -> [(neighbour, edge)] in the schema's declared
    direction; inv[...] is the reverse, needed because M4 traverses
    PARTICIPATES_IN forwards then backwards.
    """

    def __init__(self) -> None:
        self.nodes: dict[str, dict] = {}
        self.out: dict[tuple[str, str], list[tuple[str, dict]]] = defaultdict(list)
        self.inv: dict[tuple[str, str], list[tuple[str, dict]]] = defaultdict(list)
        self.degree: Counter = Counter()

    @classmethod
    def load(cls, release: Path, skip_predicates: set[str] | None = None) -> Graph:
        skip = skip_predicates or set()
        g = cls()
        for line in (Path(release) / "nodes.jsonl").read_text().splitlines():
            if line.strip():
                n = json.loads(line)
                g.nodes[n["id"]] = n
        for line in (Path(release) / "edges.jsonl").read_text().splitlines():
            if not line.strip():
                continue
            e = json.loads(line)
            p = e["predicate"]
            # Degree counts ALL edges, including held-out ones: degree is a
            # property of the graph, and excluding label edges from it would
            # make the degree baseline easier to beat than it really is.
            g.degree[e["subject"]] += 1
            g.degree[e["object"]] += 1
            if p in skip:
                continue
            g.out[(p, e["subject"])].append((e["object"], e))
            g.inv[(p, e["object"])].append((e["subject"], e))
        return g

    def klass(self, nid: str) -> str | None:
        n = self.nodes.get(nid)
        return n["class"] if n else None

    def is_viral(self, nid: str) -> bool:
        n = self.nodes.get(nid)
        return bool(n and n["properties"].get("is_viral"))

    def drugs(self) -> list[str]:
        return [n for n, d in self.nodes.items() if d["class"] == "SmallMolecule"]

    def viruses(self) -> list[str]:
        return [n for n, d in self.nodes.items() if d["class"] == "OrganismTaxon"]


def _qualifier_ok(edge: dict, predicate: str, constraints: dict | None) -> bool:
    if not constraints:
        return True
    for key, want in constraints.items():
        pred, _, qual = key.partition(".")
        if pred != predicate:
            continue
        if (edge.get("qualifiers") or {}).get(qual) != want:
            return False
    return True


def _hub_set(g: Graph, limits: dict) -> set[str]:
    """Class-scoped hub exclusion, per the schema.

    A global cut would exclude the virus nodes and the viral replicase, which
    is where every metapath ENDS.
    """
    cut = limits.get("hub_degree_cut", 460)
    classes = set(limits.get("hub_exclude_classes") or ["Protein"])
    include_viral = bool(limits.get("hub_exclude_viral", False))
    return {n for n, d in g.degree.items()
            if d > cut and g.klass(n) in classes
            and (include_viral or not g.is_viral(n))}


def _oversized_pathways(g: Graph, limits: dict) -> set[str]:
    """Pathways too general to be a mechanism.

    A GO term with 1,389 members generates ~1.9M protein pairs in M4, so every
    drug touching any transcription regulator looks adjacent to every host
    factor.
    """
    cap = limits.get("max_pathway_members")
    if not cap:
        return set()
    members: Counter = Counter()
    for (pred, _node), targets in g.out.items():
        if pred == "PARTICIPATES_IN":
            for t, _e in targets:
                members[t] += 1
    return {p for p, n in members.items() if n > cap}


def walk_metapath(g: Graph, start: str, path: list, constraints: dict | None,
                  blocked: set[str], w: float = DEFAULT_W,
                  max_paths: int = 20000) -> dict[str, float]:
    """DWPC from `start` along a metapath, returning {endpoint: score}.

    Each frontier entry carries the accumulated degree-product weight. Nodes
    already on the path are not revisited: a path that loops back through its
    own intermediate is not evidence.
    """
    frontier: list[tuple[str, float, frozenset]] = [(start, 1.0, frozenset({start}))]
    steps = [(path[i], path[i + 1]) for i in range(1, len(path) - 1, 2)]

    for step_i, (predicate, want_class) in enumerate(steps):
        nxt: list[tuple[str, float, frozenset]] = []
        last = step_i == len(steps) - 1
        for node, weight, visited in frontier:
            candidates = g.out.get((predicate, node), [])
            reverse = not candidates and bool(g.inv.get((predicate, node)))
            if reverse:
                candidates = g.inv.get((predicate, node), [])
            for nbr, edge in candidates:
                if nbr in visited or g.klass(nbr) != want_class:
                    continue
                if not _qualifier_ok(edge, predicate, constraints):
                    continue
                # Hub and oversized-pathway exclusion applies to INTERMEDIATE
                # nodes only. Endpoints are the answer, not a shortcut.
                if not last and nbr in blocked:
                    continue
                d = g.degree.get(nbr, 1) or 1
                nxt.append((nbr, weight * d ** -w, visited | {nbr}))
            if len(nxt) > max_paths:
                break
        frontier = nxt
        if not frontier:
            return {}

    out: dict[str, float] = defaultdict(float)
    for node, weight, _v in frontier:
        out[node] += weight
    return dict(out)


def dwpc_scores(g: Graph, metapaths: dict, limits: dict,
                w: float = DEFAULT_W) -> dict[str, dict[str, dict[str, float]]]:
    """{metapath: {drug: {virus: score}}} for every drug with any path."""
    blocked = _hub_set(g, limits) | _oversized_pathways(g, limits)
    results: dict[str, dict[str, dict[str, float]]] = {}
    for name, mp in metapaths.items():
        per_drug: dict[str, dict[str, float]] = {}
        for drug in g.drugs():
            hits = walk_metapath(g, drug, mp["path"], mp.get("constraints"),
                                 blocked if mp.get("exclude_hubs") else set(), w)
            if hits:
                per_drug[drug] = hits
        results[name] = per_drug
    return results


def combine(dwpc: dict[str, dict[str, dict[str, float]]], virus: str,
            normalize: bool = True) -> dict[str, float]:
    """Combine DWPC across metapaths for one virus.

    RAW SUMMATION IS WRONG and was the first version of this function. The
    metapaths produce incomparable scales: M1 yields one specific path
    (drug -> nsp5 -> gene -> virus, score ~0.3) while M4 yields hundreds of
    diffuse ones summing above 2.0. Summing them lets the highest-variance
    route win by construction.

    Measured on SARS-CoV-2, every metapath carries real signal -- lift over
    its own pool's base rate: M4 2.84, M1 2.64, M3 2.18, M5 2.03, M2 1.82 --
    but M1's pool has a 15.5% positive rate against M4's 6.3%, so M1 is the
    more informative channel while producing the smaller numbers. Under raw
    summation the top 25 became a pure kinase-inhibitor list, nirmatrelvir
    fell to rank 2,973, and hydroxychloroquine outranked it.

    Percentile rank within each metapath makes the scales comparable without
    introducing a fitted weight, so this stays a no-training baseline. A drug
    absent from a metapath's pool contributes nothing rather than a zero:
    having no route of one kind is not evidence against the routes it does
    have.
    """
    if not normalize:
        total: dict[str, float] = defaultdict(float)
        for per_drug in dwpc.values():
            for drug, by_virus in per_drug.items():
                if virus in by_virus:
                    total[drug] += by_virus[virus]
        return dict(total)

    # MAX, not sum. Summing percentiles rewards APPEARING IN MANY POOLS rather
    # than scoring well in any: a promiscuous kinase inhibitor present in all
    # five metapaths reaches ~4.0 by stacking five mediocre percentiles, while
    # nirmatrelvir -- which has one clean M1 route and no polypharmacology --
    # is capped at 1.0. Measured: under sum-of-percentiles chloroquine ranked
    # 142 and nirmatrelvir 1,428, and ranks 6-25 were tied within 0.08.
    #
    # Max means a drug is ranked by its BEST mechanistic route. M1's pool has a
    # 15.5% positive rate and MRR 1.0, so one strong direct-acting path is
    # better evidence than five diffuse host-directed ones.
    best: dict[str, float] = defaultdict(float)
    for per_drug in dwpc.values():
        pool = {drug: by[virus] for drug, by in per_drug.items() if virus in by}
        if not pool:
            continue
        # A single-occupant pool has no defined percentile spread, but the drug
        # IS the top of that pool -- dropping it would silently discard every
        # drug whose only route is a rare metapath.
        if len(pool) == 1:
            drug = next(iter(pool))
            best[drug] = max(best[drug], 1.0)
            continue
        ordered = sorted(pool.items(), key=lambda kv: kv[1])
        n = len(ordered) - 1
        for i, (drug, _v) in enumerate(ordered):
            pct = i / n
            best[drug] = max(best[drug], pct)
    return dict(best)


def degree_ranking(g: Graph) -> dict[str, float]:
    """The null hypothesis: rank drugs by degree, ignoring the virus entirely."""
    return {d: float(g.degree.get(d, 0)) for d in g.drugs()}


def _bfs_distances(g: Graph, sources: set[str], predicate: str,
                   blocked: set[str], max_hops: int = 4) -> dict[str, int]:
    dist = {s: 0 for s in sources}
    frontier = set(sources)
    for hop in range(1, max_hops + 1):
        nxt: set[str] = set()
        for node in frontier:
            for nbr, _e in (g.out.get((predicate, node), [])
                            + g.inv.get((predicate, node), [])):
                if nbr in dist or nbr in blocked:
                    continue
                dist[nbr] = hop
                nxt.add(nbr)
        frontier = nxt
        if not frontier:
            break
    return dist


def network_proximity(g: Graph, virus: str, limits: dict,
                      n_random: int = 100, seed: int = 0) -> dict[str, float]:
    """Z-scored closest distance from a drug's targets to the virus's host
    factors, against degree-matched random target sets.

    The z-score is what makes this a real baseline rather than a restatement of
    degree: a drug whose targets are close to the host factors only because it
    has many well-connected targets scores near zero.
    """
    rng = random.Random(seed)
    blocked = _hub_set(g, limits)

    factors = {n for n, _e in g.inv.get(("HOST_FACTOR_FOR", virus), [])
               if (_e.get("qualifiers") or {}).get("direction") == "dependency"}
    if not factors:
        return {}
    dist = _bfs_distances(g, factors, "PHYSICALLY_INTERACTS_WITH", blocked)

    proteins = [n for n, d in g.nodes.items()
                if d["class"] == "Protein" and not d["properties"].get("is_viral")]
    by_bin: dict[int, list[str]] = defaultdict(list)
    for p in proteins:
        by_bin[int(math.log2(g.degree.get(p, 1) or 1))].append(p)

    def closest(targets: list[str]) -> float | None:
        ds = [dist[t] for t in targets if t in dist]
        return sum(ds) / len(ds) if ds else None

    drug_targets: dict[str, list[str]] = defaultdict(list)
    for (pred, drug), edges in g.out.items():
        if pred == "TARGETS":
            drug_targets[drug].extend(t for t, _e in edges)

    out: dict[str, float] = {}
    for drug, targets in drug_targets.items():
        obs = closest(targets)
        if obs is None:
            continue
        null: list[float] = []
        for _ in range(n_random):
            sample = []
            for t in targets:
                bucket = by_bin[int(math.log2(g.degree.get(t, 1) or 1))] or proteins
                sample.append(rng.choice(bucket))
            v = closest(sample)
            if v is not None:
                null.append(v)
        if len(null) < 10:
            continue
        mu = sum(null) / len(null)
        sd = (sum((x - mu) ** 2 for x in null) / len(null)) ** 0.5
        # Negated so that higher is better, consistent with the other scorers.
        out[drug] = -((obs - mu) / sd) if sd > 0 else 0.0
    return out


def rank(scores: dict[str, float]) -> list[tuple[str, float]]:
    return sorted(scores.items(), key=lambda kv: -kv[1])


def hits_at_k(ranked: list[tuple[str, float]], positives: set[str], k: int) -> int:
    return sum(1 for d, _s in ranked[:k] if d in positives)


def mrr(ranked: list[tuple[str, float]], positives: set[str]) -> float:
    for i, (d, _s) in enumerate(ranked, 1):
        if d in positives:
            return 1.0 / i
    return 0.0
