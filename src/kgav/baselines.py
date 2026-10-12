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

from kgav.schema import parse_hop

DEFAULT_W = 0.4


def label_publication_index(release: Path, label_predicates: set[str]
                            ) -> dict[str, set[str]]:
    """subject -> the publications its LABEL edges cite.

    Feeds Graph.load(drop_pubs_shared_with=...). Measured on v0.1: of 1,070
    compounds carrying both an INHIBITS edge and a
    HAS_ANTIVIRAL_ACTIVITY_AGAINST label where both sides cite a publication,
    1,015 -- 94.9% -- share at least one PMID.

    That is not bad data. A paper reporting an Mpro IC50 almost always reports
    a cell-based EC50 in the same table, and ChEMBL files them as two
    activities. But it means M1 traverses a measurement taken alongside the
    label it is scored against, so its AUC is partly retrieval. Neither the
    compound-level split nor the temporal split removes it: both facts hang
    off one compound and carry one date.
    """
    idx: dict[str, set[str]] = defaultdict(set)
    for line in (Path(release) / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        if e["predicate"] in label_predicates:
            pubs = {p for p in (e.get("publications") or []) if p}
            if pubs:
                idx[e["subject"]] |= pubs
    return dict(idx)


class Graph:
    """Adjacency built for metapath walking.

    out[(predicate, node)] -> [(neighbour, edge)] in the schema's declared
    direction; inv[...] is the reverse, needed because M4 traverses
    PARTICIPATES_IN forwards then backwards.

    Which of the two a hop reads is decided by neighbours(), from the schema's
    `symmetric` declaration and the hop's own direction marker -- never from
    what a particular node happens to have.
    """

    def __init__(self, symmetry: dict[str, bool] | None = None) -> None:
        self.nodes: dict[str, dict] = {}
        self.out: dict[tuple[str, str], list[tuple[str, dict]]] = defaultdict(list)
        self.inv: dict[tuple[str, str], list[tuple[str, dict]]] = defaultdict(list)
        self.degree: Counter = Counter()
        # predicate -> symmetric, straight from the schema. A predicate absent
        # here has no declared symmetry and MUST NOT be traversed: see
        # neighbours().
        self.symmetry: dict[str, bool] = dict(symmetry or {})
        # load-time accounting: edges seen, skipped as held-out, dropped as
        # publication-coincident with a label. Empty unless load() filled it.
        self.stats: Counter = Counter()

    def neighbours(self, predicate: str, node: str, reverse: bool = False
                   ) -> list[tuple[str, dict]]:
        """Adjacent nodes for one hop, in the direction the SCHEMA declares.

        A symmetric predicate reads both adjacency maps, because which side a
        pair was stored on is an ingest artifact -- STRING canonicalises each
        interaction to one row ordered by accession. A directed predicate
        reads exactly the map the hop asks for.

        This replaces `follow outgoing edges, and fall back to incoming only
        when there are none`, which truncated every node holding edges in both
        directions -- 9,016 host proteins in v0.1.

        The union does not dedup, because on v0.1 it cannot double-count:
        zero symmetric pairs are stored in BOTH orientations (checked for
        PHYSICALLY_INTERACTS_WITH and CHEMICALLY_SIMILAR_TO, 0 of 150,432 and
        0 of 1,583). That is an invariant of the ingests -- STRING is
        canonicalised on read and the similarity layer is bipartite -- not of
        this function. An ingest that emitted both orientations of one fact
        would inflate its DWPC weight twofold.
        """
        symmetric = self.symmetry.get(predicate)
        if symmetric is None:
            raise ValueError(
                f"predicate {predicate!r} has no declared symmetry, so its "
                f"traversal direction is unknown. Declare `symmetric` on it in "
                f"kg_schema.yaml and load the graph with Schema.symmetry().")
        if symmetric:
            return self.out.get((predicate, node), []) + self.inv.get((predicate, node), [])
        return (self.inv if reverse else self.out).get((predicate, node), [])

    @classmethod
    def load(cls, release: Path, skip_predicates: set[str] | None = None,
             symmetry: dict[str, bool] | None = None,
             drop_pubs_shared_with: dict[str, set[str]] | None = None) -> Graph:
        """Build the traversal graph.

        drop_pubs_shared_with maps a node to the publications its LABEL edges
        cite (see label_publication_index). An edge whose SUBJECT appears there
        and which cites any of the same publications is withheld from
        traversal: the measurement it records came out of the same paper as the
        answer, so a path across it is retrieval rather than prediction.

        Keyed on the subject deliberately. Compound A's label PMID appearing on
        compound B's edge is not leakage for B -- B is scored against B's own
        label -- so a global PMID blocklist would discard sound evidence. The
        edge is also only dropped from TRAVERSAL; degree still counts it, for
        the same reason held-out predicates do.
        """
        skip = skip_predicates or set()
        shared = drop_pubs_shared_with or {}
        g = cls(symmetry)
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
            g.stats["edges_seen"] += 1
            if p in skip:
                g.stats["skipped_held_out"] += 1
                continue
            label_pubs = shared.get(e["subject"])
            if label_pubs and label_pubs & {q for q in (e.get("publications") or []) if q}:
                g.stats["dropped_publication_coincident"] += 1
                g.stats[f"dropped_{p}"] += 1
                continue
            g.out[(p, e["subject"])].append((e["object"], e))
            g.inv[(p, e["object"])].append((e["subject"], e))
            g.stats["traversable"] += 1
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


def qualifier_ok(edge: dict, predicate: str, constraints: dict | None) -> bool:
    """Whether one edge satisfies a metapath's qualifier constraints.

    PUBLIC because the agent needs it too. src/kgav/agent/tools.py enumerates
    M1 and M2 as its own loops and applied no constraints to M1 and a
    HAND-COPIED one to M2, so a constraint added to the schema reached the
    scorer and not the agent -- which is how the agent came to explain routes
    through edges the ranking refuses to walk. One implementation, read from
    the schema by both.

    A MISSING QUALIFIER IS A FAILED MATCH. That is why the schema requires
    target_resolution on INHIBITS: optional, an edge emitted without it would
    vanish from the direct-acting channel in silence.
    """
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
    steps = [(parse_hop(path[i]), path[i + 1]) for i in range(1, len(path) - 1, 2)]

    for step_i, ((predicate, reverse), want_class) in enumerate(steps):
        nxt: list[tuple[str, float, frozenset]] = []
        last = step_i == len(steps) - 1
        for node, weight, visited in frontier:
            candidates = g.neighbours(predicate, node, reverse)
            for nbr, edge in candidates:
                if nbr in visited or g.klass(nbr) != want_class:
                    continue
                if not qualifier_ok(edge, predicate, constraints):
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


def metapath_reach(g: Graph, metapaths: dict,
                   dwpc: dict[str, dict[str, dict[str, float]]]) -> dict[str, dict]:
    """Per metapath: whether it reached anything, and why it could not.

    EVERY DECLARED METAPATH IS REPORTED, INCLUDING THE ONES THAT REACH NOTHING.
    That is the whole point. M8 was declared in schema 0.8.0, described there as
    "THE ROUTE v2 EXISTS FOR", and has produced zero paths on every run since:
    its TargetClass nodes and MEMBER_OF_CLASS edges are not in the assembled
    release at all, and the FOLD_SIMILAR_TO edges the fold ingest writes are
    self-loops that walk_metapath discards as already-visited.

    None of that was recorded anywhere. compare_metapaths printed "no paths" to
    stdout and then omitted the metapath from its results file; hard_negatives
    skipped it with `if not sc: continue` and printed nothing. A metapath that
    silently contributes nothing is indistinguishable from one that was never
    declared, so the absence has to be a reported value rather than a missing
    key.

    `missing_predicates` is the diagnosis: a hop whose predicate has no edges in
    this graph explains the zero without anyone having to reconstruct it. The
    project's own stated integrity check is comparing counts across entities
    that should be similar, and this is that check applied to metapaths.
    """
    present = {k[0] for k in g.out} | {k[0] for k in g.inv}
    out: dict[str, dict] = {}
    for name, mp in sorted(metapaths.items()):
        path = list(mp.get("path") or [])
        hops = [parse_hop(path[i])[0] for i in range(1, len(path) - 1, 2)]
        drugs = dwpc.get(name) or {}
        viruses = {v for by in drugs.values() for v in by}
        out[name] = {
            "drugs_reached": len(drugs),
            "viruses_reached": sorted(viruses),
            "predicates": hops,
            # Deduped, order preserved. M8 traverses MEMBER_OF_CLASS twice --
            # forward then reverse -- so the raw list printed it twice and read
            # as two separate problems.
            "missing_predicates": list(dict.fromkeys(
                p for p in hops if p not in present)),
        }
    return out


# Potency fields in FIXED PRECEDENCE, never chosen by value. selectivity.py
# records what choosing by value costs: "Taking the worst ratio defamed
# nirmatrelvir (SI 1.8) and remdesivir (0.9). Taking the best exonerated
# chloroquine (17.7) and hydroxychloroquine (21.2). Both were the same error
# -- selecting a study by its answer." So the FIELD is picked by this order
# and only then is a value read. ec50_nm is absent deliberately: it is a
# cellular readout, and chembl.py files cellular rows on the organism rather
# than the protein.
POTENCY_PRECEDENCE = ("ic50_nm", "ki_nm", "kd_nm")

# A relation in this set is a BOUND, not a measurement. chembl.py's rule:
# "IC50 > 10000 nM is a NEGATIVE result -- the compound was measured and
# found inactive -- and coercing it to a potency number inverts its meaning."
CENSORED_RELATIONS = frozenset({">", "<", ">=", "<=", "~"})


def potency_key(edge: dict) -> tuple[float, str] | None:
    """(pX, field) for one INHIBITS edge, or None if it carries no measurement.

    pX is the standard 9 - log10(nM), a unit conversion rather than a fitted
    constant: 75 nM becomes 7.125, 7,280 nM becomes 5.138. Higher is more
    potent.

    None is returned for a censored bound and for an edge with no potency
    field at all, and the two are NOT distinguished here because neither can
    order anything. A caller that needs the difference counts the relations
    itself.
    """
    q = edge.get("qualifiers") or {}
    if q.get("relation") in CENSORED_RELATIONS:
        return None
    for field in POTENCY_PRECEDENCE:
        v = q.get(field)
        if v is None:
            continue
        try:
            v = float(v)
        except (TypeError, ValueError):
            continue
        if v <= 0:                      # the schema permits 0; log10 does not
            continue
        return 9.0 - math.log10(v), field
    return None


def potency_tiebreak(g: Graph, virus: str, path: list,
                     constraints: dict | None = None
                     ) -> tuple[dict[str, float], Counter]:
    """{drug: pX} over the INHIBITS edges this metapath actually walks.

    Only edges to a protein with an onward route to `virus` count, and the
    metapath's own qualifier constraints are honoured, so this keys on the
    evidence the channel uses rather than on every INHIBITS edge in the
    graph.

    WHERE SEVERAL EXACT VALUES EXIST THE WORST IS KEPT -- the least potent,
    the lowest pX. That is selectivity.py's rule for the same situation ("among
    several exact measurements the WORST is kept, since a compound shown
    cytotoxic in a real measurement is cytotoxic"), and it is the direction
    that cannot flatter a candidate.

    Returns the map and a Counter of which field each drug was keyed on, so
    prediction P5 -- that a drug ordered on its Ki is not strictly comparable
    to one ordered on its IC50 -- can be sized rather than assumed small.
    """
    want = (constraints or {}).get("INHIBITS.target_resolution")
    reaching = set()
    for prot in g.nodes:
        if g.klass(prot) != "Protein" or not g.is_viral(prot):
            continue
        for gene, _e in g.neighbours("ENCODED_BY", prot):
            for tax, _e2 in g.neighbours("BELONGS_TO", gene):
                if tax == virus:
                    reaching.add(prot)
    keys: dict[str, float] = {}
    fields: Counter = Counter()
    for drug in g.drugs():
        best: tuple[float, str] | None = None
        for nbr, edge in g.neighbours("INHIBITS", drug):
            if nbr not in reaching:
                continue
            if want is not None and (edge.get("qualifiers") or {}).get(
                    "target_resolution") != want:
                continue
            got = potency_key(edge)
            if got is None:
                continue
            if best is None or got[0] < best[0]:   # worst = least potent
                best = got
        if best is not None:
            keys[drug] = best[0]
            fields[best[1]] += 1
    return keys, fields


def channel_percentiles(pool: dict[str, float],
                        tiebreak: dict[str, float] | None = None
                        ) -> dict[str, float]:
    """Percentile rank within one channel's pool, ties taking their midrank.

    THE ONE IMPLEMENTATION. This arithmetic used to live inside combine, and
    scripts/diagnose_controls.py carried a second copy so it could show the
    percentile beside the combined score. When combine moved to midranks the
    copy did not, and the diagnostic printed a percentile of 0.0182 next to a
    combined score of 0.4041 for the same drug on a single-channel route --
    two numbers that cannot both be right, from a script written to make the
    scoring checkable. Exported here so there is nothing to keep in sync.

    TIES TAKE THEIR MIDRANK. `pct = i / n` read the POSITIONAL index, so
    drugs holding the identical raw DWPC received different percentiles
    according to where sorted() happened to place them -- which, for equal
    keys, is dict insertion order. Measured on v0.1 M1 for SARS-CoV-2:
    nirmatrelvir and chloroquine both score 0.000242, both by the one path
    through nsp5, and were handed 0.0180 and 0.0020. nirmatrelvir's reported
    rank of 6,432 was therefore not a measurement of anything; reinserting
    the same drugs in another order would have moved it.

    The blocks are wide because DWPC is TOPOLOGY ONLY: every drug whose sole
    route to this virus is nsp5 traverses the same nodes and so earns the
    same number, potency notwithstanding. That is a separate question from
    this one; what is fixed here is that the arithmetic must not invent an
    order it did not measure.

    A single-occupant pool has no defined spread, but the drug IS the top of
    it -- returning 0.0 would silently discard every drug whose only route is
    a rare metapath.

    Equality is exact, not within a tolerance. These ties arise from
    identical topology -- the same product of the same degrees -- so they are
    bit-identical by construction, and a tolerance would instead merge scores
    that genuinely differ while adding a fitted parameter to a no-training
    baseline. Scores separated only by floating-point association are left
    separate, which is the conservative direction: this never merges what is
    not identical.
    """
    if not pool:
        return {}
    ordered = sorted(pool.items(), key=lambda kv: kv[1])
    n = len(ordered) - 1
    if n <= 0:
        return {ordered[0][0]: 1.0}
    out: dict[str, float] = {}
    i = 0
    while i < len(ordered):
        j = i
        while j + 1 < len(ordered) and ordered[j + 1][1] == ordered[i][1]:
            j += 1
        mid = ((i + j) / 2) / n              # midrank of positions i..j
        block = [d for d, _v in ordered[i:j + 1]]
        keyed = [d for d in block if tiebreak and d in tiebreak]
        if len(keyed) < 2:
            # Nothing to order by -- one value cannot rank a block, and no
            # value certainly cannot. Everyone keeps the midrank.
            for drug in block:
                out[drug] = mid
        else:
            # A DRUG WITH NO KEY KEEPS THE SCORE IT ALREADY HAS: the block's
            # midrank. So this is a strict refinement -- no imputation
            # decides anything, and "unmeasured" is never read as "weak",
            # which is the error chembl.py refuses when it declines to
            # coerce a censored bound into a number.
            keyed_set = set(keyed)
            for drug in block:
                if drug not in keyed_set:
                    out[drug] = mid
            # The keyed drugs spread across THE BLOCK'S OWN SPAN, i..j, so
            # the block occupies exactly the percentile range it did before
            # and nothing outside it moves. Ties in the key itself take their
            # own midrank, by the same rule one level up.
            keyed.sort(key=lambda d: tiebreak[d])
            span, last = j - i, len(keyed) - 1
            r = 0
            while r <= last:
                t = r
                while t + 1 <= last and tiebreak[keyed[t + 1]] == tiebreak[keyed[r]]:
                    t += 1
                pos = i + span * ((r + t) / 2) / last
                for drug in keyed[r:t + 1]:
                    out[drug] = pos / n
                r = t + 1
        i = j + 1
    return out


def combine(dwpc: dict[str, dict[str, dict[str, float]]], virus: str,
            normalize: bool = True, calibration=None,
            tiebreak: dict[str, dict[str, float]] | None = None
            ) -> dict[str, float]:
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
    # MAX IS BIASED UPWARD BY THE NUMBER OF CHANNELS, which is harmless while
    # every channel carries signal and corrosive once some do not. The maximum
    # of k independent uniform percentiles has expectation k/(k+1): over eight
    # channels that is 0.89 for any drug present in all of them, so the
    # ordering reverts to counting pools -- the very failure the switch from
    # sum to max was made to fix -- and a drug at the top of a pure-noise pool
    # scores 1.0 outright. Measured on the single-screen protocol, COMBINED
    # lands at 0.516 while the best single channel reaches 0.523.
    #
    # `calibration` (kgav.calibration.Calibration) restricts the max to
    # channels whose measured 95% interval EXCLUDES 0.5, and weights each by
    # its margin over chance. With none qualifying it returns {} rather than a
    # ranking: refusing to order candidates is a result, and an order built
    # from eight chance-level channels is not one. Omit it and the uncalibrated
    # behaviour above is unchanged, so every published number reproduces.
    weights: dict[str, float] | None = None
    if calibration is not None:
        weights = calibration.weights()
        if not weights:
            return {}

    best: dict[str, float] = defaultdict(float)
    for name, per_drug in dwpc.items():
        if weights is not None and name not in weights:
            continue
        pool = {drug: by[virus] for drug, by in per_drug.items() if virus in by}
        if not pool:
            continue
        # The single-occupant case lives in channel_percentiles, which returns
        # 1.0 for it: that drug IS the top of its pool, and dropping it would
        # silently discard every drug whose only route is a rare metapath.
        w = weights[name] if weights is not None else 1.0
        # PER CHANNEL, because a tiebreak is a claim about one kind of
        # evidence. Potency against a viral protein orders M1; it says
        # nothing about M2-M6, which are host channels, and nothing about M7,
        # which borrows a NEIGHBOUR's target so the drug's own potency is not
        # what that path asserts.
        keys = (tiebreak or {}).get(name)
        for drug, pct in channel_percentiles(pool, keys).items():
            best[drug] = max(best[drug], pct * w)
    return dict(best)


MIN_TRUSTWORTHY_EXPECTATION = 1.0


def lift_pvalue(pool: int, positives: int, k: int, hits: int) -> float | None:
    """P(hits >= observed) under random ranking: the exact hypergeometric tail.

    THE LIFT EQUIVALENT OF auc_interval. Every AUC in this project carries a
    Hanley-McNeil interval, and lift carried only the `trustworthy` heuristic
    -- expectation >= 1 and pool > k -- which answers "could this be read at
    all", not "could this have happened by chance". On the 2021 temporal split
    the two disagree:

      SARS-CoV-2 M1   pool 296   8 pos   6 hits   exp 2.7   lift 2.22  p=0.020
      MERS-CoV   M4   pool 2,682 5 pos   3 hits   exp 0.2   lift 16.1  p=0.0007

    M4 is flagged untrustworthy because its expectation is 0.19, yet 3 hits
    from 5 positives in a 2,682 pool is the least likely row in the run. The
    heuristic is a guard against unreadable figures, not a test, and it was
    being read as one.

    Drawing the top k from a pool of `pool` containing `positives` is sampling
    WITHOUT replacement, so the null is hypergeometric, not Poisson. Computed
    in log space: pool 4,316 choose 100 is a 240-digit integer and the ratio
    is all that matters.

    Returns None when the question is empty (no positives, or k covers the
    whole pool so every draw is identical).
    """
    if positives <= 0 or pool <= 0 or k <= 0 or k >= pool:
        return None
    upper = min(positives, k)
    if hits > upper:
        return 0.0
    def _lc(n: int, r: int) -> float:
        if r < 0 or r > n:
            return float("-inf")
        return (math.lgamma(n + 1) - math.lgamma(r + 1)
                - math.lgamma(n - r + 1))
    denom = _lc(pool, k)
    terms = [_lc(positives, i) + _lc(pool - positives, k - i) - denom
             for i in range(max(hits, 0), upper + 1)]
    terms = [x for x in terms if x != float("-inf")]
    if not terms:
        return 0.0
    m = max(terms)
    return min(1.0, math.exp(m) * sum(math.exp(x - m) for x in terms))


DEGREE_KEY = "degree"
SIGNIFICANCE = 0.05


def degree_relative(per: dict[str, dict], degree_key: str = DEGREE_KEY) -> None:
    """Annotate every row with its lift relative to the degree baseline.

    WHY A P-VALUE IS NOT ENOUGH. lift_pvalue's null is RANDOM RANKING, so it
    answers "could this ordering have arisen by chance" and nothing else. It
    cannot see a confound that inflates the scorer and the labels together --
    which is the confound this project spent its whole evaluation on: 94.9% of
    compounds share a PMID between their INHIBITS edge and their antiviral
    label, and degree IS study volume.

    So a p-value makes the UNCONTROLLED protocol look full of discoveries. On
    v0.1, compare_metapaths reports p < 0.0001 for SARS-CoV M1 -- and p <
    0.0001 for degree on the same labels, at 9.51 lift against M1's 2.10. In
    every cross-sectional virus where a channel cleared 0.05, degree cleared
    it too; there is not one exception. The comparison that carries
    information is channel against DEGREE, not channel against random.

    Degree's own p-value is the diagnostic. If degree is significant, the
    protocol is measuring how much a compound has been studied, and a channel
    beating random on those labels says nothing about mechanism.
    """
    deg = per.get(degree_key)
    base = (deg or {}).get("lift")
    for name, m in per.items():
        if name == degree_key:
            m["lift_vs_degree"] = None
            m["beats_degree"] = None
            m["vs_degree_readable"] = None
            continue
        if base is None:
            # No baseline row at all, so there is nothing to compare against.
            m["lift_vs_degree"] = None
            m["beats_degree"] = None
            m["vs_degree_readable"] = None
            continue
        # A RATIO IS ONLY EVIDENCE WHEN THE DENOMINATOR IS. Two cases broke
        # the first version of this column, both visible in the v0.1 run:
        #
        #  base == 0  HCoV-229E temporal computed_only: degree lift 0.00 and
        #             M2 lift 4.65. The ratio is undefined, `not base` sent it
        #             down the no-baseline branch, and the one row where a
        #             channel beat degree outright printed as "-" --
        #             indistinguishable from the baseline row itself.
        #             beats_degree carries it instead of a ratio, and NOT
        #             float("inf"), which json.dumps writes as `Infinity` and
        #             would make every results file invalid JSON.
        #
        #  base < 1   SARS-CoV-2 cross-sectional: degree lift 0.45, p=0.99 --
        #             WORSE than random. M4's 4.44x against it reads as
        #             decisive and is 4.44 * 0.45 = 1.99, p=0.14. Beating a
        #             sub-random baseline is not an achievement, so the ratio
        #             is marked unreadable rather than shown bare.
        m["beats_degree"] = m["lift"] > base
        m["lift_vs_degree"] = (m["lift"] / base) if base else None
        # Readable means the DENOMINATOR is a real effect. When degree is
        # itself at chance the ratio divides by noise, so it describes rather
        # than tests -- SARS-CoV-2 cross-sectional reads M4 at 4.44x a degree
        # lift of 0.45 (p=0.99), which is M4 lift 1.99 at p=0.14. My first
        # threshold here was `base >= 1.0`, which wrongly flagged the clean
        # SARS-CoV-2 temporal row where degree sits at 0.93, p=0.63: that is
        # degree AT CHANCE, which is the condition a real result needs, not a
        # defect in it.
        dp = (deg or {}).get("p_value")
        m["vs_degree_readable"] = dp is not None and dp < SIGNIFICANCE


def confound_warning(per: dict[str, dict], degree_key: str = DEGREE_KEY,
                     alpha: float = SIGNIFICANCE) -> str | None:
    """The one-line read of a results table, or None when it needs no warning.

    Returns a warning when the degree baseline is itself significant, because
    that is the signature of a protocol measuring study volume rather than
    biology.
    """
    deg = per.get(degree_key)
    if not deg:
        return None
    dp, dl = deg.get("p_value"), deg.get("lift")
    # BOTH are required to read the baseline. Without dl there is nothing to
    # compare against; without dp there is no way to say whether the baseline
    # is at chance, so neither verdict can be given. This also removes a
    # latent TypeError: the at-chance branch formats dp, and lift_pvalue
    # returns None whenever a row has no positives or k covers the pool.
    if dl is None or dp is None:
        return None

    # THE BASELINE'S OWN HITS MUST BE DETERMINED, because every verdict below
    # is a comparison against it -- and so is the lift_vs_degree ratio in the
    # `vs deg` column, which divides by this row.
    #
    # THIS IS THE THIRD TIME THIS FUNCTION HAS READ ONE FIELD AND IGNORED THE
    # ONE BESIDE IT. It once read p and ignored `trustworthy`, and promoted
    # HCoV-OC43 M4 on a single positive; the comment below records that. It
    # then read p and `trustworthy` and ignored `hits_determined`, which is
    # how it came to print "BASELINE AT CHANCE: degree lift 0.97 (p=0.5909)
    # is indistinguishable from random here, so M1's lift 2.49 (p=0.0033) is
    # not explained by study volume" for SARS-CoV-2 train2021 -- where M1's
    # Hits@100 was 0-8 of 8 positives and degree's was 3-5. A sentence
    # comparing two intervals as though they were points.
    if deg.get("hits_determined") is False:
        return (f"BASELINE NOT DETERMINED: the degree row's Hits@k is one "
                f"draw from [{deg.get('hits_lo')}, {deg.get('hits_hi')}] "
                f"because a tie block straddles the cut "
                f"({deg.get('tied_at_cut')} tied, "
                f"{deg.get('tied_at_cut_positives')} of them positive). Its "
                f"lift {dl:.2f} and p={dp:.4f} are not measurements, so no "
                f"comparison against this baseline can be read -- including "
                f"every ratio in the vs-deg column, which divides by it.")

    others = [(m["lift"], m.get("p_value"), n, m.get("trustworthy"),
               m.get("hits_determined", True))
              for n, m in per.items()
              if n != degree_key and m.get("lift") is not None]
    better = sorted(x for x in others if x[0] > dl)
    if better:
        lift, _p, name, tw, det = better[-1]
        if det is False:
            tw = False          # an undetermined row is not a clean exceedance
        mark = "" if tw else " -- but that row is flagged (!)"
        tail = f"{name} exceeds it at {lift:.2f}{mark}"
    else:
        tail = "no channel exceeds it"

    if dp is not None and dp < alpha:
        return (f"CONFOUNDED: the degree baseline is itself significant "
                f"(lift {dl:.2f}, p={dp:.4f}). These labels track how much a "
                f"compound has been studied, so a channel beating RANDOM says "
                f"nothing about mechanism -- only beating DEGREE does. "
                f"{tail}.")

    # THE GOOD CASE, AND IT DESERVES SAYING OUT LOUD. Degree at chance is the
    # condition every real result in this project depends on, and it held in
    # exactly one place: SARS-CoV-2 under the temporal protocol. Printing it
    # only as an absence of warning left the reader to notice it.
    # TRUSTWORTHY IS REQUIRED, NOT JUST SIGNIFICANT. This branch read p and
    # ignored the field beside it, and promoted HCoV-OC43 temporal M4: pool
    # 2,545, ONE positive, lift 25.45, p=0.0393 -- and 100/2545 = 0.0393
    # exactly, so that p IS the sampling fraction for a single compound. The
    # row already carried (!) for expectation < 1. Requiring True explicitly
    # rather than `is not False`, because a row with no verdict has not earned
    # a promotion either.
    # hits_determined IS REQUIRED HERE TOO, for the same reason `trustworthy`
    # is: a channel whose Hits@k is an interval has no p-value to promote.
    sig = sorted((p_, lift, n) for lift, p_, n, tw, det in others
                 if p_ is not None and p_ < alpha and tw is True
                 and det is True)
    if sig:
        p_, lift, n = sig[0]
        return (f"BASELINE AT CHANCE: degree lift {dl:.2f} (p={dp:.4f}) is "
                f"indistinguishable from random here, so {n}'s lift {lift:.2f} "
                f"(p={p_:.4f}) is not explained by study volume. A channel's p "
                f"only means what it appears to mean under this condition.")
    return None


def hits_at_k_range(ranked: list[tuple[str, float]], positives: set[str],
                    k: int) -> tuple[int, int, int, int]:
    """The range Hits@k can take over tie orderings, and the tie at the cut.

    rank() is `sorted(scores.items(), key=lambda kv: -kv[1])` with NO
    tiebreak, and Python's sort is stable, so equal scores come out in dict
    insertion order. hits_at_k then reads ranked[:k]. Where a tie block
    STRADDLES position k, which of its members fall inside the cut is decided
    by insertion order rather than by any score -- so Hits@k, MRR and lift
    are reproducible for one dict and are not properties of the graph.

    Measured on v0.1 for SARS-CoV-2, M1, with 2,020 measured actives:

        k     Hits@k range   tied at cut   of which positive
        10         1 -  2              4                   1
        50         2 - 12             98                  10
       100         2 - 12             98                  10

    A 98-wide block of identical raw DWPC spans positions ~30 to ~128 and
    holds 10 positives, so both cuts land inside it and the published figure
    is one draw from [2, 12]. COMBINED is determined at all three k on the
    same release, because the blocks at its cuts hold no positives -- so this
    is not a reason to distrust every number, which is exactly why it is
    measured per row rather than asserted.

    Returns (low, high, contested, contested_positives). low is the count if
    every tied non-positive is placed first, high if every tied positive is.
    low == high means the cut lands in no tie and the figure is a
    measurement. A block of ONE is not a tie and reports 0 contested: that
    drug is simply the kth.
    """
    if len(ranked) <= k:
        n = sum(1 for d, _v in ranked if d in positives)
        return n, n, 0, 0
    cut = ranked[k - 1][1]
    inside = [d for d, v in ranked if v > cut]
    contested = [d for d, v in ranked if v == cut]
    c_pos = sum(1 for d in contested if d in positives)
    fixed = sum(1 for d in inside if d in positives)
    slots = k - len(inside)
    low = fixed + max(0, slots - (len(contested) - c_pos))
    high = fixed + min(slots, c_pos)
    if len(contested) == 1:
        return low, high, 0, 0
    return low, high, len(contested), c_pos


def lift_row(scores: dict[str, float], pos: set[str], k: int = 100) -> dict:
    """Hits@k against the pool's own random expectation.

    WHEN THE POOL IS NO LARGER THAN k, LIFT CARRIES NO INFORMATION. ranked[:k]
    is then the whole pool, so hits == every positive in it, and

        lift = p / (k * p / pool) = pool / k

    exactly -- a function of the pool size and nothing else. Measured on v0.1:
    MERS M1 pool 24 lift 0.24, HCoV-229E M1 pool 17 lift 0.17, MERS M7 pool 14
    lift 0.14, HCoV-229E M7 pool 7 lift 0.07, SARS-CoV M7 pool 100 lift 1.00.
    Every one is pool/k to two decimals, and none was flagged, because the
    guard tested only expectation >= 1 and these expectations run to 57.

    Those values then fed the M6-vs-M7 verdict, so "MERS-CoV: M6 0.88 vs M7
    0.14, M6 better" compared a real lift against a pool-size artifact.
    """
    r = rank(scores)
    pool = len(r)
    p = len(pos & set(scores))
    exp = k * p / pool if pool else 0.0
    h = hits_at_k(r, pos, k)
    # THE RANGE, REPORTED THE SAME WAY `lift_is_pool_artifact` IS: alongside
    # the figure rather than folded into it, because the figure is still the
    # one that was computed and the range says whether it was determined.
    # Where a tie block straddles k the published hit count is one draw from
    # [lo, hi] -- on v0.1 M1 that is [2, 12] at both k=50 and k=100 -- and a
    # lift computed from it inherits the same spread.
    lo, hi, n_tied, n_tied_pos = hits_at_k_range(r, pos, k)
    return {"pool": pool, "pos": p, "hits": h, "expected": exp,
            "lift": (h / exp) if exp else 0.0, "mrr": mrr(r, pos),
            "hits_lo": lo, "hits_hi": hi,
            "hits_determined": lo == hi,
            "tied_at_cut": n_tied, "tied_at_cut_positives": n_tied_pos,
            "lift_lo": (lo / exp) if exp else 0.0,
            "lift_hi": (hi / exp) if exp else 0.0,
            # BOTH conditions. Expectation >= 1 says the hit count is not a
            # coin flip; pool > k says the ranking was actually exercised.
            "trustworthy": exp >= MIN_TRUSTWORTHY_EXPECTATION and pool > k,
            "lift_is_pool_artifact": pool <= k,
            # Reported alongside `trustworthy`, not folded into it: the
            # heuristic says whether the figure is readable, the p-value says
            # whether it is surprising, and collapsing them hid MERS M4.
            "p_value": lift_pvalue(pool, p, k, h)}


def hits_cell(m: dict, width: int = 7) -> str:
    """The hits column for an evaluated row: the figure, or the range.

    Lives here rather than in each script because this project has twice now
    had one calculation in two files and watched them drift -- the percentile
    arithmetic in combine and diagnose_controls, and `trustworthy` in
    baselines and temporal_split (whose own header still records the two
    definitions that once disagreed).

    A determined row prints its hit count. An undetermined one prints
    "lo-hi", because there is no single number to print: the figure that was
    computed is one draw from that interval, decided by dict insertion order
    rather than by any score.
    """
    if m.get("hits_determined", True):
        return f"{m['hits']:>{width}}"
    return f"{m['hits_lo']}-{m['hits_hi']}".rjust(width)


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
            # Reads the same declared symmetry walk_metapath does. These two
            # used to disagree -- this one unioned both directions while
            # walk_metapath guessed -- with nothing authoritative to check
            # either against.
            for nbr, _e in g.neighbours(predicate, node):
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
