"""Day 8: assemble the per-layer releases into one graph.

Per-layer validation cannot see cross-layer problems. This pass exists to catch
the ones that only appear once the layers sit together:

  NODE CONFLICTS. The same id defined in two layers with DIFFERENT property
  values. Merging silently would pick a winner by dict order; a SmallMolecule
  appearing twice with different SMILES is a normalisation bug, not something
  to resolve quietly. Properties merge by a fixed precedence, and every
  genuine disagreement is reported.

  DANGLING REFERENCES BECOME FATAL. The vhppi and orcs layers write edges-only
  releases, so per-layer dangling was expected. After assembly it is an error:
  an edge pointing at a node that exists nowhere is a join that failed.

  DUPLICATE FACTS FROM DIFFERENT SOURCES. The same (subject, predicate, object)
  asserted by two sources is one fact with two witnesses, not two edges.
  Stacking them would double its weight in every path count. Merged edges carry
  both sources, the union of publications, and -- critically -- the EARLIEST
  first_asserted_date. A temporal split must use when a fact was first
  established, not when the second source repeated it.

  THE REAL DEGREE DISTRIBUTION. hub_percentile was calibrated against the host
  layer alone. With compounds and activity edges added the distribution shifts,
  so the assembled graph is where the hub cut actually gets set.
"""
from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any

# (predicate, qualifiers) -> qualifiers with every schema-derived qualifier
# recomputed. Schema.derive_qualifiers satisfies it; passed as a callable so
# this module keeps knowing nothing about the schema file.
DeriveFn = Callable[[str, Mapping[str, Any]], dict]

# Later layers do not overwrite earlier ones on conflict. The spine is
# hand-curated and authoritative for viral entities; the host layer is derived
# from reference proteomes; chembl mints compounds and knows least about
# anything it did not create.
# Earlier wins a qualifier conflict (see the H2 note in the README). "ncats"
# sits last because its position is immaterial: HAS_ANTIVIRAL_ACTIVITY_AGAINST
# declares identity_includes_source, so a panel measurement never shares a key
# with a chembl one and never competes for its qualifiers.
LAYER_PRECEDENCE = ["spine", "host", "vhppi", "orcs", "chembl", "selectivity",
                    "hosttargets", "similarity", "ncats"]


def _read_jsonl(path: Path) -> list[dict]:
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def edge_key(e: dict, identity: dict[str, list[str]] | None = None,
             source_identity: set[str] | None = None) -> tuple:
    """Identity of a fact.

    (subject, predicate, object) alone is WRONG for predicates with a
    direction-bearing qualifier: a dependency edge and a restriction edge for
    the same (protein, virus) are different facts, not two witnesses to one.
    320 pairs in the ORCS layer carry both, and merging on the triple alone let
    whichever source was read first decide the direction.
    """
    base = (e["subject"], e["predicate"], e["object"])
    if source_identity and e["predicate"] in source_identity:
        # A SECOND EXPERIMENT, NOT A SECOND WITNESS. Without this, a panel
        # inactive and a chembl active for one compound collapse into a single
        # edge and precedence picks the qualifiers -- so one real measurement
        # is gone before labels.build_labels can see the two disagree and drop
        # the compound as ambiguous. The safeguard exists; this is what makes
        # it reachable.
        base = base + (e.get("primary_knowledge_source"),)
    if not identity:
        return base
    quals = identity.get(e["predicate"]) or []
    if not quals:
        return base
    eq = e.get("qualifiers") or {}
    return base + tuple(eq.get(q) for q in quals)


# Width of the stable edge id. 64 bits over 371,117 edges puts the birthday
# collision probability near 4e-9, and assembly checks for one anyway rather
# than trusting the arithmetic.
FACT_ID_HEX = 16


def fact_id(key: tuple) -> str:
    """A stable, content-derived id for one fact.

    DERIVED FROM edge_key, WHICH IS THE POINT. The id therefore inherits
    exactly the identity semantics the assembler already uses: two witnesses
    to one fact get one id, a dependency and a restriction edge for the same
    pair get two, and a panel measurement gets a different id from a ChEMBL
    one because HAS_ANTIVIRAL_ACTIVITY_AGAINST is identified partly by source.

    It replaces a LINE INDEX. agent/tools.py reads `edge_id` and falls back to
    f"E{i}" -- the position of the edge in edges.jsonl -- and nothing has ever
    written the field, so every citation the agent has ever emitted was a line
    number. Those renumber whenever the release is rebuilt, which means a
    citation in a saved answer silently comes to point at a different fact. An
    agent that cites has to cite something that holds still.

    Shaped "E" + hex so agent/verify.py's EDGE_ID_RE keeps matching, and so
    old-style E<digits> ids still parse.
    """
    payload = "\x1f".join("" if k is None else str(k) for k in key)
    return "E" + hashlib.blake2b(payload.encode("utf-8"),
                                 digest_size=FACT_ID_HEX // 2).hexdigest()


# The project's marker for "this fact has no assertion date", not a date in
# 1970. temporal.edge_year and labels.build_labels both already read it that
# way; _earliest did not.
PLACEHOLDER_DATE_PREFIX = "1970"


def _earliest(a: str | None, b: str | None) -> str | None:
    """Earliest REAL date, falling back to the placeholder only if that is all
    there is.

    A plain min() over ISO strings makes "1970-01-01" win every comparison,
    so merging a dated edge with an undated witness turned a dated fact into
    an undated one. temporal.edge_year then returns None for it, and under the
    default undated_policy="include" it is kept in EVERY training graph
    regardless of its real date -- a temporal leak created by the merge, in
    the module whose docstring says the earliest date is kept precisely so the
    split can use when a fact was first established.
    """
    dates = [d for d in (a, b) if d]
    if not dates:
        return None
    real = [d for d in dates if not d.startswith(PLACEHOLDER_DATE_PREFIX)]
    return min(real) if real else min(dates)


def merge_edges(a: dict, b: dict, derive: DeriveFn | None = None) -> dict:
    """One fact, two witnesses.

    support_count records HOW MANY independent assertions back the fact. A
    (protein, virus) dependency pair supported by fifteen CRISPR screens is
    much stronger evidence than one supported by a single screen, and without
    this the two are indistinguishable after merging -- so DWPC would weight a
    lone noisy hit exactly like ACE2.
    """
    out = dict(a)
    out["support_count"] = a.get("support_count", 1) + b.get("support_count", 1)
    sources = {a.get("primary_knowledge_source"), b.get("primary_knowledge_source")}
    sources.discard(None)
    out["primary_knowledge_source"] = "|".join(sorted(sources))
    pubs = set(a.get("publications") or []) | set(b.get("publications") or [])
    if pubs:
        out["publications"] = sorted(pubs)
    # Strongest evidence wins: tier 1 beats tier 4.
    out["evidence_tier"] = min(a.get("evidence_tier", 4), b.get("evidence_tier", 4))
    out["first_asserted_date"] = _earliest(a.get("first_asserted_date"),
                                           b.get("first_asserted_date"))
    # PRECEDENCE DECIDES CONFLICTS, BUT IT CANNOT DECIDE A DERIVED FIELD.
    # `a` wins here, which is right for a measured value: the earlier layer is
    # the more authoritative source. It was catastrophic for `unquantified`,
    # where chembl's "no CC50 exists" beat the selectivity layer's CC50 on
    # 1,264 edges -- the losing layer was the one holding the evidence. Any
    # qualifier the schema declares as derived is recomputed from the merged
    # qualifiers, so no layer's local view can survive the merge.
    merged_quals = dict(b.get("qualifiers") or {})
    merged_quals.update(a.get("qualifiers") or {})
    out["qualifiers"] = derive(a["predicate"], merged_quals) if derive else merged_quals
    return out


class Assembly:
    def __init__(self, identity: dict[str, list[str]] | None = None,
                 derive: DeriveFn | None = None,
                 source_identity: set[str] | None = None) -> None:
        self.identity = identity or {}
        self.source_identity = source_identity or set()
        # Without a derive function the assembled graph carries whatever the
        # layers wrote. That is why the driver always passes one: a derived
        # qualifier stated by a layer is a claim no layer can substantiate.
        self.derive = derive
        self.nodes: dict[str, dict] = {}
        self.node_layer: dict[str, str] = {}
        self.edges: dict[tuple, dict] = {}
        self.stats: Counter = Counter()
        self.conflicts: list[str] = []
        # id -> the key it was derived from, so a collision is DETECTED rather
        # than assumed away. Two facts sharing an id would make the agent's
        # citations ambiguous in a way no downstream check could see.
        self.fact_ids: dict[str, tuple] = {}

    def add_layer(self, name: str, release_dir: Path) -> Counter:
        local: Counter = Counter()
        for n in _read_jsonl(Path(release_dir) / "nodes.jsonl"):
            local["nodes"] += 1
            nid = n["id"]
            if nid not in self.nodes:
                self.nodes[nid] = n
                self.node_layer[nid] = name
                continue
            existing = self.nodes[nid]
            if existing["class"] != n["class"]:
                self.conflicts.append(
                    f"{nid}: class {existing['class']} ({self.node_layer[nid]}) "
                    f"vs {n['class']} ({name})")
                self.stats["class_conflict"] += 1
                continue
            self._merge_properties(nid, existing, n, name)
            local["node_merged"] += 1

        for e in _read_jsonl(Path(release_dir) / "edges.jsonl"):
            local["edges"] += 1
            k = edge_key(e, self.identity, self.source_identity)
            fid = fact_id(k)
            seen_key = self.fact_ids.setdefault(fid, k)
            if seen_key != k:
                self.conflicts.append(
                    f"fact id {fid} derived from two different keys: "
                    f"{seen_key} and {k}")
                self.stats["fact_id_collision"] += 1
            e["edge_id"] = fid
            if k in self.edges:
                self.edges[k] = merge_edges(self.edges[k], e, self.derive)
                self.stats["edge_merged"] += 1
            else:
                if self.derive:
                    # Unmerged edges get it too: 304 of the selectivity
                    # layer's 1,568 edges had no chembl edge to merge with,
                    # and they are no more entitled to state the field.
                    e["qualifiers"] = self.derive(e["predicate"],
                                                  e.get("qualifiers") or {})
                self.edges[k] = e
        return local

    def _merge_properties(self, nid: str, existing: dict, incoming: dict, layer: str) -> None:
        """Fill gaps from the later layer; never overwrite, report real clashes."""
        ep, ip = existing.setdefault("properties", {}), incoming.get("properties") or {}
        for k, v in ip.items():
            if v is None:
                continue
            if k not in ep or ep[k] is None:
                ep[k] = v
            elif ep[k] != v:
                self.stats[f"property_conflict:{k}"] += 1
                if len(self.conflicts) < 40:
                    self.conflicts.append(
                        f"{nid}.{k}: {str(ep[k])[:34]!r} ({self.node_layer[nid]}) "
                        f"vs {str(v)[:34]!r} ({layer})")

    def dangling(self) -> list[tuple]:
        return [k for k in self.edges
                if k[0] not in self.nodes or k[2] not in self.nodes]

    def degree(self) -> Counter:
        d: Counter = Counter()
        # Keys may carry identity qualifiers beyond (subject, predicate,
        # object), so index rather than unpack.
        for k in self.edges:
            d[k[0]] += 1
            d[k[2]] += 1
        return d

    def hub_threshold(self, percentile: float) -> tuple[int, list[str]]:
        """Degree cut at the given percentile, plus the nodes above it."""
        deg = self.degree()
        if not deg:
            return 0, []
        ordered = sorted(deg.values())
        idx = min(int(len(ordered) * percentile / 100), len(ordered) - 1)
        cut = ordered[idx]
        return cut, [n for n, k in deg.items() if k > cut]

    def by_class(self) -> Counter:
        return Counter(n["class"] for n in self.nodes.values())

    def by_predicate(self) -> Counter:
        return Counter(k[1] for k in self.edges)

    def write(self, out_dir: Path) -> None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "nodes.jsonl").write_text(
            "\n".join(json.dumps(n) for n in self.nodes.values()))
        (out_dir / "edges.jsonl").write_text(
            "\n".join(json.dumps(e) for e in self.edges.values()))

    def manifest(self, layers: dict[str, Any], schema_version: str) -> dict:
        deg = self.degree()
        ordered = sorted(deg.values()) or [0]
        return {
            "schema_version": schema_version,
            "layers": layers,
            "nodes": len(self.nodes),
            "edges": len(self.edges),
            "nodes_by_class": dict(self.by_class()),
            "edges_by_predicate": dict(self.by_predicate()),
            "degree": {
                "median": ordered[len(ordered) // 2],
                "p99": ordered[int(len(ordered) * 0.99)],
                "p99_9": ordered[int(len(ordered) * 0.999)],
                "max": ordered[-1],
            },
            "merges": {k: v for k, v in self.stats.items()},
        }


def assemble(layers: dict[str, Path],
             identity: dict[str, list[str]] | None = None,
             derive: DeriveFn | None = None,
             source_identity: set[str] | None = None) -> Assembly:
    """Layers are added in LAYER_PRECEDENCE order regardless of dict order."""
    a = Assembly(identity, derive, source_identity)
    for name in LAYER_PRECEDENCE:
        if name in layers:
            a.add_layer(name, layers[name])
    for name, path in layers.items():
        if name not in LAYER_PRECEDENCE:
            a.add_layer(name, path)
    return a


def orphan_nodes(a: Assembly) -> Counter:
    """Nodes with no edges at all, counted by class.

    Not an error -- the host proteome contributes many proteins with no
    coronavirus evidence -- but a large orphan count in a class that should be
    connected (SmallMolecule, for instance) means a join failed.
    """
    touched: set[str] = set()
    for k in a.edges:
        touched.add(k[0])
        touched.add(k[2])
    return Counter(n["class"] for nid, n in a.nodes.items() if nid not in touched)


def connectivity_report(a: Assembly) -> dict[str, Counter]:
    """Per-class edge-type breakdown, for spotting a layer that failed to join."""
    out: dict[str, Counter] = defaultdict(Counter)
    for k in a.edges:
        sc = a.nodes.get(k[0], {}).get("class", "?")
        oc = a.nodes.get(k[2], {}).get("class", "?")
        out[k[1]][f"{sc}->{oc}"] += 1
    return out
