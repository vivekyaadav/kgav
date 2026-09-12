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

import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

# Later layers do not overwrite earlier ones on conflict. The spine is
# hand-curated and authoritative for viral entities; the host layer is derived
# from reference proteomes; chembl mints compounds and knows least about
# anything it did not create.
LAYER_PRECEDENCE = ["spine", "host", "vhppi", "orcs", "chembl"]


def _read_jsonl(path: Path) -> list[dict]:
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def edge_key(e: dict, identity: dict[str, list[str]] | None = None) -> tuple:
    """Identity of a fact.

    (subject, predicate, object) alone is WRONG for predicates with a
    direction-bearing qualifier: a dependency edge and a restriction edge for
    the same (protein, virus) are different facts, not two witnesses to one.
    320 pairs in the ORCS layer carry both, and merging on the triple alone let
    whichever source was read first decide the direction.
    """
    base = (e["subject"], e["predicate"], e["object"])
    if not identity:
        return base
    quals = identity.get(e["predicate"]) or []
    if not quals:
        return base
    eq = e.get("qualifiers") or {}
    return base + tuple(eq.get(q) for q in quals)


def _earliest(a: str | None, b: str | None) -> str | None:
    dates = [d for d in (a, b) if d]
    return min(dates) if dates else None


def merge_edges(a: dict, b: dict) -> dict:
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
    merged_quals = dict(b.get("qualifiers") or {})
    merged_quals.update(a.get("qualifiers") or {})
    out["qualifiers"] = merged_quals
    return out


class Assembly:
    def __init__(self, identity: dict[str, list[str]] | None = None) -> None:
        self.identity = identity or {}
        self.nodes: dict[str, dict] = {}
        self.node_layer: dict[str, str] = {}
        self.edges: dict[tuple, dict] = {}
        self.stats: Counter = Counter()
        self.conflicts: list[str] = []

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
            k = edge_key(e, self.identity)
            if k in self.edges:
                self.edges[k] = merge_edges(self.edges[k], e)
                self.stats["edge_merged"] += 1
            else:
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
             identity: dict[str, list[str]] | None = None) -> Assembly:
    """Layers are added in LAYER_PRECEDENCE order regardless of dict order."""
    a = Assembly(identity)
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
