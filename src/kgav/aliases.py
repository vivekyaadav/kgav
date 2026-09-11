"""Identifier alias resolution.

The rule learned three times over during the spine build: an identifier the
graph chooses not to keep as its own node must still RESOLVE. Sources cite the
ids they cite. A dropped identifier fails silently -- no error, no missing-data
signal, just edges that never join.

Every ingest downstream of the spine loads the SAME_AS edges and resolves
through them before joining on any protein id.
"""
from __future__ import annotations

import json
from pathlib import Path


class AliasIndex:
    """Resolve an id to its canonical node, following SAME_AS chains."""

    def __init__(self) -> None:
        self.alias: dict[str, str] = {}
        self.nodes: dict[str, dict] = {}

    @classmethod
    def from_releases(cls, *release_dirs: Path) -> AliasIndex:
        ix = cls()
        for d in release_dirs:
            nf, ef = Path(d) / "nodes.jsonl", Path(d) / "edges.jsonl"
            if nf.exists():
                for line in nf.read_text().splitlines():
                    if line.strip():
                        n = json.loads(line)
                        ix.nodes[n["id"]] = n
            if ef.exists():
                for line in ef.read_text().splitlines():
                    if not line.strip():
                        continue
                    e = json.loads(line)
                    if e["predicate"] == "SAME_AS":
                        ix.alias[e["subject"]] = e["object"]
        return ix

    def resolve(self, node_id: str, *, max_hops: int = 8) -> str | None:
        """Canonical id, or None if it does not exist in the loaded releases.

        Follows SAME_AS chains with a hop cap: a cycle in the alias graph is a
        data bug, and looping forever hides it. Returning None on an unresolved
        id is deliberate -- callers must count misses, not silently skip them.
        """
        seen = {node_id}
        cur = node_id
        for _ in range(max_hops):
            nxt = self.alias.get(cur)
            if nxt is None:
                break
            if nxt in seen:
                return None          # cycle
            seen.add(nxt)
            cur = nxt
        return cur if cur in self.nodes else None

    def klass(self, node_id: str) -> str | None:
        n = self.nodes.get(node_id)
        return n["class"] if n else None

    def is_viral(self, node_id: str) -> bool | None:
        n = self.nodes.get(node_id)
        if not n or n["class"] != "Protein":
            return None
        return n["properties"].get("is_viral")

    def is_alias(self, node_id: str) -> bool:
        return node_id in self.alias

    def __len__(self) -> int:
        return len(self.nodes)
