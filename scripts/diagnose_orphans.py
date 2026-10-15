#!/usr/bin/env python3
"""Name every VIRAL orphan in a release, and say which failure it is.

WHY THIS EXISTS. assemble_release has always printed "Protein 634" under
"nodes with no edges", with a sentence naming two very different causes in
one breath: the host proteome is ingested whole so host orphans are
expected, but "a large orphan count in a class that should be connected
means a join failed". Splitting on is_viral turned that into:

    Protein               634
      host                                633
      viral NCBITaxon:2697049               1   <- a join failed

One SARS-CoV-2 protein, in the graph, touched by nothing. That is a real
defect and the split is as far as a per-class counter can take it. Naming
the node, and separating WHY it is unreferenced, needs the layers.

THREE DISTINCT CAUSES, and the fix differs for each:

  DROPPED AT ASSEMBLY  some layer emitted an edge naming this id and the
                       assembled graph has none. Precedence, merging or
                       an edge-key collision ate it. The layer is the
                       evidence and the bug is in assembly.

  DISCARDED IDENTIFIER some other node carries the same mature peptide or
                       sequence hash AND has edges. The evidence went to a
                       sibling id. This repo's README calls this "the
                       dominant bug class" and counts four instances.

  NEVER REFERENCED     no layer emitted any edge for it and no sibling
                       holds its evidence. The node is genuinely isolated,
                       and the question becomes whether it should exist --
                       the spine writes ENCODED_BY only when a chain has a
                       gene symbol, so a symbol-less chain with no PPI and
                       no inhibitor orphans legitimately.

Reports which, per node. Decides nothing.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.assemble import LAYER_PRECEDENCE  # noqa: E402
from kgav.schema import print_capped  # noqa: E402


def read_jsonl(path: Path) -> list[dict]:
    if not path.exists():
        return []
    return [json.loads(ln) for ln in path.read_text().splitlines() if ln.strip()]


def untouched(nodes: list[dict], edges: list[dict]) -> list[dict]:
    """Nodes no assembled edge names, as subject or object."""
    seen: set[str] = set()
    for e in edges:
        seen.add(e["subject"])
        seen.add(e["object"])
    return [n for n in nodes if n["id"] not in seen]


def layer_references(releases: Path, prefix: str = "v0.1-") -> dict[str, Counter]:
    """id -> Counter of layer name, over every layer's OWN edges.

    Read from the per-layer files rather than the assembled graph, because
    the whole question is whether something was present before assembly and
    absent after.
    """
    out: dict[str, Counter] = defaultdict(Counter)
    for name in LAYER_PRECEDENCE:
        d = releases / f"{prefix}{name}"
        for e in read_jsonl(d / "edges.jsonl"):
            out[e["subject"]][name] += 1
            out[e["object"]][name] += 1
    return out


#: Identity keys, strongest first, with what a match actually proves.
#:
#: sequence_hash IS NOT ALWAYS A HASH OF THE SEQUENCE. ingest_spine builds
#: it two ways: _sha(seq) on an entry-level protein, and _sha(f"{acc}:{cid}")
#: on a CHAIN. The second is a hash of the identifier, so it is injective by
#: construction and two distinct chains can never share it. A chain match is
#: therefore real evidence; the absence of one proves nothing.
#:
#: mature_peptide is a human-readable chain description, not an identity.
IDENTITY_KEYS = {
    "sequence_hash": "same sequence hash",
    "mature_peptide": "same chain label",
}


def siblings(node: dict, nodes: list[dict], connected: set[str]) -> list[tuple]:
    """Other nodes of the same class AND TAXON sharing an identity key.

    THE TAXON CLAUSE IS THE WHOLE CORRECTION. Without it this reported
    SARS-CoV-2's orphaned pp1a as a DISCARDED IDENTIFIER on the strength of
    six matches -- every one of them a different coronavirus whose pp1a
    chain carries the same label, "Replicase polyprotein 1a". ingest_spine
    dedups chain labels in `seen_chains`, which is LOCAL to ingest_virus, so
    the same label across viruses is expected and means nothing. Six
    different proteins, reported as one filed twice.

    A label is only identity-bearing inside one taxon, and there the spine's
    own dedup should already have collapsed it -- so a same-taxon label
    match is itself a finding about that dedup.
    """
    props = node.get("properties") or {}
    tax = props.get("taxon_id")
    keys = [k for k in IDENTITY_KEYS if props.get(k)]
    found = []
    for other in nodes:
        if other["id"] == node["id"] or other.get("class") != node.get("class"):
            continue
        op = other.get("properties") or {}
        if op.get("taxon_id") != tax:
            continue
        for k in keys:
            if op.get(k) == props[k]:
                found.append((k, other["id"], other["id"] in connected))
    return found


def classify(node: dict, refs: Counter, sibs: list[tuple]) -> str:
    if refs:
        return "DROPPED AT ASSEMBLY"
    if any(connected for _k, _i, connected in sibs):
        return "DISCARDED IDENTIFIER"
    return "NEVER REFERENCED"


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", type=Path, default=root / "data/releases/v0.1")
    ap.add_argument("--releases", type=Path, default=root / "data/releases",
                    help="directory holding the per-layer v0.1-* releases")
    ap.add_argument("--all", action="store_true",
                    help="report host orphans too, not just viral ones")
    args = ap.parse_args()

    nodes_f = args.release / "nodes.jsonl"
    if not nodes_f.exists():
        print(f"no assembled release at {args.release}; run "
              f"scripts/assemble_release.py", file=sys.stderr)
        return 2
    nodes = read_jsonl(nodes_f)
    edges = read_jsonl(args.release / "edges.jsonl")
    orphans = untouched(nodes, edges)
    connected = {n["id"] for n in nodes} - {n["id"] for n in orphans}

    by_kind: Counter = Counter()
    for n in orphans:
        p = n.get("properties") or {}
        by_kind["viral" if p.get("is_viral") else "host/other"] += 1
    print(f"{len(nodes):,} nodes, {len(edges):,} edges, "
          f"{len(orphans):,} with no edge "
          f"({by_kind['viral']} viral, {by_kind['host/other']} host/other)")

    refs = layer_references(args.releases)
    if not refs:
        print(f"\nno per-layer releases found under {args.releases}; the "
              f"DROPPED-AT-ASSEMBLY case cannot be distinguished from "
              f"NEVER REFERENCED without them.", file=sys.stderr)

    targets = orphans if args.all else [
        n for n in orphans if (n.get("properties") or {}).get("is_viral")]
    if not targets:
        print("\nNo viral orphan. Every orphan above is the expected case: "
              "the host\nproteome is ingested whole and carries entries with "
              "no viral evidence.")
        return 0

    print(f"\n{len(targets)} orphan(s) to account for:\n")
    verdicts: Counter = Counter()
    for n in sorted(targets, key=lambda n: n["id"]):
        props = n.get("properties") or {}
        sibs = siblings(n, nodes, connected)
        verdict = classify(n, refs.get(n["id"], Counter()), sibs)
        verdicts[verdict] += 1
        print(f"  {n['id']}   [{n.get('class')}]   {verdict}")
        for k, v in sorted(props.items()):
            print(f"      {k:<18} {v}")
        r = refs.get(n["id"], Counter())
        if r:
            print("      referenced by these layers' own edges, and by no "
                  "assembled edge:")
            for layer, cnt in r.most_common():
                print(f"        {layer:<14} {cnt:>6,} edge(s)")
            print("      -> assembly dropped them. Precedence, merging or an "
                  "edge-key collision.")
        elif sibs:
            print("      another node of the SAME TAXON carries the same "
                  "identity property:")
            for k, oid, conn in sibs:
                print(f"        {IDENTITY_KEYS[k]:<20} {oid}   "
                      f"{'HAS edges' if conn else 'also orphaned'}")
            print("      -> the evidence is filed under a second id for the "
                  "same protein.")
        else:
            why = ("no gene_symbol, so the spine writes no ENCODED_BY"
                   if not props.get("gene_symbol")
                   else "has a gene_symbol, so even ENCODED_BY is missing")
            print(f"      no layer emitted any edge naming it; {why}.")
            if str(n["id"]).startswith("UniProtKB:PRO_"):
                print("      NOTE this is a chain node, whose sequence_hash "
                      "is _sha(accession:chain_id)\n      -- a hash of the "
                      "identifier. It cannot match another chain's, so it "
                      "rules\n      nothing out. Only the label was "
                      "available to compare, and a label is\n      not "
                      "identity.")
        print()

    print("verdicts:")
    for v, n in verdicts.most_common():
        print(f"  {v:<22} {n}")
    print("\nDROPPED AT ASSEMBLY is a bug in assembly. DISCARDED IDENTIFIER "
          "is a bug in\nthe ingest that minted the second id. NEVER "
          "REFERENCED is a curation question:\nshould the spine build a "
          "chain nothing can connect to?")

    if not args.all and by_kind["host/other"]:
        print(f"\n{by_kind['host/other']:,} host/other orphans not listed "
              f"(--all to see them). A sample:")
        print_capped([n["id"] for n in orphans
                      if not (n.get("properties") or {}).get("is_viral")],
                     5, unit="more host/other orphan(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
