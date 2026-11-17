#!/usr/bin/env python3
"""What this release IS, in the figures a later rebuild will move.

WHY THIS EXISTS. Fix-list L26 records the cost of not having it. The L27
rebuild moved the drugged proteome from 1,946 to 1,939 human proteins, the
superseded release had been rebuilt IN PLACE, and so *which seven proteins
left* was not answerable from this repository at all -- the old operand
survived only as a number written in prose. A rebuild that overwrites its
own input destroys the only thing that could explain it.

So this is run BEFORE a rebuild and again AFTER, and the two are diffed.
It captures operands -- the sets and counts a study is computed over --
rather than results, which `scripts/compare_results.py` already covers.

    python scripts/snapshot_release_operands.py --out before.json
    ... rebuild ...
    python scripts/snapshot_release_operands.py --out after.json
    python scripts/snapshot_release_operands.py --diff before.json after.json

MEMBERS, NOT ONLY COUNTS, for the label sets. That is the whole lesson of
L26: a count tells you something moved and a set tells you what.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.labels import build_labels                      # noqa: E402
from kgav.pipeline import release_dir                     # noqa: E402
from kgav.schema import print_capped                      # noqa: E402

#: Label sets are built per SOURCE as well as pooled, because the
#: plate-matched cells are single-source by construction and a pooled count
#: cannot show one of them moving.
SOURCES = ("infores:chembl", "infores:ncats-opendata")


def _digest(ids) -> str:
    return hashlib.sha256("\n".join(sorted(ids)).encode()).hexdigest()


def label_cells(release: Path, sources: set[str] | None) -> dict:
    out = {}
    labels = build_labels(release, sources=sources)
    for virus in sorted(set(labels.positives) | set(labels.negatives)):
        pos = sorted(labels.positives.get(virus, set()))
        neg = sorted(labels.negatives.get(virus, set()))
        if not pos and not neg:
            continue
        out[virus] = {"n_pos": len(pos), "n_neg": len(neg),
                      "pos_sha256": _digest(pos), "neg_sha256": _digest(neg),
                      "positives": pos}
    return out


def snapshot(release: Path) -> dict:
    manifest = json.loads((release / "MANIFEST.json").read_text())
    by_pred_source: Counter = Counter()
    assay_types: Counter = Counter()
    for line in (release / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        by_pred_source[(e["predicate"], e["primary_knowledge_source"])] += 1
        at = (e.get("qualifiers") or {}).get("assay_type")
        if at:
            assay_types[at] += 1
    cells = {"pooled": label_cells(release, None)}
    for s in SOURCES:
        cells[s] = label_cells(release, {s})
    return {
        "release": manifest.get("release"),
        "schema_version": manifest.get("schema_version"),
        "nodes": manifest.get("nodes"),
        "edges": manifest.get("edges"),
        "nodes_by_class": manifest.get("nodes_by_class"),
        "edges_by_predicate": manifest.get("edges_by_predicate"),
        "layers": {k: v.get("edges_sha") for k, v in
                   (manifest.get("layers") or {}).items()},
        "edges_by_predicate_and_source": {
            f"{p}|{s}": n for (p, s), n in sorted(by_pred_source.items())},
        "assay_types": dict(sorted(assay_types.items())),
        "label_cells": cells,
    }


def _flat(d: dict, prefix: str = "") -> dict:
    """Leaves, excluding the member lists -- those are diffed as SETS."""
    out = {}
    for k, v in d.items():
        key = f"{prefix}{k}"
        if k == "positives":
            continue
        if isinstance(v, dict):
            out.update(_flat(v, key + "."))
        else:
            out[key] = v
    return out


def diff(before: dict, after: dict) -> int:
    a, b = _flat(before), _flat(after)
    moved = sorted(k for k in set(a) | set(b) if a.get(k) != b.get(k))
    print(f"{len(moved)} value(s) moved")
    print_capped(moved, 60,
                 render=lambda k: f"{k}\n        {a.get(k)!r}  ->  {b.get(k)!r}",
                 indent="      ")

    # THE PART A COUNT CANNOT GIVE: which compounds entered or left each
    # positive set. This is the answer L26 could not produce.
    print("\npositive-set membership, per cell:")
    rows = 0
    for scope in sorted(set(before["label_cells"]) | set(after["label_cells"])):
        bb = before["label_cells"].get(scope, {})
        aa = after["label_cells"].get(scope, {})
        for virus in sorted(set(bb) | set(aa)):
            p0 = set(bb.get(virus, {}).get("positives", []))
            p1 = set(aa.get(virus, {}).get("positives", []))
            if p0 == p1:
                continue
            rows += 1
            print(f"  {scope} {virus}: {len(p0):,} -> {len(p1):,}  "
                  f"(+{len(p1 - p0):,} / -{len(p0 - p1):,})")
            for label, s in (("left ", p0 - p1), ("joined", p1 - p0)):
                if s:
                    print_capped(sorted(s), 6, render=str,
                                 indent=f"      {label} ")
    if not rows:
        print("  no positive set changed membership")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--diff", nargs=2, type=Path, default=None,
                    help="two snapshot files, before and after")
    args = ap.parse_args(argv)

    if args.diff:
        return diff(json.loads(args.diff[0].read_text()),
                    json.loads(args.diff[1].read_text()))

    release = args.release or release_dir(ROOT)
    if not (release / "MANIFEST.json").exists():
        print(f"no release at {release}")
        return 2
    snap = snapshot(release)
    text = json.dumps(snap, indent=2, sort_keys=True)
    if args.out:
        args.out.write_text(text + "\n")
        print(f"wrote {args.out}")
    print(f"  release        {snap['release']}")
    print(f"  schema         {snap['schema_version']}")
    print(f"  nodes/edges    {snap['nodes']:,} / {snap['edges']:,}")
    print(f"  assay_types    {snap['assay_types']}")
    for scope, cells in snap["label_cells"].items():
        tot_p = sum(c["n_pos"] for c in cells.values())
        tot_n = sum(c["n_neg"] for c in cells.values())
        print(f"  {scope:<26} {len(cells):>2} cell(s), "
              f"{tot_p:,} positives / {tot_n:,} negatives")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
