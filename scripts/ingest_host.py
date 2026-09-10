#!/usr/bin/env python3
"""Day 4 driver: build the host layer release."""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.emit import Emit
from kgav.host import (
    EXPERIMENTAL_MIN,
    build_ensp_map,
    ingest_proteome,
    ingest_string,
    load_human_proteome,
)

PROV = {
    "proteome_source": "infores:uniprot",
    "pathway_source": "infores:uniprot",
    "ppi_source": "infores:string",
    # STRING v12.0 human links file release date, from the download timestamp.
    "string_release_date": "2023-08-28",
}


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--host-dir", type=Path, default=root / "data" / "raw" / "host")
    ap.add_argument("--out", type=Path, default=root / "data" / "releases" / "v0.1-host")
    ap.add_argument("--experimental-min", type=int, default=EXPERIMENTAL_MIN)
    args = ap.parse_args()

    em = Emit()

    print("loading human proteome")
    entries = load_human_proteome(args.host_dir / "human_proteome.json.gz")
    print(f"  {len(entries):,} reviewed entries")

    s = ingest_proteome(em, entries, PROV)
    print(f"  proteins={s['proteins']:,} genes={s['genes']:,} "
          f"reactome={s['reactome_edges']:,} go={s['go_edges']:,}")

    accessions = {n["id"].split(":", 1)[1] for n in em.nodes.values() if n["class"] == "Protein"}
    print("\nmapping STRING identifiers")
    mapping = build_ensp_map(args.host_dir / "9606.protein.aliases.v12.0.txt.gz",
                             accessions, em.notes)
    print(f"  {len(mapping):,} ENSP -> UniProt, "
          f"{em.notes['ensp_ambiguous_after_filter']:,} ambiguous and dropped")

    print(f"\nstreaming STRING links (experimental >= {args.experimental_min})")
    t = ingest_string(em, args.host_dir / "9606.protein.links.detailed.v12.0.txt.gz",
                      mapping, PROV, args.experimental_min)
    print(f"  {t['rows']:,} rows -> {t['passed_filter']:,} passed filter -> {t['edges']:,} edges")
    print(f"  dropped: {t['unmapped_endpoint']:,} unmapped, "
          f"{t['reciprocal_duplicate']:,} reciprocal, {t['self_loop']:,} self-loop")

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in em.nodes.values()))
    (args.out / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in em.edges))

    print(f"\n{len(em.nodes):,} nodes, {len(em.edges):,} edges -> {args.out}")
    for c, n in Counter(x["class"] for x in em.nodes.values()).most_common():
        print(f"  {c:<16} {n:>8,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
