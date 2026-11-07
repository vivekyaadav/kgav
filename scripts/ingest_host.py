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
    # WHAT THE GO ASPECT FILTER COSTS (fix-list E7). Only P: is a pathway,
    # and F:/C: annotations were discarded with nothing recorded -- so the
    # share of UniProt's GO data this layer uses was unknown, and a release
    # that stopped prefixing the aspect would drop all of it invisibly.
    aspects = {k.removeprefix("go_aspect_"): v for k, v in s.items()
               if k.startswith("go_aspect_")}
    if aspects:
        census = ", ".join(f"{k}={v:,}" for k, v in sorted(aspects.items()))
        print(f"  GO annotations by aspect: {census}")
        print(f"    {s['go_dropped_not_a_process']:,} dropped as not a "
              f"biological process")
        if s["go_aspect_unprefixed"]:
            print(f"    ! {s['go_aspect_unprefixed']:,} carry NO aspect "
                  f"prefix; the reader expects 'P:'/'F:'/'C:' and this "
                  f"release does not use it")
    # The two halves of the Gene gate, apart: a missing xref is an
    # ID-mapping gap, a missing name is a TrEMBL entry.
    if s["no_gene_id_xref"] or s["no_gene_symbol"]:
        print(f"  no Gene node for {s['no_gene_id_xref']:,} protein(s) with "
              f"no GeneID xref and {s['no_gene_symbol']:,} with no symbol")

    accessions = {n["id"].split(":", 1)[1] for n in em.nodes.values() if n["class"] == "Protein"}
    print("\nmapping STRING identifiers")
    mapping = build_ensp_map(args.host_dir / "9606.protein.aliases.v12.0.txt.gz",
                             accessions, em.notes)
    n = em.notes
    print(f"  {len(mapping):,} ENSP -> UniProt, "
          f"{n['ensp_ambiguous_after_filter']:,} ambiguous and dropped")
    print(f"  alias rows read {n['alias_rows']:,}:")
    for k, label in (("alias_kept", "kept"),
                     ("alias_dropped_other_source", "not a UniProt_AC row"),
                     ("alias_dropped_outside_proteome",
                      "accession outside the reviewed proteome"),
                     ("alias_dropped_short_row", "fewer than 3 columns")):
        if n[k]:
            print(f"    {n[k]:>9,}  {label}")
    acc = (n["alias_kept"] + n["alias_dropped_other_source"]
           + n["alias_dropped_outside_proteome"] + n["alias_dropped_short_row"])
    if acc != n["alias_rows"]:
        print(f"    ! ACCOUNTING DOES NOT CLOSE: {n['alias_rows']:,} read, "
              f"{acc:,} accounted for")
        return 1

    print(f"\nstreaming STRING links (experimental >= {args.experimental_min})")
    t = ingest_string(em, args.host_dir / "9606.protein.links.detailed.v12.0.txt.gz",
                      mapping, PROV, args.experimental_min)
    print(f"  {t['rows']:,} rows -> {t['passed_filter']:,} passed filter -> {t['edges']:,} edges")
    print(f"  dropped: {t['below_experimental_min']:,} below "
          f"experimental {args.experimental_min}, "
          f"{t['unmapped_endpoint']:,} unmapped, "
          f"{t['reciprocal_duplicate']:,} reciprocal, "
          f"{t['self_loop']:,} self-loop")
    # THE CUTOFF DECIDES THE SIZE OF THREE QUARTERS OF THE GRAPH, and was
    # only derivable here as rows - passed_filter. src/kgav/host.py's
    # docstring justifies it with four row counts that no code in this
    # repository computes; this is the one figure that is actually measured.
    link_acc = (t["below_experimental_min"] + t["unmapped_endpoint"]
                + t["self_loop"] + t["reciprocal_duplicate"] + t["edges"])
    if link_acc != t["rows"]:
        print(f"  ! ACCOUNTING DOES NOT CLOSE: {t['rows']:,} rows read, "
              f"{link_acc:,} accounted for. A row is being dropped by a "
              f"gate with no counter.")
        return 1
    print(f"  accounting closes: every one of {t['rows']:,} link rows is "
          f"either dropped by a named gate or emitted")

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in em.nodes.values()))
    (args.out / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in em.edges))

    print(f"\n{len(em.nodes):,} nodes, {len(em.edges):,} edges -> {args.out}")
    for c, n in Counter(x["class"] for x in em.nodes.values()).most_common():
        print(f"  {c:<16} {n:>8,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
