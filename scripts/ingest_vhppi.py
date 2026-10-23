#!/usr/bin/env python3
"""Day 5 driver: virus-host PPI bridge layer."""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.aliases import AliasIndex
from kgav.emit import Emit
from kgav.virus_register import declared_taxa
from kgav.vhppi import ingest_mitab


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=root / "config" / "viruses.yaml")
    ap.add_argument("--vhppi-dir", type=Path, default=root / "data" / "raw" / "vhppi")
    ap.add_argument("--spine", type=Path, default=root / "data" / "releases" / "v0.1-spine")
    ap.add_argument("--host", type=Path, default=root / "data" / "releases" / "v0.1-host")
    ap.add_argument("--out", type=Path, default=root / "data" / "releases" / "v0.1-vhppi")
    ap.add_argument("--date", default="2024-01-01",
                    help="fallback first_asserted_date when a row has no PMID year")
    args = ap.parse_args()

    print("loading spine + host as the resolution index")
    index = AliasIndex.from_releases(args.spine, args.host)
    print(f"  {len(index):,} nodes, {len(index.alias):,} aliases")

    # ONE MAP, FROM THE REGISTER. virus_register.declared_taxa exists for
    # exactly this and its docstring says so -- "the same shape the ingests
    # build ad hoc from isolate_taxa, derived here so chain_source cannot
    # drift out of the map" -- and nothing called it: four scripts built the
    # map by hand and the canonical function was reached only by its own test.
    # Five implementations of one mapping, which is how AUC 0.826, the
    # coverage floor and lift's trustworthy rule each ended up with two.
    #
    # declared_taxa also folds in chain_source.taxon, which the hand-built
    # copies omitted. Identical on the current register -- verified before
    # this change -- because validate()'s UNDECLARED_CHAIN_TAXON already
    # forces every differing chain taxon into isolate_taxa. So this is
    # behaviour-neutral today and stops being luck tomorrow.
    taxon_map = declared_taxa(
        yaml.safe_load(args.config.read_text()))
    print(f"  {len(taxon_map)} viral taxids mapped to "
          f"{len(set(taxon_map.values()))} species")

    # TWO SOURCES, ONE LAYER. This globbed vhn_*.tab27 only, and
    # data/raw/vhppi also held intact.zip -- 1.35 GB, dated 2026-01-14,
    # never read, because the pattern did not match it. The seven vhn files
    # are one coronavirus each and no flavivirus file was ever fetched,
    # which is why Zika produced zero edges while being in the taxon map
    # since 2026-10-12. The broader source sat beside the narrower one and
    # nothing said so.
    #
    # They share a predicate and a layer, and assemble_release already
    # merges one fact from two sources into one edge carrying both
    # ("infores:string|infores:virhostnet"). What must NOT be shared is the
    # label: an IntAct row recorded as infores:virhostnet is a false
    # provenance claim, and the tier ablation reads these labels.
    groups = [("infores:virhostnet", sorted(args.vhppi_dir.glob("vhn_*.tab27"))),
              ("infores:intact", sorted(args.vhppi_dir.glob("intact_*.tab27")))]
    files = [f for _src, fs in groups for f in fs]
    if not files:
        raise SystemExit(
            f"no vhn_*.tab27 or intact_*.tab27 files in {args.vhppi_dir}.\n"
            f"If intact.zip is present, scripts/filter_intact.py cuts the "
            f"viral rows out of it.")
    print(f"\nparsing {len(files)} MITAB files")
    for src, fs in groups:
        if fs:
            print(f"  {src:<24} {len(fs)} file(s): "
                  f"{', '.join(f.name for f in fs)}")

    em = Emit()
    s = Counter()
    for src, fs in groups:
        if not fs:
            continue
        part = ingest_mitab(em, fs, index, taxon_map,
                            source=src, default_date=args.date)
        print(f"\n  {src}: {part['rows']:,} rows -> {part['edges']:,} edges")
        s += part

    print(f"  {s['rows']:,} rows -> {s['edges']:,} edges")
    print(f"  {s['direct']:,} direct interaction, {s['co_complex']:,} co-complex/"
          f"association ({s['binary_method']:,} from a binary assay)")
    print("\n  dropped:")
    for k in ("unparseable", "no_detection_method", "unresolved_id", "self_interaction",
              "not_a_protein", "host_host", "viral_viral", "unmapped_viral_taxon",
              "duplicate"):
        if s[k]:
            print(f"    {k:<24} {s[k]:>8,}")

    # This release is edges-only: every endpoint lives in the spine or host
    # release. validate_release.py would report all of them DANGLING against an
    # empty node file, so validation runs here against the real index instead.
    from kgav.schema import load_schema, report_violations
    violations = load_schema().validate_batch(index.nodes.values(), em.edges)
    if report_violations(violations):
        return 1

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "nodes.jsonl").write_text("")
    (args.out / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in em.edges))
    print(f"\n{len(em.edges):,} edges -> {args.out}")

    per_virus: Counter = Counter()
    for e in em.edges:
        n = index.nodes.get(e["object"])
        if n:
            per_virus[n["properties"].get("taxon_id", "?")] += 1
    print("\nedges per virus:")
    for tid, n in per_virus.most_common():
        label = (index.nodes.get(tid, {}).get("properties", {}) or {}).get("label", tid)
        print(f"  {label:<14} {n:>7,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
