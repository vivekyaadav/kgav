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



def coverage_lines(per_virus, taxon_map: dict, nodes: dict) -> list[str]:
    """Edges per virus over EVERY declared-and-built virus, zeros included.

    This was per_virus.most_common(), so a virus with zero edges did not
    appear at all -- and the comment in main records that exactly this
    happened: no flavivirus file was ever fetched, Zika sat in the taxon
    map producing nothing, and the build said so nowhere. The absence was
    discoverable only by noticing which labels were missing from a list
    that never claimed to be complete.

    Same correction pathway_landscape already carries for its own coverage
    table. A zero is a finding about the SOURCE rather than about the
    virus, and it is the finding that decides what to fetch next.

    A pure function returning lines, so the reporting can be tested without
    a spine release on disk.
    """
    declared = sorted(set(taxon_map.values()))
    built = [t for t in declared if f"NCBITaxon:{t}" in nodes]

    def label(tid):
        return ((nodes.get(tid, {}) or {}).get("properties") or {}).get(
            "label", tid)

    header = ("edges per virus, over every virus the register declares "
              "and this release built:")
    out = ["", header]
    zero = 0
    for taxon in built:
        tid = f"NCBITaxon:{taxon}"
        n = per_virus.get(tid, 0)
        if not n:
            zero += 1
        flag = "   <- NO virus-host PPI from either source" if not n else ""
        out.append(f"  {label(tid):<14} {n:>7,}{flag}")
    if zero:
        out.append("")
        out.append(f"  {zero} of {len(built)} built viruses have NO "
                   f"virus-host PPI edge from VirHostNet or IntAct. That is "
                   f"a gap in what was FETCHED as much as in what was "
                   f"curated: the vhn_*.tab27 files are one coronavirus "
                   f"each, and a taxon with no file cannot produce an edge "
                   f"however well the source covers it.")
    unbuilt = len(declared) - len(built)
    if unbuilt:
        out.append(f"  {unbuilt} further declared taxa are not in this "
                   f"release and were not counted.")
    return out


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
    for line in coverage_lines(per_virus, taxon_map, index.nodes):
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
