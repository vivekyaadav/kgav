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
from kgav.vhppi import ingest_mitab


def build_taxon_map(config: Path) -> dict[str, str]:
    """Every viral taxid that may appear -> the species taxid the graph uses.

    VirHostNet files MERS under the HCoV-EMC/2012 isolate and HKU1 under a
    strain taxon. Querying the species node returns zero rows for both. Without
    this map those interactions look like an unknown organism and vanish.
    """
    cfg = yaml.safe_load(config.read_text())
    m: dict[str, str] = {}
    for v in cfg["viruses"]:
        m[str(v["taxon"])] = str(v["taxon"])
        for iso in v.get("isolate_taxa") or []:
            m[str(iso)] = str(v["taxon"])
    return m


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

    taxon_map = build_taxon_map(args.config)
    print(f"  {len(taxon_map)} viral taxids mapped to "
          f"{len(set(taxon_map.values()))} species")

    files = sorted(args.vhppi_dir.glob("vhn_*.tab27"))
    if not files:
        raise SystemExit(f"no vhn_*.tab27 files in {args.vhppi_dir}")
    print(f"\nparsing {len(files)} MITAB files")

    em = Emit()
    s = ingest_mitab(em, files, index, taxon_map,
                     source="infores:virhostnet", default_date=args.date)

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
    from kgav.schema import load_schema
    violations = load_schema().validate_batch(index.nodes.values(), em.edges)
    if violations:
        print(f"\nSCHEMA FAIL — {len(violations)} violations")
        for v in violations[:8]:
            print(f"  {v}")
        return 1
    print("\nschema: PASS")

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
