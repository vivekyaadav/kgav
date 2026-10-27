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
from kgav.vhppi import ingest_mitab, slug
from kgav.virus_register import declared_taxa

#: Gates a zero-edge virus is reported against, worst-first. Not every gate
#: in classify_row: host_host and not_a_protein cannot explain a vhn file,
#: and no_detection_method is reported only when it actually fires.
JOIN_GATES = ("unresolved_id", "unmapped_viral_taxon", "viral_viral",
              "unparseable", "no_detection_method", "self_interaction",
              "duplicate")


def rows_per_taxon(per_file: dict, viruses: list, taxon_map: dict) -> dict:
    """Attribute each raw file's gate tally to the virus it was fetched for.

    The sidecar fetch_vhppi writes names the virus outright, so it is read
    first and the filename is the fallback for the hand-placed coronavirus
    files that predate the fetcher and have no sidecar.

    This exists because a row that dies at `unresolved_id` carries no taxon
    -- resolution runs before the taxon map -- so the file is the only thing
    that can say which virus lost it. Attribution is per FILE, which is
    exactly as precise as the fetcher's one-file-per-virus layout and no
    more: an intact_*.tab27 holding many viruses attributes to none of them,
    and is reported as unattributed rather than spread by guesswork.
    """
    by_stem = {slug(v["name"]): str(v["taxon"]) for v in viruses}
    out: dict[str, Counter] = {}
    unattributed: Counter = Counter()
    for path, c in per_file.items():
        taxon = None
        side = Path(str(path) + ".source.json")
        if side.exists():
            try:
                taxon = str(json.loads(side.read_text()).get("taxon") or "")
            except (ValueError, OSError):
                taxon = None
        if not taxon:
            stem = Path(path).stem
            taxon = by_stem.get(stem.removeprefix("vhn_"))
        taxon = taxon_map.get(taxon, taxon) if taxon else None
        if not taxon:
            unattributed += c
            continue
        out.setdefault(taxon, Counter())
        out[taxon] += c
    if unattributed:
        out["?"] = unattributed
    return out


def coverage_lines(per_virus, taxon_map: dict, nodes: dict,
                   fetched: dict | None = None) -> list[str]:
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
    nothing_fetched, all_dropped, unknown = [], [], []
    for taxon in built:
        tid = f"NCBITaxon:{taxon}"
        n = per_virus.get(tid, 0)
        flag = ""
        if not n:
            rows = (fetched or {}).get(taxon, Counter()).get("rows", 0)
            if fetched is None:
                unknown.append(taxon)
                flag = "   <- no edge"
            elif not rows:
                nothing_fetched.append(taxon)
                flag = "   <- no edge, and no row was fetched for it"
            else:
                all_dropped.append((taxon, rows))
                flag = f"   <- no edge from {rows:,} fetched row(s)"
        out.append(f"  {label(tid):<14} {n:>7,}{flag}")

    zero = len(nothing_fetched) + len(all_dropped) + len(unknown)
    if zero:
        out.append("")
        out.append(f"  {zero} of {len(built)} built viruses have NO "
                   f"virus-host PPI edge from VirHostNet or IntAct.")
    # THE CAUSE IS COMPUTED, NOT ASSERTED. This said "the vhn_*.tab27 files
    # are one coronavirus each" -- true when written, false the moment the
    # flaviviruses were fetched, and still printed under five viruses whose
    # files were sitting in the directory. A sentence naming a cause has to
    # be derived from the same run that prints it or it becomes a confident
    # wrong answer the moment the thing it describes is fixed.
    if nothing_fetched:
        out.append(f"    nothing fetched ({len(nothing_fetched)}): "
                   + ", ".join(label(f"NCBITaxon:{t}") for t in nothing_fetched))
        out.append("      A taxon with no file cannot produce an edge "
                   "however well the source covers it. This is a FETCH gap "
                   "and scripts/fetch_vhppi.py is what closes it.")
    if all_dropped:
        out.append(f"    fetched but nothing survived ({len(all_dropped)}):")
        for taxon, rows in all_dropped:
            c = (fetched or {}).get(taxon, Counter())
            gates = [(g, c[g]) for g in JOIN_GATES if c[g]]
            why = ", ".join(f"{g} {n:,}" for g, n in gates) or "no gate recorded"
            out.append(f"      {label(f'NCBITaxon:{taxon}'):<14} "
                       f"{rows:,} row(s): {why}")
        out.append("      This is a JOIN gap, not a fetch gap: the rows are "
                   "on disk and the ingest drops them. Fetching again "
                   "changes nothing.")
    if unknown:
        out.append(f"    {len(unknown)} with no per-file tally available, so "
                   f"this run cannot say whether they were fetched.")
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
    per_file: dict = {}
    for src, fs in groups:
        if not fs:
            continue
        part = ingest_mitab(em, fs, index, taxon_map,
                            source=src, default_date=args.date,
                            per_file=per_file)
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
    fetched = rows_per_taxon(per_file, yaml.safe_load(
        args.config.read_text())["viruses"], taxon_map)
    for line in coverage_lines(per_virus, taxon_map, index.nodes, fetched):
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
