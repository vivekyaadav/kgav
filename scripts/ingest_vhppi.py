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
from kgav.vhppi import file_source, file_taxon, ingest_mitab, slug
from kgav.virus_register import declared_taxa, ppi_accession_map

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
        taxon = file_taxon(path, by_stem, taxon_map)
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
    appear at all. Before 2026-10-27 the flavivirus files had not been
    fetched, Zika sat in the taxon map producing nothing, and the build
    said so nowhere. The absence was discoverable only by noticing which
    labels were missing from a list that never claimed to be complete.

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
    register = yaml.safe_load(args.config.read_text())
    taxon_map = declared_taxa(register)
    print(f"  {len(taxon_map)} viral taxids mapped to "
          f"{len(set(taxon_map.values()))} species")

    # ATTRIBUTION, NOT IDENTITY. Some flavivirus rows are curated against an
    # unreviewed strain-specific polyprotein; the register says which
    # reference polyprotein each is attributed to and on what basis, and the
    # edge records the accession actually measured.
    acc_map = ppi_accession_map(register)
    if acc_map:
        print(f"  {len(acc_map)} measured accession(s) attributed to "
              f"{len(set(acc_map.values()))} reference polyprotein(s)")

    # TWO SOURCES, ONE LAYER. This globbed vhn_*.tab27 only, so
    # intact.zip -- 1.35 GB, in the same directory -- went unread until
    # 2026-10-23: the pattern did not match it. The seven vhn files were
    # one coronavirus each and the flaviviruses were missing entirely,
    # which is why Zika produced zero edges while sitting in the taxon
    # map. Both gaps closed 2026-10-27; the per-source lines this prints
    # are the current answer, and no count is repeated here to go stale.
    #
    # They share a predicate and a layer, and assemble_release already
    # merges one fact from two sources into one edge carrying both
    # ("infores:string|infores:virhostnet"). What must NOT be shared is the
    # label: an IntAct row recorded as infores:virhostnet is a false
    # provenance claim, and the tier ablation reads these labels.
    files = sorted(args.vhppi_dir.glob("*.tab27"))
    if not files:
        raise SystemExit(
            f"no *.tab27 files in {args.vhppi_dir}.\n"
            f"scripts/fetch_vhppi.py fetches from a PSICQUIC service, and "
            f"if intact.zip is present scripts/filter_intact.py cuts the "
            f"viral rows out of it.")

    # EVERY FILE, GROUPED BY WHAT ITS SIDECAR SAYS. This globbed vhn_* and
    # intact_* and labelled by the pattern, which is correct exactly while
    # no third service is ever fetched -- and fetch_vhppi takes --base, so
    # a MINT file would have landed as vhn_*.tab27 and been recorded as
    # VirHostNet. A glob cannot know where a file came from; the sidecar
    # the fetcher wrote can, and a file predating sidecars is covered by
    # the declared legacy prefixes.
    grouped: dict[str, list] = {}
    unlabelled = []
    for f in files:
        src = file_source(f)
        if not src:
            unlabelled.append(f)
            continue
        grouped.setdefault(src, []).append(f)
    if unlabelled:
        # REFUSED, NOT GUESSED. Picking a source here is the exact false
        # provenance claim this grouping exists to prevent, and the tier
        # ablation reads these labels as fact.
        raise SystemExit(
            "no source recorded for: "
            + ", ".join(f.name for f in unlabelled)
            + "\nEach raw file needs a <file>.source.json naming the "
              "service, or a filename prefix in vhppi.LEGACY_PREFIX_SOURCE."
              "\nRefusing rather than guessing: a row labelled with the "
              "wrong service is read as fact by the tier ablation.")
    groups = sorted(grouped.items())
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
                            per_file=per_file, accession_map=acc_map)
        print(f"\n  {src}: {part['rows']:,} rows -> {part['edges']:,} edges")
        s += part

    print(f"  {s['rows']:,} rows -> {s['edges']:,} edges")
    print(f"  {s['direct']:,} direct interaction, {s['co_complex']:,} co-complex/"
          f"association ({s['binary_method']:,} from a binary assay)")
    if s["attributed"]:
        # SAID OUT LOUD, EVERY RUN. An edge that exists only because the
        # register attributed one accession to another is a curated claim,
        # and a curated claim that shows up as an ordinary edge is how a
        # decision stops being reviewable.
        print(f"\n  {s['attributed']:,} edge(s) carry measured_accession: "
              f"attributed to a reference polyprotein by the register's "
              f"ppi_accession_map, at polyprotein resolution")
    by_curator = {k[4:]: v for k, v in s.items() if k.startswith("src:")}
    if by_curator:
        # SAID EVERY RUN. The file a row arrived in is not the authority
        # that curated it: IntAct's release carries IMEx partner
        # records, and 1,171 of its 12,496 viral rows are MINT, DIP,
        # HPIDb or MatrixDB. Those were filed as infores:intact until
        # 2026-10-28, and the tier ablation reads these labels.
        print("\n  edges by CURATOR, read from MITAB column 12 rather than "
              "from the file:")
        for k, v in sorted(by_curator.items(), key=lambda kv: -kv[1]):
            print(f"    {k:<28} {v:>8,}")
        if s["source_db_absent"]:
            print(f"    {'(no source column)':<28} "
                  f"{s['source_db_absent']:>8,}  fell back to the file's "
                  f"service")
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
    fetched = rows_per_taxon(per_file, register["viruses"], taxon_map)
    for line in coverage_lines(per_virus, taxon_map, index.nodes, fetched):
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
