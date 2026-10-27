#!/usr/bin/env python3
"""Cut the viral rows out of IntAct's full release, without unpacking it.

    python scripts/filter_intact.py            # writes intact_viral.tab27
    python scripts/filter_intact.py --inspect  # counts per taxon, writes nothing

WHY THIS EXISTS. data/raw/vhppi holds seven vhn_<taxid>.tab27 files, one per
coronavirus, fetched 2026-09-11 -- and intact.zip, 1.35 GB, dated 2026-01-14,
which nothing has ever read. ingest_vhppi globs `vhn_*.tab27`, so the larger
and broader source sat beside the narrower one, unused and unmentioned, while
the flavivirus PPI gap was recorded as a missing download.

That gap is real for VirHostNet: no flavivirus file was ever fetched, which
is why Zika produced zero edges despite being in the taxon map since
2026-10-12. IntAct covers all species, so the rows may be here already.

STREAMED, NOT EXTRACTED. The archive expands to several GB of which a tiny
fraction is viral. zipfile.open gives a file object over the compressed
member, so this reads it in one pass and never writes the whole thing to
disk.

THE FILTER IS THE REGISTER'S OWN TAXON MAP. A row is kept when either
interactor's taxid is a taxon config/viruses.yaml declares -- species or
isolate, via virus_register.declared_taxa, the same map ingest_vhppi
resolves against. Nothing is matched by name, and a taxid the register does
not declare is not kept on the grounds that it looks viral.

WHAT THIS DOES NOT DO. It does not decide whether a row is host-viral, drop
self-interactions, or read detection methods. src/kgav/vhppi.py does all of
that and is not duplicated here; this only makes the file small enough to
hand it.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import zipfile
from collections import Counter
from contextlib import ExitStack
from datetime import UTC, datetime
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.virus_register import declared_taxa

TAXID_RE = re.compile(rb"taxid:(-?\d+)")
# MITAB 2.7 column positions, the same ones src/kgav/vhppi.py reads.
COL_TAX_A, COL_TAX_B = 9, 10
MIN_COLUMNS = 15


#: A member whose name marks it as IntAct's NEGATIVE set. Those rows record
#: pairs an experiment showed do NOT interact, in the same MITAB columns as
#: the positives, so nothing downstream can tell them apart once they are in
#: a .tab27: ingest_vhppi would emit each one as a PHYSICALLY_INTERACTS_WITH
#: edge asserting the opposite of what was measured.
NEGATIVE_MARK = "negative"


def is_negative_member(name: str) -> bool:
    return NEGATIVE_MARK in Path(name).name.lower()


def check_member(name: str, *, allow_negative: bool = False) -> str:
    """The member to read, or exit rather than invert the evidence.

    IntAct's release holds intact.txt AND intact_negative.txt. The inferring
    branch below refuses to choose between them, which is what kept the
    negative set out on the --inspect run that found them -- but --member
    takes whatever it is handed, and the write step is where that matters.

    Refused by NAME rather than by content because the MITAB columns this
    script reads are identical in both files; by the time a row is here
    there is nothing left to distinguish it by. --allow-negative exists so
    the refusal can be overridden deliberately by someone who means to
    study the negative set, never by a default.
    """
    if is_negative_member(name) and not allow_negative:
        raise SystemExit(
            f"{name} is IntAct's NEGATIVE interaction set: rows recording "
            f"that a pair does NOT interact.\n"
            f"ingest_vhppi cannot tell them from positives -- it would emit "
            f"each as a PHYSICALLY_INTERACTS_WITH edge\n"
            f"asserting the opposite of the experiment. Use the positive "
            f"member (intact.txt), or --allow-negative if\n"
            f"you mean to filter the negative set for its own sake.")
    return name


def member_name(zf: zipfile.ZipFile) -> str:
    """The one data member. IntAct ships intact.txt; be explicit if not."""
    names = [n for n in zf.namelist() if n.lower().endswith((".txt", ".tsv"))]
    if len(names) != 1:
        raise SystemExit(
            f"expected exactly one .txt/.tsv member, found {names or 'none'}. "
            f"Name it with --member rather than letting this guess.")
    return names[0]


def taxids(line: bytes) -> set[bytes]:
    f = line.split(b"\t")
    if len(f) < MIN_COLUMNS:
        return set()
    out = set()
    for col in (COL_TAX_A, COL_TAX_B):
        out |= set(TAXID_RE.findall(f[col]))
    return out


def scan(path: Path, wanted: set[bytes], member: str | None,
         out: Path | None, allow_negative: bool = False) -> Counter:
    """One pass. Writes matching lines to `out`, or counts only if None.

    WRITES TO A TEMPORARY AND RENAMES. Opening `out` directly left an empty
    file behind when nothing matched, under a run that printed "Nothing
    written" -- and an empty .tab27 is worse than none, because ingest_vhppi
    globs it, parses zero edges, and reports a source that was consulted and
    had nothing to say. A half-written file from an interrupted run has the
    same shape, and the rename makes both impossible.
    """
    stats: Counter = Counter()
    tmp = out.with_suffix(out.suffix + ".partial") if out is not None else None
    with ExitStack() as stack:
        sink = stack.enter_context(open(tmp, "wb")) if tmp is not None else None
        zf = stack.enter_context(zipfile.ZipFile(path))
        chosen = check_member(member or member_name(zf),
                              allow_negative=allow_negative)
        fh = stack.enter_context(zf.open(chosen))
        for line in fh:                      # streamed; never read whole
            stats["rows"] += 1
            if line.startswith(b"#"):
                stats["header"] += 1
                continue
            hit = taxids(line) & wanted
            if not hit:
                continue
            stats["kept"] += 1
            for t in sorted(hit):
                stats[f"taxon:{t.decode()}"] += 1
            if sink is not None:
                sink.write(line)
    if tmp is not None:
        if stats["kept"]:
            tmp.replace(out)
        else:
            tmp.unlink(missing_ok=True)
    return stats


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--zip", type=Path,
                    default=ROOT / "data/raw/vhppi/intact.zip")
    ap.add_argument("--member", default=None,
                    help="member to read; inferred when there is only one")
    ap.add_argument("--allow-negative", action="store_true",
                    help="read a member whose name marks it as IntAct's "
                         "NEGATIVE set. Refused by default: those rows say "
                         "a pair does NOT interact, and nothing downstream "
                         "can tell them from positives.")
    ap.add_argument("--config", type=Path, default=ROOT / "config/viruses.yaml")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "data/raw/vhppi/intact_viral.tab27")
    ap.add_argument("--inspect", action="store_true",
                    help="report what would be kept; write nothing")
    args = ap.parse_args(argv)

    if not args.zip.exists():
        print(f"no archive at {args.zip}")
        print("IntAct's full release is at "
              "https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/")
        return 1

    cfg = yaml.safe_load(args.config.read_text())
    taxon_map = declared_taxa(cfg)
    wanted = {t.encode() for t in taxon_map}
    species = {v["taxon"]: v["name"] for v in cfg["viruses"]}
    print(f"{args.zip.name} ({args.zip.stat().st_size / 1e9:.2f} GB)")

    # RESOLVED HERE, NOT INSIDE THE SCAN. The sidecar recorded
    # `args.member or "inferred"`, so a run that inferred the member left a
    # provenance record that does not say which file it read -- and this
    # archive holds two, one of which asserts the opposite of the other.
    # Resolving first also puts the negative-set refusal before the pass
    # rather than after the zip is open.
    with zipfile.ZipFile(args.zip) as zf:
        member = check_member(args.member or member_name(zf),
                              allow_negative=args.allow_negative)
    note = " (NEGATIVE SET, by explicit request)" if \
        is_negative_member(member) else ""
    print(f"reading member {member}{note}")
    print(f"filtering to {len(wanted)} declared viral taxids "
          f"({len(set(taxon_map.values()))} species)\n")

    stats = scan(args.zip, wanted, member,
                 None if args.inspect else args.out,
                 allow_negative=args.allow_negative)

    print(f"{stats['rows']:,} rows read, {stats['header']} header, "
          f"{stats['kept']:,} kept\n")
    rows = sorted(((k.split(":", 1)[1], v) for k, v in stats.items()
                   if k.startswith("taxon:")), key=lambda r: -r[1])
    if rows:
        print("  rows per declared taxid (a row naming two is counted twice)")
        for tax, n in rows:
            sp = taxon_map.get(tax, tax)
            label = species.get(sp, "")
            via = "" if sp == tax else f"  -> {sp}"
            print(f"    {tax:<10} {n:>8,}  {label}{via}")
    else:
        print("  NO ROW names any declared viral taxon.")
        print("  Either this archive is not what it is taken for, or IntAct "
              "files these\n  interactions under taxids the register does "
              "not declare. Check one\n  known pair by hand before "
              "concluding the data is absent.")

    if args.inspect:
        print("\n--inspect: nothing written.")
        return 0
    if not stats["kept"]:
        print("\nNothing written: an empty file would parse to zero edges and "
              "read as\na source that was consulted and had nothing.")
        return 1

    digest = hashlib.sha256(args.out.read_bytes()).hexdigest()[:16]
    side = args.out.with_suffix(args.out.suffix + ".source.json")
    side.write_text(json.dumps({
        "derived_from": str(args.zip.name),
        "archive_bytes": args.zip.stat().st_size,
        "member": member,
        "negative_set": is_negative_member(member),
        "rows_kept": stats["kept"],
        "rows_read": stats["rows"],
        "sha256_16": digest,
        "filtered_on": sorted(t.decode() for t in wanted),
        "written": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }, indent=2) + "\n")
    print(f"\nwrote {args.out} ({args.out.stat().st_size:,} bytes) "
          f"+ {side.name}")
    print("ingest_vhppi reads it as infores:intact, beside the "
          "vhn_*.tab27 files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
