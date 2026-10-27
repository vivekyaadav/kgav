#!/usr/bin/env python
"""Fetch virus-host interactions per virus from a PSICQUIC service.

    python scripts/fetch_vhppi.py --base <restUrl> --dry-run   # counts only
    python scripts/fetch_vhppi.py --base <restUrl>             # write files

WHY THIS EXISTS (HISTORY). data/raw/vhppi held seven vhn_*.tab27 files,
one coronavirus each, and the flaviviruses were missing entirely --
nine of sixteen built viruses with no virus-host PPI at all, and 9 edges
between them. The files were placed by hand from a host that has since
moved, so there was nothing to re-run and nothing recording where they
came from. This script closed that: the flaviviruses were fetched
2026-10-27 and every file it writes carries a sidecar. Ask the directory
what is there now rather than this paragraph.

THE BASE URL IS REQUIRED AND NOT DEFAULTED. VirHostNet moved from
virhostnet.prabi.fr to virhostnet.univ-lyon1.fr and distributes through
PSICQUIC; a URL written into this file from memory would be the assumption
this project refuses, and a stale one would fail in a way that looks like
an empty database. Get the current endpoint from the PSICQUIC registry,
and note `protocol=rest`:

    curl -sS "https://www.ebi.ac.uk/Tools/webservices/psicquic/registry/\\
registry?action=ACTIVE&format=txt&protocol=rest" | grep -i virhostnet

WITHOUT protocol=rest THE REGISTRY RETURNS THE SOAP ENDPOINT, and this
docstring said so for exactly as long as it took someone to follow it. The
two differ by a path segment and the SOAP one answers a REST query with
HTTP 500 and a soap:Fault reading "No binding operation info while
invoking unknown method" -- which reads as a broken service rather than a
wrong URL, and nearly sent this layer to a different data source.

The registry's REST entry ends at .../webservices/current/search/ and
--base takes everything before /webservices, so VirHostNet's is
http://virhostnet.prabi.fr:9090/psicquic -- the host that serves the data,
even though the web interface has moved to virhostnet.univ-lyon1.fr.

It is recorded in every sidecar this writes, so a release can always say
which service answered.

ISOLATE TAXA ARE QUERIED TOO. The register exists because sources curate
under a strain rather than the species -- seven of nine flaviviruses are
curated under a strain at UniProt, and DENV-2's chains live under 31634.
Querying the species alone would return nothing for exactly the viruses
this script was written to reach, and the emptiness would read as an
absent interactome rather than a wrong query.

A VIRUS WITH NO ROWS IS REPORTED, NOT OMITTED. That is the defect this
layer has produced twice: a per-virus table built from what was found
cannot distinguish a virus with no data from a virus nobody asked about.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import yaml

from kgav.vhppi import LEGACY_PREFIX_SOURCE, slug

TIMEOUT = 180
# PSICQUIC's own REST shape, not a VirHostNet invention:
#   <base>/webservices/current/search/query/<MIQL>?format=<fmt>
# Spec: psicquic.github.io/PsicquicSpec_1_4_Rest.html
QUERY_PATH = "webservices/current/search/query"


def prefix_for(source: str) -> str:
    """The filename stem prefix for a service.

    The two services already on disk keep the prefixes their files
    already carry -- renaming them would orphan the sidecars and the
    provenance those record. Anything else derives its prefix from the
    infores id, so a new service cannot collide with them by accident.
    """
    for pfx, infores in LEGACY_PREFIX_SOURCE.items():
        if infores == source:
            return pfx
    return re.sub(r"[^a-z0-9]+", "", source.removeprefix("infores:").lower())


def taxids_for(entry: dict) -> list[str]:
    """Species taxon plus every isolate the register declares.

    chain_source.taxon is folded in for the same reason declared_taxa folds
    it in: a curation that lives under the chain taxon is still this virus.
    """
    out = [str(entry["taxon"])]
    out += [str(t) for t in (entry.get("isolate_taxa") or [])]
    cs = (entry.get("chain_source") or {}).get("taxon")
    if cs:
        out.append(str(cs))
    seen, uniq = set(), []
    for t in out:
        if t not in seen:
            seen.add(t)
            uniq.append(t)
    return uniq


def miql(taxids: list[str]) -> str:
    """Either side of the interaction, over every taxid for this virus."""
    parts = [f"taxidA:{t} OR taxidB:{t}" for t in taxids]
    return " OR ".join(f"({p})" for p in parts)


def url_for(base: str, taxids: list[str], fmt: str) -> str:
    q = urllib.parse.quote(miql(taxids), safe="")
    return f"{base.rstrip('/')}/{QUERY_PATH}/{q}?format={fmt}"


def fetch(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "kgav-fetch"})
    with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
        return r.read()


def sidecar_path(target: Path) -> Path:
    return target.with_suffix(target.suffix + ".source.json")


def write_with_provenance(target: Path, raw: bytes, url: str, entry: dict,
                          taxids: list[str], n_rows: int,
                          source: str = "") -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(raw)
    sidecar_path(target).write_text(json.dumps({
        "url": url,
        # WHICH SERVICE ANSWERED. ingest_vhppi reads this to label the
        # edge; without it the label comes from the filename, and a file
        # from one service in another's naming scheme is a false
        # provenance claim the tier ablation then reads as fact.
        "source": source,
        "retrieved": datetime.now(UTC).strftime("%Y-%m-%d"),
        "sha256": hashlib.sha256(raw).hexdigest(),
        "bytes": len(raw),
        "virus": entry["name"],
        "taxon": str(entry["taxon"]),
        "taxids_queried": taxids,
        "n_rows": n_rows,
        "fetched_by": "scripts/fetch_vhppi.py",
        "consumed_by": "scripts/ingest_vhppi.py",
    }, indent=2, sort_keys=True) + "\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--base", required=True,
                    help="PSICQUIC restUrl, from the registry. Required: a "
                         "URL defaulted from memory is the assumption this "
                         "project refuses.")
    ap.add_argument("--source", required=True,
                    help="infores id of the service at --base, e.g. "
                         "infores:virhostnet. Required for the same reason "
                         "--base is: the file is labelled with this, and a "
                         "service guessed from a filename is a false "
                         "provenance claim.")
    ap.add_argument("--config", type=Path, default=ROOT / "config/viruses.yaml")
    ap.add_argument("--out", type=Path, default=ROOT / "data/raw/vhppi")
    ap.add_argument("--family", default=None,
                    help="only viruses whose register family matches, e.g. "
                         "Flaviviridae")
    ap.add_argument("--format", default="tab25",
                    help="PSICQUIC format. The site serves PSI-MITAB 2.5 and "
                         "vhppi.parse_row reads columns 0-14, which 2.5 has.")
    ap.add_argument("--dry-run", action="store_true",
                    help="ask each virus for a COUNT and write nothing")
    ap.add_argument("--overwrite", action="store_true",
                    help="replace a file that already exists")
    args = ap.parse_args(argv)

    viruses = yaml.safe_load(args.config.read_text())["viruses"]
    if args.family:
        viruses = [v for v in viruses
                   if (v.get("family") or "") == args.family]
    if not viruses:
        print(f"no virus in the register matches --family {args.family!r}")
        return 2

    print(f"{len(viruses)} virus(es) from {args.config.name} against "
          f"{args.base}")
    print(f"  labelling every file {args.source}")
    empty, failed, wrote = [], [], []
    for entry in viruses:
        taxids = taxids_for(entry)
        name = entry["name"]
        count_url = url_for(args.base, taxids, "count")
        try:
            n = int(fetch(count_url).decode().strip() or 0)
        except (urllib.error.URLError, ValueError, OSError) as exc:
            print(f"  {name:<14} QUERY FAILED: {exc}")
            failed.append(name)
            continue
        taxa_note = f" ({len(taxids)} taxids)" if len(taxids) > 1 else ""
        print(f"  {name:<14} {n:>7,} interaction(s){taxa_note}")
        if not n:
            empty.append(name)
            continue
        if args.dry_run:
            continue
        target = args.out / f"{prefix_for(args.source)}_{slug(name)}.tab27"
        if target.exists() and not args.overwrite:
            print(f"    exists, left alone: {target.name} (--overwrite to "
                  f"replace)")
            continue
        data_url = url_for(args.base, taxids, args.format)
        try:
            raw = fetch(data_url)
        except (urllib.error.URLError, OSError) as exc:
            print(f"    DOWNLOAD FAILED: {exc}")
            failed.append(name)
            continue
        rows = sum(1 for ln in raw.decode("utf-8", "replace").splitlines()
                   if ln.strip())
        if rows != n:
            # Not fatal: PSICQUIC counts interactions and a MITAB line is a
            # binary pair, so the two can legitimately differ. Reported
            # because a silent mismatch is how a truncated download looks
            # exactly like a small one.
            print(f"    note: {rows:,} row(s) written against a count of "
                  f"{n:,} -- check before trusting either")
        write_with_provenance(target, raw, data_url, entry, taxids, rows,
                              source=args.source)
        print(f"    wrote {target.name}  {len(raw):,} bytes, {rows:,} rows")
        wrote.append(name)

    print()
    if wrote:
        print(f"  {len(wrote)} file(s) written to {args.out}")
    if empty:
        print(f"  {len(empty)} virus(es) returned NO interaction: "
              + ", ".join(empty))
        print("  That is this service reporting nothing for them, not a "
              "failed fetch. It is a source limit and belongs in the "
              "methods, not a retry.")
    if failed:
        print(f"  {len(failed)} virus(es) FAILED to query: "
              + ", ".join(failed))
        print("  A failure is not an empty result; nothing was written for "
              "these and the gap is still open.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
