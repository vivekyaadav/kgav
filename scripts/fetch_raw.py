#!/usr/bin/env python
"""Fetch the UniProt files ingest_spine reads, and record where they came from.

WHY THIS EXISTS

Nothing in this repository said where data/raw/proteomes came from. There was
no fetcher, no documented URL, no retrieval date and no checksum -- the whole
release rested on files placed by hand, and `data/raw/*` is gitignored, so a
fresh clone could not rebuild v0.1 at all. The first visible consequence was
the DENV-2/Zika promotion: the register declared two viruses, the proteomes
were simply absent, ingest_spine recorded `proteome_missing` and returned
early, and the build carried on to invalidate eleven results files for two
orphan nodes. A missing file was indistinguishable from a file nobody had
fetched yet, because there was no step that was supposed to fetch it.

WHAT IT VERIFIES, AND WHY EACH CHECK IS HERE

Downloading is the easy half. Each check below exists because the thing it
checks has already gone wrong once in this project:

  shape          load_proteome does `json.loads(raw).get("results", [])` and
                 returns [] for anything else -- including an error page --
                 so a 200 response carrying a JSON error object would ingest
                 as "a virus with no proteins" and look exactly like a
                 missing file.

  chain_source   config/viruses.yaml pins, per virus, the accession the mature
                 peptides really come from, that accession's own taxon, and
                 how many Chain features it carries, with a `verified:` date.
                 That date was a human claim nothing re-checked. This compares
                 the claim against the downloaded file and fails on a
                 mismatch, so the register's curation is testable rather than
                 asserted. (Checked 2026-10-11: Q32ZE1 and P29990 both still
                 carry 13 chains, matching `verified: 2026-09-17`.)

  the MERS trap  Swiss-Prot sometimes curates a reference ISOLATE rather than
                 the species, leaving the species proteome TrEMBL-only and
                 therefore free of any Chain annotation. A proteome that
                 downloads cleanly and contains zero chains builds proteins
                 with no mature peptides, and every ChEMBL row then lands on
                 the polyprotein. DENV-2 is the live case: P29990 carries 13
                 chains and has NO proteome cross-reference at all, so no
                 proteome download produces it and extra_accessions is the
                 only route. This refuses to write a proteome whose
                 chain_source accession is missing from it, and names
                 extra_accessions as the fix.

Provenance goes in a sidecar `<name>.source.json`: the exact URL, the UTC
date, a sha256, and the counts the checks were made against. `--check` then
re-validates what is already on disk WITHOUT network, which is how a
hand-placed file with no provenance gets found.

USAGE

  python scripts/fetch_raw.py --virus Zika            # proteome (+ EXTRA)
  python scripts/fetch_raw.py --virus DENV-2
  python scripts/fetch_raw.py --all-built             # every virus in viruses:
  python scripts/fetch_raw.py --virus Zika --dry-run  # print URLs, fetch none
  python scripts/fetch_raw.py --check                 # audit files on disk

This script writes to data/raw/proteomes and nowhere else. It never touches
data/releases, so it cannot invalidate a results file.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "viruses.yaml"
PROTEOMES = ROOT / "data" / "raw" / "proteomes"

STREAM = "https://rest.uniprot.org/uniprotkb/stream"
TIMEOUT = 120

# The `stream` endpoint rather than `search`: search paginates and silently
# truncates at its default size, which for a proteome means "some of the
# proteins". stream returns the whole result set in one body. Verified
# 2026-10-11 that it returns {"results": [...]} -- the shape load_proteome
# expects -- for (proteome:UP000054557), with 2 entries: Q32ZE1 and P0DXN9.
def _url(query: str) -> str:
    return f"{STREAM}?" + urllib.parse.urlencode(
        {"query": query, "format": "json", "compressed": "true"})


def proteome_url(proteome_id: str) -> str:
    return _url(f"(proteome:{proteome_id})")


def accession_url(accessions: list[str]) -> str:
    return _url(" OR ".join(f"accession:{a}" for a in accessions))


def fetch_bytes(url: str) -> bytes:
    req = urllib.request.Request(url, headers={
        "User-Agent": "kgav/0.1 (+https://github.com/vivekyaadav/kgav)"})
    with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
        return r.read()


# --------------------------------------------------------------------------
# validation
# --------------------------------------------------------------------------
def _decompress(raw: bytes) -> bytes:
    return gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw


def chain_count(entry: dict) -> int:
    return len([f for f in entry.get("features", []) if f.get("type") == "Chain"])


def inspect(raw: bytes) -> dict:
    """Parse and describe a downloaded body, or raise ValueError saying why not.

    Raises rather than returning a flag: a body that cannot be described is
    one ingest_spine would read as an empty proteome, and that failure mode --
    indistinguishable from a missing file -- is the one this script exists to
    remove.
    """
    try:
        doc = json.loads(_decompress(raw))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as e:
        raise ValueError(f"not gzipped JSON: {e}") from e
    if not isinstance(doc, dict) or "results" not in doc:
        keys = sorted(doc)[:6] if isinstance(doc, dict) else type(doc).__name__
        raise ValueError(
            f"no top-level 'results' key (got {keys}). load_proteome would "
            f"read this as an empty proteome")
    entries = doc["results"]
    if not entries:
        raise ValueError("'results' is empty -- the query matched nothing")
    return {
        "n_entries": len(entries),
        "accessions": [e.get("primaryAccession") for e in entries],
        "chains": {e.get("primaryAccession"): chain_count(e) for e in entries},
        "total_chains": sum(chain_count(e) for e in entries),
    }


def check_against_register(report: dict, virus: dict, *, is_extra: bool) -> list[str]:
    """Compare a downloaded file against what config/viruses.yaml claims.

    Returns human-readable problems; an empty list means the register's claim
    held. The caller refuses to write when this is non-empty.
    """
    cs = virus.get("chain_source") or {}
    acc, claimed = cs.get("accession"), cs.get("chains")
    extra = virus.get("extra_accessions") or []
    problems: list[str] = []

    def _count(where: str) -> list[str]:
        """Compare the pinned chain count against the file that supplies it."""
        if not acc or claimed is None or report["chains"].get(acc) == claimed:
            return []
        return [f"{acc} carries {report['chains'].get(acc)} Chain features in "
                f"the {where}, register claims {claimed} "
                f"(verified {cs.get('verified')}) -- UniProt changed, or the "
                f"claim was wrong"]

    # WHICH FILE HOLDS THE CHAINS IS THE WHOLE QUESTION.
    #
    # ingest_spine does `entries = curated` (scripts/ingest_spine.py:221):
    # curated accessions REPLACE the proteome rather than augmenting it. So
    # for a virus declaring extra_accessions, the proteome's CONTENT is never
    # read -- only its existence matters, because the pfile.exists() guard
    # returns early before the curated branch is reached.
    #
    # The first run of --check got this wrong and reported MERS-CoV as
    # DISAGREEING: chain_source pins K9N7C7, the proteome holds eleven
    # unrelated K0B*/K4L* accessions, and that is exactly how MERS is designed
    # -- K9N7C7 is filed under isolate taxon 1263720 and arrives via
    # EXTRA_1335626.json.gz. Looking for the chains in the proteome asked the
    # wrong file. A check that cries wolf on a correctly curated virus gets
    # switched off, which would have cost more than it saved.
    replaced = bool(extra)

    if is_extra:
        missing = [a for a in extra if a not in report["accessions"]]
        if missing:
            problems.append(
                f"extra_accessions {missing} not in the response -- the "
                f"accession query returned {report['accessions']}")
        if acc in extra:
            problems += _count("curated file")
        if report["total_chains"] == 0:
            problems.append(
                "the curated file has zero Chain features, so it replaces the "
                "proteome with entries that carry no mature peptides and every "
                "ChEMBL row lands on the polyprotein")
        return problems

    if replaced:
        # Gate token only: it must exist, and nothing about its content is
        # consumed. Checking it would be checking a file the build discards.
        return problems

    if acc and acc not in report["accessions"]:
        problems.append(
            f"chain_source pins {acc} but the proteome contains "
            f"{report['accessions']} -- add extra_accessions: [{acc}] to "
            f"{virus['name']} and re-run, which is the MERS mechanism")
    else:
        problems += _count("proteome")

    if report["total_chains"] == 0:
        problems.append(
            "the proteome has zero Chain features, so it builds proteins with "
            "no mature peptides and every ChEMBL row lands on the "
            "polyprotein (the MERS trap). A curated accession is needed")
    return problems


# --------------------------------------------------------------------------
# writing
# --------------------------------------------------------------------------
def sidecar_path(target: Path) -> Path:
    return target.with_suffix(target.suffix + ".source.json")


def write_with_provenance(target: Path, raw: bytes, url: str, report: dict,
                          virus_name: str) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(raw)
    sidecar_path(target).write_text(json.dumps({
        "url": url,
        "retrieved": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
        "sha256": hashlib.sha256(raw).hexdigest(),
        "bytes": len(raw),
        "virus": virus_name,
        "n_entries": report["n_entries"],
        "accessions": report["accessions"],
        "chains": report["chains"],
        "fetched_by": "scripts/fetch_raw.py",
    }, indent=2, sort_keys=True) + "\n")


def load_register(path: Path = CONFIG) -> dict:
    return yaml.safe_load(path.read_text())


def find_virus(cfg: dict, name: str) -> tuple[dict, str]:
    for section in ("viruses", "surveyed"):
        for v in cfg.get(section) or []:
            if v["name"].lower() == name.lower():
                return v, section
    names = [v["name"] for s in ("viruses", "surveyed") for v in (cfg.get(s) or [])]
    raise SystemExit(f"no virus named {name!r} in the register. Known: {names}")


def targets_for(virus: dict) -> list[tuple[Path, str, bool]]:
    """(target path, url, is_extra) for one virus, proteome first."""
    out = [(PROTEOMES / f"{virus['proteome_id']}.json.gz",
            proteome_url(virus["proteome_id"]), False)]
    if virus.get("extra_accessions"):
        out.append((PROTEOMES / f"EXTRA_{virus['taxon']}.json.gz",
                    accession_url(virus["extra_accessions"]), True))
    return out


def fetch_virus(virus: dict, *, force: bool, dry_run: bool,
                fetch=fetch_bytes) -> int:
    """Fetch one virus's files. Returns the number of problems found.

    A FAILED PROTEOME IS NOT ALWAYS A PROBLEM. When the register declares
    extra_accessions, those curated entries REPLACE the proteome in
    ingest_spine -- the proteome file is not read at all. UniProt retired
    the TrEMBL-only flavivirus species proteomes, so six entries now fetch
    their EXTRA file successfully and their proteome not at all, which is
    the intended outcome.

    This counted every failure, so that run printed "6 problem(s)" and
    "Nothing was written for a refused file" after writing six files
    correctly -- telling the reader the ingest was broken when it had done
    exactly what it should. It also used one word for two causes: REFUSED
    means the content disagreed with the register, FAILED means the query
    returned nothing, and only the first is a reason to distrust a source.
    """
    problems = 0
    outcomes: dict[bool, str] = {}
    substituted = bool(virus.get("extra_accessions"))
    for target, url, is_extra in targets_for(virus):
        label = f"{virus['name']} {'EXTRA' if is_extra else 'proteome'}"
        print(f"\n  {label}: {target.name}")
        print(f"    {url}")
        if dry_run:
            continue
        if target.exists() and not force:
            print("    already on disk; --force to refetch")
            continue
        try:
            raw = fetch(url)
            report = inspect(raw)
        except Exception as e:                       # noqa: BLE001 - reported
            print(f"    FAILED: {e}")
            outcomes[is_extra] = "failed"
            if not (substituted and not is_extra):
                problems += 1
            continue

        found = check_against_register(report, virus, is_extra=is_extra)
        print(f"    {report['n_entries']} entries, "
              f"{report['total_chains']} Chain features: {report['chains']}")
        if found:
            # Refuse to write. A file that disagrees with the register is
            # worse than no file: ingest_spine would consume it silently.
            for p in found:
                print(f"    REFUSED: {p}")
            problems += len(found)
            outcomes[is_extra] = "refused"
            continue
        write_with_provenance(target, raw, url, report, virus["name"])
        print(f"    wrote {len(raw):,} bytes + {sidecar_path(target).name}")
        outcomes[is_extra] = "written"

    # SAY SO, rather than leaving a FAILED line as the last word on a virus
    # that is in fact ready. The EXTRA file is what ingest_spine reads.
    if substituted and outcomes.get(False) == "failed":
        extra_state = outcomes.get(True, "already on disk")
        if extra_state in ("written", "already on disk"):
            print(f"    -> {virus['name']} is READY: the proteome is not "
                  f"read for this entry;\n       extra_accessions "
                  f"{virus['extra_accessions']} replace it, and that file "
                  f"is {extra_state}.")
        else:
            print(f"    -> {virus['name']} is NOT ready: the proteome "
                  f"failed and its replacement\n       "
                  f"({'EXTRA_' + virus['taxon']}.json.gz) is {extra_state}.")
    return problems


def _provenance_note(side: Path) -> str:
    """'' when really fetched, a note otherwise."""
    if not side.exists():
        return "  NO PROVENANCE"
    try:
        rec = json.loads(side.read_text())
    except (OSError, ValueError):
        return "  UNREADABLE PROVENANCE"
    if rec.get("provenance") == "reconstructed" or not rec.get("url"):
        return "  RECONSTRUCTED"
    return ""


def check_on_disk(cfg: dict) -> int:
    """Re-validate every file on disk against the register. No network."""
    problems = 0
    seen: set[Path] = set()
    print(f"auditing {PROTEOMES}\n")
    for section in ("viruses", "surveyed"):
        for v in cfg.get(section) or []:
            # THE SAME SUBSTITUTION fetch_virus ALREADY KNOWS ABOUT, and
            # this did not. When extra_accessions are declared and their
            # EXTRA file is present, ingest_spine never reads the proteome
            # -- UniProt retired seven of the nine flavivirus species
            # proteomes on or before 2026-10-23 and the register answers
            # that by naming the strain accession instead. Reported as
            # MISSING, those seven failed `make check-raw` for a build
            # that is correct, which is a guard crying wolf about its own
            # repository.
            xfile = PROTEOMES / f"EXTRA_{v['taxon']}.json.gz"
            substitutes = bool(v.get("extra_accessions")) and xfile.exists()
            for target, _url, is_extra in targets_for(v):
                seen.add(target)
                tag = f"{v['name']} ({section})"
                if not target.exists():
                    if substitutes and not is_extra:
                        print(f"  n/a       {target.name:<28} {tag}  "
                              f"not read: extra_accessions replace it")
                    elif section == "viruses":
                        print(f"  MISSING   {target.name:<28} {tag}")
                        problems += 1
                    continue
                try:
                    report = inspect(target.read_bytes())
                except ValueError as e:
                    print(f"  CORRUPT   {target.name:<28} {tag}: {e}")
                    problems += 1
                    continue
                found = check_against_register(report, v, is_extra=is_extra)
                # RECONSTRUCTED IS NOT FETCHED. record_raw_provenance writes
                # a sidecar for a hand-placed file carrying a checksum and
                # no URL; testing only that a sidecar exists would report
                # that as full provenance and lose the distinction the
                # reconstruction was careful to record.
                prov = _provenance_note(sidecar_path(target))
                if found:
                    print(f"  DISAGREES {target.name:<28} {tag}{prov}")
                    for p in found:
                        print(f"            {p}")
                    problems += len(found)
                else:
                    print(f"  ok        {target.name:<28} {tag}{prov}")
                if prov == "  NO PROVENANCE":
                    print("            placed by hand; no URL, date or "
                          "checksum. `make record-provenance` records a "
                          "checksum;\n            refetching with --force "
                          "records the URL too")
                elif prov == "  RECONSTRUCTED":
                    print("            checksum recorded after the fact; "
                          "the URL and date are not recoverable")
    stray = sorted(p for p in PROTEOMES.glob("*.json.gz") if p not in seen) \
        if PROTEOMES.exists() else []
    for p in stray:
        print(f"  STRAY     {p.name:<28} no register entry names this file")
    return problems


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--virus", action="append", default=[],
                    help="virus name from config/viruses.yaml (repeatable)")
    ap.add_argument("--all-built", action="store_true",
                    help="every virus in `viruses:`")
    ap.add_argument("--check", action="store_true",
                    help="validate files already on disk; no network")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the URLs and fetch nothing")
    ap.add_argument("--force", action="store_true",
                    help="refetch files that already exist")
    ap.add_argument("--config", type=Path, default=CONFIG)
    args = ap.parse_args(argv)

    cfg = load_register(args.config)
    if args.check:
        problems = check_on_disk(cfg)
        print(f"\n{problems} problem(s)")
        return 1 if problems else 0

    if args.all_built:
        viruses = list(cfg["viruses"])
    elif args.virus:
        viruses = [find_virus(cfg, n)[0] for n in args.virus]
    else:
        ap.error("one of --virus, --all-built or --check is required")

    print(f"{'would fetch' if args.dry_run else 'fetching'} "
          f"{len(viruses)} virus(es) into {PROTEOMES}")
    problems = sum(fetch_virus(v, force=args.force, dry_run=args.dry_run)
                   for v in viruses)
    if args.dry_run:
        print("\ndry run: nothing fetched, nothing written")
        return 0
    print(f"\n{problems} problem(s)")
    if problems:
        print("A REFUSED file is one whose contents disagree with the "
              "register; nothing was\nwritten for it. A FAILED fetch for an "
              "entry whose extra_accessions replace the\nproteome is not "
              "counted here, because that file is never read.\n")
        print("A proteome that "
              "disagrees with the register would be consumed silently.")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
