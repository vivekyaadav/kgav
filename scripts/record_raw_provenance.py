#!/usr/bin/env python3
"""Stamp a raw input that was placed by hand, and re-verify it later.

    python scripts/record_raw_provenance.py              # stamp data/raw/host
    python scripts/record_raw_provenance.py --check      # verify, write nothing
    python scripts/record_raw_provenance.py --dir data/raw/ncats

WHAT THIS IS FOR. `make data-timeline` reports the host layer as five
files, 328,846,356 bytes, `mtime only` -- the human proteome and STRING
v12.0, feeding 157,211 PARTICIPATES_IN and 141,019
PHYSICALLY_INTERACTS_WITH edges, over three-quarters of the graph, with no
recorded URL, release or checksum. fetch_raw.py writes a .source.json
sidecar for every proteome it fetches; nothing does it for these, because
they were downloaded by hand before that existed.

IT DOES NOT INVENT A DOWNLOAD. The URL and the retrieval date of a file
placed by hand are not recoverable, and writing a plausible value for
either would be a fabricated provenance record -- worse than none, because
a reader cannot tell it from a real one. The sidecar this writes states
`url: null`, `retrieved: null` and `provenance: "reconstructed"`, and
names the date the reconstruction happened.

WHAT IT DOES ESTABLISH, which is the part that matters going forward:

    sha256    so drift is detectable from now on. Today a file could be
              replaced with a different release and nothing in the repo
              would notice.
    bytes     a cheap second check, and what data_timeline already shows.
    version   read from the filename where the source puts it there --
              STRING names its release in 9606.protein.links.detailed.
              v12.0.txt.gz. That is a fact about the file, not a guess.

A FETCHER IS A DIFFERENT JOB. Making a fresh clone able to download these
with real provenance means writing the URLs and validating the responses,
the way fetch_raw.py does for proteomes. This script makes the files that
are already here auditable; it does not make them reproducible.

IT NEVER OVERWRITES AN EXISTING SIDECAR. A record written by fetch_raw
carries a real URL and a real retrieval date. Replacing it with a
reconstruction would lose exactly the information that distinguishes them.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Versions a source writes into its own filename. A match is a fact about
# the file; no match leaves the field null rather than guessing.
VERSION_PATTERNS = (
    re.compile(r"\b(v\d+\.\d+)\b"),        # STRING: ...links.detailed.v12.0...
    re.compile(r"\b(chembl_\d+)\b"),       # ChEMBL names its release
)

# What reads each file, so a sidecar says why the file is here. Keyed by
# filename; anything unlisted is still stamped, with consumer null.
CONSUMED_BY = {
    "human_proteome.json.gz": "scripts/ingest_host.py",
    "9606.protein.aliases.v12.0.txt.gz": "scripts/ingest_host.py",
    "9606.protein.links.detailed.v12.0.txt.gz": "scripts/ingest_host.py",
}

SKIP_SUFFIXES = (".source.json", ".partial")


def sidecar_path(target: Path) -> Path:
    return target.with_suffix(target.suffix + ".source.json")


def version_of(name: str) -> str | None:
    for pat in VERSION_PATTERNS:
        m = pat.search(name)
        if m:
            return m.group(1)
    return None


def sha256(path: Path, chunk: int = 1 << 20) -> str:
    """Streamed. These files run to hundreds of megabytes."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def candidates(directory: Path) -> list[Path]:
    return sorted(p for p in directory.iterdir()
                  if p.is_file() and not p.name.endswith(SKIP_SUFFIXES))


def record(path: Path, today: str) -> dict:
    return {
        "url": None,
        "retrieved": None,
        "provenance": "reconstructed",
        "recorded": today,
        "recorded_by": "scripts/record_raw_provenance.py",
        "note": ("Placed by hand before provenance was recorded. The URL "
                 "and retrieval date are not recoverable and are left null "
                 "rather than guessed; the checksum below is what makes "
                 "later drift detectable."),
        "sha256": sha256(path),
        "bytes": path.stat().st_size,
        "version_in_filename": version_of(path.name),
        "consumed_by": CONSUMED_BY.get(path.name),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dir", type=Path, default=ROOT / "data/raw/host")
    ap.add_argument("--check", action="store_true",
                    help="re-hash every file against its sidecar; write "
                         "nothing")
    args = ap.parse_args(argv)

    if not args.dir.is_dir():
        print(f"no directory at {args.dir}")
        return 2
    files = candidates(args.dir)
    if not files:
        print(f"no files in {args.dir}")
        return 2

    today = datetime.now(UTC).strftime("%Y-%m-%d")
    print(f"{args.dir}  ({len(files)} file(s))\n")
    wrote = kept = drifted = unstamped = 0

    for p in files:
        side = sidecar_path(p)
        size = f"{p.stat().st_size:,}"
        if args.check:
            if not side.exists():
                unstamped += 1
                print(f"  UNSTAMPED {p.name:<44}{size:>16}")
                continue
            try:
                rec = json.loads(side.read_text())
            except ValueError as e:
                drifted += 1
                print(f"  UNREADABLE {p.name:<43}  {e}")
                continue
            claimed = rec.get("sha256")
            if not claimed:
                unstamped += 1
                print(f"  NO CHECKSUM {p.name:<42}{size:>16}")
                continue
            actual = sha256(p)
            if actual == claimed:
                print(f"  ok        {p.name:<44}{size:>16}  "
                      f"{rec.get('provenance', 'fetched')}")
            else:
                drifted += 1
                print(f"  DRIFTED   {p.name:<44}{size:>16}")
                print(f"            sidecar {claimed[:16]}, file "
                      f"{actual[:16]}")
            continue

        if side.exists():
            # A record from fetch_raw carries a real URL and retrieval
            # date. Replacing it with a reconstruction would destroy the
            # one thing that tells the two apart.
            kept += 1
            print(f"  kept      {p.name:<44}{size:>16}  already stamped")
            continue
        rec = record(p, today)
        side.write_text(json.dumps(rec, indent=2, sort_keys=True) + "\n")
        wrote += 1
        ver = rec["version_in_filename"] or "-"
        print(f"  stamped   {p.name:<44}{size:>16}  {ver}")

    print()
    if args.check:
        print(f"{len(files) - drifted - unstamped} verified, {drifted} "
              f"drifted, {unstamped} unstamped")
        if drifted:
            print("\nA DRIFTED file is not the one the release was built "
                  "from. Rebuild, or\nrestore the file the checksum names.")
        if unstamped:
            print("\nAn UNSTAMPED file has no checksum, so a change to it "
                  "cannot be seen.\nRun this without --check to record one.")
        return 1 if (drifted or unstamped) else 0

    print(f"{wrote} stamped, {kept} already had a record")
    if wrote:
        print("\nThese say `provenance: reconstructed` and leave url and "
              "retrieved null:\nthe files were placed by hand and those are "
              "not recoverable. What they\nestablish is a checksum, so the "
              "next change to one of them is visible.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
