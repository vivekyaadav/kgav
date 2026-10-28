#!/usr/bin/env python3
"""Check every MONDO term the register declares against the ontology.

    python scripts/verify_mondo.py
    python scripts/verify_mondo.py --mondo data/raw/mondo/mondo.json

WHY. config/viruses.yaml declares a disease per virus and a
`mondo_source` saying on what basis. Three of those were curated by
hand in 2026-09 and nothing re-read them: a MONDO term can be
deprecated between releases, an xref can become claimed by two terms,
and a typo in a CURIE produces a Disease node that joins to nothing.
ingest_spine emits a CAUSES edge straight from the declared value, so
an unchecked term becomes an unchecked edge.

normalize_bio.MondoResolver existed for exactly this and had NO caller.
Two of the three curated terms recorded a LABEL MATCH as their basis --
"severe acute respiratory syndrome; sole exact match" -- which is the
one method that module bans, in its own words: "'influenza' matches the
disease, the virus genus, and a handful of drug trade names, and a
label-matched graph is one where nobody can tell which." Both terms
happened to be right. The recorded reason for believing them was not.

WHAT IS CHECKED, per virus
  the term exists and is not deprecated
  the xref named in mondo_source resolves, by xref, to that same term
  a term declared with no xref basis is reported, not accepted

WHAT IS NOT. Whether the disease is the RIGHT one for the virus is a
curation judgement and stays with the person who made it; this only
proves the identifier means what the register says it means. The
seasonal coronaviruses share `common cold` and the four dengue
serotypes share `dengue disease`, both because MONDO has no finer term
-- the notes say so and this script does not second-guess them.

Read-only. Writes nothing. Exit 1 if any declared term fails.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav.normalize_bio import MondoResolver

#: A CURIE in mondo_source is the stated basis. Anything with a prefix
#: and a local id counts, because MONDO carries xrefs from many
#: vocabularies and pinning this to DOID would reject a legitimate one.
XREF_RE = re.compile(r"\b([A-Za-z][A-Za-z0-9.]*:[A-Za-z0-9._-]+)\b")
MONDO_RE = re.compile(r"^MONDO:\d+$")

DOWNLOAD = ("https://github.com/monarch-initiative/mondo/releases/latest/"
            "download/mondo.json")


def stated_xrefs(source: str) -> list[str]:
    """CURIEs named in mondo_source, excluding MONDO terms themselves.

    A mondo_source citing the MONDO id it is justifying would be
    circular, so those are not accepted as a basis.
    """
    return [x for x in XREF_RE.findall(source or "")
            if not MONDO_RE.match(x)]


def check(entry: dict, r: MondoResolver) -> tuple[str, str]:
    """(status, detail) for one virus. status "" means it passed."""
    term = entry.get("mondo")
    source = entry.get("mondo_source") or ""
    if not term:
        return "NO TERM", "mondo is null"
    res = r.resolve(term)
    if res.status == "deleted":
        return "DEPRECATED", f"{term} is deprecated in this MONDO release"
    if not res.ok:
        return "UNKNOWN", f"{term} is not a term in this MONDO release"
    label = r.labels.get(term, "")
    xrefs = stated_xrefs(source)
    if not xrefs:
        # A LABEL MATCH IS NOT A BASIS. MondoResolver refuses to resolve
        # by label and this refuses to accept one as justification.
        return "NO XREF BASIS", (f"{term} ({label}) is declared with no "
                                 f"cross-reference in mondo_source")
    for xref in xrefs:
        got = r.resolve(xref)
        if got.ok and got.curie == term:
            return "", f"{term}  {label}   <- {xref}"
        if got.status == "demerged":
            return "DEMERGED", f"{xref}: {got.note}"
    return "XREF DISAGREES", (
        f"{term} ({label}) is declared, but "
        + "; ".join(f"{x} resolves to {r.resolve(x).curie}" for x in xrefs))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--config", type=Path, default=ROOT / "config/viruses.yaml")
    ap.add_argument("--mondo", type=Path,
                    default=ROOT / "data/raw/mondo/mondo.json")
    args = ap.parse_args(argv)

    if not args.mondo.exists():
        print(f"no MONDO ontology at {args.mondo}")
        print(f"  curl -sSL -o {args.mondo} \\\n    {DOWNLOAD}")
        return 2

    r = MondoResolver(args.mondo)
    print(f"{len(r.xref_to_mondo):,} xrefs over {len(r.labels):,} live terms")
    doc = yaml.safe_load(args.config.read_text())
    viruses = doc.get("viruses") or []

    failed, untermed = [], []
    print(f"\n    {'virus':<14} {'verdict':<15} detail")
    for v in viruses:
        status, detail = check(v, r)
        print(f"    {v['name']:<14} {status or 'ok':<15} {detail}")
        if status == "NO TERM":
            untermed.append(v["name"])
        elif status:
            failed.append((v["name"], status, detail))

    print()
    if untermed:
        # NOT A FAILURE. A virus with no curated disease is a gap in the
        # register, and the agent cannot be asked about it by disease
        # name -- but an absent term is honest where a guessed one is
        # not, so this is reported and does not fail the check.
        print(f"  {len(untermed)} virus(es) carry no MONDO term: "
              + ", ".join(untermed))
        print("  They cannot be reached from a disease name. That is a "
              "curation gap,\n  not an error, and guessing a term would be "
              "worse than the gap.")
    if failed:
        print(f"  {len(failed)} DECLARED TERM(S) FAILED:")
        for name, status, detail in failed:
            print(f"    {name:<14} {status}: {detail}")
        print("\n  ingest_spine emits a CAUSES edge straight from these "
              "values, so a\n  failure here is an edge asserting something "
              "the ontology does not.")
        return 1
    print(f"  all {len(viruses) - len(untermed)} declared term(s) verified "
          f"by xref")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
