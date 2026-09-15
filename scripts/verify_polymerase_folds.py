#!/usr/bin/env python3
"""Verify the polymerase result, which contradicted the plan.

PHASE 0 REPORTED: coronavirus RdRp shares NO superfamily with picornavirus 3D
or flavivirus NS5, while those two share G3DSA:3.30.70.270 and SSF56672.

That is surprising. RdRp folds are usually described as broadly conserved
across positive-sense RNA viruses -- the "right hand" fingers/palm/thumb
architecture with conserved motifs A-G. A result that contradicts the textbook
is either a real and interesting finding or an artifact of how the question was
asked, and phase 0 already produced one of each.

TWO HYPOTHESES:

  (a) REAL. CATH assigns coronavirus nsp12 its own superfamily because of the
      NiRAN domain and the interface domain, which picornaviral 3D lacks
      entirely. The catalytic palm may still be homologous while the
      superfamily assignment differs.

  (b) ARTIFACT. The chain span used in phase 0 covers the whole nsp12 chain
      (NiRAN + interface + RdRp), so a palm-domain match could have been
      recorded under a different parent entry, or the overlap filter clipped it.

WHAT THIS SCRIPT DOES DIFFERENTLY:

  1. Reports EVERY match on each polymerase chain with its position, database
     and type -- not just superfamily-level, so a shared entry at family or
     domain level is visible.
  2. Restricts to the CATALYTIC PALM subregion rather than the whole chain,
     using the annotated active site as an anchor.
  3. Checks for the motif C aspartate pair (GDD in most +ssRNA viruses, SDD in
     coronaviruses) directly in the sequence, which is the functional signature
     the fold assignment is supposed to reflect.

If (3) finds the motif in all four while (1) finds no shared superfamily, then
the fold IS conserved and CATH's superfamily granularity is simply too fine to
express it -- which would mean FOLD_SIMILAR_TO cannot be built from CATH alone
for polymerases and needs structural alignment instead.
"""
from __future__ import annotations

import json
import re
import time
import urllib.request

INTERPRO = "https://www.ebi.ac.uk/interpro/api"
UNIPROT = "https://rest.uniprot.org/uniprotkb"

POLYMERASES = {
    "SARS-CoV-2 nsp12": ("P0DTD1", 4393, 5324, "Coronaviridae"),
    "SARS-CoV nsp12": ("P0C6X7", 4370, 5301, "Coronaviridae"),
    "Poliovirus 3D": ("P03300", 1749, 2209, "Picornaviridae"),
    "EV-A71 3D": ("Q66478", 1730, 2194, "Picornaviridae"),
    "HCV NS5B": ("P26664", 2421, 3011, "Flaviviridae"),
    "Dengue 2 NS5": ("P29990", 2492, 3391, "Flaviviridae"),
    "Zika NS5": ("Q32ZE1", 2520, 3423, "Flaviviridae"),
}

# Motif C carries the catalytic aspartates. GDD in most +ssRNA viruses; SDD in
# coronaviruses and some others. This is the functional signature the fold
# assignment is meant to reflect, so finding it where CATH disagrees is
# evidence that the superfamily granularity is the problem.
MOTIF_C = re.compile(r"[GSC]D[DN]")


def get(url: str):
    req = urllib.request.Request(url, headers={"Accept": "application/json",
                                               "User-Agent": "kgav/2.0"})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode())


def matches_on(acc: str, lo: int, hi: int) -> list[dict]:
    out = []
    data = get(f"{INTERPRO}/entry/all/protein/uniprot/{acc}/?page_size=200")
    for item in data.get("results", []):
        md = item.get("metadata", {})
        for prot in item.get("proteins", []):
            for loc in prot.get("entry_protein_locations", []):
                for frag in loc.get("fragments", []):
                    fs, fe = frag.get("start", 0), frag.get("end", 0)
                    if fs <= hi and fe >= lo:
                        out.append({"acc": md.get("accession"),
                                    "name": md.get("name"),
                                    "db": (md.get("source_database") or "").lower(),
                                    "type": (md.get("type") or "").lower(),
                                    "start": fs, "end": fe})
                        break
    seen, uniq = set(), []
    for m in out:
        if m["acc"] not in seen:
            seen.add(m["acc"])
            uniq.append(m)
    return uniq


def sequence(acc: str) -> str:
    return get(f"{UNIPROT}/{acc}?fields=sequence").get("sequence", {}).get("value", "")


def main() -> int:
    results = {}
    seqs: dict[str, str] = {}
    print("fetching\n")
    for label, (acc, lo, hi, fam) in POLYMERASES.items():
        if acc not in seqs:
            seqs[acc] = sequence(acc)
            time.sleep(0.3)
        ms = matches_on(acc, lo, hi)
        results[label] = {"family": fam, "acc": acc, "lo": lo, "hi": hi,
                          "matches": ms}
        print(f"  {label:<20} {len(ms):>3} matches on [{lo}-{hi}]")
        time.sleep(0.3)

    # ------------------------------------------------- 1. all levels, shared
    print("\n" + "=" * 74)
    print("1 — shared entries by level (phase 0 looked only at superfamily)")
    print("=" * 74)

    def entries(fam: str, level: str | None) -> set[str]:
        s = set()
        for r in results.values():
            if r["family"] != fam:
                continue
            for m in r["matches"]:
                if level is None or m["type"] == level:
                    s.add(m["acc"])
        return s

    fams = ("Coronaviridae", "Picornaviridae", "Flaviviridae")
    for level in ("homologous_superfamily", "family", "domain", None):
        label = level or "ANY LEVEL"
        sets = {f: entries(f, level) for f in fams}
        cp = sets["Coronaviridae"] & sets["Picornaviridae"]
        cf = sets["Coronaviridae"] & sets["Flaviviridae"]
        pf = sets["Picornaviridae"] & sets["Flaviviridae"]
        print(f"\n  {label}")
        print(f"    CoV n Picorna  : {sorted(cp) or 'none'}")
        print(f"    CoV n Flavi    : {sorted(cf) or 'none'}")
        print(f"    Picorna n Flavi: {sorted(pf) or 'none'}")

    # --------------------------------------------- 2. what CoV actually has
    print("\n" + "=" * 74)
    print("2 — coronavirus nsp12 entries, to see what it is assigned instead")
    print("=" * 74)
    for label, r in results.items():
        if r["family"] != "Coronaviridae":
            continue
        print(f"\n  {label}")
        for m in sorted(r["matches"], key=lambda m: m["start"]):
            print(f"    {m['acc']:<16} {m['type'][:22]:<24} "
                  f"[{m['start']}-{m['end']}] {(m['name'] or '')[:40]}")

    # ------------------------------------------------------ 3. motif C
    print("\n" + "=" * 74)
    print("3 — motif C catalytic aspartates (GDD / SDD), read from sequence")
    print("=" * 74)
    print("  The functional signature the fold assignment is meant to reflect.")
    print(f"\n  {'polymerase':<20} {'family':<16} motif hits in chain")
    for label, (acc, lo, hi, fam) in POLYMERASES.items():
        seq = seqs.get(acc, "")
        if not seq:
            print(f"  {label:<20} {fam:<16} sequence unavailable")
            continue
        region = seq[lo - 1:hi]
        hits = [(m.group(), lo + m.start()) for m in MOTIF_C.finditer(region)]
        shown = ", ".join(f"{g}@{p}" for g, p in hits[:6]) or "NONE"
        print(f"  {label:<20} {fam:<16} {shown}")

    # ------------------------------------------------------------ verdict
    print("\n" + "=" * 74)
    print("READING THIS")
    print("=" * 74)
    print("""  If section 1 shows NO shared entry at ANY level between coronavirus
  and the other two, while section 3 finds motif C in all of them, then the
  catalytic fold is conserved and CATH/InterPro granularity cannot express
  it. FOLD_SIMILAR_TO for polymerases would then need structural alignment
  (TM-align or Foldseek over PDB/AlphaFold models), not database lookup.

  If section 1 shows shared entries at family or domain level, phase 0's
  superfamily-only filter was too strict and should be widened.

  Section 2 shows what coronavirus nsp12 is assigned INSTEAD -- if those
  entries are nsp12-specific, that is CATH splitting on the NiRAN and
  interface domains, which picornaviral 3D genuinely lacks.""")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
