#!/usr/bin/env python3
"""Phase 0 — verify fold classifications AT CHAIN LEVEL.

WHY THIS IS THE SECOND VERSION. The first queried polyprotein accessions and
compared every domain in the entire replicase. It reported that coronavirus and
flavivirus proteases "share" IPR027417 -- which is P-loop NTPase, i.e. the
HELICASE, present in both polyproteins and nothing to do with either protease.
The comparison was never about proteases at all.

That is the polyprotein problem from v1 chapter 9, recurring in the first script
of v2. Coronavirus replicase is one UniProt entry containing nsp1-16; flavivirus
and picornavirus genomes are likewise single polyproteins. Any question about a
specific protein must be asked of a specific CHAIN.

WHAT CHANGED. For each target we now:
  1. fetch the PRO_ chain feature and its coordinates from UniProt
  2. fetch InterPro matches WITH their positions
  3. keep only matches whose span overlaps the chain
  4. compare folds and catalytic types on that restricted set

ONE RESULT FROM THE FIRST RUN SURVIVES AND MATTERS. SSF50494 -- the trypsin-like
serine protease superfamily -- appeared for both coronavirus and flavivirus
proteases. That is very likely REAL rather than an artifact: coronavirus nsp5 is
chymotrypsin-like in fold, chymotrypsin and trypsin share that superfamily, and
flavivirus NS3 is trypsin-like.

If chain-level analysis confirms it, the plan's premise needs restating:

    The FOLD IS SHARED across all three families.
    What differs is the CATALYTIC RESIDUE -- cysteine vs serine.

That makes the catalytic-type gate the single thing separating them, which is a
stronger justification for the gate than the plan had, and it means fold
similarity ALONE would have recommended cysteine-protease inhibitors for serine
proteases. Exactly the failure the gate exists to prevent.

EXIT CODES
  0  premise confirmed, in whichever form the data supports
  1  premise refuted -- phases 5-10 must not be built on it
  2  insufficient data to decide (chains or annotations missing)
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, field

INTERPRO = "https://www.ebi.ac.uk/interpro/api"
UNIPROT = "https://rest.uniprot.org/uniprotkb"

TARGETS = {
    # --- Coronaviridae 3C-like / main protease (nsp5). Expected CYSTEINE.
    "SARS-CoV-2 3CLpro": {
        "acc": "P0DTD1", "family": "Coronaviridae", "role": "protease",
        "expect": "cysteine",
        "chain": ("3C-like proteinase", "nsp5", "3C-like protease")},
    "SARS-CoV 3CLpro": {
        "acc": "P0C6X7", "family": "Coronaviridae", "role": "protease",
        "expect": "cysteine", "chain": ("3C-like proteinase", "nsp5")},
    "MERS-CoV 3CLpro": {
        "acc": "K9N7C7", "family": "Coronaviridae", "role": "protease",
        "expect": "cysteine", "chain": ("3C-like proteinase", "nsp5")},

    # --- Picornaviridae 3C protease. Expected CYSTEINE. POSITIVE CONTROL.
    "EV-A71 3C": {
        "acc": "Q66478", "family": "Picornaviridae", "role": "protease",
        "expect": "cysteine", "chain": ("Picornain 3C", "Protease 3C", "3C")},
    "Poliovirus 3C": {
        "acc": "P03300", "family": "Picornaviridae", "role": "protease",
        "expect": "cysteine", "chain": ("Picornain 3C", "Protease 3C")},
    "HAV 3C": {
        "acc": "P08617", "family": "Picornaviridae", "role": "protease",
        "expect": "cysteine", "chain": ("Picornain 3C", "Protease 3C")},

    # --- Flaviviridae NS3 protease. Expected SERINE. NEGATIVE CONTROL.
    "Dengue 2 NS3": {
        "acc": "P29990", "family": "Flaviviridae", "role": "protease",
        "expect": "serine", "chain": ("Serine protease NS3", "NS3")},
    "Zika NS3": {
        "acc": "Q32ZE1", "family": "Flaviviridae", "role": "protease",
        "expect": "serine", "chain": ("Serine protease NS3", "NS3")},
    "West Nile NS3": {
        "acc": "P06935", "family": "Flaviviridae", "role": "protease",
        "expect": "serine", "chain": ("Serine protease NS3", "NS3")},
    "HCV NS3": {
        "acc": "P26664", "family": "Flaviviridae", "role": "protease",
        "expect": "serine", "chain": ("Serine protease/helicase NS3", "NS3")},

    # --- Polymerases. Second prediction: RdRp folds cross families more
    #     broadly than protease folds do.
    "SARS-CoV-2 RdRp": {
        "acc": "P0DTD1", "family": "Coronaviridae", "role": "polymerase", "expect": None,
        "chain": ("RNA-directed RNA polymerase", "nsp12")},
    "HCV NS5B": {
        "acc": "P26664", "family": "Flaviviridae", "role": "polymerase", "expect": None,
        "chain": ("RNA-directed RNA polymerase", "NS5B")},
    "Dengue 2 NS5": {
        "acc": "P29990", "family": "Flaviviridae", "role": "polymerase", "expect": None,
        "chain": ("RNA-directed RNA polymerase NS5", "NS5")},
    "Poliovirus 3D": {
        "acc": "P03300", "family": "Picornaviridae", "role": "polymerase", "expect": None,
        "chain": ("RNA-directed RNA polymerase", "3D")},
}

# Superfamily-level databases only. Family-level entries are virus-specific by
# construction and can never cross a family boundary, so including them would
# guarantee a null result and make the test unfalsifiable.
SUPERFAMILY_DBS = {"cathgene3d", "ssf", "superfamily"}

CYS_MARKERS = ("cysteine peptidase", "cysteine protease", "peptidase c3",
               "peptidase c30", "picornain", "thiol protease", "3c-like")
SER_MARKERS = ("serine peptidase", "serine protease", "peptidase s7",
               "peptidase s29", "peptidase s31", "trypsin-like serine",
               "chymotrypsin-like serine")


@dataclass
class Chain:
    label: str
    accession: str
    family: str
    role: str
    expect: str | None
    chain_name: str = ""
    start: int = 0
    end: int = 0
    matches: list[dict] = field(default_factory=list)
    note: str = ""

    @property
    def found(self) -> bool:
        return self.end > self.start > 0

    def superfamilies(self) -> set[str]:
        return {m["accession"] for m in self.matches
                if (m.get("db") or "").lower() in SUPERFAMILY_DBS}

    def catalytic(self) -> str | None:
        text = " ".join((m.get("name") or "") for m in self.matches).lower()
        cys = any(k in text for k in CYS_MARKERS)
        ser = any(k in text for k in SER_MARKERS)
        if cys and not ser:
            return "cysteine"
        if ser and not cys:
            return "serine"
        if cys and ser:
            return "ambiguous"
        return None


def _get(url: str, retries: int = 3) -> dict | None:
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url, headers={"Accept": "application/json",
                              "User-Agent": "kgav-phase0/2.0"})
            with urllib.request.urlopen(req, timeout=90) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if attempt == retries - 1:
                raise
        except Exception:
            if attempt == retries - 1:
                raise
        time.sleep(1.5 * (attempt + 1))
    return None


def find_chain(acc: str, patterns: tuple[str, ...]) -> tuple[str, int, int, str]:
    """Locate a PRO_ chain by description.

    Patterns are tried in order, longest-intent first, because chain names nest:
    "NS3" is a substring of "Serine protease/helicase NS3". Fifth instance of
    this hazard in the project.
    """
    data = _get(f"{UNIPROT}/{acc}?fields=ft_chain,ft_peptide")
    if not data:
        return "", 0, 0, "accession not found"
    feats = [f for f in data.get("features", [])
             if f.get("type") in ("Chain", "Peptide")]
    if not feats:
        return "", 0, 0, "entry has no chain features (TrEMBL-only?)"
    lowered = [(f, (f.get("description") or "").lower()) for f in feats]
    for pat in patterns:
        p = pat.lower()
        for f, desc in lowered:
            if p in desc:
                loc = f.get("location", {})
                s = (loc.get("start") or {}).get("value")
                e = (loc.get("end") or {}).get("value")
                if s and e:
                    return f.get("description", pat), int(s), int(e), ""
    names = "; ".join(sorted({d for _f, d in lowered})[:6])
    return "", 0, 0, f"no chain matched {patterns}; entry has: {names}"


def fetch_matches(acc: str, start: int, end: int) -> list[dict]:
    """InterPro matches overlapping [start, end] on this accession."""
    out: list[dict] = []
    data = _get(f"{INTERPRO}/entry/all/protein/uniprot/{acc}/?page_size=200")
    if not data:
        return out
    for item in data.get("results", []):
        md = item.get("metadata", {})
        for prot in item.get("proteins", []):
            for loc in prot.get("entry_protein_locations", []):
                hit = False
                for frag in loc.get("fragments", []):
                    fs, fe = frag.get("start", 0), frag.get("end", 0)
                    if fs <= end and fe >= start:
                        out.append({"accession": md.get("accession"),
                                    "name": md.get("name"),
                                    "db": md.get("source_database"),
                                    "type": (md.get("type") or "").lower(),
                                    "start": fs, "end": fe})
                        hit = True
                        break
                if hit:
                    break
    seen, uniq = set(), []
    for m in out:
        if m["accession"] not in seen:
            seen.add(m["accession"])
            uniq.append(m)
    return uniq


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data/raw/folds/phase0_chains.json")
    ap.add_argument("--offline", action="store_true")
    args = ap.parse_args()
    out = pathlib.Path(args.out)

    if args.offline and out.exists():
        chains = {k: Chain(**v) for k, v in json.loads(out.read_text()).items()}
    else:
        chains = {}
        for label, spec in TARGETS.items():
            print(f"{label:<22} ", end="", flush=True)
            c = Chain(label=label, accession=spec["acc"], family=spec["family"],
                      role=spec["role"], expect=spec["expect"])
            name, s, e, note = find_chain(spec["acc"], spec["chain"])
            c.chain_name, c.start, c.end, c.note = name, s, e, note
            if c.found:
                c.matches = fetch_matches(spec["acc"], s, e)
                print(f"{name[:38]:<40} [{s}-{e}] {len(c.matches)} matches")
            else:
                print(f"NOT FOUND — {note[:72]}")
            chains[label] = c
            time.sleep(0.4)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps({k: v.__dict__ for k, v in chains.items()},
                                  indent=2))
        print(f"\nwrote {out}")

    missing = [c for c in chains.values() if not c.found]
    usable = [c for c in chains.values() if c.found]
    if missing:
        print(f"\n{len(missing)} chain(s) not located:")
        for c in missing:
            print(f"  {c.label}: {c.note[:110]}")

    print("\n" + "=" * 74)
    print("CHECK 1 — catalytic type, from CHAIN-level annotation only")
    print("=" * 74)
    print(f"{'chain':<22} {'family':<16} {'expect':<9} {'observed':<10} ok")
    cat_fail = 0
    for c in usable:
        if c.role != "protease":
            continue
        obs = c.catalytic()
        ok = obs == c.expect
        cat_fail += 0 if ok else 1
        print(f"{c.label:<22} {c.family:<16} {c.expect or '-':<9} "
              f"{obs or 'unknown':<10} {'yes' if ok else 'NO'}")

    def sfs(family: str, role: str) -> set[str]:
        s: set[str] = set()
        for c in usable:
            if c.family == family and c.role == role:
                s |= c.superfamilies()
        return s

    cov, pic, fla = (sfs(f, "protease")
                     for f in ("Coronaviridae", "Picornaviridae", "Flaviviridae"))
    print("\n" + "=" * 74)
    print("CHECK 2 — protease fold, superfamily level, chain-restricted")
    print("=" * 74)
    for name, s in (("Coronaviridae", cov), ("Picornaviridae", pic),
                    ("Flaviviridae", fla)):
        print(f"  {name:<16}: {sorted(s) or 'none'}")
    shared_pic, shared_fla = cov & pic, cov & fla
    print(f"\n  CoV n Picorna: {sorted(shared_pic) or 'NONE'}")
    print(f"  CoV n Flavi  : {sorted(shared_fla) or 'none'}")

    rcov, rpic, rfla = (sfs(f, "polymerase")
                        for f in ("Coronaviridae", "Picornaviridae",
                                  "Flaviviridae"))
    print("\n" + "=" * 74)
    print("CHECK 3 — polymerase fold (expect broader sharing than protease)")
    print("=" * 74)
    print(f"  CoV n Picorna  : {sorted(rcov & rpic) or 'none'}")
    print(f"  CoV n Flavi    : {sorted(rcov & rfla) or 'none'}")
    print(f"  Picorna n Flavi: {sorted(rpic & rfla) or 'none'}")

    print("\n" + "=" * 74)
    print("VERDICT")
    print("=" * 74)
    if len(usable) < 6:
        print(f"INSUFFICIENT DATA: only {len(usable)} of {len(TARGETS)} chains "
              "resolved.\nFix the chain patterns above before drawing any "
              "conclusion.")
        return 2

    cat_separates = cat_fail == 0

    if shared_pic and shared_fla and cat_separates:
        print("PREMISE RESTATED, AND THE GATE IS NECESSARY.")
        print(f"  All three families share protease superfamilies: "
              f"{sorted(shared_pic & shared_fla) or sorted(shared_pic)}")
        print("  Catalytic type SEPARATES them cleanly (cysteine vs serine).")
        print()
        print("  Consequence: fold similarity ALONE would recommend")
        print("  cysteine-protease inhibitors for serine proteases. The")
        print("  catalytic-type gate is not a refinement -- it is the only")
        print("  thing standing between the breadth score and that error.")
        print()
        print("  The falsifiable prediction holds, but its MECHANISM is")
        print("  catalytic chemistry, not fold divergence. Restate it that way")
        print("  in the plan, then proceed to phase 1.")
        return 0

    if shared_pic and not shared_fla and cat_separates:
        print("PREMISE CONFIRMED AS ORIGINALLY STATED.")
        print(f"  CoV and Picorna proteases share {sorted(shared_pic)}")
        print("  CoV and Flavi proteases share no superfamily.")
        print("  Catalytic types are as expected. Proceed to phase 1.")
        return 0

    print("PREMISE NOT CONFIRMED:")
    if not shared_pic:
        print("  - CoV and Picorna proteases share NO superfamily: the")
        print("    positive control fails, so fold conservation does not")
        print("    predict the transfer the plan is built on.")
    if not cat_separates:
        print(f"  - {cat_fail} protease(s) do not show the expected catalytic")
        print("    type at chain level, so the gate cannot be trusted.")
    print("\nPhases 5-10 must not be built on the conservation score.")
    print("This is the outcome the phase exists to detect.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
