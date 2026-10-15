#!/usr/bin/env python3
"""Phase 2 driver — build the fold layer across three virus families.

WHAT IS INGESTED AND WHAT IS SCORED ARE DIFFERENT SETS.

Everything below is ingested: proteases, polymerases, helicases,
methyltransferases. Only proteases and polymerases are marked score-eligible.

The reason is sequencing rather than squeamishness. Phase 4 -- the HCV
retrospective -- tests the conservation premise on PROTEASES. If it fails, the
premise is wrong and phases 5-10 should not be built. Adding weak-discrimination
targets to the score before that test is building on an untested foundation, and
the v1 result already showed what diffuse evidence does: specific mechanistic
routes reached AUC 0.826 while diffuse host-directed ones sat below chance.

Helicases and methyltransferases land in the gate's MOTIF mode, which has no
catalytic discriminator and is always marked weak. They also have far fewer
known inhibitors, so they add targets without adding compounds to transfer
from. Ingesting them now costs an afternoon; scoring on them now costs the
validity of phase 4.

THE NS3 SPLIT. Flavivirus NS3 is a protease AND a helicase in one chain --
UniProt describes it as "Serine protease/helicase NS3" and phase 0 saw 22-26
InterPro matches on it, against 7-9 for a coronavirus or picornavirus protease.
Treating it as one target lets the protease comparison inherit helicase fold
matches, which is the polyprotein problem one level down.

THIS PARAGRAPH USED TO END "Both are declared separately below, with the
boundary taken from the domain annotation rather than assumed." Neither
clause was ever true. There is one NS3 spec per flavivirus, not two, and
until this commit nothing in the codebase read a Domain annotation at all --
grep for ft_domain returned zero. The NS3pro chains resolve to 618-631
residues, the whole protease-AND-helicase product, which is why they carry
22-26 classes against 7-12 for a chain that is one enzyme. Third instance of
a comment outliving its code, after selectivity.py's "WORST is kept" and
catalytic_machinery's "last catalytic residue" paragraph -- and the worst of
the three, because this one described a design that was planned and never
built, in the present tense.

WHAT IS TRUE NOW: ChainSpec takes an optional `domain`, which narrows the
analysis window -- the region the class fetch and the catalytic rule read --
to a named UniProt Domain inside the matched chain, while the chain's PRO_ id
stays the node identity. SARS-CoV-2 PLpro uses it. The six NS3pro chains do
not yet: narrowing them needs the domain's name as UniProt spells it, which
is measured by the run (see the analysis-window report) rather than recalled.
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kgav.emit import Emit
from kgav.folds import (DERIVED_FIELDS, ChainSpec, ResolvedChain, gate,
                        resolve)
from kgav.schema import load_schema, print_capped, report_violations

# Roles whose comparisons may enter a score. Others are ingested and stored but
# excluded until phase 4 validates the premise.
SCORE_ELIGIBLE_ROLES = {"protease", "polymerase"}

# Chain patterns are MOST SPECIFIC FIRST. "NS3" is a substring of
# "Serine protease/helicase NS3"; "nsp1" is a substring of "nsp10". Sixth
# instance of this hazard in the project.
#
# ---------------------------------------------------------------------------
# THE FLAVIVIRUS ACCESSIONS BELOW ARE STRAIN-TAXON ACCESSIONS, AND THE
# virus_taxon FIELD BESIDE THEM IS THE SPECIES. That reassignment is
# deliberate -- cross-viral metapaths must resolve against one taxon per virus
# -- but it was undocumented, which is how this class of bug survives.
# Measured against UniProt 2026-09-17:
#
#   P29990  Dengue 2   filed under 31634  (strain Thailand/16681/1984)
#                      declared here as   NCBITaxon:11060  (species)
#   P27395  JEV        filed under 11073  (strain SA-14)
#                      declared here as   NCBITaxon:11072  (species)
#   P03314  Yellow fev filed under 11090  (strain 17D vaccine)
#                      declared here as   NCBITaxon:11089  (species)
#   P26664  HCV        filed under 11104  (genotype 1a isolate 1)
#                      declared here as   NCBITaxon:11103  (species)
#
# Q32ZE1 (Zika, 64320) and P06935 (West Nile, 11082) are the only two curated
# AT the species taxon, so they need no reassignment.
#
# Querying the species taxon for any of the four above returns no curated
# polyprotein and therefore no chains -- the MERS trap, which cost this
# project Mpro and RdRp once already. Every strain taxon here is now also
# registered in config/viruses.yaml under `surveyed:`, with chain_source
# naming the accession and the taxon it really lives under, and
# kgav.virus_register.validate enforces that the pairing is declared rather
# than discovered. When this list moves into the config mechanism, that
# register is the source to read; until then the two must agree.
# ---------------------------------------------------------------------------
SPECS: list[ChainSpec] = [
    # ---------------------------------------------------- Coronaviridae
    ChainSpec("SARS-CoV-2 3CLpro", "P0DTD1",
              ("3C-like proteinase", "nsp5"), "Coronaviridae", "protease",
              "NCBITaxon:2697049"),
    ChainSpec("SARS-CoV 3CLpro", "P0C6X7",
              ("3C-like proteinase", "nsp5"), "Coronaviridae", "protease",
              "NCBITaxon:694009"),
    ChainSpec("MERS-CoV 3CLpro", "K9N7C7",
              ("3C-like proteinase", "nsp5"), "Coronaviridae", "protease",
              "NCBITaxon:1335626"),
    ChainSpec("HCoV-229E 3CLpro", "P0C6X1",
              ("3C-like proteinase", "nsp5"), "Coronaviridae", "protease",
              "NCBITaxon:11137"),
    # THE CHAIN IS nsp3, 819-2763 -- nine domains and 1,945 residues, of
    # which the protease is one. Measured consequences of analysing all of
    # it: 60 InterPro classes against 7-12 for a chain that is only a
    # protease, and an active-site list carrying N1062 and G1070, the
    # macrodomain's ADP-ribosylhydrolase sites, 600 residues away. The
    # catalytic rule returned cysteine anyway, but only because neither of
    # those is a C or an S -- a right answer that did not depend on the
    # window being right.
    #
    # "Peptidase C16" 1634-1898 is UniProt's own annotation on P0DTD1,
    # read off the entry rather than recalled.
    ChainSpec("SARS-CoV-2 PLpro", "P0DTD1",
              ("Papain-like protease", "papain-like proteinase"),
              "Coronaviridae", "protease", "NCBITaxon:2697049",
              domain="Peptidase C16"),
    # ALL SEVEN POLYMERASES NARROW TO "RdRp catalytic", and that name was
    # chosen because it is the ONLY domain UniProt annotates on every one of
    # them -- not because it is the prettiest boundary. The 2026-10-15 run
    # printed each entry's domains; intersecting the seven name sets leaves
    # exactly {"RdRp catalytic"}. The broad coronavirus domain, "Nsp12
    # RNA-dependent RNA polymerase" (568 aa), has no counterpart in the
    # flaviviruses or picornaviruses, so using it would compare a 568 aa
    # window against a 116 aa one -- unlike regions, scored as if alike.
    #
    #   SARS-CoV-2 RdRp  5004-5166  163 aa   sites S5151 D5152 D5153 inside
    #   SARS-CoV RdRp    4981-5143  163 aa   sites S5128 D5129 D5130 inside
    #   HCV NS5B         2634-2752  119 aa   no annotated sites
    #   Dengue 2 NS5     3020-3169  150 aa   its 4 sites are OUTSIDE, below
    #   Zika NS5         3045-3195  151 aa   its 4 sites are OUTSIDE, below
    #   Poliovirus 3D    1975-2090  116 aa   no annotated sites
    #   EV-A71 3D        1958-2073  116 aa   no annotated sites
    #
    # What the coronavirus windows drop is NiRAN, a separate
    # nucleotidyltransferase with its own druggable site, and the Nsp12
    # interface region.
    ChainSpec("SARS-CoV-2 RdRp", "P0DTD1",
              ("RNA-directed RNA polymerase", "nsp12"), "Coronaviridae",
              "polymerase", "NCBITaxon:2697049", domain="RdRp catalytic"),
    ChainSpec("SARS-CoV RdRp", "P0C6X7",
              ("RNA-directed RNA polymerase", "nsp12"), "Coronaviridae",
              "polymerase", "NCBITaxon:694009", domain="RdRp catalytic"),
    ChainSpec("SARS-CoV-2 helicase", "P0DTD1",
              ("Helicase", "nsp13"), "Coronaviridae", "helicase",
              "NCBITaxon:2697049"),

    # ----------------------------------------------------- Flaviviridae
    # NS3 IS ONE CHAIN HOLDING TWO OR THREE ENZYMES, and the spec matched all
    # of it: 617-631 residues and 22-26 InterPro classes, against 7-12 for a
    # chain that is one enzyme. The comment here used to say the protease and
    # helicase were "declared separately"; there is one spec per flavivirus,
    # below, and there always was.
    #
    # The WINDOW narrows instead, so the chain's PRO_ id stays the node
    # identity. Domain names are UniProt's, read off the 2026-10-15 run:
    #
    #   Dengue 2      Peptidase S7   1476-1653   + Helicase ATP-binding, C-terminal
    #   Zika          Peptidase S7   1499-1676   + Helicase ATP-binding, C-terminal
    #   West Nile     Peptidase S7   1502-1679   + Helicase ATP-binding, C-terminal
    #   HCV           Peptidase S29  1027-1208   + Helicase ATP-binding
    #   Yellow fever  Peptidase S7   1485-1665   + Helicase ATP-binding, C-terminal
    #   JEV           Peptidase S7   1505-1682   + Helicase ATP-binding, C-terminal
    #
    # CHECKED BEFORE CHANGING: all 18 catalytic residues -- six triads -- lie
    # inside the protease domain, so no catalytic type can move. What moves is
    # the class list, which is the point: the gate compares on classes.
    ChainSpec("Dengue 2 NS3pro", "P29990",
              ("Serine protease NS3", "NS3"), "Flaviviridae", "protease",
              "NCBITaxon:11060", domain="Peptidase S7"),
    ChainSpec("Zika NS3pro", "Q32ZE1",
              ("Serine protease NS3", "NS3"), "Flaviviridae", "protease",
              "NCBITaxon:64320", domain="Peptidase S7"),
    ChainSpec("West Nile NS3pro", "P06935",
              ("Serine protease/Helicase NS3", "NS3"), "Flaviviridae",
              "protease", "NCBITaxon:11082", domain="Peptidase S7"),
    # HCV is MEROPS S29, the five flaviviruses are S7. Different families,
    # both chymotrypsin-fold serine proteases -- not a typo to normalise.
    ChainSpec("HCV NS3pro", "P26664",
              ("Serine protease/helicase NS3", "NS3"), "Flaviviridae",
              "protease", "NCBITaxon:11103", domain="Peptidase S29"),
    ChainSpec("Yellow fever NS3pro", "P03314",
              ("Serine protease NS3", "NS3"), "Flaviviridae", "protease",
              "NCBITaxon:11089", domain="Peptidase S7"),
    ChainSpec("JEV NS3pro", "P27395",
              ("Serine protease NS3", "NS3"), "Flaviviridae", "protease",
              "NCBITaxon:11072", domain="Peptidase S7"),
    ChainSpec("HCV NS5B", "P26664",
              ("RNA-directed RNA polymerase", "NS5B"), "Flaviviridae",
              "polymerase", "NCBITaxon:11103", domain="RdRp catalytic"),
    # NS5 LOSES ALL FOUR OF ITS ANNOTATED ACTIVE SITES, DELIBERATELY. The
    # run placed every one of them -- K2552, D2637, K2672, E2708 for Dengue,
    # K2577, D2662, K2698, E2734 for Zika -- inside "mRNA cap 0-1 NS5-type
    # MT". They are the METHYLTRANSFERASE's, not the polymerase's. A
    # polymerase window that excludes a capping enzyme's residues is more
    # correct, not less, and nothing changes in the gate: NS5 reads "none"
    # with sites nobody could use and "none" with no sites at all.
    ChainSpec("Dengue 2 NS5", "P29990",
              ("RNA-directed RNA polymerase NS5", "NS5"), "Flaviviridae",
              "polymerase", "NCBITaxon:11060", domain="RdRp catalytic"),
    ChainSpec("Zika NS5", "Q32ZE1",
              ("RNA-directed RNA polymerase NS5", "NS5"), "Flaviviridae",
              "polymerase", "NCBITaxon:64320", domain="RdRp catalytic"),

    # --------------------------------------------------- Picornaviridae
    ChainSpec("EV-A71 3C", "Q66478",
              ("Picornain 3C", "Protease 3C"), "Picornaviridae", "protease",
              "NCBITaxon:39054"),
    ChainSpec("Poliovirus 3C", "P03300",
              ("Picornain 3C", "Protease 3C"), "Picornaviridae", "protease",
              "NCBITaxon:12080"),
    ChainSpec("HAV 3C", "P08617",
              ("Picornain 3C", "Protease 3C"), "Picornaviridae", "protease",
              "NCBITaxon:12092"),
    ChainSpec("CVB3 3C", "P03313",
              ("Picornain 3C", "Protease 3C"), "Picornaviridae", "protease",
              "NCBITaxon:12072"),
    ChainSpec("Rhinovirus A 3C", "P07210",
              ("Picornain 3C", "Protease 3C"), "Picornaviridae", "protease",
              "NCBITaxon:147711"),
    ChainSpec("Poliovirus 3D", "P03300",
              ("RNA-directed RNA polymerase", "3D"), "Picornaviridae",
              "polymerase", "NCBITaxon:12080", domain="RdRp catalytic"),
    ChainSpec("EV-A71 3D", "Q66478",
              ("RNA-directed RNA polymerase", "3D"), "Picornaviridae",
              "polymerase", "NCBITaxon:39054", domain="RdRp catalytic"),
]


# --------------------------------------------------------------- the cache
# THE CACHE HOLDS NETWORK PAYLOAD, NEVER A VERDICT.
#
# It used to hold the whole ResolvedChain, derived fields included, and
# `--offline` rebuilt the object straight from JSON. `catalytic_machinery`
# ran only on the network path, so a change to the rule had no effect
# offline: the run that was supposed to show the widened composition rule
# printed "composition 14, none 9, annotated 1" -- identical to the run
# before it, because it WAS the run before it, replayed.
#
# The guard below catches a field the cache LACKS. Nothing caught a field
# the cache HAD and should not have, which is the worse of the two: a
# missing field at least defaults to something visibly wrong, while a stale
# verdict reads as a fresh one.
#
# Writer and reader are module-level so a test can round-trip through the
# code the script actually runs, rather than a copy of it.
def cache_payload(chains: list) -> list[dict]:
    """What gets written: inputs only, verdicts dropped."""
    return [{**{k: v for k, v in c.__dict__.items()
                if k != "spec" and k not in DERIVED_FIELDS},
             "spec": c.spec.__dict__} for c in chains]


def missing_inputs(raw: list[dict]) -> list[str]:
    """Input fields ResolvedChain carries that this cache file lacks.

    Derived from the dataclass, not from a list someone has to remember to
    extend: adding an input field invalidates an old cache automatically.
    DERIVED_FIELDS are excluded because they are no longer stored, so their
    absence is correct rather than stale.
    """
    if not raw:
        return []
    required = ({f.name for f in dataclasses.fields(ResolvedChain)}
                - {"spec"} - DERIVED_FIELDS)
    return sorted(required - set().union(*(set(d) for d in raw)))


def load_cache(raw: list[dict]) -> list:
    """Rebuild chains, recomputing every derived field.

    A stored verdict is dropped rather than trusted, so a cache written by
    an older rule cannot survive into a run of the newer one.
    """
    chains = []
    for d in raw:
        d = dict(d)
        spec = ChainSpec(**d.pop("spec"))
        c = ResolvedChain(spec=spec,
                          **{k: v for k, v in d.items()
                             if k not in DERIVED_FIELDS})
        c.redecide()
        chains.append(c)
    return chains


def window_label(c) -> str:
    """How a chain's analysed region is printed.

    The per-chain line used to print the CHAIN's coordinates beside the
    WINDOW's class count: "[1476-2093] ... 11 classes" when the 11 came
    from 1476-1653. Both numbers are shown when they differ, because the
    chain is what supplies the node id and the window is what was read.
    """
    lo, hi = c.window
    if c.window_is_chain:
        return f"[{lo}-{hi}]"
    return f"[{lo}-{hi} of {c.start}-{c.end}]"


def active_site_rows(chains) -> list[tuple[str, str, str, str]]:
    """(label, role, residue, domain) for every score-eligible active site.

    EVERY ONE. The filter used to require `c.domains` as well, which
    silently dropped HAV 3C -- three catalytic residues, zero annotated
    Domains -- from the one report whose job is to say where each active
    site sits. A table that omits the case it cannot explain is the failure
    this project keeps finding: correct about what it looks at, silent
    about what it does not.

    All containing domains are named, not the first: UniProt nests, and
    SARS-CoV-2 RdRp's S5151 is in both "Nsp12 RNA-dependent RNA polymerase"
    and "RdRp catalytic".
    """
    rows = []
    for c in chains:
        if c.spec.role not in SCORE_ELIGIBLE_ROLES or not c.active_sites:
            continue
        for pos, aa, _d in sorted(tuple(x) for x in c.active_sites):
            where = [dd for ds, de, dd in (c.domains or []) if ds <= pos <= de]
            if where:
                dom = " + ".join(where)
            elif not c.domains:
                dom = "(entry annotates no Domain)"
            else:
                dom = "(no annotated domain)"
            rows.append((c.spec.label, c.spec.role, f"{aa}{pos}", dom))
    return rows


def allowed_pairs(results) -> list:
    """The transfers that were made."""
    return [r for r in results if r[2].allowed]


def refused_with_shared_fold(results) -> list:
    """Refusals where the two chains DO share a class.

    These carry the gate's justification: a pair the fold calls alike and
    the nucleophile calls different. A refusal with no shared class is the
    uninteresting kind -- nothing wanted to join them in the first place.
    """
    return [r for r in results
            if not r[2].allowed and r[2].mode == "nucleophile" and r[2].shared]


def class_census(pairs) -> dict:
    """mode -> (number of pairs, Counter of class -> pairs carrying it).

    Counts whatever it is handed, so the same arithmetic serves the
    transfers and the refusals instead of being written twice.

    A class carried by EVERY pair in a bucket is a constant. It cannot
    distinguish one pair from another, so it is not evidence about any
    particular transfer -- only evidence that both ends are the kind of
    enzyme the mode is for. IPR007094 is "RNA-directed RNA polymerase,
    catalytic domain": two polymerases sharing it says they are both
    polymerases.

    Counted rather than asserted. The header of src/kgav/folds.py asserted
    for months that all three protease families share SSF50494, and the
    first run that printed every shared class instead of the first two
    showed SARS-CoV-2 and SARS-CoV 3CLpro do not carry it at all.
    """
    out: dict = {}
    for _a, _b, g in pairs:
        n, c = out.setdefault(g.mode, [0, Counter()])
        out[g.mode][0] = n + 1
        for acc in g.shared:
            c[acc] += 1
    return {m: (n, c) for m, (n, c) in out.items()}


def class_coverage(pairs) -> dict:
    """mode -> {class: [labels of participating chains that LACK it]}.

    THE COUNTS ALONE INVITE A WRONG READING, and I made it one message
    after writing them. Seeing "IPR043502 10/16" I called the class
    discriminating. It is not. Every count below the maximum in the
    2026-10-15 run is a hole in annotation coverage:

        SSF50494          10/20 allowed   absent from both SARS 3CLpros
        G3DSA:2.40.10.10  18/54 refused   absent from 4 of the 6 NS3pros
        SSF56672           6/16 allowed   absent from both coronavirus RdRps
        IPR001205          4/16 allowed   absent from the 3 flavivirus NS5s

    A member database either has a hit on an entry or it does not, and a
    pair count falls mechanically with every chain missing the class. Read
    as conservation it says the opposite of what it means: SSF50494 at
    10/20 looks like a class that distinguishes some coronavirus proteases
    from others, when it is a class SUPERFAMILY did not annotate on two
    entries.

    A count is evidence about chemistry only where coverage is complete.
    """
    seen: dict = {}
    for a, b, g in pairs:
        m = seen.setdefault(g.mode, {"chains": {}, "classes": set()})
        for c in (a, b):
            m["chains"][c.spec.label] = {x["accession"] for x in c.classes}
        m["classes"] |= set(g.shared)
    return {mode: {acc: sorted(lab for lab, cl in v["chains"].items()
                               if acc not in cl)
                   for acc in sorted(v["classes"])}
            for mode, v in seen.items()}


def print_census(title: str, census: dict, coverage: dict | None = None,
                 *, note: str = "") -> None:
    """Each count with the reason it is not the maximum, where there is one.

    A bare count reads as conservation and is usually coverage. Naming the
    chains that lack the class puts the cause beside the number instead of
    leaving the reader -- me, last time -- to supply a flattering one.
    """
    if not census:
        return
    coverage = coverage or {}
    print(f"\n  {title}")
    for mode in sorted(census):
        n, counts = census[mode]
        universal = [a for a, k in counts.items() if k == n]
        print(f"    {mode} mode — {n} pair(s), {len(counts)} distinct class(es)")
        for acc, k in counts.most_common():
            if k == n:
                why = "   <- in every pair, so it ranks nothing"
            else:
                lacking = coverage.get(mode, {}).get(acc, [])
                why = ("   <- absent from " + ", ".join(lacking)
                       if lacking else "")
            print(f"      {acc:<22} {k:>3}/{n}{why}")
        if universal and len(universal) == len(counts):
            print(f"      EVERY shared class in {mode} mode is universal: all "
                  f"{n} pairs carry\n      exactly the same fold evidence, so "
                  f"the fold cannot prefer one over\n      another. "
                  f"Eligibility, not similarity.")
        partial = [a for a, k in counts.items() if k < n]
        if partial and all(coverage.get(mode, {}).get(a) for a in partial):
            print("      Every count below the maximum here is explained by a "
                  "chain lacking the\n      class, which is annotation "
                  "coverage and not conservation. None of\n      these spreads "
                  "is evidence about chemistry.")
    if note:
        print(note)


def walkable_transfers(results, spine_ids) -> tuple[list, list]:
    """(walkable, unwalkable) allowed transfers.

    A FOLD_SIMILAR_TO edge is Protein->Protein, so it needs BOTH endpoints
    in the spine. A transfer whose comparator the spine never built is not
    a weak hop -- it is no hop. Separated here rather than counted with the
    rest, because "36 allowed" as a reach figure would be unreproducible:
    M8 walking this graph would find 2.
    """
    walkable, unwalkable = [], []
    for a, b, g in allowed_pairs(results):
        both = a.node_id in spine_ids and b.node_id in spine_ids
        (walkable if both else unwalkable).append((a, b, g))
    return walkable, unwalkable


def blocking_chains(unwalkable, spine_ids) -> Counter:
    """How many transfers each un-built comparator costs."""
    out: Counter = Counter()
    for a, b, _g in unwalkable:
        for c in (a, b):
            if c.node_id not in spine_ids:
                out[c.spec.label] += 1
    return out


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=root / "data/releases/v0.2-folds")
    ap.add_argument("--cache", type=Path,
                    default=root / "data/raw/folds/phase2_chains.json")
    ap.add_argument("--spine", type=Path,
                    default=root / "data/releases/v0.1-spine",
                    help="checked for the Protein nodes this layer references")
    ap.add_argument("--offline", action="store_true")
    args = ap.parse_args()

    # ------------------------------------------------------------- resolve
    if args.offline and args.cache.exists():
        raw = json.loads(args.cache.read_text())
        missing = missing_inputs(raw)
        if missing:
            print(f"{args.cache} predates {', '.join(missing)}.\n"
                  f"  ResolvedChain would default those, and a default is "
                  f"not a measurement -- the\n  last time this happened the "
                  f"report said 'none 24' beside a table saying\n  "
                  f"'cysteine 4'. Re-run without --offline to refresh it.",
                  file=sys.stderr)
            return 2
        chains = load_cache(raw)
        print(f"loaded {len(chains)} chains from cache")
    else:
        chains = []
        for spec in SPECS:
            print(f"{spec.label:<22} ", end="", flush=True)
            r = resolve(spec)
            if r.ok:
                # THE WINDOW, not the chain. This line printed the chain's
                # coordinates beside the WINDOW's class count, so after
                # narrowing it read "[1476-2093] ... 11 classes" when the 11
                # came from 1476-1653. Two numbers from different regions,
                # side by side, which is the defect this whole sequence has
                # been about. The chain is still shown when they differ,
                # because the chain is what supplies the node id.
                print(f"{window_label(r)} {r.catalytic_type:<9} "
                      f"{r.nucleophile or '-':<8} {len(r.classes)} classes")
            else:
                print(f"NOT FOUND — {r.note[:64]}")
            chains.append(r)
            time.sleep(0.4)
        args.cache.parent.mkdir(parents=True, exist_ok=True)
        args.cache.write_text(json.dumps(cache_payload(chains), indent=2))
        print(f"\ncached -> {args.cache}")

    usable = [c for c in chains if c.ok]
    print(f"\n{len(usable)}/{len(chains)} chains resolved")

    # ------------------------------------------------- the analysis window
    # WHAT EACH CHAIN'S CLASSES AND CATALYTIC TYPE WERE READ FROM. Every
    # verdict in this run is a function of this window, and until now it
    # was never printed -- so a chain spanning nine domains and a chain
    # spanning one looked identical in the output.
    #
    # THE SUSPECT CONDITION IS "MORE THAN ONE ANNOTATED DOMAIN IN THE
    # WINDOW", not a size or class-count threshold. The first version of
    # this flagged a protease whose window was twice the median protease
    # window, and it fired on nothing: six of the twelve proteases are the
    # un-narrowed NS3 chains, so they ARE the median and masked themselves.
    # A statistic computed over a population that is mostly outliers does
    # not find the outliers. Domain count is a fact about the chain and
    # cannot be diluted that way.
    print("\nanalysis window — what the classes and catalytic type were "
          "read from:")
    # 'dom' is domains-in-window / domains-in-chain, because for a
    # narrowed chain those differ and only the first drives the flag.
    print(f"  {'chain':<22} {'window':<14} {'aa':>5} {'cls':>4} {'dom':>6}  "
          f"basis")
    for c in usable:
        lo, hi = c.window
        all_d = c.domains or []
        inside = [d for d in all_d if d[0] >= lo and d[1] <= hi]
        flag = ("   <- window holds %d annotated domains" % len(inside)
                if c.spec.role in SCORE_ELIGIBLE_ROLES and len(inside) > 1
                else "")
        print(f"  {c.spec.label:<22} {f'{lo}-{hi}':<14} {hi - lo + 1:>5} "
              f"{len(c.classes):>4} {f'{len(inside)}/{len(all_d)}':>6}  "
              f"{c.window_note or 'whole chain'}{flag}")

    # EVERY score-eligible chain, not only the crowded ones. The previous
    # version printed this block only for windows holding more than one
    # domain, so a single-domain chain never showed its domain's NAME or
    # SPAN -- which is exactly what the polymerase narrowing needs, and why
    # that decision could not be made from the last run: the coronavirus and
    # flavivirus polymerases were listed and the picornavirus 3Ds were not.
    #
    # A chain is also listed when its window is barely covered by any
    # annotated domain. HAV 3C has NO Domain feature at all, and the one
    # report that accounts for active sites dropped it silently for exactly
    # that reason.
    print("\n  domains UniProt annotates in each score-eligible window —")
    print("  the only place a narrowing may be read from:")
    for c in usable:
        if c.spec.role not in SCORE_ELIGIBLE_ROLES:
            continue
        lo, hi = c.window
        inside = [d for d in (c.domains or []) if d[0] >= lo and d[1] <= hi]
        covered = len({p for ds, de, _d in inside
                       for p in range(max(ds, lo), min(de, hi) + 1)})
        pct = 100 * covered / (hi - lo + 1)
        note = ""
        if not inside:
            note = "   <- entry annotates no Domain here"
        elif len(inside) > 1:
            note = f"   <- {len(inside)} enzymes in one window"
        elif pct < 90:
            note = f"   <- domain covers only {pct:.0f}% of the window"
        print(f"    {c.spec.label:<22} [{lo}-{hi}]  {pct:>3.0f}% annotated"
              f"{note}")
        for ds, de, dd in inside:
            print(f"      {ds:>6}-{de:<6} {de - ds + 1:>4} aa  {dd}")

    # WHICH DOMAIN EACH ACTIVE SITE IS IN. A chain's catalytic type is read
    # from its active sites, and over an un-narrowed window those sites can
    # belong to a different enzyme. PLpro's did: N1062 and G1070 are the
    # macrodomain's, 600 residues from the protease.
    #
    # No heuristic here, just the fact, because the last heuristic in this
    # report fired on nothing. A chain labelled polymerase whose every site
    # sits in a methyltransferase is visible without being flagged.
    rows = active_site_rows(usable)
    off = [r for r in rows if r[3].startswith("(")]
    if rows:
        print("\nactive sites by annotated domain — the catalytic type is "
              "read from these:")
        last = None
        for label, role, res, dom in rows:
            head = (f"  {label:<22} {role:<11}" if label != last
                    else " " * 36)
            print(f"{head} {res:<8} {dom}")
            last = label
        if off:
            print(f"\n  {len(off)} site(s) lie in no annotated domain.")

    refused_narrowing = [c for c in chains
                         if c.spec.domain and c.window_is_chain]
    if refused_narrowing:
        print("\n  WINDOW NOT NARROWED THOUGH ONE WAS ASKED FOR — the whole")
        print("  chain was analysed, which is the condition this field exists")
        print("  to prevent. Fix the domain name or drop the request:")
        for c in refused_narrowing:
            print(f"    {c.spec.label:<22} {c.window_note}")

    # A PROTEASE WITH NO CATALYTIC MACHINERY IS A RESOLUTION FAILURE, not a
    # biological fact, and nothing said so.
    #
    # SARS-CoV-2 PLpro was the one this caught, and it is now empty: the
    # widened composition rule reads it cysteine/C1674. Kept because the
    # check is about the class of failure, not that one chain -- and
    # because what it caught was real. This is the polyprotein problem from
    # phase0's header recurring one level down: "Any question about a
    # specific protein must be asked of a specific CHAIN" -- and the chain
    # can itself be multi-domain, which is what ChainSpec.domain addresses.
    # Reported rather than dropped: silently excluding a target is worse.
    miscalled = [c for c in usable
                 if c.spec.role == "protease" and c.catalytic_type == "none"]
    if miscalled:
        print(f"\n  {len(miscalled)} chain(s) declared a protease but "
              f"resolved no catalytic type:")
        for c in miscalled:
            print(f"    {c.spec.label:<22} {c.spec.accession} "
                  f"[{c.start}-{c.end}] {c.end - c.start + 1:,} aa, "
                  f"{len(c.classes)} classes")
        print("    A protease has catalytic machinery by definition, so this "
              "is the chain match\n    being wrong, not the enzyme. Such a "
              "chain cannot be gated on its nucleophile\n    and contributes "
              "nothing to M8.")

    # WHICH RULE DECIDED EACH TYPE. catalytic_machinery has two paths: an
    # explicitly annotated nucleophile, and a composition fallback whose
    # test is "the LAST active site is C or S" -- the chymotrypsin-like
    # assumption the annotated path was already fixed to stop making.
    # Widening the fallback would re-decide every chain that reaches it,
    # and until now nothing said which chains those are.
    by_basis: Counter = Counter(c.catalytic_basis for c in usable)
    print(f"\n  how each catalytic type was decided: "
          + ", ".join(f"{k} {v}" for k, v in by_basis.most_common()))
    fallback = [c for c in usable if c.catalytic_basis == "composition"]
    if fallback:
        print("    decided by composition, so a change to that rule would "
              "move these:")
        for c in fallback:
            print(f"      {c.spec.label:<22} {c.catalytic_type:<9} "
                  f"{', '.join(c.catalytic_residues[:4])}")
    unread = [c for c in usable
              if c.catalytic_basis == "none" and c.catalytic_residues]
    if unread:
        print("    HAVE active sites but no rule could read them:")
        for c in unread:
            print(f"      {c.spec.label:<22} {c.spec.role:<10} "
                  f"{', '.join(c.catalytic_residues[:5])}")
    unkeyed = [c for c in chains
               if c.end > c.start > 0 and not c.chain_pro_id]
    if unkeyed:
        print(f"  {len(unkeyed)} matched a chain but carry no PRO_ id, so "
              f"they have no stable node key:")
        print_capped(unkeyed, 6,
                     lambda c: f"{c.spec.label:<22} {c.note[:60]}",
                     unit="more unkeyed chain(s)")

    # DOES THIS LAYER JOIN THE SPINE? Until now it could not: node_id was
    # built from coordinates while the spine keys the same chain on its PRO_
    # id, so every MEMBER_OF_CLASS edge referenced a node that does not
    # exist. That is almost certainly why `folds` is absent from
    # LAYER_PRECEDENCE -- assemble_release returns 1 on a dangling edge.
    #
    # Reported rather than assumed, because the ids matching in form does
    # not mean the spine built the chain: a virus outside the register, or a
    # chain the spine skipped, still dangles.
    spine_ids: set[str] = set()
    if args.spine and (args.spine / "nodes.jsonl").exists():
        spine_ids = {json.loads(line)["id"]
                     for line in (args.spine / "nodes.jsonl").read_text().splitlines()
                     if line.strip()}
        hit = [c for c in usable if c.node_id in spine_ids]
        miss = [c for c in usable if c.node_id not in spine_ids]
        print(f"\njoin against {args.spine.name}: {len(hit)}/{len(usable)} "
              f"resolved chains exist as Protein nodes")
        # ALL OF THEM. This was miss[:10] while the line beneath it said
        # "14 chain(s) would produce DANGLING edges", so four -- CVB3 3C,
        # Rhinovirus A 3C, Poliovirus 3D, EV-A71 3D -- were counted and not
        # named. In the one report that gates LAYER_PRECEDENCE.
        for c in miss:
            print(f"    MISSING  {c.node_id:<30} {c.spec.label} "
                  f"({c.spec.family}, taxon {c.spec.virus_taxon})")
        if miss:
            print(f"  {len(miss)} chain(s) would produce DANGLING edges. The "
                  f"layer must not be added to\n  LAYER_PRECEDENCE until this "
                  f"is zero, or until those chains are dropped.")

    else:
        print(f"\nno spine release at {args.spine}, so the join is unchecked")

    # --------------------------------------------------- catalytic summary
    print("\ncatalytic type by family and role:")
    grid: Counter = Counter()
    for c in usable:
        grid[(c.spec.family, c.spec.role, c.catalytic_type)] += 1
    for (fam, role, ct), n in sorted(grid.items()):
        flag = "" if role in SCORE_ELIGIBLE_ROLES else "   (ingested, not scored)"
        print(f"  {fam:<16} {role:<14} {ct:<10} {n}{flag}")

    # ----------------------------------------------------- pairwise gate
    print("\npairwise comparison, score-eligible roles only:")
    eligible = [c for c in usable if c.spec.role in SCORE_ELIGIBLE_ROLES]
    results = []
    for i, a in enumerate(eligible):
        for b in eligible[i + 1:]:
            if a.spec.family == b.spec.family:
                continue                     # cross-family only
            g = gate(a, b)
            results.append((a, b, g))

    allowed = [r for r in results if r[2].allowed]
    refused = [r for r in results if not r[2].allowed]
    print(f"  {len(results)} cross-family pairs | {len(allowed)} allowed | "
          f"{len(refused)} refused")

    by_mode: Counter = Counter(g.mode for _a, _b, g in results)
    for mode, n in by_mode.most_common():
        print(f"    {mode:<12} {n}")

    # ALL OF THEM, not the first 14. The cap was [:14] with no "and N more",
    # and sorted strong-first -- so the printed list was 14 strong protease
    # rows and the 22 hidden ones included EVERY weak polymerase pair. The
    # comment below then said all 36 were cysteine protease pairs, and the
    # truncation was what made that look true. Measured: 20 protease pairs
    # and 16 polymerase pairs, and 16 is also the motif-mode total, so every
    # cross-family polymerase pair is allowed. These rows are the layer's
    # actual output -- the FOLD_SIMILAR_TO edges it would write -- and the
    # one list that must not be abridged.
    print(f"\n  allowed transfers, all {len(allowed)} (strong first):")
    for a, b, g in sorted(allowed, key=lambda r: (r[2].strength != "strong",
                                                  a.spec.label, b.spec.label)):
        # EVERY shared class, not the first two. The row cap was fixed one
        # commit ago and this one was left: all 20 protease rows printed
        # "G3DSA:2.40.10.10, IPR009003" and nothing said whether that was
        # the whole evidence or the start of it.
        print(f"    [{g.strength:<6}] {a.spec.label:<20} -> {b.spec.label:<20} "
              f"{', '.join(g.shared)}")

    # WHICH CHAINS HAVE ANY PARTNER AT ALL. M8's entire claim is reaching a
    # virus with no activity data of its own, so a chain with zero allowed
    # transfers is a virus M8 cannot help however the path is wired.
    #
    # Reported per chain because the aggregate hides it. THIS COMMENT USED
    # TO SAY "on v0.1 every one of them is a cysteine protease pair", which
    # the same run refuted one screen below: the partners add to 40 for
    # proteases and 32 for polymerases, so 20 pairs and 16 pairs. It read as
    # true only because the transfer list above was cut at 14 rows, all of
    # them strong protease pairs. Fourth comment outliving its code.
    #
    # The part that IS true and matters: comparison is cross-family only
    # (same-family pairs are skipped above), and the catalytic table shows
    # serine proteases exist in ONE family here -- so a flavivirus protease
    # has no cross-family partner that is not cysteine, and the gate
    # correctly refuses all of
    # them. That is structural, not a data gap.
    partners: Counter = Counter()
    for a, b, g in results:
        if g.allowed:
            partners[a.spec.label] += 1
            partners[b.spec.label] += 1
    stranded = [c for c in eligible if not partners[c.spec.label]]
    print(f"\n  allowed partners per chain "
          f"({len(eligible) - len(stranded)}/{len(eligible)} have at least one):")
    for c in sorted(eligible, key=lambda c: -partners[c.spec.label]):
        mark = "   <- M8 cannot reach this one" if not partners[c.spec.label] else ""
        print(f"    {c.spec.label:<22} {c.catalytic_type:<9} "
              f"{partners[c.spec.label]:>3}{mark}")

    # WHAT THE SHARED CLASSES ACTUALLY DISTINGUISH.
    #
    # A class carried by every pair in a bucket is a constant, and a
    # constant cannot rank anything. If every shared class a pair has is
    # universal, that pair's fold evidence is identical to every other
    # pair's -- the transfer is eligible, but the fold says nothing about
    # WHICH comparator to prefer. That matters directly for M8: DWPC over
    # FOLD_SIMILAR_TO edges carrying identical evidence weights them
    # identically, so the hop would be a filter, not a score.
    #
    # THE REFUSALS ARE CENSUSED TOO, and they are where the gate's
    # justification lives. "Fold alone cannot separate a cysteine protease
    # from a serine one" is a claim about pairs the gate REFUSES, so the
    # class that is universal among those is the one doing the work. The
    # module header named SSF50494 for years; this is what checks it.
    taken, spurned = allowed_pairs(results), refused_with_shared_fold(results)
    print_census("shared classes among ALLOWED transfers, by how many "
                 "pairs carry them:",
                 class_census(taken), class_coverage(taken))
    print_census("shared classes among REFUSED-DESPITE-SHARED-FOLD pairs — "
                 "the gate's own\n  justification, which rests on whatever "
                 "is universal here:",
                 class_census(spurned), class_coverage(spurned))

    # WHAT M8 COULD ACTUALLY WALK.
    #
    # 36 allowed transfers is the fold layer's answer. It is not M8's:
    # the agreed design writes one Protein->Protein FOLD_SIMILAR_TO per
    # evaluated pair, and an edge needs BOTH endpoints in the spine.
    # A transfer whose comparator the spine never built is not a weak
    # hop, it is no hop -- and counting it with the rest is how "36
    # allowed" would become a reach figure nobody could reproduce.
    if spine_ids:
        walkable, unwalkable = walkable_transfers(results, spine_ids)
        print(f"\n  of {len(walkable) + len(unwalkable)} allowed transfers, "
              f"{len(walkable)} have BOTH ends in the spine and could be "
              f"written\n  as Protein->Protein FOLD_SIMILAR_TO edges:")
        for a, b, g in sorted(walkable, key=lambda r: (r[0].spec.label,
                                                       r[1].spec.label)):
            print(f"    [{g.strength:<6}] {a.spec.label:<20} -> "
                  f"{b.spec.label}")
        if not walkable:
            print("    (none)")
        blocked = blocking_chains(unwalkable, spine_ids)
        if blocked:
            print(f"\n  the other {len(unwalkable)} are blocked by comparator "
                  f"chains the spine does not build:")
            for label, n in blocked.most_common():
                print(f"    {label:<22} blocks {n:>2} transfer(s)")
            strengths = Counter(g.strength for _a, _b, g in unwalkable)
            print("    by strength: " + ", ".join(
                f"{k} {v}" for k, v in sorted(strengths.items())))

    # The refusals that matter: same fold, different nucleophile.
    key_refusals = refused_with_shared_fold(results)
    if key_refusals:
        print(f"\n  refused despite a SHARED FOLD ({len(key_refusals)}) — the "
              f"gate doing its job:")
        for a, b, g in sorted(key_refusals,
                              key=lambda r: (r[0].spec.label, r[1].spec.label)):
            print(f"    {a.spec.label:<20} -> {b.spec.label:<20} "
                  f"shares {', '.join(g.shared)} but "
                  f"{a.catalytic_type} != {b.catalytic_type}")

    # ------------------------------------------------------------- write
    em = Emit()
    from kgav.folds import emit_classes
    stats = emit_classes(em, usable, "infores:interpro")
    # NO FOLD_SIMILAR_TO EDGES ARE EMITTED, AND THAT IS DELIBERATE UNTIL THE
    # PATH SHAPE IS SETTLED.
    #
    # This loop used to write `em.edge(cid, "FOLD_SIMILAR_TO", cid, ...)` --
    # subject and object the SAME node, because a shared class is one node, not
    # two. M8 requires two DISTINCT TargetClass nodes joined by that predicate,
    # and walk_metapath discards a neighbour already in `visited`, so every
    # self-loop was skipped and M8 produced zero paths on every run since 0.8.0.
    #
    # The self-loop is not a typo to swap for a pair of ids. It is the data
    # model disagreeing with the metapath: when two chains SHARE a class, the
    # structure that connects them is MEMBER_OF_CLASS in both directions
    # through that one node, and no similarity edge is needed or meaningful.
    # M8's declared path asks for a hop between two classes that this ingest
    # has no principled pair to supply. Inventing one -- every class of A
    # crossed with every class of B -- would mint thousands of similarity
    # assertions nobody measured, which is the failure this project exists to
    # avoid. See the note in the schema on M8.
    #
    # The gate itself is correct and is reported above; what is unresolved is
    # how to write its verdict into the graph. Until that is decided the layer
    # emits classes and membership only, and metapath_reach reports M8 as
    # reaching nothing rather than letting it look scored.
    if any(g.allowed for _a, _b, g in results):
        n_allowed = len(allowed_pairs(results))
        # "n allowed transfers are NOT written" implied that writing them
        # would give M8 n hops. It would give it the walkable subset and
        # n - that dangling edges, which is the self-loop mistake again in
        # a different shape: edges that validate and cannot be walked.
        n_walk = (len(walkable_transfers(results, spine_ids)[0])
                  if spine_ids else None)
        reach = (f"{n_walk} of them have both endpoints in the spine, so "
                 f"writing all\n  {n_allowed} would give M8 {n_walk} hop(s) "
                 f"and {n_allowed - n_walk} dangling edge(s)."
                 if n_walk is not None else
                 "No spine was checked, so how many could be written is "
                 "unknown.")
        print(f"\n  NOTE: none of the {n_allowed} allowed transfers is written "
              f"as a\n  FOLD_SIMILAR_TO edge, so M8 has no similarity hop and "
              f"reports zero reach.\n  {reach}\n  See the comment at this line "
              f"before changing it -- the previous version\n  emitted "
              f"self-loops, which validate cleanly and are silently "
              f"unwalkable.")

    print(f"\n{len(em.nodes):,} TargetClass nodes | {len(em.edges):,} edges")
    for k, v in stats.most_common():
        print(f"  {k}: {v:,}")

    # VALIDATED AGAINST THE SPINE TOO, because this layer emits TargetClass
    # nodes and MEMBER_OF_CLASS edges whose SUBJECTS are Proteins living in
    # the spine. Checking the edges against em.nodes alone reported all 402
    # as DANGLING on every run, including the ten chains the spine really
    # does hold -- a validator that fails identically whether or not the
    # data is right, which is the condition it exists to detect.
    #
    # ingest_vhppi.py solved exactly this and said so: "This release is
    # edges-only: every endpoint lives in the spine or host release ...
    # validation runs here against the real index instead." Ninth instance
    # in this project of a solved problem not reaching its second site.
    known = list(em.nodes.values())
    if args.spine and (args.spine / "nodes.jsonl").exists():
        known += [json.loads(line)
                  for line in (args.spine / "nodes.jsonl").read_text().splitlines()
                  if line.strip()]
    else:
        print(f"\nno spine at {args.spine}: every Protein endpoint will read "
              f"DANGLING, which says\n  nothing about this layer")
    violations = load_schema().validate_batch(known, em.edges)
    if report_violations(violations, note=(
            "\n  A DANGLING subject here is a chain the spine does not "
            "build -- see the join\n  report above. Those are comparator "
            "chains for viruses outside the register,\n  not a defect in "
            "this layer.")):
        return 1
    print("\nschema: PASS")

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in em.nodes.values()))
    (args.out / "edges.jsonl").write_text(
        "\n".join(json.dumps(e) for e in em.edges))
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
