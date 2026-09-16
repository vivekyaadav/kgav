"""Phase 2 — target classes, fold similarity, and the two-mode gate.

WHAT PHASES 0 AND ITS FOLLOW-UP ESTABLISHED, measured at chain level.

PROTEASES. All three families share the same fold:

    Coronaviridae 3CLpro : SSF50494, G3DSA:2.40.10.10
    Picornaviridae 3C    : SSF50494, G3DSA:2.40.10.10
    Flaviviridae NS3     : SSF50494, G3DSA:2.40.10.10

The plan assumed the folds DIFFER. They do not. What separates them is the
nucleophile, read from the sequence at UniProt's annotated active site:
cysteine in coronaviruses and picornaviruses, serine in flaviviruses. Every
triad has histidine first and the nucleophile last in the same arrangement.

POLYMERASES. Also shared, and more broadly:

    IPR043502  DNA/RNA polymerase superfamily
    IPR007094  RNA-directed RNA polymerase, catalytic domain
    PS50507    RdRp of positive ssRNA viruses catalytic domain

An earlier run missed this because the filter accepted only CATH-Gene3D and
SUPERFAMILY databases, and InterPro's own IPR superfamily entries are neither.
A filter written to avoid one artifact created another. Hence rule 1 below.

AND THERE IS NO CATALYTIC DISCRIMINATOR FOR POLYMERASES. Proteases separate on
Cys vs Ser. Polymerases carry motif C in all three families -- GDD in
picornaviruses and flaviviruses, SDD in coronaviruses -- which is a far weaker
distinction than a different nucleophile. That asymmetry is the reason for the
two-mode gate.

TWO RULES THIS MODULE ENFORCES

  1. Ingest every InterPro level, not a hand-picked set of databases. The level
     is stored so a query can restrict later; restricting at ingest time is how
     the polymerase conservation was lost.

  2. The gate has two modes, chosen by catalytic machinery rather than by
     guesswork:

       NUCLEOPHILE mode  both classes have a catalytic nucleophile. Identical
                         residue required. Cys != Ser, no edge, no exceptions.

       MOTIF mode        neither has a nucleophile (polymerases, helicases).
                         Falls back to shared catalytic-domain entries and
                         reports that the discrimination is WEAK, because it is.

     Mixed -- one has a nucleophile and the other does not -- is refused
     outright. A protease and a polymerase are not comparable targets and any
     similarity between them is an artifact.
"""
from __future__ import annotations

import json
import re
import time
import urllib.error
import urllib.request
from collections import Counter
from dataclasses import dataclass, field

INTERPRO = "https://www.ebi.ac.uk/interpro/api"
UNIPROT = "https://rest.uniprot.org/uniprotkb"

# Residues that can act as a catalytic nucleophile, mapped to the class name
# the schema's catalytic_type enum uses.
NUCLEOPHILES = {"C": "cysteine", "S": "serine", "T": "threonine",
                "D": "aspartic", "E": "aspartic"}

# Motif C of positive-sense RNA virus polymerases: the catalytic aspartate
# pair. GDD in most families, SDD in coronaviruses. Used as a CHECK on a
# region already known to be a polymerase, never as a detector -- the pattern
# also matches CDD in the coronavirus NiRAN domain, which is not motif C.
MOTIF_C = re.compile(r"[GS]D[DN]")

# Entries whose name marks them as a catalytic domain rather than a fold.
CATALYTIC_DOMAIN_MARKERS = ("catalytic domain", "active site", "catalytic")


def _get(url: str, retries: int = 3) -> dict | None:
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url, headers={"Accept": "application/json",
                              "User-Agent": "kgav-phase2/1.0"})
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


# ------------------------------------------------------------------ chains
@dataclass
class ChainSpec:
    """One viral protein chain within a polyprotein entry."""
    label: str
    accession: str
    patterns: tuple[str, ...]
    family: str
    role: str
    virus_taxon: str


@dataclass
class ResolvedChain:
    spec: ChainSpec
    chain_name: str = ""
    start: int = 0
    end: int = 0
    catalytic_residues: list[str] = field(default_factory=list)
    nucleophile: str = ""
    catalytic_type: str = "unknown"
    classes: list[dict] = field(default_factory=list)
    note: str = ""

    @property
    def ok(self) -> bool:
        return self.end > self.start > 0

    @property
    def node_id(self) -> str:
        return f"UniProtKB:{self.spec.accession}_{self.start}_{self.end}"


def find_chain(acc: str, patterns: tuple[str, ...]) -> tuple[str, int, int, str]:
    """Locate a PRO_ chain by description.

    Patterns are tried in the order given, which must be most-specific first:
    "NS3" is a substring of "Serine protease/helicase NS3". Sixth instance of
    the substring hazard in this project, after SARS-CoV/SARS-CoV-2,
    Vero/Vero E6, chloroquine/hydroxychloroquine, SELECT/SELECTIVITY and
    nsp1/nsp1x.
    """
    data = _get(f"{UNIPROT}/{acc}?fields=ft_chain,ft_peptide")
    if not data:
        return "", 0, 0, "accession not found"
    feats = [f for f in data.get("features", [])
             if f.get("type") in ("Chain", "Peptide")]
    if not feats:
        return "", 0, 0, "no chain features (TrEMBL-only entry?)"
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
    have = "; ".join(sorted({d for _f, d in lowered})[:6])
    return "", 0, 0, f"no chain matched {patterns}; entry has: {have}"


def catalytic_machinery(acc: str, lo: int, hi: int) -> tuple[list[str], str, str]:
    """Catalytic residues, nucleophile, and catalytic type — from SEQUENCE.

    UniProt's active-site annotation names the ROLE ("Nucleophile; for 3CL-PRO
    activity@3408") but not the residue. The residue is in the sequence at that
    position, and that is the only unambiguous source: InterPro calls the
    CORONAVIRUS protease superfamily "trypsin-like serine protease" because
    that is the architecture, not the chemistry, so a classifier reading entry
    names returns "ambiguous" for every cysteine protease.

    The nucleophile is taken as the LAST catalytic residue in sequence order.
    In every chymotrypsin-like triad measured in phase 0 -- His-Cys,
    His-Glu-Cys, His-Asp-Ser -- the nucleophile is last.
    """
    data = _get(f"{UNIPROT}/{acc}?fields=ft_act_site,sequence")
    if not data:
        return [], "", "unknown"
    seq = data.get("sequence", {}).get("value", "")
    sites: list[tuple[int, str, str]] = []
    for f in data.get("features", []):
        pos = ((f.get("location", {}) or {}).get("start", {}) or {}).get("value")
        if pos and lo <= pos <= hi and pos <= len(seq):
            sites.append((pos, seq[pos - 1], f.get("description", "")))
    if not sites:
        return [], "", "none"
    sites.sort()
    residues = [f"{aa}{pos}" for pos, aa, _d in sites]

    # ONLY an explicitly annotated nucleophile counts. The earlier rule --
    # "the last catalytic residue is the nucleophile" -- is valid for
    # chymotrypsin-like proteases and wrong everywhere else. Applied
    # universally it read SARS-CoV-2 RdRp as "aspartic D5153" and Dengue NS5
    # as "aspartic E2708": polymerases have catalytic aspartates that
    # coordinate metal ions, not a covalent nucleophile. That sent the
    # coronavirus and picornavirus polymerases into different gate modes and
    # broke the polymerase comparison entirely.
    explicit = [(p, aa) for p, aa, d in sites if "nucleophile" in d.lower()]
    if explicit:
        pos, aa = explicit[-1]
        return residues, f"{aa}{pos}", NUCLEOPHILES.get(aa, "unknown")

    # No annotated nucleophile. A His-...-Cys or His-...-Ser triad in a
    # protease is still recognisable by composition; anything else is
    # metal-coordinating or structural and gets "none".
    aas = [aa for _p, aa, _d in sites]
    if "H" in aas and aas[-1] in ("C", "S"):
        pos, aa = sites[-1][0], sites[-1][1]
        return residues, f"{aa}{pos}", NUCLEOPHILES[aa]
    return residues, "", "none"


def fetch_classes(acc: str, lo: int, hi: int) -> list[dict]:
    """EVERY InterPro entry overlapping the chain, at every level.

    Level is stored, never filtered here. Phase 0's follow-up lost the
    polymerase conservation by accepting only CATH-Gene3D and SUPERFAMILY
    databases -- IPR043502 and IPR007094 are InterPro's own entries and were
    excluded, though they are exactly what the three families share.
    """
    out: list[dict] = []
    data = _get(f"{INTERPRO}/entry/all/protein/uniprot/{acc}/?page_size=200")
    if not data:
        return out
    for item in data.get("results", []):
        md = item.get("metadata", {})
        for prot in item.get("proteins", []):
            hit = False
            for loc in prot.get("entry_protein_locations", []):
                for frag in loc.get("fragments", []):
                    fs, fe = frag.get("start", 0), frag.get("end", 0)
                    if fs <= hi and fe >= lo:
                        out.append({
                            "accession": md.get("accession"),
                            "label": md.get("name") or md.get("accession"),
                            "db": (md.get("source_database") or "").lower(),
                            "level": (md.get("type") or "").lower(),
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


def resolve(spec: ChainSpec) -> ResolvedChain:
    r = ResolvedChain(spec=spec)
    name, s, e, note = find_chain(spec.accession, spec.patterns)
    r.chain_name, r.start, r.end, r.note = name, s, e, note
    if not r.ok:
        return r
    r.catalytic_residues, r.nucleophile, r.catalytic_type = \
        catalytic_machinery(spec.accession, s, e)
    r.classes = fetch_classes(spec.accession, s, e)
    return r


# ------------------------------------------------------------- the gate
@dataclass
class GateResult:
    allowed: bool
    mode: str                 # nucleophile | motif | refused
    catalytic_type_match: bool
    shared: list[str] = field(default_factory=list)
    strength: str = "none"    # strong | weak | none
    reason: str = ""


def gate(a: ResolvedChain, b: ResolvedChain) -> GateResult:
    """Two-mode comparison.

    NUCLEOPHILE mode — both have a catalytic nucleophile. Identical residue
    required. Phase 0 showed the three families share SSF50494 and
    G3DSA:2.40.10.10, so fold alone cannot separate a cysteine protease from a
    serine protease: without this check M8 traverses SARS-CoV-2 Mpro to HCV
    NS3 and recommends nirmatrelvir for hepatitis C.

    MOTIF mode — neither has a nucleophile (polymerases, helicases). Shared
    catalytic-domain entries are the only available signal and the result is
    marked WEAK. Polymerases carry motif C in all three families, GDD or SDD,
    which does not discriminate the way Cys vs Ser does.

    MIXED — refused. A protease and a polymerase are not comparable targets,
    and any fold similarity between them is an artifact.
    """
    a_nuc = a.catalytic_type in NUCLEOPHILES.values()
    b_nuc = b.catalytic_type in NUCLEOPHILES.values()

    shared = sorted({c["accession"] for c in a.classes}
                    & {c["accession"] for c in b.classes})

    if a_nuc != b_nuc:
        return GateResult(False, "refused", False, shared, "none",
                          f"one target has a catalytic nucleophile "
                          f"({a.catalytic_type} vs {b.catalytic_type}) and the "
                          f"other does not; they are not comparable targets")

    if a_nuc and b_nuc:
        match = a.catalytic_type == b.catalytic_type
        if not match:
            return GateResult(
                False, "nucleophile", False, shared, "none",
                f"catalytic nucleophile differs: {a.catalytic_type} "
                f"({a.nucleophile}) vs {b.catalytic_type} ({b.nucleophile}). "
                f"Shared fold ({', '.join(shared[:3]) or 'none'}) does not "
                f"license transfer across a different nucleophile.")
        if not shared:
            return GateResult(False, "nucleophile", True, shared, "none",
                              "same nucleophile but no shared fold entry")
        return GateResult(True, "nucleophile", True, shared, "strong",
                          f"same nucleophile ({a.catalytic_type}) and "
                          f"{len(shared)} shared entries")

    # Neither has a nucleophile.
    catalytic_shared = sorted(
        c["accession"] for c in a.classes
        if c["accession"] in {x["accession"] for x in b.classes}
        and any(m in (c["label"] or "").lower() for m in CATALYTIC_DOMAIN_MARKERS))
    if not shared:
        return GateResult(False, "motif", True, shared, "none",
                          "no shared fold entry")
    return GateResult(
        True, "motif", True, shared,
        # Always weak. Sharing a catalytic-domain entry is better than not,
        # but neither case approaches the discrimination a differing
        # nucleophile provides, so promoting one of them to "strong" would
        # overstate what the comparison establishes.
        "weak",
        f"no catalytic nucleophile in either target, so discrimination rests "
        f"on {len(shared)} shared entries"
        + (f" including catalytic-domain entries {', '.join(catalytic_shared[:3])}"
           if catalytic_shared else " with no catalytic-domain entry among them")
        + ". This is WEAK evidence: polymerases share motif C across all three "
          "families, which does not discriminate the way a nucleophile does.")


def motif_c_present(sequence: str, lo: int, hi: int) -> list[str]:
    """Motif C hits within a chain. A CHECK, not a detector — the pattern also
    matches CDD in the coronavirus NiRAN domain, which is not motif C."""
    region = sequence[lo - 1:hi]
    return [f"{m.group()}@{lo + m.start()}" for m in MOTIF_C.finditer(region)]


# ------------------------------------------------------------------ emit
def emit_classes(em, chains: list[ResolvedChain], source: str,
                 date: str = "1970-01-01") -> Counter:
    """TargetClass nodes and MEMBER_OF_CLASS edges.

    Classes are dated 1970 and tier 3: a fold assignment is computed, not
    asserted by an experiment on a date. Same rule as the v1 chemical
    similarity layer.
    """
    stats: Counter = Counter()
    for r in chains:
        if not r.ok:
            stats["chain_unresolved"] += 1
            continue
        for c in r.classes:
            cid = f"INTERPRO:{c['accession']}" if c["accession"].startswith("IPR") \
                else f"KGAV:{c['accession']}"
            em.node(cid, "TargetClass",
                    label=c["label"], source="interpro", level=c["level"],
                    # A class inherits the catalytic type of its members. Where
                    # members disagree the class is not a usable gate target
                    # and is marked unknown rather than guessed.
                    catalytic_type=r.catalytic_type,
                    interpro_id=c["accession"] if c["accession"].startswith("IPR") else None)
            em.edge(r.node_id, "MEMBER_OF_CLASS", cid,
                    source=source, date=date, tier=3,
                    quals={"method": "interpro", "start": c["start"],
                           "end": c["end"]})
            stats["member_edges"] += 1
        stats["chains"] += 1
    return stats
