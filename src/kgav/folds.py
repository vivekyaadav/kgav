"""Phase 2 — target classes, fold similarity, and the two-mode gate.

WHAT PHASES 0 AND ITS FOLLOW-UP ESTABLISHED, measured at chain level.

PROTEASES. All three families share a fold, and the plan's assumption that
they differ was wrong. What separates them is the nucleophile, read from the
sequence at UniProt's annotated active site: cysteine in coronaviruses and
picornaviruses, serine in flaviviruses.

THE TABLE THAT USED TO BE HERE WAS WRONG, AND IT WAS THE PROJECT'S
FOUNDATIONAL MEASUREMENT. It read:

    Coronaviridae 3CLpro : SSF50494, G3DSA:2.40.10.10
    Picornaviridae 3C    : SSF50494, G3DSA:2.40.10.10
    Flaviviridae NS3     : SSF50494, G3DSA:2.40.10.10

The 2026-10-15 run, printing every shared class instead of the first two,
shows SSF50494 on MERS-CoV and HCoV-229E 3CLpro and on all five 3Cs and all
six NS3pros -- but NOT on SARS-CoV-2 or SARS-CoV 3CLpro, the two chains this
project exists to serve. Deducible without a new fetch, because a shared set
is an intersection: MERS and 229E each share SSF50494 with EV-A71 3C and the
two SARS 3CLpros do not, so EV-A71 3C has it and those two do not. The
census counted it at exactly 10/20 allowed pairs, which is MERS's five plus
229E's five.

G3DSA:2.40.10.10 is likewise absent from four of the six NS3pros.

WHAT IS ACTUALLY UNIVERSAL is IPR009003, Peptidase S1 PA clan: every
cross-family protease pair carries it, allowed and refused alike. That, not
SSF50494, is the class that makes fold alone unable to separate a cysteine
protease from a serine one -- which is the danger the gate exists for, so
the gate's justification is intact and only its cited evidence was wrong.

Triad ORDER is not universal either. "Every triad has histidine first and
the nucleophile last" was in this header until SARS-CoV-2 PLpro was read
correctly: C1674, H1835, D1849, a papain-like triad with the nucleophile
FIRST. The composition rule was widened for exactly that chain, and this
sentence survived the commit that disproved it.

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
# The same pattern ingest_spine._chain_id applies, for the same reason:
# a chain without a PRO_ id cannot be given a stable node key.
PRO_RE = re.compile(r"^PRO_\d+$")

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
    #: UniProt Domain description to narrow the ANALYSIS WINDOW to, within
    #: the matched chain. Empty means the whole chain, which is right for a
    #: chain that is one enzyme and wrong for a chain that is several.
    #:
    #: THE CHAIN STAYS THE IDENTITY. Only the window moves: the node is
    #: still keyed on the chain's PRO_ id, because that is what the spine
    #: has. What narrows is the region the class fetch and the catalytic
    #: rule look at -- the two things a neighbouring domain corrupts.
    domain: str = ""


#: Fields of ResolvedChain that are DERIVED by `decide_machinery` from
#: `active_sites`, rather than fetched. They are never written to the chain
#: cache and always recomputed on load.
#:
#: WHY: `--offline` rebuilt ResolvedChain straight from JSON, so
#: `catalytic_machinery` ran only on the network path. Widening the
#: composition rule therefore changed nothing offline -- the run reported
#: "composition 14, none 9, annotated 1" identically before and after,
#: because it was replaying the verdicts of the run before it. A cached
#: verdict looks exactly like a computed one, which is the same failure the
#: cache-completeness guard below already names for defaulted fields, one
#: level up: that guard catches a field that is ABSENT, not one that is
#: present and stale.
DERIVED_FIELDS = frozenset({"catalytic_residues", "nucleophile",
                            "catalytic_type", "catalytic_basis"})


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
    #: Which rule decided catalytic_type: "annotated", "composition" or
    #: "none". A type nobody annotated and one nobody could infer are
    #: different situations and the value alone cannot tell them apart.
    catalytic_basis: str = "none"
    #: The chain's UniProt featureId, e.g. PRO_0000449623. THE JOIN to the
    #: spine, which keys the same protein on it. Empty means unresolvable.
    chain_pro_id: str = ""
    #: THE INPUT the catalytic rule reads: [(position, residue, description)],
    #: straight off UniProt. Stored so the four DERIVED_FIELDS above can be
    #: recomputed instead of replayed -- see DERIVED_FIELDS.
    #:
    #: None and [] are different and must stay different: None means the
    #: fetch never happened or failed (catalytic_type "unknown"), [] means it
    #: happened and the chain has no annotated active site ("none").
    active_sites: list | None = None
    #: The region the class fetch and the catalytic rule actually read.
    #: 0 means "not narrowed", and `window` falls back to the whole chain.
    window_start: int = 0
    window_end: int = 0
    #: Why the window is what it is. Non-empty and window_start == 0 means a
    #: narrowing was ASKED FOR and did not happen -- which must be visible,
    #: because silently analysing the whole chain is the bug being fixed.
    window_note: str = ""
    #: Every Domain UniProt annotates inside the chain, as
    #: [(start, end, description)]. Carried so the report can show what a
    #: chain COULD be narrowed to without anyone naming domains from memory.
    domains: list | None = None

    @property
    def window(self) -> tuple[int, int]:
        """The coordinates analysed, which are not always the chain's.

        SARS-CoV-2 PLpro's chain is nsp3, 819-2763: nine domains, 1,945
        residues, 60 InterPro classes against 7-12 for a protease that is
        only a protease. Its active-site list carried N1062 and G1070 --
        the MACRODOMAIN's ADP-ribosylhydrolase sites, 600 residues from the
        protease -- and the catalytic rule read cysteine only because
        neither of those happens to be a C or an S.
        """
        return (self.window_start or self.start,
                self.window_end or self.end)

    @property
    def window_is_chain(self) -> bool:
        return not self.window_start

    def redecide(self) -> None:
        """Re-run the rule over `active_sites`, discarding any stored verdict.

        Called on every cache load. A cached verdict is indistinguishable
        from a freshly computed one, so the only safe thing to do with one
        is not keep it.
        """
        (self.catalytic_residues, self.nucleophile,
         self.catalytic_type, self.catalytic_basis) = \
            decide_machinery(self.active_sites)

    @property
    def ok(self) -> bool:
        # The PRO_ id is part of being resolved, not an extra. A chain with
        # coordinates but no id produced a node_id that joined to nothing,
        # which is the failure this field exists to end -- and `ok` is what
        # every caller filters on.
        return self.end > self.start > 0 and bool(self.chain_pro_id)

    @property
    def node_id(self) -> str:
        """The spine's id for this chain, not a coordinate-derived one."""
        if not self.chain_pro_id:
            raise ValueError(
                f"{self.spec.label}: no PRO_ id, so no stable node key. "
                "Callers must filter on .ok before reading node_id.")
        return f"UniProtKB:{self.chain_pro_id}"


def find_chain(acc: str,
               patterns: tuple[str, ...]) -> tuple[str, int, int, str, str]:
    """Locate a PRO_ chain by description. Returns its featureId too.

    Patterns are tried in the order given, which must be most-specific first:
    "NS3" is a substring of "Serine protease/helicase NS3". Sixth instance of
    the substring hazard in this project, after SARS-CoV/SARS-CoV-2,
    Vero/Vero E6, chloroquine/hydroxychloroquine, SELECT/SELECTIVITY and
    nsp1/nsp1x.

    THE featureId IS THE JOIN, and dropping it is why this layer has never
    reached the graph. ResolvedChain.node_id was built from coordinates --
    UniProtKB:P0DTD1_3264_3569 -- while ingest_spine keys the same chain on
    its PRO_ id, UniProtKB:PRO_0000449623 (see its _chain_id, which reads
    exactly this field). Two id schemes for one protein, so every
    MEMBER_OF_CLASS edge pointed at a node that does not exist, assembling
    the layer would have produced dangling edges, and assemble_release
    returns 1 on dangling. The response already carried the key:

        {'type': 'Chain', 'featureId': 'PRO_0000449618',
         'description': 'Replicase polyprotein 1ab'}

    verified against rest.uniprot.org for P0DTD1 on 2026-10-15 with this
    exact field selection, because a field present in a proteome download is
    not evidence about a single-entry query.
    """
    data = _get(f"{UNIPROT}/{acc}?fields=ft_chain,ft_peptide")
    if not data:
        return "", 0, 0, "accession not found", ""
    feats = [f for f in data.get("features", [])
             if f.get("type") in ("Chain", "Peptide")]
    if not feats:
        return "", 0, 0, "no chain features (TrEMBL-only entry?)", ""
    lowered = [(f, (f.get("description") or "").lower()) for f in feats]
    for pat in patterns:
        p = pat.lower()
        for f, desc in lowered:
            if p in desc:
                loc = f.get("location", {})
                s = (loc.get("start") or {}).get("value")
                e = (loc.get("end") or {}).get("value")
                if s and e:
                    # Same rule as the spine: a chain with no PRO_ id has no
                    # stable key, so it is reported rather than given a
                    # coordinate-based one that joins to nothing.
                    fid = f.get("featureId") or ""
                    if not PRO_RE.match(fid):
                        return (f.get("description", pat), int(s), int(e),
                                f"chain matched but carries no PRO_ id "
                                f"(featureId {fid!r}); no stable node key", "")
                    return f.get("description", pat), int(s), int(e), "", fid
    have = "; ".join(sorted({d for _f, d in lowered})[:6])
    return "", 0, 0, f"no chain matched {patterns}; entry has: {have}", ""


def fetch_domains(acc: str, lo: int,
                  hi: int) -> list[tuple[int, str]] | None:
    """Every UniProt Domain inside [lo, hi] as (start, end, description).

    Fetched for EVERY chain, not just the ones being narrowed, so the report
    can show what each chain could be narrowed to. The alternative is
    someone naming domains from memory, and the boundary is the whole point.

    None if the fetch failed.
    """
    data = _get(f"{UNIPROT}/{acc}?fields=ft_domain")
    if not data:
        return None
    out = []
    for f in data.get("features", []):
        if f.get("type") != "Domain":
            continue
        loc = f.get("location", {}) or {}
        s = (loc.get("start", {}) or {}).get("value")
        e = (loc.get("end", {}) or {}).get("value")
        if s and e and s >= lo and e <= hi:
            out.append((int(s), int(e), f.get("description", "")))
    return sorted(out)


def pick_domain(domains: list | None,
                description: str) -> tuple[int, int, str]:
    """The named domain's span, or (0, 0, why not).

    NO SILENT FALLBACK. A requested narrowing that cannot be performed
    returns a reason, and the caller must surface it: quietly analysing the
    whole chain instead is exactly the behaviour this is here to end.
    """
    if domains is None:
        return 0, 0, "domain fetch failed; window not narrowed"
    want = description.strip().lower()
    hits = [(s, e, d) for s, e, d in domains
            if d.strip().lower() == want]
    if not hits:
        # NEAR MISSES FIRST. Sorted alphabetically and cut at 8, this list
        # dropped "Peptidase C16" from the message for a request of
        # "Peptidase C30" -- hiding the one name that would have fixed the
        # typo. The point of the message is to be actionable.
        words = set(want.split())
        names = sorted({d for _s, _e, d in domains},
                       key=lambda d: (-len(words & set(d.lower().split())), d))
        have = ", ".join(names[:8]) or "none"
        if len(names) > 8:
            have += f", +{len(names) - 8} more"
        return 0, 0, f"no Domain matched {description!r}; entry has: {have}"
    if len(hits) > 1:
        return 0, 0, (f"{len(hits)} Domains match {description!r}; "
                      f"ambiguous, so the window was not narrowed")
    s, e, _d = hits[0]
    return s, e, ""


def fetch_sites(acc: str, lo: int,
                hi: int) -> list[tuple[int, str, str]] | None:
    """NETWORK HALF. The active sites in [lo, hi] as (position, residue,
    description).

    UniProt's active-site annotation names the ROLE ("Nucleophile; for
    3CL-PRO activity@3408") but not the residue. The residue is in the
    sequence at that position, and that is the only unambiguous source:
    InterPro calls the CORONAVIRUS protease superfamily "trypsin-like serine
    protease" because that is the architecture, not the chemistry, so a
    classifier reading entry names returns "ambiguous" for every cysteine
    protease. So the sequence lookup happens here, and what comes out is a
    payload the rule can be re-run over without the network.

    Returns None if the fetch failed, [] if it succeeded and found nothing.
    """
    data = _get(f"{UNIPROT}/{acc}?fields=ft_act_site,sequence")
    if not data:
        return None
    seq = data.get("sequence", {}).get("value", "")
    sites: list[tuple[int, str, str]] = []
    for f in data.get("features", []):
        pos = ((f.get("location", {}) or {}).get("start", {}) or {}).get("value")
        if pos and lo <= pos <= hi and pos <= len(seq):
            sites.append((pos, seq[pos - 1], f.get("description", "")))
    sites.sort()
    return sites


def decide_machinery(
        sites: list | None) -> tuple[list[str], str, str, str]:
    """PURE HALF — the rule. No network, so it can be re-run over a cache.

    THE SPLIT EXISTS BECAUSE this function's verdict used to be what got
    cached. `--offline` then replayed decisions made by whatever version of
    this rule happened to be installed when the cache was written, and a
    change here had no effect until someone ran against the network. The
    cache now holds `sites`, this runs on every load, and the rule in force
    is always the rule in the file.

    `sites` is [(position, residue, description)]; None means the fetch
    failed. The two produce different types -- "unknown" against "none" --
    and conflating them would report a network outage as a chain with no
    catalytic machinery.

    THIS PARAGRAPH USED TO SAY "the nucleophile is taken as the LAST
    catalytic residue in sequence order", and the code below stopped doing
    that -- see the comment on `explicit`, which records reading SARS-CoV-2
    RdRp as "aspartic D5153" because of it. The rule was removed from the
    code and left in the docstring, where it described behaviour that no
    longer existed. Second time in this project after selectivity.py's
    "WORST is kept"; a comment that outlives its code is a claim nobody can
    check against anything.

    What the code does: an explicitly annotated nucleophile decides, and
    nothing else does. Failing that, composition is consulted -- see the
    fallback, whose own assumption is narrower than it looks.

    Returns (residues, nucleophile, catalytic_type, basis), where basis is
    "annotated", "composition" or "none" and says WHICH of those decided.
    Without it the per-chain report cannot distinguish a type nobody
    annotated from one nobody could infer.
    """
    if sites is None:
        return [], "", "unknown", "none"
    if not sites:
        return [], "", "none", "none"
    # JSON round-trips a tuple as a list, so a cache load hands this lists.
    sites = sorted(tuple(s) for s in sites)
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
        return (residues, f"{aa}{pos}",
                NUCLEOPHILES.get(aa, "unknown"), "annotated")

    # No annotated nucleophile. A His-...-Cys or His-...-Ser triad in a
    # protease is still recognisable by composition; anything else is
    # metal-coordinating or structural and gets "none".
    aas = [aa for _p, aa, _d in sites]
    # POSITION-INDEPENDENT NOW. The test was "the LAST residue is C or S"
    # -- the chymotrypsin-like assumption the explicit branch above was
    # fixed to stop making, still in force down here. A papain-like
    # protease runs Cys-His-Asp with the nucleophile FIRST, so SARS-CoV-2
    # PLpro (C1674, H1835, D1849, none annotated) read "none" and the gate
    # could not evaluate a cysteine protease.
    #
    # THE RULE: a histidine must be present, and exactly ONE site may be C
    # or S. That one is the nucleophile wherever it sits. Both conditions
    # carry weight, and the first is what stops this re-introducing the bug
    # the explicit branch fixed:
    #
    #   the histidine   SARS-CoV-2 RdRp's sites are S5151, D5152, D5153 --
    #                   a serine, but no His, so this refuses it. Without
    #                   that test it would read "serine" for a polymerase,
    #                   one step from the "aspartic D5153" error.
    #   exactly one     two candidates mean the set spans more than one
    #                   active site and nothing here can say which belongs
    #                   to the protease. Refusing is sayable; guessing is
    #                   not.
    #
    # MEASURED BEFORE CHANGING, against all 19 chains that have sites. The
    # 14 this branch already decided come out bit-identical -- every one of
    # them does carry its nucleophile last -- the four polymerases stay
    # "none", and exactly one chain moves: PLpro, "none" -> cysteine.
    cands = [(p, aa) for p, aa, _d in sites if aa in ("C", "S")]
    if "H" in aas and len(cands) == 1:
        pos, aa = cands[0]
        return residues, f"{aa}{pos}", NUCLEOPHILES[aa], "composition"
    return residues, "", "none", "none"


def catalytic_machinery(acc: str, lo: int,
                        hi: int) -> tuple[list[str], str, str, str]:
    """Fetch then decide. The whole operation, for callers with a network.

    Kept so there is still one name for "work out this chain's catalytic
    machinery", but it is now composition, not implementation: the rule
    lives in `decide_machinery` alone and nothing may copy it.
    """
    return decide_machinery(fetch_sites(acc, lo, hi))


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
    name, s, e, note, pro = find_chain(spec.accession, spec.patterns)
    r.chain_name, r.start, r.end, r.note, r.chain_pro_id = name, s, e, note, pro
    if not r.ok:
        return r
    # The window is settled BEFORE anything is read, because everything
    # downstream -- active sites, InterPro classes, and therefore the
    # catalytic type and every gate verdict -- is a function of it.
    r.domains = fetch_domains(spec.accession, s, e)
    if spec.domain:
        ws, we, why = pick_domain(r.domains, spec.domain)
        r.window_start, r.window_end, r.window_note = ws, we, why
        if not why:
            r.window_note = f"narrowed to {spec.domain} ({ws}-{we})"

    lo, hi = r.window
    # Fetch the sites, store them, then decide from the stored copy -- the
    # same call the cache path makes, so offline and online cannot diverge.
    r.active_sites = fetch_sites(spec.accession, lo, hi)
    r.redecide()
    r.classes = fetch_classes(spec.accession, lo, hi)
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
    required, because fold alone cannot separate a cysteine protease from a
    serine one: without this check M8 traverses SARS-CoV-2 Mpro to HCV NS3
    and recommends nirmatrelvir for hepatitis C.

    THE CLASS THAT MAKES THAT TRUE IS IPR009003, Peptidase S1 PA clan,
    carried by every cross-family protease pair. This docstring used to
    cite SSF50494 and G3DSA:2.40.10.10 instead; the run shows SARS-CoV-2
    3CLpro carries neither SSF50494 nor, against four of the six NS3pros,
    a shared G3DSA entry. The Mpro-to-NS3 traversal is still real and
    still refused -- the two do share IPR009003 and G3DSA:2.40.10.10 --
    but the evidence named here was not the evidence doing the work.

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

    # ROLE FIRST, because the paragraph above claims this rule and the code
    # did not implement it. MIXED was decided purely by whether a nucleophile
    # was FOUND, which is a property of the sequence scan, not of the target.
    # SARS-CoV-2 PLpro is the live case: its spec matches UniProt's
    # "Papain-like protease nsp3", a 1,945-residue multi-domain chain
    # (819-2763, 60 InterPro classes against 7-26 for every other chain), so
    # no single nucleophile resolves and catalytic_type reads "none". A
    # protease then entered MOTIF mode alongside the polymerases and could
    # have produced an allowed weak transfer between a protease and a
    # polymerase -- the one comparison this docstring says is an artifact.
    #
    # On v0.1 it did not fire, because PLpro shares no catalytic-domain entry
    # with any cross-family polymerase and scored 0 allowed partners. That is
    # the data, not the gate. A gate that holds by luck is the kind this
    # project keeps finding.
    if a.spec.role != b.spec.role:
        return GateResult(False, "refused", False, shared, "none",
                          f"different target roles ({a.spec.role} vs "
                          f"{b.spec.role}); not comparable targets, and any "
                          f"fold similarity between them is an artifact "
                          f"whatever the catalytic type says")

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

    def _cid(c: dict) -> str:
        return (f"INTERPRO:{c['accession']}" if c["accession"].startswith("IPR")
                else f"KGAV:{c['accession']}")

    # RESOLVED BEFORE ANYTHING IS EMITTED, because Emit.node merges with
    # `properties.update(props)` -- last write wins, silently, recording only
    # a node_merged counter. So the comment below USED TO BE FALSE: a class
    # holding a cysteine protease and a serine protease took whichever chain
    # happened to be written last, and reversing the chain order flipped it
    # from "serine" to "cysteine".
    #
    # That is the fifth insertion-order dependence found in this project, and
    # the worst placed: catalytic_type is the field M8's gate reads. The
    # schema warns that without that gate M8 "traverses from SARS-CoV-2 Mpro
    # to HCV NS3 -- same superfamily, opposite nucleophile -- and recommends
    # nirmatrelvir for hepatitis C", and those two land in the same InterPro
    # superfamily. A gate whose value depends on emission order is not a gate.
    #
    # "unknown" is a declared value of the catalytic_type enum, so a mixed
    # class is sayable rather than guessable.
    seen: dict[str, set[str]] = {}
    for r in chains:
        if r.ok:
            for c in r.classes:
                seen.setdefault(_cid(c), set()).add(r.catalytic_type)
    resolved = {cid: (next(iter(ts)) if len(ts) == 1 else "unknown")
                for cid, ts in seen.items()}
    stats["classes_mixed_catalytic_type"] = sum(
        1 for ts in seen.values() if len(ts) > 1)

    for r in chains:
        if not r.ok:
            stats["chain_unresolved"] += 1
            continue
        for c in r.classes:
            cid = _cid(c)
            em.node(cid, "TargetClass",
                    label=c["label"], source="interpro", level=c["level"],
                    # A class inherits the catalytic type of its members.
                    # Where members disagree the class is not a usable gate
                    # target and is marked unknown rather than guessed. This
                    # now reflects the code, which it did not before.
                    catalytic_type=resolved[cid],
                    interpro_id=c["accession"] if c["accession"].startswith("IPR") else None)
            em.edge(r.node_id, "MEMBER_OF_CLASS", cid,
                    source=source, date=date, tier=3,
                    quals={"method": "interpro", "start": c["start"],
                           "end": c["end"]})
            stats["member_edges"] += 1
        stats["chains"] += 1
    return stats
