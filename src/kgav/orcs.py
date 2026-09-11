"""Day 6: host dependency and restriction factors from BioGRID ORCS.

The densest repurposing signal in the graph, and the only layer where getting a
sign backwards produces a confidently wrong recommendation rather than a missing
one.

WHY DIRECTION IS THE WHOLE PROBLEM

    dependency factor knocked out -> virus cannot replicate -> cell survives
    restriction factor knocked out -> virus replicates better -> cell dies

Both register as "hit". Topologically identical, biologically opposite. Inhibit
a restriction factor with a drug and you HELP the virus. Nothing downstream
catches the error: the graph is structurally perfect either way, the metapaths
resolve, the validator passes.

HOW DIRECTION IS DETERMINED, in order of confidence:

  1. SCREEN_RATIONALE is one-sided ("Increased resistance to virus" ->
     dependency; "Decreased resistance to virus" -> restriction). Screen-level,
     unambiguous, 41 of 96 coronavirus screens.
  2. SCREEN_RATIONALE reports both directions AND Score.1 is a signed effect
     (Z-score, Log2FC, Beta). Direction is sign(Score.1) per gene. 35 screens.
     Verified against screen 1379 (Hoffmann 2020, Huh-7.5, HCoV-OC43): positive
     Z gives RAB7A, VPS11, VPS39, ATP6AP1, ATP6V1A -- late endosome, HOPS, and
     v-ATPase, i.e. the acidification-dependent entry route OC43 requires.
     Negative Z gives TBK1, the central interferon-induction kinase. So
     positive = dependency, negative = restriction.
  3. Neither applies (one-sided p-value such as MAGeCK pos score) -> the screen
     needs manual curation from its paper and is SKIPPED unless an override is
     supplied.

  Screens that are not virus-resistance screens at all -- cell-essential gene
  screens, spike-binding screens, frameshifting reporters -- are excluded
  entirely. A cell-essential gene is not a host factor; ingesting it would put
  ribosomal proteins on every path.

SIGN CALIBRATION, PER SCREEN

There is no global sign convention. On the real data 21 screens score positive
= dependency and 11 score positive = restriction. A single global rule would
have inverted the biology for a third of them, silently.

So each screen is asked which convention it uses, via genes whose polarity is
not in doubt: the type-I interferon axis (restriction only) and established
entry/trafficking factors (dependency only). Both sets are needed, because the
interferon axis only registers as a hit in cells with an intact interferon
response and Huh-7.5 -- common in these screens -- has none.

Calibration runs for EVERY screen, including one-sided ones: a one-sided
SCREEN_RATIONALE describes the screen's primary selection, not the polarity of
every hit in it. A screen whose anchors contradict its rationale is rejected.

KNOWN LIMITATION. Screen-level calibration cannot separate the two tails WITHIN
a screen. Genes in the interferon pathway remain contested across screens
(IFNAR1 2:3, IRF9 2:2). Direction for the dependency side is well supported by
independent controls in neither anchor set -- LY6E 1:7 restriction, NPC1 8:0,
TMEM106B 7:0, PIK3C3 4:0, SCAP 6:1, SREBF2 3:0 dependency -- but
RESTRICTION-FACTOR FINDINGS ARE NOT SUPPORTED IN v1. Resolving this needs
per-screen score semantics read from each paper's methods.
"""
from __future__ import annotations

import re
from collections import Counter, defaultdict
from pathlib import Path

# Type-I interferon induction and signalling. A knockout of any of these makes
# viral infection WORSE; none can legitimately be a dependency factor for a
# coronavirus in a survival screen.
RESTRICTION_SENTINELS = {
    "TBK1", "IFNAR1", "IFNAR2", "STAT1", "STAT2", "IRF3", "IRF9",
    "MAVS", "JAK1", "TYK2", "IFIH1", "DDX58", "RNASEL", "OAS1",
}

DEPENDENCY_ANCHORS = {
    "ACE2", "TMPRSS2", "CTSL", "ANPEP", "DPP4", "NPC1",
    "RAB7A", "RAB2A", "VPS11", "VPS16", "VPS18", "VPS33A", "VPS39", "VPS41",
    "CCDC22", "COMMD3", "VPS35", "VPS26A", "VPS29",
    "ATP6AP1", "ATP6AP2", "ATP6V1A", "ATP6V0C", "ATP6V1B2",
    "TMEM106B", "EXT1", "EXTL3", "B3GALT6", "SLC35B2", "B4GALT7",
}

SIGNED_SCORES = {"Z-score", "Log2FC", "Beta Score", "LFC", "log2 fold change"}

RATIONALE_DEPENDENCY = "increased resistance to virus"
RATIONALE_RESTRICTION = "decreased resistance to virus"
RATIONALE_BOTH = "increased/decreased resistance to virus"

METHODOLOGY_TO_SCREEN_TYPE = {
    "knockout": "CRISPRko",
    "activation": "CRISPRa",
    "inhibition": "CRISPRi",
    "interference": "CRISPRi",
    "knockdown": "RNAi",
}

VIRUS_PATTERNS = [
    (re.compile(r"sars-?cov-?2|severe acute respiratory syndrome coronavirus 2", re.IGNORECASE), "2697049"),
    (re.compile(r"\bmers\b|middle east respiratory", re.IGNORECASE), "1335626"),
    (re.compile(r"sars-?cov\b(?!-?2)|sars coronavirus", re.IGNORECASE), "694009"),
    (re.compile(r"229e", re.IGNORECASE), "11137"),
    (re.compile(r"nl63", re.IGNORECASE), "277944"),
    (re.compile(r"oc43", re.IGNORECASE), "31631"),
    (re.compile(r"hku1", re.IGNORECASE), "290028"),
]


def virus_taxon(screen: dict) -> str | None:
    """Map a screen's condition text to one of the seven species taxa.

    SARS-CoV is matched only after SARS-CoV-2 has been ruled out: the substring
    'SARS-CoV' occurs inside 'SARS-CoV-2', and pattern order is the only thing
    keeping every SARS-CoV-2 screen from being filed under SARS-CoV.
    """
    text = " ".join(str(screen.get(f, "")) for f in
                    ("CONDITION_NAME", "PHENOTYPE", "NOTES", "SCREEN_NAME"))
    for pattern, taxon in VIRUS_PATTERNS:
        if pattern.search(text):
            return taxon
    return None


def screen_type(screen: dict) -> str:
    m = str(screen.get("METHODOLOGY", "")).lower()
    for key, val in METHODOLOGY_TO_SCREEN_TYPE.items():
        if key in m:
            return val
    return "other"


def polarity_mode(screen: dict) -> tuple[str, str | None]:
    """Return (mode, fixed_direction).

    mode is 'screen_level' (direction fixed for every hit), 'signed' (direction
    from sign(Score.1) per gene), or 'unresolvable'.
    """
    rationale = str(screen.get("SCREEN_RATIONALE", "")).strip().lower()
    if rationale == RATIONALE_DEPENDENCY:
        return "screen_level", "dependency"
    if rationale == RATIONALE_RESTRICTION:
        return "screen_level", "restriction"
    if rationale == RATIONALE_BOTH:
        if str(screen.get("SCORE.1_TYPE", "")).strip() in SIGNED_SCORES:
            return "signed", None
        return "unresolvable", None
    return "unresolvable", None


def is_virus_resistance_screen(screen: dict) -> bool:
    """Exclude cell-essential, spike-binding and reporter screens.

    A cell-essential gene is not a host factor -- it is a gene the cell needs
    with or without the virus. Ingesting those would put ribosomal proteins
    into every path.
    """
    rationale = str(screen.get("SCREEN_RATIONALE", "")).strip().lower()
    return "resistance to virus" in rationale


def gene_direction(score1: str | float, fixed: str | None,
                   multiplier: int = 1) -> str | None:
    """multiplier comes from calibrate_sign: +1 if positive means dependency in
    this screen, -1 if positive means restriction."""
    if fixed:
        return fixed
    try:
        v = float(score1) * multiplier
    except (TypeError, ValueError):
        return None
    if v == 0:
        return None
    return "dependency" if v > 0 else "restriction"


def calibrate_sign(hits: list[dict]) -> tuple[int, list[str]]:
    """Infer a screen's sign convention from its own anchor-gene hits.

    A global rule does not work. Screen 1379 (Hoffmann 2020) scores positive Z
    for RAB7A/VPS11/ATP6AP1 -- dependency. Screen 2361 (Le Pen 2024) states
    outright that z-score >= 2 means antiviral -- positive is RESTRICTION. On
    the real data 21 screens use one convention and 11 the other.

    So each screen is asked which convention it uses, via genes whose polarity
    is not in doubt: the type-I interferon axis (restriction only) and the
    established entry/trafficking factors (dependency only). Both sets are
    needed -- the interferon axis only registers as a hit in cells with an
    intact interferon response, and Huh-7.5, used in many screens, has none.

    Returns (multiplier, anchors_seen). +1 when positive means dependency, -1
    when positive means restriction, 0 when the screen cannot be calibrated.
    An uncalibrated screen is skipped, never given a guessed convention.
    """
    votes: list[int] = []
    seen: list[str] = []
    for g in hits:
        sym = (g.get("OFFICIAL_SYMBOL") or "").upper()
        if sym in RESTRICTION_SENTINELS:
            anchor = -1
        elif sym in DEPENDENCY_ANCHORS:
            anchor = 1
        else:
            continue
        try:
            v = float(g.get("SCORE.1"))
        except (TypeError, ValueError):
            continue
        if v == 0:
            continue
        votes.append(anchor if v > 0 else -anchor)
        seen.append(sym)
    if not votes:
        return 0, []
    agree = sum(1 for x in votes if x > 0)
    if agree * 2 == len(votes):
        return 0, seen
    return (1 if agree * 2 > len(votes) else -1), seen


AUTHOR_YEAR_RE = re.compile(r"\((\d{4})\)")


def screen_date(screen: dict) -> str:
    """first_asserted_date from the AUTHOR field, e.g. "Wang T (2014)".

    ORCS exposes no publication date field, and every temporal split on Day 12
    depends on this being real rather than a constant. A screen with no
    parseable year is dated 1970 so it sorts before any split point and can be
    found by query, rather than being given a plausible-looking guess.
    """
    m = AUTHOR_YEAR_RE.search(str(screen.get("AUTHOR", "")))
    return f"{m.group(1)}-01-01" if m else "1970-01-01"


def build_gene_to_proteins(edges_paths) -> dict[str, list[str]]:
    """NCBIGene CURIE -> host protein ids, from the Day 4 ENCODED_BY edges.

    ORCS keys on NCBI Gene; HOST_FACTOR_FOR takes a Protein subject. One gene
    can carry several reviewed proteins, and all of them get the edge: the
    screen knocked out the gene, so the claim applies to every product.
    """
    import json
    out: dict[str, list[str]] = defaultdict(list)
    for path in edges_paths:
        p = Path(path)
        if not p.exists():
            continue
        for line in p.read_text().splitlines():
            if not line.strip():
                continue
            e = json.loads(line)
            if e.get("predicate") == "ENCODED_BY" and e["object"].startswith("NCBIGene:"):
                out[e["object"]].append(e["subject"])
    return out


def ingest_screens(em, screens: list[dict], gene_rows: dict[str, list[dict]],
                   gene_to_protein: dict[str, list[str]], source: str,
                   overrides: dict[str, str] | None = None) -> Counter:
    """Emit HOST_FACTOR_FOR edges from per-screen gene-level results."""
    stats: Counter = Counter()
    overrides = overrides or {}
    seen: set[tuple] = set()

    for screen in screens:
        sid = str(screen["SCREEN_ID"])
        taxon = virus_taxon(screen)
        if taxon is None:
            stats["no_virus_match"] += 1
            continue
        if not is_virus_resistance_screen(screen):
            stats["not_a_resistance_screen"] += 1
            continue

        mode, fixed = polarity_mode(screen)
        if mode == "unresolvable":
            if sid in overrides:
                mode, fixed = "screen_level", overrides[sid]
                stats["override_applied"] += 1
            else:
                stats["polarity_unresolvable"] += 1
                continue

        rows = gene_rows.get(sid)
        if not rows:
            stats["no_gene_data"] += 1
            continue

        hits = [g for g in rows if str(g.get("HIT", "")).upper() == "YES"]
        if not hits:
            stats["no_hits"] += 1
            continue

        # Calibration runs for EVERY screen. A one-sided SCREEN_RATIONALE
        # describes the screen's primary selection; it does not guarantee every
        # hit shares that polarity. Screen 1704 ("Increased resistance to
        # virus") contains IFNAR1, STAT1 and IRF9 -- none can be dependency
        # factors -- and trusting the rationale labelled all three dependency
        # without ever reading a score.
        multiplier = 1
        if mode == "screen_level":
            _m, anchors = calibrate_sign(hits)
            if anchors and _m == -1:
                stats["rationale_contradicted_by_anchors"] += 1
                em.notes[f"rationale_contradicted:screen_{sid}"] += 1
                continue
        if mode == "signed":
            multiplier, _anchors = calibrate_sign(hits)
            if multiplier == 0:
                stats["uncalibrated_sign"] += 1
                em.notes[f"uncalibrated:screen_{sid}"] += 1
                continue
            stats["positive_is_dependency" if multiplier > 0
                  else "positive_is_restriction"] += 1

        cell_line = str(screen.get("CELL_LINE", "")).strip() or "unknown"
        stype = screen_type(screen)
        pmid = str(screen.get("SOURCE_ID", "")).strip()

        for g in hits:
            direction = gene_direction(g.get("SCORE.1"), fixed, multiplier)
            if direction is None:
                stats["no_direction"] += 1
                continue
            gene_curie = f"NCBIGene:{g['IDENTIFIER_ID']}"
            proteins = gene_to_protein.get(gene_curie)
            if not proteins:
                stats["gene_not_in_graph"] += 1
                continue
            try:
                effect = float(g.get("SCORE.1"))
            except (TypeError, ValueError):
                effect = None

            for pid in proteins:
                key = (pid, taxon, sid, direction)
                if key in seen:
                    stats["duplicate"] += 1
                    continue
                seen.add(key)
                quals = {"direction": direction, "screen_type": stype,
                         "cell_line": cell_line}
                if effect is not None:
                    quals["effect_size"] = effect
                em.edge(pid, "HOST_FACTOR_FOR", f"NCBITaxon:{taxon}",
                        source=source, date=screen_date(screen), tier=1, quals=quals,
                        pmids=[f"PMID:{pmid}"] if pmid.isdigit() else None)
                stats["edges"] += 1
                stats[direction] += 1
        stats["screens_used"] += 1

    return stats
