"""Day 5: virus-host physical interactions from PSI-MITAB 2.7.

This is the layer that makes the graph antiviral rather than generic. It is the
only bridge between the viral spine and the host interactome, so metapath M3
(drug -> host protein -> viral protein -> virus) exists entirely because of it.

WHAT THE EVIDENCE ACTUALLY SAYS. The SARS-CoV-2 detection-method histogram is
dominated by proximity-dependent biotin identification (MI:1314) and affinity
capture (MI:0400) -- roughly 6,600 of ~7,000 rows. BioID labels anything within
about 10 nm, so those edges mean "in the same complex neighbourhood", not
"binds". Binary methods (two-hybrid, pull-down) are a rounding error by
comparison.

That distinction is preserved rather than averaged away:

    detection_method  PSI-MI code from column 7 -- the assay
    interaction_type  PSI-MI code from column 12 -- MI:0407 direct interaction
                      vs MI:0915 physical association (co-complex)
    source_score      VirHostNet miscore from column 15

All of these edges are evidence_tier 1. Tier means HOW THE FACT WAS ESTABLISHED
(curated experimental / curated inferred / computational / text-mined), not how
good the assay was. BioID is a curated experimental result; demoting it to tier
2 would put it alongside GO's curated-inferred annotations, which it has nothing
in common with, and would make Day 12's tier ablation measure an incoherent
mixture. Filtering to direct interactions is a query against interaction_type,
not a re-ingest.

Consequence to carry into the write-up: paths through this layer mean the drug
target sits in the same complex neighbourhood as a viral protein. That is
weaker than binding, and "physically interacts with" should not be allowed to
imply more than the assay supports.
"""
from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

MI_RE = re.compile(r'psi-mi:"?(MI:\d+)"?')
PRO_RE = re.compile(r"(PRO_\d+)")
ACC_RE = re.compile(r"uniprotkb:([A-Z0-9]+)(?:-PRO_\d+)?", re.IGNORECASE)
TAXID_RE = re.compile(r"taxid:(-?\d+)")
PMID_RE = re.compile(r"pubmed:(\d+)")
SCORE_RE = re.compile(r"miscore:([0-9.]+)")

# MI:0407 "direct interaction" asserts a binary contact. MI:0915 "physical
# association" and MI:0914 "association" only assert co-membership of a complex.
DIRECT_INTERACTION = "MI:0407"
# Assays that read out a binary contact rather than co-purification/proximity.
BINARY_METHODS = {"MI:0018", "MI:0096", "MI:0900", "MI:0397", "MI:0055"}


def _first(pattern: re.Pattern, text: str) -> str | None:
    m = pattern.search(text or "")
    return m.group(1) if m else None


def parse_identifier(col: str) -> str | None:
    """MITAB identifier column -> a graph node id.

    A chain reference like uniprotkb:P0DTD1-PRO_0000449624 resolves to the
    CHAIN, not the parent: the source is telling you which mature peptide it
    measured, and collapsing to the polyprotein throws that away.
    """
    if not col or col == "-":
        return None
    chain = _first(PRO_RE, col)
    if chain:
        return f"UniProtKB:{chain}"
    acc = _first(ACC_RE, col)
    return f"UniProtKB:{acc}" if acc else None


def parse_row(line: str) -> dict | None:
    f = line.rstrip("\n").split("\t")
    if len(f) < 15:
        return None
    return {
        "a": parse_identifier(f[0]),
        "b": parse_identifier(f[1]),
        "detection_method": _first(MI_RE, f[6]),
        "pmid": _first(PMID_RE, f[8]),
        "taxid_a": _first(TAXID_RE, f[9]),
        "taxid_b": _first(TAXID_RE, f[10]),
        "interaction_type": _first(MI_RE, f[11]),
        "score": _first(SCORE_RE, f[14]),
    }


@dataclass(frozen=True)
class RowVerdict:
    """Where one MITAB row dies, or that it survives to be emitted.

    WHY THIS IS A FUNCTION RATHER THAN A COMMENT. The gates run in a fixed
    order and a row is counted at the FIRST one it fails, so a reader who
    wants to know what `unresolved_id` contains cannot simply re-test
    resolution: rows that would also fail it may have died at `unparseable`
    or `no_detection_method` first and been counted there.

    diagnose_flavivirus_join.py already carries a paragraph about this --
    "ingest_mitab checks the taxon map LAST, after unparseable,
    no_detection_method, unresolved_id, self_interaction and viral_viral" --
    written after an earlier version of that script asserted rows were lost
    at a gate they never reached. A second implementation of this order is
    how that claim became wrong; there is one implementation now.

    `gate` is None for a row that survives every gate here. Deduplication is
    NOT a gate: it needs state across rows, so ingest_mitab keeps it.
    """
    gate: str | None
    row: dict | None = None
    a: str | None = None
    b: str | None = None
    host: str | None = None
    viral: str | None = None
    viral_tax: str | None = None
    #: For an unresolved_id row: which side(s) failed to resolve, and the
    #: taxid column of each, verbatim from the file.
    unresolved: tuple[tuple[str, str | None, str | None], ...] = ()


def classify_row(line: str, index, taxon_map: dict[str, str]) -> RowVerdict:
    """Run one MITAB line through the ingest's gates, in the ingest's order."""
    r = parse_row(line)
    if not r or not r["a"] or not r["b"]:
        return RowVerdict("unparseable", r)
    if not r["detection_method"]:
        return RowVerdict("no_detection_method", r)

    a, b = index.resolve(r["a"]), index.resolve(r["b"])
    if a is None or b is None:
        # Recorded per side, because "which endpoint failed" is the question
        # the counter alone cannot answer: an unresolved HOST id is lost
        # interactome, an unresolved VIRAL id is a virus this build does not
        # carry, and the single tally reads the same either way.
        bad = tuple((raw, tax, side)
                    for raw, tax, side, got in
                    ((r["a"], r["taxid_a"], "a", a), (r["b"], r["taxid_b"], "b", b))
                    if got is None)
        return RowVerdict("unresolved_id", r, a, b, unresolved=bad)
    if a == b:
        return RowVerdict("self_interaction", r, a, b)

    va, vb = index.is_viral(a), index.is_viral(b)
    if va is None or vb is None:
        return RowVerdict("not_a_protein", r, a, b)
    if va == vb:
        return RowVerdict("viral_viral" if va else "host_host", r, a, b)

    host, viral = (b, a) if va else (a, b)
    viral_tax = r["taxid_a"] if va else r["taxid_b"]
    if viral_tax not in taxon_map:
        return RowVerdict("unmapped_viral_taxon", r, a, b, host, viral, viral_tax)
    return RowVerdict(None, r, a, b, host, viral, viral_tax)


def ingest_mitab(em, paths: list[Path], index, taxon_map: dict[str, str],
                 source: str, default_date: str) -> Counter:
    """Emit host<->viral PHYSICALLY_INTERACTS_WITH edges.

    taxon_map maps every viral taxid that may appear (species AND isolate) to
    the species taxid the graph uses. Swiss-Prot and VirHostNet both file
    several coronaviruses under isolate taxa; without the map those rows look
    like an unknown organism and vanish.
    """
    stats: Counter = Counter()
    seen: set[tuple] = set()

    for path in paths:
        for line in Path(path).read_text(errors="replace").splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            stats["rows"] += 1
            # Every gate, in order, lives in classify_row -- including the
            # ones whose comments used to live here: detection_method is a
            # required qualifier, so an edge without it is dropped rather
            # than defaulted; and viral-viral and host-host pairs are real
            # but out of scope for a BRIDGE layer (host-host belongs to the
            # STRING pass).
            v = classify_row(line, index, taxon_map)
            if v.gate:
                stats[v.gate] += 1
                continue
            r, host, viral = v.row, v.host, v.viral

            key = (host, viral, r["detection_method"], r["pmid"])
            if key in seen:
                stats["duplicate"] += 1
                continue
            seen.add(key)

            # No host_taxon qualifier: it would duplicate taxon_id on the
            # host protein node, and a fact stored twice is a fact that can
            # disagree with itself.
            quals = {"detection_method": r["detection_method"]}
            # interaction_type is NOT stored. VirHostNet normalises every row
            # to MI:0915 "physical association" regardless of assay -- all
            # 8,252 SARS-CoV-2 rows, including all 153 two-hybrid ones. A
            # qualifier with one value across every edge is noise that looks
            # like signal. The assay distinction survives in detection_method.

            em.edge(host, "PHYSICALLY_INTERACTS_WITH", viral,
                    source=source, date=default_date, tier=1, quals=quals,
                    pmids=[f"PMID:{r['pmid']}"] if r["pmid"] else None,
                    score=float(r["score"]) if r["score"] else None)
            stats["edges"] += 1
            stats["direct" if r["interaction_type"] == DIRECT_INTERACTION else "co_complex"] += 1
            if r["detection_method"] in BINARY_METHODS:
                stats["binary_method"] += 1

    return stats
