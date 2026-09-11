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
            r = parse_row(line)
            if not r or not r["a"] or not r["b"]:
                stats["unparseable"] += 1
                continue
            if not r["detection_method"]:
                # Required qualifier; an edge without it cannot be filtered
                # later, so it is dropped rather than defaulted.
                stats["no_detection_method"] += 1
                continue

            a, b = index.resolve(r["a"]), index.resolve(r["b"])
            if a is None or b is None:
                stats["unresolved_id"] += 1
                continue
            if a == b:
                stats["self_interaction"] += 1
                continue

            va, vb = index.is_viral(a), index.is_viral(b)
            if va is None or vb is None:
                stats["not_a_protein"] += 1
                continue
            if va == vb:
                # Viral-viral and host-host pairs are real, but this ingest is
                # the BRIDGE layer. Host-host belongs to the STRING pass;
                # viral-viral is out of scope for v1.
                stats["viral_viral" if va else "host_host"] += 1
                continue

            host, viral = (b, a) if va else (a, b)
            viral_tax = r["taxid_a"] if va else r["taxid_b"]
            if viral_tax not in taxon_map:
                stats["unmapped_viral_taxon"] += 1
                continue

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
