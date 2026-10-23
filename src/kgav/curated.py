"""Host factors from published screens that BioGRID ORCS does not carry.

WHY THIS LAYER EXISTS. The pathway landscape needs one virus-specific
input: HOST_FACTOR_FOR edges. ORCS supplies them for five coronaviruses
and, once its flavivirus patterns are rebuilt, for DENV-2, Zika and West
Nile. It carries nothing for DENV-1, DENV-3, DENV-4, JEV or TBEV.

That gap is mostly not ORCS's. Surveying the flavivirus CRISPR literature
returns screens for DENV, ZIKV, WNV and YFV and none for JEV or TBEV, so
for two of the nine surveyed flaviviruses the data does not exist at
genome scale anywhere. For the rest it exists in supplementary tables,
which is a curation problem rather than an engineering one.

So this reads a table a human filled in, with one row per hit, and emits
the same HOST_FACTOR_FOR edges ORCS does. Adding a paper is then a
transcription task, not a new ingest.

THE DIRECTION IS STATED, NEVER INFERRED. ORCS needs calibrate_sign
because screen score conventions disagree -- the same positive Z means
dependency in one screen and restriction in another. A curated table has
no such excuse: the curator read the paper and knows which way it ran. A
row with no direction is REFUSED rather than defaulted, because
src/kgav/orcs.py's own header records that a sign error here produces a
confidently wrong recommendation rather than a missing one.

TIER 2, AND SEPARATELY SOURCED. These edges carry
`infores:curated-screen`, not `infores:biogrid-orcs`, and evidence tier 2
rather than 1. They are hand-transcribed from a PDF or a spreadsheet by
one person without a second reader, which is a weaker provenance than a
database export, and anything downstream must be able to tell them apart
and discount them.
"""
from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

CURATED_SOURCE = "infores:curated-screen"
CURATED_TIER = 2
DIRECTIONS = ("dependency", "restriction")
REQUIRED = ("ncbi_gene_id", "taxon", "direction", "pmid")
OPTIONAL = ("gene_symbol", "cell_line", "screen_type", "effect_size", "date")


class CuratedRowError(ValueError):
    """A row that cannot become an edge, with the reason and the row number."""


def read_table(path: Path) -> list[dict]:
    """Rows of a curated screen table, as written.

    Tab-separated because gene symbols and cell-line names contain commas
    often enough that a CSV round-trips badly, and a mis-split row here
    becomes a wrong edge rather than a parse error.
    """
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in REQUIRED if c not in (reader.fieldnames or ())]
        if missing:
            raise CuratedRowError(
                f"{path.name}: missing required column(s) {missing}; "
                f"required are {list(REQUIRED)}")
        return [dict(r) for r in reader]


def validate_row(row: dict, n: int) -> dict:
    """One row, checked and normalised, or CuratedRowError naming the row.

    Every check here is a thing that would otherwise become a silently
    wrong edge: a blank direction defaulting to dependency, a taxon that
    is not a taxon, a PMID that is a paper title.
    """
    for col in REQUIRED:
        if not str(row.get(col, "")).strip():
            raise CuratedRowError(f"row {n}: {col} is empty")
    direction = str(row["direction"]).strip().lower()
    if direction not in DIRECTIONS:
        raise CuratedRowError(
            f"row {n}: direction {direction!r} is not one of {DIRECTIONS}. "
            f"It is stated by the curator and never inferred.")
    gene = str(row["ncbi_gene_id"]).strip()
    if not gene.isdigit():
        raise CuratedRowError(
            f"row {n}: ncbi_gene_id {gene!r} is not numeric. A gene SYMBOL "
            f"will not map; resolve it to an NCBI Gene id first.")
    taxon = str(row["taxon"]).strip()
    if not taxon.isdigit():
        raise CuratedRowError(
            f"row {n}: taxon {taxon!r} is not a bare NCBI taxon id")
    pmid = str(row["pmid"]).strip().removeprefix("PMID:")
    if not pmid.isdigit():
        raise CuratedRowError(f"row {n}: pmid {pmid!r} is not numeric")
    effect = str(row.get("effect_size", "")).strip()
    try:
        effect_val = float(effect) if effect else None
    except ValueError:
        raise CuratedRowError(f"row {n}: effect_size {effect!r} is not a number")
    return {
        "gene": f"NCBIGene:{gene}",
        "taxon": taxon,
        "direction": direction,
        "pmid": f"PMID:{pmid}",
        "cell_line": str(row.get("cell_line", "")).strip() or "unknown",
        "screen_type": str(row.get("screen_type", "")).strip() or "unknown",
        "effect_size": effect_val,
        "date": str(row.get("date", "")).strip() or None,
    }


def resolve_symbols(symbols, gene_nodes: dict[str, str]
                    ) -> tuple[dict[str, str], list[str], dict[str, list[str]]]:
    """Gene SYMBOLS -> NCBIGene ids, from the graph's own Gene nodes.

    WHY THIS IS NOT A LOOKUP TABLE SOMEWHERE ELSE. validate_row refuses a
    symbol because build_gene_to_proteins keys on NCBIGene CURIEs, so the
    curator has to supply ids -- and typing a few hundred ids off a
    supplementary table by hand is the single likeliest way a wrong edge
    enters this layer. The host layer already carries `symbol` on every
    Gene node, so the ids come from the graph the edges will land in,
    exactly, with nothing to look up and nothing to guess.

    Returns (resolved, unresolved, ambiguous). An UNRESOLVED symbol is
    worth as much as a resolved one: a screen hit whose gene is not in the
    host layer can never become an edge, and knowing that before
    transcribing the row saves the work. Matching is case-insensitive;
    ALIASES ARE NOT HANDLED -- a screen writing HRD1 where UniProt writes
    SYVN1 comes back unresolved, and that is reported rather than guessed.
    """
    by_symbol: dict[str, list[str]] = {}
    for gid, sym in gene_nodes.items():
        if sym:
            by_symbol.setdefault(sym.strip().upper(), []).append(gid)
    resolved, unresolved, ambiguous = {}, [], {}
    for raw in symbols:
        s = str(raw).strip()
        if not s:
            continue
        hits = by_symbol.get(s.upper())
        if not hits:
            unresolved.append(s)
        elif len(hits) > 1:
            ambiguous[s] = sorted(hits)
        else:
            resolved[s] = hits[0]
    return resolved, unresolved, ambiguous


def ingest_curated(em, rows: list[dict], gene_to_protein: dict[str, list[str]],
                   known_taxa: set[str] | None = None,
                   default_date: str = "1970-01-01") -> Counter:
    """Emit HOST_FACTOR_FOR edges for validated rows.

    `known_taxa` is the set of taxa the spine actually built. A row for a
    taxon with no node would otherwise emit an edge pointing at nothing,
    which validates as a dangling reference much later and reads as a
    graph fault rather than a curation one.
    """
    stats: Counter = Counter()
    seen = set()
    for row in rows:
        if known_taxa is not None and row["taxon"] not in known_taxa:
            stats["taxon_not_built"] += 1
            continue
        proteins = gene_to_protein.get(row["gene"])
        if not proteins:
            stats["gene_not_in_graph"] += 1
            continue
        for pid in proteins:
            key = (pid, row["taxon"], row["pmid"], row["direction"])
            if key in seen:
                stats["duplicate"] += 1
                continue
            seen.add(key)
            quals = {"direction": row["direction"],
                     "screen_type": row["screen_type"],
                     "cell_line": row["cell_line"]}
            if row["effect_size"] is not None:
                quals["effect_size"] = row["effect_size"]
            em.edge(pid, "HOST_FACTOR_FOR", f"NCBITaxon:{row['taxon']}",
                    source=CURATED_SOURCE, date=row["date"] or default_date,
                    tier=CURATED_TIER, quals=quals, pmids=[row["pmid"]])
            stats["edges"] += 1
            stats[row["direction"]] += 1
    return stats
