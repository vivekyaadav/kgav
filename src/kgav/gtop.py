"""Drug -> human protein TARGETS from the IUPHAR/BPS Guide to PHARMACOLOGY.

WHY THIS SOURCE AND NOT ANOTHER. Two gaps in the graph are measured, and
this closes the one that matters more.

  DIRECTION. 12,155 of 14,381 TARGETS edges state no direction, because
  ChEMBL's bulk activity table records potency and only drug_mechanism
  carries action_type, for a minority of approved drugs. So "this drug
  targets a dependency factor" cannot become "inhibiting it should block
  infection". GtP curates `action` on essentially every interaction, which
  is the field that gap is made of.

  COVERAGE. 18,665 of 21,667 small molecules carry no human target. GtP
  adds clinical and approved ligands, so the census's approved-plus-
  clinical count is the ceiling on what it can reach -- not the whole gap,
  most of which is bulk screening compounds no curated source covers.

REDISTRIBUTABLE, WHICH DRUGBANK IS NOT. A graph published beside a paper
has to ship, and DrugBank's academic licence restricts that while GtP's
downloads do not. Confirm the licence text against the release you
actually download; this module states the reason for the choice, not a
warranty about a file it has never seen.

COMPOUND IDENTITY IS THE JOIN, AND IT IS NOT FREE. The graph keys small
molecules on INCHIKEY:<standard key>. interactions.csv carries GtP ligand
ids and a PubChem SID, not an InChIKey, so ligands.csv (or
ligand_id_mapping.csv) has to supply it. A ligand with no InChIKey is
SKIPPED AND COUNTED: without one it cannot be identified against a
compound already in the graph, and minting a new node on a name would
create a duplicate of a drug that is already there under its key.

BUILT AGAINST THE PUBLISHED COLUMN LIST, NOT AGAINST A FILE. The schema
comes from the GtP download documentation. Real files drift from docs, so
every reader checks its columns and refuses by name rather than reading a
shifted field as a direction.
"""
from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

GTOP_SOURCE = "infores:guidetopharmacology"
GTOP_TIER = 1
HUMAN = "Human"

# GtP `action` -> the schema's target_direction enum, which
# src/kgav/host_targets.py already populates from ChEMBL action_type. One
# vocabulary, so a direction means the same thing whichever source wrote
# it.
#
# WHAT IS DELIBERATELY ABSENT. GtP also uses Binding, Unknown, None and
# several context-dependent terms. They map to "unknown" by falling
# through, because a binder is not an inhibitor and recording it as one
# would manufacture the mechanistic claim this layer exists to supply
# honestly.
ACTION_TO_DIRECTION = {
    "INHIBITOR": "inhibitor",
    "ANTAGONIST": "inhibitor",
    "CHANNEL BLOCKER": "inhibitor",
    "GATING INHIBITOR": "inhibitor",
    "NEGATIVE": "inhibitor",
    "IRREVERSIBLE INHIBITION": "inhibitor",
    "ANTISENSE INHIBITION": "inhibitor",
    "AGONIST": "agonist",
    "PARTIAL AGONIST": "agonist",
    "FULL AGONIST": "agonist",
    "ACTIVATOR": "agonist",
    "POSITIVE": "agonist",
    "BIASED AGONIST": "agonist",
    "ALLOSTERIC MODULATOR": "modulator",
    "ALLOSTERIC ENHANCER": "modulator",
    "MODULATOR": "modulator",
    "MIXED": "modulator",
}

INTERACTION_COLUMNS = ("ligand_id", "target_uniprot", "target_species",
                       "action", "pubmed_id")
LIGAND_COLUMNS = ("ligand_id",)
INCHIKEY_COLUMNS = ("inchikey", "InChIKey", "InChIKey ")


class GtopFormatError(ValueError):
    """A file whose columns are not what the documentation describes."""


def _rows(path: Path, required: tuple[str, ...]) -> list[dict]:
    """CSV rows, after checking the columns the reader will use exist.

    GtP files lead with a '#' comment line holding the version, so that is
    skipped before the header is read.
    """
    with open(path, newline="", encoding="utf-8-sig") as fh:
        lines = [l for l in fh if not l.startswith("#")]
    reader = csv.DictReader(lines)
    have = set(reader.fieldnames or ())
    missing = [c for c in required if c not in have]
    if missing:
        raise GtopFormatError(
            f"{path.name}: missing column(s) {missing}. The reader was built "
            f"against the published column list; this file has "
            f"{sorted(have)[:8]}... Check the release rather than letting a "
            f"shifted field be read as a direction.")
    return [dict(r) for r in reader]


def ligand_inchikeys(path: Path) -> dict[str, str]:
    """{GtP ligand_id: InChIKey}, from ligands.csv or ligand_id_mapping.csv.

    Whichever file carries the key, the column name has varied across
    releases, so several spellings are accepted and the absence of all of
    them is an error naming what was looked for.
    """
    rows = _rows(path, LIGAND_COLUMNS)
    col = next((c for c in INCHIKEY_COLUMNS
                if rows and c in rows[0]), None)
    if col is None:
        raise GtopFormatError(
            f"{path.name}: no InChIKey column; looked for "
            f"{list(INCHIKEY_COLUMNS)}. Without one no ligand can be "
            f"matched to a compound already in the graph.")
    out = {}
    for r in rows:
        key = str(r.get(col) or "").strip()
        lid = str(r.get("ligand_id") or "").strip()
        if lid and key:
            out[lid] = key
    return out


def direction_of(action: str) -> str:
    return ACTION_TO_DIRECTION.get((action or "").strip().upper(), "unknown")


def ingest_interactions(em, rows: list[dict], inchikeys: dict[str, str],
                        known_proteins: set[str] | None = None,
                        known_compounds: set[str] | None = None,
                        default_date: str = "1970-01-01") -> Counter:
    """Emit TARGETS edges for human interactions whose ligand has a key.

    `known_compounds`, when given, restricts to molecules already in the
    graph. Minting new ones here is possible and deliberately not done by
    default: this layer exists to add DIRECTION and reach to compounds the
    graph has, and a new node keyed on an InChIKey the chemistry layer
    never normalised would sit beside its own duplicate.
    """
    stats: Counter = Counter()
    seen = set()
    for r in rows:
        if (r.get("target_species") or "").strip() != HUMAN:
            stats["not_human"] += 1
            continue
        acc = (r.get("target_uniprot") or "").strip()
        if not acc:
            stats["no_uniprot"] += 1
            continue
        # A pipe-separated accession list means a complex, not one protein.
        if "|" in acc:
            stats["multi_accession"] += 1
            continue
        pid = f"UniProtKB:{acc}"
        if known_proteins is not None and pid not in known_proteins:
            stats["protein_not_in_graph"] += 1
            continue
        lid = (r.get("ligand_id") or "").strip()
        key = inchikeys.get(lid)
        if not key:
            stats["ligand_without_inchikey"] += 1
            continue
        cid = f"INCHIKEY:{key}"
        if known_compounds is not None and cid not in known_compounds:
            stats["compound_not_in_graph"] += 1
            continue
        direction = direction_of(r.get("action"))
        if (cid, pid) in seen:
            stats["duplicate"] += 1
            continue
        seen.add((cid, pid))
        quals = {"direction": direction}
        affinity = (r.get("original_affinity_median_nm") or "").strip()
        try:
            if affinity:
                quals["affinity_nm"] = float(affinity)
        except ValueError:
            stats["unparseable_affinity"] += 1
        pmids = [f"PMID:{p.strip()}"
                 for p in str(r.get("pubmed_id") or "").split("|")
                 if p.strip().isdigit()]
        em.edge(cid, "TARGETS", pid, source=GTOP_SOURCE, date=default_date,
                tier=GTOP_TIER, quals=quals, pmids=pmids or None)
        stats["edges"] += 1
        stats[f"direction_{direction}"] += 1
    return stats
