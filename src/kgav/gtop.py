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

# THE COLUMNS, AS THE DOWNLOAD ACTUALLY SPELLS THEM.
#
# These were wrong until 2026-11-01 and so was every fixture in
# tests/test_gtop.py. The reader was built from the published column list
# -- snake_case names like `target_uniprot` -- and the tests were written
# to match the reader, so 286 lines of them passed against a format the
# site has never served. The first real download refused on the first
# column check, which is the only reason it was caught.
#
# Spelled out rather than normalised by lower-casing, because
# interactions.csv carries BOTH `Type` and `Ligand Type` and they are
# different things: `Type` is the interaction, `Ligand Type` is the
# molecule class. A reader that matched loosely could take either.
C_LIGAND_ID = "Ligand ID"
C_UNIPROT = "Target UniProt ID"
C_SPECIES = "Target Species"
C_TYPE = "Type"
C_ACTION = "Action"
C_PMID = "PubMed ID"
C_AFFINITY = "Original Affinity Median nm"

INTERACTION_COLUMNS = (C_LIGAND_ID, C_UNIPROT, C_SPECIES, C_TYPE,
                       C_ACTION, C_PMID)
LIGAND_COLUMNS = (C_LIGAND_ID,)
INCHIKEY_COLUMNS = ("InChIKey", "inchikey", "InChIKey ")

# GtoP `Type` -> the schema's target_direction enum, which
# src/kgav/host_targets.py already populates from ChEMBL action_type. One
# vocabulary, so a direction means the same thing whichever source wrote
# it.
#
# `Type`, NOT `Action`. This map was always the Type vocabulary --
# Inhibitor, Activator, Channel blocker, Gating inhibitor and Allosteric
# modulator are Type values, none of them Action values -- while the
# reader passed it `Action`, whose largest value in 2026.3.1 is
# `Inhibition` at 8,159 rows and does not appear here. Every one of those
# would have fallen through to `unknown`: the layer would have run,
# reported success, and supplied almost no inhibitor directions, which is
# the one direction an antiviral repurposing graph is for.
#
# WHAT IS DELIBERATELY ABSENT. None, Antibody, Fusion protein and
# Subunit-specific fall through to "unknown" -- 1,600 rows of 24,968 --
# because none of them states a direction. An antibody may antagonise or
# may only bind, and recording a guess would manufacture the mechanistic
# claim this layer exists to supply honestly.
TYPE_TO_DIRECTION = {
    "INHIBITOR": "inhibitor",
    "ANTAGONIST": "inhibitor",
    "CHANNEL BLOCKER": "inhibitor",
    "GATING INHIBITOR": "inhibitor",
    "AGONIST": "agonist",
    "ACTIVATOR": "agonist",
    "ALLOSTERIC MODULATOR": "modulator",
}

# `Action` ASSIGNS NOTHING. It is read for exactly one purpose: to notice
# the source contradicting itself, and to refuse a direction when it does.
#
# In 2026.3.1 that happens 136 times -- 77 rows reading Type=Antagonist
# with Action=Partial agonist, and 59 reading Type=Agonist with
# Action=Inverse agonist. Taking `Type` there would assert the opposite of
# what the same row's other column says, so those get "unknown".
#
# ALLOSTERIC MODULATORS ARE NOT A CONTRADICTION and are exempt below. For
# them the two columns are complementary rather than competing: `Type`
# records that the compound binds an allosteric site, `Action` records
# which way it pushes -- Positive 364 rows, Negative 355. Decided
# 2026-11-01 to record all of them as `modulator` and keep `Action` out of
# it: a positive allosteric modulator has no intrinsic activity of its
# own, so calling it an agonist would overstate what GtoP asserts. The
# cost is real and recorded here -- M2, M4 and M5 see `modulator` and
# cannot tell a PAM from a NAM.
ACTION_TO_DIRECTION = {
    "INHIBITION": "inhibitor",
    "ANTAGONIST": "inhibitor",
    "NEGATIVE": "inhibitor",
    "IRREVERSIBLE INHIBITION": "inhibitor",
    "PORE BLOCKER": "inhibitor",
    "FEEDBACK INHIBITION": "inhibitor",
    "VOLTAGE-DEPENDENT INHIBITION": "inhibitor",
    "INVERSE AGONIST": "inhibitor",
    "AGONIST": "agonist",
    "FULL AGONIST": "agonist",
    "PARTIAL AGONIST": "agonist",
    "BIASED AGONIST": "agonist",
    "IRREVERSIBLE AGONIST": "agonist",
    "ACTIVATION": "agonist",
    "POSITIVE": "agonist",
    "POTENTIATION": "agonist",
    "MIXED": "modulator",
}

#: The first line of a GtoP export, which is a comment and is QUOTED.
COMMENT_PREFIXES = ("#", '"#')


class GtopFormatError(ValueError):
    """A file whose columns are not what the documentation describes."""


def _version_and_body(path: Path) -> tuple[str | None, list[str]]:
    """The leading version comment, and the lines after it.

    GtoPdb writes the version as the first line and QUOTES it:

        "# GtoPdb Version: 2026.3.1 - published: 2026-10-12"

    The previous reader filtered `if not l.startswith("#")`, which that
    line does not -- it starts with a double quote -- so DictReader read
    the comment as the header and every column check failed against it.
    tests/test_gtop.py wrote its fixture as an UNQUOTED `#GtP v2026.1`, so
    the test for this passed while the real form could not be read at all.

    ONLY THE FIRST LINE is treated as a comment; a filter applied to every
    line would also drop any data row that happened to begin that way. The
    version is returned rather than discarded because DATA_LICENSES.md
    requires the release recorded in the MANIFEST, and this is the release.
    """
    with open(path, newline="", encoding="utf-8-sig") as fh:
        lines = fh.readlines()
    if lines and lines[0].lstrip().startswith(COMMENT_PREFIXES):
        return lines[0].strip().strip('"').lstrip("#").strip(), lines[1:]
    return None, lines


def file_version(path: Path) -> str | None:
    """The GtoPdb release string a download declares, or None."""
    return _version_and_body(path)[0]


def _rows(path: Path, required: tuple[str, ...]) -> list[dict]:
    """CSV rows, after checking the columns the reader will use exist."""
    reader = csv.DictReader(_version_and_body(path)[1])
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
        lid = str(r.get(C_LIGAND_ID) or "").strip()
        if lid and key:
            out[lid] = key
    return out


def direction_of(type_value: str, action_value: str = "") -> str:
    """The schema direction a row asserts, or "unknown".

    `Type` decides. `Action` is consulted only to refuse: where the two
    columns map to opposite directions the row contradicts itself, and
    "unknown" is the honest reading of it rather than whichever column is
    checked first.

    Allosteric modulators are exempt from that test by design. There
    `Action` is mechanism detail, not a competing claim, and all of them
    record `modulator` -- see ACTION_TO_DIRECTION's note for why a
    positive allosteric modulator is not filed as an agonist.
    """
    t = TYPE_TO_DIRECTION.get((type_value or "").strip().upper())
    if t is None:
        return "unknown"
    if t == "modulator":
        return "modulator"
    a = ACTION_TO_DIRECTION.get((action_value or "").strip().upper())
    return "unknown" if a is not None and a != t else t


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
        if (r.get(C_SPECIES) or "").strip() != HUMAN:
            stats["not_human"] += 1
            continue
        acc = (r.get(C_UNIPROT) or "").strip()
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
        lid = (r.get(C_LIGAND_ID) or "").strip()
        key = inchikeys.get(lid)
        if not key:
            stats["ligand_without_inchikey"] += 1
            continue
        cid = f"INCHIKEY:{key}"
        if known_compounds is not None and cid not in known_compounds:
            stats["compound_not_in_graph"] += 1
            continue
        direction = direction_of(r.get(C_TYPE), r.get(C_ACTION))
        if direction == "unknown" and TYPE_TO_DIRECTION.get(
                (r.get(C_TYPE) or "").strip().upper()) is not None:
            # Type mapped, and Action disagreed with it. Counted so the
            # driver can report how often the source contradicts itself
            # rather than hiding it inside direction_unknown.
            stats["type_action_conflict"] += 1
        if (cid, pid) in seen:
            stats["duplicate"] += 1
            continue
        seen.add((cid, pid))
        quals = {"direction": direction}
        affinity = (r.get(C_AFFINITY) or "").strip()
        try:
            if affinity:
                quals["affinity_nm"] = float(affinity)
        except ValueError:
            stats["unparseable_affinity"] += 1
        pmids = [f"PMID:{p.strip()}"
                 for p in str(r.get(C_PMID) or "").split("|")
                 if p.strip().isdigit()]
        em.edge(cid, "TARGETS", pid, source=GTOP_SOURCE, date=default_date,
                tier=GTOP_TIER, quals=quals, pmids=pmids or None)
        stats["edges"] += 1
        stats[f"direction_{direction}"] += 1
    return stats
