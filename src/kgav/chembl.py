"""Day 7: chemistry from ChEMBL 37.

Three edge types come out of this layer:

    SmallMolecule -INHIBITS->                     Protein (viral)
    SmallMolecule -TARGETS->                      Protein (host)
    SmallMolecule -HAS_ANTIVIRAL_ACTIVITY_AGAINST-> OrganismTaxon

WHAT THE DATA ACTUALLY SUPPORTS, established before writing any of this:

  CYTOTOXICITY IS ABSENT. Two CC50 records across 21,900 coronavirus
  activities. Selectivity index is therefore uncomputable for essentially
  every compound, and the graph cannot distinguish a genuine antiviral from
  something that merely kills the cell. This matters most for the 2,865
  organism-level EC50s, which are cell-based assays where cytotoxicity is
  exactly the confounder. Every activity edge carries unquantified=true.

  MOST VALUES ARE CENSORED. Only 6,600 of 17,366 IC50 records have
  standard_relation '='. The rest are '>' or '<'. "IC50 > 10000 nM" is a
  NEGATIVE result -- the compound was measured and found inactive -- and
  coercing it to a potency number inverts its meaning. Censored records are
  kept with their relation, never silently numeric. They are also the best
  hard negatives available for Day 12, better than mining failed trials,
  because they are measurements rather than inferences about intent.

  CHEMBL HAS NO CHAIN-LEVEL CORONAVIRUS TARGETS. Every activity maps to the
  whole replicase polyprotein (P0DTD1 and friends), so nirmatrelvir's Mpro
  data and remdesivir's RdRp data would land on the same node. The assay
  DESCRIPTION carries what the target record does not: 4,865 of 6,720
  protein-level activities (72%) name a specific nsp. Those resolve to the
  chain; the rest attach to the polyprotein with
  target_resolution=parent_unresolved so an unresolved edge is never mistaken
  for a measured one.

  NSP3 IS TWO DRUGGABLE SITES. PLpro and the ADP-ribose macrodomain live on
  the same protein -- 623 and 1,479 activities respectively. A macrodomain
  binder says nothing about protease inhibition, so the site goes in a
  `domain` qualifier rather than being collapsed into one undifferentiated
  nsp3 target.
"""
from __future__ import annotations

import re
import sqlite3
from collections import Counter
from pathlib import Path

# Description patterns -> (nsp family, domain). Order matters: the macrodomain
# patterns must be tested before the generic nsp3/PLpro ones, since a
# macrodomain assay description also mentions nsp3.
NSP_PATTERNS: list[tuple[re.Pattern, str, str | None]] = [
    (re.compile(r"macrodomain|mac1\b|adp[- ]ribose", re.IGNORECASE), "nsp3", "macrodomain"),
    (re.compile(r"papain|pl[- ]?pro\b|pl protease|plpro", re.IGNORECASE), "nsp3", "PLpro"),
    (re.compile(r"3cl|3c-like|3-cl|main protease|mpro|m-protease|m protease", re.IGNORECASE), "nsp5", None),
    (re.compile(r"rdrp|rna[- ]dependent rna pol|rna[- ]directed rna pol|nsp12", re.IGNORECASE), "nsp12", None),
    (re.compile(r"helicase|nsp13", re.IGNORECASE), "nsp13", None),
    (re.compile(r"exoribonuclease|nsp14", re.IGNORECASE), "nsp14", None),
    (re.compile(r"endoribonuclease|nendou|nsp15", re.IGNORECASE), "nsp15", None),
    (re.compile(r"2'-o-methyl|nsp16", re.IGNORECASE), "nsp16", None),
    (re.compile(r"nsp1\b", re.IGNORECASE), "nsp1", None),
]

QUANT_TYPES = {"IC50": "ic50_nm", "EC50": "ec50_nm", "Ki": "ki_nm", "Kd": "kd_nm"}
VALID_RELATIONS = {"=", ">", "<", ">=", "<=", "~"}


def resolve_nsp(description: str) -> tuple[str | None, str | None]:
    """(nsp family, domain) named by an assay description, or (None, None)."""
    for pattern, nsp, domain in NSP_PATTERNS:
        if pattern.search(description or ""):
            return nsp, domain
    return None, None


def to_nm(value: float, units: str | None) -> float | None:
    """Normalise to nanomolar, or None when the unit is not a concentration.

    ug.mL-1 is deliberately NOT converted: that needs the molecular weight, and
    silently guessing one produces a number that looks measured.
    """
    if value is None:
        return None
    u = (units or "").strip()
    if u == "nM":
        return float(value)
    if u == "uM":
        return float(value) * 1000
    if u == "pM":
        return float(value) / 1000
    if u == "mM":
        return float(value) * 1_000_000
    return None


def build_viral_target_index(index, taxon_map: dict[str, str]) -> dict[str, dict[str, str]]:
    """taxon -> {nsp family: chain node id}, from the spine's protein_family."""
    out: dict[str, dict[str, str]] = {}
    for n in index.nodes.values():
        if n["class"] != "Protein":
            continue
        p = n["properties"]
        if not p.get("is_viral") or not p.get("protein_family"):
            continue
        tax = (p.get("taxon_id") or "").split(":")[-1]
        if tax in taxon_map:
            out.setdefault(taxon_map[tax], {})[p["protein_family"]] = n["id"]
    return out


def query_activities(db: Path, taxa: list[str]) -> list[dict]:
    """Every usable coronavirus activity, with its assay description."""
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    placeholders = ",".join("?" * len(taxa))
    rows = con.execute(f"""
        SELECT md.chembl_id AS compound, cs_struct.canonical_smiles AS smiles,
               cs_struct.standard_inchi_key AS inchikey,
               md.max_phase, md.pref_name AS drug_name,
               t.chembl_id AS target, t.pref_name AS target_name,
               t.target_type, t.tax_id,
               comp.accession,
               act.standard_type, act.standard_value, act.standard_units,
               act.standard_relation,
               a.description, a.assay_type, a.chembl_id AS assay,
               d.year, d.pubmed_id
        FROM activities act
        JOIN assays a          ON act.assay_id = a.assay_id
        JOIN target_dictionary t ON a.tid = t.tid
        JOIN molecule_dictionary md ON act.molregno = md.molregno
        LEFT JOIN compound_structures cs_struct ON md.molregno = cs_struct.molregno
        LEFT JOIN docs d       ON a.doc_id = d.doc_id
        LEFT JOIN target_components tc ON t.tid = tc.tid
        LEFT JOIN component_sequences comp ON tc.component_id = comp.component_id
        WHERE t.tax_id IN ({placeholders})
          AND act.standard_value IS NOT NULL
          AND act.standard_type IN ('IC50','EC50','Ki','Kd')
    """, taxa).fetchall()
    con.close()
    return [dict(r) for r in rows]


def ingest_activities(em, rows: list[dict], viral_targets: dict[str, dict[str, str]],
                      taxon_map: dict[str, str], source: str) -> Counter:
    stats: Counter = Counter()
    seen: set[tuple] = set()

    for r in rows:
        stats["rows"] += 1
        inchikey = r.get("inchikey")
        if not inchikey:
            stats["no_inchikey"] += 1
            continue
        tax = taxon_map.get(str(r["tax_id"]))
        if tax is None:
            stats["unmapped_taxon"] += 1
            continue

        qual_field = QUANT_TYPES.get(r["standard_type"])
        nm = to_nm(r["standard_value"], r["standard_units"])
        if nm is None:
            stats["unconvertible_units"] += 1
            continue
        relation = r["standard_relation"] if r["standard_relation"] in VALID_RELATIONS else "="
        if relation != "=":
            stats["censored"] += 1

        drug = f"INCHIKEY:{inchikey}"
        em.node(drug, "SmallMolecule",
                smiles=r.get("smiles") or "",
                inchikey_skel=inchikey[:14],
                chembl_id=r["compound"],
                max_phase=int(r["max_phase"]) if r.get("max_phase") is not None else None,
                is_approved=bool(r.get("max_phase") == 4),
                salt_collapsed=False, stereo_collapsed=False)

        date = f"{int(r['year'])}-01-01" if r.get("year") else "1970-01-01"
        pmids = [f"PMID:{r['pubmed_id']}"] if r.get("pubmed_id") else None
        quals = {qual_field: nm, "relation": relation,
                 # CC50 is absent from ChEMBL's coronavirus data (2 records in
                 # 21,900), so no activity here has a selectivity index.
                 "unquantified": True}

        if r["target_type"] == "ORGANISM":
            # Ki and Kd are binding constants against a purified enzyme; they
            # are not meaningful against a whole organism in a cell-based
            # assay. Two such records exist in ChEMBL 37 and are mis-entered.
            # Dropped rather than accommodated -- widening the schema to admit
            # them would let a real error through unnoticed later.
            if r["standard_type"] in ("Ki", "Kd"):
                stats["binding_constant_on_organism"] += 1
                continue
            quals["assay_type"] = "cell_based_antiviral"
            key = (drug, "HAS_ANTIVIRAL_ACTIVITY_AGAINST", f"NCBITaxon:{tax}",
                   r["assay"], r["standard_type"])
            if key in seen:
                stats["duplicate"] += 1
                continue
            seen.add(key)
            em.edge(drug, "HAS_ANTIVIRAL_ACTIVITY_AGAINST", f"NCBITaxon:{tax}",
                    source=source, date=date, tier=1, quals=quals, pmids=pmids)
            stats["organism_edges"] += 1
            continue

        # Multi-component targets (CHEMBL4888460 "Spike glycoprotein/ACE2",
        # CHEMBL6067603 "cereblon/Replicase polyprotein 1ab") map to BOTH a
        # viral and a host accession, so the join returns whichever row comes
        # first. A spike-ACE2 disruption assay is neither INHIBITS on a viral
        # protein nor TARGETS on a host one -- it measures complex disruption.
        # Skipped in v1 rather than forced into a predicate that overstates it.
        if "/" in (r.get("target_name") or ""):
            stats["multicomponent_target"] += 1
            continue

        nsp, domain = resolve_nsp(r.get("description") or "")
        target_node = viral_targets.get(tax, {}).get(nsp) if nsp else None
        if target_node:
            quals["target_resolution"] = "chain_from_assay_description"
            stats[f"resolved_{nsp}"] += 1
        else:
            acc = r.get("accession")
            if not acc:
                stats["no_target_node"] += 1
                continue
            target_node = f"UniProtKB:{acc}"
            quals["target_resolution"] = "parent_unresolved"
            stats["unresolved"] += 1
        if domain:
            quals["domain"] = domain

        quals["assay_type"] = "biochemical"
        key = (drug, "INHIBITS", target_node, r["assay"], r["standard_type"])
        if key in seen:
            stats["duplicate"] += 1
            continue
        seen.add(key)
        em.edge(drug, "INHIBITS", target_node,
                source=source, date=date, tier=1, quals=quals, pmids=pmids)
        stats["protein_edges"] += 1

    return stats
