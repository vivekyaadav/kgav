"""Day 7: chemistry from ChEMBL 37.

Three edge types come out of this layer:

    SmallMolecule -INHIBITS->                     Protein (viral)
    SmallMolecule -MEASURED_INACTIVE_AGAINST->    Protein (viral)
    SmallMolecule -TARGETS->                      Protein (host)
    SmallMolecule -HAS_ANTIVIRAL_ACTIVITY_AGAINST-> OrganismTaxon

  A viral-protein row becomes INHIBITS or MEASURED_INACTIVE_AGAINST according
  to labels.classify -- the function that also decides the labels, so the
  features and the labels cannot disagree about what "acts on" means.

WHAT THE DATA ACTUALLY SUPPORTS, established before writing any of this:

  CYTOTOXICITY IS ABSENT. Two CC50 records across 21,900 coronavirus
  activities. Selectivity index is therefore uncomputable for essentially
  every compound, and the graph cannot distinguish a genuine antiviral from
  something that merely kills the cell. This matters most for the 2,865
  organism-level EC50s, which are cell-based assays where cytotoxicity is
  exactly the confounder. (This layer used to record that as
  unquantified=true on every edge it wrote. It no longer states the field at
  all: the selectivity layer supplies CC50 for 1,568 of these pairs from
  same-document evidence, so the claim was false by the time the graph was
  assembled. The field is derived after merging -- see the schema.)

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

from kgav.labels import classify

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


# WHAT WAS THE MEASUREMENT MADE ON. ChEMBL's assay_type code does not answer
# this: measured on the 7,677 coronavirus protein-target rows, 444 B-coded
# ("binding") rows describe a cellular readout, including "Inhibition of eGFP
# fused SARS-CoV2 main protease transfected in HEK293T/17 cells ... by
# fluorescence based flow cytometry". The description answers it.
#
# This matters because the predicate depends on it. A pseudovirus entry assay
# against spike does not show a compound binds spike; it shows the compound
# blocks spike-mediated entry. Chloroquine's "IC50 160 nM on Spike
# glycoprotein" is exactly that -- "spike glycoprotein S in SARS-CoV-2
# pseudovirus infected in human 293T/ACE2 cells assessed as inhibition of
# viral infection" -- and it blocks endosomal acidification, not spike. The
# same compound's real spike measurement, an ACE2-interaction ELISA, is
# 7,000 nM. Filed on the protein, the surrogate wins the best-potency
# selection by 44x and becomes the headline the agent reports.
SURROGATE_PATTERNS = (
    r"pseudo[- ]?typed", r"pseudovirus", r"pseudo[- ]?particle",
    r"\breplicon\b", r"\bluciferase\b", r"\beGFP\b", r"\bGFP\b",
    r"reporter", r"transfect",
)
LIVE_VIRUS_PATTERNS = (
    r"cytopath", r"\bCPE\b", r"\binfect", r"viral load", r"virus yield",
    r"\bvirus titer\b", r"\btiter\b",
)
PLAQUE_PATTERNS = (r"plaque",)
# A cell line name alone makes a readout cellular. Kept separate from the
# above because it says WHERE, not WHAT.
CELL_LINE_PATTERNS = (
    r"\bvero\b", r"calu-?3", r"\bhuh-?7", r"\bcaco-?2", r"\bhek ?293",
    r"\b293t", r"\bA549\b", r"\bcells?\b",
)
# Protein-level, even when the material came from cells. A lysate's
# description necessarily says "cell", so without these it reads as cellular.
BIOCHEMICAL_PATTERNS = (
    r"recombinant", r"purified", r"enzymatic", r"\bFRET\b", r"fluorogenic",
    r"\bSPR\b", r"thermal shift", r"\bITC\b", r"cell[- ]free", r"lysate",
    r"\bextract\b", r"\bin vitro\b", r"\bELISA\b",
)

_SURROGATE = re.compile("|".join(SURROGATE_PATTERNS), re.I)
_LIVE = re.compile("|".join(LIVE_VIRUS_PATTERNS), re.I)
_PLAQUE = re.compile("|".join(PLAQUE_PATTERNS), re.I)
_CELL = re.compile("|".join(CELL_LINE_PATTERNS), re.I)
_BIOCHEM = re.compile("|".join(BIOCHEMICAL_PATTERNS), re.I)


def assay_readout(description: str, assay_type_code: str | None = None) -> tuple[str, str]:
    """(where, assay_type) for one assay description.

    `where` is "cells", "protein" or "unclear", and decides whether the row
    describes the protein it is filed against. `assay_type` is one of the
    schema's enum values -- which already contains reporter,
    plaque_reduction and cell_based_antiviral. The vocabulary was there; the
    ingest wrote "biochemical" on every protein row regardless.

    "unclear" is returned when the text says both, as in a recombinant
    protease assayed in a cell lysate. 537 rows read that way and they are
    left on the protein: moving a measurement requires knowing it was
    cellular, and these are the ones a human has to read.
    """
    d = description or ""
    bio = bool(_BIOCHEM.search(d))
    surrogate = bool(_SURROGATE.search(d))
    live = bool(_LIVE.search(d))
    cellular = surrogate or live or bool(_CELL.search(d))

    if cellular and bio:
        return "unclear", "other"
    if cellular:
        if _PLAQUE.search(d):
            return "cells", "plaque_reduction"
        if surrogate:
            return "cells", "reporter"
        return "cells", "cell_based_antiviral"
    if bio:
        return "protein", "biochemical"
    # Nothing either way: 1,947 rows, all B-coded. Trust the code, which is
    # the weaker signal but the only one left, and say "binding" rather than
    # "biochemical" because that is what B means.
    return ("protein", "binding" if (assay_type_code or "B").upper() == "B"
            else "biochemical")


# EVERY ROW'S FATE. Kept here, beside the code that increments the counters,
# because they lived in the driver and drifted from it: the accounting check
# failed by exactly the 1,238 rows a new counter tracked, the layer was never
# written, and the reassembly that followed silently used the previous one.
#
# EMITTED and DROPPED must partition every row. BREAKDOWN_PREFIXES name
# counters that DESCRIBE rows instead of deciding their fate -- they are
# incremented before dedup and classify, so they are a superset of what is
# emitted and must not enter the sum.
# Qualifiers that mean something on an ORGANISM edge. An ALLOWLIST, not a
# copy of the protein row's dict: copying carried target_resolution and domain
# -- which describe how the PROTEIN was resolved out of the assay text and say
# nothing once the edge is on the organism -- and produced 2,110
# UNDECLARED_QUALIFIER violations, so the layer was not written.
#
# The schema is deliberately not widened to admit them. nominal_protein_target
# already records WHICH protein the source named, which is the auditable fact;
# how its accession was derived belongs to a protein edge. Declaring
# protein-resolution fields on the organism predicate to keep a detail is the
# widening this module already refuses elsewhere, on the grounds that it lets
# a real error through unnoticed later.
ORGANISM_QUALIFIERS = ("ec50_nm", "ic50_nm", "cc50_nm", "relation",
                       "assay_type", "cell_line", "selectivity_index",
                       "nominal_protein_target")

EMITTED_COUNTERS = ("protein_edges", "measured_inactive_edges",
                    "organism_edges", "reattributed_to_organism")
DROPPED_COUNTERS = ("no_inchikey", "unmapped_taxon", "unconvertible_units",
                    "no_target_node", "duplicate", "multicomponent_target",
                    "binding_constant_on_organism",
                    "undecidable_protein_measurement")
BREAKDOWN_COUNTERS = ("rows", "censored", "unresolved")
BREAKDOWN_PREFIXES = ("resolved_", "readout_", "assay_type_")


def accounting(stats) -> tuple[int, int, list[str]]:
    """(emitted, dropped, counters whose fate is unstated).

    An unaccounted row is a row whose fate nobody can state, which is how
    1,433 pairs spent four months asserting inhibition they had no
    measurement for.
    """
    emitted = sum(stats[k] for k in EMITTED_COUNTERS)
    dropped = sum(stats[k] for k in DROPPED_COUNTERS)
    unlisted = sorted(
        set(stats) - set(EMITTED_COUNTERS) - set(DROPPED_COUNTERS)
        - set(BREAKDOWN_COUNTERS)
        - {k for k in stats if k.startswith(BREAKDOWN_PREFIXES)})
    return emitted, dropped, unlisted


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
        # `unquantified` is NOT set here. This layer once wrote true on every
        # edge, reasoning that CC50 is absent from ChEMBL's coronavirus
        # activities (2 records in 21,900) so no selectivity index is
        # computable. The reasoning was about THIS layer; the field is about
        # the assembled edge, and the selectivity layer later supplied the
        # CC50 that made the claim false. Derived at assembly instead.
        quals = {qual_field: nm, "relation": relation}

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

        # THE MEASUREMENT DECIDES THE PREDICATE, AND SO DOES WHAT IT WAS MADE
        # ON. The rule below already stopped a non-active measurement becoming
        # INHIBITS; this is the same rule one level deeper. 1,884 of 7,677
        # protein-target rows (24.5%) describe a CELLULAR readout, so filing
        # them on the protein asserts a direct interaction the assay never
        # measured -- and M1, M7 and M8 walk those edges, which means the
        # "direct-acting" channel was never purely direct-acting.
        #
        # A cellular row is emitted on the ORGANISM instead, which is where
        # the measurement points: it is an antiviral activity measurement that
        # happens to name a protein as the assay's nominal target. That makes
        # it a held-out label rather than traversable evidence, which is
        # correct -- a cell-based EC50 is exactly what this project evaluates
        # against.
        where, assay_type = assay_readout(r.get("description") or "",
                                          r.get("assay_type"))
        quals["assay_type"] = assay_type
        stats[f"readout_{where}"] += 1
        stats[f"assay_type_{assay_type}"] += 1

        if where == "cells":
            # SAME RULE AS THE ORGANISM BRANCH ABOVE. A Ki or Kd is a binding
            # constant against a purified enzyme and is not meaningful against
            # a whole organism -- and a row reporting one while describing a
            # cellular readout contradicts itself twice over. The existing
            # branch drops these; this path bypassed that rule, which is how a
            # reattribution would have carried ki_nm onto an organism edge.
            if r["standard_type"] in ("Ki", "Kd"):
                stats["binding_constant_on_organism"] += 1
                continue
            # The nominal protein target is kept so the reattribution is
            # auditable rather than silent: it is why this row exists at all.
            quals["nominal_protein_target"] = target_node
            org_quals = {k: v for k, v in quals.items()
                         if k in ORGANISM_QUALIFIERS}
            org_key = (drug, "HAS_ANTIVIRAL_ACTIVITY_AGAINST",
                       f"NCBITaxon:{tax}", r["assay"], r["standard_type"])
            if org_key in seen:
                stats["duplicate"] += 1
                continue
            seen.add(org_key)
            em.edge(drug, "HAS_ANTIVIRAL_ACTIVITY_AGAINST", f"NCBITaxon:{tax}",
                    source=source, date=date, tier=1, quals=org_quals,
                    pmids=pmids)
            stats["reattributed_to_organism"] += 1
            continue

        # Every row used to become an
        # INHIBITS edge regardless of what it measured, so 1,433 of 5,486
        # pairs asserted inhibition with no active measurement behind them --
        # 26% of the direct-acting channel that M1, M7 and M8 walk. classify()
        # is the same function that decides labels, so "acts on" means one
        # thing in the features and in the labels.
        verdict = classify({"qualifiers": quals})
        if verdict == "active":
            predicate = "INHIBITS"
        elif verdict == "inactive":
            predicate = "MEASURED_INACTIVE_AGAINST"
        else:
            # Neither: a bound too weak to prove inactivity ("IC50 > 100 nM"
            # sits below the concentration a hit would be pursued at). Counted
            # rather than forced to a side, the same rule labels.py applies.
            stats["undecidable_protein_measurement"] += 1
            continue

        key = (drug, predicate, target_node, r["assay"], r["standard_type"])
        if key in seen:
            stats["duplicate"] += 1
            continue
        seen.add(key)
        em.edge(drug, predicate, target_node,
                source=source, date=date, tier=1, quals=quals, pmids=pmids)
        stats["protein_edges" if predicate == "INHIBITS" else "measured_inactive_edges"] += 1

    return stats
