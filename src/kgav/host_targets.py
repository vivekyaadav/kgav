"""Drug -> human protein TARGETS edges from ChEMBL.

WHY THIS LAYER EXISTS. After Day 7 the graph had 7,380 HOST_FACTOR_FOR edges
and ZERO TARGETS edges, because Day 7 only ingested activities against
coronavirus targets. That leaves five of six metapaths with nothing to
traverse:

    M2  Drug -TARGETS-> HostProtein -HOST_FACTOR_FOR-> Virus
    M3  Drug -TARGETS-> HostProtein -PHYS_INTERACTS-> ViralProtein -> Virus
    M4  Drug -TARGETS-> HostProtein -PARTICIPATES_IN-> Pathway -> ... -> Virus
    M5  Drug -TARGETS-> HostProtein -PHYS_INTERACTS-> HostProtein -> Virus
    M6  Drug -CHEM_SIMILAR-> Drug -TARGETS-> HostProtein -> Virus

Only M1 (direct-acting) worked. The entire host-directed channel -- the reason
the CRISPR and interactome layers exist -- had no entry point.

SCOPE. All compounds at max_phase >= 2 with pChEMBL >= 6 against a human
SINGLE PROTEIN target: 3,667 drugs, 63,826 activities, 1,670 targets. 504 of
those targets are already coronavirus dependency factors, which is what makes
M2 produce real single-hop paths.

Deliberately NOT restricted to compounds already in the graph. Repurposing
means proposing drugs nobody tested against the virus; limiting the pool to
compounds someone already assayed defeats the purpose.

EVERY EDGE CARRIES ITS PUBLICATIONS, as of 2026-11-07. It carried none
before: `d.year` was selected and `d.pubmed_id` was not, so all 18,954
TARGETS edges in v0.1 are unciteable -- the only predicate in the graph
with zero publication coverage, while 86.1% of the rows behind them
(50,627 of 58,817, over 10,737 distinct documents) name a PMID in the
same `docs` row the date already came from. An agent asked for evidence
on a host-directed route could show a potency and no source for it.

THE PMIDS ARE UNIONED ACROSS THE ROWS THAT COLLAPSE INTO ONE EDGE. The
key is (compound, target, assay type) and about three source rows share
one, so taking the first row's PMID would report one paper for an edge
measured in several. Everything else about the edge still comes from the
first row seen, unchanged: this commit adds publications and moves
nothing else.

DIRECTION IS 'unknown' BY DEFAULT, and that is not laziness. Inhibiting a
dependency factor should block infection; ACTIVATING one should not. ChEMBL's
bulk activity table records potency, not whether the compound is an agonist or
an antagonist -- only drug_mechanism carries action_type, and only for a
minority of approved drugs. Asserting 'inhibitor' everywhere would manufacture
a mechanistic claim the data does not make, and every M2 path would silently
inherit it.
"""
from __future__ import annotations

import sqlite3
from collections import Counter
from pathlib import Path

from kgav.chembl import CHEMBL_UNITS
from kgav.units import to_nm as _to_nm

# ChEMBL action_type -> the schema's target_direction enum.
ACTION_TO_DIRECTION = {
    "INHIBITOR": "inhibitor",
    "ANTAGONIST": "inhibitor",
    "BLOCKER": "inhibitor",
    "NEGATIVE ALLOSTERIC MODULATOR": "inhibitor",
    "NEGATIVE MODULATOR": "inhibitor",
    "PARTIAL AGONIST": "agonist",
    "AGONIST": "agonist",
    "ACTIVATOR": "agonist",
    "POSITIVE ALLOSTERIC MODULATOR": "agonist",
    "POSITIVE MODULATOR": "agonist",
    "MODULATOR": "modulator",
    "ALLOSTERIC MODULATOR": "modulator",
}

QUANT_FIELD = {"IC50": "ic50_nm", "EC50": "ec50_nm", "Ki": "ki_nm", "Kd": "kd_nm"}
VALID_RELATIONS = {"=", ">", "<", ">=", "<=", "~"}

DEFAULT_MIN_PCHEMBL = 6.0
DEFAULT_MIN_PHASE = 2


def load_mechanisms(db: Path) -> dict[tuple[str, str], str]:
    """(compound chembl_id, accession) -> action_type, from drug_mechanism.

    This is the only place ChEMBL states whether a drug inhibits or activates
    its target. It covers a minority of approved drugs; everything else stays
    'unknown'.
    """
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    rows = con.execute("""
        SELECT md.chembl_id, cs.accession, dm.action_type
        FROM drug_mechanism dm
        JOIN molecule_dictionary md ON dm.molregno = md.molregno
        JOIN target_dictionary t ON dm.tid = t.tid
        JOIN target_components tc ON t.tid = tc.tid
        JOIN component_sequences cs ON tc.component_id = cs.component_id
        WHERE dm.action_type IS NOT NULL AND cs.accession IS NOT NULL
    """).fetchall()
    con.close()
    return {(c, a): act for c, a, act in rows}


def query_host_activities(db: Path, min_pchembl: float = DEFAULT_MIN_PCHEMBL,
                          min_phase: int = DEFAULT_MIN_PHASE) -> list[dict]:
    """Human single-protein potency rows for compounds at phase >= min_phase.

    ORDER BY IS LOAD BEARING, not tidiness. About three rows share each
    (compound, target, assay type) key and the first one seen sets the
    edge's value, relation, direction and date while the rest only add
    citations. Unordered, which row that is comes from the query planner,
    so two rebuilds of the same database could give one edge two different
    potencies and two different dates -- and the second would land in a
    different half of a temporal split. activity_id is ChEMBL's own
    insertion order and is stable across its releases for existing rows.
    """
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    rows = con.execute("""
        SELECT md.chembl_id AS compound, md.max_phase,
               st.canonical_smiles AS smiles, st.standard_inchi_key AS inchikey,
               cs.accession,
               act.standard_type, act.standard_value, act.standard_units,
               act.standard_relation, act.pchembl_value,
               d.year, d.pubmed_id
        FROM activities act
        JOIN assays a            ON act.assay_id = a.assay_id
        JOIN target_dictionary t ON a.tid = t.tid
        JOIN target_components tc ON t.tid = tc.tid
        JOIN component_sequences cs ON tc.component_id = cs.component_id
        JOIN molecule_dictionary md ON act.molregno = md.molregno
        LEFT JOIN compound_structures st ON md.molregno = st.molregno
        LEFT JOIN docs d         ON a.doc_id = d.doc_id
        WHERE t.tax_id = 9606
          AND t.target_type = 'SINGLE PROTEIN'
          AND act.pchembl_value >= ?
          AND md.max_phase >= ?
          AND act.standard_type IN ('IC50','EC50','Ki','Kd')
        ORDER BY act.activity_id
    """, (min_pchembl, min_phase)).fetchall()
    con.close()
    return [dict(r) for r in rows]


def to_nm(value: float | None, units: str | None) -> float | None:
    """Nanomolar, through the one implementation in kgav.units.

    This was a verbatim-looking copy of chembl.to_nm that was not one: it
    multiplied by 0.001 for pM where chembl divided by 1000, which differs
    in the last bit for 13.75% of values. Harmless at pM scale, and the
    reason to notice is that two copies believed identical had never been
    compared -- by then their unit sets had already diverged.
    """
    return _to_nm(value, units, allow=CHEMBL_UNITS)


def ingest_host_targets(em, rows: list[dict], mechanisms: dict[tuple[str, str], str],
                        host_proteins: set[str], source: str) -> Counter:
    """Emit TARGETS edges, keeping only targets the graph already contains.

    A target absent from the host layer would dangle after assembly, and an
    activity against a protein the graph does not model contributes nothing to
    any metapath -- so it is counted and dropped rather than minting a node
    this pass has no proteome data for.
    """
    stats: Counter = Counter()
    # key -> the edge's arguments, with publications accumulated. A set
    # rather than a counter: two rows citing the same paper are one piece
    # of evidence, and sorted() below makes the output order stable so a
    # rebuild that changes nothing produces an identical edges.jsonl.
    pending: dict[tuple, dict] = {}

    for r in rows:
        stats["rows"] += 1
        inchikey, acc = r.get("inchikey"), r.get("accession")
        if not inchikey:
            stats["no_inchikey"] += 1
            continue
        target = f"UniProtKB:{acc}"
        if target not in host_proteins:
            stats["target_not_in_graph"] += 1
            continue

        nm = to_nm(r["standard_value"], r["standard_units"])
        if nm is None:
            stats["unconvertible_units"] += 1
            continue
        relation = r["standard_relation"] if r["standard_relation"] in VALID_RELATIONS else "="
        if relation != "=":
            # pChEMBL >= 6 already excludes most censored rows, but a '<'
            # relation on a potent value is a real lower bound worth keeping.
            stats["censored"] += 1

        drug = f"INCHIKEY:{inchikey}"
        phase = r.get("max_phase")
        em.node(drug, "SmallMolecule",
                smiles=r.get("smiles") or "",
                inchikey_skel=inchikey[:14],
                chembl_id=r["compound"],
                max_phase=int(phase) if phase is not None else None,
                is_approved=bool(phase == 4),
                # Must match the Day 7 chembl layer exactly, or the assembler
                # reports property conflicts on thousands of shared compounds.
                salt_collapsed=False, stereo_collapsed=False)

        action = mechanisms.get((r["compound"], acc))
        direction = ACTION_TO_DIRECTION.get((action or "").upper(), "unknown")

        pmid = r.get("pubmed_id")
        cite = f"PMID:{int(pmid)}" if pmid else None

        key = (drug, target, r["standard_type"])
        prev = pending.get(key)
        if prev is not None:
            # NOT `continue`. A duplicate row is a second measurement of the
            # same pair and it has its own paper; the edge keeps the first
            # row's value, relation, direction and date, and gains the
            # citation. Dropping the row outright is what left edges
            # measured in four papers citing one.
            stats["duplicate"] += 1
            if cite and cite not in prev["pmids"]:
                prev["pmids"].add(cite)
                stats["duplicate_added_publication"] += 1
            continue

        pending[key] = {
            "drug": drug, "target": target,
            "date": (f"{int(r['year'])}-01-01" if r.get("year")
                     else "1970-01-01"),
            "quals": {QUANT_FIELD[r["standard_type"]]: nm,
                      "relation": relation,
                      "direction": direction,
                      "assay_type": "binding"},
            "pmids": {cite} if cite else set(),
        }
        stats["edges"] += 1

    for rec in pending.values():
        # PER EDGE, NOT PER ROW. Counted in the row loop until 2026-11-07,
        # where three rows that collapse into one edge each added one --
        # so the driver printed a direction census summing to the ~58,000
        # source rows directly beneath "18,954 TARGETS edges", and
        # test_direction_defaults_to_unknown_not_inhibitor only held
        # because its fixture had no duplicates to expose it.
        stats[f"direction_{rec['quals']['direction']}"] += 1
        if rec["pmids"]:
            stats["with_publication"] += 1
        else:
            stats["no_publication"] += 1
        em.edge(rec["drug"], "TARGETS", rec["target"],
                source=source, date=rec["date"], tier=1,
                quals=rec["quals"],
                pmids=sorted(rec["pmids"]) or None)

    return stats


def host_protein_ids(index) -> set[str]:
    return {nid for nid, n in index.nodes.items()
            if n["class"] == "Protein" and not n["properties"].get("is_viral")}
