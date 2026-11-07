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

assay_type IS READ, NOT ASSERTED, as of 2026-11-07. All 18,954 edges in
v0.1 say `assay_type: "binding"` -- a constant at the emit site -- while
ChEMBL classifies the rows behind them as 52,485 binding, 4,936
functional, 1,301 ADMET and 95 toxicity. See CHEMBL_ASSAY_TYPE below for
the mapping and for the 1,301 ADMET rows, which are CYP and hERG
liability panels that have been feeding M2-M6 as target engagement.

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
from kgav.temporal import PLACEHOLDER_DATE
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

# ChEMBL's OWN assay classification -> the schema's assay_type enum.
#
# All 18,954 TARGETS edges in v0.1 say `assay_type: "binding"`, a constant
# written at the emit site. ChEMBL classifies the 58,817 rows behind them
# as 52,485 B, 4,936 F, 1,301 A and 95 T, so the constant is wrong on
# 6,332 of them. src/kgav/chembl.py found the identical bug in the
# antiviral layer -- "the ingest wrote 'biochemical' on every protein row
# regardless" -- and the fix there was to stop asserting and start
# reading. Same fix, same reason.
#
# NOT DERIVED FROM standard_type, which is what the fix list proposed:
# 1,868 EC50 rows are B-coded and 1,426 IC50 rows are F-coded, so the
# readout name does not determine the assay kind and a mapping from it
# would be wrong in both directions at once.
#
# F -> biochemical IS THE ONE JUDGEMENT HERE. ChEMBL's F is a functional
# (activity) assay, cell-based or enzymatic; `biochemical` is this
# schema's only in-vitro activity member, because the vocabulary was
# written around antiviral assays and its cellular members all say
# "antiviral", which a CYP3A4 measurement is not. `other` would be the
# alternative and would throw away the 4,936/52,485 split. The gap is
# real and is recorded in the fix list rather than papered over by
# picking a wrong-but-specific value.
#
# A -> other is NOT a tidy-up, it is a finding. 1,301 of these rows are
# ADMET liability panels: CYP3A4 (247), CYP2D6 (125), CYP2C9 (92), hERG
# (64 as T), CYP2C19, CYP2J2, CYP2B6, aldehyde oxidase. Every one became
# a TARGETS edge indistinguishable from therapeutic target engagement and
# fed M2-M6. They are counted and named at ingest now, not dropped:
# inhibiting CYP3A4 at 100 nM is a true statement about the compound, and
# whether a repurposing route may walk through it is a judgement for the
# person reading the census, not for this mapping.
CHEMBL_ASSAY_TYPE = {
    "B": "binding",
    "F": "biochemical",
    "T": "cytotoxicity",
    "A": "other",
}

# ChEMBL's letter codes, for the census the driver prints.
CHEMBL_ASSAY_LABEL = {
    "B": "binding", "F": "functional", "A": "ADMET", "T": "toxicity",
}
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


def gate_census(db: Path, min_pchembl: float = DEFAULT_MIN_PCHEMBL,
                min_phase: int = DEFAULT_MIN_PHASE) -> Counter:
    """What query_host_activities' WHERE clause discards, by reason.

    FOUR GATES IN ONE WHERE CLAUSE, none countable from inside it
    (fix-list E7): non-human target, target not a SINGLE PROTEIN, pChEMBL
    below the floor, and max_phase below the floor. The last is fix-list
    B1's subject -- a proposal to drop it entirely -- and nothing in this
    repository could say what it costs.

    The gates are counted in the order the WHERE clause applies them
    conceptually, so each number is "rows this gate was the FIRST to
    reject". A row failing two gates is attributed to one of them, which
    is what makes the four add up to the total.

    MEASURED ON chembl_37, 2026-11-07, and recorded here because this is a
    full scan of 24.5 million activity rows and takes 76 seconds against
    an 8-second ingest -- so the driver runs it behind --gate-census
    rather than on every rebuild:

        rows scanned                       24,527,044
        dropped: not human                 10,573,458
        dropped: target not SINGLE PROTEIN  5,654,400
        dropped: standard_type not a potency 5,530,961
        dropped: pchembl below 6            1,404,641
        dropped: max_phase below 2          1,304,767   <- fix-list B1
        kept                                   58,817

    `kept` is exactly what query_host_activities returns, and the six
    numbers sum to the scan, so the census and the WHERE clause cannot
    disagree about what the gates mean without this failing to add up.
    """
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    out: Counter = Counter()
    rows = con.execute("""
        SELECT t.tax_id = 9606                        AS human,
               t.target_type = 'SINGLE PROTEIN'       AS single,
               act.standard_type IN ('IC50','EC50','Ki','Kd') AS potency,
               act.pchembl_value IS NOT NULL AND act.pchembl_value >= ?
                                                      AS pchembl_ok,
               md.max_phase IS NOT NULL AND md.max_phase >= ?
                                                      AS phase_ok,
               COUNT(*) AS n
        FROM activities act
        JOIN assays a            ON act.assay_id = a.assay_id
        JOIN target_dictionary t ON a.tid = t.tid
        JOIN molecule_dictionary md ON act.molregno = md.molregno
        GROUP BY 1, 2, 3, 4, 5
    """, (min_pchembl, min_phase)).fetchall()
    con.close()
    for human, single, potency, pchembl_ok, phase_ok, n in rows:
        out["rows_scanned"] += n
        if not human:
            out["dropped_not_human"] += n
        elif not single:
            out["dropped_target_not_single_protein"] += n
        elif not potency:
            out["dropped_standard_type_not_a_potency"] += n
        elif not pchembl_ok:
            out[f"dropped_pchembl_below_{min_pchembl:g}"] += n
        elif not phase_ok:
            # FIX-LIST B1'S NUMBER. The proposal is to drop this filter;
            # this is what it currently withholds.
            out[f"dropped_max_phase_below_{min_phase:g}"] += n
        else:
            out["kept"] += n
    return out


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
               a.assay_type AS chembl_assay_type,
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
        # ROW-level, named so: the per-EDGE census is in the emit loop. The
        # direction counters conflated the two until 2026-11-07 and printed
        # a census summing to the source rows beneath an edge count.
        #
        # Counted BEFORE any gate, so this census covers exactly the
        # population the ADMET measurement was taken over and
        # sum(row_chembl_assay_type_*) == rows is an assertable partition.
        code = (r.get("chembl_assay_type") or "").upper()
        assay_type = CHEMBL_ASSAY_TYPE.get(code, "other")
        stats[f"row_chembl_assay_type_{code or 'none'}"] += 1
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
        # THE COERCION IS COUNTED, which it was not. An absent or unexpected
        # standard_relation became "=" silently and was then indistinguishable
        # from a value ChEMBL actually stated as exact. Measured on chembl_37:
        # ZERO of the 58,817 rows need it -- every one states "=" -- so the
        # count is 0 today and the branch is a guard against a future release,
        # not a description of this one. Worth the two lines precisely because
        # a silent 0 and a silent 40,000 looked the same from here.
        given = r["standard_relation"]
        if given in VALID_RELATIONS:
            relation = given
        else:
            relation = "="
            stats["relation_imputed"] += 1
        if relation != "=":
            # pChEMBL >= 6 already excludes most censored rows, but a '<'
            # relation on a potent value is a real lower bound worth keeping.
            # Also zero on chembl_37, for the same reason.
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
            # A YEAR IS WRITTEN AS A YEAR (fix-list C10).
            "date": (str(int(r["year"])) if r.get("year")
                     else PLACEHOLDER_DATE),
            "quals": {QUANT_FIELD[r["standard_type"]]: nm,
                      "relation": relation,
                      "direction": direction,
                      "assay_type": assay_type},
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
        stats[f"assay_type_{rec['quals']['assay_type']}"] += 1
        # THE PLACEHOLDER, COUNTED (fix-list C10). "undated" is not one
        # fact: 14,888 edges across seven sources carry it, and each is a
        # different gap with a different remedy.
        stats["undated_placeholder" if rec["date"] == PLACEHOLDER_DATE
              else "dated_year_only"] += 1
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
