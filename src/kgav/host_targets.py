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
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    rows = con.execute("""
        SELECT md.chembl_id AS compound, md.max_phase,
               st.canonical_smiles AS smiles, st.standard_inchi_key AS inchikey,
               cs.accession,
               act.standard_type, act.standard_value, act.standard_units,
               act.standard_relation, act.pchembl_value,
               d.year
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
    """, (min_pchembl, min_phase)).fetchall()
    con.close()
    return [dict(r) for r in rows]


def to_nm(value: float | None, units: str | None) -> float | None:
    if value is None:
        return None
    factors = {"nM": 1.0, "uM": 1000.0, "pM": 0.001, "mM": 1_000_000.0}
    f = factors.get((units or "").strip())
    return float(value) * f if f else None


def ingest_host_targets(em, rows: list[dict], mechanisms: dict[tuple[str, str], str],
                        host_proteins: set[str], source: str) -> Counter:
    """Emit TARGETS edges, keeping only targets the graph already contains.

    A target absent from the host layer would dangle after assembly, and an
    activity against a protein the graph does not model contributes nothing to
    any metapath -- so it is counted and dropped rather than minting a node
    this pass has no proteome data for.
    """
    stats: Counter = Counter()
    seen: set[tuple] = set()

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
        if direction != "unknown":
            stats[f"direction_{direction}"] += 1
        else:
            stats["direction_unknown"] += 1

        key = (drug, target, r["standard_type"])
        if key in seen:
            stats["duplicate"] += 1
            continue
        seen.add(key)

        em.edge(drug, "TARGETS", target,
                source=source,
                date=f"{int(r['year'])}-01-01" if r.get("year") else "1970-01-01",
                tier=1,
                quals={QUANT_FIELD[r["standard_type"]]: nm,
                       "relation": relation,
                       "direction": direction,
                       "assay_type": "binding"})
        stats["edges"] += 1

    return stats


def host_protein_ids(index) -> set[str]:
    return {nid for nid, n in index.nodes.items()
            if n["class"] == "Protein" and not n["properties"].get("is_viral")}
