"""Drug -> Disease edges from ChEMBL's drug_indication table.

THE TWO PREDICATES THIS FILLS WERE DECLARED AND EMPTY. The schema marks
`TREATS` and `IN_TRIAL_FOR` `evaluation_label: true`, and
doc2/kg_release_census.txt printed the consequence on every release:

    2 of 3 evaluation label predicate(s) carry no edges, so every
    evaluation rests on HAS_ANTIVIRAL_ACTIVITY_AGAINST alone.

So every AUC in this repository was one predicate derived from in-vitro
concentration-response measurements at one threshold. A scorer that
recovered the idiosyncrasies of that screening and nothing else would
pass all of them. This layer adds a label of a different KIND -- did a
compound reach a clinical trial, or an approval, for the disease this
virus causes -- registered in advance by doc2/registrations/prereg_clinical_labels.md.

NO NEW COMPOUNDS, AND THAT IS A REGISTERED REFUSAL rather than a
convenience. The registration forbids this layer from adding a
SmallMolecule: a new node would change every other population in the
repository -- the scorer's candidate sets, the metapath reach, the
degree baseline -- and would make `make auc-delta` uninterpretable for
the release that introduced the label it is supposed to be an
independent check on. A compound ChEMBL indicates and this graph does
not already hold is DROPPED and counted.

THE JOIN IS DERIVED, NOT TYPED. Register MONDO id -> mondo.json xrefs ->
MeSH / EFO / DOID -> drug_indication. An earlier hand-written version of
this mapping put the COVID-19 MeSH term `D000086382` under BOTH
MONDO:0100096 and MONDO:0100116, because MERS has no MeSH xref and the
nearest-looking term was reached for; that attributed 66 COVID-19 trials
to MERS. Deriving the xrefs from the ontology cannot make that mistake,
and `ambiguous_xref` below refuses outright if two of the registered
diseases ever claim one identifier.

WHAT CHEMBL DOES NOT CARRY, stated because the schema requires both
fields and silence would look like data:

  status      `IN_TRIAL_FOR.status` is required and its enum has
              `unknown` for exactly this. ChEMBL records that an
              indication was studied at a phase, not whether the trial
              recruited, completed or was terminated. Every edge here
              carries `unknown`, so a consumer filtering for
              `completed` gets nothing rather than something invented.

  dates       there is no usable per-row date, which is why the
              registration states that NO TEMPORAL SPLIT is available
              and that a positive result cannot license a predictive
              claim. first_asserted_date is the epoch placeholder this
              project uses for a source that states none, as
              ingest_spine does for CAUSES.

  phase 4 IS NOT "WORKS".  max_phase_for_ind == 4 means approved for
              THIS indication, so it maps to TREATS with
              approval_status `approved`. Anything below it -- 0.5
              through 3 -- is IN_TRIAL_FOR. 200 of the 222 COVID-19
              in-trial compounds this layer can place were measured
              INACTIVE in vitro against SARS-CoV-2, which the
              registration counts before any ranking exists: the two
              axes are close to orthogonal and that is the point of
              adding this one.
"""
from __future__ import annotations

import json
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path

#: ChEMBL is the source of the indication; the trial registry it cites is
#: recorded per edge as nct_id. There is no separate infores for this
#: table.
INDICATION_SOURCE = "infores:chembl"

#: Curated inferred, not curated experimental: ChEMBL's curators read a
#: registry entry or a label, which is a record of an intention to treat
#: and not a measurement.
EVIDENCE_TIER = 2

#: max_phase_for_ind at or above this is an approval for the indication.
TREATS_PHASE = 4.0

#: No date in the source. The same placeholder ingest_spine uses for
#: CAUSES, so "1970-01-01" means "the source states none" consistently
#: across the graph rather than meaning it here and something else there.
NO_DATE = "1970-01-01"

#: `drug_indication.mesh_id` holds a BARE MeSH identifier (D045169), so
#: the MONDO xref prefix is stripped. `drug_indication.efo_id` holds a
#: PREFIXED one, and ChEMBL files them under several ontologies.
MESH_PREFIX = "MESH:"
EFO_PREFIXES = ("EFO:", "DOID:", "MONDO:", "Orphanet:", "HP:")


class IndicationError(RuntimeError):
    """The join cannot be performed as registered."""


def disease_xrefs(viruses: list[dict], mondo_path: Path
                  ) -> dict[str, tuple[set[str], set[str]]]:
    """{MONDO curie: (bare MeSH ids, prefixed EFO-column ids)}.

    Read from the MONDO release, for the register's own distinct ids and
    no others. The MONDO curie itself is included among the prefixed
    ids because ChEMBL files some rows under MONDO directly.

    REFUSES on an identifier claimed by two of the registered diseases.
    Attribution has to be unique or a trial lands on the wrong virus,
    which is the defect this function was written to make impossible.
    """
    wanted = {v["mondo"] for v in viruses if v.get("mondo")}
    if not wanted:
        raise IndicationError(
            "no register entry carries a MONDO id, so there is no disease "
            "to attach an indication to")
    graph = json.loads(Path(mondo_path).read_text())["graphs"][0]
    out: dict[str, tuple[set[str], set[str]]] = {
        curie: (set(), {curie}) for curie in sorted(wanted)}
    for node in graph["nodes"]:
        curie = node["id"].replace(
            "http://purl.obolibrary.org/obo/MONDO_", "MONDO:")
        if curie not in out:
            continue
        mesh, efo = out[curie]
        for xref in (node.get("meta") or {}).get("xrefs") or ():
            val = xref["val"]
            if val.startswith(MESH_PREFIX):
                mesh.add(val[len(MESH_PREFIX):])
            elif val.startswith(EFO_PREFIXES):
                efo.add(val)
    missing = [c for c, (m, e) in out.items() if not m and len(e) == 1]
    if missing:
        # Not fatal: MERS-CoV's term carries a DOID and no MeSH, and the
        # MONDO curie alone is a legitimate key. Recorded by the caller.
        pass
    claimed: dict[str, list[str]] = defaultdict(list)
    for curie, (mesh, efo) in out.items():
        for ident in [f"{MESH_PREFIX}{m}" for m in mesh] + sorted(efo):
            if ident == curie:
                continue
            claimed[ident].append(curie)
    ambiguous = {k: v for k, v in claimed.items() if len(v) > 1}
    if ambiguous:
        raise IndicationError(
            "an identifier is claimed by more than one registered disease, "
            f"so an indication cannot be attributed: {ambiguous}")
    return out


def _rows(db: Path, mesh: set[str], efo: set[str]) -> list[tuple]:
    """(inchikey, max_phase_for_ind, nct_id or None) for one disease.

    LEFT JOIN on the ClinicalTrials refs: a row with no registry
    reference still comes back, so the caller can count it as a drop
    rather than never see it. An inner join would have hidden it.
    """
    clauses, args = [], []
    if mesh:
        clauses.append(f"di.mesh_id in ({','.join('?' * len(mesh))})")
        args += sorted(mesh)
    if efo:
        clauses.append(f"di.efo_id in ({','.join('?' * len(efo))})")
        args += sorted(efo)
    if not clauses:
        return []
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    try:
        return con.execute(
            "select cs.standard_inchi_key, di.max_phase_for_ind, ir.ref_id "
            "from drug_indication di "
            "join compound_structures cs using(molregno) "
            "left join indication_refs ir "
            "  on ir.drugind_id = di.drugind_id "
            " and ir.ref_type = 'ClinicalTrials' "
            f"where {' or '.join(clauses)}", args).fetchall()
    finally:
        con.close()


def ingest_indications(db: Path, xrefs: dict[str, tuple[set[str], set[str]]],
                       compounds: set[str], diseases: set[str],
                       em, stats: Counter | None = None) -> Counter:
    """Emit TREATS / IN_TRIAL_FOR for every placeable indication.

    ONE EDGE PER (COMPOUND, DISEASE), which is what the registration
    counted. `nct_id` is not an identity qualifier, so several trial
    registrations for one pair are ONE fact with `support_count`
    witnesses; the representative id is the lowest, chosen so a rebuild
    produces the same edge rather than whichever row the database
    returned first.
    """
    s = stats if stats is not None else Counter()
    for mondo in sorted(xrefs):
        mesh, efo = xrefs[mondo]
        if mondo not in diseases:
            s["disease_not_in_graph"] += 1
            s[f"disease_absent:{mondo}"] += 1
            continue
        pairs: dict[str, dict] = defaultdict(
            lambda: {"phase": set(), "nct": set()})
        for key, phase, nct in _rows(db, mesh, efo):
            s["rows"] += 1
            if key is None:
                s["row_without_structure"] += 1
                continue
            nid = f"INCHIKEY:{key}"
            if nid not in compounds:
                s["compound_not_in_graph"] += 1
                continue
            if phase is None:
                s["row_without_phase"] += 1
                continue
            pair = pairs[nid]
            pair["phase"].add(float(phase))
            if nct:
                pair["nct"].add(nct)
        for nid, pair in sorted(pairs.items()):
            phase = max(pair["phase"])
            if len(pair["phase"]) > 1:
                # ChEMBL files one compound-disease pair at several
                # phases when distinct trials reached distinct stages.
                # The furthest is the fact about the pair; the others
                # are not contradicted, they are earlier.
                s["pair_with_several_phases"] += 1
            ncts = sorted(pair["nct"])
            if phase >= TREATS_PHASE:
                em.edge(nid, "TREATS", mondo,
                        source=INDICATION_SOURCE, date=NO_DATE,
                        tier=EVIDENCE_TIER,
                        quals={"approval_status": "approved"},
                        support=len(ncts) or None)
                s["treats"] += 1
                s[f"treats:{mondo}"] += 1
                continue
            if not ncts:
                # REQUIRED QUALIFIER, SO THE EDGE CANNOT SHIP. The
                # schema requires nct_id on IN_TRIAL_FOR and ChEMBL
                # cites no ClinicalTrials reference for this row --
                # DailyMed or an ATC code instead. Dropping it is the
                # only honest option: inventing an id or relaxing the
                # schema would both make an unverifiable trial edge
                # indistinguishable from a verifiable one.
                s["in_trial_without_nct"] += 1
                s[f"no_nct:{mondo}"] += 1
                continue
            em.edge(nid, "IN_TRIAL_FOR", mondo,
                    source=INDICATION_SOURCE, date=NO_DATE,
                    tier=EVIDENCE_TIER,
                    quals={"nct_id": ncts[0],
                           "phase": _phase_label(phase),
                           "status": "unknown"},
                    support=len(ncts))
            s["in_trial"] += 1
            s[f"in_trial:{mondo}"] += 1
            if len(ncts) > 1:
                s["pair_with_several_trials"] += 1
    s["edges"] = s["treats"] + s["in_trial"]
    return s


def _phase_label(phase: float) -> str:
    """ChEMBL's max_phase_for_ind as a string, 0.5 kept as 0.5.

    `phase` is a free string in the schema and ChEMBL uses 0.5 for
    early-phase studies, so rounding it to 0 or 1 would state something
    the source does not. 2.0 prints as "2".
    """
    return str(int(phase)) if float(phase).is_integer() else str(phase)
