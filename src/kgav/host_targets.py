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
SINGLE PROTEIN target, measured on chembl_37 at 2026-11-10: 3,022 drugs
(3,003 distinct InChIKeys), 58,817 activities, 1,646 targets. 504 of those
targets are already coronavirus dependency factors, which is what makes M2
produce real single-hop paths.

THE GATE WAS RELAXED ON 2026-11-10 AND PUT BACK THE SAME DAY, and the
reason is the one thing in this module worth reading twice. Fix-list B1
asked what the phase floor costs. Measured: relaxing it for compounds the
graph ALREADY holds -- no new nodes at all -- admits 20,512 further rows
over 1,988 compounds, 1,327 of which had no human target of any kind
before, so M2-M6 could not reach them whatever their labels said. That
took the layer to 79,329 activities / 4,990 compounds / 2,265 targets /
28,138 edges, with 778 dependency factors reached and 9,251 edges feeding
a single-hop M2 path. Dropping the floor outright instead would have
admitted 1,304,767 rows and 676,997 new compounds; `make
diagnose-phase-gate` measures all three readings in ~20s.

IT WAS DECLINED BECAUSE THE ADMISSION RULE IS DOWNSTREAM OF THE LABELS.
"Compounds the graph already has" means compounds somebody screened, and
the ones that were missing human targets are disproportionately the ones
that came back active. For DENV-2 (NCBITaxon:11060), of the compounds that
gained a first host-side route, 121 were positive and 22 negative -- 84.6%
against a 24.7% base rate, 3.43x. SARS-CoV-2 was 1.17x, so this is not
uniform and cannot be corrected by a constant. The visible effect: the
DENV-2 M3 channel went from 2 hits to 16 and p = 0.4257 to 0.000003,
clearing every correction family, while its `degree` baseline did not move
a single hit (p = 0.757 both sides) because it already reached everything.
A channel that improves only where the baseline cannot is normally
evidence the channel is doing work; here it is evidence the labels came in
with the edges. So "floor" is the default and v0.1 is built with it.

THESE THREE FIGURES READ 3,667 / 63,826 / 1,670 UNTIL 2026-11-10, and the
cause is not a mismeasurement -- re-running the query without the
`standard_type IN ('IC50','EC50','Ki','Kd')` gate returns exactly 3,667 /
63,826 / 1,670. They were right before that gate was added and nothing
re-derived them afterwards. Recorded with the date and the gates they
depend on, since the 30 GB sqlite is gitignored and no test can check them;
`make diagnose-phase-gate` re-derives them in ~20s.

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

import json
import sqlite3
from collections import Counter
from collections.abc import Iterable
from pathlib import Path

from kgav.chembl import CHEMBL_UNITS
from kgav.temporal import PLACEHOLDER_DATE
from kgav.normalize_chem import policy_key
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

#: THE PHASE GATE, AS A NAMED CHOICE (fix-list B1).
#: "floor"   -- max_phase >= DEFAULT_MIN_PHASE for everyone. How v0.1 is
#:              built, and the default: see the SCOPE note above for the
#:              label correlation that decided it.
#: "relaxed" -- the floor, OR the graph already holds the compound (option
#:              D). KEPT RUNNABLE AND NOT EVALUATED. Any figure measured
#:              under it describes a release this project does not ship,
#:              and a run must say which policy produced it.
#:              It also depends on layer order: `hosttargets` precedes
#:              `ncats`, so this policy sees 667 of the 1,327 recoverable
#:              compounds, not all of them -- see pipeline.LAYERS.
#: An unrecognised value is REFUSED rather than read as either one.
PHASE_GATE_POLICIES = ("floor", "relaxed")
DEFAULT_PHASE_GATE = "floor"


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


def known_inchikeys(releases: Path, layer_names: Iterable[str]) -> set[str]:
    """InChIKeys the already-built layers carry as SmallMolecule nodes.

    FIX-LIST B1, OPTION D. The phase gate exists to keep 679,000 bulk
    screening compounds out of the graph, and it does that job; what it also
    did was withhold human targets from compounds the graph ALREADY holds,
    because they happened to be early-phase or unphased. Those are not new
    nodes -- they are compounds that already carry antiviral labels and had
    no route into the host side, so M2-M6 could not reach them at all.

    Returned WITHOUT the `INCHIKEY:` prefix, because this set is compared
    against compound_structures.standard_inchi_key in SQL. A layer with no
    nodes.jsonl is skipped rather than refused: data/raw is populated by
    hand and an absent layer is an ordinary state of a checkout.

    **AND SINCE FIX-LIST L27 THAT COMPARISON IS NO LONGER EXACT.** Node ids
    are now the identity policy's key rather than ChEMBL's, so a compound
    whose key the policy moved -- 549 of 25,383 on v0.1 -- will not match
    its own `standard_inchi_key` in the IN-list, and the `relaxed` gate
    will quietly admit fewer compounds than it did. That is a real
    degradation and it is recorded rather than repaired here for one
    reason: `relaxed` is fix-list B1, which was DECLINED after measurement,
    and `floor` is the registered policy of
    `docs/prereg_pathway_overlap.md`. So nothing in the release depends on
    it. Anyone re-opening B1 has to map node ids back to source keys first,
    which needs the source key stored on the node and is a schema change.

    Pass `pipeline.layers_before("hosttargets")` rather than a literal list.
    """
    out: set[str] = set()
    for name in layer_names:
        f = releases / f"v0.1-{name}" / "nodes.jsonl"
        if not f.exists():
            continue
        with open(f) as fh:
            for line in fh:
                n = json.loads(line)
                if n.get("class") != "SmallMolecule":
                    continue
                nid = n.get("id", "")
                if nid.startswith("INCHIKEY:"):
                    out.add(nid[len("INCHIKEY:"):])
    return out


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
                          min_phase: int = DEFAULT_MIN_PHASE,
                          known: set[str] | None = None) -> list[dict]:
    """Human single-protein potency rows for compounds at phase >= min_phase,
    OR at any phase when `known` holds the compound's InChIKey.

    `known` comes from known_inchikeys() and is the set of compounds the
    graph already has from earlier layers (fix-list B1 option D). Omit it
    and the phase floor applies to everyone, which is how v0.1 was built.

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

    # FIX-LIST B1 OPTION D: THE GATE IS RELAXED, NOT REMOVED. With `known`
    # supplied, a row survives if its compound clears the phase floor OR the
    # graph already holds that compound. Dropping the floor outright admits
    # 1,304,767 rows and 676,997 new nodes for no gain any AUC can see; this
    # admits the 1,327 compounds that gain a first route and nothing else.
    #
    # A TEMP TABLE, NOT AN `IN (...)` LIST. The set is ~25,000 keys, which as
    # a literal list is a query SQLite has to re-parse and cannot index. The
    # temp database is separate from the main one, so this works against a
    # mode=ro connection.
    #
    # NULL-SAFE BY ACCIDENT OF OR, AND WORTH SAYING SO: for an unphased
    # compound `md.max_phase >= ?` is NULL, and NULL OR TRUE is TRUE while
    # NULL OR FALSE is NULL, which is falsy. So an unphased compound is kept
    # exactly when it is in `known`, which is the intent.
    if known:
        con.execute("CREATE TEMP TABLE _known(k TEXT PRIMARY KEY)")
        con.executemany("INSERT OR IGNORE INTO _known VALUES(?)",
                        ((k,) for k in known))
        phase_sql = ("(md.max_phase >= ? OR st.standard_inchi_key IN "
                     "(SELECT k FROM _known))")
    else:
        phase_sql = "md.max_phase >= ?"

    rows = con.execute(f"""
        SELECT md.chembl_id AS compound, md.max_phase,
               st.canonical_smiles AS smiles, st.standard_inchi_key AS inchikey,
               cs.accession,
               act.standard_type, act.standard_value, act.standard_units,
               act.standard_relation, act.pchembl_value,
               a.assay_type AS chembl_assay_type,
               a.chembl_id AS assay,
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
          AND {phase_sql}
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
        # THE SAME FIX AS chembl.py (fix-list L27). This layer read the
        # same compound_structures.standard_inchi_key and emitted it
        # verbatim, so it minted the same duplicate nodes from the other
        # side -- and it is the layer that carries TARGETS, which is the
        # first operand of the whole pathway analysis.
        inchikey, route = policy_key(r.get("smiles") or "", inchikey)
        stats[f"identity_{route}"] += 1
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
                                # THE IDENTITY POLICY, DECLARED (fix-list L24). Not a
                # per-compound measurement: normalize_chem's docstring
                # says these record WHICH KEY IS THE NODE ID, so a
                # constant is correct and both layers must state the
                # same one. salt_collapsed was False on all 25,383
                # nodes and that is wrong -- largest_fragment.choose
                # runs before the InChIKey, so `Cl.X` and `X` really
                # do share a node. stereo_collapsed stays False: the
                # id is the FULL key and L- and D-alanine differ in it.
                # tests/test_identity_policy.py derives both from what
                # normalize_chemical actually does.
                salt_collapsed=True, stereo_collapsed=False,
                tautomer_collapsed=True)

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
            # THE EARLIEST ASSERTION WINS THE DATE, because that is what the
            # field is called (fix-list P8). The edge took its date from the
            # first row by activity_id -- ChEMBL's insertion order -- which
            # is stable across rebuilds and is not the same thing as first
            # asserted. Measured on chembl_37 over 19,020 keys:
            #
            #   766  keep a year LATER than their earliest assertion,
            #        54 of them by ten years or more
            #   143  keep an UNDATED row while a dated sibling exists, so
            #        the edge gets the 1970 placeholder although the year
            #        is known
            #    71  change which side of a 2021 cutoff they fall on
            #
            # The docstring on query_host_activities already names the
            # stake -- "the second would land in a different half of a
            # temporal split" -- and ORDER BY fixed the DETERMINISM of that
            # choice without making it the right choice.
            # THE SPREAD TRAVELS WITH THE POINT (fix-list P13, option b).
            # The edge keeps the first row's value by activity_id -- a
            # deterministic, arbitrary choice -- and nothing on the edge
            # said so. Measured on chembl_37: of 19,020 keys, 7,257 have
            # more than one row, 5,517 disagree (spread median 6.67x, p90
            # 182x, max 55,625x) and 4,296 keep a non-median value; 7,116
            # of the replicate sets come from DIFFERENT assays, so they are
            # not repeats of one experiment. The point is unchanged and the
            # range is now readable beside it.
            prev["values"].append(nm)
            year = str(int(r["year"])) if r.get("year") else None
            if year and (prev["date"] == PLACEHOLDER_DATE
                         or year < prev["date"]):
                # The placeholder is compared by identity, not by value:
                # "2014" < "1970-01-01" is False lexically, so a `<` alone
                # would let the placeholder beat every real year.
                stats["date_moved_earlier" if prev["date"] != PLACEHOLDER_DATE
                      else "date_recovered_from_a_dated_duplicate"] += 1
                prev["date"] = year
            continue

        pending[key] = {
            "drug": drug, "target": target,
            # A YEAR IS WRITTEN AS A YEAR (fix-list C10).
            "date": (str(int(r["year"])) if r.get("year")
                     else PLACEHOLDER_DATE),
            "quals": {QUANT_FIELD[r["standard_type"]]: nm,
                      "relation": relation,
                      "direction": direction,
                      "assay_type": assay_type,
                      # WHICH ASSAY (fix-list V6). Already in the dedup
                      # key's sibling fields and on no edge.
                      **({"assay_id": str(r["assay"])}
                         if r.get("assay") else {})},
            "pmids": {cite} if cite else set(),
            # Every converted value for this key, in activity_id order, so
            # values[0] IS the potency the edge carries.
            "values": [nm],
            "potency_field": QUANT_FIELD[r["standard_type"]],
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
        # THE THREE SPREAD QUALIFIERS (fix-list P13, option b). Emitted on
        # EVERY edge, including the 11,763 keys with a single row, where
        # they read potency_rows=1 and min==max==the value. Emitting them only
        # for replicates would make their absence mean two things -- one
        # measurement, or a layer that predates this -- and an absent
        # qualifier is the ambiguity this project keeps paying for.
        values = rec["values"]
        # NAMED, so the range cannot be read against the wrong number.
        # 3,801 of 14,381 (drug, target) pairs carry rows of more than one
        # standard_type and merge on (subject, predicate, object), so an
        # unnamed range would describe one quantity on an edge stating two
        # -- which is fix-list A4 reintroduced under new qualifier names.
        rec["quals"]["potency_spread_field"] = rec["potency_field"]
        rec["quals"]["potency_rows"] = len(values)
        rec["quals"]["potency_min_nm"] = float(min(values))
        rec["quals"]["potency_max_nm"] = float(max(values))
        if len(values) > 1:
            stats["edges_collapsing_several_rows"] += 1
            if min(values) != max(values):
                stats["edges_whose_rows_disagree"] += 1
        em.edge(rec["drug"], "TARGETS", rec["target"],
                source=source, date=rec["date"], tier=1,
                quals=rec["quals"],
                pmids=sorted(rec["pmids"]) or None)

    return stats


def host_protein_ids(index) -> set[str]:
    return {nid for nid, n in index.nodes.items()
            if n["class"] == "Protein" and not n["properties"].get("is_viral")}
