"""Drug -> human protein TARGETS from the IUPHAR/BPS Guide to PHARMACOLOGY.

WHY THIS SOURCE AND NOT ANOTHER. Two gaps in the graph are measured, and
this closes the one that matters more.

  DIRECTION. Before this layer landed, 12,155 of 14,381 TARGETS edges
  stated no direction -- 84.5% -- because ChEMBL's bulk activity table
  records potency and only drug_mechanism carries action_type, for a
  minority of approved drugs. So "this drug targets a dependency factor"
  could not become "inhibiting it should block infection". GtP curates
  `action` on essentially every interaction, which is the field that gap
  is made of.

  IT WORKED, AND THE FIGURE ABOVE IS NOT v0.1's. On the assembled v0.1,
  which CONTAINS this layer, 10,175 of 17,587 TARGETS edges state no
  direction -- 57.8%, down from 84.5%. The directed count went 2,226 ->
  7,587, a 3.4x rise, by two routes: new directed edges, and existing
  ChEMBL edges gaining an action through the qualifier ladder. The
  pre-gtop numbers were quoted as "v0.1's" in seven places across this
  module, assemble.py, guards.py, ingest_gtop.py and two documents --
  correct numbers attributed to the wrong release, which is harder to
  notice than a wrong number.

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
import re
from collections import Counter
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path

from kgav.assemble import POTENCY_FIELD_ORDER, resolve_qualifier
from kgav.normalize_bio import UniProtResolver

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
C_CHEMBL_ID = "ChEMBL ID"
C_AFF_UNITS = "Original Affinity Units"
C_AFF_LOW = "Original Affinity Low nm"
C_AFF_HIGH = "Original Affinity High nm"

#: GtoP's affinity unit -> the TARGETS qualifier the SCHEMA declares for it.
#:
#: There is no untyped `affinity` field and there should not be: TARGETS
#: declares kd_nm, ic50_nm, ki_nm and ec50_nm separately because which
#: measurement a number is forms part of the claim. This layer invented
#: `affinity_nm` instead, which nothing noticed until the layer first
#: produced edges -- 4,270 UNDECLARED_QUALIFIER violations, 8 seconds into
#: a phase-4 assembly, 9 minutes into a rebuild.
#:
#: AN UNRECOGNISED UNIT EMITS NO POTENCY AT ALL, and is counted. This table
#: was written against GtoP's documented vocabulary rather than against a
#: measured one, which is the mistake that produced this module's whole
#: history -- so it fails to silence rather than to a guess. Extend it from
#: `affinity_unit_unmapped` in the driver's report, not from documentation.
AFFINITY_UNIT_TO_QUALIFIER = {
    "PKI": "ki_nm", "KI": "ki_nm",
    "PIC50": "ic50_nm", "IC50": "ic50_nm",
    "PKD": "kd_nm", "KD": "kd_nm",
    "PEC50": "ec50_nm", "EC50": "ec50_nm",
}

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

#: The InChIKey's first block: the connectivity skeleton, without the
#: stereo/salt/protonation blocks that follow it.
SKELETON_LEN = 14


@dataclass(frozen=True)
class CompoundIndex:
    """How a GtoP ligand finds its node in the graph, and by which route.

    THE EXACT KEY IS NOT THE ONLY HONEST JOIN. The chemistry layer salt-
    and stereo-collapses InChIKeys, so a GtoP ligand can BE a compound the
    graph already holds, under a different form of the same key. Measured
    on 2026.3.1 against v0.1: of 9,827 human interactions whose exact key
    is absent, 661 resolve by the 14-character skeleton and 24 more by
    ChEMBL id.

    AND THE ALTERNATIVE IS NOT MINTING. A new node keyed on a raw GtoP
    InChIKey would sit beside its own collapsed duplicate, splitting one
    chemical's degree, its DWPC and its labels across two ids. The
    remaining 9,142 are genuinely absent compounds, and adding them is a
    different decision with a different cost -- 9,142 unlabelled molecules
    enlarge every `pool` in `lift = hits * pool / (k * positives)` while
    adding nothing to any numerator.

    RESOLUTION IS ORDERED, most exact first, and the route is returned so
    the driver can report how many edges each one earned. A join nobody
    can see the provenance of is a join nobody can audit.
    """

    exact: frozenset[str]
    by_skeleton: Mapping[str, str] = field(default_factory=dict)
    by_chembl: Mapping[str, str] = field(default_factory=dict)

    def resolve(self, inchikey: str,
                chembl_id: str = "") -> tuple[str | None, str]:
        """(compound id, route) for one ligand, or (None, "absent")."""
        key = (inchikey or "").strip()
        if not key:
            return None, "absent"
        cid = f"INCHIKEY:{key}"
        if cid in self.exact:
            return cid, "exact"
        alt = self.by_skeleton.get(key[:SKELETON_LEN])
        if alt:
            return alt, "skeleton"
        alt = self.by_chembl.get((chembl_id or "").strip()) if chembl_id else None
        if alt:
            return alt, "chembl_id"
        return None, "absent"

    @classmethod
    def from_nodes(cls, nodes: Iterable[dict]) -> CompoundIndex:
        """Built from a release's SmallMolecule nodes.

        A skeleton shared by two nodes is NOT indexed. The graph holds
        24,051 skeletons over 25,383 compounds, so some are shared -- and
        a skeleton that names two compounds names neither. Resolving it to
        whichever node was read first would attach an interaction to an
        arbitrary one of two molecules.
        """
        exact: set[str] = set()
        skel: dict[str, str | None] = {}
        chembl: dict[str, str | None] = {}
        for n in nodes:
            if n.get("class") != "SmallMolecule":
                continue
            exact.add(n["id"])
            props = n.get("properties") or {}
            k = props.get("inchikey_skel")
            if k:
                skel[k] = None if k in skel else n["id"]
            c = props.get("chembl_id")
            if c:
                c = str(c)
                chembl[c] = None if c in chembl else n["id"]
        return cls(frozenset(exact),
                   {k: v for k, v in skel.items() if v},
                   {k: v for k, v in chembl.items() if v})


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


# "GtoPdb Version: 2026.3.1 - published: 2026-10-12" -> the two facts in it.
# Both groups optional and matched independently: a future download that
# drops the date must still yield the release rather than nothing.
_RELEASE_RE = re.compile(r"Version:\s*([0-9][0-9.]*)", re.I)
_PUBLISHED_RE = re.compile(r"published:\s*(\d{4}-\d{2}-\d{2})", re.I)


def parse_version(decl: str | None) -> dict[str, str | None]:
    """The declaration split into {raw, release, published}.

    THE PUBLISHED DATE IS THE ONE FIGURE OF ITS KIND IN THIS PROJECT.
    Every other source's "release date" in the repo is a download mtime --
    `string_release_date` is literally a file timestamp. GtoPdb states its
    own publication date in the file, so this is the real thing, and it is
    kept as a separate field rather than left inside a sentence where only
    a human reading the MANIFEST could use it.

    NOTHING IS GUESSED. A declaration that states a release and no date
    gives `published: None`; a file with no declaration at all gives all
    three None. A fabricated date in a provenance record is worse than an
    absent one, because a reader cannot tell it from a real one.
    """
    if not decl:
        return {"raw": None, "release": None, "published": None}
    rel = _RELEASE_RE.search(decl)
    pub = _PUBLISHED_RE.search(decl)
    return {"raw": decl,
            "release": rel.group(1) if rel else None,
            "published": pub.group(1) if pub else None}


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


def ligand_chembl_ids(path: Path) -> dict[str, str]:
    """{GtP ligand_id: ChEMBL id}, for ligands the InChIKey join misses.

    A second route, not a replacement: 2,298 GtoP ligands carry no
    InChIKey at all -- peptides, antibodies, biologics -- and 44 of those
    do carry a ChEMBL id. The graph holds one on 20,149 of its 25,383
    compounds, so the two spaces overlap enough to be worth asking.
    """
    out = {}
    for r in _rows(path, LIGAND_COLUMNS):
        lid = str(r.get(C_LIGAND_ID) or "").strip()
        cid = str(r.get(C_CHEMBL_ID) or "").strip()
        if lid and cid:
            out[lid] = cid
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


def _attach_spread(quals: dict, values: dict[str, list[float]],
                   stats: Counter) -> None:
    """State the spread behind this edge's potency (fix-list P14).

    WHY IT IS COUNTED AND NOT READ OFF THE SOURCE. GtoPdb exports
    `Original Affinity Low nm` and `Original Affinity High nm` beside the
    median, so the obvious implementation is to carry the source's own
    interval. Measured on 2026.3.1: of the 14,066 rows that state a usable
    median, **0** state a low or a high. The 940 rows that do state one
    state no median (counted as `affinity_range_without_median`), so there
    is no row anywhere in this export where the two coexist. The range
    therefore comes from collapsing the duplicate rows of one pair -- the
    same thing host_targets counts -- which is what makes `potency_rows` mean
    ONE thing across both layers: source records collapsed into the named
    potency, with n=1 implying min == max == the point.

    THE POINT IS THE BOTTOM OF THE RANGE HERE, because the accumulation
    above keeps the most potent value per type. host_targets' point is its
    first row's and sits anywhere inside. Both satisfy the only invariant
    that matters downstream and that assemble enforces: the range contains
    the qualifier it names.

    ONE RANGE SLOT, SEVERAL MEASUREMENTS. A pair can carry a Ki and an
    IC50; the edge can describe the spread of one. The field with the most
    collapsed rows is named, since a field with one row has no spread to
    report, and ties go to POTENCY_FIELD_ORDER so the choice cannot depend
    on dict order. What the other field's spread would have said is lost,
    and `spread_not_recorded_for_another_measured_field` says how often.
    """
    if not values:
        return
    named = max(values, key=lambda q: (len(values[q]),
                                       -POTENCY_FIELD_ORDER.index(q)))
    rows = values[named]
    quals["potency_spread_field"] = named
    quals["potency_rows"] = len(rows)
    quals["potency_min_nm"] = float(min(rows))
    quals["potency_max_nm"] = float(max(rows))
    if len(rows) > 1:
        stats["edges_collapsing_several_rows"] += 1
        if min(rows) != max(rows):
            stats["edges_whose_rows_disagree"] += 1
    if len(values) > 1:
        stats["edge_measures_several_fields_one_spread_kept"] += 1
        for q, seen in values.items():
            if q != named and len(seen) > 1 and min(seen) != max(seen):
                stats["spread_not_recorded_for_another_measured_field"] += 1


def ingest_interactions(em, rows: list[dict], inchikeys: dict[str, str],
                        known_proteins: set[str] | None = None,
                        compounds: CompoundIndex | None = None,
                        ligand_chembl: Mapping[str, str] | None = None,
                        default_date: str = "1970-01-01",
                        accessions: UniProtResolver | None = None) -> Counter:
    """Emit TARGETS edges for human interactions whose ligand has a key.

    `compounds`, when given, restricts to molecules already in the graph
    and records which join route found each one. Minting new ones is
    possible and deliberately not the default: this layer exists to add
    DIRECTION and reach to compounds the graph HAS, and a new node keyed
    on an InChIKey the chemistry layer never normalised would sit beside
    its own duplicate. Passing None mints.
    """
    stats: Counter = Counter()
    pending: dict[tuple[str, str], dict] = {}
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
            # WHY IT DID NOT JOIN, NOT JUST THAT IT DID NOT (fix-list P3).
            # One counter stood here, and it could not tell "this protein is
            # not in our proteome" from "UniProt renamed it and we never
            # noticed". normalize_bio.UniProtResolver returns exactly that
            # distinction -- primary / secondary / demerged / deleted /
            # malformed -- and was reached by nothing but its own test,
            # while data/raw/refs/sec_ac.txt (1,159,114 accessions) and
            # delac_sp.txt (1,912) were downloaded on every build and read
            # by no code at all.
            #
            # AND THE RECOVERY IS COUNTED PER ROUTE, which is the rule the
            # compound join below already states: "a join that recovers
            # edges by a looser key has to be auditable". A secondary
            # accession whose primary IS in the graph is a row this layer
            # was dropping for a reason that is not about the protein.
            #
            # ON THE RELEASE OF 2026-11-09 THIS RECOVERS NOTHING: all
            # 19,313 human single-accession rows join directly, and 0 of
            # the graph's 20,888 UniProtKB nodes carry a secondary or
            # deleted accession. That is what makes it safe to add now --
            # the behaviour is unchanged and the blind spot is named -- and
            # it is why the counters matter more than the recovery: a
            # future GtoP release citing a renamed accession becomes a
            # line of output instead of a silent subtraction.
            got = accessions.resolve(acc) if accessions is not None else None
            if got is not None:
                stats[f"not_in_graph_accession_{got.status}"] += 1
                if got.curie and got.curie in known_proteins:
                    stats["recovered_via_secondary_accession"] += 1
                    pid = got.curie
            if pid not in known_proteins:
                stats["protein_not_in_graph"] += 1
                continue
        lid = (r.get(C_LIGAND_ID) or "").strip()
        key = inchikeys.get(lid)
        if not key:
            stats["ligand_without_inchikey"] += 1
            continue
        if compounds is None:
            cid = f"INCHIKEY:{key}"
        else:
            cid, route = compounds.resolve(
                key, (ligand_chembl or {}).get(lid, ""))
            if cid is None:
                stats["compound_not_in_graph"] += 1
                continue
            # Counted per route, because a join that recovers edges by a
            # looser key has to be auditable: 661 by skeleton and 24 by
            # ChEMBL id on 2026.3.1, and if either number ever dwarfs the
            # exact matches the looser key is matching things it should not.
            if route != "exact":
                stats[f"compound_by_{route}"] += 1
        direction = direction_of(r.get(C_TYPE), r.get(C_ACTION))
        if direction == "unknown" and TYPE_TO_DIRECTION.get(
                (r.get(C_TYPE) or "").strip().upper()) is not None:
            # Type mapped, and Action disagreed with it. Counted so the
            # driver can report how often the source contradicts itself
            # rather than hiding it inside direction_unknown.
            stats["type_action_conflict"] += 1
        # The value is already in nM; the UNIT says which measurement it
        # is, and so which declared qualifier it belongs under. `relation`
        # is deliberately not carried: it is one field on the edge, an edge
        # can gather several potency types from several rows, and attaching
        # whichever row's relation won would be a line-number claim again.
        potency: dict[str, float] = {}
        raw_affinity = (r.get(C_AFFINITY) or "").strip()
        if not raw_affinity:
            # 10,852 of 24,968 rows in the 2026.3.1 export report no
            # MEDIAN affinity. Counted so that "this edge has no potency"
            # can be told from "this edge had one and it could not be read"
            # -- the distinction the unmapped-unit counter above exists
            # for. THE COUNTER IS NARROWER THAN THAT FIGURE: it is reached
            # only by rows that survived the species, accession and
            # compound-join filters above, which is 3,043 of the 10,852.
            stats["no_affinity_reported"] += 1
            # AND 940 OF THE 10,852 DO STATE AN AFFINITY (260 of the 3,043
            # that reach here), as a low and/or a high with no median
            # between them (fix-list P14). This layer still carries no
            # number for them: the qualifier the schema declares is a
            # point, and picking one out of a bare interval would invent a
            # measurement -- a range without a point is also exactly what
            # assemble drops, since potency_spread_field has to name a
            # qualifier that is there. Counted rather than folded into the
            # line above, which read as if GtoPdb had said nothing about
            # these rows.
            if ((r.get(C_AFF_LOW) or "").strip()
                    or (r.get(C_AFF_HIGH) or "").strip()):
                stats["affinity_range_without_median"] += 1
        if raw_affinity:
            qual = AFFINITY_UNIT_TO_QUALIFIER.get(
                (r.get(C_AFF_UNITS) or "").strip().upper())
            if qual is None:
                # NAMED, NOT JUST COUNTED. The edge still ships -- with a
                # direction and no number -- so it is indistinguishable
                # from one where GtoPdb reported no affinity at all, and
                # those are different facts. 50 rows on release 2026.3.1:
                # 43 are KB (the antagonist dissociation constant from a
                # functional assay, a real measure the schema has no
                # qualifier for), 6 are "-" and 1 is blank. Counting the
                # unit means a release that introduces a type shows up as
                # a new line rather than as a slightly smaller number.
                stats["affinity_unit_unmapped"] += 1
                unit = (r.get(C_AFF_UNITS) or "").strip().upper() or "(blank)"
                stats[f"affinity_unit:{unit}"] += 1
            else:
                try:
                    potency[qual] = float(raw_affinity)
                except ValueError:
                    stats["unparseable_affinity"] += 1
        pmids = {f"PMID:{p.strip()}"
                 for p in str(r.get(C_PMID) or "").split("|")
                 if p.strip().isdigit()}

        # ACCUMULATED, NOT FIRST-WINS. GtoP states a pair more than once --
        # 987 times on 2026.3.1, different papers for one compound and one
        # target -- and `if (cid, pid) in seen: continue` kept whichever row
        # the CSV happened to list first. On 61 of those pairs the rows
        # DISAGREE about direction, so the edge's mechanism was decided by a
        # line number: 24 pairs read `agonist` or `inhibitor` depending only
        # on file order.
        #
        # Reconciled through assemble.resolve_qualifier, the same function
        # that reconciles two LAYERS, because the question is identical and
        # the answer must not depend on where it is asked. It is also
        # order-independent, which is what makes the layer reproducible: the
        # rows arrive in whatever order the export lists them.
        prev = pending.get((cid, pid))
        if prev is None:
            pending[(cid, pid)] = {"direction": direction,
                                   "potency": potency, "pmids": pmids,
                                   # Every value seen per measurement type,
                                   # so the spread behind the point below
                                   # can be stated (fix-list P14).
                                   "values": {q: [v]
                                              for q, v in potency.items()}}
            continue
        stats["duplicate"] += 1
        settled = resolve_qualifier("direction", prev["direction"], direction,
                                    stats)
        if settled != prev["direction"]:
            stats["duplicate_changed_direction"] += 1
        prev["direction"] = settled
        prev["pmids"] |= pmids
        # The MOST POTENT reported value PER MEASUREMENT TYPE, by a stated
        # rule rather than by row order. One qualifier cannot hold several
        # papers' numbers, and the publications above say which papers they
        # were. Per type, because a Ki and an IC50 are not comparable and
        # taking the smaller of the two would invent a measurement.
        for qual, value in potency.items():
            have = prev["potency"].get(qual)
            prev["potency"][qual] = value if have is None else min(have, value)
            prev["values"].setdefault(qual, []).append(value)

    for (cid, pid), rec in pending.items():
        quals = {"direction": rec["direction"], **rec["potency"]}
        _attach_spread(quals, rec["values"], stats)
        em.edge(cid, "TARGETS", pid, source=GTOP_SOURCE, date=default_date,
                tier=GTOP_TIER, quals=quals,
                pmids=sorted(rec["pmids"]) or None)
        stats["edges"] += 1
        stats[f"direction_{rec['direction']}"] += 1
    return stats
