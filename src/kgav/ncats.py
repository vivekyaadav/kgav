"""NCATS OpenData SARS-CoV-2 CPE screen: labels in the population that matters.

    SmallMolecule -HAS_ANTIVIRAL_ACTIVITY_AGAINST-> OrganismTaxon

WHY THIS LAYER EXISTS. Four evaluation protocols have now scored at chance or
below once their confounds were controlled (docs/evaluation.md). Two of the six
measured causes are properties of the LABEL SET, not of the graph:

  * Negatives cannot reach the host-directed metapaths. 1.5% of
    MEASURED_INACTIVE compounds reach a dependency host factor against 12.1%
    expected -- 7.9x under-represented. ChEMBL negatives are antiviral research
    compounds, which carry viral-target annotations and almost no host ones.
  * 94.9% of compounds carrying both an INHIBITS edge and an antiviral label
    cite the same PMID on both, because one paper reports the target IC50 and
    the cell-based EC50 together.

A PANEL fixes both at once, and this one in particular:

  * It reports the compounds that did NOT work. 8,810 screened, 56 confirmed
    hits -- so the overwhelming majority of rows are measured inactives drawn
    from the approved-and-investigational library, which is the population
    holding the TARGETS edges.
  * It is cell-based and generated independently of ChEMBL target annotation,
    so no edge this layer's labels could leak through cites the same paper.
    --publication-disjoint therefore withholds nothing when scoring against
    these labels, which is the whole point of using them.
  * EC50 and CC50 come from the same plate, so the selectivity index is
    same-document by construction -- the strongest form this project accepts.

WHAT IT DOES NOT FIX. The screen ran in 2020, so with cutoff 2021 these labels
are PRE-cutoff and do not enter the post-2021 test set. This layer serves the
cross-sectional protocol, which the publication-disjoint result made the
primary one anyway.

IDENTITY IS THE RISK HERE, not parsing. The graph's compound keys came out of
normalize_chemical(); a key taken from this file instead would differ on salt,
stereo or tautomer and the new nodes would sit beside the existing ones without
joining to a single TARGETS edge. Nothing would error and no count would look
wrong -- the README calls a discarded identifier the dominant bug class, found
four times, silent every time. So: structures go through the project's own
normalizer, and the driver refuses to write a release whose join rate against
the existing graph is implausible.
"""
from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

ANTIVIRAL_PREDICATE = "HAS_ANTIVIRAL_ACTIVITY_AGAINST"

# Chen CZ et al., Front Pharmacol 2021;11:592737. Vero E6, MOI 0.002, 72h.
CPE_PUBLICATION = "PMID:33708112"
SOURCE = "infores:ncats-opendata"
# Vero E6, the line this screen used. Cellosaurus: Vero E6 is CVCL_0574, which
# is NOT Vero (CVCL_0059) -- a distinct subclone selected for ACE2 expression.
VERO_E6 = "CVCL:0574"
SARS_COV_2 = "NCBITaxon:2697049"

# First public assertion: bioRxiv 2020.08.18.255877. Not the portal's file
# mtime (2023-05-22), which is a re-export date and would put a 2020 screen
# after both PPI layers. Overridable; never inferred.
DEFAULT_DATE = "2020-08-18"

# A compound is only callable INACTIVE if the assay reached a concentration at
# or above the threshold the rest of the pipeline uses for inactivity. Tested
# to 5 uM and found inactive says nothing about 10 uM, and labels.classify()
# would silently return None for it; counting it as a negative anyway would
# invent evidence. Such rows are reported as `below_threshold`, not dropped
# quietly.
DEFAULT_INACTIVE_AT_NM = 10_000.0


class ColumnError(ValueError):
    """Raised when the header cannot be mapped. Carries the whole header.

    Deliberately fatal. Guessing a column here produces a release that
    validates, counts plausibly, and encodes the wrong measurement.
    """


# field -> header spellings seen across NCATS ODP exports, lowercased and
# stripped of punctuation by _key(). Extend rather than loosen the matcher.
COLUMN_CANDIDATES: dict[str, tuple[str, ...]] = {
    "smiles": ("smiles", "structure", "canonicalsmiles"),
    "sample": ("ncgcsid", "sid", "samplename", "sample", "ncgcid", "name"),
    "activity": ("activity", "activityclass", "outcome", "calledactivity"),
    "ac50": ("ac50", "ac50m", "ac50um", "ec50", "ec50um", "ec50m",
             "cpeac50", "efficacyac50"),
    "lac50": ("lac50", "logac50", "logac50m"),
    "efficacy": ("efficacy", "efficacypct", "efficacy%", "maxresponse"),
    "curve": ("curveclass", "cclass", "curveclass2", "cc"),
    "cc50": ("cc50", "cc50m", "cc50um", "cytotoxac50", "cytotoxicityac50",
             "tox ac50", "toxac50"),
}


def _key(s: str) -> str:
    return "".join(c for c in s.strip().lower() if c.isalnum())


def resolve_columns(header: list[str], required=("smiles",)) -> dict[str, str]:
    """Map our field names onto this file's actual header.

    Returns only what it found. `required` must all be present or it raises:
    without a structure column there is nothing to normalize, so a partial
    match is not a degraded mode, it is a different file.
    """
    index = {_key(h): h for h in header}
    found: dict[str, str] = {}
    for field, cands in COLUMN_CANDIDATES.items():
        for c in cands:
            if _key(c) in index:
                found[field] = index[_key(c)]
                break
    missing = [f for f in required if f not in found]
    if missing:
        raise ColumnError(
            f"could not find {missing} in header.\n"
            f"  header seen ({len(header)} columns): {header}\n"
            f"  add the real spelling to COLUMN_CANDIDATES in src/kgav/ncats.py "
            f"rather than renaming the file's columns, so the next export of "
            f"this dataset keeps working."
        )
    return found


def read_rows(path: Path) -> tuple[list[dict], dict[str, str]]:
    """Read the TSV and resolve its header. Returns (rows, column mapping)."""
    with Path(path).open(newline="", encoding="utf-8-sig") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        if not rdr.fieldnames:
            raise ColumnError(f"{path} has no header row")
        cols = resolve_columns(list(rdr.fieldnames))
        return list(rdr), cols


def to_nm(value: str | float | None, unit: str) -> float | None:
    """Concentration to nanomolar. `unit` is one of M, uM, nM, logM.

    logM is NCATS's LAC50: log10 of molar. -5.0 is 10 uM, i.e. 10,000 nM.
    """
    if value is None or value == "":
        return None
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    if unit == "nM":
        return v
    if unit == "uM":
        return v * 1_000.0
    if unit == "M":
        return v * 1e9
    if unit == "logM":
        return (10.0 ** v) * 1e9
    raise ValueError(f"unknown concentration unit {unit!r}")


def unit_for(column_name: str, default: str = "uM") -> str:
    """Infer the unit from the column's own name, since NCATS varies it.

    Explicit beats inferred: the driver takes --ac50-units, and this only
    supplies its default.
    """
    k = _key(column_name)
    if k.startswith("lac50") or k.startswith("logac50"):
        return "logM"
    if k.endswith("m") and not k.endswith("um") and not k.endswith("nm"):
        return "M"
    if k.endswith("nm"):
        return "nM"
    if k.endswith("um"):
        return "uM"
    return default


# Values in an NCATS `activity`/outcome column meaning the compound was tested
# and did nothing. Anything unrecognised is counted, never assumed either way.
INACTIVE_WORDS = {"inactive", "inconclusive inactive", "not active", "inert"}
ACTIVE_WORDS = {"active", "activeagonist", "activeantagonist", "activeinhibitor",
                "inconclusive active", "activator", "inhibitor"}


def classify_row(row: dict, cols: dict[str, str], *, ac50_units: str,
                 cc50_units: str, max_conc_nm: float,
                 inactive_at_nm: float = DEFAULT_INACTIVE_AT_NM,
                 ) -> tuple[str | None, dict, str]:
    """One screen row -> (verdict, qualifiers, reason).

    verdict is "active", "inactive" or None. The qualifiers are written so that
    labels.classify() reaches the SAME verdict from the emitted edge alone:
    an active carries relation "=" with its fitted EC50, an inactive carries
    relation ">" with the highest concentration tested. That symmetry is the
    contract -- this function must not be the only place the call is recorded.
    """
    ac50_col = cols.get("ac50") or cols.get("lac50")
    ac50 = to_nm(row.get(ac50_col), ac50_units) if ac50_col else None
    cc50 = to_nm(row.get(cols["cc50"]), cc50_units) if "cc50" in cols else None

    quals: dict = {"assay_type": "cell_based_antiviral", "cell_line": VERO_E6,
                   "virus_strain": "USA-WA1/2020"}

    called = _key(row.get(cols["activity"], "")) if "activity" in cols else ""
    inactive_called = called in {_key(w) for w in INACTIVE_WORDS}
    active_called = called in {_key(w) for w in ACTIVE_WORDS}

    # An AC50 present and finite is the strongest statement the row makes, and
    # it outranks a text outcome column: the outcome is a curve-class heuristic
    # applied by the depositor, the AC50 is the fit.
    if ac50 is not None and ac50 > 0:
        quals["ec50_nm"] = ac50
        quals["relation"] = "="
        if cc50 is not None and cc50 > 0:
            quals["cc50_nm"] = cc50
            quals["selectivity_index"] = cc50 / ac50
            quals["selectivity_verified"] = True
            quals["selectivity_n_studies"] = 1
        if "efficacy" in cols:
            try:
                quals["source_efficacy_pct"] = float(row[cols["efficacy"]])
            except (TypeError, ValueError, KeyError):
                pass
        return "active", quals, "fitted_ac50"

    if inactive_called:
        # No fitted AC50 and the depositor called it inactive: the compound was
        # tested to max_conc and did not act. Recorded as a CENSORED value at
        # that concentration, which is exactly what was observed, and is the
        # form labels.classify() already understands.
        if max_conc_nm < inactive_at_nm:
            return None, quals, "below_threshold"
        quals["ec50_nm"] = max_conc_nm
        quals["relation"] = ">"
        if cc50 is not None and cc50 > 0:
            quals["cc50_nm"] = cc50
        return "inactive", quals, "called_inactive"

    if active_called:
        # Called active but no usable AC50. Asserting activity with no value is
        # what `unquantified` is for, and that field is DERIVED at assembly
        # (see chembl.py) -- so emit the edge without a value and let assembly
        # set it. Do not invent a potency.
        quals["relation"] = "="
        return "active", quals, "called_active_unquantified"

    return None, quals, "no_call"


def ingest_cpe(em, rows: list[dict], cols: dict[str, str], normalizer, *,
               existing_compounds: frozenset[str] | set[str] = frozenset(),
               taxon: str = SARS_COV_2, source: str = SOURCE,
               date: str = DEFAULT_DATE, ac50_units: str = "uM",
               cc50_units: str = "uM", max_conc_nm: float = 46_000.0,
               inactive_at_nm: float = DEFAULT_INACTIVE_AT_NM) -> Counter:
    """Emit one antiviral-activity edge per usable row.

    `normalizer` is injected so tests need no rdkit and so the function cannot
    quietly fall back to the file's own identifiers: it takes a structure
    string and returns an object with .curie and .inchikey, or raises.

    `existing_compounds` is why this layer does not become the authority on
    `is_approved`. That property is required on SmallMolecule and boolean, with
    no unknown, and this file does not say whether a compound is approved. A
    layer that asserts it anyway is the H2 bug again -- a true statement about
    one layer that becomes false after the merge, resolved by layer precedence
    rather than by evidence. So a compound the graph already holds gets NO node
    from here, only its label edge: the layer that knows keeps the property.
    Genuinely new compounds get is_approved False, which is a floor rather than
    a claim -- no metapath or filter reads the field, so an understated value
    costs nothing while an overstated one would be a false assertion.
    """
    stats: Counter = Counter()
    seen: set[str] = set()

    for row in rows:
        stats["rows"] += 1
        struct = (row.get(cols["smiles"]) or "").strip()
        if not struct:
            stats["no_structure"] += 1
            continue
        try:
            cid = normalizer(struct)
        except Exception:
            stats["unnormalizable"] += 1
            continue

        verdict, quals, reason = classify_row(
            row, cols, ac50_units=ac50_units, cc50_units=cc50_units,
            max_conc_nm=max_conc_nm, inactive_at_nm=inactive_at_nm)
        stats[f"reason_{reason}"] += 1
        if verdict is None:
            continue

        drug = cid.curie
        if drug in seen:
            # One structure, two plate records. Keeping both would double its
            # DWPC weight and let one compound vote twice in the AUC.
            stats["duplicate_structure"] += 1
            continue
        seen.add(drug)

        if drug in existing_compounds:
            stats["joined_existing_compound"] += 1
        else:
            em.node(drug, "SmallMolecule",
                    smiles=struct,
                    inchikey_skel=cid.inchikey[:14],
                    # See the docstring: a floor, not a claim.
                    is_approved=False,
                    salt_collapsed=False, stereo_collapsed=False)
            stats["new_compound_node"] += 1
        em.edge(drug, ANTIVIRAL_PREDICATE, taxon,
                source=source, date=date, tier=1, quals=quals,
                pmids=[CPE_PUBLICATION])
        stats[verdict] += 1
        if quals.get("selectivity_verified"):
            stats["with_same_plate_si"] += 1

    return stats


def join_report(labelled: set[str], existing: set[str]) -> dict:
    """How much of this layer lands on compounds the graph already holds.

    `labelled` is every compound this layer wrote an edge for -- taken from the
    edges, not from em.nodes, because a compound that joined deliberately gets
    no node from this layer (see ingest_cpe).

    The number that matters is not how many nodes were created but how many
    joined, because a node that joins nothing contributes no path. Reported so
    a bad normalization shows up as a count rather than as a null result three
    scripts later -- cross-entity count comparison is this project's primary
    detector and it only works if the count is printed.
    """
    joined = labelled & existing
    return {"emitted": len(labelled), "already_in_graph": len(joined),
            "new_compounds": len(labelled - existing),
            "join_rate": (len(joined) / len(labelled)) if labelled else 0.0}
