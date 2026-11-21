"""doc2/registrations/prereg_clinical_labels.md holds only operands, and they are real.

A registration signed before the measurement is worth exactly what its
operands are worth. This document counts compounds, diseases and class
sizes from data/raw and from the release BEFORE any ranking exists, and
every one of those counts is a number a later reader would otherwise
have to take on trust -- including the two the design turns on: that
SARS and COVID-19 are one cohort rather than two, and that 200 of the
222 in-trial compounds in the labelled population were measured
inactive.

The last test is the one that matters most: it fails if the document
ever acquires an AUC while no results file exists to carry it.
"""
from __future__ import annotations

import json
import re
import sqlite3
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav import virus_register as vr  # noqa: E402
from kgav.labels import build_labels  # noqa: E402
from kgav.pipeline import release_dir  # noqa: E402

PREREG = ROOT / "doc2" / "registrations" / "prereg_clinical_labels.md"
CHEMBL = (ROOT / "data" / "raw" / "chembl" / "chembl_37" / "chembl_37_sqlite"
          / "chembl_37.db")
MONDO = ROOT / "data" / "raw" / "mondo" / "mondo.json"

#: The xrefs the document's join uses, disease by disease. Written here
#: rather than derived so that a change to the join is a change to this
#: list and shows up as one -- the ingest will import THIS mapping, not
#: a second copy of it.
XREFS = {
    "MONDO:0002282": (["D014901"], ["DOID:2366"]),
    "MONDO:0005091": (["D045169"], ["EFO:0000694"]),
    "MONDO:0005502": (["D003715"], ["EFO:0005547"]),
    "MONDO:0005709": (["D003139"], ["EFO:0007214"]),
    "MONDO:0017572": (["D004675"], ["EFO:1001309"]),
    "MONDO:0018661": (["D000071243"], ["DOID:0060478"]),
    "MONDO:0019209": (["D004672"], ["EFO:0007332"]),
    "MONDO:0020502": (["D015004"], ["DOID:9682"]),
    "MONDO:0100096": (["C000657245", "D000086382"], ["MONDO:0100096"]),
    "MONDO:0100116": ([], ["DOID:0080642", "MONDO:0100116"]),
}

SARS_COV_2 = "NCBITaxon:2697049"


@pytest.fixture(scope="module")
def prereg() -> str:
    return PREREG.read_text()


@pytest.fixture(scope="module")
def indicated() -> dict[str, dict[str, float]]:
    """{mondo: {inchikey: max_phase}} straight from ChEMBL 37."""
    if not CHEMBL.exists():
        pytest.skip(f"{CHEMBL} absent")
    con = sqlite3.connect(f"file:{CHEMBL}?mode=ro", uri=True)
    out: dict[str, dict[str, float]] = {}
    for mondo, (mesh, other) in XREFS.items():
        clauses, args = [], []
        if mesh:
            clauses.append(f"di.mesh_id in ({','.join('?' * len(mesh))})")
            args += mesh
        if other:
            clauses.append(f"di.efo_id in ({','.join('?' * len(other))})")
            args += other
        if not clauses:
            out[mondo] = {}
            continue
        rows = con.execute(
            "select cs.standard_inchi_key, max(di.max_phase_for_ind) "
            "from drug_indication di "
            "join molecule_dictionary md using(molregno) "
            "join compound_structures cs using(molregno) "
            f"where {' or '.join(clauses)} "
            "group by cs.standard_inchi_key", args).fetchall()
        out[mondo] = {f"INCHIKEY:{k}": p for k, p in rows}
    con.close()
    return out


@pytest.fixture(scope="module")
def small_molecules() -> set[str]:
    nodes = release_dir(ROOT) / "nodes.jsonl"
    return {json.loads(ln)["id"] for ln in nodes.open()
            if '"SmallMolecule"' in ln}


@pytest.fixture(scope="module")
def labels():
    return build_labels(release_dir(ROOT))


def test_the_ten_mondo_ids_are_the_registers_own(prereg):
    """The document's disease list is the register's, not a selection.

    A registration that quietly dropped a disease with no indications
    would read as a complete join over a smaller world.
    """
    registered = {v["mondo"] for v in vr.load()["viruses"] if v.get("mondo")}
    assert set(XREFS) == registered, (
        "the join's disease list has drifted from config: "
        f"missing {registered - set(XREFS)}, extra {set(XREFS) - registered}")
    for mondo in registered:
        assert mondo in prereg, f"{mondo} is not named in the registration"


def test_every_disease_is_in_the_table_including_the_ones_that_yield_nothing(
        prereg, indicated, small_molecules):
    """Five of ten yield little or nothing, and all ten are printed.

    A table of the three that work is the selection failure this project
    is built against, so the two zero rows are asserted explicitly.
    """
    block = prereg.split("MONDO          disease")[1].split("1,358")[0]
    for mondo in XREFS:
        assert mondo in block, f"{mondo} missing from the yield table"
    for empty in ("MONDO:0017572", "MONDO:0020502"):
        assert not indicated[empty], (
            f"{empty} now has indications; the table's zero row is wrong")
        row = [ln for ln in block.splitlines() if empty in ln][0]
        assert row.split()[-3:] == ["0", "0", "0"], row


def test_sars_and_covid_are_one_cohort_which_is_why_the_family_is_one(
        prereg, indicated):
    """The load-bearing operand. If these two sets ever separate, the
    family-of-one alpha in this registration is wrong and the document
    has to be amended rather than reinterpreted."""
    cov = set(indicated["MONDO:0100096"])
    sars = set(indicated["MONDO:0005091"])
    shared = len(cov & sars)
    assert shared / min(len(cov), len(sars)) >= 0.90, (
        f"SARS and COVID-19 now share only {shared} of "
        f"{min(len(cov), len(sars))}; they are no longer one cell")
    assert f"share **{shared} of {max(len(cov), len(sars))}**" in prereg, (
        f"the registration does not state the measured overlap {shared} of "
        f"{max(len(cov), len(sars))}")
    assert "Family size 1, alpha = 0.05" in prereg.replace("\n", " ")


def test_the_population_split_is_the_releases_and_carries_its_direction(
        prereg, indicated, labels, small_molecules):
    """22 active / 200 inactive / 138 unlabelled, each asserted with its
    class name rather than as a bare integer.

    A bare count matches any document that happens to contain it, which
    is how an earlier publication-split test in this repository passed
    against a figure from a different table.
    """
    cov = {k for k in indicated["MONDO:0100096"] if k in small_molecules}
    pos = set(labels.positives[SARS_COV_2])
    neg = set(labels.negatives[SARS_COV_2])
    # The structure-bearing cohort, before the release's node filter: the
    # document's 360 is "with a structure in ChEMBL", not "in the graph".
    cohort = set(indicated["MONDO:0100096"])
    active, inactive = len(cohort & pos), len(cohort & neg)
    unlabelled = len(cohort - pos - neg)

    for label, n in (("in-vitro ACTIVE   vs SARS-CoV-2", active),
                     ("in-vitro INACTIVE vs SARS-CoV-2", inactive),
                     ("not labelled at all", unlabelled)):
        assert re.search(rf"{re.escape(label)}\s+{n}\b", prereg), (
            f"the registration does not state {n} against {label!r}")
    assert inactive > active, (
        "the registration's premise is that the clinical and in-vitro axes "
        "are near-orthogonal; they no longer are")

    in_trial = len(cohort & (pos | neg))
    screened = len((pos | neg) - cohort)
    assert f"{in_trial} positives" in prereg
    assert f"{screened:,} negatives" in prereg
    assert cov, "no indicated compound is a node, so the cell is empty"


def test_the_untestable_stratum_is_counted_so_it_cannot_become_testable(
        prereg, indicated, labels):
    """P4: 13 of 222 carry a viral-target annotation.

    Registered as not a cell. Without the count on record, a later run
    could introduce it as though 13 had been a surprise.
    """
    nodes = release_dir(ROOT) / "nodes.jsonl"
    viral = {json.loads(ln)["id"] for ln in nodes.open()
             if '"Protein"' in ln
             and json.loads(ln)["properties"].get("taxon_id") == SARS_COV_2}
    assert viral, "no SARS-CoV-2 proteins in the release"
    annotated = set()
    for ln in (release_dir(ROOT) / "edges.jsonl").open():
        if "INHIBITS" not in ln and "TARGETS" not in ln:
            continue
        e = json.loads(ln)
        if e["predicate"] in ("INHIBITS", "TARGETS") and e["object"] in viral:
            annotated.add(e["subject"])
    pop = set(labels.positives[SARS_COV_2]) | set(labels.negatives[SARS_COV_2])
    cohort = set(indicated["MONDO:0100096"]) & pop
    n = len(cohort & annotated)
    assert f"**{n}** carry an `INHIBITS`" in prereg, (
        f"the registration does not state {n} for P4")
    assert n < 30, (
        f"{n} compounds now carry a viral-target annotation; P4 may have "
        "become testable, which requires an amendment, not a reading")
    assert f"{len(viral)} SARS-CoV-2 proteins" in prereg


def test_the_document_states_whichever_is_true_of_the_two_predicates(prereg):
    """WRITTEN TO FAIL WHEN THE INGEST LANDED, AND IT DID.

    The first version of this test asserted both predicates were still
    empty, with a note that the failure is what would tell a reader the
    document had gone historical. The ingest ran the same day and the
    census went 0/0 -> 568/9, so the test fired and this is the amended
    form: it requires the document to carry the amendment, with the
    counts the release actually holds, and still requires the signature
    above it to say the registration preceded the ingest.
    """
    census = (ROOT / "doc2" / "kg_release_census.txt").read_text()
    trial = re.search(r"IN_TRIAL_FOR\s+L\s+([\d,]+)", census)
    treats = re.search(r"TREATS\s+L\s+([\d,]+)", census)
    assert trial and treats, "the census no longer reports both predicates"
    n_trial = int(trial.group(1).replace(",", ""))
    n_treats = int(treats.group(1).replace(",", ""))

    assert "SIGNED 2026-11-17, BEFORE THE INGEST AND BEFORE ANY AUC" in prereg
    if n_trial == 0 and n_treats == 0:
        assert "# First amendment" not in prereg, (
            "the document claims the ingest ran but the release holds no "
            "indication edge")
        return
    assert "# First amendment" in prereg, (
        f"the release holds {n_trial} IN_TRIAL_FOR and {n_treats} TREATS, so "
        "this registration reads as a plan for work already done")
    amendment = prereg.split("# First amendment")[1]
    assert f"{n_trial}    IN_TRIAL_FOR" in amendment, (
        f"the amendment does not state the {n_trial} edges that landed")
    assert f"{n_treats}         {n_treats}    TREATS" in amendment, (
        f"the amendment does not state the {n_treats} TREATS edges")
    # The gap between what was registered and what shipped has to be
    # named, not absorbed: one edge was refused by required_qualifiers.
    assert "in_trial_without_nct" in amendment


def test_the_registration_claims_no_auc_while_no_results_file_exists(prereg):
    """The whole value of a signature.

    Scans the document for anything shaped like a measured AUC. The
    operand counts above are integers and percentages, so a four-decimal
    0.x figure can only be a score -- and none may appear until
    data/results carries a file to regenerate it from.
    """
    existing = sorted((ROOT / "data" / "results").glob("*clinical*"))
    body = re.sub(r"`[^`]*`", "", prereg)
    body = body.replace("0.05", "").replace("0.5–3", "")
    scores = re.findall(r"\b0\.\d{3,4}\b", body)
    if not existing:
        assert not scores, (
            f"the registration quotes {scores} but no clinical-label results "
            "file exists to regenerate them from")
    else:
        assert scores, (
            f"{[p.name for p in existing]} exists, so the registration's "
            "verdict section should now quote its figures")


# ------------------------------------------------- the 2026-11-18 amendment
def _amendment() -> str:
    """The appended correction, or a failure saying it is missing."""
    doc = PREREG.read_text()
    marker = "# Amendment, 2026-11-18"
    assert marker in doc, (
        "this document asserts a reason for declining the rebuild that was "
        "measured to be false; the correction must be appended")
    return doc[doc.index(marker):]


def test_the_falsified_sentence_is_corrected_and_not_deleted():
    """A REGISTRATION IS APPENDED TO, NEVER REWRITTEN.

    Both halves matter. If the original sentence were edited away, the
    record would no longer show what was signed; if the amendment were
    missing, the document would assert something the repository
    disproves. So the wrong sentence must still be there AND be
    answered below it.
    """
    doc = PREREG.read_text()
    wrong = "would re-fetch BioGRID ORCS over the network"
    assert wrong in doc, "the original claim was rewritten instead of amended"
    tail = _amendment()
    assert wrong in tail, "the amendment does not quote what it corrects"
    assert "cache-first" in tail


def test_fetch_screen_really_is_cache_first():
    """THE CLAIM, DERIVED FROM THE CODE RATHER THAN FROM THE AMENDMENT.

    A document correcting a document proves nothing. `fetch_screen` must
    return a cached screen before it builds a URL, or the original
    sentence was right after all and this amendment is the error.
    """
    src = (ROOT / "scripts" / "ingest_orcs.py").read_text()
    body = src[src.index("def fetch_screen("):]
    body = body[:body.index("\ndef ", 1)] if "\ndef " in body[1:] else body
    read = body.index("read_text()")
    url = body.index("urlopen")
    assert read < url, (
        "fetch_screen reaches the network before reading its cache, so the "
        "amendment's cache-first claim is wrong")
    assert "if f.exists()" in body


def test_the_amendment_does_not_claim_the_decision_was_wrong():
    """The decision held; only its stated reason was false. An amendment
    that overcorrected into 'the deferral was a mistake' would misstate
    the record in the other direction."""
    tail = _amendment()
    assert "THE DECISION WAS STILL RIGHT" in tail


def test_the_rebuild_figure_is_the_measured_one():
    """The amendment reports the rebuild it was declining. The runtime
    and the two shas are facts, and the shas must match the release the
    results on disk were measured against."""
    tail = _amendment()
    # release_fingerprint, NOT MANIFEST.json. The first version of this
    # test read manifest["nodes_sha"], which does not exist -- the
    # manifest carries `release`, `code` and the census, and the shas
    # come from kgav.provenance. `if sha:` then skipped both keys and
    # the test asserted NOTHING while reading green. A check that cannot
    # fail is worse than no check: it occupies the place where one
    # should be.
    from kgav.provenance import release_fingerprint

    release = ROOT / "data" / "releases" / "v0.1"
    if not (release / "nodes.jsonl").exists():
        pytest.skip("no assembled release")
    fp = release_fingerprint(release)
    for key in ("nodes_sha", "edges_sha"):
        sha = fp[key]
        assert sha, f"release_fingerprint returned no {key}"
        assert sha in tail, (
            f"the amendment does not name the release's {key} ({sha})")
