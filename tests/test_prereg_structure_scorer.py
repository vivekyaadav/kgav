"""The structure-scorer registration, checked against the release it names.

WHY THIS FILE EXISTS. docs/prereg_structure_scorer.md fixes a protocol
before the measurement, and every number in it is a property of the release
or of the labels -- class sizes, SMILES coverage, which cells are readable,
the corrected alpha. A registration whose figures were typed by hand drifts,
and this project has already paid for that twice in one document: the
pathway-overlap registration maintained a release-identification count by
hand and it read "eight" when git said eleven, with a per-date breakdown
that matched neither.

So the figures are re-derived here from the same functions the run will use.

AND THE FIRST DRAFT OF THAT DOCUMENT NEEDED IT. Written by hand it quoted
the POOLED labelled-compound count (18,998) for a ChEMBL-only study
(13,417), the wrong edge count for the release it names, a skeleton-group
figure from the pooled population, and four virus names that were simply
wrong -- 11082 as a flavivirus when it is West Nile, 11072 as an
enterovirus when it is JEV, 11084 as a morbillivirus when it is TBEV, and a
West Nile PMID figure attributed to HCoV-NL63. Every one of those is the
kind of error a derived figure cannot make.

WHAT THIS FILE DELIBERATELY DOES NOT DO. It computes no score and no AUC.
The registration is unsigned-in-spirit the moment a test asserts an outcome,
and the scorer does not exist yet.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

PREREG = ROOT / "docs/prereg_structure_scorer.md"


@pytest.fixture(scope="module")
def prereg() -> str:
    assert PREREG.exists(), "the registration is missing"
    return " ".join(PREREG.read_text().split())


@pytest.fixture(scope="module")
def design():
    """The registered population's design facts, from the release."""
    from kgav.labels import build_labels
    from kgav.pipeline import release_dir
    import diagnose_structure_scorer_design as D

    rd = release_dir(ROOT)
    if not (rd / "nodes.jsonl").exists():
        pytest.skip("no assembled release on disk")
    props = D.small_molecules(rd)
    labels = build_labels(rd, sources={"infores:chembl"})
    return D, props, labels


def test_it_says_it_was_signed_before_the_code_and_the_run(prereg):
    """The claim that distinguishes this registration from the pathway one.

    THE FIRST VERSION OF THIS TEST DEMANDED THE WRONG THING. It required the
    phrase to DISAPPEAR once scripts/structure_scorer.py existed -- as
    though the claim expired. It does not: "signed before the code" is a
    permanent fact about 2026-11-11, and editing a registration to match
    what happened later destroys the only thing it is for, which is what
    prereg_pathway_overlap.md says at length about its own release
    identification.

    What a reader is owed once the code exists is the ORDERING, dated, so
    the claim can be checked rather than taken. So the phrase stays and the
    document must additionally record when the code was written.
    """
    assert "SIGNED 2026-11-11, BEFORE THE CODE AND BEFORE THE RUN" in prereg
    if (ROOT / "scripts/structure_scorer.py").exists():
        assert "## When the code and the run happened" in prereg, (
            "scripts/structure_scorer.py exists and the registration does "
            "not record when it was written; the ordering claim is then "
            "unverifiable")
        assert "`scripts/structure_scorer.py` was written on" in prereg, (
            "the registration must date the code it was signed before")


def test_the_quoted_release_is_the_one_on_disk(prereg):
    """DERIVED FROM THE CENSUS. The pathway registration's equivalent test
    exists because a hand-typed release identification went stale, and the
    first draft of THIS document named 392,363 edges when the release has
    392,375."""
    census = ROOT / "docs/kg_release_census.txt"
    if not census.exists():
        pytest.skip("census absent")
    body = census.read_text()
    m = re.search(r"([\d,]+) nodes\s+([\d,]+) edges", body)
    assert m, "the census no longer states a node and edge count"
    for token in (f"{m.group(1)} nodes", f"{m.group(2)} edges"):
        assert token in prereg, (
            f"the census says {token} and the registration does not quote it")
    for sha in sorted(set(re.findall(r"\b[0-9a-f]{16}\b", body))):
        assert sha in prereg, f"the census states {sha} and the prereg does not"


def test_the_premise_holds_every_labelled_compound_has_a_smiles(prereg, design):
    """THE PREMISE OF THE WHOLE STUDY. A structure scorer's floor is 0.5
    only because its coverage gap is zero. If one labelled compound has no
    SMILES the registered P3 is already false."""
    D, props, labels = design
    labelled = set()
    for v in set(labels.positives) | set(labels.negatives):
        labelled |= labels.positives.get(v, set())
        labelled |= labels.negatives.get(v, set())
    with_smiles = {i for i in labelled
                   if (props.get(i, {}).get("smiles") or "")}
    assert with_smiles == labelled, (
        f"{len(labelled - with_smiles):,} labelled compound(s) carry no "
        f"SMILES; the registration's premise of complete coverage is false")
    assert f"{len(labelled):,} of {len(labelled):,}" in prereg, (
        f"the registered population holds {len(labelled):,} labelled "
        f"compounds and the registration does not quote that figure")


def test_the_registered_cells_are_the_readable_ones(prereg, design):
    """The correction family is the readable cells, so the count and the
    membership both have to be the release's and not a reading of it."""
    from kgav.calibration import MIN_REACHED_POSITIVES
    D, props, labels = design
    viruses = sorted(set(labels.positives) | set(labels.negatives))
    readable, unreadable = [], []
    for v in viruses:
        pos = labels.positives.get(v, set())
        neg = labels.negatives.get(v, set())
        (readable if (len(pos) >= MIN_REACHED_POSITIVES
                      and len(neg) >= MIN_REACHED_POSITIVES)
         else unreadable).append((v, len(pos), len(neg)))

    m = re.search(r"alpha = 0\.05/(\d+) = ([\d.]+)", prereg)
    assert m, "the registration no longer states its corrected alpha"
    family = int(m.group(1))
    assert family == len(readable), (
        f"{len(readable)} cells are readable on this release and the "
        f"registration fixes the family at {family}")
    stated_alpha = float(m.group(2))
    assert abs(stated_alpha - 0.05 / family) < 5e-7, (
        f"0.05/{family} is {0.05 / family:.6f}, not {stated_alpha}")

    # and every excluded cell must be named with its counts
    names = D.virus_names()
    for v, npos, nneg in unreadable:
        n = names.get(v, v)
        assert n in prereg, (
            f"{n} is not readable ({npos} positives, {nneg} negatives) and "
            f"the registration does not name it; a cell dropped silently is "
            f"the powerless row this project has relabelled twice")


def test_every_registered_class_size_is_the_releases(prereg, design):
    """THE POWER TABLE, CELL BY CELL. Fourteen rows of positives against
    negatives, each a count the registration will be read against."""
    from kgav.calibration import MIN_REACHED_POSITIVES
    D, props, labels = design
    names = D.virus_names()
    checked = 0
    for v in sorted(set(labels.positives) | set(labels.negatives)):
        pos = len(labels.positives.get(v, set()))
        neg = len(labels.negatives.get(v, set()))
        if pos < MIN_REACHED_POSITIVES or neg < MIN_REACHED_POSITIVES:
            continue
        pair = f"{pos:,} / {neg:,}"
        assert pair in prereg, (
            f"{names.get(v, v)} is {pair} on this release and the "
            f"registration's table does not carry that pair")
        checked += 1
    assert checked >= 14, f"only {checked} cells checked"


def test_the_minimum_detectable_auc_is_recomputed(prereg, design):
    """Each cell's registered detectable floor, from the project's own
    interval rather than from the document."""
    from kgav.calibration import MIN_REACHED_POSITIVES
    D, props, labels = design
    for v in sorted(set(labels.positives) | set(labels.negatives)):
        pos = len(labels.positives.get(v, set()))
        neg = len(labels.negatives.get(v, set()))
        if pos < MIN_REACHED_POSITIVES or neg < MIN_REACHED_POSITIVES:
            continue
        det = D.detectable_auc(pos, neg)
        assert det is not None, (v, pos, neg)
        assert f"{det:.4f}" in prereg, (
            f"{D.virus_names().get(v, v)} at {pos}/{neg} can first clear "
            f"chance at {det:.4f} and the registration does not state it")


def test_the_corrected_z_matches_the_family(prereg):
    """z is what the intervals will be built from, so it is derived from the
    registered family rather than carried beside it."""
    import math
    m = re.search(r"alpha = 0\.05/(\d+) = [\d.]+, \*\*z = ([\d.]+)\*\*",
                  prereg)
    assert m, "the registration no longer states alpha and z together"
    family, stated = int(m.group(1)), float(m.group(2))
    alpha = 0.05 / family
    lo, hi = 0.0, 10.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if math.erfc(mid / math.sqrt(2)) > alpha:
            lo = mid
        else:
            hi = mid
    assert abs(stated - (lo + hi) / 2) < 5e-5, (
        f"two-sided z for alpha={alpha:.6f} is {(lo + hi) / 2:.4f}, not "
        f"{stated}")


def test_the_secondary_cell_is_the_panels(prereg):
    """The NCATS matched cell carries the registration's load-bearing
    prediction P5, so its class sizes are derived too.

    AND SOURCE ALONE NO LONGER NAMES IT. Until 2026-11-17 this test asked
    for `sources={"infores:ncats-opendata"}` because that source published
    one screen, and the registration's sentence said "one cell". The source
    now publishes four screens; the same expression gives three cells and a
    SARS-CoV-2 figure of 149 / 7,474 that POOLS two protocols -- which is
    the confound the registration's population decision refuses one level
    up. The cell is identified by source AND protocol, by the same constant
    the scorer runs, so the test and the run cannot drift apart.
    """
    from kgav.labels import build_labels
    from kgav.pipeline import release_dir
    import structure_scorer as S
    rd = release_dir(ROOT)
    if not (rd / "nodes.jsonl").exists():
        pytest.skip("no assembled release on disk")
    spec = [(src, at) for key, src, at in S.POPULATIONS
            if key == "secondary_ncats"]
    assert len(spec) == 1, "the scorer no longer runs one secondary_ncats"
    source, assay_types = spec[0]
    assert assay_types, (
        "the secondary population is pinned by source alone; that source "
        "publishes more than one protocol and the cell would be pooled")
    panel = build_labels(rd, sources={source}, assay_types=assay_types)
    cells = {v: (len(panel.positives.get(v, set())),
                 len(panel.negatives.get(v, set())))
             for v in set(panel.positives) | set(panel.negatives)}
    readable = {v: c for v, c in cells.items() if c[0] >= 10 and c[1] >= 10}
    assert len(readable) == 1, (
        f"the registration calls the panel one cell and it is {len(readable)}")
    (pos, neg), = readable.values()
    assert f"{pos} positives against {neg:,} negatives" in prereg, (
        f"the panel cell is {pos} against {neg:,} and the registration does "
        f"not state it")

    # AND THE POOLED FIGURE MUST BE ON THE RECORD, because the expression
    # the registration was signed with still runs and still returns it.
    pooled = build_labels(rd, sources={source})
    for virus, (p, n) in sorted(
            ((v, (len(pooled.positives.get(v, set())),
                  len(pooled.negatives.get(v, set()))))
             for v in set(pooled.positives) | set(pooled.negatives)),
            key=lambda r: -r[1][0]):
        assert f"{p:,} / {n:,}" in prereg, (
            f"source alone gives {virus} as {p:,} / {n:,} and the "
            f"registration does not state what the signed expression now "
            f"returns")


def test_the_fingerprint_is_imported_not_restated(prereg):
    """The registration promises the scorer reuses the similarity layer's
    fingerprint. Read the constants, not the sentence."""
    from kgav.similarity import FP_BITS, FP_RADIUS
    assert f"radius {FP_RADIUS}, {FP_BITS} bits" in prereg, (
        f"kgav.similarity builds ECFP at radius {FP_RADIUS} over {FP_BITS} "
        f"bits and the registration states something else")


def test_the_diagnostic_computes_no_similarity(prereg):
    """THE BOUNDARY, ENFORCED. The design diagnostic exists so the facts
    above could be gathered without producing the study's result. If it ever
    computes a Tanimoto, the separation between the two classes' score
    distributions IS the registered outcome, and gathering it before the run
    would make this registration a description."""
    src = (ROOT / "scripts/diagnose_structure_scorer_design.py").read_text()
    code = "\n".join(ln for ln in src.splitlines()
                     if not ln.lstrip().startswith("#"))
    body = code.split('"""', 2)[-1]          # past the module docstring
    for banned in ("TanimotoSimilarity", "BulkTanimoto", "GetMorganGenerator",
                   "fingerprints(", "DataStructs"):
        assert banned not in body, (
            f"the design diagnostic calls {banned}; it is registered as "
            f"computing no score, and a similarity distribution over the "
            f"label classes is the outcome itself")


def test_the_pooled_figure_it_quotes_beside_the_registered_one(prereg, design):
    """THE FIGURE I FABRICATED. Updating this registration after the L27
    rebuild, I replaced the pooled count `18,998 of 18,998` with
    `18,914 of 18,914` -- a number I did not measure and which was simply
    wrong; the pooled population holds 18,900. Every other figure in the
    document was pinned by a test and this one was not, so it passed.

    The registration quotes the pooled count only to say that the
    REGISTERED population is the ChEMBL one, so the contrast is the whole
    point of the sentence. It is now derived too.
    """
    from kgav.labels import build_labels
    from kgav.pipeline import release_dir
    rd = release_dir(ROOT)
    if not (rd / "nodes.jsonl").exists():
        pytest.skip("no assembled release on disk")
    pooled = build_labels(rd)
    ids = set()
    for v in set(pooled.positives) | set(pooled.negatives):
        ids |= pooled.positives.get(v, set()) | pooled.negatives.get(v, set())
    assert f"{len(ids):,} of {len(ids):,} pooled" in prereg, (
        f"the pooled population holds {len(ids):,} labelled compounds and "
        f"the registration quotes a different figure beside the registered "
        f"one")


def test_the_amendments_publication_split_is_the_releases(prereg):
    """THE AMENDMENT'S CLAIM THAT P5 STILL REFUSES, derived.

    The 2026-11-17 amendment re-derives P5's refusal rather than inheriting
    it, because the merge changed the publications on that cell's edges.
    The refusal survives only if EVERY edge in the cell still shares a PMID
    with every other -- which is a property of the release, so the figures
    the amendment states for it are read off the release here.
    """
    import collections
    import json
    from kgav.labels import ACTIVITY_PREDICATE, edge_sources
    from kgav.pipeline import release_dir
    import structure_scorer as S
    rd = release_dir(ROOT)
    if not (rd / "edges.jsonl").exists():
        pytest.skip("no assembled release on disk")

    plate = [(src, at) for _k, src, at in S.POPULATIONS if at]
    sources = {src for src, _ in plate}
    groups: dict[tuple, collections.Counter] = collections.defaultdict(
        collections.Counter)
    for line in (rd / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        if e["predicate"] != ACTIVITY_PREDICATE:
            continue
        if not (edge_sources(e) & sources):
            continue
        at = (e.get("qualifiers") or {}).get("assay_type")
        pubs = tuple(sorted(p for p in (e.get("publications") or []) if p))
        groups[(at, e["object"])][pubs] += 1

    # THE REGISTERED SECONDARY CELL carries P5, so its split is stated.
    sec = [(src, at) for key, src, at in S.POPULATIONS
           if key == "secondary_ncats"]
    (sec_src, sec_at), = sec
    sec_keys = [k for k in groups if k[0] in sec_at]
    assert len(sec_keys) == 1, f"{len(sec_keys)} cells for the secondary"
    counts = groups[sec_keys[0]]
    assert len(counts) > 1, (
        "the secondary cell's edges all cite the same publication list; "
        "the amendment re-derives P5 from a SPLIT and there is none, so "
        "this test is checking nothing")
    plain = prereg.replace("**", "")   # emphasis is not data
    for pubs, n in sorted(counts.items(), key=lambda r: -r[1]):
        # Counted WITH its publication list, because a bare "2,622" could
        # be satisfied by any four-digit figure elsewhere in the document
        # -- and a bare "1" by every document ever written. The first
        # draft of this test asserted the bare count and passed against a
        # DENV-1 class size.
        cited = " and ".join(f"`{p}`" for p in pubs)
        claim = f"{n:,} CPE-typed | {cited}"
        assert claim in plain, (
            f"the secondary cell has {n:,} edges citing {list(pubs)} and "
            f"the amendment's table does not carry the row "
            f"'{claim}'")

    # EVERY PLATE CELL: the refusal is total disjointness being impossible,
    # so what matters is that one publication is shared by all of its edges.
    for (at, virus), counts in sorted(groups.items(),
                                      key=lambda r: -sum(r[1].values())):
        shared = set.intersection(*(set(p) for p in counts))
        assert shared, (
            f"{at} / {virus}: the edges share NO publication, so a "
            f"publication-disjoint reference set is no longer empty for "
            f"this cell and the amendment's claim that it still refuses is "
            f"false")
        for pmid in sorted(shared):
            assert pmid in prereg, (
                f"every edge of {at} / {virus} cites {pmid}, which is why "
                f"that cell refuses, and the amendment does not name it")
