"""The entry-assay registration and its two amendments, against the release.

WHY THIS FILE EXISTS. doc2/registrations/prereg_ncats_entry_assays.md was signed before
the ingest code existed, and both of its amendments record a prediction
that FAILED. A registration whose failures are recorded by hand is a
registration that can be quietly improved later, and the whole point of
signing one is that it cannot be.

So every figure the second amendment states is re-derived here from the
release and the project's own functions: the cell sizes, the pooled figure
the "no pooling" decision forbids, the disjointness that justifies refusing
a cell, the plate family size, the corrected z, the detectable floors, and
the compound-by-compound accounting of why P1c missed by five.

WHAT THIS FILE DELIBERATELY DOES NOT DO. It computes no AUC. The
registration is unsigned-in-spirit the moment a test asserts an outcome.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from statistics import NormalDist

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

PREREG = ROOT / "doc2/registrations/prereg_ncats_entry_assays.md"
SOURCE = "infores:ncats-opendata"
ARCHIVE = ROOT / "data/archive/v0.1-pre-entry-assays-2026-11-17/release"
SNAPSHOT = ROOT / "doc2/records/release_operands_2026-11-17_pre_entry_assays.json"


@pytest.fixture(scope="module")
def prereg() -> str:
    assert PREREG.exists(), "the registration is missing"
    return " ".join(PREREG.read_text().split())


@pytest.fixture(scope="module")
def amendment2() -> str:
    """The SECOND amendment alone.

    Scoped because the document deliberately keeps the superseded family of
    four in place -- a registration edited to match what happened is worth
    nothing -- so a search over the whole text finds 0.05/4 first and the
    check passes against the figure it was written to catch.
    """
    text = PREREG.read_text()
    marker = "# Second amendment,"
    assert marker in text, "the second amendment is gone"
    return " ".join(text[text.index(marker):].split())


@pytest.fixture(scope="module")
def release() -> Path:
    from kgav.pipeline import release_dir
    rd = release_dir(ROOT)
    if not (rd / "edges.jsonl").exists():
        pytest.skip("no assembled release on disk")
    return rd


@pytest.fixture(scope="module")
def cells(release):
    """Every plate cell, keyed by the protocol that produced it."""
    from kgav.labels import build_labels
    import structure_scorer as S

    out = {}
    for key, source, assay_types in S.POPULATIONS:
        if assay_types is None:          # the ChEMBL population is not a plate
            continue
        labels = build_labels(release, sources={source},
                              assay_types=assay_types)
        refused = S.complement_cells(release, source, assay_types, labels)
        for virus in sorted(set(labels.positives) | set(labels.negatives)):
            out[(key, virus)] = {
                "pos": set(labels.positives.get(virus, set())),
                "neg": set(labels.negatives.get(virus, set())),
                "refused": virus in refused,
            }
    return out


def test_every_cell_size_in_the_second_amendment_is_the_releases(prereg,
                                                                 cells):
    """Four rows of positives against negatives, including the refused one.

    Typed by hand these drift, and the first amendment's own P1c drifted by
    five within a day of being written."""
    checked = 0
    for (key, virus), c in sorted(cells.items()):
        pair = f"{len(c['pos']):,} | {len(c['neg']):,}"
        stripped = f"~~{len(c['pos']):,}~~ | ~~{len(c['neg']):,}~~"
        assert pair in prereg or stripped in prereg, (
            f"{key} {virus} is {len(c['pos']):,} / {len(c['neg']):,} on this "
            f"release and the amendment's table carries neither that pair "
            f"nor a struck-through form of it")
        checked += 1
    assert checked == 4, f"{checked} plate cells, expected 4"


def test_the_pooled_figure_the_no_pooling_decision_forbids_is_quoted(prereg,
                                                                     release):
    """The amendment claims that without the assay_type filter the cell is
    the pooled one. That is a property of the release, so it is derived --
    and it is the figure a reader needs in order to see what the filter
    bought."""
    from kgav.labels import build_labels
    pooled = build_labels(release, sources={SOURCE})
    merged = [(v, len(pooled.positives.get(v, set())),
               len(pooled.negatives.get(v, set())))
              for v in sorted(set(pooled.positives) | set(pooled.negatives))]
    # the virus screened under two protocols is the one whose pooled cell
    # exceeds either of its per-protocol cells -- derived, not named
    biggest = max(merged, key=lambda r: r[1])
    assert f"{biggest[1]:,} / {biggest[2]:,}" in prereg, (
        f"pooled, this release gives {biggest[1]:,} / {biggest[2]:,} for "
        f"{biggest[0]} and the amendment does not quote it")


def test_the_refused_cell_shares_no_compound_with_the_protocol_that_won(
        prereg, cells, release):
    """THE REFUSAL'S WHOLE JUSTIFICATION, and it is a set property.

    If the refused cell shared even one compound with the other protocol's
    cell for the same virus, "complement" would be the wrong word and the
    refusal would need a different reason. Zero is not "few"."""
    refused = [(k, c) for k, c in cells.items() if c["refused"]]
    assert len(refused) == 1, (
        f"{len(refused)} cells refused as complements; the amendment "
        f"registers exactly one")
    (key, virus), cell = refused[0]
    mine = cell["pos"] | cell["neg"]
    siblings = [c for (k, v), c in cells.items()
                if v == virus and k != key]
    assert siblings, f"nothing for {virus} to be the complement OF"
    theirs = set().union(*(s["pos"] | s["neg"] for s in siblings))
    assert not (mine & theirs), (
        f"{len(mine & theirs)} compounds are in both protocols for {virus}; "
        f"the refusal is stated as total disjointness and is not")
    assert f"0 of {len(mine):,}" in prereg, (
        f"the overlap is 0 of {len(mine):,} and the amendment does not say "
        f"so in those terms")
    assert len(theirs) > len(mine), (
        "the refused side must be the smaller one, or the rule refuses the "
        "wrong cell")

    # AND THE REFUSAL MUST BE DOING WORK. If this cell were too thin to
    # score, the size gate would already have excluded it and the
    # complement rule would be decoration -- which is the shape of a
    # check that cannot fail.
    from kgav.calibration import MIN_REACHED_POSITIVES
    assert (len(cell["pos"]) >= MIN_REACHED_POSITIVES
            and len(cell["neg"]) >= MIN_REACHED_POSITIVES), (
        f"the refused cell is {len(cell['pos'])}/{len(cell['neg'])} and "
        f"would have been excluded by the size gate anyway; the complement "
        f"refusal is then untested")


def test_the_plate_family_is_the_cells_that_survived_and_not_the_signed_four(
        amendment2, cells):
    """The document was signed for four cells at 0.05/4. Three exist. The
    family size is therefore a derived count, and the alpha must follow it
    rather than the signature."""
    from kgav.calibration import MIN_REACHED_POSITIVES
    import structure_scorer as S

    readable = [k for k, c in cells.items()
                if not c["refused"]
                and len(c["pos"]) >= MIN_REACHED_POSITIVES
                and len(c["neg"]) >= MIN_REACHED_POSITIVES]
    assert len(readable) == S.PLATE_FAMILY, (
        f"{len(readable)} plate cells are readable and the scorer fixes "
        f"PLATE_FAMILY at {S.PLATE_FAMILY}")
    m = re.search(r"alpha = 0\.05/(\d+) = \*{0,2}([\d.]+)\*{0,2}, "
                  r"\*{0,2}z = \*{0,2}([\d.]+)", amendment2)
    assert m, "the amendment no longer states its corrected alpha and z"
    assert int(m.group(1)) == len(readable), (
        f"the amendment divides by {m.group(1)} and {len(readable)} cells "
        f"survive")
    assert abs(float(m.group(2)) - 0.05 / len(readable)) < 5e-7
    expected_z = NormalDist().inv_cdf(1 - (0.05 / len(readable)) / 2)
    assert abs(float(m.group(3)) - expected_z) < 5e-5, (
        f"two-sided z for alpha {0.05 / len(readable):.6f} is "
        f"{expected_z:.4f}, not {m.group(3)}")
    assert abs(S.PLATE_Z_CORRECTED - expected_z) < 5e-5


def test_every_detectable_floor_in_the_amendment_is_recomputed(prereg, cells):
    """Both columns, because the amendment states both and says they
    disagree. A document carrying two methods must carry them right."""
    from diagnose_structure_scorer_design import detectable_auc
    from kgav.labels import auc_interval
    import structure_scorer as S

    z = S.PLATE_Z_CORRECTED
    for (key, virus), c in sorted(cells.items()):
        p, n = len(c["pos"]), len(c["neg"])
        lo, _hi = auc_interval(0.5, p, n, z=z)
        at_half = round(0.5 + (0.5 - lo), 4)
        assert f"{at_half:.4f}" in prereg, (
            f"{key} {virus} ({p}/{n}) has a half-width floor of "
            f"{at_half:.4f} at z={z} and the amendment does not state it")
        if c["refused"]:
            continue
        walked = detectable_auc(p, n, z)
        assert walked is not None and f"{walked:.4f}" in prereg, (
            f"{key} {virus} has detectable_auc {walked} and the amendment "
            f"does not state it")


def test_the_chembl_cells_did_not_move_and_the_amendment_counts_them(
        prereg, release):
    """H2's label route, which the amendment reports as closed. Closed by
    comparison against the committed snapshot, not by assertion."""
    import hashlib
    from kgav.labels import build_labels
    if not SNAPSHOT.exists():
        pytest.skip("no pre-entry snapshot committed")
    before = json.loads(SNAPSHOT.read_text())["label_cells"]["infores:chembl"]
    now = build_labels(release, sources={"infores:chembl"})

    def digest(ids) -> str:
        return hashlib.sha256("\n".join(sorted(ids)).encode()).hexdigest()

    same, moved = 0, []
    tot_pos = tot_neg = 0
    for virus, snap in before.items():
        pos = now.positives.get(virus, set())
        neg = now.negatives.get(virus, set())
        tot_pos += len(pos)
        tot_neg += len(neg)
        if (digest(pos) == snap["pos_sha256"]
                and digest(neg) == snap["neg_sha256"]):
            same += 1
        else:
            moved.append(virus)
    assert f"**{same}** ChEMBL cells" in prereg, (
        f"{same} cells are byte-identical and the amendment does not say "
        f"{same}")
    assert f"{tot_pos:,} positives and {tot_neg:,} negatives" in prereg
    if moved:
        # NOT A FAILURE BY ITSELF -- a finding the amendment must carry.
        for virus in moved:
            assert virus in prereg, (
                f"{virus}'s ChEMBL cell moved in the rebuild and the "
                f"amendment reports the label route as unmoved")


def test_the_p1c_accounting_adds_up_against_the_archived_release(prereg,
                                                                 release):
    """FIVE, AND WHY FIVE. The first amendment predicted ten negatives would
    leave the CPE cell and five did. The second amendment explains the gap
    as 4 weak-hit removals + 6 structures that kept an independent inactive
    row + 1 edge lost to the merge, and every one of those is countable."""
    from kgav.labels import (ACTIVITY_PREDICATE, DEFAULT_ACTIVE_BELOW_NM,
                             build_labels, edge_sources)
    if not (ARCHIVE / "edges.jsonl").exists():
        pytest.skip("no archived pre-entry release on disk")

    old = build_labels(ARCHIVE, sources={SOURCE})
    new = build_labels(release, sources={SOURCE},
                       assay_types={"cell_based_antiviral"})
    virus = max(old.positives, key=lambda v: len(old.negatives.get(v, set())))
    o_neg = set(old.negatives[virus])
    n_neg = set(new.negatives[virus])
    left, joined = o_neg - n_neg, n_neg - o_neg
    assert not joined, f"{len(joined)} compounds JOINED the negative class"

    # the weak hits, read off the archived edges rather than named
    weak = set()
    for line in (ARCHIVE / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        if e["predicate"] != ACTIVITY_PREDICATE or e["object"] != virus:
            continue
        if SOURCE not in edge_sources(e):
            continue
        q = e.get("qualifiers") or {}
        if q.get("relation") == "=" and (q.get("ec50_nm") or 0) > \
                DEFAULT_ACTIVE_BELOW_NM:
            weak.add(e["subject"])

    gone = weak - n_neg
    kept = weak & n_neg
    displaced = left - weak
    assert len(gone) + len(displaced) == len(left)
    for n, phrase in ((len(gone), "**%d** of the ten weak-hit"),
                      (len(kept), "**%d** remain negatives"),
                      (len(displaced), "**%d** further compound")):
        assert (phrase % n) in prereg, (
            f"the accounting is {len(gone)} gone, {len(kept)} kept, "
            f"{len(displaced)} displaced, and the amendment does not state "
            f"{n} where it says '{phrase % n}'")
    assert f"{len(new.positives[virus]):,} / {len(n_neg):,}" in prereg
