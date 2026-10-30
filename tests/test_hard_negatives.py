"""Gate for the hard-negative driver: the two label sides are filtered by one
function, the metapath reach survives into the results, and the no-op control
is actually compared against something.
"""
import importlib.util
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _mod():
    spec = importlib.util.spec_from_file_location(
        "hard_negatives", ROOT / "scripts" / "hard_negatives.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


class _Labels:
    def __init__(self, pos, neg):
        self.positives, self.negatives = pos, neg


V = "NCBITaxon:2697049"


def _labels():
    return _Labels({V: {"A", "B", "C", "D"}}, {V: {"W", "X", "Y", "Z"}})


# si: A and W are selective, B and X have an index but are cytotoxic,
# C, D, Y, Z have no paired cytotoxicity data at all.
SI = {("A", V): 50.0, ("B", V): 2.0, ("W", V): 80.0, ("X", V): 1.5}


def test_default_filters_positives_only_which_is_not_like_for_like():
    """The published 0.826 is computed this way, and the asymmetry is real.

    Requiring a verified index of a POSITIVE requires paired same-document
    CC50/EC50, which selects for well-characterised compounds -- the ones most
    likely to carry the target annotations M1 reads. Measured on v0.1 the
    positive set halves 1,726 -> 808 and M1 rises 0.729 -> 0.823, while
    cov_pos goes 51.5% -> 70.3% and cov_neg stays at 5.3%.
    """
    m = _mod()
    lab = _labels()
    args = SimpleNamespace(selectivity="selective-only", si_threshold=10.0,
                           selectivity_applies="positives")
    m._filter_labels(lab, SI, args)
    assert lab.positives[V] == {"A"}          # only the selective active
    assert lab.negatives[V] == {"W", "X", "Y", "Z"}   # untouched


def test_both_applies_the_same_requirement_to_negatives():
    """The control. A compound must clear the same evidence bar to be a
    negative as to be a positive."""
    m = _mod()
    lab = _labels()
    args = SimpleNamespace(selectivity="selective-only", si_threshold=10.0,
                           selectivity_applies="both")
    m._filter_labels(lab, SI, args)
    assert lab.positives[V] == {"A"}
    assert lab.negatives[V] == {"W"}


def test_verified_only_keeps_cytotoxic_compounds_on_both_sides():
    """verified-only isolates 'has matched data' from 'passed the threshold',
    so a cytotoxic compound with an index stays in."""
    m = _mod()
    lab = _labels()
    args = SimpleNamespace(selectivity="verified-only", si_threshold=10.0,
                           selectivity_applies="both")
    m._filter_labels(lab, SI, args)
    assert lab.positives[V] == {"A", "B"}
    assert lab.negatives[V] == {"W", "X"}


def test_one_helper_serves_both_the_cross_sectional_and_temporal_sides():
    """The temporal block held a second hand-written copy of the loop, which
    is how the two sides came to differ without anyone deciding they should."""
    src = (ROOT / "scripts" / "hard_negatives.py").read_text()
    assert src.count("_filter_labels(") >= 3      # definition + both call sites
    assert "si_t.get(" not in src                 # the duplicated loop is gone


def test_cross_sectional_reach_is_not_discarded_before_it_is_written():
    """Batch 2a reset cross_reach to {} one line AFTER run() returned it, so
    the cross-sectional metapath reach never reached the results file -- the
    exact thing 2a was written to record."""
    src = (ROOT / "scripts" / "hard_negatives.py").read_text()
    assigned = src.index("cross, cross_reach = run(")
    reinit = src.find("cross_reach: dict = {}")
    assert reinit == -1 or reinit < assigned, \
        "cross_reach is re-initialised after run() set it"


# --------------------------------------------------------- the no-op control
#
# docs/evaluation.md states that `--selectivity all --selectivity-applies
# both` "reproduces the row-2 figure, so the flag does nothing it should
# not." The Reproducing block ran it and NOTHING COMPARED THE RESULT. Every
# hard_negatives.py variant writes to one hard_negatives.json, so by the
# time the control ran, row 2's file had been overwritten two steps earlier
# -- and the control's own output was overwritten one step later by row 5.
# About two and a half minutes per rebuild producing a file nobody read,
# reported in the document in the past tense as a check that had passed.
#
# It runs here instead, on every commit rather than once per rebuild, on a
# fixture small enough to be cheap and large enough to answer the question.
#
# TWO ARMS, and the second one is the point. That two documents agree proves
# nothing unless the fixture could have shown them disagreeing: the smoke
# release yields `cross_sectional: []` and would pass the agreement arm by
# comparing two empty documents, which is the fixture-too-small failure this
# suite keeps finding. So the same pair of settings must DIVERGE under
# `--selectivity selective-only`, where the filter actually runs.

VIRUS = "NCBITaxon:2697049"
VIRUS_LABEL = "Test virus"
PROTEINS = [f"UniProtKB:P{i}" for i in range(1, 5)]

#: 12 selective actives and 4 cytotoxic ones; 12 selective inactives and 6
#: with no verified index. Every arm leaves both classes at or above
#: MIN_CLASS_SIZE, so no arm is skipped for want of labels -- a skip would
#: make the two documents agree by both being empty.
COMPOUNDS: list[tuple[str, float | None]] = (
    [(f"INCHIKEY:A{i:02d}", 50.0) for i in range(1, 13)]
    + [(f"INCHIKEY:A{i:02d}", 2.0) for i in range(13, 17)]
    + [(f"INCHIKEY:N{i:02d}", 80.0) for i in range(1, 13)]
    + [(f"INCHIKEY:N{i:02d}", None) for i in range(13, 19)]
)


def _edge(subject, predicate, obj, qualifiers):
    return {"subject": subject, "predicate": predicate, "object": obj,
            "qualifiers": qualifiers, "primary_knowledge_source": "infores:test",
            "evidence_tier": 1, "first_asserted_date": "2020-06-01"}


def _build_release(out: Path) -> Path:
    nodes = [{"id": VIRUS, "class": "OrganismTaxon",
              "properties": {"label": VIRUS_LABEL}}]
    nodes += [{"id": p, "class": "Protein",
               "properties": {"taxon_id": "NCBITaxon:9606", "is_viral": False}}
              for p in PROTEINS]
    nodes += [{"id": d, "class": "SmallMolecule", "properties": {}}
              for d, _ in COMPOUNDS]

    # Two of the four proteins are host dependencies, so M2 reaches some
    # compounds and not others. A channel that reaches everything or nothing
    # would make coverage constant and the comparison less sensitive.
    edges = [_edge(p, "HOST_FACTOR_FOR", VIRUS, {"direction": "dependency"})
             for p in PROTEINS[:2]]
    for k, (drug, si) in enumerate(COMPOUNDS):
        targets = [PROTEINS[k % 4]]
        if k % 3 == 0:                      # varied degree: the degree
            targets.append(PROTEINS[(k + 1) % 4])   # baseline is not constant
        edges += [_edge(drug, "TARGETS", p, {}) for p in targets]
        active = drug.startswith("INCHIKEY:A")
        q = {"assay_type": "cell_based_antiviral",
             "ec50_nm": 200.0 if active else 50_000.0,
             "relation": "=" if active else ">"}
        if si is not None:
            q["selectivity_index"] = si
        edges.append(_edge(drug, "HAS_ANTIVIRAL_ACTIVITY_AGAINST", VIRUS, q))

    out.mkdir(parents=True, exist_ok=True)
    (out / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n, sort_keys=True) for n in nodes))
    (out / "edges.jsonl").write_text(
        "\n".join(json.dumps(e, sort_keys=True) for e in edges))
    return out


def _run(release: Path, results: Path, selectivity: str, applies: str) -> dict:
    """The real entry point, as the block invokes it -- not main() imported.

    main() calls parse_args() with no argv, so importing it would test a
    different surface from the one the document runs.
    """
    proc = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "hard_negatives.py"),
         "--release", str(release),
         # Absent on purpose: the temporal column needs a stamped training
         # graph, both arms are denied one identically, and the skip and its
         # reason are recorded in the document being compared.
         "--train", str(release.parent / "no-training-graph"),
         "--results", str(results),
         "--selectivity", selectivity,
         "--selectivity-applies", applies],
        cwd=ROOT, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    doc = json.loads((results / "hard_negatives.json").read_text())
    # The one field that moves between two runs of the SAME configuration.
    # Dropped by name rather than by dropping `provenance`, so a difference
    # in release fingerprint or code commit would still show up.
    doc["provenance"].pop("generated_utc", None)
    return doc


@pytest.fixture(scope="module")
def arms(tmp_path_factory) -> dict[tuple[str, str], dict]:
    release = _build_release(tmp_path_factory.mktemp("hn_release") / "release")
    out = tmp_path_factory.mktemp("hn_results")
    return {(f, a): _run(release, out / f"{f}_{a}", f, a)
            for f in ("all", "selective-only")
            for a in ("positives", "both")}


def test_the_fixture_actually_evaluates_something(arms):
    """GUARDS EVERY COMPARISON BELOW. If the virus were skipped for want of
    labels, or no channel reached a compound, the documents would agree by
    recording nothing -- a pass that means the opposite of what it looks
    like."""
    per_virus = arms[("all", "positives")]["cross_sectional"]
    assert list(per_virus) == [VIRUS_LABEL], per_virus
    per = per_virus[VIRUS_LABEL]
    assert {"M2", "COMBINED", "degree"} <= set(per), sorted(per)
    floor = _mod().MIN_CLASS_SIZE
    assert per["degree"]["n_pos"] >= floor
    assert per["degree"]["n_neg"] >= floor
    assert 0 < per["M2"]["coverage_pos"] < 1, "M2 reaches all or none"


def test_the_fixture_can_tell_the_two_settings_apart(arms):
    """THE CONTROL FOR THE CONTROL. Where the filter actually runs,
    --selectivity-applies both must change the evaluation. Without this arm
    the agreement below would be satisfied by any fixture too small to
    disagree."""
    a = arms[("selective-only", "positives")]
    b = arms[("selective-only", "both")]
    assert a["cross_sectional"] != b["cross_sectional"]
    pa = a["cross_sectional"][VIRUS_LABEL]["degree"]
    pb = b["cross_sectional"][VIRUS_LABEL]["degree"]
    assert pa["n_pos"] == pb["n_pos"], "the positives side is not the one 'both' adds"
    assert pb["n_neg"] < pa["n_neg"], "the negatives were never filtered"


def test_selectivity_applies_is_a_no_op_when_the_filter_is_all(arms):
    """The claim docs/evaluation.md makes. `all` short-circuits the filter
    block entirely, so the flag is recorded in the output and otherwise
    unread."""
    a, b = arms[("all", "positives")], arms[("all", "both")]
    differing = {k for k in set(a) | set(b) if a.get(k) != b.get(k)}
    # EXACTLY the recorded flag, not 'no evaluation difference'. Requiring
    # that one field to differ is what proves the two subprocesses really ran
    # different configurations; without it this assertion would be satisfied
    # just as well by comparing a document with itself.
    assert differing == {"selectivity_applies"}, sorted(differing)
    assert (a["selectivity_applies"], b["selectivity_applies"]) \
        == ("positives", "both")
