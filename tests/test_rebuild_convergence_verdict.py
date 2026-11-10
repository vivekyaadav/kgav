"""The rebuild must not report convergence it has not measured (L20).

WHY. `make rebuild` assembles at phase 2 from the PREVIOUS run's phase-3
directories, so what `ingest_similarity` reads at phase 3 contains `gtop`'s
output from the last rebuild -- and `gtop` is also phase 3, contributing
TARGETS. `ingest_similarity` partitions compounds into query and reference
sets by their TARGETS / INHIBITS / MEASURED_INACTIVE_AGAINST subjects, so
the two layers are mutually dependent through the assembly between them and
one pass of 3-then-4 is not a fixed point whenever that subject set moves.

MEASURED, NOT ARGUED. Reverting fix-list B1 -- a phase-1 change -- went

    af01e12ce8a86891 -> 8402bbaf9ed5c531     one extra CHEMICALLY_SIMILAR_TO
                     -> f8a1ae20bfd89b37     converged, byte-identical

and `rebuild_release.py` printed "The release CHANGED this run" BOTH TIMES.
`staleness_verdict` had four branches -- no release, none after, unchanged,
changed -- and no fifth for "changed because it has not converged yet",
which is the one an operator most needs told apart from the others.

This is L20 option (c): say so, and ask for another run. Options (a) loop
until the fingerprint repeats and (b) build phase 3 from the phase-1
directories remain open, and this test does not assert either.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import rebuild_release as R  # noqa: E402
from kgav.similarity import ANNOTATING, annotating_subjects_sha  # noqa: E402

FP_A = {"id": "v0.1", "nodes_sha": "n1", "edges_sha": "af01e12ce8a86891"}
FP_B = {"id": "v0.1", "nodes_sha": "n2", "edges_sha": "8402bbaf9ed5c531"}


def test_an_unchanged_release_still_reports_convergence():
    """The branch that matters most often: a no-op rebuild. It must not
    acquire a scary new caveat, which is the risk in adding a fifth case."""
    out = R.staleness_verdict(FP_A, FP_A, "aaaa", "aaaa")
    assert "did not change the release" in out
    assert "FIXED POINT" not in out


def test_a_changed_release_with_a_settled_partition_says_only_stale():
    """Nine of ten rebuilds in this repository's history: the graph moved,
    but not the set similarity partitions on, so phase 3 is settled."""
    out = R.staleness_verdict(FP_A, FP_B, "aaaa", "aaaa")
    assert "The release CHANGED this run" in out
    assert "MAY NOT BE THE FIXED POINT" not in out


def test_a_moved_annotating_set_demands_another_run():
    """THE B1 REVERT, in the shape the run actually had. This is the case
    that printed the same sentence as a converged run."""
    out = R.staleness_verdict(FP_A, FP_B, "aaaa", "bbbb")
    assert "The release CHANGED this run" in out
    assert "MAY NOT BE THE FIXED POINT -- RUN THE REBUILD AGAIN" in out
    assert "aaaa -> bbbb" in out
    assert "TARGETS / INHIBITS / MEASURED_INACTIVE_AGAINST" in out
    assert "L20" in out


def test_the_branch_is_silent_when_the_digests_are_unknown():
    """A caller that cannot measure the set -- no release on disk either
    side -- must get the old four-branch behaviour, not a warning derived
    from two Nones."""
    assert "FIXED POINT" not in R.staleness_verdict(FP_A, FP_B, None, None)
    assert "FIXED POINT" not in R.staleness_verdict(FP_A, FP_B, "aaaa", None)
    assert "FIXED POINT" not in R.staleness_verdict(FP_A, FP_B, None, "bbbb")


def test_the_predicate_list_has_one_owner():
    """ingest_similarity.py used to define this tuple as a local inside
    main(). The rebuild now asks the same question, and a second copy would
    be a second definition of the feedback loop -- the defect that produced
    L20 in the first place was a docstring disagreeing with the code."""
    src = (ROOT / "scripts/ingest_similarity.py").read_text()
    assert "ANNOTATING = (" not in src, (
        "ingest_similarity.py defines ANNOTATING again; import it from "
        "kgav.similarity instead")
    assert "ANNOTATING" in src, "the ingest must still use the tuple"
    assert ANNOTATING == ("TARGETS", "INHIBITS", "MEASURED_INACTIVE_AGAINST")


def test_the_digest_is_over_subjects_and_ignores_edge_count(tmp_path):
    """It must track WHICH compounds are annotated, not how many edges say
    so: a second TARGETS edge for a compound already annotated does not
    move similarity's partition and must not demand a rerun."""
    import json
    a = tmp_path / "a"; a.mkdir()
    b = tmp_path / "b"; b.mkdir()
    rows_a = [{"subject": "C1", "predicate": "TARGETS", "object": "P1"}]
    rows_b = rows_a + [{"subject": "C1", "predicate": "TARGETS", "object": "P2"}]
    (a / "edges.jsonl").write_text("".join(json.dumps(r) + "\n" for r in rows_a))
    (b / "edges.jsonl").write_text("".join(json.dumps(r) + "\n" for r in rows_b))
    assert annotating_subjects_sha(a) == annotating_subjects_sha(b)

    c = tmp_path / "c"; c.mkdir()
    rows_c = rows_a + [{"subject": "C2", "predicate": "TARGETS", "object": "P1"}]
    (c / "edges.jsonl").write_text("".join(json.dumps(r) + "\n" for r in rows_c))
    assert annotating_subjects_sha(a) != annotating_subjects_sha(c)


def test_a_predicate_outside_ANNOTATING_does_not_move_the_digest(tmp_path):
    """PARTICIPATES_IN is 157,211 of the 392,363 edges and has nothing to
    do with the partition. If it moved the digest, nearly every rebuild
    would demand a second run and the warning would be ignored."""
    import json
    a = tmp_path / "a"; a.mkdir()
    b = tmp_path / "b"; b.mkdir()
    base = [{"subject": "C1", "predicate": "TARGETS", "object": "P1"}]
    (a / "edges.jsonl").write_text("".join(json.dumps(r) + "\n" for r in base))
    (b / "edges.jsonl").write_text("".join(
        json.dumps(r) + "\n" for r in base
        + [{"subject": "P9", "predicate": "PARTICIPATES_IN", "object": "PW1"}]))
    assert annotating_subjects_sha(a) == annotating_subjects_sha(b)


def test_the_digest_is_none_without_a_release(tmp_path):
    assert annotating_subjects_sha(tmp_path / "absent") is None
