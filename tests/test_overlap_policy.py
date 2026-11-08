"""The cross-cutoff overlap policy, re-derived (fix-list C8).

TWO DEFECTS, ONE SENTENCE APART. `scripts/hard_negatives.py` built its
temporal test set with `labels.build_labels(year_from=cutoff+1)` and applied
no overlap policy at all, while the results file it wrote reported
`train_split: {"overlap_policy": "drop-from-test"}` -- a field copied from
the training graph's stamp, which records how `temporal_split.py` cut the
GRAPH and says nothing about these labels. The file stated a policy no code
path applied, which a reader cannot distinguish from one that was applied.

And `temporal.build_split`'s comment justified the policy with a measurement
that is false on this release: "All 19 agreed on polarity across both
sources, which is why dropping them costs no information". Re-derived here,
the overlap is larger than 19 and nearly half of it DISAGREES across the
cutoff -- in one direction, inactive before and active after.

The code was right and its stated reason was not, which is the harder of
the two to notice: `drop_already_known` matches positives against train
positives and negatives against train negatives, so a compound whose answer
CHANGED survives the filter. Those survivors are the hardest prospective
cases the split contains. Whether they should survive is a protocol
question and not one to settle after seeing the AUCs, so this file pins the
behaviour and the census, not a preference.

Every figure in that comment is parsed out of the comment and compared with
one derived from the release, for the reason test_landscape.py states: a
check written with the same literals it is checking cannot survive the
event it exists for.
"""
from __future__ import annotations

import inspect
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kgav import temporal  # noqa: E402
from kgav.labels import build_labels  # noqa: E402
from kgav.temporal import drop_already_known  # noqa: E402

RELEASE = ROOT / "data" / "releases" / "v0.1"
CUTOFF = 2021


@pytest.fixture(scope="module")
def sides():
    if not (RELEASE / "edges.jsonl").exists():
        pytest.skip("release v0.1 is not built in this checkout")
    return (build_labels(RELEASE, year_to=CUTOFF),
            build_labels(RELEASE, year_from=CUTOFF + 1))


@pytest.fixture(scope="module")
def census(sides):
    """(virus, compound) -> (polarity before, polarity after).

    A compound can appear in both the positive and the negative set of one
    side -- different cell lines, strains, labs -- so a side's polarity is a
    SET, not a flag. labels.build_labels keeps that disagreement rather than
    forcing it to one side, and collapsing it here would invent the very
    label the loader refuses to invent.
    """
    pre, post = sides
    out: dict[tuple[str, str], tuple[frozenset, frozenset]] = {}
    for virus in sorted(set(post.positives) | set(post.negatives)):
        before = {"active": pre.positives.get(virus, set()),
                  "inactive": pre.negatives.get(virus, set())}
        after = {"active": post.positives.get(virus, set()),
                 "inactive": post.negatives.get(virus, set())}
        seen = before["active"] | before["inactive"]
        for compound in sorted(after["active"] | after["inactive"]):
            if compound not in seen:
                continue
            out[(virus, compound)] = (
                frozenset(k for k, v in before.items() if compound in v),
                frozenset(k for k, v in after.items() if compound in v))
    return out


def agreeing(census):
    return {k for k, (b, a) in census.items() if b == a}


def flipping(census):
    return {k for k, (b, a) in census.items() if b != a}


def test_the_overlap_is_not_empty(census):
    """Everything below is vacuous if nothing crosses the cutoff, which is
    how a policy test passes while the policy does nothing."""
    assert census, (
        f"no compound in the post-{CUTOFF} test set was measured before "
        f"{CUTOFF}; this file would then pass without testing anything")


def test_the_filter_drops_the_repeats_and_keeps_the_changed_answers(sides, census):
    """The property the comment's reasoning rests on, stated as a property.

    Not "it drops 18": it drops EXACTLY the pairs whose answer did not
    change, and keeps exactly the ones whose answer did.
    """
    pre, post = sides
    test_pos = {v: set(s) for v, s in post.positives.items()}
    test_neg = {v: set(s) for v, s in post.negatives.items()}
    before = {(v, c) for v, s in test_pos.items() for c in s}
    before |= {(v, c) for v, s in test_neg.items() for c in s}

    drop_already_known(test_pos, test_neg, pre.positives, pre.negatives, CUTOFF)

    after = {(v, c) for v, s in test_pos.items() for c in s}
    after |= {(v, c) for v, s in test_neg.items() for c in s}
    dropped = before - after

    # An agreeing pair may hold BOTH polarities on both sides, in which case
    # each polarity is dropped on its own side and the pair leaves entirely.
    assert dropped == agreeing(census), (
        f"dropped but did not agree: {sorted(dropped - agreeing(census))}; "
        f"agreed but was not dropped: {sorted(agreeing(census) - dropped)}")
    assert flipping(census) <= after, (
        f"compounds whose answer changed across {CUTOFF} were dropped: "
        f"{sorted(flipping(census) - after)}")


def test_build_split_drops_what_the_helper_drops(census):
    """The caller that cuts the graph, measured against the same census.

    build_split re-reads the release and partitions the edges, so it can
    disagree with the helper for reasons the helper cannot see -- the
    ambiguity rule runs first there, and the undated policy decides which
    edges exist at all. This is the check that it does not.
    """
    import kgav.schema as schema_mod
    sch = schema_mod.Schema.load(ROOT / "schema" / "kg_schema.yaml")
    split = temporal.build_split(
        RELEASE, CUTOFF, set(sch.held_out_predicates()),
        undated_policy="include",
        label_predicates=set(sch.evaluation_label_predicates()))
    dropped = (split.stats.get("test_positive_already_known_dropped", 0)
               + split.stats.get("test_negative_already_known_dropped", 0))
    # Per POLARITY, not per pair: a pair holding both polarities on both
    # sides is counted once on each side.
    expected = sum(len(b & a) for b, a in
                   (census[k] for k in agreeing(census)))
    assert dropped == expected, (
        f"build_split dropped {dropped} already-known label(s) and the "
        f"census derived from build_labels implies {expected}")


def test_the_two_policies_are_distinguishable_on_this_release(sides, census):
    """`overlap_policy` is recorded in the results file, and a recorded
    setting that changes nothing is a setting nobody can audit. This asserts
    the two values actually produce different test sets here, so the field
    in hard_negatives.json carries information."""
    pre, post = sides
    test_pos = {v: set(s) for v, s in post.positives.items()}
    test_neg = {v: set(s) for v, s in post.negatives.items()}
    stats, messages = drop_already_known(
        test_pos, test_neg, pre.positives, pre.negatives, CUTOFF)
    assert sum(stats.values()) > 0 and messages, (
        "the helper reported nothing to drop, so a results file recording "
        "either policy would look the same")


def _comment() -> str:
    """build_split's source with comment markers stripped and runs of
    whitespace collapsed, so a claim spanning six `#` lines reads as one
    sentence."""
    src = inspect.getsource(temporal.build_split)
    flat = " ".join(re.sub(r"^\s*#\s?", "", ln) for ln in src.splitlines())
    return " ".join(flat.split())


def test_the_comments_own_figures_are_the_releases(census):
    """The three numbers the comment states, parsed out of the comment.

    THE POINT OF PARSING THEM. The previous version of that comment stated
    19 and "all agreed"; both were wrong and nothing could tell, because the
    only record of the measurement was the sentence making the claim. Now a
    release that moves these counts fails here and names the comment.
    """
    text = _comment()
    m = re.search(r"cutoff 2021: (\d+) \(virus, compound\) pairs overlap the "
                  r"cutoff, (\d+) agree on polarity and (\d+) DISAGREE", text)
    assert m, ("build_split's comment no longer states the overlap census in "
               "the form this test reads; update both together")
    stated = tuple(int(g) for g in m.groups())
    derived = (len(census), len(agreeing(census)), len(flipping(census)))
    assert stated == derived, (
        f"the comment states {stated} (overlap, agreeing, flipping) and the "
        f"release gives {derived}")
    share = re.search(r"false for (\d+)% of the overlap", text)
    assert share, "the comment no longer states the share that disagrees"
    assert int(share.group(1)) == round(100 * derived[2] / derived[0]), (
        f"the comment says {share.group(1)}% of the overlap disagrees and "
        f"the release gives {round(100 * derived[2] / derived[0])}%")


def test_the_direction_of_every_disagreement_is_the_one_stated(census):
    """The comment says the disagreements run one way only. If a release
    ever produces an active-then-inactive pair, dropping it is a different
    decision -- a withdrawn activity, not a missed one -- and the sentence
    has to be rewritten rather than kept."""
    wrong_way = {k: (sorted(b), sorted(a)) for k, (b, a) in census.items()
                 if b != a and "active" not in a}
    text = _comment()
    assert "inactive before and active after" in text
    assert not wrong_way, (
        f"the comment says every disagreement is inactive-then-active and "
        f"these are not: {wrong_way}")


def test_the_false_justification_is_not_reinstated():
    """A WITHDRAWN CLAIM MAY BE QUOTED; IT MAY NOT BE ASSERTED. The comment
    quotes the old sentence in order to retract it, so a flat substring
    search calls the retraction the defect -- which is what it did in
    test_temporal_table.py and test_release_plan_figures.py before this.
    Quoted spans are stripped, and the retraction is then required to
    survive."""
    text = _comment()
    unquoted = re.sub(r'"[^"]*"', " ", text)
    assert "costs no information" not in unquoted, (
        "build_split's comment asserts, outside a quotation, that dropping "
        "the overlap costs no information; it does not")
    assert "costs no information" in text, (
        "the retraction has been deleted rather than the claim corrected; "
        "the comment should record what it used to assert")


DOC = ROOT / "docs" / "evaluation.md"
MARKER = "<!-- FIGURE: overlap_census -->"


def _paragraph() -> str:
    text = DOC.read_text()
    assert MARKER in text, (
        "the overlap census lost its machine-readable marker, so nothing "
        "re-derives the figures in that paragraph")
    return " ".join(text.split(MARKER, 1)[1].split("\n\n")[0].split())


def test_the_documents_census_is_the_releases(census):
    """`evaluation.md` stated 19 and "19/19 by source"; the release says
    something else. Parsed from the paragraph so a release that moves these
    counts fails here rather than in a reader's head."""
    para = _paragraph()
    m = re.search(r"the overlap is \*\*(\d+)\*\* pairs.*?\*\*(\d+)\*\* agree "
                  r"on polarity.*?\*\*(\d+) disagree", para)
    assert m, ("the census paragraph no longer states the three counts in "
               "the form this test reads; update both together")
    stated = tuple(int(g) for g in m.groups())
    derived = (len(census), len(agreeing(census)), len(flipping(census)))
    assert stated == derived, (
        f"evaluation.md states {stated} (overlap, agreeing, disagreeing) "
        f"and the release gives {derived}")


def test_the_documents_label_split_is_the_releases(sides, census):
    """The "5 positive and 13 negative labels" half of the sentence, which
    counts LABELS and not pairs -- a pair holding both polarities on both
    sides is dropped once on each side, so the two totals need not agree
    and the document must not conflate them."""
    para = _paragraph()
    m = re.search(r"as (\d+) positive and (\d+) negative labels", para)
    assert m, "the paragraph no longer states the per-polarity split"
    pre, post = sides
    test_pos = {v: set(s) for v, s in post.positives.items()}
    test_neg = {v: set(s) for v, s in post.negatives.items()}
    stats, _ = drop_already_known(test_pos, test_neg,
                                  pre.positives, pre.negatives, CUTOFF)
    derived = (stats.get("test_positive_already_known_dropped", 0),
               stats.get("test_negative_already_known_dropped", 0))
    assert tuple(int(g) for g in m.groups()) == derived, (
        f"evaluation.md states {m.groups()} positive/negative labels "
        f"dropped and the release gives {derived}")


def test_the_document_does_not_reinstate_the_old_census():
    """The retraction quotes the withdrawn sentence, so quoted spans are
    stripped first and the retraction is then required to survive."""
    flat = " ".join(DOC.read_text().split())
    unquoted = re.sub(r'"[^"]*"', " ", flat)
    assert "19/19 by source" not in unquoted, (
        "evaluation.md asserts the 19/19 split outside a quotation; the "
        "overlap is not 19 and it is not evenly split")
    assert "19/19 by source" in flat, (
        "the retraction of that figure has been deleted rather than the "
        "figure corrected; the document should record what it used to say")


RESULTS = ROOT / "data" / "results"


@pytest.mark.parametrize("name", ["hard_negatives.json",
                                  "hard_negatives_row1.json",
                                  "hard_negatives_row2.json",
                                  "hard_negatives_row3.json",
                                  "hard_negatives_row4.json"])
def test_each_results_file_records_the_policy_it_applied(name):
    """THE DEFECT THIS EXISTS FOR, and it is not a missing feature.

    These files have been reporting `train_split: {"overlap_policy":
    "drop-from-test"}` -- read off the TRAINING GRAPH's stamp, which
    describes how temporal_split.py cut the graph -- beside a temporal
    column built by labels.build_labels, which applied no overlap policy of
    any kind. The field was true of a different artifact, so the file looked
    audited and was not, and no reader could tell.

    `temporal_notes.overlap_policy` is what the column actually applied, and
    `already_known_dropped` is what it cost. A file whose temporal column
    ran without either is the old defect, so it fails here rather than
    skipping.
    """
    import json
    path = RESULTS / name
    if not path.exists():
        pytest.skip(f"{name} is not built in this checkout")
    doc = json.loads(path.read_text())
    if doc.get("temporal_skipped"):
        pytest.skip(f"{name}'s temporal column had no training graph")
    notes = doc.get("temporal_notes") or {}
    assert notes.get("overlap_policy"), (
        f"{name}'s temporal column ran and the file does not say which "
        f"overlap policy it applied, while train_split reports "
        f"{(doc.get('train_split') or {}).get('overlap_policy')!r} for the "
        f"graph it was cut from")
    assert "already_known_dropped" in notes, (
        f"{name} names a policy without recording what it dropped, which "
        f"cannot be distinguished from a policy that dropped nothing")
    if notes["overlap_policy"] == "drop-from-test" and not notes.get("empty_label_set"):
        assert sum(notes["already_known_dropped"].values()) > 0, (
            f"{name} applied drop-from-test to a non-empty test set and "
            f"dropped nothing; on this release the overlap is not empty")
