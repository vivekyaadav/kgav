"""Temporal split: train on what was known, test on what came next.

WHY THIS IS THE EVALUATION THAT MATTERS. Every metric so far is measured
against compounds people chose to screen, so it partly measures literature
recovery rather than prediction. The clearest symptom was M7's pool: a 59%
positive rate, meaning the compounds chemical similarity reached were
overwhelmingly ones already screened. A temporal split removes that by
construction -- a compound first assayed in 2023 was not in the 2021 graph's
pool at all.

CUTOFF. End of 2021, chosen from the label distribution rather than convention:
5,181 HAS_ANTIVIRAL_ACTIVITY_AGAINST edges are dated 2007-2021 and 2,027 are
dated 2022-2025, so both sides are large enough to measure. The 2021 spike
(4,541 labels) is the COVID screening wave, which makes the test side the
steadier post-pandemic tail -- arguably the right shape, since "what gets
screened next" is the real question.

THE TRAP THIS MODULE EXISTS TO AVOID. Filtering only the LABELS leaks the
future. If post-cutoff INHIBITS and TARGETS edges stay in the graph, a compound
whose nsp5 activity was published in 2023 gets an M1 path built from 2023
evidence to predict a 2023 label. The training graph must be filtered on the
same date as the labels, and the audits below check that it was.

UNDATED EDGES. 1,583 CHEMICALLY_SIMILAR_TO edges are dated 1970 deliberately --
a computed edge has no assertion date -- so they appear in every training
graph. That is correct but means the split cannot test whether similarity helps
on genuinely future compounds. 2,209 TARGETS edges are undated for a different
reason: ChEMBL had no year for the source document. Those are a real leak of
unknown size, so `undated_policy` makes the choice explicit and both settings
are worth reporting.
"""
from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from kgav.labels import (
    DEFAULT_ACTIVE_BELOW_NM,
    DEFAULT_INACTIVE_ABOVE_NM,
    classify,
)

PLACEHOLDER_DATE = "1970-01-01"
COMPUTED_SOURCES = {"infores:kgav-computed"}


def edge_year(edge: dict) -> int | None:
    """None when the date is the 1970 placeholder rather than a real year."""
    d = edge.get("first_asserted_date") or ""
    if not d or d.startswith("1970"):
        return None
    try:
        return int(d[:4])
    except ValueError:
        return None


@dataclass
class Split:
    cutoff: int
    train_edges: list[dict] = field(default_factory=list)
    # POSITIVES. Named *_labels for continuity with the callers that predate
    # the negative sets.
    test_labels: dict[str, set[str]] = field(default_factory=dict)
    train_labels: dict[str, set[str]] = field(default_factory=dict)
    # MEASURED negatives: compounds assayed against the virus and found
    # inactive, per labels.classify. Previously these were dropped as
    # "label_censored" -- 4,685 measurements discarded -- which left the
    # temporal evaluation with no negatives of its own.
    test_negatives: dict[str, set[str]] = field(default_factory=dict)
    train_negatives: dict[str, set[str]] = field(default_factory=dict)
    stats: Counter = field(default_factory=Counter)
    violations: list[str] = field(default_factory=list)


def build_split(release: Path, cutoff: int, held_out: set[str],
                undated_policy: str = "include",
                inactive_above_nm: float = DEFAULT_INACTIVE_ABOVE_NM,
                active_below_nm: float = DEFAULT_ACTIVE_BELOW_NM) -> Split:
    """Partition a release at `cutoff`.

    undated_policy:
      include   keep undated non-computed edges in the training graph
      exclude   drop them (conservative: an unknown date may be post-cutoff)
      computed_only  keep only undated edges from a computed source

    LABEL POLARITY COMES FROM labels.classify, not from the relation alone.
    This module used to keep every exact-relation label as a positive and
    discard every censored one. Both halves were wrong on the same data: 578
    of 2,826 exact labels sit above 10 uM and are measurements of INACTIVITY,
    while 4,685 censored rows are the best negatives the project has. A
    positive set containing measured-inactive compounds makes every AUC below
    it uninterpretable.
    """
    if undated_policy not in {"include", "exclude", "computed_only"}:
        raise ValueError(f"unknown undated_policy {undated_policy!r}")

    s = Split(cutoff=cutoff)
    # virus -> compounds, kept per side so the ambiguity rule can be applied
    # within each of them before the sets are handed out.
    pos: dict[str, dict[str, set[str]]] = {"train": {}, "test": {}}
    neg: dict[str, dict[str, set[str]]] = {"train": {}, "test": {}}
    for line in (Path(release) / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        pred = e["predicate"]
        year = edge_year(e)
        is_label = pred in held_out

        if is_label:
            # "EC50 > 10000 nM" is a measurement of INACTIVITY: a negative,
            # not a discard. "EC50 = 50000 nM" is one too, despite its exact
            # relation. classify() is the single place that judgement lives.
            verdict = classify(e, inactive_above_nm, active_below_nm)
            if verdict is None:
                s.stats["label_undecidable"] += 1
                continue
            if year is None:
                s.stats["label_undated"] += 1
                continue
            side = "train" if year <= cutoff else "test"
            bucket = pos[side] if verdict == "active" else neg[side]
            bucket.setdefault(e["object"], set()).add(e["subject"])
            s.stats[f"label_{side}_{'active' if verdict == 'active' else 'inactive'}"] += 1
            continue

        if year is None:
            computed = e.get("primary_knowledge_source") in COMPUTED_SOURCES
            if undated_policy == "exclude":
                s.stats["undated_dropped"] += 1
                continue
            if undated_policy == "computed_only" and not computed:
                s.stats["undated_dropped"] += 1
                continue
            s.stats["undated_kept_computed" if computed else "undated_kept"] += 1
            s.train_edges.append(e)
            continue

        if year <= cutoff:
            s.train_edges.append(e)
            s.stats["edge_train"] += 1
        else:
            s.stats["edge_dropped_future"] += 1

    # labels.py's rule: a compound measured both active and inactive against
    # the same virus is a real disagreement -- different cell lines, strains,
    # labs -- and forcing it to one side manufactures a label the data does
    # not support. Applied WITHIN each side, because each is its own label set
    # used for its own question; disagreement ACROSS the cutoff is a different
    # thing and the L5 audit already reports it.
    for side, positives, negatives in (("train", s.train_labels, s.train_negatives),
                                       ("test", s.test_labels, s.test_negatives)):
        # sorted for the same reason as labels.build_labels: this order
        # reaches the audit list and TEST_LABELS.json, so leaving it to set
        # iteration made two runs over one release differ textually.
        for virus in sorted(set(pos[side]) | set(neg[side])):
            a = pos[side].get(virus, set())
            i = neg[side].get(virus, set())
            both = a & i
            if both:
                s.stats[f"{side}_ambiguous_dropped"] += len(both)
            positives[virus] = a - both
            negatives[virus] = i - both
    return s


def audit(split: Split, release: Path, held_out: set[str]) -> list[str]:
    """The leakage checks. Each returns a message when it FAILS."""
    problems: list[str] = []

    # L3 temporal bleed: no dated training edge may postdate the cutoff.
    bled = [e for e in split.train_edges
            if (y := edge_year(e)) is not None and y > split.cutoff]
    if bled:
        problems.append(f"L3 temporal bleed: {len(bled):,} training edges postdate "
                        f"{split.cutoff} (e.g. {bled[0]['predicate']} "
                        f"{bled[0]['first_asserted_date']})")

    # Held-out predicates must not be traversable at all.
    leaked = [e for e in split.train_edges if e["predicate"] in held_out]
    if leaked:
        problems.append(f"label leak: {len(leaked):,} held-out edges are in the "
                        f"training graph")

    # A test-set compound whose only evidence is post-cutoff is unpredictable
    # by construction. Not a leak, but it caps achievable recall, so it must be
    # reported rather than silently depressing every metric.
    train_subjects = {e["subject"] for e in split.train_edges}
    train_subjects |= {e["object"] for e in split.train_edges}
    for virus, drugs in split.test_labels.items():
        unreachable = len(drugs - train_subjects)
        if unreachable:
            problems.append(f"ceiling: {unreachable:,}/{len(drugs):,} test compounds "
                            f"for {virus} have no pre-{split.cutoff} evidence")

    # Test labels must not also be train labels for the same virus: a compound
    # measured both before and after the cutoff is already known. Checked for
    # negatives too -- a compound already known inactive is no more a
    # prospective test case than one already known active.
    for kind, test, train in (("positive", split.test_labels, split.train_labels),
                              ("negative", split.test_negatives, split.train_negatives)):
        for virus, drugs in test.items():
            overlap = drugs & train.get(virus, set())
            if overlap:
                problems.append(f"L5 label overlap: {len(overlap):,} {kind} compounds "
                                f"for {virus} appear in BOTH train and test labels")
    return problems


def write_split(split: Split, release: Path, out: Path) -> None:
    """Write the training graph, keeping only nodes its edges reach."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    keep = {e["subject"] for e in split.train_edges} | \
           {e["object"] for e in split.train_edges}
    nodes = []
    for line in (Path(release) / "nodes.jsonl").read_text().splitlines():
        if line.strip():
            n = json.loads(line)
            if n["id"] in keep:
                nodes.append(n)
    (out / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in nodes))
    (out / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in split.train_edges))
    # Negatives are emitted alongside the positives. A consumer that reads
    # only "test" gets the same positive set as before; one that wants a
    # hard-negative evaluation no longer has to rebuild the labels itself and
    # risk using a different rule from the one the split was built with.
    (out / "TEST_LABELS.json").write_text(json.dumps(
        {"cutoff": split.cutoff,
         "test": {v: sorted(d) for v, d in split.test_labels.items()},
         "train": {v: sorted(d) for v, d in split.train_labels.items()},
         "test_negatives": {v: sorted(d) for v, d in split.test_negatives.items()},
         "train_negatives": {v: sorted(d) for v, d in split.train_negatives.items()}},
        indent=2))
