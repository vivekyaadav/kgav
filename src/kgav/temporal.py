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
    test_labels: dict[str, set[str]] = field(default_factory=dict)
    train_labels: dict[str, set[str]] = field(default_factory=dict)
    stats: Counter = field(default_factory=Counter)
    violations: list[str] = field(default_factory=list)


def build_split(release: Path, cutoff: int, held_out: set[str],
                undated_policy: str = "include") -> Split:
    """Partition a release at `cutoff`.

    undated_policy:
      include   keep undated non-computed edges in the training graph
      exclude   drop them (conservative: an unknown date may be post-cutoff)
      computed_only  keep only undated edges from a computed source
    """
    if undated_policy not in {"include", "exclude", "computed_only"}:
        raise ValueError(f"unknown undated_policy {undated_policy!r}")

    s = Split(cutoff=cutoff)
    for line in (Path(release) / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        pred = e["predicate"]
        year = edge_year(e)
        is_label = pred in held_out

        if is_label:
            # Only exact relations are usable labels: "EC50 > 10000 nM" is a
            # measurement of INACTIVITY and would invert the target.
            if (e.get("qualifiers") or {}).get("relation") not in (None, "="):
                s.stats["label_censored"] += 1
                continue
            if year is None:
                s.stats["label_undated"] += 1
                continue
            bucket = s.train_labels if year <= cutoff else s.test_labels
            bucket.setdefault(e["object"], set()).add(e["subject"])
            s.stats["label_train" if year <= cutoff else "label_test"] += 1
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
    # measured both before and after the cutoff is already known.
    for virus, drugs in split.test_labels.items():
        overlap = drugs & split.train_labels.get(virus, set())
        if overlap:
            problems.append(f"L5 label overlap: {len(overlap):,} compounds for "
                            f"{virus} appear in BOTH train and test labels")
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
    (out / "TEST_LABELS.json").write_text(json.dumps(
        {"cutoff": split.cutoff,
         "test": {v: sorted(d) for v, d in split.test_labels.items()},
         "train": {v: sorted(d) for v, d in split.train_labels.items()}}, indent=2))
