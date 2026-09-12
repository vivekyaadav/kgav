"""Label sets with MEASURED negatives, not untested compounds.

WHAT WAS WRONG BEFORE. Every evaluation so far treated the positives as
compounds measured active and everything else in the pool as negative. But most
of that pool was never tested against the virus. So "lift 1.40" meant "ranks
actives above UNTESTED compounds" -- which a graph could achieve merely by
separating well-studied chemical matter from obscure chemical matter, with no
antiviral insight at all.

ChEMBL's censored values fix this. A row reading "EC50 > 10000 nM" is a
compound that WAS assayed and found inactive. There are 11,260 such records
across the coronavirus targets, and the temporal split alone discarded 4,748 of
them. They are the best negatives available anywhere in this project -- better
than mining failed clinical trials, because they are measurements rather than
inferences about why a programme stopped.

THE POTENCY THRESHOLD MATTERS. "EC50 > 100 nM" is not an inactive compound; it
is a weak lower bound that may sit below the assay's top concentration. Only a
censored value above `inactive_above_nm` counts as evidence of inactivity, and
the default (10,000 nM = 10 uM) is the concentration above which a cell-based
antiviral hit is conventionally not pursued.

AMBIGUOUS COMPOUNDS ARE EXCLUDED, NOT ASSIGNED. A compound with both an active
and an inactive measurement against the same virus is real -- different cell
lines, different strains, different labs -- and forcing it to one side would
manufacture a label the data does not support. It is dropped and counted.
"""
from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

ACTIVITY_PREDICATE = "HAS_ANTIVIRAL_ACTIVITY_AGAINST"
POTENCY_FIELDS = ("ec50_nm", "ic50_nm")
DEFAULT_INACTIVE_ABOVE_NM = 10_000.0
DEFAULT_ACTIVE_BELOW_NM = 10_000.0


@dataclass
class Labels:
    """Per-virus label sets. `negatives` are MEASURED inactive."""
    positives: dict[str, set[str]] = field(default_factory=dict)
    negatives: dict[str, set[str]] = field(default_factory=dict)
    stats: Counter = field(default_factory=Counter)

    def evaluable(self, virus: str) -> set[str]:
        """The only compounds a hard-negative evaluation may rank."""
        return self.positives.get(virus, set()) | self.negatives.get(virus, set())


def _potency(edge: dict) -> float | None:
    q = edge.get("qualifiers") or {}
    for f in POTENCY_FIELDS:
        v = q.get(f)
        if v is not None:
            return float(v)
    return None


def classify(edge: dict, inactive_above_nm: float = DEFAULT_INACTIVE_ABOVE_NM,
             active_below_nm: float = DEFAULT_ACTIVE_BELOW_NM) -> str | None:
    """'active', 'inactive', or None when the measurement decides nothing."""
    q = edge.get("qualifiers") or {}
    relation = q.get("relation") or "="
    value = _potency(edge)
    if value is None:
        return None

    if relation in ("=", "~"):
        return "active" if value <= active_below_nm else "inactive"
    if relation in (">", ">="):
        # A lower bound only proves inactivity if the bound itself is already
        # past the point where a hit would be pursued.
        return "inactive" if value >= inactive_above_nm else None
    if relation in ("<", "<="):
        return "active" if value <= active_below_nm else None
    return None


def build_labels(release: Path, year_from: int | None = None,
                 year_to: int | None = None,
                 inactive_above_nm: float = DEFAULT_INACTIVE_ABOVE_NM,
                 active_below_nm: float = DEFAULT_ACTIVE_BELOW_NM) -> Labels:
    """Label sets from the activity edges, optionally restricted by year.

    year_from/year_to are inclusive and let the same function produce a
    temporal test set (2022-2025) or a cross-sectional one (no bounds).
    """
    out = Labels()
    active: dict[str, set[str]] = {}
    inactive: dict[str, set[str]] = {}

    for line in (Path(release) / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        if e["predicate"] != ACTIVITY_PREDICATE:
            continue
        date = e.get("first_asserted_date") or ""
        year = int(date[:4]) if date[:4].isdigit() and not date.startswith("1970") else None
        if year_from is not None and (year is None or year < year_from):
            out.stats["out_of_window"] += 1
            continue
        if year_to is not None and (year is None or year > year_to):
            out.stats["out_of_window"] += 1
            continue

        verdict = classify(e, inactive_above_nm, active_below_nm)
        if verdict is None:
            out.stats["undecidable"] += 1
            continue
        bucket = active if verdict == "active" else inactive
        bucket.setdefault(e["object"], set()).add(e["subject"])
        out.stats[verdict] += 1

    for virus in set(active) | set(inactive):
        a = active.get(virus, set())
        i = inactive.get(virus, set())
        both = a & i
        if both:
            # Real disagreement: different cell lines, strains, labs. Forcing a
            # side would manufacture a label the data does not support.
            out.stats["ambiguous_dropped"] += len(both)
        out.positives[virus] = a - both
        out.negatives[virus] = i - both
    return out


def evaluate_against_negatives(scores: dict[str, float], pos: set[str],
                               neg: set[str], impute_missing: bool = True) -> dict:
    """Rank only measured compounds and report AUC plus enrichment.

    ROC AUC is the right summary here, not Hits@k: with a measured negative set
    the question is separation between two known classes, and AUC is
    prevalence-independent so it stays comparable across viruses whose base
    rates differ by two orders of magnitude. 0.5 is chance.
    """
    # IMPUTATION IS NOT OPTIONAL IN PRACTICE. A metapath scorer only reaches
    # compounds that have paths, and measured-INACTIVE compounds frequently
    # have none -- in the fixture M1's pool contained actives only. Evaluating
    # each scorer on just the compounds it happens to reach lets it choose its
    # own test set, which is how a scorer with 3% coverage posts a perfect AUC.
    # A compound with no path has no evidence, which IS a prediction of "not
    # active", so it takes the floor score.
    floor = min(scores.values(), default=0.0) - 1.0 if impute_missing else None
    def _get(d: str) -> float | None:
        if d in scores:
            return scores[d]
        return floor
    ranked_pos = [v for v in (_get(d) for d in pos) if v is not None]
    ranked_neg = [v for v in (_get(d) for d in neg) if v is not None]
    n_p, n_n = len(ranked_pos), len(ranked_neg)
    if not n_p or not n_n:
        return {"n_pos": n_p, "n_neg": n_n, "auc": None, "evaluable": False}

    # Rank-sum AUC with ties counted as half, which matters because DWPC
    # produces many exact ties at percentile boundaries.
    wins = 0.0
    neg_sorted = sorted(ranked_neg)
    import bisect
    for p in ranked_pos:
        lower = bisect.bisect_left(neg_sorted, p)
        equal = bisect.bisect_right(neg_sorted, p) - lower
        wins += lower + 0.5 * equal
    auc = wins / (n_p * n_n)
    return {"n_pos": n_p, "n_neg": n_n, "auc": auc, "evaluable": True,
            "imputed": bool(impute_missing),
            # Coverage is now reported separately from n_pos/n_neg: it is what
            # fraction the scorer actually REACHED, and a high AUC on low
            # coverage means the floor is doing the work, not the paths.
            "coverage_pos": sum(1 for d in pos if d in scores) / max(len(pos), 1),
            "coverage_neg": sum(1 for d in neg if d in scores) / max(len(neg), 1)}
