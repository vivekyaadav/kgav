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
censored value above `inactive_above_nm` counts as evidence of inactivity.

WHERE 10 uM COMES FROM -- the inactive floor is MEASURED, not conventional.

The original justification here was that 10 uM "is the concentration above
which a cell-based antiviral hit is conventionally not pursued". That is an
appeal to convention with no citation behind it, and it was the weakest
sentence in this module. It is not needed, because a censored row states
the answer: "EC50 > 20000" is a screen reporting where it STOPPED -- it
tested to 20 uM and saw nothing. The distribution of those bounds is
therefore a direct measurement of the assay ceilings in this corpus.

censored_bound_census() below computes it. On data/releases/v0.1, over the
18,074 censored HAS_ANTIVIRAL_ACTIVITY_AGAINST rows:

    bound >= 10 uM    17,734   98.1%
    bound <  10 uM       340    1.9%   <- the rows the floor must refuse
    distinct bounds      191
    mode              20 uM    9,924 rows (54.9%)
    next              99.5 uM  2,549 | 100 uM 2,276 | 10 uM 1,009

So 10 uM admits 98.1% of the censored evidence while rejecting exactly the
weak bounds the rule exists to reject. The corpus would support anywhere in
10-20 uM and 10 is the conservative end of that range. tests/test_labels.py
re-derives these figures from the release, so they cannot rot into a claim
with no way back to the data.

THE ACTIVE CUT IS A DIFFERENT QUANTITY, AND ONLY HALF DERIVED.

The two constants are not independent. A censored row "> X" proves
inactivity only if X >= active_below_nm: below that, the compound's true
value could lie between X and the active cut, which would make it ACTIVE.
So `inactive_above_nm >= active_below_nm` is a CONSISTENCY requirement, not
a convention -- violate it and the two classes overlap, as

    classify({">": 2000}, inactive_above=1000, active_below=10000)
      -> "inactive", while "= 5000" -> "active"

shows: one compound can be both. Equality is the tightest consistent
choice, because it admits the most censored evidence without overlap, and
that is a real argument for the two defaults matching. test_labels.py
pins the inequality.

What has NO derivation is the VALUE of active_below_nm. The commit that
introduced both (88074a6, 2026-09-12) reasons only about the censored
floor. "Above what bound can an inactive call be believed" is a fact about
assay ceilings and is measured above; "what counts as a hit" is a decision
about the task, and no distribution in this corpus answers it.

It is also only half of the convention it invokes. Published antiviral
screening criteria pair the potency cut with a SELECTIVITY requirement --
EC50 < 10 uM AND SI >= 5 to 10, SI being CC50/EC50 -- and this module
applies the potency clause alone. Measured on v0.1, of 5,465 positive
(compound, virus) pairs: 2,167 (39.7%) pass SI >= 5, 392 of the 2,559 that
have a selectivity index FAIL it, 65 of those have SI < 1 (cytotoxic BELOW
the effective concentration), and 2,906 (53.2%) have no cytotoxicity data
at all. scripts/hard_negatives.py --selectivity defaults to "all", so the
published numbers apply no cytotoxicity control.

Two consequences are left as declared open choices rather than changed
here, because they move every published figure:

  - the two cuts are TIED by default, so there is no undecided band: 9,999
    nM is active and 10,001 nM is inactive. classify() takes them as
    independent parameters and will accept a gap; only the CLI ties them,
    and --active-below-nm does not exist.
  - 56.5% of positives at 10 uM are not positive at 1 uM, so the majority
    of the positive set is weak-potency chemical matter grouped with
    sub-nanomolar compounds under one label.

AMBIGUOUS COMPOUNDS ARE EXCLUDED, NOT ASSIGNED. A compound with both an active
and an inactive measurement against the same virus is real -- different cell
lines, different strains, different labs -- and forcing it to one side would
manufacture a label the data does not support. It is dropped and counted.
"""
from __future__ import annotations

import json
import math
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

# provenance imports nothing from this package, so this is a plain
# top-level import and not a cycle dodged with a function-local one.
from kgav.provenance import check_file

ACTIVITY_PREDICATE = "HAS_ANTIVIRAL_ACTIVITY_AGAINST"
# Binding constants belong here. Reading only EC50/IC50 meant classify()
# returned None for all 667 Ki/Kd measurements in the chembl layer -- not
# because anyone judged them undecidable, but because the function never
# looked. 504 of them are sub-10uM binders, the strongest direct-acting
# evidence in that layer, and routing edges by this verdict would have
# deleted every one.
#
# LABEL-NEUTRAL, measured before changing it: no HAS_ANTIVIRAL_ACTIVITY_AGAINST
# edge in either layer carries ki_nm or kd_nm, because chembl.py drops binding
# constants on ORGANISM targets as mis-entered. So this widens what counts as
# "acts on" at the PROTEIN level and moves no organism-level label.
POTENCY_FIELDS = ("ec50_nm", "ic50_nm", "ki_nm", "kd_nm")
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


def censored_bound_census(release: Path,
                          predicate: str = ACTIVITY_PREDICATE) -> Counter:
    """{bound in nM: rows} over the censored ('>' / '>=') activity rows.

    WHAT THIS IS FOR. `inactive_above_nm` is the one threshold in this
    project that does NOT need a convention, because the data states it. A
    row reading "EC50 > 20000" is a screen reporting where it stopped: it
    tested to 20 uM and saw nothing. So the distribution of those bounds is
    a measurement of the assay ceilings actually present in the corpus, and
    the floor can be read off it instead of asserted.

    Scoped to the predicate build_labels reads, because that is the
    population the threshold decides. A census over a wider set would be a
    true statement about something else.

    Kept here rather than in a script so the docstring's figures can be
    re-derived by the test that checks them; a number in a comment with no
    way back to the data is the thing this repository keeps finding.
    """
    out: Counter = Counter()
    for line in (Path(release) / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        if e["predicate"] != predicate:
            continue
        q = e.get("qualifiers") or {}
        if (q.get("relation") or "=") not in (">", ">="):
            continue
        v = _potency(e)
        if v is not None:
            out[v] += 1
    return out


def edge_sources(edge: dict) -> set[str]:
    """Knowledge sources behind one edge.

    merge_edges joins them with "|" when two layers witness one fact, so a
    merged edge reads "infores:string|infores:virhostnet". Membership has to
    be tested per component or a merged edge matches no source at all.
    """
    raw = edge.get("primary_knowledge_source") or ""
    return {s for s in raw.split("|") if s}


def build_labels(release: Path, year_from: int | None = None,
                 year_to: int | None = None,
                 inactive_above_nm: float = DEFAULT_INACTIVE_ABOVE_NM,
                 active_below_nm: float = DEFAULT_ACTIVE_BELOW_NM,
                 sources: set[str] | None = None) -> Labels:
    """Label sets from the activity edges, optionally restricted by year.

    year_from/year_to are inclusive and let the same function produce a
    temporal test set (2022-2025) or a cross-sectional one (no bounds).

    `sources` restricts the labels to measurements from given knowledge
    sources, which is how a SINGLE-SCREEN evaluation is built.

    Why that matters more than it sounds. Pooling ChEMBL and panel labels
    reintroduces the confound the panel was ingested to remove: ChEMBL actives
    are antiviral research compounds carrying viral-target annotations, panel
    inactives are library compounds carrying host-target annotations, and a
    scorer can separate the two classes by recognising which collection a
    compound came from. Restricted to one screen, both classes come off the
    same plates under one protocol over one concentration range, so the
    populations are matched BY CONSTRUCTION rather than by a filter --
    which is what --selectivity-applies both was approximating.
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
        if sources is not None and not (edge_sources(e) & sources):
            out.stats["other_source"] += 1
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

    # SORTED, because a set of strings iterates in a different order in every
    # process (str hashing is randomised per interpreter). That order reaches
    # the results files through the label dicts, so two runs over the SAME
    # release produced files that differed textually while every number in
    # them was identical -- which is precisely the signal the provenance
    # stamps exist to make trustworthy.
    for virus in sorted(set(active) | set(inactive)):
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


def auc_interval(auc: float, n_pos: int, n_neg: int,
                 z: float = 1.96) -> tuple[float, float]:
    """95% interval for an AUC, Hanley & McNeil (1982).

    Assumes exponential score distributions, which is conservative for the
    small-n cases that matter here. Clamped to [0, 1] because the normal
    approximation runs outside the range when n is tiny, and an AUC of 1.04
    reported to a reader is worse than a wide interval.
    """
    if n_pos < 1 or n_neg < 1:
        return (0.0, 1.0)
    q1 = auc / (2.0 - auc)
    q2 = 2.0 * auc * auc / (1.0 + auc)
    var = (auc * (1 - auc) + (n_pos - 1) * (q1 - auc * auc)
           + (n_neg - 1) * (q2 - auc * auc)) / (n_pos * n_neg)
    se = math.sqrt(var) if var > 0 else 0.0
    return (max(0.0, auc - z * se), min(1.0, auc + z * se))


def auc_ceiling(reached_pos: int, reached_neg: int,
                n_pos: int, n_neg: int) -> float:
    """The highest AUC a scorer with this coverage could attain, ranking
    everything it reaches PERFECTLY.

    WHY THIS IS NOT A DETAIL. Unreached compounds take the imputed floor,
    so a scorer that reaches few positives cannot separate the classes
    however good its ranking is -- the floor does the arithmetic. On
    release v0.1, SARS-CoV-2, every channel sits within 0.04 of this
    ceiling and three sit ON it:

        M1  reaches  3 of 89   measured 0.5135   ceiling 0.5136
        M8  reaches  1 of 89   measured 0.506    ceiling 0.506
        M7  reaches  0 of 89   measured 0.4973   ceiling 0.4973

    M1 achieved the best score it was capable of, and that score is
    0.514. M7's ceiling is BELOW chance: it reaches 39 negatives and no
    positives, so those negatives necessarily outrank the floor-tied
    positives. "Slightly below chance" is arithmetic there, not a result.

    AND THE EXISTING FLOOR DOES NOT CATCH THIS. MIN_REACHED_POSITIVES
    gates on the COUNT of positives reached, and the count does not
    determine the ceiling -- both coverages do. M2 reaches 11, clears
    the floor of 10, and has a ceiling of 0.5125: it was admitted to the
    evaluable set while being incapable of demonstrating separation.

    THE BEST CASE, spelled out, because the formula is only obvious once:
    every reached positive outranks every negative; the unreached
    positives sit at the floor, tied with the unreached negatives (half
    credit each) and below the reached negatives (no credit).

    Reported, not enforced. Whether this replaces the count as the
    readability gate is a decision about what the evaluation claims, and
    it belongs in the prereg rather than in a helper.
    """
    if n_pos <= 0 or n_neg <= 0:
        return 0.5
    k = max(0, min(reached_pos, n_pos))
    j = max(0, min(reached_neg, n_neg))
    return (k * n_neg + 0.5 * (n_pos - k) * (n_neg - j)) / (n_pos * n_neg)


def measured_reach(results: Path, release: Path) -> dict[str, int] | None:
    """Held-out positives each virus's metapaths actually reach.

    `{virus label: COMBINED.pos}`, or None when that cannot be known.

    WHY THIS IS NOT DERIVED FROM THE GRAPH. Whether a metapath reaches a
    labelled compound is the result of walking it against a label set,
    which is what compare_metapaths.py does. Edge counts cannot answer
    it: DENV-1 has 349 labels and 53 virus-host interactions and reaches
    none of the 349.

    ONLY IF THE FILE IS CURRENT. provenance.check_file is the same audit
    `make check-results` runs, so a comparison measured on a graph that
    has since moved is reported as unmeasured rather than quoted here.
    Returning None for "cannot be known" rather than 0 for "reaches
    nothing" is the whole point: those are opposite claims, and this
    project has twice shipped a figure from a previous graph.

    LIVED IN scripts/audit_build.py UNTIL 2026-11-02, which is why the
    agent could not express the tier this measures. The audit report
    prints four per-virus tiers and says "THE AGENT ENFORCES THIS"; the
    agent could see three of them, because the one fact separating
    ROUTES REACH NO LABEL from `evaluable` was computed in a script it
    does not import. Five viruses -- DENV-1, DENV-3, HCoV-NL63, JEV,
    TBEV -- were therefore ranked with a caveat saying their reliability
    was "unknown", when it is known and is zero. One owner now, read by
    the audit and by guards.VirusRoutes alike.
    """
    # A DIRECTORY, AND THIS REFUSES ANYTHING ELSE. None means "this cannot
    # be known", and every branch below returns it for a reason a reader of
    # the output can act on: the file is absent, stale, or unreadable. A
    # caller that passed a FILE here -- scripts/kgav_ask.py did, for weeks,
    # because --results names a results file -- got the same None, so the
    # miswiring was indistinguishable from an honest absence and the agent
    # reported every virus's reach as unmeasured. That is a bug in the
    # caller, not a fact about the data, and it is now loud.
    results = Path(results)
    if results.exists() and not results.is_dir():
        raise ValueError(
            f"measured_reach wants the results DIRECTORY, not a file in it: "
            f"got {results}. Pass {results.parent} instead.")
    path = results / "metapath_comparison.json"
    if not path.exists():
        return None
    status, _msg = check_file(path, release)
    if status != "ok":
        return None
    try:
        per = json.loads(path.read_text())["per_virus"]
    except (OSError, ValueError, KeyError):
        return None
    out = {}
    for virus, channels in per.items():
        pos = (channels.get("COMBINED") or {}).get("pos")
        if isinstance(pos, int):
            out[virus] = pos
    return out or None


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
    lo, hi = auc_interval(auc, n_p, n_n)
    reached_p = sum(1 for d in pos if d in scores)
    reached_n = sum(1 for d in neg if d in scores)
    # AGAINST THE COUNTS ACTUALLY RANKED, not the class sizes. Under
    # imputation these are the same (every compound is ranked, reached or
    # floored). Without it, unreached compounds leave the evaluation
    # entirely, so nothing is holding the scorer down and a perfect
    # ranking reaches 1.0 -- passing len(pos)/len(neg) there would
    # impose an imputation ceiling on a cell with no imputation. Found
    # by a surviving mutant, which turned out to be the correct code.
    ceiling = auc_ceiling(reached_p, reached_n, n_p, n_n)
    return {"n_pos": n_p, "n_neg": n_n, "auc": auc, "evaluable": True,
            # THE INTERVAL TRAVELS WITH THE NUMBER. A point AUC says nothing
            # about whether the scorer discriminates: 0.523 on 89 positives
            # spans [0.462, 0.584] and is indistinguishable from chance, while
            # 0.523 on 10,000 would not be. Everything that later decides
            # whether to trust a channel reads these two fields, so they
            # cannot be left for a reader to compute.
            "auc_ci_lo": lo, "auc_ci_hi": hi,
            "imputed": bool(impute_missing),
            # Coverage is now reported separately from n_pos/n_neg: it is what
            # fraction the scorer actually REACHED, and a high AUC on low
            # coverage means the floor is doing the work, not the paths.
            "coverage_pos": reached_p / max(len(pos), 1),
            "coverage_neg": reached_n / max(len(neg), 1),
            # WHAT THIS CELL COULD HAVE SCORED AT BEST. An AUC near 0.5
            # means one of two completely different things -- the scorer
            # ranked badly, or it could not have ranked well -- and
            # nothing in this dict distinguished them. See auc_ceiling.
            "auc_ceiling": ceiling,
            "auc_headroom": ceiling - 0.5}
