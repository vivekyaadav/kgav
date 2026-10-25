"""What the pathway route is measured to be worth.

THE CLAIM THIS EXISTS TO STOP. An agent that orders compounds and mentions a
route through a host pathway has implied, without saying it, that the route
bears on the order. On this release it does not. pathway_landscape.py measures
the overlap between the drugged human proteome and each virus's CRISPR
dependency factors twice -- naive, and against a background sampled to match
the observed degree distribution -- and the matched test is the one that
counts, because well-studied proteins are both likelier to be drugged and
likelier to be called dependency factors.

WHY A FILE AND NOT A SENTENCE. The result used to exist only as terminal
output from a script that wrote nothing, so the agent could not cite it and
nothing re-ran it. A sentence in a docstring would have been a figure with no
record of what produced it -- the exact objection calibration.py raises about
hardcoded weights. This loads a results file carrying its own release,
permutation count and seed, so a claim about the pathway route can always be
traced to the run that measured it, and a run on a different release can be
refused rather than quoted.

THE TEST IS ONE-SIDED, and saying so is the whole of what this module gets
right that the first version got wrong. empirical_p measures the ENRICHMENT
tail only: the probability that a degree-matched random protein set hits as
many dependency factors as the drugged proteome did. A ratio below 1 is
therefore a description, not a result. Nothing here has tested depletion, and
a summary that prints "ratios 0.00-0.94" beside "no enrichment" invites a
reader to take the low end as a finding when no test was run in that
direction.

That matters on the measured release rather than hypothetically. Of the eight
taxa with dependency factors, six cluster at 0.86-0.94 and two sit far below:
one with 18 dependency factors against a matched expectation of 3.1 and no
overlap at all, and one with 50 against 7.5 and two. Read as Poisson tails
those are nominally 0.045 and 0.020, and neither survives eight taxa' worth of
multiple testing -- but this module cannot say even that, because the run it
loads never computed a lower tail. It says the direction is untested and names
the two taxa whose ratio is least readable, which is what it knows.

WHAT MIN_EXPECTED_HITS GATES, after that correction. It was introduced as a
power threshold for depletion, which was answering a question the module never
asks -- an internal inconsistency that the real numbers exposed: it classified
all eight taxa as "testable" and so reported a 0.00 ratio with no qualification
at all. It now gates whether the RATIO can be read as anything, not whether the
taxon was tested. Enrichment is tested for every taxon with a dependency
factor, however small, because a small taxon can still return a significant
excess; it is the low end of the ratio that goes uninterpretable first.

Enrichment is reported wherever it is found, whatever the expectation, because
a significant result is significant and an agent built on "nothing is enriched"
must say so the moment that stops being true. That mirrors the discriminates
branch in guards._measured_at_chance, for the same reason.

FORWARD COMPATIBLE. If a results file carries `p_matched_less`, the depletion
direction was measured and this module reports it. Until the writer produces
one, depletion is stated as untested rather than inferred from the ratio.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

# Below this matched expectation the RATIO carries no information in the low
# direction, so it is reported but not read. Under a Poisson background with
# mean m the most extreme downward observation is zero hits, with probability
# exp(-m); exp(-3) = 0.0498 and exp(-2.9) = 0.055, so below m = 3 even a total
# absence of overlap is an ordinary draw. Derived, not chosen: a threshold
# picked for looking round would be a fitted parameter with no record of what
# fitted it.
#
# This does NOT gate whether a taxon was tested. The run measures the
# enrichment tail, and a taxon with a small expectation can still return a
# significant excess -- NCBITaxon:11082 would need 8 hits against 3.1 expected.
MIN_EXPECTED_HITS = 3.0

# An empirical p from N permutations cannot resolve below 1/(N+1). At 200 the
# finest p is 0.005, which is enough for a 0.05 test but leaves a "p = 0.0000"
# unable to distinguish 0.004 from 0.0000001. Runs coarser than this are
# refused rather than quoted: the agent states these figures to a user.
MIN_PERMUTATIONS = 200

ALPHA = 0.05


@dataclass(frozen=True)
class TaxonOverlap:
    """One taxon's measured overlap, and whether it could be measured."""
    virus: str
    n_dependency: int
    n_observed: int
    expected: float
    p_naive: float
    p_matched: float
    # The DEPLETION tail, when the run computed one. None means the direction
    # was not measured, which is not the same as measured and null, and the
    # distinction is the whole correction this class carries.
    p_matched_less: float | None = None
    # The threshold this taxon's verdicts are judged against, set by the
    # loader from how many tests the run performed in total. Carried per-row
    # rather than read from the module constant so a verdict can never be
    # computed against a different threshold from the one reported next to
    # it.
    alpha: float = ALPHA
    n_tests: int = 1

    @property
    def ratio(self) -> float:
        """Observed over matched expectation, derived rather than stored.

        A stored ratio is a third number that can disagree with the two it
        comes from, and the disagreement would be invisible. It is also not
        always finite: a taxon whose degree-matched background is expected to
        hit nothing has an infinite ratio if anything at all is observed, and
        JSON has no way to write that. Deriving it here keeps one source of
        truth and keeps the results file plain.
        """
        if self.expected > 0:
            return self.n_observed / self.expected
        return float("inf") if self.n_observed else 0.0

    @property
    def ratio_readable(self) -> bool:
        """Whether the LOW end of this taxon's ratio carries information.

        Not whether the taxon was tested -- it was, in the enrichment
        direction, which is the direction the run measures. Below the
        threshold a ratio of 0.00 is an ordinary draw from the background and
        says nothing, so it is printed and not interpreted.
        """
        return self.expected >= MIN_EXPECTED_HITS

    @property
    def enriched(self) -> bool:
        """Clears the degree-matched null in the enrichment direction.

        Against the CORRECTED threshold. See Landscape.alpha: a run over
        eight taxa in two directions performs sixteen tests, and an
        uncorrected 0.05 would manufacture roughly one finding from noise
        every other run -- in a caveat whose entire job is to stop the agent
        claiming more than was measured.
        """
        return self.p_matched < self.alpha and self.ratio > 1.0

    @property
    def depletion_tested(self) -> bool:
        """Whether the run computed a lower tail at all."""
        return self.p_matched_less is not None

    @property
    def depleted(self) -> bool:
        """Clears the null downward. False when the direction was not tested.

        False here means "not shown to be depleted", never "shown not to be".
        `depletion_tested` is what separates the two, and every sentence this
        module emits about a sub-1 ratio consults it.
        """
        return (self.depletion_tested
                and self.p_matched_less < self.alpha and self.ratio < 1.0)

    @property
    def at_null(self) -> bool:
        """Tested for enrichment and not enriched.

        Every taxon with a dependency factor qualifies: the run measures the
        enrichment tail for all of them. This used to require
        `ratio_readable`, which silently dropped the two smallest taxa out of
        the count of what had been tested, while still printing their ratios
        in the range.
        """
        return not self.enriched

    def describe(self) -> str:
        if self.enriched:
            return (f"{self.virus}: {self.n_observed:,} of "
                    f"{self.n_dependency:,} dependency factors are drugged "
                    f"against {self.expected:.1f} expected under a "
                    f"degree-matched background — ratio {self.ratio:.2f}, "
                    f"matched p {self.p_matched:.4f}, which CLEARS the null")
        if self.depleted:
            return (f"{self.virus}: {self.n_observed:,} of "
                    f"{self.n_dependency:,} dependency factors are drugged "
                    f"against {self.expected:.1f} expected — ratio "
                    f"{self.ratio:.2f}, lower-tail matched p "
                    f"{self.p_matched_less:.4f}, DEPLETED")
        head = (f"{self.virus}: {self.n_observed:,} of {self.n_dependency:,} "
                f"dependency factors are drugged against {self.expected:.1f} "
                f"expected under a degree-matched background — ratio "
                f"{self.ratio:.2f}, matched p {self.p_matched:.2f}, not "
                f"enriched")
        if self.ratio < 1.0 and not self.depletion_tested:
            head += ("; the run measured the enrichment tail only, so the "
                     "sub-1 ratio is not a measured depletion")
        elif self.ratio < 1.0:
            # Tested downward and found nothing is a result, not a reason to
            # go quiet. Saying only "not enriched" beside a ratio of 0.27
            # leaves the reader to supply the other half themselves.
            head += (f"; tested downward too and not depleted, lower-tail p "
                     f"{self.p_matched_less:.4f} against a threshold of "
                     f"{self.alpha:.4f}")
        if not self.ratio_readable:
            head += (f"; and with only {self.expected:.1f} expected, below "
                     f"{MIN_EXPECTED_HITS:.0f}, even a zero overlap would be "
                     f"an ordinary draw, so the ratio is not readable "
                     f"downward")
        return head


@dataclass
class Landscape:
    """Every taxon's measured pathway overlap, with the run that measured it."""
    release: str = ""
    permutations: int = 0
    seed: int = 0
    measured_on: str = ""
    taxa: dict[str, TaxonOverlap] = field(default_factory=dict)
    uncovered: list[str] = field(default_factory=list)
    source_file: str = ""
    refused_reason: str = ""

    @property
    def usable(self) -> bool:
        return not self.refused_reason

    @property
    def protocol(self) -> str:
        return (f"a degree-matched permutation null, "
                f"{self.permutations:,} permutations, seed {self.seed}, "
                f"significance at {self.alpha:.4f} "
                f"({ALPHA} Bonferroni-corrected for {self.n_tests} tests)")

    @property
    def n_tests(self) -> int:
        """Tests performed: one per taxon per measured direction."""
        rows = self._sorted()
        return max(1, len(rows) * (2 if self.depletion_tested() else 1))

    @property
    def alpha(self) -> float:
        """The corrected threshold every verdict here is judged against."""
        return ALPHA / self.n_tests

    def for_virus(self, virus: str) -> TaxonOverlap | None:
        """This taxon's overlap, or None. Never another taxon's."""
        if not self.usable:
            return None
        return self.taxa.get(virus)

    def tested(self) -> list[TaxonOverlap]:
        """Every taxon the run measured: all of them, in the enrichment tail."""
        return self._sorted()

    def unreadable_ratio(self) -> list[TaxonOverlap]:
        """Taxa whose ratio says nothing downward, named so a reader can
        discount the low end of the range rather than read a finding into it."""
        return [t for t in self._sorted() if not t.ratio_readable]

    def enriched(self) -> list[TaxonOverlap]:
        return [t for t in self._sorted() if t.enriched]

    def depleted(self) -> list[TaxonOverlap]:
        return [t for t in self._sorted() if t.depleted]

    def depletion_tested(self) -> bool:
        """Whether the run measured the lower tail for every taxon it covers.

        All-or-nothing on purpose. A partially measured direction cannot
        support "none is depleted", and the honest report of a mixed file is
        that the direction is untested.
        """
        rows = self._sorted()
        return bool(rows) and all(t.depletion_tested for t in rows)

    def at_null(self) -> list[TaxonOverlap]:
        return [t for t in self._sorted() if t.at_null]

    def _sorted(self) -> list[TaxonOverlap]:
        if not self.usable:
            return []
        return [self.taxa[k] for k in sorted(self.taxa)]

    def ratio_range(self) -> tuple[float, float] | None:
        """Range over every tested taxon, or None if there are none.

        Over ALL of them, not the subset whose ratio is readable downward.
        Restricting the range silently dropped taxa that the run did test,
        and the summary names the unreadable ones instead -- a reader is
        owed the real spread plus the reason its low end is soft, not a
        tidier number.
        """
        rs = [t.ratio for t in self._sorted()]
        return (min(rs), max(rs)) if rs else None

    def summary(self) -> str:
        """What the pathway route is worth on this release, in one sentence.

        Every branch is reachable, including the one that says the route works.
        A module whose only output is "no enrichment" would keep saying it
        after the finding changed.
        """
        if not self.usable:
            return (f"the pathway route has not been measured on this "
                    f"release ({self.refused_reason})")
        if not self.taxa:
            return ("no taxon in this release has a dependency factor, so "
                    "there is no overlap to measure and the pathway route is "
                    "untested rather than null")
        tested = self.tested()
        enriched = self.enriched()

        # WHY THE LOW END OF THE RANGE IS QUALIFIED RATHER THAN TRIMMED. On
        # the measured release the range runs 0.00-0.94 and the two lowest
        # values are the two smallest taxa. Dropping them would report a
        # tidier 0.86-0.94 and hide that two taxa were measured; printing
        # 0.00 bare invites it to be read as depletion, which nothing tested.
        # So both are said.
        soft = self.unreadable_ratio()
        tail = ""
        if soft:
            tail += (f" The low end of that range is not a finding: "
                     f"{len(soft)} taxon(s) have a matched expectation under "
                     f"{MIN_EXPECTED_HITS:.0f}, where even a zero overlap is "
                     f"an ordinary draw "
                     + ", ".join(f"{t.virus} ({t.expected:.1f} expected)"
                                 for t in soft) + ".")
        if not self.depletion_tested():
            tail += (" Depletion was not tested at all: the run measures the "
                     "enrichment tail only, so no ratio below 1 here is a "
                     "measured result in that direction.")
        elif self.depleted():
            tail += (" In the other direction " + ", ".join(
                f"{t.virus} (ratio {t.ratio:.2f}, lower-tail p "
                f"{t.p_matched_less:.4f})" for t in self.depleted())
                + " clears the null downward: fewer dependency factors are "
                  "drugged than chance predicts.")
        else:
            # A MEASURED NULL IN THE OTHER DIRECTION IS A RESULT, and this
            # branch is the whole return on measuring the lower tail. Without
            # it the summary fell silent about depletion the moment it became
            # measurable: the "not tested" clause switched off and the
            # "clears downward" clause never fired, so the reader learned
            # less from a run that had done more work.
            down = min(self._sorted(), key=lambda t: t.p_matched_less)
            tail += (f" Nor depleted in any: the strongest lower tail is "
                     f"{down.virus} at p {down.p_matched_less:.4f}, above "
                     f"the corrected threshold. The ratios sit below 1 "
                     f"throughout and none of them is a measured depletion.")

        if enriched:
            names = ", ".join(f"{t.virus} (ratio {t.ratio:.2f}, p "
                              f"{t.p_matched:.4f})" for t in enriched)
            return (f"measured under {self.protocol}, the drugged proteome IS "
                    f"enriched for dependency factors in {len(enriched)} of "
                    f"{len(tested)} taxa: {names}.{tail}")
        lo, hi = self.ratio_range()
        # The SMALLEST matched p: the taxon that comes closest to clearing the
        # null, which is the only one whose figure supports the "none below
        # ALPHA" beside it. Quoting the largest would put a number next to a
        # claim it cannot carry. "the strongest matched p" was the original
        # wording and reads either way, so the sentence names the comparison.
        closest = min(tested, key=lambda t: t.p_matched)
        return (f"measured under {self.protocol}, the drugged proteome is not "
                f"enriched for dependency factors in any of the "
                f"{len(tested)} taxa with dependency factors — ratios "
                f"{lo:.2f}–{hi:.2f}, and the closest any taxon comes to "
                f"clearing the null is {closest.virus} at matched p "
                f"{closest.p_matched:.2f}, none below the corrected "
                f"{self.alpha:.4f}.{tail}")

    def report(self) -> list[str]:
        lines = [f"pathway landscape: {self.release or 'unknown release'}"
                 + (f"  (from {self.source_file})" if self.source_file else "")]
        if self.refused_reason:
            lines.append(f"  REFUSED: {self.refused_reason}")
            return lines
        for t in self._sorted():
            lines.append(f"  {t.describe()}")
        for v in self.uncovered:
            lines.append(f"  {v}: no dependency factor in this release, so no "
                         f"overlap to measure")
        lines.append(f"  {self.summary()}")
        return lines


def _overlap(virus: str, row: dict, alpha: float = ALPHA,
             n_tests: int = 1) -> TaxonOverlap:
    return TaxonOverlap(
        alpha=alpha,
        n_tests=n_tests,
        virus=virus,
        n_dependency=int(row.get("n_dependency") or 0),
        n_observed=int(row.get("n_observed") or 0),
        expected=float(row.get("expected") or 0.0),
        p_naive=float(row.get("p_naive") if row.get("p_naive") is not None else 1.0),
        p_matched=float(row.get("p_matched") if row.get("p_matched") is not None else 1.0),
        # Absent means the run did not measure the lower tail. It is NOT
        # defaulted to 1.0 the way the others are: a 1.0 would read as
        # "measured, and nowhere near depleted", which is a claim about a
        # test that was never run.
        p_matched_less=(float(row["p_matched_less"])
                        if row.get("p_matched_less") is not None else None),
    )


def load(results: Path, release: str | None = None) -> Landscape:
    """Read a landscape results file, refusing one that cannot be quoted.

    `release` is the release the agent is serving. A landscape measured on a
    different one is refused rather than quoted: the mis-keyed calibration bug
    in agent/tools.py was exactly this failure, where a figure measured on one
    thing was stated about another, and it is the most confident sentence
    either module can emit.
    """
    path = Path(results)
    src = path.name
    if not path.exists():
        return Landscape(source_file=str(path),
                         refused_reason=f"no results file at {path}")
    try:
        doc = json.loads(path.read_text())
    except (json.JSONDecodeError, OSError) as exc:
        return Landscape(source_file=src,
                         refused_reason=f"{src} could not be read: {exc}")

    measured_release = str(doc.get("release") or "")
    perms = int(doc.get("permutations") or 0)
    base = Landscape(
        release=measured_release,
        permutations=perms,
        seed=int(doc.get("seed") or 0),
        measured_on=str(doc.get("measured_on") or ""),
        uncovered=list(doc.get("uncovered") or []),
        source_file=src,
    )
    if release is not None and measured_release and measured_release != release:
        base.refused_reason = (
            f"measured on release {measured_release!r}, but the graph being "
            f"served is {release!r}; a landscape describes the release it was "
            f"run on and nothing else")
        return base
    if doc.get("partial"):
        why = str(doc.get("partial_reason") or "it covers a subset of taxa")
        base.refused_reason = (
            f"{src} is a partial run ({why}); every statement this module "
            f"makes is about all taxa in the release, and a subset cannot "
            f"support one")
        return base
    if perms <= 0:
        base.refused_reason = ("no permutation count is recorded, so the "
                               "matched null cannot be interpreted")
        return base
    if perms < MIN_PERMUTATIONS:
        base.refused_reason = (
            f"{perms:,} permutations cannot resolve a p below "
            f"{1 / (perms + 1):.3f}, which is too coarse to quote")
        return base

    rows = {v: r for v, r in (doc.get("taxa") or {}).items()
            if isinstance(r, dict)}

    # HOW MANY TESTS THIS RUN ACTUALLY PERFORMED, which is what the verdict
    # threshold has to answer to. One per taxon per measured direction: on
    # v0.1 that is eight taxa in two directions, sixteen tests. At an
    # uncorrected 0.05 roughly one in every other run would turn up a
    # "finding" from noise, and this module exists to stop the agent
    # claiming more than was measured -- producing the false claim itself
    # would be the same failure wearing a lab coat.
    #
    # Bonferroni rather than anything cleverer: the taxa are not independent
    # (they share the drugged proteome and much of the protein universe), and
    # a procedure assuming independence would understate the correction. It
    # is conservative for a reader asking about ONE virus, and the raw p is
    # printed beside every verdict so nothing is hidden by it.
    directions = 2 if rows and all(
        r.get("p_matched_less") is not None for r in rows.values()) else 1
    n_tests = max(1, len(rows) * directions)
    alpha = ALPHA / n_tests

    base.taxa = {v: _overlap(v, r, alpha=alpha, n_tests=n_tests)
                 for v, r in rows.items()}
    return base
