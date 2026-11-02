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
import math
from dataclasses import dataclass, field
from pathlib import Path

ALPHA = 0.05


def min_expected_hits(alpha: float) -> float:
    """The matched expectation below which the RATIO carries no information
    in the low direction, for a test judged at `alpha`.

    Under a Poisson background with mean m the most extreme downward
    observation is zero hits, with probability exp(-m). So the ratio is
    readable downward only where exp(-m) < alpha, i.e. m > -ln(alpha).
    Derived, not chosen: a threshold picked for looking round would be a
    fitted parameter with no record of what fitted it.

    IT TAKES ALPHA AS AN ARGUMENT, and it used to be a module constant of
    3.0 -- because exp(-3) = 0.0498 < 0.05. That is the UNCORRECTED
    alpha. Every other verdict in this module is judged against
    ALPHA / n_tests, which on release v0.1 is 0.05/44 = 0.001136, and the
    threshold against that is -ln(0.001136) = 6.78. So one taxon was
    passing as readable-downward while a zero overlap could never have
    cleared the bar its own row quotes: DENV-2, 6.22 expected.

    TaxonOverlap already carried `alpha` per row, with a comment saying
    it is held there "so a verdict can never be computed against a
    different threshold from the one reported next to it" -- and
    ratio_readable reached past it for the global. The data was right and
    the one property that needed it looked somewhere else.

    (The external audit named DENV-2, DENV-3 and TBEV. DENV-3 and TBEV
    have no dependency factors, so they have no row in this analysis and
    were never tested in either direction. One taxon, not three.)
    """
    if not 0.0 < alpha < 1.0:
        raise ValueError(f"alpha must be in (0, 1), got {alpha!r}")
    return -math.log(alpha)


def resolvable_steps(permutations: int, alpha: float) -> int:
    """How many attainable p-values sit strictly below `alpha`.

    An empirical p from N permutations takes values k/(N+1), so the
    number below alpha is floor(alpha * (N + 1)) less the k = 0 case that
    cannot occur. ZERO means the test cannot reject at any observation
    whatsoever, however extreme.

    THE SAME ERROR AS ABOVE, IN THE SAME MODULE. MIN_PERMUTATIONS = 200
    was justified as "the finest p is 0.005, which is enough for a 0.05
    test" -- the uncorrected alpha again. At the corrected 0.001136 a
    200-permutation run resolves ZERO steps below the threshold and the
    old floor accepted it. On release v0.1, with 2,000 permutations:

        p floor         1/2001   = 0.000500
        attainable      0.000500, 0.001000          <- below alpha
        next            0.001499                    <- above alpha
        corrected alpha 0.001136

    Two values. A test with two resolvable steps below its own threshold
    can report "nothing was significant" honestly -- failing to reject is
    still failing to reject -- but it cannot distinguish a true p of
    0.0012 from one of 0.0005, and its negative result is weaker than the
    four decimal places it is quoted to.
    """
    if permutations <= 0:
        return 0
    return max(0, int(math.floor(alpha * (permutations + 1) + 1e-12)))


def minimum_detectable_ratio(expected: float, alpha: float,
                             permutations: int = 0) -> tuple[int, float, str]:
    """`(smallest detectable count, that count over expected, what binds)`.

    THE TABLE A NULL RESULT CANNOT BE PUBLISHED WITHOUT, and this project
    had it nowhere. "Nothing survived correction" is only a finding if
    the reader knows what WOULD have survived: a row expected to hit 0.84
    pathways needs six hits -- a 7.1x excess -- before it can clear the
    corrected threshold, and failing to find that is not evidence against
    a real effect of 2x.

    Poisson upper tail, which is the same approximation
    min_expected_hits already uses for the lower one, computed with an
    iterative pmf because m = 425 overflows a factorial.

    TWO THINGS CAN BIND, and reporting only the first would hide the
    second. The Poisson tail gives the smallest k whose p falls below
    alpha; the permutation grid cannot express a p below 1/(N+1) at all.
    Where the floor is the coarser of the two, no k is detectable and
    "floor" is returned -- a row whose minimum detectable effect is
    infinite because of how the run was configured, not because of the
    biology.
    """
    if not 0.0 < alpha < 1.0:
        raise ValueError(f"alpha must be in (0, 1), got {alpha!r}")
    if permutations > 0 and 1.0 / (permutations + 1) >= alpha:
        return 0, float("inf"), "floor"
    if expected <= 0:
        return 1, float("inf"), "empty-background"
    pmf = math.exp(-expected)
    cum = pmf
    k = 0
    while k < 100_000:
        k += 1
        if 1.0 - cum < alpha:
            return k, k / expected, "poisson"
        pmf *= expected / k
        cum += pmf
    raise AssertionError(f"no detectable count for expected={expected}")


#: Resolvable steps below alpha under which a run is reported as
#: underpowered rather than quoted plainly. Ten is a convention and is
#: labelled as one: the arithmetic only gives the hard floor, which is
#: one step, and a test sitting at the hard floor is one permutation draw
#: away from being unable to reject at all.
MIN_RESOLVABLE_STEPS = 10

# An absolute sanity floor, independent of alpha: below this a run is
# refused whatever its family size. Kept because it catches a truncated or
# aborted run, which the alpha-derived check would pass for a large alpha.
MIN_PERMUTATIONS = 200


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
    #: Carried for the same reason as `alpha`: the minimum detectable
    #: effect depends on the permutation grid as well as the threshold,
    #: and a row computing it from a module default would report a
    #: sensitivity the run never had.
    n_permutations: int = 0

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

        AGAINST THIS ROW'S OWN ALPHA. It was the module constant, derived
        from the uncorrected 0.05 -- see min_expected_hits.
        """
        return self.expected >= self.min_expected_hits

    @property
    def min_expected_hits(self) -> float:
        """The readability threshold for THIS row's corrected alpha."""
        return min_expected_hits(self.alpha)

    @property
    def depletion_can_reject(self) -> bool:
        """Whether a DOWNWARD test on this row could reject at all.

        The most extreme downward observation is zero, so the test has a
        non-empty rejection region only where exp(-expected) < alpha.
        Below that no observation whatsoever clears the bar -- the test
        is in the correction family and cannot contribute a finding to
        it, which makes every other test in the family harder to pass
        for nothing.
        """
        return self.expected > self.min_expected_hits

    @property
    def minimum_detectable_ratio(self) -> float:
        """The smallest UPWARD excess this row could have detected."""
        return minimum_detectable_ratio(
            self.expected, self.alpha, self.n_permutations)[1]

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
                     f"{self.min_expected_hits:.1f}, even a zero overlap would be "
                     f"an ordinary draw, so the ratio is not readable "
                     f"downward")
        return head


@dataclass(frozen=True)
class FamilyRow:
    """One test in the correction family, for POWER reporting only.

    Deliberately thin and deliberately not a TaxonOverlap: it carries no
    p-value and no verdict, because its only job is to answer "could
    this test have rejected, and at what effect size". A record that
    could be mistaken for a result would be read as one.
    """
    channel: str
    virus: str
    expected: float
    n_observed: int
    n_dependency: int
    alpha: float
    n_permutations: int = 0

    @property
    def depletion_can_reject(self) -> bool:
        return self.expected > min_expected_hits(self.alpha)

    @property
    def minimum_detectable(self) -> tuple[int, float, str]:
        return minimum_detectable_ratio(self.expected, self.alpha,
                                        self.n_permutations)


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
    #: Tests the run performed in OTHER channels of the same file. The
    #: correction answers to the family of tests performed, not to the
    #: subset this module reports -- see n_tests.
    other_channel_tests: int = 0
    #: Every row in the correction family, both channels, for POWER
    #: reporting only. See the loader: these are not results and never
    #: reach for_virus.
    family: list[FamilyRow] = field(default_factory=list)

    @property
    def usable(self) -> bool:
        return not self.refused_reason

    @property
    def protocol(self) -> str:
        base = (f"a degree-matched permutation null, "
                f"{self.permutations:,} permutations, seed {self.seed}, "
                f"significance at {self.alpha:.4f} "
                f"({ALPHA} Bonferroni-corrected for {self.n_tests} tests)")
        # THE RESOLUTION, BESIDE THE THRESHOLD. Quoting a four-decimal
        # alpha next to a permutation count that cannot resolve it is
        # the false precision this whole module is written against.
        if self.underpowered:
            base += (f". Only {self.resolvable_steps} attainable p-value(s) "
                     f"lie below that threshold at {self.permutations:,} "
                     f"permutations, so a negative result here means the "
                     f"test could not reject, not that the effect is small; "
                     f"{int(math.ceil(MIN_RESOLVABLE_STEPS / self.alpha)) - 1:,} "
                     f"permutations would give {MIN_RESOLVABLE_STEPS}")
        return base

    @property
    def n_tests(self) -> int:
        """Every test the RUN performed, not every test this module reports.

        One per taxon per measured direction, plus the tests any other
        channel in the same file performed. The interactome channel
        landed 2026-10-28 and measures fourteen taxa in two directions:
        counting only this module's eight would leave the threshold at
        0.0031 when the family is forty-four tests and the bar is
        0.0011. The agent would then have called a p of 0.0035
        significant on a run where nothing was.

        A correction answers to what was DONE, not to what gets quoted.
        Leaving the other channel out because the agent cannot cite it
        would be choosing the tests after seeing them.
        """
        rows = self._sorted()
        mine = len(rows) * (2 if self.depletion_tested() else 1)
        return max(1, mine + max(0, self.other_channel_tests))

    @property
    def alpha(self) -> float:
        """The corrected threshold every verdict here is judged against."""
        return ALPHA / self.n_tests

    @property
    def min_expected_hits(self) -> float:
        """The downward-readability threshold for this run's own alpha."""
        return min_expected_hits(self.alpha)

    @property
    def resolvable_steps(self) -> int:
        """Attainable p-values strictly below this run's corrected alpha."""
        return resolvable_steps(self.permutations, self.alpha)

    def powerless_depletion(self) -> list[FamilyRow]:
        """Family rows whose DOWNWARD test has an empty rejection region.

        On release v0.1 this is six of twenty-two, and the external
        audit's count is exactly right:

            interactome  JEV          0.84 expected
            interactome  West Nile    0.88
            pathway      West Nile    2.53
            interactome  DENV-3       4.44
            interactome  TBEV         4.60
            pathway      DENV-2       6.22

        WHAT IS NOT DONE HERE. These six are still in the family, so the
        threshold every other test is judged against is 0.05/44 rather
        than 0.05/38 -- six tests that cannot contribute a finding make
        the sixteen that can harder to pass. Dropping them is defensible:
        emptiness depends on the matched EXPECTATION, which is fixed by
        the degree-matched design and knowable before any outcome is
        looked at, so it is not choosing tests after seeing them.

        It is still not a change to make in passing. This project's own
        standing rule is that the correction family is fixed in the
        preregistration before a p-value is in hand, and there is no
        preregistration yet. So the rows are NAMED and the family is
        left alone, and the decision belongs in that document.
        """
        return [r for r in self.family if not r.depletion_can_reject]

    @property
    def underpowered(self) -> bool:
        """Whether the permutation count can resolve this run's own alpha.

        SEPARATE FROM `usable`, deliberately. A run with too few
        resolvable steps has still measured something: failing to reject
        is failing to reject, and "nothing survived correction" is the
        finding this project is publishing. What it cannot do is
        distinguish a p just under the threshold from one far below it,
        so the negative result is weaker than its four decimal places
        suggest. Refusing the file would delete a true finding to punish
        its precision; reporting the resolution lets a reader weigh it.

        A run resolving ZERO steps is different and IS refused by the
        loader: it could not have rejected at any observation, so
        quoting it as a null test asserts a test that never happened.
        """
        return self.resolvable_steps < MIN_RESOLVABLE_STEPS

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
                     f"{self.min_expected_hits:.1f}, where even a zero overlap is "
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
             n_tests: int = 1, n_permutations: int = 0) -> TaxonOverlap:
    return TaxonOverlap(
        alpha=alpha,
        n_tests=n_tests,
        n_permutations=n_permutations,
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

    # TESTS THE RUN PERFORMED ELSEWHERE IN THE SAME FILE. The writer
    # records its own total, because only the writer knows how many
    # channels it ran -- a reader that inferred the count from keys it
    # recognises would silently miss a channel added later, which is
    # exactly how this module came to be correcting for eight tests on a
    # run of forty-four.
    #
    # Taken as a FLOOR, never as the answer: max() against what can be
    # seen here means a file cannot weaken the threshold by
    # under-reporting, only strengthen it by reporting more than this
    # module can account for.
    mine = len(rows) * directions
    declared = int(doc.get("n_tests") or 0)
    base.other_channel_tests = max(0, declared - mine)
    # COMPUTED HERE, not read back off base.n_tests: that property counts
    # base.taxa, which is still empty on this line. Reading it gave every
    # TaxonOverlap alpha 0.05 and n_tests 1 -- an UNCORRECTED threshold,
    # the exact failure this function exists to prevent, introduced while
    # tightening it.
    n_tests = max(1, mine + base.other_channel_tests)
    alpha = ALPHA / n_tests

    # A RUN THAT CANNOT REJECT AT ANY OBSERVATION IS NOT A TEST.
    # MIN_PERMUTATIONS above is an absolute floor derived against the
    # UNCORRECTED 0.05 -- "at 200 the finest p is 0.005, which is enough
    # for a 0.05 test". Against a corrected 0.001136 a 200-permutation
    # run resolves ZERO attainable p-values below the threshold, and that
    # floor accepted it. This check is the one that answers to the alpha
    # the verdicts are actually judged against, so it has to sit here,
    # after n_tests is known, rather than beside the constant.
    steps = resolvable_steps(perms, alpha)
    if steps < 1:
        base.refused_reason = (
            f"{perms:,} permutations cannot resolve the corrected threshold "
            f"of {alpha:.6f} ({ALPHA} over {n_tests} tests): the finest "
            f"attainable p is {1 / (perms + 1):.6f} and no attainable value "
            f"lies below it, so this run could not have rejected at any "
            f"observation. "
            f"{int(math.ceil(MIN_RESOLVABLE_STEPS / alpha)) - 1:,} "
            f"permutations would give {MIN_RESOLVABLE_STEPS} steps")
        return base

    base.taxa = {v: _overlap(v, r, alpha=alpha, n_tests=n_tests,
                             n_permutations=perms)
                 for v, r in rows.items()}
    # EVERY ROW IN THE FAMILY, both channels, for power reporting only.
    #
    # The module corrects alpha for the interactome channel's fourteen
    # taxa -- see n_tests -- and then loads none of them, so fourteen of
    # the twenty-two rows paying the multiple-testing price were
    # invisible to the module charging it. They are NOT in `taxa` and
    # never reach for_virus: an interactome overlap reported as a pathway
    # overlap would be a wrong answer, which is worse than a missing one.
    # They are here so the power table can cover the real family.
    base.family = []
    for channel, src in (("pathway", rows),
                         ("interactome", doc.get("interactome") or {})):
        for v, r in sorted((src or {}).items()):
            if not isinstance(r, dict):
                continue
            base.family.append(FamilyRow(
                channel=channel, virus=v,
                expected=float(r.get("expected") or 0.0),
                n_observed=int(r.get("n_observed") or 0),
                n_dependency=int(r.get("n_dependency") or 0),
                alpha=alpha, n_permutations=perms))
    return base
