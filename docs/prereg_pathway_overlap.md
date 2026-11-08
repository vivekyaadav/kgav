# Pre-registration: the pathway and interactome overlap, re-run

**UNSIGNED DRAFT.** The values under *Fixed in advance* are proposals with
their arithmetic attached. This is not a registration until it is signed and
committed, and the re-run must not be executed before then. Two decisions
are flagged **DECISION** and belong to the author, not to this draft: the
permutation count, and whether the six empty-rejection-region tests stay in
the correction family. A third candidate -- whether `auc_ceiling` replaces
`MIN_REACHED_POSITIVES` as the readability gate -- belongs to the
hard-negatives analysis and not to this one, and is left out rather than
listed here for completeness.

Date drafted: 2026-11-02. Release: v0.1 (78,789 nodes / 392,674 edges,
`nodes_sha` 3356fae68afdf18d, `edges_sha` f3bb418a89779e94). Scripts:
`scripts/pathway_landscape.py`, read by `src/kgav/landscape.py`.

**The release identification was updated on 2026-11-08 and the analysis was
not re-registered**, because the rebuild of 2026-11-07 did not touch it.
`data/results/pathway_landscape.json` was regenerated against the new
release and differs from the committed version in **four fields, all
provenance**: the code commit, the dirty flag, the timestamp and
`edges_sha`. Every measured figure -- every p-value, every overlap ratio,
every degree-matched null -- is byte-identical. The two edges this rebuild
added are `HAS_ANTIVIRAL_ACTIVITY_AGAINST`, which this analysis does not
read. A registration that quoted a release no longer on disk would be
unverifiable, which is what `tests/test_landscape.py::
test_the_quoted_release_matches_the_census` exists to prevent; it caught
this.

## Why this is a re-registration, and what that costs

The other four pre-registrations in this directory open by saying they were
written **before** the code, "because a prediction made afterwards is a
description". This one cannot say that. The analysis ran on 2026-10-28 and
its figures are quoted in `docs/scope_and_limits.md` and in the paper.

**So the 2026-10-28 run was unregistered and its figures are exploratory.**
Writing registered predictions for a measurement already in hand would be
the exact failure this directory's format exists to prevent, and no
paragraph here retroactively fixes that. What is genuinely still ahead is a
**re-run**, forced by two defects found on 2026-11-02 that change every
p-value in the file, and whose parameters are genuinely undecided. This
document registers that re-run.

The published null does not change status because of this. It was
exploratory before and it is exploratory now; what changes is that the
re-run will not be.

## The defects forcing the re-run

**One generator consumed in row order.** `matched_null` already documents
three bugs of this kind -- unsorted pools, unsorted bin order, degree ties
broken by set iteration -- and each was the same sentence: the seed stopped
determining the draw. All three were fixed and a single run-wide generator
survived them. Measured on the real function with `Random(0)` and identical
data, one taxon's null total:

    rows A, B        A = 1209
    rows C, A, B     A = 1165      a virus added BEFORE it
    rows B, A        A = 1165      the same two, other order

Adding a taxon to `config/viruses.yaml`, reordering the register, or running
one channel alone moves every downstream p-value. The curated layer will add
five viruses with dependency factors, so this is not hypothetical.

**Two resolvable steps below the threshold.** At 2,000 permutations an
empirical p takes values k/2001, so below the corrected 0.001136 there are
exactly two attainable values -- 0.000500 and 0.001000 -- and the third,
0.001499, is above it. The test could distinguish "significant" from "not"
and nothing finer.

Both are fixed in code. The figures in the file still come from the old
behaviour, and `landscape.load` reports that: a file with no `seed_scheme`
key now describes itself as shared-generator in the protocol string.

## The measurement

Unchanged, and stated here because a re-run that quietly changes the
estimand is a different study. For each (channel, taxon): the overlap
between the drugged proteome and that taxon's factor set, against a
**degree-matched permutation null** that draws, from each degree bin, as
many proteins as the drugged set has in that bin.

    pathway channel       factors = CRISPR dependency factors
    interactome channel   factors = host proteins a viral protein touches

**The estimand is the overlap count**, compared against its matched
expectation. Both tails from the same null.

## Fixed in advance

- **Permutations: 100,000.** **DECISION**. The hard floor is one attainable
  value below alpha, which needs N+1 > 1/alpha = 880. Ten steps needs
  N >= 8,799. 100,000 gives 113 steps and costs 50x the current runtime.
  The proposal is 100,000 rather than 8,799 because the family size will
  change when the curated layer lands, and a count chosen to just clear
  today's alpha has to be revisited then.

- **Correction family: every test the run performs, both channels, both
  directions.** On v0.1 that is 22 rows x 2 = 44, so alpha = 0.05/44 =
  0.001136. Computed by the script, which is the only thing that knows how
  many channels ran; `landscape.load` takes the declared count as a floor
  and may raise it, never lower it.

- **The six empty-rejection-region tests stay in the family.** **DECISION**.
  Six of the 22 downward tests cannot reject at any observation: the most
  extreme downward value is zero hits with probability exp(-m), and where
  that exceeds alpha -- m < 6.78 -- no observation clears the bar.

      interactome  JEV          0.83 expected
      interactome  West Nile    0.87
      pathway      West Nile    2.55
      interactome  DENV-3       4.41
      interactome  TBEV         4.59
      pathway      DENV-2       6.40

  Dropping them would give 0.05/38 = 0.001316 and is defensible: emptiness
  depends on the matched expectation, which the degree-matched design fixes
  before any outcome is looked at, so it is not choosing tests after seeing
  them. The proposal is to KEEP them, for two reasons. The threshold moves
  16%, from 0.001136 to 0.001316, and the smallest p in the family is
  0.003998 -- 3.0x the looser bar -- so nothing turns on it. And a family
  defined as "the tests that could have rejected" has to be recomputed
  whenever an expectation moves, which makes the correction depend on the
  release in a way a reader cannot check by counting.

- **Both procedures reported: Bonferroni and Benjamini-Hochberg at
  q = 0.05.** On the 2026-11-07 figures both reject nothing, and the
  smallest p is 0.003998 against a Bonferroni bar of 0.001136 and a BH
  critical value of 0.001136 at k=1. Bonferroni over 22 correlated rows is
  over-conservative and that objection is worth answering with a number
  rather than a sentence. Primary is Bonferroni; BH is reported beside it.

- **One generator per row, seeded `Random(f"{seed}|{channel}|{virus}")`,
  default seed 0.** A string seed because `Random` hashes a `str` through
  sha512 of its bytes rather than through `hash()`, so `PYTHONHASHSEED`
  cannot move it -- the property the three earlier fixes were about,
  asserted in subprocesses rather than reasoned about. The scheme is
  recorded in the output file, because a file written under either scheme
  is structurally identical and differs in every p-value.

- **Downward readability gate: `expected > -ln(alpha)`, derived from the
  run's own alpha.** It was a constant 3.0, from `exp(-3) = 0.0498 < 0.05`
  -- the uncorrected alpha -- while every verdict in the module is judged
  against alpha/n_tests. Derived per row, so a change in family size moves
  it automatically.

- **Primary direction: enrichment. Secondary: depletion.** Both are
  computed from the same null and both are in the correction family. The
  asymmetry is in interpretation only: the project's claim is about
  pathway-mediated routes existing, and depletion is reported because
  without it a low ratio cannot be told from an ordinary low draw.

- **Degree basis: `--degree-from study`**, the default, not the assembled
  graph. The graph basis is circular -- it derives degree from the same
  TARGETS and virus-host edges that define the sets being compared -- and
  biases every ratio downward.

- **No floor on taxon size for being tested.** A taxon with a small
  expectation can still return a significant excess, and excluding it would
  be a power decision disguised as a coverage decision. The minimum
  detectable effect is reported per row instead, in
  `docs/scope_and_limits.md`.

## Registered predictions

**P1. No test clears the corrected alpha, under Bonferroni or under BH.**
The null this document points at, stated so that failing to reject it is a
result rather than an absence.

**P2. Raising the permutation count will not produce a rejection.** The
smallest p in the 2026-11-07 family is 0.003998, which is 3.5x the
Bonferroni threshold. Resolution was never what blocked a finding, and this
prediction says so in advance: if a rejection appears at 100,000
permutations it did not come from finer resolution and something else
changed.

**P3. Per-row seeding will move individual p-values and not the verdict.**
Every row's draws change, so every p changes. The ordering of ratios and
the set of rejections should not. A verdict that flips on the seeding scheme
alone would mean the result was a draw rather than a measurement.

**P4. The six empty-rejection-region rows stay empty.** Emptiness is a
property of alpha and the expectation, not of N, so no permutation count
fixes it. If the curated layer later gives DENV-3, TBEV or JEV real
dependency-factor sets, their expectations rise and some of the six leave
that state -- which is the only thing that can.

## What would change the conclusion

A test clearing the corrected alpha under both procedures, on a row whose
minimum detectable effect is below the observed ratio. That would be the
first pathway-level finding in this project and would justify a per-pathway
test rather than a pooled overlap.

## What this cannot do

- It does not test pathway ENRICHMENT in the gene-set sense. The statistic
  is an overlap between two flat protein sets against a degree-matched
  background. `scripts/pathway_landscape.py` computes no per-pathway
  quantity, and the claim in the paper is being renamed accordingly (P1 in
  the audit fix list) rather than defended here.
- It does not establish mechanism for any row it cannot reject, and seven of
  22 rows need a 2.4x excess or more before they could.
- It does not settle `MIN_REACHED_POSITIVES` or `MIN_CLASS_SIZE`. Both are
  10; `MIN_REACHED_POSITIVES` justifies itself by matching
  `MIN_CLASS_SIZE`, and `MIN_CLASS_SIZE` at `scripts/hard_negatives.py:29`
  carries no comment at all. The pair is circular and neither is used by
  this analysis, so it is named here and left to its own registration.
- It does not license pooling the two channels. They ask different
  questions of the same virus -- what it NEEDS against what it TOUCHES --
  and are reported side by side.

## Implementation constraints

No scipy. The Poisson tail behind the minimum-detectable-effect table is
hand-written with an iterative pmf, because m = 425 overflows a factorial,
and is checked against a direct summation at sizes where the factorial
still fits.

100,000 permutations over 22 rows is 50x the current run. If that is not
affordable the fallback is 8,799, which clears ten steps at today's family
size and must be revisited when the family changes -- and the fallback has
to be chosen HERE, not after seeing how long the first attempt takes.
