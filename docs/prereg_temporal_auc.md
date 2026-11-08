# Pre-registration: a tie-safe metric for the temporal protocol

Written before any code. Date: 2026-10-14. Release: v0.1 / v0.1-train2021,
schema v0.12.0.

The two thresholds in `docs/prereg_m1_potency.md` were both asserted without
the arithmetic and both turned out unreachable — P2 needed a fitted exponent,
and P1 demanded a tie block smaller than the unkeyed population the design
guarantees. Every number below is computed first and the working is shown.

## The defect, already measured

`rank()` sorts by score with no tiebreak and `hits_at_k` reads `ranked[:k]`,
so where a tie block straddles k the figure is decided by dict insertion
order. On the SARS-CoV-2 train2021 split:

    M1       284   8   0-8   2.8   2.49  0.0033  ~tied
    degree  4,514 232   3-5   5.1   0.97  0.5909  ~tied

Not one of M1's 8 positives scores strictly above the cut; all eight sit in
one tie block straddling position 100, with at least 8 non-positives. The
hit count ranges over the whole 0 to 8 and the lift over 0.00 to 2.86.

## What the existing code already does right

**The AUC is tie-safe and has been all along.** `labels.evaluate_against_
negatives` computes a rank-sum AUC "with ties counted as half, which matters
because DWPC produces many exact ties at percentile boundaries". Tied pairs
contribute 0.5, which is the correct treatment and makes the statistic
invariant to the order equal scores happen to be stored in.

So this is not a project-wide defect. It is confined to the metrics that
read a cut: Hits@k, MRR and the lift computed from them.

**And the temporal protocol holds everything it needs to use the safe one
and does not.** `split.test_negatives` is populated — the run reports 3,630
test inactives — and `scripts/temporal_split.py` passes only `test_pos` to
`lift_row`. It computes a tie-sensitive lift against positives alone while
holding a measured negative set it never touches and a tie-safe AUC function
it never calls.

## The change

1. The temporal protocol reports the tie-safe AUC with its Hanley-McNeil
   interval, from `evaluate_against_negatives`, alongside the existing lift.
2. `lift_row` additionally reports `hits_expected`, the hit count under
   RANDOM tie-breaking: `fixed + slots * c_pos / contested`, the
   hypergeometric mean. It is a deterministic function of the data, so it
   does not move when a dict is rebuilt, and it is the honest centre of the
   interval already reported. **It is never to be reported without
   `hits_lo`/`hits_hi` beside it**, because a point estimate at the centre of
   [0, 8] would read as a measurement.

Nothing is invented here: both pieces exist and the first is simply called.

## The arithmetic, done before the thresholds

**What E[hits] can be for M1.** The printed range pins `fixed = 0` and
`c_pos = 8` and constrains `8 <= slots <= contested - 8`; it does not pin
`contested`. Across that feasible space:

    contested   slots   E[hits]   lift
           20      12      4.80   1.70
           40      32      6.40   2.27
           60      52      6.93   2.46
          120     100      6.67   2.37
          184     100      4.35   1.54
    ... and with slots at its floor of 8, E[hits] falls as low as 0.33

E[hits] can land anywhere from 0.33 to 7.41 against an expectation of 2.82.
Those two are exact rather than sampled: the maximum is at contested=108
with slots at its cap of 100, and the minimum at contested=192 with slots at
its floor of 8. An earlier pass read the maximum off a coarse grid and got
7.29, which is the kind of small wrongness this section exists to prevent.

**The reported figure was 6, so E[hits] landing near it is a likely
outcome and not a vindication** — and it would not rescue the result either,
because the interval is still [0, 8]. Any prediction of the form "the lift
falls" would therefore be predicting something the arithmetic says is
unknowable from here. None is made.

**What AUC the interval needs to clear chance**, by Hanley-McNeil:

    n_pos   n_neg   min AUC   CI at that AUC
        8      50     0.713   [0.501, 0.925]
        8   3,630     0.704   [0.501, 0.907]
       48      50     0.612   [0.500, 0.724]
       48     200     0.593   [0.500, 0.686]
       48   3,630     0.586   [0.501, 0.671]

n_pos binds; past a few hundred negatives n_neg barely moves the bar. So the
bar is known before the run: **about 0.59 for COMBINED at its 48 positives,
and about 0.70 for M1 at its 8.** For reference the cross-sectional protocol
measured its best channels at 0.507 and 0.516 and the degree baseline at
0.582 — so COMBINED's bar sits just above the point estimate of the null
hypothesis, which is the right place for it and is not a bar this project
has ever cleared.

## Predictions

**P1 — invariance (certain, and the point of the change).** The temporal AUC
is identical when the same release is loaded with its edges in a different
order, while Hits@k on the same rows is not. If the AUC moves, the
implementation is wrong rather than the metric.

**P2 — M1 STAYS NOT EVALUABLE, stated in advance so the result cannot be
read as a disappointment.** It reaches 8 post-2021 positives against
`MIN_REACHED_POSITIVES = 10`. Adding an AUC does not change which compounds
a channel reaches, so M1 will be refused by the same gate, and even were it
evaluated its bar would be 0.704. **The headline is not rescued by this
change and nothing here should be read as trying to rescue it.**

**P3 — COMBINED becomes evaluable.** It reaches 48 post-2021 positives,
above the gate. This predicts only that a number appears, not what it is.

**P4 — determined rows are unchanged.** Where the cut lands in no tie,
`hits_expected == hits` exactly, so every determined figure in every results
file keeps its value. A change that moved determined numbers would be
changing the measurement rather than the handling of ties.

**NOT PREDICTED: the value of COMBINED's AUC.** The bar is 0.586-0.612 and
whether it is cleared is the question the measurement answers. Writing a
direction here would be the error the amendment to the previous
pre-registration was written to record.

## What cannot be concluded, whatever the numbers

- Nothing about M1, which stays not evaluable.
- A COMBINED AUC clearing its interval would be the first evaluable
  non-null in this project, and it would still carry the degree baseline's
  question: degree must be computed on the same split and compared, because
  `docs/evaluation.md` records that in every cross-sectional virus where a
  channel cleared p<0.05 the degree baseline cleared it too, 4 of 4.
- An AUC below the bar does not retire the channels; it says this release
  and this split cannot evaluate them, which is what the hard-negative
  protocol already says.

## Decision rule

Keep the AUC reporting unconditionally: a tie-safe metric replacing a
tie-sensitive one is correct independently of what it reads. Keep
`hits_expected` only if P4 holds exactly. Claim a result for COMBINED only
if its interval clears 0.5 AND the degree baseline on the same split does
not.

---

# Verdict — the tie-safe metric finds no result either

2026-10-14, measured on v0.1-train2021, both undated policies.

    SARS-CoV-2, undated=include
      scorer       hits        lift        p              AUC [95% CI]
      M1            0-8        2.49   0.0033       not evaluable (8)  ~tied
      M4              2        4.16   0.0812       0.504 [0.48, 0.53]
      M7             26        0.85      —         0.506 [0.48, 0.53]
      COMBINED        5        3.36   0.0154       0.512 [0.49, 0.54]
      degree        3-5        0.97   0.5909       0.504 [0.48, 0.53]  ~tied

    undated=exclude
      COMBINED        2        2.99   0.1428       0.505 [0.48, 0.53]
      degree          2        2.72   0.1667       0.506 [0.48, 0.53]

**P1 held** — the AUC is invariant to edge order by construction and the
implementation reproduces it.

**P2 held, as written in advance.** M1 reports "not evaluable (8)" under
both policies. The headline was not rescued and nothing here tried to.

**P3 held.** COMBINED became evaluable at 48 reached positives.

**The question that was deliberately not predicted has an answer: 0.512,
interval [0.49, 0.54], against a bar of 0.586.** It does not clear, and it
misses by 0.074. Under `exclude` it reads 0.505. No result is claimed.

## What this settles, and it is more than a null

The same scorer, on the same split, reads lift 3.36 at p=0.0154 and AUC
0.512. Those are not in tension: the lift reads a cut and the AUC does not.
**COMBINED's lift was the strongest surviving row after the determinacy
finding, and the tie-safe metric says it separates nothing.**

And the condition written in so it could not be skipped is satisfied in the
cleanest possible way. **Degree on the same split reads 0.504 [0.48, 0.53]
— at chance.** So the temporal protocol does what it was built to do: it
removes the study-volume confound that `docs/evaluation.md` records in four
of four cross-sectional viruses. Under that clean condition, no channel
separates either.

That is a better negative than "everything is confounded". It says the
protocol works and the signal is not there, which is a finding about the
graph rather than about the evaluation.

Across every virus in the run: either too few reached positives to
evaluate, or an interval containing 0.5, or a degree baseline at or above
the channel — SARS-CoV COMBINED 0.616 against degree 0.779, MERS COMBINED
0.548 against degree 0.660. **No channel on any protocol on any virus has
an AUC whose interval excludes chance.** The only scorer that ever did is
degree, cross-sectionally, and degree is the null hypothesis.

## Decision

Per the rule above: the AUC reporting is kept unconditionally, because a
tie-safe metric replacing a tie-sensitive one is correct regardless of what
it reads. `hits_expected` is kept, P4 having held exactly. **No result is
claimed for COMBINED.**

## One defect this run exposed

The `vs deg` column printed 2.55, 2.51 and 4.27 directly above a warning
reading "no comparison against this baseline can be read -- including every
ratio in the vs-deg column, which divides by it". Printing a figure and then
telling the reader to ignore it is the defect rather than the fix. The cell
now reads `n/d` when the degree row's own Hits@k is an interval.

---

# Amendment, 2026-11-08: the bar was computed on the wrong n

The verdict above is kept as written, because it records what was measured
on 2026-10-14 and the sequence is the result. Three of its statements are
no longer true of the release, and one of them was never true.

## The bar

**"At COMBINED's 48 reached positives a Hanley-McNeil interval excludes 0.5
only from 0.586 upward" used the wrong n, and the arithmetic above it is
right.** `evaluate_against_negatives` ranks the whole labelled test set and
imputes the floor to every compound a scorer cannot reach, so the AUC's n
is the test set and is identical for all nine scorers. The results file has
recorded it all along, in `auc_n_pos` and `auc_n_neg`: **1,571 and 979** on
`undated=include`. The reach count — 48 then, 96 now — is the denominator
of the LIFT, not of the AUC.

At n=1,571/979 the interval half-width is 0.023 and the bar is **0.523**.
The table in "What AUC the interval needs to clear chance" is correct for
every row it contains; the row the run needed was not in it, and the
conclusion drawn from the table therefore understated the bar by 0.063.

## What that changes

| 2026-10-14 | 2026-11-08 | |
|---|---|---|
| COMBINED 0.512 [0.49, 0.54] | **0.526 [0.503, 0.549]** | now excludes 0.5 |
| degree 0.504 [0.48, 0.53] | **0.555 [0.533, 0.578]** | now excludes 0.5 |
| COMBINED lift 3.36, p=0.0154 | 2.39, p=0.0574 | no longer significant |
| degree lift 0.97, p=0.5909 | 0.44, p=0.9723, hits 3–4 | now undetermined |

**The decision rule settles it without an amendment, which is the whole
argument for having written one.** The rule was: *claim a result for
COMBINED only if its interval clears 0.5 AND the degree baseline on the
same split does not.* The first condition is now met. The second is not —
degree clears by a wider margin, from a ranking with no imputation in it at
all. **No result is claimed for COMBINED.** Nothing here is a judgement
made after seeing the number.

## Two statements that were never true

**"Degree on the same split reads 0.504 — at chance"** is wrong on this
release and was wrong on the committed results of 2026-11-07 as well:
degree read 0.556 [0.534, 0.579] there. The sentence it licensed — "the
temporal protocol removes the study-volume confound" — is true of the LIFT,
where degree falls from 10.02 (p<.0001) cross-sectionally to 0.44
(p=0.9723), and false of the AUC. The two are not in tension for the same
reason the verdict already gives for COMBINED: the lift reads a cut and the
AUC does not. The claim must name which metric it is about.

**"The only scorer that ever did is degree, cross-sectionally"** is wrong in
both halves. Seven (run, virus, scorer) cells on this release have an AUC
interval excluding 0.5, five of them temporal and two on SARS-CoV rather
than SARS-CoV-2, and `tests/test_evaluation_table.py` now derives that
census from every results file. The claim that no INDIVIDUAL reasoning
channel clears chance anywhere survives and is asserted by the same test.

**"Four of four cross-sectional viruses"** is recorded as D2 in
`docs/audit_2026-11-03_fixes.md` and is corrected in `docs/evaluation.md`:
five cross-sectional viruses have a channel at p<0.05 and degree clears it
in one of them. The sentence counted the viruses where DEGREE was
significant.

## 0.526 IS the imputation floor

COMBINED reaches 96 of 1,575 test positives, 6.1%, and 8 of 984 test
negatives, 0.8% (the figures of the leaking graph; see the next
amendment). The other 93.9% of the positives carry the imputed floor.
An AUC computed that way exceeds 0.5 whenever a scorer reaches a larger
fraction of the positives than of the negatives, with no ranking skill
inside the reached set at all: if reached compounds all outrank unreached
ones and are internally unordered, the AUC is exactly
`0.5 + (coverage_pos - coverage_neg) / 2`.

    coverage gap            6.1% - 0.8%  =  5.26 points
    floor  0.5 + 0.0526                  =  0.5263
    measured AUC                            0.5263
    ceiling, ranking the reached set
      PERFECTLY                             0.5266

**The floor equals the measured AUC to four decimal places, and the
ceiling is three ten-thousandths above it.** COMBINED could not have
scored above 0.5266 however well it ranked, and it scored what its reach
alone produces. The interval [0.503, 0.549] excludes 0.5 because n is
1,571/979, not because anything was separated. `kgav.labels.auc_floor` and
`auc_ceiling` compute both figures and
`tests/test_auc_decomposition.py::test_the_one_chance_clearing_ensemble_cell_is_its_own_floor`
fails if a later release pulls the three numbers apart.

**The decision rule's first condition should never have been read as met.**
It said *claim a result for COMBINED only if its interval clears 0.5 AND
the degree baseline on the same split does not*. The interval clears 0.5,
and what clears it is the imputation. The rule still reaches the right
answer here, for the second reason rather than the first, but a bar stated
on an imputed AUC without a floor beside it can be cleared by a scorer that
reaches nothing well. **Any successor to this registration must state the
floor and the ceiling before the run, not only the interval bar.** That is
the one amendment to the methodology this run earns.

**This section said the opposite on 2026-11-08**, under the heading "What
is still not measured": that `temporal_split.py` "writes `auc_n_pos` and
`auc_n_neg` but not `coverage_pos` and `coverage_neg`", so the row "cannot
be distinguished from its own imputation floor" and the measurement was "a
change to the script". The script records them as `auc_coverage_pos` and
`auc_coverage_neg`. The claim was wrong by a prefix, and the thing it
called unmeasurable is the sharpest number in this file. Recorded rather
than deleted, because a registration that quietly drops a wrong statement
about its own data cannot be audited on how it found it.

---

# Amendment, 2026-11-09: the result was a leak, and the floor caught it first

The 2026-11-08 amendment above is kept as written. One day later the
similarity layer was found to be leaking post-cutoff compounds into the
2021 training graph (fix-list C4), and with that fixed **no scorer on any
undated policy has an AUC interval excluding 0.5.**

## What leaked

Every `CHEMICALLY_SIMILAR_TO` edge carries the placeholder date, because
the layer is computed rather than curated, so all 3,149 of them entered
the 2021 graph under `include`. The layer was computed over the whole
release, so those edges join compounds nothing had measured until 2024.

An undated computed edge is no older than what it was computed over.
`kgav.temporal.derived_year` now dates one by the later of its two
endpoints' first appearances -- earliest per node, because a compound
exists from its first appearance anywhere; latest across the pair, because
the pair does not exist until both do. Measured on v0.1:

    computed edges, undated                 3,173
      datable from their endpoints          2,864
      of those, endpoint after the cutoff   1,412   <- removed
      endpoints nothing can date              309   <- policy still decides

So **44% of the computed layer in the training graph was post-cutoff.**

## What it was holding up

    SARS-CoV-2 test compounds reachable   231 (14.7%)  ->  64 (4.1%)
    evaluable viruses, include             7           ->  3
    COMBINED AUC                           0.526 [0.503, 0.549]
                                                       ->  0.509 [0.486, 0.532]
    degree AUC                             0.555 [0.533, 0.578]
                                                       ->  0.516 [0.493, 0.539]
    cells with an interval excluding 0.5   5           ->  0

`exclude` moved the other way, 1 evaluable virus to 3, because 1,452
computed edges are now genuinely dated pre-cutoff instead of being
discarded as dateless. A fix that only ever removed evidence would be
suspicious; this one removes post-cutoff edges and admits pre-cutoff ones,
which is what dating them means.

## The decomposition got there first, by a different route

The 2026-11-08 amendment had already established that COMBINED's 0.526 was
its own imputation floor to four decimals, with a ceiling of 0.5266 --
that it cleared the bar on a figure its reach produced by itself. The leak
fix now puts it below the bar. **Two independent findings, one day apart,
and either alone retires the result.** That is the argument for reporting a
floor beside every imputed AUC: the decomposition did not need to know why
the coverage gap was there.

Across all 69 temporal AUC cells on this release, 53 equal their floor to
four decimal places, the median deviation is exactly zero, the largest is
0.0067, and the largest floor-to-ceiling span any scorer had is 0.0320.
**The tie-safe AUC on this protocol cannot separate**, whatever scorer is
applied to it: a channel could rank its reached set perfectly and still
read 0.509.

## The methodology amendment stands and gets sharper

The 2026-11-08 amendment required that any successor registration state
the floor and the ceiling before the run, not only the interval bar. That
is now the binding requirement rather than a suggestion, and it is not
enough on its own. **A successor must also state the minimum reach at
which the AUC will be read**, because a bar on an imputed AUC is cleared
by coverage asymmetry alone and the floor only reveals that afterwards.
`MIN_REACHED_POSITIVES = 10` gates the LIFT and gates the AUC only through
`auc_evaluable`; on this split COMBINED passes it at 32 reached positives
and still has 0.0000 of span between floor and ceiling.

## Decision

Unchanged, and now for the simplest possible reason. **No result is
claimed for COMBINED**: its interval contains 0.5.

---

# Amendment, 2026-11-09 (second): the test-set sizes moved under L3

The figures above were measured before fix-list L3, which made a label
edge whose potency fields DISAGREE decide nothing. 173 such edges left the
label set, so the temporal test set is slightly smaller:

    auc_n_pos  1,580 -> 1,571        auc_n_neg  982 -> 979

The bar is unchanged at **0.523** -- the Hanley-McNeil half-width at
1,571/979 is still 0.023 -- and so is every conclusion above. COMBINED
reads 0.5087 against a floor of 0.5087, degree 0.5163 against 0.5163, and
neither interval excludes 0.5. Recorded because an n that moves without
being re-derived is how the 0.586 bar came to be computed on the wrong
one in the first place.
