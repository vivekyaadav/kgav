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
