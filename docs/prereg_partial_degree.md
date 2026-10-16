# Pre-registration: does any channel carry signal beyond degree

Written **before** the code, and committed before the measurement, because a
prediction made afterwards is a description. If the result contradicts what
is written here, this file stays as it is and the verdict goes beneath it.

Date: 2026-10-16. Release: v0.1 (74,935 nodes / 374,792 edges, schema
v0.13.0). Scorings: `baselines.SCORING_LEGACY` and `SCORING_PDP`.

## The defect

**Every verdict in `docs/evaluation.md` compares a channel against a degree
baseline, and both scorings of that channel are degree with a sign.**

Measured by `scripts/diagnose_scorer.py` on v0.1, Spearman between a
channel's score and the compound's degree in the layer the channel walks:

    legacy    +0.3 to +0.7 on most cells; M4 reaches +0.831
    pdp       negative on 25 of 36 channel-virus cells, to -0.846
    agreement the global-degree and layer-degree columns agree in SIGN
              36 times out of 36, so the flip is not an artifact of which
              degree it is measured against

Four cells read exactly **-1.000**: M8 on SARS-CoV-2, Zika and SARS-CoV, and
M1 on HCoV-229E. Rank for rank, those channels ARE the compound's edge
count.

A comparison of two quantities that are each a monotone function of a third
cannot say which one explains the labels. `confound_warning()` reports
whether the degree baseline is itself significant, which is necessary and
not sufficient: it cannot separate a channel that adds mechanism from a
channel that restates degree, because it never puts both in one model.

## The measurement

For each (virus, channel) cell, a logistic regression on the channel's own
evaluated pool:

    label  ~  b_c * channel  +  b_d * degree  +  intercept

fitted twice, once per scoring. **`b_c` is the estimand.** It is what the
channel contributes with degree already accounted for, which is the question
every table in `evaluation.md` has been asking by hand.

This is Rephetio's design reduced to two features. That model carried "a
prior probability of treatment, node degrees for 14 metaedges, and DWPCs for
123 metapaths" over 755 positives; this project has 0 to 904 positives per
virus and 8 channels, so the full model is not yet fittable and the pairwise
one is. The ordering is deliberate: if no channel clears this, there is
nothing for a larger model to find.

### Fixed in advance

- **Features are midrank percentiles**, from `baselines.channel_percentiles`
  rather than a second copy of its arithmetic, then centred and scaled.
  Raw DWPC is skewed and tie-heavy enough that a z-score is dominated by its
  tail, and midranks are already the form `combine` uses.
- **Degree is the same global node degree `degree_ranking` returns**, the
  one the published tables compare against, percentiled the same way.
- **L2 penalty, lambda = 1.0**, on the two coefficients and not the
  intercept. Ridge, not lasso: selection over two features at these counts
  is noise, and Rephetio could afford elastic net at 755 positives.
- **Test: likelihood-ratio**, the two-feature model against degree-only,
  both penalised identically, one degree of freedom. Not a Wald test on
  `b_c`, whose standard error under penalisation is approximate.
- **alpha = 0.05, Bonferroni-corrected across every cell actually run**,
  both scorings counted. Raw p is reported beside the corrected one.
- **Floor: 10 positives AND 10 negatives in the pool**, matching
  `calibration.MIN_REACHED_POSITIVES`. A cell below it is NOT RUN and says
  so; it is not a null.
- **Primary is cross-sectional**, for sample size. Secondary is the temporal
  `computed_only` split, which is publication-disjoint by construction.

### Why cross-sectional is defensible here and is not elsewhere

92.2% of compounds share a PMID between their `INHIBITS` edge and their
label, which is why raw cross-sectional lift means little. Degree is that
confound's best available proxy, and this model controls for it. So `b_c`
is more robust to the shared-PMID problem than the lift it replaces --
**partially, not fully**: a shared publication is not only degree, and
nothing here makes the cross-sectional protocol clean.

## Registered predictions

**P1. No cell clears the corrected alpha under both scorings with the same
sign.** This is the null the whole document points at, stated so that
failing to reject it is a result rather than an absence.

**P2. Where a cell IS significant under one scoring, it will not be sign-
consistent across the two.** This is the falsification test and the reason
both scorings are run. A channel carrying real mechanism should contribute
in the same direction however the degree damping is parameterised; a
"signal" that reverses when `w` moves from 0 to 0.4 on the source term was
degree in one direction or the other.

**P3. M8 everywhere and M1 on HCoV-229E return NOT RUN or null.** Their pdp
score is a perfect monotone function of the covariate, so the two features
are collinear and `b_c` is unidentifiable. They are expected to fail the
floor first -- 11, 29, 11 and 17 drugs reached -- which is itself the point.

**P4. M4 produces the largest `|b_d|` of any channel.** It is the most
degree-correlated channel on every virus that reaches it, under both
scorings and both degrees.

## What would change the conclusion

A cell significant under both scorings, same sign, surviving Bonferroni,
with the floor met. That would be the first evidence in this project of a
channel carrying something degree does not, and it would justify the full
regularized model.

## What this cannot do

- It does not rank drugs and will not be used to. The output is a
  coefficient and a p-value per cell.
- It does not establish mechanism. Controlling for degree removes one
  confound, named and measured; the label-provenance confound is only
  partly proxied by it.
- It does not choose a scoring. If P2 holds, neither scoring is preferable
  and the question of which to default to stays open.
- With 2 features and a floor of 10 positives the events-per-variable is 5,
  which is low. A null here is weak evidence of absence, and the verdict
  must say so rather than reading P1's confirmation as a finding.

## Implementation constraints

`pyproject` declares `pyyaml`, `pandas`, `pyarrow`. There is no scipy and no
sklearn, so the fit and the chi-square tail are hand-written, as
`labels.auc_interval` and `phase4_hcv_retrospective.fisher_one_sided`
already are. Both are to be checked against sklearn's
`LogisticRegression(penalty="l2", C=1/lambda)` and `scipy.stats.chi2` over
randomised inputs, with the agreement recorded in the test file and the
reference values pasted rather than recalled.

Separation is expected on small cells. The penalty bounds the coefficients,
and a cell whose fit does not converge is reported as NOT FITTED rather than
given a coefficient.

# Verdict — P1 not rejected, P2 held, and a defect in the registered statistic

Run 2026-10-17 on v0.1 by `scripts/partial_degree.py`, both scorings, both
readings of the undefined "negative".

## P1 is not rejected

    negatives     cells  tested  significant after Bonferroni
    measured         78      16  0
    all              78      26  1   (SARS-CoV-2 M1, legacy only)

**No cell clears the corrected alpha under both scorings with the same
sign, under either reading.** P1 stands.

## P2 held, on the one cell that could test it

SARS-CoV-2 M1 under the all-negatives reading is significant in the legacy
scoring (`b_c` = −0.201, p = 0.0006, corrected 0.0168) and absent in pdp
(`b_c` = −0.022, p = 0.7033). A contribution that survives one
parameterisation of the degree damping and vanishes under the other is the
signature P2 named in advance.

## The closest thing to a signal, which does not clear the bar

SARS-CoV-2 M4, measured negatives:

    scoring              b_c      b_d       LR        p    p_adj
    entered-node-global  +0.401  -0.295   4.126   0.0422   0.6758
    path-degree-product  +0.403  -0.215   6.765   0.0093   0.1488

Same sign, near-identical magnitude, significant uncorrected under both,
and the only cell in the table with that property. It fails Bonferroni and
is therefore **not a result**; it is the one cell worth powering up. Note
that M4 is also the most degree-correlated channel in
`diagnose_scorer.py` -- and here, with degree in the model, `b_c` is
positive while `b_d` goes negative.

## The reading of "negative" decides what the model measures

`b_d`, SARS-CoV-2, by how a negative is defined:

    channel   measured negatives   unlabelled counted as negative
    M1        -0.393  -0.383       +1.972  +1.956
    M4        -0.295  -0.215       +0.560  +0.406
    M5        -0.218  -0.203       +0.325  +0.284

Negative on 5 of 7 channels under measured negatives; positive on 6 of 7
under the other, reaching +2.23 on SARS-CoV. **Counting unlabelled
compounds as inactive turns the model into a predictor of which compounds
got measured at all**, which is the 92.2% shared-PMID confound in another
costume. The GAP resolution taken in the script -- measured negatives as
primary -- was argued from principle before the run and is now supported by
the run.

## A defect in the registered statistic, found at implementation

**The registered likelihood ratio can be negative.** MERS-CoV M4 under pdp
returned −0.031 on this release.

The registration says the two models are "both penalised identically",
meaning the same lambda -- but over different numbers of coefficients. The
full model maximises `loglik - lambda/2 * (b_c^2 + b_d^2)`, the reduced one
`loglik - lambda/2 * b_d^2`, so admitting `b_c` lets the fit trade `b_d`
away against the penalty and settle at a lower DATA log-likelihood than the
reduced model reached. Reproduced in 146 of 300 synthetic null cells shaped
like the real ones (n ≈ 3,000, 15 positives, features correlated at 0.95),
most negative −0.601.

**It is conservative.** Those cells get p = 1, so the defect can only fail
to reject. P1 is a null and its non-rejection survives. What it costs is
power, precisely where channel and degree are collinear -- which is every
cell worth testing, so the null is weaker than the corrected p-values make
it look.

The fix belongs in a new registration rather than here: take the ratio of
the **penalised objectives**, which cannot go negative because the full
model evaluated at `b_c = 0` is exactly the reduced model.

## What this does and does not license

It does not license the full regularized model. P1 predicted that and P1
held; the criterion registered for proceeding -- a cell significant under
both scorings with the same sign -- was not met.

It does not establish that the channels are empty. Sixty-two of 78 cells
never ran, at 2 features and a floor of 10 the events-per-variable is 5,
and the statistic is conservative by the defect above. Three independent
reasons the instrument is blunt, all of them stated before the result
except the third.

The one actionable item is SARS-CoV-2 M4 under measured negatives, and the
action is more labels rather than more modelling.

## Addendum — the verdict's own recommendation does not survive measurement

The verdict above says the one actionable item is SARS-CoV-2 M4 and that
"the action is more labels rather than more modelling". Sized by
`scripts/label_headroom.py` on 2026-10-17, that recommendation is wrong,
and wrong in a way that matters more than the cell it was about.

**Required labelled-set size for 80% power at the corrected alpha, against
the most any panel over this graph could supply** -- the channel's reach,
since a compound no metapath arrives at cannot enter the regression however
well it is screened:

    cell                ceiling   needed   ratio
    SARS-CoV-2 M4         2,982    3,869    1.30x
    SARS-CoV-2 M5         2,192    6,504    2.97x
    SARS-CoV-2 M1         3,021   17,748    5.87x
    SARS-CoV   M7            93    1,038   11.16x
    SARS-CoV-2 M2         1,133   24,962   22.03x
    SARS-CoV-2 M7           232   11,117   47.92x
    SARS-CoV-2 M6           232   11,182   48.20x
    SARS-CoV-2 M3           836   55,025   65.82x

**Zero of eight.** Labelling every compound every channel already reaches
would not bring one cell to 80% power. The binding constraint is REACH, not
labels, and a second screening panel helps only insofar as it brings
compounds the graph can already route to.

M4 is the near miss at 1.30x, short by 887 compounds. That is the one place
where a modest extension of reach -- not of labels -- would change the
answer.

**The scoring choice is worth real power.** The required count is lower
under pdp on four of five cells where both fit, because pdp decorrelates
the channel from degree and the `(1 - rho^2)` term rewards that directly:

    SARS-CoV-2 M4   rho +0.775 -> +0.575   needed  6,534 -> 3,869
    SARS-CoV-2 M2   rho +0.565 -> -0.032   needed  398,194 -> 24,962

A sixteen-fold difference in required labels from the scoring alone. That
is the strongest argument yet for switching the default, and it is still
not a reason to switch it without a registration.

**Caveat, and it is the load-bearing one.** Every `needed` figure is
derived from a coefficient estimated on a cell that was not significant.
An effect size taken from a null result is the weakest possible input to a
power calculation, and the column should be read as an order of magnitude.
What survives that caveat is the ordering and the scale: the shortfalls run
from 1.3x to 66x, and no plausible revision of a coefficient moves a 48x
cell into reach.
