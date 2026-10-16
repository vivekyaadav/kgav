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
