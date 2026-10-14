# Pre-registration: letting potency into the direct-acting channel

Written **before** the code, and committed before the measurement, because a
prediction made afterwards is a description. If the result contradicts what
is written here, this file stays as it is and the verdict goes beneath it.

Date: 2026-10-12. Release: v0.1 (74,872 nodes / 374,686 edges, schema
v0.12.0). Virus: NCBITaxon:2697049.

## The defect

**M1's DWPC is a property of the target, not of the drug.**

`walk_metapath` starts its frontier at weight 1.0 and multiplies by each
neighbour's `degree ** -w`. Nothing drug-specific ever enters. For
`M1 = [SmallMolecule, INHIBITS, Protein, ENCODED_BY, Gene, BELONGS_TO,
OrganismTaxon]` the score is therefore a function of which viral proteins the
drug reaches and of nothing else. Two drugs hitting the same protein are
identical by construction, whatever their potency.

Measured on v0.1, per channel, by `scripts/diagnose_controls.py`:

    chan     drugs  distinct  largest tie  % of pool  % in any tie
    M1       3,020        12        2,441      80.8%         99.9%
    M2       1,133       348           58       5.1%         78.6%
    M3         836       345           81       9.7%         69.1%
    M4       2,982     1,835           28       0.9%         50.5%
    M5       2,192     1,202           37       1.7%         56.0%
    M6         232       133            7       3.0%         63.8%
    M7         232        97           14       6.0%         76.7%

**Twelve distinct values across 3,020 drugs.** The degeneracy is specific to
M1: M4 and M5 resolve 1,835 and 1,202 distinct values over comparable pools.
M1 is not ranking, it is grouping by target, and 80.8% of the pool sits in
one group — every drug whose only route to this virus is nsp5
(`UniProtKB:PRO_0000449623`), all of them scoring 0.000242.

nirmatrelvir and chloroquine are both in that group, with byte-identical raw
DWPC, through the same node. The channel cannot distinguish the approved
antiviral from the cell-line artifact.

And the degeneracy reaches the metrics. `rank()` has no tiebreak and
`hits_at_k` reads `ranked[:k]`, so where a block straddles k the figure is
decided by dict insertion order:

    ranking       k    Hits@k range  tied at cut  of which pos
    COMBINED     10           2 - 2            4             0
    COMBINED     50           4 - 4            0             0
    COMBINED    100           7 - 7            4             0
    M1           10           1 - 2            4             1   NOT DETERMINED
    M1           50          2 - 12           98            10   NOT DETERMINED
    M1          100          2 - 12           98            10   NOT DETERMINED

M1's Hits@50 and Hits@100 are each one draw from [2, 12]. COMBINED is
determined at all three cuts, so this is a property of M1 rather than of the
evaluation.

## The proposed change

Let the measured potency on an `INHIBITS` edge modulate that edge's
contribution to the DWPC weight, so that two drugs reaching the same protein
are separated by what was measured about them.

Nothing is designed here. The design follows the precondition below.

## Precondition, to be measured first

**How many M1-traversable `INHIBITS` edges carry an exact potency at all.**

If most do not, potency weighting mostly imputes, and an imputed value that
decides a ranking is worse than the tie it replaces. The measurement must
report, over edges `M1` actually walks for this virus:

- count with an exact `ic50_nm`, `ki_nm` or `kd_nm` (`relation` of `=`)
- count with a censored relation (`>`, `<`, `>=`, `<=`, `~`)
- count with none of them
- the same three broken down for the 2,441-drug nsp5 block specifically

**If fewer than half of the nsp5 block's edges carry an exact value, this
change is abandoned**, and the honest statement is that the graph records
which protein a compound hits and not how well.

Censored values are not coerced. `chembl.py` already refuses that —
"IC50 > 10000 nM is a NEGATIVE result ... coercing it to a potency number
inverts its meaning" — and a censored edge must either keep the tie or leave
the channel, never be read as a number.

## Predictions

Stated so they can fail.

**P1 — tie structure (near-certain; a precondition check, not evidence).**
M1 distinct values rise from 12 to over 1,000, and the largest tie block
falls from 80.8% of the pool to under 5%. If this fails, the change did not
take effect.

**P2 — the decisive directional prediction.** nirmatrelvir and chloroquine
separate in M1: nirmatrelvir's M1 percentile above 0.90, chloroquine's below
0.50. They currently share 0.4041. Nirmatrelvir's nsp5 IC50 is low
nanomolar; chloroquine's nsp5 record is not. **If they do not separate, the
change has failed**, whatever happens to P1 — a less degenerate ranking that
still cannot tell those two apart has bought nothing.

**P3 — the confound that would invalidate P2.** Potency weighting must not
be a popularity term in disguise. Well-studied compounds have more
measurements, more of them exact, and better-characterised potency. If M1's
improvement tracks the degree baseline, we have re-derived degree. The
degree baseline is already the only scorer on this release whose interval
clears chance (0.582 [0.520, 0.645]) and it is the null hypothesis.
**Required: the rank correlation between the new M1 score and node degree
must not exceed that between the old M1 score and node degree by more than
0.1.** If it does, the change is reverted regardless of P1 and P2.

**P4 — a selection effect to state, not to test.** A compound with a strong
measured IC50 against a viral protein is more likely to have been carried
forward into the cell-based assay that produces the label. So potency
correlates with the *existence* of a label independently of whether the
compound works. `publication_disjoint` removes the shared-publication case
and does not remove this one. Any improvement in a label-based metric is
therefore an upper bound on the real improvement, and this paragraph exists
so that cannot be forgotten when the numbers look better.

## What this change CANNOT be validated by, on this release

This is the part most likely to be misread later, so it is stated before the
result exists.

**Neither protocol can evaluate M1.** `docs/evaluation.md` records that M1
reaches **3 of 89** measured actives under the hard-negative protocol,
against `MIN_REACHED_POSITIVES = 10`. Potency weighting reorders the drugs
M1 reaches; it does not change *which* drugs it reaches. So `n_reached_pos`
stays at 3 and **M1 remains NOT EVALUABLE after this change.** The temporal
protocol reaches 8 post-cutoff positives, also below 10.

Therefore:

- **No AUC claim may be made for M1 after this change**, improved or not.
  Its interval will still describe the imputed floor.
- The headline figure — SARS-CoV-2 M1 temporal lift 2.22, p=0.0196 — rests
  on **8 positives and 6 hits**, and `hits_at_k_range` must be run on that
  split before and after. If that cut lands in a tie block, the 2.22 was
  never determined either, and that is a finding about the existing result
  rather than about this change.
- A calibrated ranking will still return `{}`. This change cannot make M1
  earn a place in a ranking on v0.1, and anything that appears to is a bug.

**What can be concluded** is confined to: the tie structure (P1), the
control separation (P2), the degree confound (P3), and whether the headline
lift was determined. Nothing about discrimination.

## Decision rule

Abandon if the precondition fails. Revert if P3 fails. Keep only if P1 and
P2 both hold, and in that case claim the control separation and the tie
structure, and claim nothing about AUC or about a shortlist.

---

# Amendment 1 — the design, and a corrected P2

2026-10-12, after the precondition and **before any implementation or any
measurement of P1, P2 or P3.** The original text above is unedited; this is
appended, dated and reasoned, because a pre-registration that gets quietly
rewritten is not one.

## The precondition was met

    LARGEST TIE BLOCK (2,441 drugs)
                     edges   % edges     drugs   % drugs
      exact          2,084     85.4%     2,084     85.4%
      censored         349     14.3%       349     14.3%
      none               8      0.3%         8      0.3%

85.4% against the 50% the text above required. The controls carry usable
numbers on the same node, `UniProtKB:PRO_0000449623`:

    nirmatrelvir   ic50_nm=75    ki_nm=3   kd_nm=190
    chloroquine    ic50_nm=7280
    remdesivir     ic50_nm=1640  (nsp12, a different group)

A 97-fold potency ratio between the approved antiviral and the artifact,
both currently scoring 0.4041.

## P2 AS WRITTEN IS UNREACHABLE, and the reason is arithmetic

P2 required nirmatrelvir's M1 percentile above 0.90. To pass 90% of 3,020
drugs it must leave its own 2,441-drug block AND overtake roughly 277 drugs
in higher-scoring groups. The nearest of those is the nsp12 group, whose
DWPC is **3.15x** nirmatrelvir's on the degree term alone — confirmed
against the release, 0.000762 / 0.000242 = 3.149.

Potency in log space does not span that. With the standard pX = 9 -
log10(nM), nirmatrelvir is 7.125 and remdesivir 5.785: a factor of 1.232.
Closing 3.15x from 1.232x requires **pX ** 5.5**, an exponent chosen to
produce the wanted answer — a fitted parameter in a baseline whose entire
claim is that it has none.

So P2's threshold was mis-specified, and it was mis-specified in a
particular direction: it demanded that potency overrule the degree
correction ACROSS targets, which is a different proposition from the defect
being fixed and one this change does not argue for. Nothing in the
measurement caused this correction; the arithmetic above was available when
the original was written, and was not done.

## The design, chosen before the numbers exist

**Potency breaks exact ties within M1 and does nothing else.** The M1 score
remains a function of the target; a second key orders drugs inside a block
of identical raw DWPC.

- **Nothing outside a tie block moves.** Block boundaries come from the raw
  DWPC sort, which is untouched, so the block occupies the same span of
  percentiles and every cross-target ordering — remdesivir above the nsp5
  group included — is preserved exactly. COMBINED does not change.
- **Path specificity and measured affinity stay separate quantities.** This
  repo refuses such collapses elsewhere: `domain` is recorded rather than
  folded into nsp3, and `unquantified`'s three coexisting definitions are
  documented as a defect. Multiplying a degree penalty by a binding constant
  would make one number out of two different claims.
- **The quantity is chosen by a FIXED PRECEDENCE, never by value:**
  `ic50_nm`, then `ki_nm`, then `kd_nm`. `selectivity.py` records what
  choosing by value costs — "Taking the worst ratio defamed nirmatrelvir (SI
  1.8) ... Taking the best exonerated chloroquine (17.7). Both were the same
  error -- selecting a study by its answer." Where several exact values of
  the selected field exist, the WORST (least potent) is kept, which is that
  module's own rule.
- **Censored bounds are not coverage.** A relation of >, <, >=, <= or ~ is a
  bound, and on this project's rule a bound in the wrong direction is a
  measurement of inactivity. Censored and absent edges are not ordered.
- **An unmeasured drug keeps the score it has today**: the block's midrank.
  That is its current value, so the change is a strict refinement — no
  imputation decides anything, and "unmeasured" is never read as "weak".
- **M1 only.** M7 borrows a neighbour's target, so its drug's own potency is
  not what that path asserts, and its largest tie is 14 drugs (6% of pool).
  M2–M6 are host channels where potency against a viral protein is
  irrelevant.

## P2, corrected

**P2' — the controls separate within the block.** nirmatrelvir lands in the
**top decile of the tie block** by potency rank and chloroquine in the
**bottom half**, and their M1 percentiles differ by **at least 0.30**. They
currently differ by 0.0000.

The ceiling is stated so it cannot be mistaken for a shortfall later: the
block spans roughly percentile 0.002 to 0.806, so **nirmatrelvir cannot
exceed ~0.81 by construction** and no result above that should be expected
or believed.

If 75 nM against 7,280 nM does not place those two in those regions, the
potency data is internally inconsistent and the change has failed.

## P1 and P3 are unchanged, with one clarification

P1's "distinct values past 1,000" now refers to M1's distinct values after
tie-breaking, which is the same quantity. P3's degree-correlation bound of
0.1 is unchanged and is the condition most likely to bite: potency exists
preferentially for well-studied compounds, and 14.3% of the block is
censored while 0.3% has nothing — both populations skew.

## A new caveat, stated with P4

**P5 — mixed quantity types are not strictly comparable.** A drug ordered on
its Ki is being compared against one ordered on its IC50, and an IC50 is
assay-dependent in a way a Ki is not. The precedence makes the choice
reproducible, not commensurable. The implementation must therefore report
how many drugs were ordered on each field, so the size of this problem is
visible rather than assumed small.

## Decision rule, restated

Keep only if P1 and P2' both hold and P3 does not trip. Claim the control
separation and the tie structure. Claim nothing about AUC, nothing about a
shortlist, and nothing about M1 having become evaluable — it reaches 3 of 89
measured actives before this change and will reach 3 after.

---

# Verdict — NOT ADOPTED

2026-10-12, measured by `scripts/diagnose_m1_potency.py` on v0.1,
NCBITaxon:2697049. The decision rule above requires P1 and P2' both to hold.
Both fail, so **the tie-break is not adopted**: it stays implemented and
off, `combine` is called without it, and no published figure moves.

    P1 FAILS | P2' FAILS | P3 HOLDS

    distinct M1 scores        12  ->   1,385      predicted past 1,000   MET
    largest tie block      80.8%  ->    11.8%     predicted under 5%     NOT MET

    nirmatrelvir   pct 0.4041 -> 0.5855   potency rank 72.4th pctile of block
    chloroquine    pct 0.4041 -> 0.0295   potency rank  3.6th pctile of block
    gap +0.5560

    rank corr with degree  -0.0048 -> +0.0785   delta +0.0833  (bound 0.10)

    ordered on ic50_nm  2,489  93.5%
    ordered on ki_nm      172   6.5%
    ordered on kd_nm        1   0.0%

## P1 failed on a threshold that was impossible, and that is my error

The residual largest block is 356 drugs. The block holds 349 censored edges
and 8 with no potency field: 357 drugs, 11.8% of the pool. **The residual
block IS the unkeyed population.** By the design fixed in Amendment 1 an
unmeasured drug keeps the block's midrank, so those 357 must stay tied and
the largest block cannot fall below 11.8% — while P1 demanded under 5%.

P1 was written before the design existed. Amendment 1 then carried it
forward with the single remark that it "now refers to M1's distinct values
after tie-breaking, which is the same quantity", and did not check its second
clause against the design being adopted in the same document. That is the
same failure as the original P2: a threshold asserted without the arithmetic.
Twice, in two consecutive amendments.

The mechanism did what it was built to do — 12 distinct values became 1,385,
and every drug that remains tied is one the graph has no exact measurement
for. That is not a pass. P1 as written is failed, and a threshold corrected
after seeing the number it scores is not a threshold.

## P2' failed on one of its three criteria

    gap >= 0.30                     +0.5560         PASS
    chloroquine in the bottom half  3.6th pctile    PASS
    nirmatrelvir in the top decile  72.4th pctile   FAIL

The separation is large and in the right direction: chloroquine moves from
0.4041 to 0.0295, near the bottom of the nsp5 block, and the gap is 0.5560
against a required 0.30. What fails is the claim that nirmatrelvir should sit
in the top decile of nsp5 inhibitors by potency.

**It does not, and it probably should not.** Its chosen field is ic50_nm = 75
nM, and roughly 580 of the 2,084 exactly-measured nsp5 binders in ChEMBL are
more potent than that. Nirmatrelvir did not reach the clinic by being the
most potent Mpro inhibitor ever measured; it reached it on oral
bioavailability and a ritonavir boost that compensates for rapid CYP3A4
metabolism. The top-decile criterion encoded an assumption — that clinical
success implies top-decile biochemical potency — which this measurement
gives no reason to believe and some reason to doubt.

That is an explanation, not a pass. The prediction was mine, it was wrong
about the world, and the rule it feeds says do not adopt.

## The precedence question, raised with its timing declared

nirmatrelvir carries `ic50_nm=75`, `ki_nm=3` and `kd_nm=190`. The fixed
precedence reads ic50 and lands at the 72.4th percentile. Reading Ki would
read 3 nM and rank it far higher, and there is a real pharmacological
argument for Ki: for a competitive inhibitor it is a concentration-
independent constant, while an IC50 depends on substrate concentration and
assay conditions and is therefore less comparable across sources.

**That argument is being raised after seeing that it would flip the
verdict**, which is exactly the position Amendment 1 quoted selectivity.py
against — "selecting a study by its answer". It is recorded here and not
acted on. P5 also shows the question is narrow: 93.5% of keyed drugs are
ordered on ic50 and 6.5% on ki, so the precedence changes almost nothing
except for the one compound whose number has now been seen.

If the precedence is to change it must be decided on pharmacological grounds
and then tested against evidence this run has not already shown.

## P3 held, with 83% of its budget spent

Rank correlation with node degree moved from -0.0048 to +0.0785, a delta of
+0.0833 against a bound of 0.10. It held, and it nearly did not. Potency
exists preferentially for well-studied compounds and this is that effect,
measured: four fifths of the allowance, from a change that only reorders
within blocks. Any future design that reorders ACROSS blocks should expect
to breach it.

## What would be evidence, if this is pursued

Not a re-score of this run. The controls in the block are n=2, and both
their numbers have now been seen, so nothing further can be concluded from
them.

The available independent test uses the held-out labels: within the
2,441-drug block, do the measured ACTIVES sit higher by potency rank than
the measured INACTIVES? There are 2,020 compounds measured active against
this virus and 10,956 measured inactive, so the n is large, the comparison
does not depend on anyone's beliefs about two compounds, and it asks the
question that matters — whether potency rank carries information about
antiviral activity at all, inside the population the tie-break reorders.

That test must be pre-registered before it is run, in a new section, with
its own threshold and the arithmetic done first.

---

# Addendum — the design cited a stale comment

2026-10-14, after the verdict above. Recorded because the design was
pre-registered and the ground under part of it moved.

Amendment 1 justified keeping the WORST of several exact values like this:

> The quantity is chosen by a FIXED PRECEDENCE, never by value.
> `selectivity.py` records what choosing by value costs [...] Where several
> exact values of the selected field exist, the WORST is kept, which is that
> module's own rule.

**It is not that module's rule.** `selectivity.py` takes the MEDIAN of its
exact measurements, and has since a fix recorded in its own emit loop: the
worst of 34 remdesivir measurements is an outlier that reported a licensed
antiviral as cytotoxic at SI 0.9. `tests/test_selectivity.py` pins the
median. What I quoted was a selection comment that had been superseded in
the same function and never corrected, with "collect all; median chosen
below" on the line directly beneath it.

The fixed-precedence half of the rule stands, and it is the half that
mattered: choosing the FIELD by value is the error that module documents,
and the precedence avoids it.

The aggregation half is now stated as this project's own choice rather than
as an inherited one. The situations differ: assembly merges edges on
(subject, predicate, object), so repeated measurements against one protein
are already one edge, and what `potency_tiebreak` aggregates is a drug's
DISTINCT viral targets. Taking the weakest is conservative; it is not
selectivity.py's median rule and should not have claimed to be.

Nothing is re-measured, because the verdict is unchanged: fix 3 was not
adopted, and `combine` is still called without a tiebreak. The stale comment
is corrected at its source and a test refuses its return.
