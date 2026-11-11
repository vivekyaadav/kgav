# Pre-registration: the same strength, symmetrically

**SIGNED 2026-11-11, BEFORE THE CODE AND BEFORE THE RUN.** No score under
this protocol has been computed. `scripts/structure_scorer.py` implements
two controls, `publication` and `knn`, and nothing in it reads a per-cell
k; `tests/test_prereg_strength_matched.py` fails if this document still
claims to predate its code once that changes.

**NOT BLIND, AND FOR THE SECOND TIME.** The design responds to the verdict
of `docs/prereg_structure_scorer_symmetric.md`, run earlier today, which
held all six of its predictions and then named this study as the thing that
would decide whether to believe them. Fixing the parameters before the run
is what the signature buys; it does not make this independent.

Date: 2026-11-11. Release: v0.1 (78,789 nodes / 392,375 edges, `nodes_sha`
c800b9005c8b074c, `edges_sha` ad7121bc2195a93a). Prior runs:
`data/results/structure_scorer.json` (publication control),
`data/results/structure_scorer_symmetric.json` (k=10).

## The confound this exists to remove

The symmetric control at k=10 raised ten of fourteen cells to clear the
corrected alpha above chance, none below, against five above and two below
under publication-disjointness. Its own verdict named the reason that might
be an artifact rather than a finding: **the two controls differ in strength
as well as in symmetry.** Publication-disjointness removed between 3.16 and
26.97 references per positive depending on the cell; (3') removed exactly
10. So in eleven of fourteen cells the new control was also *weaker*, and a
weaker control admits more leak.

The prior verdict bounded that confound and could not eliminate it: in three
cells the old control removed fewer than 10 per positive, so (3') was
strictly stronger there, and the AUC rose in all three — by +0.0292, +0.0496
and +0.0600, which are small beside JEV's +0.6550. **So the confound cannot
explain the direction and may explain much of the magnitude.** This study
removes it.

## The measurement

Identical to `docs/prereg_structure_scorer_symmetric.md` in every respect
— ECFP4 radius 2 over 2048 bits, maximum Tanimoto to admissible actives,
the tie-safe rank-sum AUC, the floor/ceiling/skill decomposition,
ChEMBL-only per virus, exclusions (1) self and (2) `inchikey_skel` applied
first — **except that k is per cell**:

    k_cell = ceil(max(pub_per_positive, pub_per_negative))

where both quantities are the prior run's measured exclusion rates. The
control therefore removes, from **every** candidate of **both** classes, at
least as many references as publication-disjointness removed from either
class in that cell. It is symmetric by construction and no weaker than the
control it is compared against.

### Fixed in advance

- **The k values, derived from the prior run and listed here so they cannot
  drift.** `refs` is the reference count after self-exclusion, `left` what
  survives (3'), and `at k=10` the figure this study is testing:

    virus            refs  pub/pos  pub/neg  k_cell   left  at k=10
    SARS-CoV-2      1,785    20.63    32.23      33  1,752   0.8286
    DENV-2            995    26.97     2.12      27    968   0.7882
    HCoV-OC43         163    25.05    26.99      27    136   0.5677
    DENV-3            212    23.66     1.94      24    188   0.6496
    DENV-1            194    13.68     2.75      14    180   0.6066
    Yellow fever      134    12.13     2.99      13    121   0.8388
    MERS-CoV          132    10.36    12.61      13    119   0.4827
    DENV-4            152    11.93     1.98      12    140   0.5868
    JEV                18    11.16     2.00      12      6   0.8002
    TBEV               71    10.31     5.91      11     60   0.7462
    HCoV-229E         119    10.90     3.49      11    108   0.7114
    SARS-CoV          142     8.97     6.41       9    133   0.5785
    Zika              747     6.86     0.84       7    740   0.6812
    West Nile         260     3.16     2.18       4    256   0.6600

  Every cell is feasible; the tightest is JEV, which keeps 6 references.
  **That feasibility is the reason `ceil(max(...))` is registered rather
  than something sterner** — `ceil(pub_per_positive + pub_per_negative)`
  would put JEV at 14 of 18 and is admissible, but leaves 4, and a cell
  reduced to single-figure references is measuring its own sparsity.

- **The correction family is unchanged: the 14 readable ChEMBL cells,
  two-sided, alpha = 0.05/14 = 0.003571, z = 2.9137.** The k values differ
  per cell and the family does not, because the cells are the tests.

- **k = 10 is retained as a registered companion column**, not a second
  family member, so every cell's movement is readable in one table.

- **THE MATCHED NCATS CELL CANNOT BE STRENGTH-MATCHED, AND THAT IS A
  FINDING ABOUT THE PRIOR CONTROL RATHER THAN A GAP HERE.**
  Publication-disjointness removed **88.00 references per positive and
  89.00 per negative there, out of 88** — every admissible reference for
  every candidate, which is why that cell was vacuous in the first study.
  `ceil(max(...))` is 89 and leaves -1. So no symmetric control can match
  that strength, because that strength was total. The cell is registered to
  run at **k = 10 and k = 33**: 10 for continuity with the symmetric study,
  33 because it equals the ChEMBL SARS-CoV-2 cell's `k_cell` and makes the
  source-matched and plate-matched populations comparable at an identical
  handicap. **Neither is strength-matched and the output must say so.**

- **Degree is recomputed and must be unchanged**, as in both prior runs; a
  discrepancy is an implementation fault, not a result.

## Registered predictions

**P1. SARS-CoV-2 at k_cell = 33 reads at or above 0.6346.** The same
threshold as the prior two studies — half the margin its
publication-disjoint figure held above chance — and the same load-bearing
role. At k = 33 the control removes more from both classes than
publication-disjointness removed from either, so **if P1 holds, leniency is
ruled out as the explanation of the symmetric study's result.** If it fails,
that study was measuring its own weakness and its verdict must be withdrawn.

**P2. Every candidate in every cell loses exactly `k_cell` references to
(3'), plus its skeleton group.** The invariant, checked where it can be
observed: `excluded_knn_pos == k_cell * candidates_pos` and likewise for
negatives, with `knn_short` zero. If it fails nothing else is
interpretable.

**P3. AUC at `k_cell` is at or below AUC at k = 10 in at least 10 of the 14
cells.** `k_cell` is >= 10 in eleven cells and < 10 in three (SARS-CoV at 9,
Zika at 7, West Nile at 4), so most cells get a sterner control and should
not improve. The three where it is weaker are named now so they cannot be
read later as evidence of anything.

**P4. No cell clears the corrected alpha BELOW chance.** This is the
sharpest test of the asymmetry mechanism. Strength rises in eleven cells
while remaining symmetric; if anti-ranking returns anyway, then strength and
not asymmetry produced the publication control's below-chance cells, and the
symmetric study's central claim is wrong.

**P5. At least 8 of the 14 cells clear the corrected alpha above chance.**
**This prediction is what discriminates the two explanations.** Under the
asymmetry account the count should stay near the symmetric study's ten, since
symmetry is preserved; under the leniency account it should fall back toward
the publication control's five, since strength is restored. Registering 8
puts the threshold between them.

**P6. The cell ordering is stable: Spearman between the `k_cell` AUC and
the k = 10 AUC exceeds 0.544 over the 14 cells.** If the ranking of viruses
survives a change in the control, it is a property of the data; if it does
not, every per-virus statement in all three verdicts is about the control.

## What would change the conclusion

P1 and P5 together, with P4 holding. That would mean the structure signal
survives a control that is symmetric and at least as strong as the one which
produced the original below-chance readings — at which point the symmetric
study's verdict stands and the collection-confound bound of roughly 0.22 AUC
from its matched cell is the main remaining caveat. P1 holding while P5 falls
to five would mean SARS-CoV-2 is real and most other cells were leniency.

## What this cannot do

- **It cannot make either follow-up blind**, and this is now the second
  registration in a chain responding to its own predecessor's result.
- **It cannot strength-match the plate-matched population**, for the reason
  given above: that control removed everything there. The 0.22 AUC
  collection-confound bound therefore rests on a k = 10 comparison and is
  re-measured here at k = 33 rather than confirmed.
- **Equal in count is still not equal in effect.** A candidate with forty
  close analogs loses its 33 closest and keeps seven; one with 15 loses all
  but none. The claim is the weaker one.
- **It still is not a mechanism**, and everything the two earlier
  registrations say under this heading stands.
- **It does not explain MERS-CoV**, which fell from 0.6677 to 0.4827 between
  the two runs and which no registration has accounted for. Its `k_cell` is
  13 against 41 negatives, the **second-smallest** negative class in the
  family after JEV's 27. The first draft of this line called it the
  smallest, which is wrong and was the convenient reading: JEV is smaller
  and moved the opposite way, from 0.1452 to 0.8002, so small-sample noise
  does not on its own account for MERS-CoV's direction.

## Implementation constraints

No new computation: the per-cell similarity matrix is already built and (3')
is a partial sort of each row, so a per-cell k costs exactly what a fixed k
costs. The k values are read from the prior run's counters rather than
recomputed, which is what makes them auditable; a test re-derives all
fourteen from `data/results/structure_scorer.json`.

**The run must print `k_cell` beside every row**, because P3 and P5 are
unreadable without knowing which cells got a sterner control and which got a
weaker one.

## What is still the owner's to change — but only now

Every value above is fixed and this document is signed. The three most
contestable, with the reasoning already given: `ceil(max(...))` rather than
a sterner combination; keeping the family at 14 cells rather than adding the
plate-matched one; and running the plate-matched cell at k = 10 and 33
rather than declaring it unmeasurable.

---

# Verdict — leniency is ruled out, one prediction failed, and the plate-matched evidence did not survive

Run 2026-11-11, `make structure-scorer-strength`, written to
`data/results/structure_scorer_strength_matched.json`. Release v0.1,
`edges_sha` ad7121bc2195a93a. The registration above is unchanged.

| | prediction | outcome |
|---|---|---|
| **P1** | SARS-CoV-2 at k_cell=33 reads >= 0.6346 | **HELD — 0.7841** |
| **P2** | every candidate loses exactly k_cell | **HELD in all 14 cells and the matched one** |
| **P3** | AUC at k_cell <= AUC at k=10 in >= 10 of 14 | **FAILED — 8 of 14** |
| **P4** | no cell clears the corrected alpha below chance | **HELD — zero** |
| **P5** | at least 8 of 14 clear above chance | **HELD — 10, the symmetric study's count, not the publication control's 5** |
| **P6** | cell ordering stable, Spearman > 0.544 | **HELD — +0.912** |

## The numbers

    virus             pos     neg   k      AUC      corrected 95%   skill    k=10   prior     deg
    Yellow fever      135     882  13  0.8117^ [0.7443, 0.8791]  0.6234  0.8388  0.5545  0.7264
    JEV                19      27  12  0.7885^ [0.5810, 0.9960]  0.5770  0.8002  0.1452  0.6374
    SARS-CoV-2      1,786   5,362  33  0.7841^ [0.7641, 0.8041]  0.5682  0.8286  0.7693  0.5131
    TBEV               72      65  11  0.7385^ [0.6159, 0.8611]  0.4769  0.7462  0.5764  0.4203
    West Nile         261     311   4  0.7079^ [0.6440, 0.7717]  0.4157  0.6600  0.6308  0.5657
    Zika              748   2,577   7  0.7032^ [0.6695, 0.7369]  0.4065  0.6812  0.6316  0.5534
    DENV-2            996   3,043  27  0.6919^ [0.6621, 0.7216]  0.3837  0.7882  0.5735  0.4016
    HCoV-229E         120      88  11  0.6757^ [0.5685, 0.7829]  0.3514  0.7114  0.5951  0.5908
    DENV-3            213      95  24  0.6653^ [0.5728, 0.7579]  0.3307  0.6496  0.3861  0.1715
    DENV-1            195     123  14  0.5967^ [0.5034, 0.6901]  0.1935  0.6066  0.4489  0.3144
    HCoV-OC43         164      70  27  0.5858  [0.4704, 0.7012]  0.1716  0.5677  0.5396  0.5345
    SARS-CoV          143      96   9  0.5825  [0.4744, 0.6907]  0.1651  0.5785  0.5186  0.6812
    MERS-CoV          133      41  13  0.4893  [0.3381, 0.6405] -0.0215  0.4827  0.6677  0.7085
    DENV-4            153     109  12  0.4792  [0.3735, 0.5850] -0.0416  0.5868  0.4644  0.2302

    ^ clears the corrected alpha above chance.  No cell clears it below.

## P1 and P5 together: the symmetric study was not measuring its own leniency

This was the question the registration existed to answer. At k_cell = 33 the
control removes 33 references from **every** candidate of **both** classes,
against publication-disjointness removing 20.63 per positive and 32.23 per
negative in that cell — so it is strictly sterner on both sides and
symmetric. SARS-CoV-2 reads **0.7841 [0.7641, 0.8041]**, above the 0.6346
threshold and above the publication-disjoint 0.7693 it is compared with.

**And P5 is the discriminating result.** Ten of fourteen cells clear above
chance — the symmetric study's count, not the publication control's five.
The leniency account predicted a fall toward five when strength was
restored. It did not happen. **Combined with P4, zero cells clearing below
chance, the asymmetry explanation stands and the symmetric study's verdict
is not withdrawn.**

P6 adds that the cell ordering is almost unchanged: Spearman **+0.912**
between the k_cell and k=10 readings. The per-virus pattern is a property of
the data, not of the control.

## P3 failed, and the registered reasoning behind it was wrong

P3 predicted AUC would not rise under a sterner control in at least 10 of 14
cells. It rose in six, so only 8 of 14 fell:

    West Nile      k_cell= 4   +0.0479   k < 10, named in advance
    Zika           k_cell= 7   +0.0221   k < 10, named in advance
    SARS-CoV       k_cell= 9   +0.0040   k < 10, named in advance
    HCoV-OC43      k_cell=27   +0.0181   k > 10, NOT anticipated
    DENV-3         k_cell=24   +0.0157   k > 10, NOT anticipated
    MERS-CoV       k_cell=13   +0.0066   k > 10, NOT anticipated

Three were named before the run, because their k_cell is below 10 and they
therefore get a *weaker* control than the companion column. **The other
three are the failure, and the registration's reasoning for P3 was simply
unsound.** It argued that "removing more of a candidate's nearest references
cannot add information". That is true of a single candidate's score and false
of an AUC, which compares two classes: removing more references lowers
scores in *both*, and whether the comparison improves depends on which class
loses more. Monotonicity in k was never implied by the design.

The three unanticipated rises are +0.0181, +0.0157 and +0.0066 — small, and
in cells with 70, 95 and 41 negatives. That is consistent with noise and is
not evidence of anything; it is stated as the reason P3's threshold should
never have been set from the symmetric study's incidental 10-of-14, which is
what calibrating a prediction on a prior result's margin amounts to.

## The result that weakens the previous verdict

**The plate-matched NCATS cell no longer clears the corrected alpha.** At
k=33 it reads **0.5898 [0.4968, 0.6827]** — the interval contains chance. At
k=10 it read 0.6088 and cleared, which is what the symmetric verdict called
"the first result in this project where a channel clears a corrected
threshold on a population matched by construction".

**That claim does not survive the stronger handicap, and it is the headline
of the previous verdict.** What survives is weaker and still worth having:
on the plate-matched population the scorer reads above chance at both k
values, with skill 0.1795 at k=33, and the interval excludes chance at k=10
but not at k=33. The honest summary is that the plate-matched evidence is
**marginal and k-dependent**, not established.

The collection-confound gap narrows slightly when both cells are read at the
same handicap: **0.1944 at k=33**, against 0.22 measured at k=10. Both are
large.

## What this licenses now

- **The structure scorer separates antiviral actives from measured inactives
  on the source-matched population, under a symmetric control at least as
  strong as the one that produced the original below-chance cells.** Ten of
  fourteen cells, zero below chance, ordering stable at +0.912. This is the
  firmest result in the project.
- **It does not license the plate-matched claim.** That cell is marginal,
  and the previous verdict's headline about it should be read as retracted to
  "above chance at both k, clearing at one of them".
- **The collection confound remains large and measured**, 0.19 to 0.22 AUC
  depending on the handicap.
- **MERS-CoV is still unexplained**, now across three runs: 0.6677
  publication, 0.4827 at k=10, 0.4893 at k_cell=13. It is the only cell that
  the publication control favoured, and no registration has accounted for
  it.
- **It still is not a mechanism.**

## What would change the conclusion

The plate-matched cell is the whole question now, and it is powered at 89
positives against 7,343 negatives — its interval is +-0.09 at the corrected
level, which is why it flips between clearing and not. **More
plate-matched positives is the only thing that settles it**, and that is a
data-acquisition problem rather than an analysis one: curating further
published single-screen panels, which `docs/v1_release_plan.md` already
lists as out of scope and large. No further control variant will help,
because the limit is the class size and not the protocol.
