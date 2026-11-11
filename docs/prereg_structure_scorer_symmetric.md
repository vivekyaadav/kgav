# Pre-registration: a control that handicaps both classes equally

**SIGNED 2026-11-11, BEFORE THE CODE AND BEFORE THE RUN.** No score under
this protocol has been computed. `scripts/structure_scorer.py` exists and
will be extended; nothing in it currently implements a k-nearest exclusion,
and `tests/test_prereg_symmetric_control.py` fails if this document claims
otherwise once it does.

**AND THIS IS NOT A BLIND REGISTRATION, WHICH IS A DIFFERENT THING FROM A
RETROSPECTIVE ONE.** The design responds to a result already in hand:
`docs/prereg_structure_scorer.md`'s verdict, run 2026-11-11, which found
that the publication-disjoint control it registered is asymmetric between
the classes and that the direction of that asymmetry tracks the sign of the
result. A follow-up motivated by a prior finding is ordinary science; what
would not be ordinary is choosing this protocol's parameters after seeing
what they produce, and that is what the signature forbids. The distinction
is stated here so no reader has to infer it.

Date: 2026-11-11. Release: v0.1 (78,500 nodes / 391,727 edges, `nodes_sha`
b0d6fc69d63942c2, `edges_sha` 545fb7d87b0aec08). Prior run:
`data/results/structure_scorer.json`.

## The defect this registers a fix for

The prior registration excluded, from each candidate's reference set, every
known active citing a publication the candidate also cites. Measured on the
run, that exclusion is **not symmetric between the classes**: a positive's
nearest known actives are its own paper-mates and are removed, while a
negative sitting in a paper containing no active loses nothing. References
excluded per candidate, positives against negatives:

    DENV-2          26.97 /   2.12  =  12.7x    AUC 0.5735   <- counterexample
    DENV-3          23.66 /   1.94  =  12.2x    AUC 0.3861
    Zika             6.86 /   0.84  =   8.2x    AUC 0.6316   <- counterexample
    DENV-4          11.93 /   1.98  =   6.0x    AUC 0.4644
    JEV             11.16 /   2.00  =   5.6x    AUC 0.1452
    DENV-1          13.68 /   2.75  =   5.0x    AUC 0.4489
    Yellow fever    12.13 /   2.99  =   4.1x    AUC 0.5545
    HCoV-229E       10.90 /   3.49  =   3.1x    AUC 0.5951
    TBEV            10.31 /   5.91  =   1.7x    AUC 0.5764
    West Nile        3.16 /   2.18  =   1.5x    AUC 0.6308
    SARS-CoV         8.97 /   6.41  =   1.4x    AUC 0.5186
    HCoV-OC43       25.05 /  26.99  =   0.9x    AUC 0.5396
    MERS-CoV        10.36 /  12.61  =   0.8x    AUC 0.6677
    SARS-CoV-2      20.63 /  32.23  =   0.6x    AUC 0.7693

**ALL FOURTEEN CELLS ARE SHOWN, AND THE FIRST DRAFT OF THIS TABLE SHOWED
SIX.** It carried the four highest ratios and the two lowest, which are
exactly the rows that support the hypothesis, and omitted DENV-2 at 12.7x
and Zika at 8.2x — both heavily asymmetric and both ABOVE chance. Selecting
table rows by whether they fit is the failure this project is built against,
and it was caught by a test requiring every readable cell's ratio to appear.
The two counterexamples are marked.

So the pattern is real but partial: the four cells where the control takes
5x or more from positives include the four lowest AUCs **and** the two
counterexamples, and the two cells where it takes more from negatives are
the two highest AUCs. Spearman between AUC and log ratio is **-0.490**,
which does NOT clear the ~0.544 critical value at n=14 — so the prior
verdict reported it as a hypothesis and declined to claim the four
below-chance cells as a finding about chemistry. **This study tests that
hypothesis, and the two counterexamples are the reason it needs testing
rather than asserting.**

**A second reason to replace the control, independent of the first.** PMID
presence is itself correlated with the label, by amounts that differ by
virus and change sign: 184 of 213 DENV-3 positives carry one against 18 of
95 negatives (86% against 19%), while SARS-CoV-2 reads 80% against 42%. A
control keyed to publication metadata therefore varies in strength with a
quantity that is not independent of the outcome. A control keyed to a count
does not.

## The measurement

Identical to the prior registration in every respect — ECFP4 at radius 2
over 2048 bits from `kgav.similarity`, maximum Tanimoto to the admissible
known actives, the tie-safe rank-sum AUC in
`labels.evaluate_against_negatives`, the `auc_floor`/`auc_ceiling`/
`auc_skill` decomposition, ChEMBL-only per virus, `MIN_REACHED_POSITIVES`
of 10 — **except for the third exclusion**, which becomes:

    3'. the candidate's k NEAREST remaining references, by Tanimoto,
        for a single k fixed below and applied identically to every
        candidate in every cell

Exclusions (1) the candidate itself and (2) every reference sharing its
`inchikey_skel` are unchanged and are applied first, so (3') removes exactly
k beyond them.

### Fixed in advance

- **k = 10 is the primary, and the arithmetic is the reason.** The smallest
  reference set in the family is JEV's: 19 positives, 18 after
  self-exclusion. k must leave a non-empty set in **every** cell or the
  study repeats the prior run's P5 failure, where a cell printed 0.5000
  because nothing had been ranked. k=10 leaves 8 admissible references in
  JEV and more everywhere else. **k=30 WAS IN THE DRAFT OF THIS DOCUMENT
  AND IS REMOVED:** it leaves JEV at -12, which empties the set and makes
  the cell vacuous. The prior run needed a run to discover that failure
  mode; this one needed only the class sizes, which is the whole argument
  for stating design arithmetic before signing.

- **The registered sensitivity curve is k ∈ {1, 3, 10}, reported for every
  cell, with k=10 carrying the headline and the correction family.** The
  curve is the estimand that actually answers the question — whether the
  separation survives removing the local neighbourhood — and a single k
  cannot show decay. k=15 is admissible by the arithmetic (3 left in JEV)
  and is excluded as too close to emptying the smallest cell to be read.

- **The publication exclusion is DROPPED from the primary, not added to.**
  Keeping both would reintroduce the asymmetry this study exists to remove.
  It is retained as a **registered secondary**, "both controls", so the two
  runs remain comparable; it is not in the correction family.

- **Correction family: the 14 readable ChEMBL cells at k=10, two-sided.**
  alpha = 0.05/14 = 0.003571, z = 2.9137. Unchanged from the prior
  registration, and for the same reason: a below-chance AUC is a finding
  here. The k=1 and k=3 curve points and the "both controls" secondary are
  reported at the same corrected level and are not added to the
  denominator.

- **The matched NCATS cell is the secondary population and is now
  EXPECTED TO BE MEASURABLE.** It has 89 positives, 88 references after
  self-exclusion, so k=10 leaves 78. The count-based control cannot empty a
  reference set the way publication-disjointness emptied all 7,432
  candidates there. This is the single largest gain of the redesign: the
  only cell that can test the collection confound becomes readable.

- **Degree is recomputed on the identical population and is unchanged by
  the control**, since it does not use a reference set. Its figures will
  therefore equal the prior run's, and a discrepancy is an implementation
  error rather than a result.

- **The residual asymmetry is bounded and will be reported, not assumed
  away.** Exclusion (2) still removes a variable number: 0% to 7% of
  positives and 0% to 14% of negatives sit in a shared-`inchikey_skel`
  group, with no group larger than 3. So total references removed is k plus
  at most 2 for a minority of candidates, against the 5x-to-12x imbalance
  being fixed. The run records `excluded_skeleton_pos` and
  `excluded_skeleton_neg` so the residual is a number rather than a claim.

## Registered predictions

**P1. SARS-CoV-2 at k=10 reads at or above 0.6346.** Derived, not chosen:
it read 0.7693 under the prior control, a margin of 0.2693 above chance, and
0.6346 is the point at which half that margin survives. Above it, the
separation is substantially not co-publication; below it, most of what the
prior run reported was the neighbourhood. **This is the prediction the study
exists to test and it is the one I expect to be closest.**

**P2. Every candidate in every cell loses exactly k references to (3'),
plus its skeleton group.** An invariant, not an outcome: the counters must
read `excluded_knn == k * candidates` in all 14 cells. If they do not, the
control is not symmetric and nothing below means anything. This replaces the
prior registration's P3 as the implementation check.

**P3. The four below-chance cells move toward 0.5, and at least three of
them end above their prior values.** DENV-1 would need +0.0511, DENV-4
+0.0356, DENV-3 +0.1139, JEV +0.3548 to reach chance. **This is the direct
test of the asymmetry hypothesis.** If they stay where they are, the
asymmetry was not the cause, and a real and unexplained anti-ranking in the
dengue serotypes becomes the finding — which would be more interesting than
the artifact.

**P4. AUC is non-increasing in k in at least 10 of the 14 cells**
(k=1 >= k=3 >= k=10, allowing ties). Removing more of a candidate's nearest
references cannot add information. Cells that violate it are reported
individually, because a non-monotone cell is either noise at small n or an
implementation fault, and the class sizes say which is plausible.

**P5. The matched NCATS cell is measurable, and reads below the ChEMBL
SARS-CoV-2 cell at the same k.** The first half is arithmetic and should
not fail; the second is the load-bearing claim the prior run could not
reach. If the matched cell reads as high as the unmatched one, the
collection confound is not driving the ChEMBL figure and the result is
stronger than this registration expects.

**P6. The leak correlation weakens.** Spearman between the k=10 AUC and the
prior run's publication-leak magnitude will be closer to zero than the
-0.684 measured within the prior run. The leak was measured against a
control that caused part of it; a symmetric control should decouple them. A
correlation that stays near -0.684 means co-publication structure reaches
the result through something other than the exclusion rule.

## What would change the conclusion

P1 and P5 together. SARS-CoV-2 above 0.6346 at k=10 **and** a matched NCATS
cell that is measurable and lower would mean the scorer separates the
classes on structure, with the collection confound bounded by a cell that
can see it. Either one alone is not enough: P1 without P5 leaves the
confound untested, and P5 without P1 says only that both cells are weak.

## What this cannot do

- **It cannot make the follow-up blind.** The design was chosen after seeing
  the first result. Fixing the parameters in advance prevents the specific
  failure of choosing them to suit an outcome; it does not turn this into an
  independent replication, and nothing here claims otherwise.
- **Equal in count is not equal in effect.** Every candidate loses exactly
  k references, but a candidate with twenty close analogs loses its ten
  closest and keeps ten, while one with two loses both. This is a strictly
  weaker claim than "equal handicap" and the registration makes only the
  weaker one.
- **It still is not a mechanism.** Everything the prior registration said
  under this heading stands: the scorer names a neighbour, not a target.
- **It does not test the temporal question**, which remains unregistered and
  must be written before it runs.
- **It cannot recover the prior run's P5.** That cell was vacuous and no
  reanalysis makes it informative; it is re-measured here, not repaired
  there.

## Implementation constraints

The full candidate-by-reference similarity matrix is already computed per
cell, so (3') is a partial sort of each candidate's row and costs nothing
new: the largest cell is 7,148 candidates against 1,786 references, 12.8M
similarities, about 51 MB as float32. The three curve points reuse one
matrix per cell. 242 GB and 32 cores available; the per-virus work is
independent.

**The run must print `excluded_knn` per class and per cell**, because P2 is
an invariant over exactly that counter, and the prior run's per-class
counters had to be added after the fact to test a hypothesis its own output
could not see.

## What is still the owner's to change — but only now

Every value above is fixed and this document is signed. Changing any of them
is legitimate until the first score under this protocol exists and
indefensible afterwards. The three most contestable, with the reasoning
already given: k=10 rather than another admissible value; dropping the
publication exclusion from the primary rather than keeping both; and
k ∈ {1, 3, 10} rather than a finer curve.

---

# Verdict — all six predictions held, and the control I was defending was destroying signal

Run 2026-11-11, `make structure-scorer-symmetric`, written to
`data/results/structure_scorer_symmetric.json`. Release v0.1, `edges_sha`
545fb7d87b0aec08. The registration above is unchanged; this is appended
beneath it.

| | prediction | outcome |
|---|---|---|
| **P1** | SARS-CoV-2 at k=10 reads >= 0.6346 | **HELD — 0.8286, and it rose rather than fell** |
| **P2** | every candidate loses exactly k | **HELD in all 14 cells, 10.00/10.00, zero short** |
| **P3** | the four below-chance cells move toward 0.5, >=3 above prior | **HELD — all four, and all four end ABOVE chance** |
| **P4** | AUC non-increasing in k in >=10 of 14 | **HELD at exactly 10 of 14** |
| **P5** | the matched cell is measurable and reads below the ChEMBL cell | **HELD, both halves** |
| **P6** | the leak correlation weakens from -0.684 | **HELD — it is +0.367** |

## The numbers

    virus             pos     neg     k=1     k=3     k=10      corrected 95%   skill    both   prior     deg
    Yellow fever      135     882  0.8282  0.8231  0.8388^ [0.7751, 0.9026]  0.6777  0.5433  0.5545  0.7264
    SARS-CoV-2      1,786   5,362  0.8677  0.8577  0.8286^ [0.8102, 0.8470]  0.6572  0.7570  0.7693  0.5131
    JEV                19      27  0.8090  0.8021  0.8002^ [0.5973, 1.0000]  0.6004  0.1608  0.1452  0.6374
    DENV-2            996   3,043  0.8474  0.8226  0.7882^ [0.7616, 0.8149]  0.5765  0.5828  0.5735  0.4016
    TBEV               72      65  0.7010  0.6864  0.7462^ [0.6250, 0.8673]  0.4923  0.7330  0.5764  0.4203
    HCoV-229E         120      88  0.6563  0.7401  0.7114^ [0.6087, 0.8140]  0.4227  0.5731  0.5951  0.5908
    Zika              748   2,577  0.7918  0.7529  0.6812^ [0.6470, 0.7154]  0.3624  0.5749  0.6316  0.5534
    West Nile         261     311  0.7319  0.7219  0.6600^ [0.5932, 0.7268]  0.3200  0.6236  0.6308  0.5657
    DENV-3            213      95  0.7710  0.7151  0.6496^ [0.5554, 0.7439]  0.2993  0.2370  0.3861  0.1715
    DENV-1            195     123  0.7270  0.6795  0.6066^ [0.5139, 0.6994]  0.2132  0.3485  0.4489  0.3144
    DENV-4            153     109  0.7386  0.6617  0.5868  [0.4842, 0.6895]  0.1737  0.3588  0.4644  0.2302
    SARS-CoV          143      96  0.6074  0.5954  0.5785  [0.4701, 0.6869]  0.1571  0.5633  0.5186  0.6812
    HCoV-OC43         164      70  0.5834  0.5838  0.5677  [0.4509, 0.6845]  0.1354  0.6081  0.5396  0.5345
    MERS-CoV          133      41  0.6197  0.5958  0.4827  [0.3312, 0.6342] -0.0347  0.6076  0.6677  0.7085

    ^ clears the corrected alpha above chance.  NO CELL CLEARS IT BELOW.

**Ten of fourteen cells clear the corrected alpha above chance and none
below**, against five above and two below under the publication-disjoint
control.

## P1 held in the direction I did not predict

P1 asked whether SARS-CoV-2 would retain half its margin: 0.6346 or above.
It reads **0.8286 [0.8102, 0.8470]**, higher than the 0.7693 it read under
the control this replaces. The registration was written expecting the
symmetric control to be a sterner test; on this cell it is a less severe one,
and the prediction was satisfied from the wrong side. That is worth stating
plainly, because a prediction satisfied in an unanticipated direction is
weaker evidence than one satisfied as reasoned.

## P3 is the result, and the 'both' column is the mechanism

All four below-chance cells not only moved toward chance but crossed it:

    JEV        0.1452 -> 0.8002    +0.6550
    DENV-3     0.3861 -> 0.6496    +0.2636
    DENV-1     0.4489 -> 0.6066    +0.1577
    DENV-4     0.4644 -> 0.5868    +0.1224

DENV-3's movement read +0.2635 in the first draft of this paragraph, because
it was computed by subtracting the two DISPLAYED four-decimal values rather
than the stored ones. The other three are unaffected at this precision. A
figure derived from a rounded figure is the defect
`scripts/diagnose_doc_figures.py` exists to catch in documents, and it
appeared here in the arithmetic instead.

**And re-adding the publication exclusion puts them back.** The `both`
column keeps both controls, and JEV returns to 0.1608, DENV-3 to 0.2370,
DENV-1 to 0.3485, DENV-4 to 0.3588 — within 0.03 of their prior values in
three of four. So the publication exclusion is what produced the
anti-ranking, not the chemistry, and the prior verdict was right to refuse
the below-chance cells as a finding.

**The direct test of the asymmetry hypothesis now clears the bar it
previously missed.** Spearman between a cell's prior exclusion asymmetry
(log of positives-per-candidate over negatives-per-candidate) and its change
in AUC is **+0.635** over 14 cells, against a ~0.544 critical value. The
prior run could only offer -0.490 against the AUC level, which did not
clear; measured against the CHANGE, the association is there.

## The confound this registration did not anticipate

**The two controls differ in strength as well as in symmetry, and this
design cannot fully separate them.** Publication-disjointness removed
between 3.16 and 26.97 references per positive depending on the cell; (3')
removes exactly 10. So in most cells the new control is also *weaker*, and
a weaker control admits more leak — which would raise AUC without any of
the symmetry argument being true. The registration fixed k on feasibility
grounds and never asked how k compared with what the old control removed.

**It is bounded, not eliminated, and here is the bound.** In three cells the
publication control removed FEWER than 10 references per positive, so (3')
is strictly stronger there — and the AUC still rose in all three:

    West Nile   publication 3.16/pos   vs knn 10   +0.0292
    Zika        publication 6.86/pos   vs knn 10   +0.0496
    SARS-CoV    publication 8.97/pos   vs knn 10   +0.0600

Three cells, and the gains are small (+0.03 to +0.06) beside the big movers.
So the strength confound cannot explain the direction, and it can plausibly
explain much of the magnitude. **The clean experiment is to set k per cell
to the number the publication control actually removed, matching strength
while keeping symmetry, and it must be registered before it runs.**

## P5 held, and it is the most informative cell in either run

The matched NCATS cell — one screen, one protocol, both classes off the same
plates — is **measurable for the first time**: 89 positives against 7,343
negatives, **0.6088 [0.5157, 0.7020]**, clearing the corrected alpha above
chance, skill 0.2177. The count-based control cannot empty a reference set,
which is exactly why it was registered as the gain.

**Both halves of P5 held: it is measurable, and 0.6088 is below the ChEMBL
cell's 0.8286 at the same k.** So the collection confound is real and large
— roughly 0.22 of AUC separates a source-matched population from a
plate-matched one — and a structure signal survives plate matching anyway.
That is the first result in this project where a channel clears a corrected
threshold on a population matched by construction rather than by a filter.

## What this licenses

- **A structure-derived scorer ranks antiviral actives above measured
  inactives, on a plate-matched population, at a corrected threshold.** That
  is a real finding and the first of its kind here.
- **It is bounded by the collection confound at roughly 0.22 AUC**, which is
  measured rather than argued, because the matched cell can now be read.
- **It still is not a mechanism.** Everything both registrations say under
  that heading stands.
- **The ChEMBL figures are inflated relative to the matched one** and should
  be quoted with the matched cell beside them, never alone.
- **MERS-CoV is the one cell that fell** (0.6677 -> 0.4827, no longer
  clearing) and the registration offers no account of it. 133 positives
  against 41 negatives, the smallest negative class in the family, and the
  publication control there was nearly symmetric (10.36/12.61). It is
  reported as unexplained.

## What would change the conclusion

The strength-matched control described above. If SARS-CoV-2 and the matched
NCATS cell both survive a control that removes as many references as
publication-disjointness did while removing the same number from each
class, the finding is secure. If they collapse, then k=10 was simply a
weaker test and this verdict is measuring its own leniency.

**IT IS NOW REGISTERED:** `docs/prereg_structure_scorer_strength_matched.md`,
signed 2026-11-11 before its code, with `k_cell = ceil(max(pub_per_positive,
pub_per_negative))` per cell -- 4 for West Nile up to 33 for SARS-CoV-2, all
fourteen feasible, the tightest keeping 6 references. **And it found one
thing this verdict got wrong and one it could not have known.** The claim
that MERS-CoV holds "the smallest negative class in the family" is false --
JEV's 27 is smaller than its 41, and JEV moved the opposite way, so
small-sample noise does not account for MERS-CoV's direction on its own;
both documents now say second-smallest. And **the plate-matched cell cannot
be strength-matched at all**: publication-disjointness removed 88.00
references per positive and 89.00 per negative there out of 88, so the
matching k is 89 and leaves -1. The 0.22 AUC collection-confound bound
therefore rests on a k=10 comparison, which the new study re-measures at
k=33 rather than confirms.
