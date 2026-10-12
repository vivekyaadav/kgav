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
clears chance (0.582 [0.520, 0.644]) and it is the null hypothesis.
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
