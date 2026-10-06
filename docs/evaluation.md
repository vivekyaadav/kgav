# Evaluation: what the metapath scores actually separate

## The claim

**Against measured-inactive compounds with the same evidence basis, direct-acting
mechanistic paths (M1) separate actives from inactives at AUC 0.660 [0.619,
0.701], against 0.528 [0.481, 0.575] for a degree-only baseline.**

The two intervals do not overlap. That is the result.

An earlier figure of 0.823 for the same metapath is reported below and is not
withdrawn, but it is not the headline, because it compares two populations with
different amounts of evidence behind them and most of the gap is that
difference rather than biology.

## Why three protocols

Measured negatives come from ChEMBL censored values — a `>` relation above
10 µM, i.e. somebody tested the compound and it did not work. That is a far
better negative than a random compound. But it is not automatically a
*comparable* negative.

A positive set filtered to compounds with a verified selectivity index
(SI = CC50/EC50 from same-document pairing) is a set of unusually
well-characterised compounds. Characterised compounds have more edges. More
edges means more reachable paths. If the negative set is not filtered the same
way, the classifier can score above chance by detecting how much work has been
done on a compound, which is not the question being asked.

`scripts/hard_negatives.py --selectivity-applies {positives,both}` exists to
make that filter symmetric and measure the difference.

## Results

Metapath M1 (direct-acting: Compound -> viral Protein) and the degree-only
baseline, across three filters. `cov` is the fraction of each class reachable
by at least one path. CIs are Hanley–McNeil.

| Filter | pos | neg | cov_pos | cov_neg | ratio | **M1** | **degree** | **M1 − degree** |
|---|---|---|---|---|---|---|---|---|
| `all` | 1,730 | 4,874 | 49.4% | 4.0% | 12.4× | 0.726 [0.711, 0.741] | 0.613 [0.597, 0.629] | **+0.113** |
| `selective-only`, positives only | 816 | 4,874 | 68.8% | 4.0% | 17.2× | 0.823 [0.805, 0.841] | 0.661 [0.639, 0.683] | **+0.162** |
| `verified-only`, **both** | 1,028 | 165 | 56.6% | 25.5% | 2.2× | **0.660 [0.619, 0.701]** | **0.528 [0.481, 0.575]** | **+0.132** |
| `selective-only`, **both** | 816 | 62 | — | — | 1.1× | 0.551 [0.479, 0.623] | 0.500 [0.425, 0.575] | +0.051 |

Three things to read off it.

**AUC tracks the coverage ratio, not the filter's stringency.** 17.2× -> 0.823,
12.4× -> 0.726, 2.2× -> 0.660, 1.1× -> 0.551. The ordering of the absolute AUCs
is the ordering of how unevenly the two classes are annotated. This is the
confound, measured.

**The degree baseline collapses when the filter is symmetric** — 0.613 to 0.528,
within touching distance of chance. The degree baseline *is* "how well-studied
is this compound", so a symmetric filter removing almost all of its signal is
the control behaving exactly as a control should. It is also direct evidence
that the asymmetric comparison was partly measuring study volume.

**M1's margin over that baseline does not collapse. It widens** — +0.113 to
+0.132. The absolute number falls; the quantity that carries the scientific
claim does not. This is why the headline is stated as a margin.

The last row is degenerate and is listed only to be ruled out:
`--selectivity selective-only --selectivity-applies both` requires a *negative*
to have SI >= 10, i.e. to be a selective antiviral, which is close to a
contradiction. It leaves 62 negatives and both scorers sit at chance. It is not
a valid control.

### Control

`--selectivity all --selectivity-applies both` must be a no-op, because `all`
short-circuits the filter block entirely. It reproduces the `all` row exactly
(0.726). The flag does nothing it should not.

The 0.726 here against 0.729 quoted in earlier runs is release drift — the
label work moved the counts from 1,726/4,890 to 1,730/4,874 — not flag
behaviour.

## Host-directed metapaths

Under the symmetric filter, every host-directed channel is at chance:

| | M2 | M3 | M4 | M5 | M6 | M7 | COMBINED |
|---|---|---|---|---|---|---|---|
| AUC | 0.493 | 0.497 | 0.481 | 0.476 | 0.496 | 0.473 | 0.630 |
| cov_pos | 4.0% | 3.1% | 7.1% | 5.9% | 1.1% | 7.6% | 71.0% |
| cov_neg | 5.5% | 3.6% | 10.9% | 10.9% | 1.8% | 13.3% | 49.7% |

This reproduces the asymmetric result and is now the second protocol to show
it, so it is not an artefact of either filter.

Note COMBINED (0.630) scores *below* M1 alone (0.660). Mixing a channel that
carries signal with six that do not dilutes it. Any ranking that pools
metapaths uniformly is currently worse than using M1 by itself.

M3 and M5 additionally cannot be evaluated in the temporal split at all: both
PPI layers are entirely post-cutoff (STRING `2023-08-28`, VirHostNet
`2024-01-01`), so no edge they traverse exists before the 2021 split point.

## What is not claimed

- No prospective validation. Every number here is retrospective.
- 165 negatives is thin. The interval is twice as wide as the `all` row's
  (0.082 vs 0.030). Non-overlapping intervals is a conservative test for a
  paired difference, so the separation is if anything understated, but the
  point estimate should not be quoted to three decimals as if it were stable.
- Restriction-factor findings remain unsupported (see README).
- M8 (cross-family fold) contributes nothing to any row: it is inert in the
  current release for three independent reasons and emits no paths.

## Reproducing

```bash
# headline
python scripts/hard_negatives.py --selectivity verified-only --selectivity-applies both

# the asymmetric comparison, for contrast
python scripts/hard_negatives.py --selectivity selective-only

# no-op control: must reproduce the --selectivity all row
python scripts/hard_negatives.py --selectivity all --selectivity-applies both
```
