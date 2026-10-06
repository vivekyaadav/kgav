# Evaluation: what the metapath scores actually separate

## The claim

**Once evidence that shares a publication with its own label is withheld, no
metapath separates measured-active from measured-inactive compounds. M1 scores
0.499 [0.452, 0.546]; the degree-only baseline scores 0.528 [0.481, 0.575].
Both intervals span 0.5.**

Three earlier figures for the same metapath — 0.823, 0.726, 0.660 — were each
produced by a protocol that failed a control the next one applied. They are
reported here in full, because the sequence is the result.

## Three nested confounds, and what each removed

Each row below adds one control to the row above it. CIs are Hanley–McNeil.

| # | Protocol | pos | neg | cov_pos | cov_neg | **M1** | degree |
|---|---|---|---|---|---|---|---|
| 1 | `selective-only`, positives only | 816 | 4,874 | 68.8% | 4.0% | 0.823 [0.805, 0.841] | 0.661 |
| 2 | `all` — no selectivity filter | 1,730 | 4,874 | 49.4% | 4.0% | 0.726 [0.711, 0.741] | 0.613 |
| 3 | `verified-only`, **both** sides | 1,028 | 165 | 56.6% | 25.5% | 0.660 [0.619, 0.701] | 0.528 |
| 4 | **+ `--publication-disjoint`** | 1,028 | 165 | **1.1%** | 1.2% | **0.499 [0.452, 0.546]** | 0.528 |

### 1 -> 2: the filter was asymmetric

Row 1 filtered positives to compounds with a verified selectivity index and
left negatives whole. Characterised compounds carry more edges, so the positive
class was 17× more reachable than the negative class. Remove the filter
entirely and 0.823 becomes 0.726 — most of the gap was annotation density.

### 2 -> 3: equalising the evidence basis

`--selectivity-applies both` requires the same paired cytotoxicity evidence of
a compound before it may be a negative. The coverage ratio falls from 12.4× to
2.2× and M1 to 0.660. The degree baseline falls from 0.613 to 0.528, which is
the control working as intended: degree *is* study volume, so equalising
characterisation should remove most of it.

At this point M1's margin over the baseline looked intact (+0.113 -> +0.132)
and the intervals did not overlap. That reading was wrong, because a fourth
confound had not been tested.

### 3 -> 4: the evidence and the answer came from one paper

Of the 1,070 compounds carrying both an `INHIBITS` edge and a
`HAS_ANTIVIRAL_ACTIVITY_AGAINST` label where both sides cite a publication,
**1,015 — 94.9% — share at least one PMID.**

A paper reporting an Mpro IC50 reports a cell-based EC50 in the same table.
ChEMBL files them as two activity records. The ingest turns one into the edge
M1 traverses and the other into the label M1 is scored against. The AUC is
then retrieval, not prediction.

Neither existing split separates them. A compound-level split cannot: both
facts hang off one compound. The temporal split cannot: one paper, one date.

`--publication-disjoint` withholds from traversal every edge citing a paper
that also reported that compound's own label. It withheld 1,150 of 3,938
`INHIBITS` edges — 29% — and **M1's reach over the evaluated positives fell
from 56.6% to 1.1%**, from roughly 582 compounds to 11.

That coverage collapse, not the AUC, is the finding. 98% of the positives M1
could reach were reachable only through same-paper evidence.

## Nothing beats the null hypothesis

SARS-CoV-2, cross-sectional, protocol 4:

| M1 | M2 | M3 | M4 | M5 | M6 | M7 | COMBINED | degree |
|---|---|---|---|---|---|---|---|---|
| 0.499 | 0.493 | 0.497 | 0.481 | 0.476 | 0.496 | 0.509 | 0.487 | 0.528 |

`baselines.py` states the bar: *"A learned scorer that does not clearly beat
all three is not a result."* The best scorer here is the degree baseline, and
it does not clear chance either.

Two further readings:

- **Seasonal coronaviruses are degree, undisguised.** HCoV-OC43 degree 0.746,
  HCoV-229E 0.699, on ~80 tested compounds each. With a pool that small,
  "has this compound been studied" is a strong predictor and no mechanism is
  involved.
- **The temporal split is worse than chance.** degree 0.439, and M1 reaches
  0.0% of post-2021 actives from the pre-2021 graph. Consistent with the
  earlier finding that 77% of post-cutoff compounds are absent from the prior
  graph entirely.

## Why it cannot currently work — all six reasons, measured

| Cause | Measurement |
|---|---|
| Evidence co-reported with labels | 94.9% shared PMIDs (1,015/1,070) |
| `INHIBITS` is effectively one protein | 2,586 of 3,938 edges on `PRO_0000449623`; 19 distinct viral proteins of 343 |
| Negatives cannot reach host-directed paths | 1.5% in M2's pool vs 12.1% expected — 7.9× under-represented |
| The graph is two graphs | 45 of 12,821 compounds have both `INHIBITS` and `TARGETS` |
| M8 inert | no `MEMBER_OF_CLASS` or `FOLD_SIMILAR_TO` edges in the release |
| M3/M5 untestable temporally | both PPI layers wholly post-cutoff (STRING 2023-08-28, VirHostNet 2024-01-01) |

Note which of these are *evaluation* faults rather than graph faults. The graph
validates, the provenance is intact, and 500 host proteins that drugs target
are dependency host factors by CRISPR screen — a real repurposing surface. No
protocol here has ever scored against it, because no labels exist in that
population.

## What is not concluded

- **Not** that knowledge-graph repurposing fails. These protocols never tested
  prediction, so they cannot have refuted it.
- **Not** that host-directed metapaths are uninformative. They were scored
  against a negative set they structurally cannot reach.
- 165 negatives is thin; the interval is 0.094 wide against 0.030 for row 2.

## What would make it testable

Labels in the population where repositioning happens. An approved-drug
antiviral screening panel (ReFRAME, Broad Drug Repurposing Hub, NCATS OpenData)
supplies actives *and* inactives among approved compounds, which is where the
`TARGETS` edges are; and because such screens are generated independently of
ChEMBL target annotations, they carry no co-reporting leakage to withhold.

That is a labels problem, not a topology problem. More interactome added to a
graph nothing can be scored against changes none of the six rows above.

## Reproducing

```bash
# row 4 -- the current claim
python scripts/hard_negatives.py --selectivity verified-only \
    --selectivity-applies both --publication-disjoint

# rows 1-3, in order
python scripts/hard_negatives.py --selectivity selective-only
python scripts/hard_negatives.py --selectivity all
python scripts/hard_negatives.py --selectivity verified-only --selectivity-applies both

# no-op control: --selectivity all short-circuits the filter, so adding
# --selectivity-applies both must reproduce row 2 exactly. It does (0.726).
python scripts/hard_negatives.py --selectivity all --selectivity-applies both
```

Diagnostics behind the six causes: `diagnose_hop_attrition.py` (hop-by-hop
survivors) and `diagnose_label_overlap.py` (population overlap, shared PMIDs).
