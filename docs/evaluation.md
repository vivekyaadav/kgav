# Evaluation: what the metapath scores actually separate

## The claim

**Of seven reasoning channels, two were ever evaluated. Both are
indistinguishable from chance. The other five reach too few measured actives
for their AUCs to describe anything but the imputed floor.**

    evaluated   M4  0.507 [0.447, 0.568]   reaches 15 of 89 actives
                M5  0.516 [0.455, 0.577]   reaches 14 of 89

    not         M3  reaches  9    M2 reaches 5    M1 reaches 3
    evaluated   M6  reaches  1    M7 reaches 0    M8 emits no paths at all

    baseline    degree 0.582 [0.520, 0.644] — the only scorer clear of 0.5,
                and it is the null hypothesis, not a reasoning channel

So a calibrated ranking returns nothing, and that is the correct output. The
graph can report evidence for a compound it is asked about; it has no measured
basis for a shortlist.

**The temporal protocol is a different story, and it is the first non-null
this project has produced.** It was blocked until 2026-10-11 and is reported
in its own section below.

Five earlier figures for M1 -- 0.823, 0.726, 0.660, 0.499, 0.513 -- are
reported below rather than withdrawn, because the sequence is the result. Each
was produced by a protocol that failed a control the next one applied.

## The sequence

Each row adds one control to the row above. CIs are Hanley-McNeil.

| # | Protocol | pos | neg | cov_pos | cov_neg | **M1** |
|---|---|---|---|---|---|---|
| 1 | `selective-only`, positives only | 816 | 4,874 | 68.8% | 4.0% | 0.823 [0.805, 0.841] |
| 2 | `all` — no selectivity filter | 1,730 | 4,874 | 49.4% | 4.0% | 0.726 [0.711, 0.741] |
| 3 | `verified-only`, **both** sides | 1,028 | 165 | 56.6% | 25.5% | 0.660 [0.619, 0.701] |
| 4 | + `--publication-disjoint` | 1,028 | 165 | 1.1% | 1.2% | 0.499 [0.452, 0.546] |
| 5 | + `--label-source` (one screen) | 89 | 7,343 | 3.4% | 0.7% | **not evaluable** |

**1 -> 2. The filter was asymmetric.** Row 1 filtered positives to compounds
with a verified selectivity index and left negatives whole, so the positive
class was 17x more reachable. Remove the filter and 0.823 becomes 0.726.

**2 -> 3. Equalising the evidence basis.** Requiring the same paired
cytotoxicity evidence of a negative drops the coverage ratio to 2.2x and M1 to
0.660. The degree baseline falls 0.613 -> 0.528, which is the control working:
degree *is* study volume.

**3 -> 4. Evidence and answer from one paper.** Of 1,070 compounds carrying
both an `INHIBITS` edge and an antiviral label where both cite a publication,
**1,015 -- 94.9% -- share a PMID**. A paper reporting an Mpro IC50 reports a
cell-based EC50 in the same table, and ChEMBL files them as two records: M1
traverses a measurement taken alongside the answer. Neither split separates
them, because both facts hang off one compound and carry one date.
Withholding that evidence took M1's reach over the evaluated positives from
56.6% to 1.1% -- 582 compounds to 11. **That coverage collapse, not the AUC,
is the finding.**

**4 -> 5. One screen, so both classes are matched by construction.** Pooling
ChEMBL and NCATS labels restores the confound by a new route: ChEMBL actives
are antiviral research compounds carrying viral-target annotations, panel
inactives are library compounds carrying host-target ones, so the classes
differ by provenance as much as by activity. Pooled, M4 reads 0.445 -- *below*
chance -- purely because the negatives are 4x more reachable than the
positives. Restricted to the NCATS CPE screen, both classes come off the same
plates under one protocol: matched structurally rather than by a filter.

### Control

`--selectivity all --selectivity-applies both` must be a no-op, because `all`
short-circuits the filter block entirely. It reproduces the row-2 figure, so
the flag does nothing it should not.

## Why the coverage floor exists

`evaluate_against_negatives` imputes a floor score for every compound a scorer
does not reach. That is deliberate and stays: a scorer evaluated only on what
it reaches chooses its own test set, which is how a channel with 3% coverage
posts a perfect AUC.

But it means a low-coverage channel's AUC is dominated by the imputation, while
its interval is computed from `n_pos` and implies far more evidence than the
channel touched. M1's 0.513 rests on **three** compounds out of 89.

The assay-readout fix is what made this visible: it removed 540 `INHIBITS`
edges, a seventh of the channel, and M1's AUC did not move by 0.001. Nothing
could move it, because almost nothing reaches it.

`MIN_REACHED_POSITIVES = 10` matches `hard_negatives.MIN_CLASS_SIZE` -- the
minimum class size that script already requires before evaluating a virus is
the minimum for evaluating one channel of it.

## The temporal protocol: the first signal above the baseline

Pre-2021 graph, post-2021 measurements. **This protocol is
publication-disjoint by construction** -- evidence is pre-cutoff and the label
post-cutoff, so the 94.9% shared-PMID confound that took cross-sectional M1 to
0.499 cannot occur. That is why it can show something the cross-sectional
protocol cannot.

It had been refusing to run. `audit()`'s L5 check found 19 SARS-CoV-2
compounds on both sides of the cutoff and aborted; the NCATS CPE layer is one
screen under one publication dated 2020, so every label it contributes routes
to train, and ChEMBL re-measured 19 of them in 2022-2024.
`diagnose_temporal_leak.py` split them 19/19 by source with no mixed case, so
they are now dropped from the TEST set -- already-known compounds are not
prospective cases -- and `--overlap-policy keep` reproduces the refusal.

### SARS-CoV-2, cutoff 2021, undated=include

231 of 1,575 post-2021 compounds exist in the pre-2021 graph: a **14.7%
ceiling on recall**. `p` is the exact hypergeometric tail, P(hits >= observed)
under random ranking -- the lift equivalent of the AUC intervals above.

| scorer | pool | pos | hits | exp | lift | p | |
|---|---|---|---|---|---|---|---|
| **M1** | 296 | 8 | 6 | 2.7 | **2.22** | **0.020** | |
| M2 | 855 | 7 | 2 | 0.8 | 2.44 | 0.193 | (!) |
| M4 | 2,701 | 13 | 2 | 0.5 | 4.16 | 0.081 | (!) |
| M6 | 186 | 4 | 3 | 2.2 | 1.40 | 0.369 | |
| M7 | 91 | 29 | 29 | 31.9 | 0.91 | — | (!) |
| COMBINED | 3,246 | 50 | 4 | 1.5 | 2.60 | 0.066 | |
| degree | 4,316 | 231 | 5 | 5.4 | 0.93 | 0.628 | |

**M1 is the first scorer in this project to beat the degree baseline on a
protocol that controls its confound.** Degree is at chance here (0.93,
p=0.63); M1 is at 2.22, p=0.020.

**What that does and does not establish.** 28 scorer/virus rows were evaluated
in this run, so a 0.05 threshold expects between one and two rows at p<0.05 by
chance alone. Bonferroni over 28 is alpha=0.0018, and **M1's 0.020 does not
survive it.** The result is suggestive and unreplicated: six hits, eight
reachable positives.

One row does survive correction -- MERS-CoV M4, 3 hits from 5 positives in a
pool of 2,682, **p=0.00048** -- and it is flagged untrustworthy, because its
expectation is 0.19 and lift 16.09 is unstable at that sample size. It is the
least likely row in the run and it rests on five compounds. A hypothesis for a
targeted test, not a finding.

### Two implementations of one rule, again

`temporal_split.py` carried its own copy of `baselines.lift_row` predating the
`and pool > k` fix, so M7 printed unflagged at pool 91 (SARS-CoV-2), 34
(SARS-CoV) and 14 (MERS-CoV) -- where `ranked[:k]` IS the whole pool and lift
is exactly pool/k. `compare_metapaths.py` flagged the same channel on the same
release. The driver now imports `lift_row`, and a test asserts it defines no
numeric scoring rule of its own.

This is the third figure in this project to have had two implementations with
only one of them fixed, after AUC 0.826 in `guards.py` and `tools.py` and the
coverage floor.

## Six measured causes

| Cause | Measurement |
|---|---|
| Evidence co-reported with labels | 94.9% shared PMIDs (1,015/1,070) |
| M1 was never purely direct-acting | 1,884 of 7,677 protein-target rows (24.5%) described a CELLULAR readout and were filed against a viral protein; now reattributed to the organism |
| `INHIBITS` is effectively two proteins | nsp5 4,745 and nsp3 2,269 of ~7,400 resolved chains |
| The graph is two compound populations | 45 of 12,821 compounds have both `INHIBITS` and `TARGETS` |
| M8 inert | no `MEMBER_OF_CLASS` or `FOLD_SIMILAR_TO` edges; `ingest_folds` writes a layer `LAYER_PRECEDENCE` does not list (strict xfail in `test_assemble`) |
| M3/M5 untestable temporally | both PPI layers wholly post-cutoff (STRING 2023-08-28, VirHostNet 2024-01-01); the protocol itself was blocked by L5 until 2026-10-11 |

Note which are *evaluation* faults rather than graph faults. The graph
validates, provenance is intact, and 500 host proteins that drugs target are
dependency host factors by CRISPR screen -- a real repurposing surface that no
protocol here has scored against, because almost no labelled compound reaches
it.

## What is not claimed

- **Not** that knowledge-graph repurposing fails. Five of seven channels were
  never evaluated; the two that were are thin at 14 and 15 reached actives.
- **Not** that the host-directed hypothesis is refuted. M4 and M5 are the
  channels that *were* tested, and both at chance on 15 compounds is weak
  evidence either way.
- **Not** that M1 works. Its temporal lift of 2.22 is nominally p=0.020 on
  six hits and does not survive correction for the 28 rows evaluated
  alongside it.
- No prospective validation. Every figure here is retrospective, the temporal
  split included: it simulates prospective use on data that already exists.

## What would strengthen it

More labelled compounds that the host-directed channels can reach. The NCATS
CPE screen supplied 89 actives and 7,343 measured inactives in the
approved-drug population, which is what made rows 4 and 5 possible at all; a
second independent panel would raise M4 and M5 from 15 reached actives to
something an interval can be read off.

## Reproducing

```bash
# row 5 -- the current protocol
python scripts/hard_negatives.py --selectivity all --publication-disjoint \
    --label-source infores:ncats-opendata
python scripts/calibrate.py          # which channels were evaluable, and why not

# rows 1-4, in order
python scripts/hard_negatives.py --selectivity selective-only
python scripts/hard_negatives.py --selectivity all
python scripts/hard_negatives.py --selectivity verified-only --selectivity-applies both
python scripts/hard_negatives.py --selectivity verified-only \
    --selectivity-applies both --publication-disjoint

# no-op control: must reproduce row 2 exactly
python scripts/hard_negatives.py --selectivity all --selectivity-applies both

# the temporal protocol
python scripts/temporal_split.py --cutoff 2021 --undated include
python scripts/temporal_split.py --cutoff 2021 --undated computed_only

# control: must refuse, naming the 19 already-known compounds
python scripts/temporal_split.py --cutoff 2021 --undated include \
    --overlap-policy keep

# every results file must describe the graph on disk
make check-results
```

Diagnostics behind the causes: `diagnose_hop_attrition.py` (hop-by-hop
survivors), `diagnose_label_overlap.py` (population overlap, shared PMIDs),
`diagnose_assay_type.py` (cellular readouts filed against proteins),
`diagnose_temporal_leak.py` (which source put a label on which side of the
cutoff).
