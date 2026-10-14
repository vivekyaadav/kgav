# Evaluation: what the metapath scores actually separate

## The claim

**Of seven reasoning channels, two were ever evaluated. Both are
indistinguishable from chance. The other five reach too few measured actives
for their AUCs to describe anything but the imputed floor.**

    evaluated   M4  0.507 [0.447, 0.568]   reaches 15 of 89 actives
                M5  0.516 [0.455, 0.577]   reaches 14 of 89

    not         M3  reaches  9    M2 reaches 5    M1 reaches 3
    evaluated   M6  reaches  1    M7 reaches 0    M8 emits no paths at all

    baseline    degree 0.582 [0.520, 0.645] — the only scorer clear of 0.5,
                and it is the null hypothesis, not a reasoning channel

So a calibrated ranking returns nothing, and that is the correct output. The
graph can report evidence for a compound it is asked about; it has no measured
basis for a shortlist.

**The temporal protocol was the first non-null this project produced. On
2026-10-12 it was found NOT DETERMINED, and on 2026-10-14 a tie-safe metric
found no result in it either.** M1's Hits@100 on the SARS-CoV-2 train2021
split can take any value from 0 to 8 of its 8 positives depending on the
order drugs were inserted into a dict. The ROC AUC does not read a cut and
is invariant to that order, and on the same split it reads **0.512
[0.49, 0.54] for COMBINED and 0.504 [0.48, 0.53] for degree** — both
containing chance. The section below reports all three numbers.

The one-line version of every table here:

> **No channel, on any protocol, on any virus, has an AUC whose interval
> excludes chance. The only scorer that ever did is degree, the null
> hypothesis, and only cross-sectionally.**
>
> In every cross-sectional virus where a reasoning channel cleared p<0.05,
> the degree baseline cleared it too -- 4 of 4, no exception. Exactly one
> virus and channel ever had a trustworthy significant result while degree
> was at chance -- SARS-CoV-2 M1 under the temporal protocol -- and its
> Hits@k is NOT DETERMINED: all 8 positives sit inside one tie block at the
> cut, so the hit count ranges over 0 to 8 and the lift over 0.00 to 2.86.
> The tie-safe AUC on the same split reads 0.512 [0.49, 0.54] for COMBINED
> against 0.504 [0.48, 0.53] for degree.
>
> **The temporal protocol therefore works and the signal is not there.** It
> removes the study-volume confound -- degree is at chance on it -- and
> under that clean condition no channel separates either. That is a
> different and better-supported statement than "everything is confounded".

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

**3 -> 4. Evidence and answer from one paper.** Of 1,066 compounds carrying
both an `INHIBITS` edge and an antiviral label, 1,012 have a publication on
both sides and **933 -- 92.2% -- share a PMID**. A paper reporting an Mpro
IC50 reports a cell-based EC50 in the same table, and ChEMBL files them as two
records: M1 traverses a measurement taken alongside the answer. Neither split
separates them, because both facts hang off one compound and carry one date.

The remaining **50 compounds (4.7%) carry no publication on the `INHIBITS`
side at all**, and those are outside the control rather than cleared by it:
`Graph.load(drop_pubs_shared_with=...)` withholds an edge that SHARES a PMID
with its label, so an edge with no PMID is never a candidate for withholding.
At 4.7% that does not move row 4, which is why the figure is stated rather
than qualified -- but the same measurement on DENV-2 reads 68%, and a virus
can look clean under this control by having nothing to be disjoint from.
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

## Why a p-value is not the test

`lift_pvalue` reports the exact hypergeometric tail: the probability of that
many hits under a **random ranking**. That is all it reports. It cannot see a
confound that inflates the scorer and the labels together -- which is the
confound this whole document is about, since 92.2% of compounds share a PMID
between their `INHIBITS` edge and their antiviral label, and **degree IS study
volume**.

So p-values make the *uncontrolled* protocol look full of discoveries.
`compare_metapaths.py` runs on unfiltered `build_labels` output, and there:

| virus | best channel | lift | p | degree | p | |
|---|---|---|---|---|---|---|
| SARS-CoV | M1 | 2.10 | <.0001 | **9.51** | **<.0001** | degree wins 4.5x |
| MERS-CoV | M4 | 11.86 | <.0001 | 10.87 | <.0001 | both significant |
| HCoV-229E | M4 | 8.48 | 0.021 | 4.48 | 0.030 | both significant |
| HCoV-NL63 | — | 0.00 | 1.000 | **36.15** | **0.027** | degree only |
| SARS-CoV-2 | M4 | 1.99 | 0.140 | 0.45 | 0.990 | neither |
| HCoV-OC43 | M4 | 7.40 | 0.129 | 2.20 | 0.230 | neither |

**Degree's own p-value is the diagnostic.** Where it is significant, the
labels track how much a compound has been studied, and a channel beating
random on those labels says nothing about mechanism -- only beating *degree*
does. `confound_warning()` prints the verdict under every table, both ways:

    CONFOUNDED            degree is itself significant -- SARS-CoV, MERS-CoV,
                          HCoV-229E, HCoV-NL63 cross-sectionally, and
                          MERS-CoV and SARS-CoV temporally
    BASELINE AT CHANCE    degree is indistinguishable from random AND a
                          TRUSTWORTHY channel is significant -- held for
                          SARS-CoV-2 temporal alone, and no longer holds
                          anywhere, because both rows are ~tied
    BASELINE NOT          the degree row's own Hits@k is an interval, so no
    DETERMINED            comparison against it can be read -- including
                          every ratio in the vs-deg column, which divides
                          by it. SARS-CoV-2 temporal, both policies
    (nothing)             neither, so there is nothing to report

Stating the good case out loud matters: printed only as an *absence* of
warning, the one condition every real result here depends on was left for the
reader to notice.

`trustworthy` is required, not just significance. The first version checked
only `p` and promoted HCoV-OC43 temporal M4: pool 2,545, **one** positive,
lift 25.45, p=0.0393 -- and 100/2,545 = 0.0393 exactly, so that p *is* the
sampling fraction for a single compound. The row already carried `(!)`.

**`hits_determined` is required too, and this function has now read one
field while ignoring the one beside it three times.** Having learned to read
`trustworthy` alongside `p`, it went on to read both while ignoring whether
the hit count was determined -- and printed "BASELINE AT CHANCE: degree lift
0.97 (p=0.5909) is indistinguishable from random here, so M1's lift 2.49
(p=0.0033) is not explained by study volume" for a split where M1's Hits@100
was 0–8 and degree's was 3–5. The one sentence in this project that promoted
a result was comparing two intervals as though they were points.

**It appears under both undated policies, and that is not a replication.** M1
reads pool 284, 8 positives, lift 2.49, p=0.0033 in `include` and `exclude`
alike, because M1 traverses `INHIBITS`, which is fully dated, so the undated
policy never touches its pool. Only degree moves (0.97 -> 2.72). Two
identical numbers from one pool are one observation -- and here they are two
identical UNDETERMINED numbers, since the same tie block straddles the cut
either way.

**The tables kept printing those points for another four commits.** The
paragraph above is right, and the fix it describes reached `confound_warning`
and the `hits` column (`hits_cell` prints `0-8` rather than a draw) — but
`lift` and `p` sat in the next two columns as point values the whole time,
because `lift_row` computes both from the same `hits_at_k` draw that
`hits_cell` refuses to print. So a row could read

    M1   206   19   7-19    9.2    1.95  < .0001

with its own hits column saying the count is unknown within [7, 19] and the
two columns beside it quoting a lift and a p computed from one number inside
it. On v0.1 that is Zika M1, and its real lift interval is
[7/9.2, 19/9.2] = **[0.76, 2.07]** — it contains 1.00, so the row is
consistent with no enrichment while reading, as printed, as the strongest
result in the regeneration: `p < .0001` against a Zika degree baseline that
is itself at chance (lift 0.85, p=0.6887), which is exactly the clean
non-null this project has been looking for. `COMBINED` is the same tie block:
[0.82, 2.24]. Both now print the interval, `p` prints `n/d`, and a note under
the table says when an interval spans 1.00. The M6-vs-M7 verdict refuses
rather than comparing two draws.

Both tables also carry a `vs deg` column, with `(~)` marking a ratio whose
denominator is itself at chance. SARS-CoV-2 cross-sectional is why: M4 reads
4.44x a degree lift of **0.45**, which is M4 lift 1.99 at p=0.14. A ratio
against a noisy denominator describes; it does not test.

Note SARS-CoV especially: M1 at p<0.0001 looks like the strongest result in
the project until you see that degree scores 9.51 on the same labels against
M1's 2.10. Reported alone, that p-value would have been the fifth wrong M1
number in this document.

## The temporal protocol: the first signal above the baseline

Pre-2021 graph, post-2021 measurements. **This protocol is
publication-disjoint by construction** -- evidence is pre-cutoff and the label
post-cutoff, so the 92.2% shared-PMID confound that took cross-sectional M1 to
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

The figures below are from 2026-10-12, after M1 was constrained to resolved
targets (`INHIBITS.target_resolution = chain_from_assay_description`), which
moved its pool from 296 to 284. The superseded run read M1 at lift 2.22,
p=0.020 on 6 hits; it is kept here rather than withdrawn, as the five earlier
AUCs are, because the sequence is the result.

| scorer | pool | pos | hits | exp | lift | p | |
|---|---|---|---|---|---|---|---|
| **M1** | 284 | 8 | **0–8** | 2.8 | 2.49 | 0.0033 | **~tied** |
| M2 | 855 | 7 | 2 | 0.8 | 2.44 | 0.193 | (!) |
| M4 | 2,701 | 13 | 2 | 0.5 | 4.16 | 0.081 | (!) |
| M6 | 186 | 4 | 3 | 2.2 | 1.40 | 0.369 | |
| M7 | 85 | 26 | 26 | 30.6 | 0.85 | — | (!) |
| COMBINED | 3,229 | 48 | 5 | 1.5 | 3.36 | 0.0154 | |
| degree | 4,514 | 232 | **3–5** | 5.1 | 0.97 | 0.5909 | **~tied** |

**M1'S RESULT IS NOT DETERMINED, AND NEITHER IS THE BASELINE IT BEAT.**
`rank()` sorts by score with no tiebreak and Python's sort is stable, so
equal scores emerge in dict insertion order, and `hits_at_k` reads
`ranked[:k]`. On this split not one of M1's 8 positives scores strictly
above the cut at k=100: all 8 sit inside a single tie block that straddles
position 100, together with at least 8 non-positives. So Hits@100 ranges
over the full 0 to 8, and the lift over 0.00 to 2.86 against an expectation
of 2.8. The reported 2.49 at p=0.0033 is one draw from that interval.

The degree baseline is undetermined too, at 3–5. So the comparison that
licensed the result -- channel significant, degree at chance -- was a
comparison of two intervals written as though they were points.
`confound_warning` printed exactly that sentence, and it now refuses: a
baseline whose own Hits@k is an interval cannot license anything, and that
includes every ratio in the `vs deg` column, which divides by it.

This is a property of M1 rather than of the protocol. M2, M4, M6, M7 and
COMBINED are all determined on the same split, and COMBINED reaches lift
3.36 at p=0.0154 with a determined hit count -- which made it the strongest
surviving row until the tie-safe metric was applied to it.

**THE TIE-SAFE AUC, ADDED 2026-10-14.** `labels.evaluate_against_negatives`
counts tied pairs as half, so it does not read a cut and is invariant to
insertion order. This protocol already held it, and held 3,630 measured
test inactives, and passed only `test_pos` to `lift_row`. Calling it:

| scorer | lift | p | AUC [95% CI] | |
|---|---|---|---|---|
| M1 | 2.49 | 0.0033 | not evaluable (8) | ~tied |
| M4 | 4.16 | 0.081 | 0.504 [0.48, 0.53] | (!) |
| M7 | 0.85 | — | 0.506 [0.48, 0.53] | (!) |
| **COMBINED** | 3.36 | 0.0154 | **0.512 [0.49, 0.54]** | |
| degree | 0.97 | 0.5909 | 0.504 [0.48, 0.53] | ~tied |

The bar was computed and written down before the run, in
`docs/prereg_m1_potency.md`'s successor `docs/prereg_temporal_auc.md`: at
COMBINED's 48 reached positives a Hanley-McNeil interval excludes 0.5 only
from **0.586** upward. COMBINED came in at 0.512 and missed by 0.074. Under
`undated=exclude` it reads 0.505 against degree's 0.506.

The same scorer reads lift 3.36 at p=0.0154 and AUC 0.512, and the two are
not in tension: the lift reads a cut and the AUC does not.

**Degree reads 0.504 on the same split, which is the point.** The temporal
protocol removes the confound it was built to remove -- degree is at chance
here, where cross-sectionally it was 0.582 and significant in four of four
viruses. Under that clean condition no channel separates either.

WHY M1 AND NOT THE OTHERS. On the full release M1 resolves **12 distinct
scores across 3,020 drugs**, with 2,441 of them (80.8%) sharing one value,
because `walk_metapath` multiplies node degrees and nothing drug-specific
ever enters: M1's score is a property of the TARGET, so every drug whose
only route is nsp5 earns the identical number. M4 and M5 resolve 1,835 and
1,202 distinct values over comparable pools. Letting measured potency break
those ties was pre-registered, implemented, measured and NOT ADOPTED --
`docs/prereg_m1_potency.md` carries the predictions and the verdict.

The same table's MERS-CoV rows do NOT qualify: M4 reaches p=0.0005 there, but
degree reaches p=0.026 on the same labels, so that virus is still inside the
confound even temporally.

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
| Evidence co-reported with labels | 92.2% shared PMIDs (933/1,012); a further 50 compounds have no `INHIBITS` publication and are outside the control |
| M1 was never purely direct-acting | 1,884 of 7,677 protein-target rows (24.5%) described a CELLULAR readout and were filed against a viral protein; now reattributed to the organism |
| `INHIBITS` is effectively two proteins | nsp5 and nsp3 dominate; `INHIBITS` reaches 18 distinct viral proteins in all |
| The graph is two compound populations | 44 of 18,075 compounds have both `INHIBITS` and `TARGETS` |
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
- **Not** that M1 works. Its temporal lift reads 2.49 at a nominal p=0.0033,
  and the hit count behind it is NOT DETERMINED -- 0 to 8 of 8 positives,
  all of them inside one tie block at the cut -- so there is no figure for
  multiple-comparison correction to be applied to. It also did not survive
  correction for the 28 rows evaluated alongside it when it read 2.22.
- **Not** that the significant cross-sectional rows are results. Degree is
  significant in every one of them, which is the signature of a protocol
  measuring study volume.
- No prospective validation. Every figure here is retrospective, the temporal
  split included: it simulates prospective use on data that already exists.

## What would strengthen it

More labelled compounds that the host-directed channels can reach. The NCATS
CPE screen supplied 89 actives and 7,343 measured inactives in the
approved-drug population, which is what made rows 4 and 5 possible at all; a
second independent panel would raise M4 and M5 from 15 reached actives to
something an interval can be read off.

## The figures in this document are not under `make check-results`

`provenance.py` stamps every results FILE with the release it was computed
from, and `check_results.py` rejects one whose release has moved. The numbers
written in this document's prose carry no such stamp, so they rot silently --
which is the stale-before-column problem one layer up from the one
`provenance.py` was written for.

It had already happened. These figures were recorded against an earlier
release and were wrong by the time anyone re-ran the diagnostic:

| was | is | why it moved |
|---|---|---|
| 94.9% shared PMIDs (1,015/1,070) | 92.2% (933/1,012) | the assay-readout fix removed 540 `INHIBITS` edges |
| 45 of 12,821 with both edges | 44 of 18,075 | same, plus a different denominator |
| `INHIBITS` reaches 19 proteins | 18 | same |

Each figure below names the script that re-derives it. Re-run that script
before quoting a number from this document, and treat a figure with no named
source as unverified.

| figure | re-derived by |
|---|---|
| shared PMIDs, no-publication population, compound overlap, proteins reached | `diagnose_label_overlap.py` |
| hop-by-hop survivors, druggable host factors | `diagnose_hop_attrition.py` |
| cellular readouts filed against proteins | `diagnose_assay_type.py` |
| which source dated a label on which side of the cutoff | `diagnose_temporal_leak.py` |
| per-channel AUC, coverage, evaluability | `hard_negatives.py` then `calibrate.py` |
| lift, p, vs-degree, the baseline read | `compare_metapaths.py`, `temporal_split.py` |
| flavivirus volumes and substrate concentration | `diagnose_flavivirus_join.py`, `diagnose_flavivirus_targets.py` |

## Reproducing

**Order matters, for two reasons.** `temporal_split.py` cuts the training
graph that `hard_negatives.py` reads, so the split comes first — and the
three `--undated` protocols write to *one* output directory, so the last one
run is the split that stays on disk. Every `hard_negatives.py` variant
likewise writes to one `hard_negatives.json`, so the last one run is the row
that stays. Both lists below therefore end on the run whose numbers the rest
of this document quotes.

```bash
# ---- after any rebuild of data/releases/v0.1, everything below is required.
# Each block regenerates files that `make check-results` audits; a block
# skipped leaves a file describing a graph that no longer exists.

# 1. the temporal split, which is also an input to step 3
python scripts/temporal_split.py --cutoff 2021 --undated exclude
python scripts/temporal_split.py --cutoff 2021 --undated computed_only
python scripts/temporal_split.py --cutoff 2021 --undated include   # LAST: the
                                                                   # split the
                                                                   # rest reads

# control: must refuse, naming the 19 already-known compounds. Writes a
# refusal record and exits 1 -- expected, and it writes no split, so it
# cannot disturb the one above.
python scripts/temporal_split.py --cutoff 2021 --undated include \
    --overlap-policy keep || true

# 2. per-virus baselines -- one file each, all eight taxa in config/viruses.yaml
for t in 2697049 694009 1335626 11137 277944 31631 290028 64320; do
    python scripts/run_baselines.py --virus "NCBITaxon:$t"
done

# 3. metapath comparison, and the hard-negative protocols
python scripts/compare_metapaths.py

# rows 1-4, in order
python scripts/hard_negatives.py --selectivity selective-only
python scripts/hard_negatives.py --selectivity all
python scripts/hard_negatives.py --selectivity verified-only --selectivity-applies both
python scripts/hard_negatives.py --selectivity verified-only \
    --selectivity-applies both --publication-disjoint

# no-op control: must reproduce row 2 exactly
python scripts/hard_negatives.py --selectivity all --selectivity-applies both

# row 5 -- the current protocol, run LAST so its numbers are the ones left
# on disk. Running it first, as this block used to, meant the no-op control
# overwrote it and hard_negatives.json held row 2 while being read as row 5.
python scripts/hard_negatives.py --selectivity all --publication-disjoint \
    --label-source infores:ncats-opendata
python scripts/calibrate.py          # which channels were evaluable, and why not

# 4. every results file, AND the derived release they depend on, must
#    describe the graph on disk
make check-results
```

`hard_negatives.py` will **skip** its temporal column, rather than compute
one, if `data/releases/v0.1-train2021` was cut from a different graph or
carries no stamp — which is why step 1 precedes step 3. The skip and its
reason are recorded in `hard_negatives.json` under `temporal_skipped`, and
the split's own policy under `train_split`.

Diagnostics behind the causes: `diagnose_hop_attrition.py` (hop-by-hop
survivors), `diagnose_label_overlap.py` (population overlap, shared PMIDs),
`diagnose_assay_type.py` (cellular readouts filed against proteins),
`diagnose_temporal_leak.py` (which source put a label on which side of the
cutoff).
