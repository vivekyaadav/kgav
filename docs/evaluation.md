# Evaluation: what the metapath scores actually separate

## The claim

**Of eight reasoning channels, four were ever evaluated. All four are
indistinguishable from chance. The other four reach too few measured actives
for their AUCs to describe anything but the imputed floor.**

The Guide to PHARMACOLOGY layer, added 2026-11-01, is why it is four rather
than two: it roughly doubled the reach of every host-directed channel — M4
from 15 positives to 32, M5 from 14 to 27 — and carried M2 and M3 over the
floor for the first time. **Twice the evidence, four channels instead of
two, and the same answer.** That is a better-supported negative result than
the one this document carried on 2026-10-31, not a weaker one.

<!-- BLOCK: claim_channels — asserted line by line by
     tests/test_claim_block.py against data/results/hard_negatives.json,
     the single-screen publication-disjoint protocol (row 5). Do not
     hand-edit; re-derive by running the Reproducing block.

     Everything here is derived: which channels are `evaluated` is
     reached >= kgav.calibration.MIN_REACHED_POSITIVES, the intervals are
     kgav.labels.auc_interval, the reach counts are
     round(coverage_pos * n_pos), and the three counts in the sentence
     above are the channel totals. "M8 emits no paths at all" was true
     once and stopped being true -- it reaches 82 drugs and one labelled
     active -- which is also why the first count is eight, not seven. -->

    evaluated   M4  0.546 [0.484, 0.608]   reaches 31 of 89 actives
                M5  0.546 [0.484, 0.607]   reaches 26 of 89
                M3  0.541 [0.479, 0.602]   reaches 16 of 89
                M2  0.507 [0.447, 0.568]   reaches 11 of 89

    not         M6  reaches  5    M1 reaches 3
    evaluated   M8  reaches  1    M7 reaches 0

    baseline    degree   0.577 [0.514, 0.639]   clear of 0.5, and it is the
                                                null hypothesis, not a
                                                reasoning channel
    ensemble    COMBINED 0.573 [0.510, 0.635]   reaches 39 of 89. Also clear
                                                of 0.5, and 0.004 below
                                                degree -- see the
                                                decomposition below before
                                                reading that as a tie

So a calibrated ranking returns nothing, and that is the correct output. The
graph can report evidence for a compound it is asked about; it has no measured
basis for a shortlist.

### Those AUCs are not comparable to each other, and the decomposition says why

**Every figure above is an IMPUTED AUC.** A compound no channel reaches is
ranked at the floor rather than dropped, so a scorer that reaches a larger
fraction of the positives than of the negatives scores above chance with no
ranking ability whatsoever. That free component is exactly computable --
`kgav.labels.auc_floor`, which reduces to `0.5 + (coverage_pos -
coverage_neg) / 2` -- and so is the best a given coverage could reach,
`auc_ceiling`. Between them they bracket what the ranking could possibly
have contributed.

<!-- TABLE: auc_decomposition — asserted cell-by-cell by
     tests/test_auc_decomposition.py against data/results/hard_negatives.json.
     floor and ceiling are kgav.labels.auc_floor and auc_ceiling called on
     the coverages the file records; skill is kgav.labels.auc_skill, which
     is (auc - floor) / (ceiling - floor) and `n/a` when a scorer had no
     span to work in. Do not hand-edit; re-derive by running the
     Reproducing block. -->

| scorer | cov_pos | cov_neg | floor | AUC | ceiling | skill |
|---|---|---|---|---|---|---|
| M1 | 3.4% | 0.7% | 0.5135 | 0.5135 | 0.5136 | 19.7% |
| M2 | 12.4% | 11.2% | 0.5056 | 0.5074 | 0.5126 | 26.3% |
| M3 | 18.0% | 10.1% | 0.5394 | 0.5405 | 0.5484 | 13.0% |
| M4 | 34.8% | 27.7% | 0.5358 | 0.5463 | 0.5840 | 21.7% |
| M5 | 29.2% | 21.0% | 0.5412 | 0.5456 | 0.5718 | 14.3% |
| M6 | 5.6% | 2.5% | 0.5156 | 0.5156 | 0.5163 | -0.9% |
| M7 | 0.0% | 0.4% | 0.4982 | 0.4982 | 0.4982 | n/a |
| M8 | 1.1% | 0.0% | 0.5056 | 0.5056 | 0.5056 | n/a |
| **COMBINED** | 43.8% | 30.9% | **0.5645** | 0.5725 | 0.6322 | **11.9%** |
| **degree** | 100.0% | 100.0% | **0.5000** | 0.5767 | 1.0000 | **15.3%** |

**Degree's floor is exactly 0.5 and COMBINED's is 0.5645.** Degree reaches
every compound in both classes, so it has no coverage asymmetry at all and
every one of its 0.0767 points above chance is ranking. COMBINED's margin
above chance is 0.0746, of which **0.0645 -- 89% -- is its coverage gap**
and 0.0081 is ranking. The two raw numbers are four thousandths apart and
they are not measuring the same thing.

**Normalised, the finding survives and gets a reason.** The fraction of the
available span each scorer used is 11.9% for COMBINED and 15.3% for degree,
and the four evaluated channels run 13.0% to 26.3%. **Every scorer on this
graph, the null hypothesis included, extracts between 12% and 27% of the
separation its own coverage allows, and none of them is distinguishable
from the others on that measure.** That is a sharper
statement than either "degree wins" or "they tie": the differences between
these AUCs are differences in reach, and the ranking quality behind them is
uniformly poor and uniformly similar.

**Three cells have a floor equal to their AUC to four decimals** -- M1, M6
and M8 -- so their entire score is coverage and `skill` is computed over a
span of 0.0001 or less, which is why M1 reads 19.7% and M6 reads -0.9%
without either figure meaning anything. `auc_skill` returns `n/a` only when
the span is exactly zero; a span this small should be read the same way, and
the `not evaluated` block above is the gate that already does so.

**And the same decomposition applies to the temporal table, where it
leaves nothing behind.** `scripts/temporal_split.py` records the coverages
as `auc_coverage_pos` and `auc_coverage_neg`, so nothing new had to be
measured.

<!-- TABLE: temporal_auc_decomposition — asserted cell-by-cell by
     tests/test_auc_decomposition.py against the temporal results files.
     The rows are exactly the (policy, virus, scorer) cells whose
     `auc_ci_lo` exceeds 0.5, derived; on this release there are none, and
     the test asserts that the document says so rather than printing a
     table. floor and ceiling are kgav.labels.auc_floor and auc_ceiling on
     the recorded coverages. Do not hand-edit. -->

**No temporal cell has an AUC interval excluding 0.5 on this release**, so
this table is empty and the paragraph below is what the decomposition says
instead. Five cells were here on 2026-11-08 and all five were the
similarity leak.

**ACROSS ALL 69 TEMPORAL AUC CELLS, 53 EQUAL THEIR OWN IMPUTATION FLOOR TO
FOUR DECIMAL PLACES.** The median deviation from the floor is **exactly
zero**, the largest anywhere is 0.0067, and the largest floor-to-ceiling
span any scorer had to work in is 0.0320. On SARS-CoV-2 under `include`:

    scorer     cov_pos  cov_neg    floor      AUC   ceiling
    COMBINED      2.0%     0.3%   0.5086   0.5086    0.5086
    degree        4.1%     0.8%   0.5162   0.5162    0.5163

Both are their floors to four decimals. **The temporal AUC is a coverage
statistic.** A scorer on this split reaches 2-4% of the test positives and
a smaller fraction of the test negatives, the gap between those two numbers
is the whole figure, and ranking quality inside the reached set has at most
0.032 of room to show itself and demonstrably uses none of it.

That is a stronger negative than "no channel separates". It says the
tie-safe AUC, on this protocol, **cannot** separate: the statistic is
dominated by imputation before any scorer is applied, so a future channel
could rank its reached set perfectly and still read 0.509.

**Before the C4 fix, COMBINED read 0.526 [0.503, 0.549] here -- an interval
excluding chance, and equal to its floor of 0.5263 to four decimals, with a
ceiling of 0.5266.** It cleared the bar on a number its reach produced by
itself. Degree read 0.555 against a floor of 0.5526, and SARS-CoV degree
read 0.8343 -- the largest AUC anywhere in this project -- against a floor
of 0.8169 from an 88.4%-vs-25.0% coverage gap. All three are gone, and they
were never results.

**One correction to record here.** The replacement text written on
2026-11-08 said degree "reaches every compound so none of its figure is
imputed". That is true CROSS-sectionally, where degree's floor is exactly
0.5, and false temporally: degree can only reach what the pre-2021 graph
holds, so it has a coverage asymmetry like any channel. And a paragraph in
this document, in `docs/prereg_temporal_auc.md` and in the fix list said
`temporal_split.py` "records `auc_n_pos` and `auc_n_neg` but not the
coverages", so the row "cannot be told apart from its own imputation
floor". The fields were there under a different prefix; the claim was wrong
by one `grep`, and the thing it called unmeasurable is the sharpest
measurement in this section.

**The temporal protocol was the first non-null this project produced. On
2026-10-12 it was found NOT DETERMINED, and on 2026-10-14 a tie-safe metric
found no result in it either.** M1's Hits@100 on the SARS-CoV-2 train2021
split can take any value from 0 to 8 of its 8 positives depending on the
order drugs were inserted into a dict. The ROC AUC does not read a cut and
is invariant to that order, and on the same split (`--undated include`, the
one the rest of the pipeline reads) it reads **0.509 [0.486, 0.532] for
COMBINED and 0.516 [0.493, 0.539] for degree** — both containing chance,
and both equal to their own imputation floor to four decimal places. The
section below reports all three numbers and the floor under each.

**THIS PAIR OF FIGURES HAS NOW BEEN WRONG, RIGHT FOR THE WRONG REASON, AND
RIGHT, IN THREE DAYS.** It read "0.512 for COMBINED and 0.504 for degree —
both containing chance" until 2026-11-08, from a release two rebuilds
earlier that nothing re-derived; the temporal table was the one results
table in this document with no test behind it. Re-derived, it became 0.526
and 0.555, NEITHER containing chance — and that reversal was itself an
artifact, because the similarity layer was leaking post-cutoff compounds
into the 2021 graph (C4 below). With that fixed it is 0.509 and 0.516, both
containing chance again. The original CONCLUSION was right the whole time
and both of its numbers were wrong; one correction replaced them with
figures that were accurate and misleading, and only the decomposition and
the leak fix together got to a statement that holds. Every figure in that
sentence is now re-derived by `tests/test_temporal_table.py`.

The one-line version of every table here:

> **No individual reasoning channel, on any protocol, on any virus, has an
> AUC whose interval excludes chance.** Every scorer that does is either
> the degree baseline or the COMBINED ensemble of all eight channels.
> There are seven such cells on this release and they are all of those
> two scorers:
>
>     run                  virus        scorer    AUC [95% CI]           floor  from reach
>     cross-sectional      SARS-CoV-2   COMBINED  0.573 [0.510, 0.635]   0.5645       89%
>                          SARS-CoV-2   degree    0.577 [0.514, 0.639]   0.5000        0%
>
> **Two cells, and the third column is how much of each margin above
> chance is the imputation floor rather than the ranking.** Degree has no
> floor at all, because it is the only scorer here that reaches every
> compound in both classes; 89% of COMBINED's is reach.
>
> **Five of these cells were temporal until 2026-11-09 and are gone.**
> They were: SARS-CoV-2 COMBINED 0.526, SARS-CoV-2 degree 0.555 and
> SARS-CoV degree 0.834 under `include`, and the latter two again under
> `computed_only`. Every one of them was an artifact of the similarity
> layer leaking post-cutoff compounds into the 2021 training graph --
> fix-list C4, measured below -- and dating a computed edge by when its
> endpoints first appear removes all five. **Nothing in the temporal
> protocol now has an AUC interval excluding chance.**
>
> Cross-sectionally degree and COMBINED are **0.004 apart**, and that
> number is not the comparison it looks like either. Normalised for reach
> they use 11.9% and 15.3% of the span available to them, which is the
> defensible version of "indistinguishable". Both decomposition tables are
> above.
>
> **THIS CLAIM IS BACK TO WHAT IT SAID, AND THE ROUND TRIP IS THE POINT.**
> For weeks it read "two scorers do, both only cross-sectionally and both
> on SARS-CoV-2". On 2026-11-08 the census was derived for the first time
> and that was **false**: seven cells cleared chance, five of them
> temporal and two on SARS-CoV, and the pre-rebuild results files carried
> every one of them, so the sentence had never described any run. The
> claim was corrected to enumerate all seven. On 2026-11-09 fixing the
> similarity leak removed the five temporal cells, and the enumeration is
> two cells, both cross-sectional, both SARS-CoV-2 -- the original words,
> true for the first time.
>
> **A sentence that is right today and was wrong yesterday, with the same
> words, is the argument for deriving it.** Nothing re-derived this census
> before 2026-11-08; a sentence enumerating the cells could say anything,
> and for two releases it said something that happened to be false and is
> now accidentally correct. `tests/test_evaluation_table.py` derives the
> whole set from every results file and asserts this block lists exactly
> it, both ways, so neither state can be reached again by assertion.
>
> That second scorer was missing from this claim until 2026-11-07, and the
> claim was wrong without it: it read "the only scorer that ever did is
> degree". COMBINED's cross-sectional interval excluded chance on the
> shipped graph too -- 0.581 [0.518, 0.643] -- at 46% coverage of the
> positives, so it was neither a new result nor a degenerate cell. It was
> simply never printed anywhere in this document, and the rule here is
> that an unprinted figure is unchecked: five test files re-derive every
> table in this document from the results, and none of them could assert a
> number the document does not contain.
>
> It does not rescue the conclusion, and it sharpens it. An ensemble of
> eight channels clears chance, and so does counting a compound's edges,
> and neither of them ranks better than the other once reach is taken out
> -- 11.9% against 15.3% of the span each had available. What the graph
> separates is well-studied chemical matter from obscure chemical matter,
> which is what degree measures directly and what COMBINED's coverage
> asymmetry measures indirectly.
>
> **Until 2026-11-08 this read "still LOSES to counting a compound's
> edges", on a degree figure of 0.616 [0.553, 0.678].** That figure was
> inflated by the baseline counting the held-out label edges themselves:
> `degree_ranking` ranked each compound by its total degree, and 57.6% of a
> labelled compound's degree is its own `HAS_ANTIVIRAL_ACTIVITY_AGAINST`
> edges. Excluding them (L7) moved degree to 0.577 and left all eight
> channels and COMBINED unchanged to three decimals, because nothing but
> the baseline read those edges. The gap between the ensemble and the null
> hypothesis was 0.041 and is 0.002: degree was never beating the graph,
> it was reading the answer sheet.
>
> **Five cross-sectional viruses have a channel at p<0.05 and the degree
> baseline clears it in one of them** -- MERS-CoV. This sentence read "4 of
> 4, no exception" until 2026-11-08, which was wrong when written and is
> recorded as D2 in `docs/audit_2026-11-03_fixes.md`: it counted the
> viruses where DEGREE was significant and called them the viruses where a
> CHANNEL was. The five are SARS-CoV-2 (M7), DENV-2 (M1, M4, M5), Zika
> (M5), West Nile (M4) and MERS-CoV (M4, M5), and in four of them degree is
> nowhere near 0.05. That is not a defence of those four cells -- every one
> of them fails `trustworthy` or `hits_determined`, which is a stronger
> objection than the confound and is why the table prints them as
> em-dashes. Only SARS-CoV-2 M7 is readable, and its degree row is
> BASELINE NOT DETERMINED.
>
> Exactly one virus and channel ever had a trustworthy significant result
> while degree was at chance -- SARS-CoV-2 M1 under the temporal protocol
> -- and its Hits@k is NOT DETERMINED: all 8 positives sit inside one tie
> block at the cut, so the hit count ranges over 0 to 8 and the lift over
> 0.00 to 2.84.
>
> **The temporal protocol removes the confound on both metrics, and no
> channel separates under it either.** Degree's Hits@100 lift there is 0.93
> at p=0.6628 against a cross-sectional 10.02 at p<.0001 for SARS-CoV, and
> its tie-safe AUC is 0.516 [0.493, 0.539], containing chance. That is the
> control working, and it is a better-supported statement than "everything
> is confounded".
>
> **It was not true of the AUC between 2026-11-08 and 2026-11-09**, and
> this document said so at length: degree read 0.555 [0.533, 0.578] there,
> excluding chance, and the paragraph argued carefully that "degree is at
> chance on it" held of the LIFT and not of the AUC, because lift reads a
> cut and AUC does not. That distinction is sound and was the wrong
> explanation. The real cause was the similarity leak (C4): with
> post-cutoff compounds removed from the 2021 graph, degree's AUC falls to
> 0.516 and contains chance, so both metrics agree after all. **A correct
> piece of reasoning about the wrong data produced a wrong conclusion**,
> which is why the figures are now derived and the leak has a test.

Five figures for M1 -- 0.779, 0.673, 0.622, and then two protocols under which
the channel reaches too few positives to be evaluated at all -- are reported
below rather than withdrawn, because the sequence is the result. Each was
produced by a protocol that failed a control the next one applied.

## The sequence

SARS-CoV-2, cross-sectional. Each row adds one control to the row above. CIs
are Hanley-McNeil. `reached` is how many of the positives M1 actually has a
path to; below `MIN_REACHED_POSITIVES = 10` the AUC describes the imputed
floor rather than the channel, and the cell reads *not evaluable*.

<!-- TABLE: hard_negatives_sequence — asserted cell-by-cell by
     tests/test_sequence_table.py against data/results/hard_negatives_row1..4.json
     and data/results/hard_negatives.json (row 5, the canonical protocol).
     Do not hand-edit; re-derive by running the Reproducing block.

     Every cell is derived: pos and neg are n_pos and n_neg, cov_pos and
     cov_neg are the coverages to one decimal, `reached` is
     round(cov_pos * pos), and the M1 cell is the AUC with the interval
     `kgav.labels.auc_interval` computes from it -- or `not evaluable`
     when `reached` is below `kgav.calibration.MIN_REACHED_POSITIVES`.

     This table went four rows stale without anything noticing, because
     every protocol wrote one hard_negatives.json and only the last
     survived. Each row now writes its own file. -->

| # | Protocol | pos | neg | cov_pos | cov_neg | reached | **M1** |
|---|---|---|---|---|---|---|---|
| 1 | `selective-only`, positives only | 784 | 10,998 | 58.0% | 1.8% | 455 | 0.781 [0.762, 0.801] |
| 2 | `all` — no selectivity filter | 1,976 | 10,998 | 36.9% | 1.8% | 730 | 0.676 [0.662, 0.689] |
| 3 | `verified-only`, **both** sides | 960 | 185 | 48.8% | 23.8% | 468 | 0.623 [0.582, 0.664] |
| 4 | + `--publication-disjoint` | 960 | 185 | 0.5% | 2.2% | 5 | **not evaluable** |
| 5 | + `--label-source` (one screen) | 89 | 7,343 | 3.4% | 0.7% | 3 | **not evaluable** |

**These figures replace an earlier set -- 0.823, 0.726, 0.660, 0.499 -- that
the graph had stopped producing.** They were measured before the release was
rebuilt, and the measured negative set for rows 1 and 2 has since more than
doubled (4,874 to 10,953). Only row 5, whose labels come from one screen
untouched by that rebuild, reproduced unchanged. The sequence survives the
correction: the AUC still falls monotonically as each control is applied, and
the collapse at 3 -> 4 is sharper than before. Nothing was checking these
cells, which is why the drift was invisible: `check-results` compares a
results *file* against the graph and the file was current, while four of
these five rows had no file left at all, every protocol having written to one
`hard_negatives.json`. Each row now writes its own, and
`tests/test_sequence_table.py` re-derives every cell from it.

**1 -> 2. The filter was asymmetric.** Row 1 filtered positives to compounds
with a verified selectivity index and left negatives whole, so the positive
class was 32x more reachable. Remove the filter and 0.779 becomes 0.673.

**2 -> 3. Equalising the evidence basis.** Requiring the same paired
cytotoxicity evidence of a negative drops the coverage ratio to 2.0x and M1 to
0.622. The degree baseline falls 0.613 -> 0.527, which is the control working:
degree *is* study volume.

**3 -> 4. Evidence and answer from one paper.** Of 1,165 compounds carrying
both an `INHIBITS` edge and an antiviral label, 1,109 have a publication on
both sides and **1,020 -- 92.0% -- share a PMID**. A paper reporting an Mpro
IC50 reports a cell-based EC50 in the same table, and ChEMBL files them as two
records: M1 traverses a measurement taken alongside the answer. Neither split
separates them, because both facts hang off one compound and carry one date.

The remaining **50 compounds (4.3%) carry no publication on the `INHIBITS`
side at all**, and those are outside the control rather than cleared by it:
`Graph.load(drop_pubs_shared_with=...)` withholds an edge that SHARES a PMID
with its label, so an edge with no PMID is never a candidate for withholding.
At 4.3% that does not move row 4, which is why the figure is stated rather
than qualified -- but the same measurement on DENV-2 reads 68%, and a virus
can look clean under this control by having nothing to be disjoint from.
Withholding that evidence took M1's reach over the evaluated positives from
48.3% to 0.5% -- 473 compounds to **5**, which is below
`MIN_REACHED_POSITIVES`, so row 4 is no longer an AUC at all. **That coverage
collapse, not the AUC, is the finding** -- and the collapse is now severe
enough that the figure it used to produce has to be withdrawn rather than
quoted. On the figures this table previously carried, row 4 reached 11
compounds and scraped past the floor; it does not any more.

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
short-circuits the filter block entirely, leaving the flag recorded in the
output and otherwise unread.

This was asserted here for weeks on the strength of a step in the Reproducing
block that **could not have checked it**. Every `hard_negatives.py` variant
writes to one `hard_negatives.json`, so when the control ran, row 2's file had
been overwritten two steps earlier, and the control's own output was
overwritten one step later by row 5. Nothing compared anything; the run cost
about two and a half minutes per rebuild and was reported here in the past
tense as a check that had passed.

The comparison now runs in `tests/test_hard_negatives.py`, on every commit
rather than once per rebuild, against a fixture built for it: both classes
clear `MIN_CLASS_SIZE` under every arm and M2 reaches some compounds and not
others. The two configurations differ in exactly one field, the recorded
`selectivity_applies` — and requiring that field to differ is what shows the
two runs were really configured differently rather than compared with
themselves.

**The second arm is the one that makes the first mean anything.** Under
`--selectivity selective-only`, where the filter does run, the same pair of
settings must *diverge*: the negatives fall from 18 to 12 and the AUC moves.
Without it, the agreement would be satisfied by any fixture too small to
disagree — the smoke release yields `cross_sectional: []` and would pass by
comparing two empty documents.

## Why the coverage floor exists

`evaluate_against_negatives` imputes a floor score for every compound a scorer
does not reach. That is deliberate and stays: a scorer evaluated only on what
it reaches chooses its own test set, which is how a channel with 3% coverage
posts a perfect AUC.

But it means a low-coverage channel's AUC is dominated by the imputation, while
its interval is computed from `n_pos` and implies far more evidence than the
channel touched. M1's 0.514 rests on **three** compounds out of 89.

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
confound this whole document is about, since 92.0% of compounds share a PMID
between their `INHIBITS` edge and their antiviral label, and **degree IS study
volume**.

So p-values make the *uncontrolled* protocol look full of discoveries.
`compare_metapaths.py` runs on unfiltered `build_labels` output, and there:

<!-- TABLE: metapath_comparison — asserted cell-by-cell by
     tests/test_evaluation_table.py against data/results/metapath_comparison.json.
     Do not hand-edit; re-derive. The rule for the channel column is the one
     this document already requires of a reader: the single metapath with the
     highest lift among those that are BOTH `trustworthy` and
     `hits_determined`, and an em-dash when no channel qualifies. The verdict
     column is the first token of `confound_warning()`.

     The `significant` column names the LOWEST-p readable channel where one
     is significant at p<0.05, because the channel column selects on LIFT and
     can pass over a significant one: on 2026-10-30 SARS-CoV-2 read M4 at
     p=0.3266 while M7 sat at p=0.000379, the only cell on this graph that
     clears any correction family. Lift and significance answer different
     questions and this table was showing one while labelling it with the
     other's number. -->

| virus | best readable channel | lift | p | significant | degree | p | baseline verdict |
|---|---|---|---|---|---|---|---|
| SARS-CoV-2 | M4 | 1.79 | 0.1850 | M7 0.0013 | 0.51 | 0.9574 | BASELINE NOT DETERMINED |
| DENV-2 | — | — | — | — | 0.76 | 0.7623 | — |
| Zika | — | — | — | — | 0.67 | 0.8039 | — |
| West Nile | — | — | — | — | 0.97 | 0.6464 | — |
| DENV-3 | — | — | — | — | 2.38 | 0.2050 | — |
| DENV-1 | — | — | — | — | 2.60 | 0.1793 | — |
| HCoV-OC43 | — | — | — | — | 1.55 | 0.4777 | — |
| DENV-4 | — | — | — | — | 3.32 | 0.1222 | — |
| SARS-CoV | M3 | 2.37 | 0.0766 | — | 10.09 | <.0001 | CONFOUNDED |
| Yellow fever | — | — | — | — | 0.00 | 1.0000 | — |
| MERS-CoV | M2 | 1.39 | 0.4374 | — | 13.36 | <.0001 | CONFOUNDED |
| HCoV-229E | — | — | — | — | 4.20 | 0.0825 | — |
| TBEV | — | — | — | — | 0.00 | 1.0000 | — |
| JEV | — | — | — | — | 0.00 | 1.0000 | — |
| HCoV-NL63 | — | — | — | — | 0.00 | 1.0000 | — |

**JEV and HCoV-NL63 were printed TWICE in this table until 2026-11-08**, and
the test that re-derives it cell by cell could not see it: `parse()` keys
rows by virus name, so a duplicated row overwrites itself and set equality
against the results file still holds. The parser now has no row to collapse,
and the duplicate is recorded here because a machine-checked table that can
carry a repeated row is a weaker guarantee than it looks.

**Twelve of fifteen viruses have no readable channel at all**, and that is
the first thing this table should say. It was four of seven when this table
covered seven viruses; adding the eight remaining flaviviruses on 2026-10-23
did not dilute the finding, it concentrated it. Only four rows name a
channel, none is significant, and their lifts run 1.13 to 2.17.

A row is unreadable for one of two recorded reasons, and on this release
both are severe.

`hits_determined: False` means the published hit count is one draw from a
tie block straddling k, so the lift computed from it is one draw too. The
tie blocks here are not marginal. **DENV-2 M1** has a nominal lift of 1.61
at p=0.000222 -- the smallest p-value anywhere in this project -- and its
Hits@100 range is `[0, 41]` out of 41 positives, with 166 of a 236-compound
pool tied at the cut and every one of the 41 positives inside that block.
The hit count could as truthfully have been zero. **SARS-CoV M1** reads lift
0.76 at p=0.8914 with range `[2, 22]` of 22, 219 of a 239 pool tied, 20
positives in the block. **Zika M1** is `[7, 28]`, **West Nile M1** is
`[0, 5]` with the entire pool of 106 tied. Four of the five largest-looking
results in the table are a tie block being reported as a measurement.

`trustworthy: False` means `exp >= 1.0` fails -- the pool is large, the
positives are few, and the expectation at k=100 is a fraction of a hit, so
a single hit produces an enormous lift. **MERS-CoV M4** at 11.86 expects
0.51 hits; **HCoV-229E M4** at 8.48 expects 0.24; **HCoV-OC43 M4** at 7.40
expects 0.14; **Zika M5** at 5.23 expects 0.57. Every one of the eye-catching
double-digit lifts in this project is of this shape.

**Seven of fifteen viruses have no positive reachable by any metapath**
(`COMBINED.pos == 0`): DENV-1, DENV-3, DENV-4, JEV, TBEV, Yellow fever and
HCoV-NL63. That is not a weak result, it is the absence of a measurement,
and it comes in two different shapes that should not be read alike.
**DENV-4 and Yellow fever reach no compound at all** -- `COMBINED` pool 0 --
so nothing was ranked. **HCoV-NL63 reaches 2,660 compounds and not one of
them carries a label**, so a ranking exists, is complete, and cannot be
scored. DENV-1 and DENV-3 reach four compounds each, JEV five, TBEV fifty.
A ranking over these is not a weak prediction; it is a prediction nothing
in this graph can mark right or wrong.

That distinction is now load-bearing outside this document. The agent's
no-route guard refuses DENV-4 and Yellow fever, which are exactly the two
with an empty pool -- but it answers for DENV-1, DENV-3, JEV, TBEV and
HCoV-NL63, which have a route and no reachable positive. Edge counts are
what that guard reads, and edge counts cannot see this.

### The one readable cell that is nominally significant, and whether it survives correction

Across all fifteen viruses and eight metapaths, **eleven cells are both
`trustworthy` and `hits_determined`** -- readable, in the sense the rest of
this document means. One of them is nominally significant at p<0.05:

**SARS-CoV-2 M7, lift 1.29, p=0.001321**, reaching 225 positives in a pool
of 477. It is the only cell in this table that clears any threshold while
being readable.

**It no longer has a baseline to be compared against.** Until 2026-11-08
this sentence ended "against a degree baseline at 0.75 (p=0.8167) that is
itself indistinguishable from random", which was the whole argument for the
cell: a readable channel clearing 0.05 where the null hypothesis does not.
Stopping the degree baseline from counting held-out label edges (L7) moved
its lift to 0.63 (p=0.9088) and -- because the compounds that lost their
only degree now sit at zero -- put its Hits@k inside a tie block, so
`hits_determined` is `False`. The verdict on that row is therefore BASELINE
NOT DETERMINED, not BASELINE AT CHANCE, and M7's advantage over degree is
not a quantity this release can state. The channel's own figure did not
move; what moved is that there is nothing to divide it by.

**IT CLEARS ONE CORRECTION FAMILY OF FOUR, DOWN FROM THREE.** On the
release of 2026-10-28 this cell read p=0.000379 and cleared every family
but the strictest. The rebuild of 2026-11-07 -- which re-dated 12,514
virus-host edges from a constant to their publication year, calibrated the
ORCS effect-size sign, and removed 434 similarity edges joining two keys
of one molecule -- moved it to p=0.001655, and at that value it clears
only the most permissive family. Nothing about M7 changed; the graph it is
measured on became more accurate, and the cell got weaker. That is the
expected direction when a confound is removed, and it is recorded here
rather than in a commit message because it is the project's one
nominally-significant result.

**DENV-2 M4 was the second such cell on 2026-10-30 and is no longer one,
and the way it left is the more instructive half of this section.** It now
reads **lift 4.11 at p=0.014623** -- a HIGHER lift and a SMALLER p-value
than the 3.75 at p=0.019033 that put it here -- and it has dropped out of
the readable set anyway. The Guide to PHARMACOLOGY layer enlarged M4's
reachable pool for DENV-2 to 2,566 compounds, which cut the expected hit
count to **0.974**, just under the `MIN_TRUSTWORTHY_EXPECTATION = 1.0`
floor in `kgav.baselines`. Four hits against an expectation below one is
an arithmetic ratio, not a measurement, so `trustworthy` is now `False`
and the cell is excluded from the eleven -- while both of the numbers a
reader would quote from it moved in the direction that looks like progress.

Its earlier history is still worth recording. Before 2026-10-30 it read
lift 1.88 at p=0.2892, computed from a walk that stopped at an
undocumented 20,000-path cap -- and which of M4's paths it had seen by
then depended on the order edges.jsonl happened to be written in. 114 of
M4's 2,983 reaching drugs truncated that way. The cap was not only
inflating figures: there it suppressed one. It has now been a non-result
for two different reasons in three days, and at no point did it clear a
correction family.

It survives correction under most reasonable families and not the
strictest. This document has never declared one for this table, which is
itself the defect; the sensitivity is given rather than a choice made after
seeing the number:

<!-- TABLE: correction_sensitivity — asserted cell-by-cell by
     tests/test_correction_families.py against the `correction` block in
     data/results/metapath_comparison.json, which
     kgav.baselines.correction_families derives. Do not hand-edit.

     ALL FOUR SIZES WERE TYPED BY HAND until 2026-11-09 and nothing could
     recompute one. Three were right. The fourth, `cells_computed` at 61,
     had a definition written down nowhere: it is every (virus, metapath)
     cell the run PRODUCED, which is not every cell with a p-value (34)
     and not every cell planned (120). Three readings, three numbers, and
     the one the table meant was unrecorded. The verdict column is derived
     too -- `clears` iff the cell's p is below that family's alpha. -->

| family | n | alpha | M7 at p=0.001321 |
|---|---|---|---|
| readable cells only | 11 | 0.00455 | clears, by a factor of 3 |
| every cell computed | 61 | 0.00082 | **fails** |
| every cell planned (15 viruses x 8 metapaths) | 120 | 0.00042 | **fails** |
| planned, plus COMBINED and degree | 150 | 0.00033 | **fails** |

Correcting over cells that could be read is the defensible choice -- an
unreadable cell yields no inference, so it is not a test -- but the family
belongs in the prereg, before the next run, not here after this one.

**The temporal protocol now reports the same four families**, from the same
function, and they are much smaller because twelve of fifteen viruses fall
below the reachability floor: **2, 17, 24 and 30**, for thresholds of
0.02500, 0.00294, 0.00208 and 0.00167. No readable cell survives any of
them, which is a different statement from the four significant rows noted
earlier -- those four clear even the strictest family here and every one of
them fails `trustworthy`. **A row can survive a multiplicity correction and
still be unreadable**, and the correction is the weaker of the two gates.

One implementation, called by both scripts, for the reason this repository
keeps rediscovering: the paragraph above has quoted a family of 28, then
51, then 23 rows on three successive releases, each time re-derived by hand
from stdout.

**What M7 is matters more than its p-value.** M7 is `SmallMolecule
-CHEMICALLY_SIMILAR_TO-> SmallMolecule -INHIBITS-> Protein -ENCODED_BY->
Gene -BELONGS_TO-> OrganismTaxon`: a compound resembling a known inhibitor
of a viral protein. That is a QSAR baseline wearing a metapath's clothes.
It predicting activity is not news, it is the premise of medicinal
chemistry, and it is the channel in this schema that depends least on the
knowledge graph being a knowledge graph. The honest statement of the best
result here is: the only reasoning channel that beats a chance baseline on
readable data is the one that reasons least.

**The table names M4 on that row, not M7.** `best_readable` picks the
highest lift among readable channels, so it prints M4 at 1.79 (p=0.1850)
while the verdict beside it is about M7 at 1.29 (p=0.0013). A reader sees a
channel named, a verdict given, and no way to tell they are about different
things. Highest lift was the wrong selector on a page whose entire argument
is that lift without a determined hit count and a credible expectation is
not a result; the column is kept, because removing it would hide the
comparison, and the row now carries the significant channel beside it.

**Degree's own p-value is the diagnostic.** Where it is significant, the
labels track how much a compound has been studied, and a channel beating
random on those labels says nothing about mechanism -- only beating *degree*
does. `confound_warning()` prints the verdict under every table, both ways:

    CONFOUNDED            degree is itself significant -- SARS-CoV and
                          MERS-CoV cross-sectionally, and the same two
                          temporally under all three undated policies.
                          EIGHT (run, virus) pairs, the largest verdict
                          class in both tables. Zika held this verdict
                          temporally until 2026-11-09 and is now below
                          the reachability floor instead
    BASELINE AT CHANCE    degree is indistinguishable from random AND a
                          TRUSTWORTHY, DETERMINED channel is significant --
                          **NO (run, virus) pair on this release**. The
                          condition needs a determined degree row, and the
                          one pair that held it on 2026-11-07 --
                          cross-sectional SARS-CoV-2 -- lost the
                          determinacy rather than the chance
    BASELINE NOT          the degree row's own Hits@k is an interval, so no
    DETERMINED            comparison against it can be read -- including
                          every ratio in the vs-deg column, which divides
                          by it. TWO pairs, both SARS-CoV-2: the
                          cross-sectional run and the temporal `include`
                          run. It was four on 2026-11-08; DENV-2 under
                          include and SARS-CoV-2 under computed_only both
                          fell out when the similarity leak was fixed --
                          DENV-2 below the reachability floor, and
                          computed_only's degree hit count off its tie
    (nothing)             neither, so there is nothing to report

**BASELINE AT CHANCE now holds nowhere, and the way it left is the finding
of the 2026-11-08 regeneration.** It held for exactly one pair --
cross-sectional SARS-CoV-2, degree lift 0.75 (p=0.8167) against M7's 1.28
(p=0.0017) -- and the obvious reading of its disappearance would be that
degree stopped being at chance there. The opposite happened: degree's lift
fell from 0.75 to 0.63 and its p rose from 0.8167 to 0.9088, further from
significance, not closer. What it lost was `hits_determined`. Excluding
held-out label edges from the baseline (L7) left every compound whose only
edges were its own labels at degree zero, which built a tie block across
the Hits@100 cut, which made degree's hit count the interval `[5, 5+]`
rather than 6. A verdict that requires a readable denominator cannot be
issued against an unreadable one, so the pair moved to BASELINE NOT
DETERMINED. **The project's one significant readable cell is unchanged and
has lost its control**, which is a worse position than the one this block
described yesterday, not a better one.

**That is the second time this verdict's own census moved without any
baseline changing its mind, and the third time it changed table.** On
2026-11-02 the Guide to PHARMACOLOGY layer's 3,601 `TARGETS` edges gave
degree more neighbours to rank, four cross-sectional degree rows landed on
a tie block of 10 compounds straddling Hits@k, and those pairs fell from AT
CHANCE to NOT DETERMINED while SARS-CoV-2 under temporal `include` moved
the other way. On 2026-11-07 the rebuild re-dated the virus-host edges,
moved those four hit counts off the tie, and the verdict came back
cross-sectionally. On 2026-11-08 L7 removed it again. No verdict here was
edited to taste; every one is `confound_warning()` read off the current
results files, and the census is re-derived by
`tests/test_evaluation_table.py` in both directions -- a pair that holds it
must be named, and a claim that one holds it must fail when none does.

**An earlier version of this block said BASELINE AT CHANCE "no longer holds
anywhere", and that sentence is now true for a reason it was not true
then.** When written it was simply wrong -- the condition held for
SARS-CoV-2 and the artifact diff proved the claim was already false rather
than overtaken: regenerating across the 2026-10-16 `NSP_PATTERNS` fix moved
SARS-CoV and nothing else, leaving that entire block byte-identical. It was
never stale; it was never right. It is right today, and the difference
between the two states is exactly why the census is derived: the words are
the same and the thing they assert is not.

Its history is worth keeping. The condition first held for SARS-CoV-2 under
the temporal `computed_only` policy -- degree lift 1.15 (p=0.4262) against
COMBINED's 3.31 (p=0.0162). That run now reads NOT DETERMINED, because
degree's Hits@100 there is `[7, 8]`: a tie block straddles the cut with one
positive inside it, so its lift 0.86 and p=0.7180 are not measurements and
nothing can be compared against them.

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
[7/13.6, 28/13.6] = **[0.52, 2.06]** — it contains 1.00, so the row is
consistent with no enrichment while reading, as printed, as the strongest
result in the regeneration: `p < .0001` against a Zika degree baseline that
is itself at chance (lift 0.66, p=0.8132), which is exactly the clean
non-null this project has been looking for. `COMBINED` is determined on this
virus and reads 2.28 at p=0.1434. Both now print the interval, `p` prints
`n/d`, and a note under the table says when an interval spans 1.00. The
M6-vs-M7 verdict refuses rather than comparing two draws.

Both tables also carry a `vs deg` column, with `(~)` marking a ratio whose
denominator is itself at chance. SARS-CoV-2 cross-sectional is why, and on
this release there is **no ratio to print**: degree's Hits@100 there is an
interval, `hits_determined` is `False`, and `lift_vs_degree` is therefore
`None` for every channel on that virus -- not a number marked unreliable
but no number at all. M4 reads lift 1.79 at p=0.1850 and M7 1.29 at
p=0.0013, and neither can be expressed as a multiple of anything. The
published row does carry a degree lift of **0.51** (p=0.9574), and that
figure is one draw from the tie block rather than a measurement, which is
precisely why nothing is divided by it. A ratio
over a denominator at chance would be a statement about how little degree
achieved on this virus rather than evidence of a channel; a ratio over a
denominator that is an interval is not defined.

**This figure has now been quotable, withdrawn, quotable and withdrawn
again across four releases.** It read 2.37x until 2026-11-02, when the
Guide to PHARMACOLOGY layer's 3,601 `TARGETS` edges put degree's own hit
count inside a tie block and the ratio stopped existing. The rebuild of
2026-11-07 returned that count to a determined 6 and the paragraph quoted
2.04x. Stopping the baseline from counting held-out label edges (L7) sent
it back to an interval on 2026-11-08 -- 5 hits inside a tie block -- and
the ratio is withdrawn again. Four states in six days for one number,
which is the reason it is derived by a test rather than typed, and the
reason `tests/test_evaluation_table.py` refuses to let this document
express a channel's lift as a multiple of degree's anywhere, in any
wording, while `vs_degree_readable` is False.

Note SARS-CoV especially. M1 there reads p<0.0001 and looks like the
strongest result in the project -- and it is not readable at all:
`hits_determined` is `False`, its Hits@100 is one draw from [2, 22] of 22
positives, so its lift is one draw from [0.22, 2.39], an interval that
contains 1.00. An earlier version of this paragraph compared M1's 2.10
against a degree lift of 9.51 as though both were measurements; they were
two draws. Reported alone, that p-value would have been the fifth wrong M1
number in this document -- and reported as a comparison, it was the sixth.

**NO POINT ESTIMATE IS QUOTED FOR IT, AND THAT IS THE POINT.** This
paragraph said "lift 2.08" until 2026-10-30, when writing the release in a
canonical order changed nothing about the graph and moved the figure to
0.83. `rank()` sorts on score with no tiebreak and Python's sort is
stable, so which of 23 tied compounds fall inside Hits@100 is decided by
the order edges.jsonl happens to be written in. Both numbers are draws
from [0.22, 2.39]; quoting either one, and then saying in the next clause
that the interval contains 1.00, is the document disagreeing with itself
in two sentences.

**On the 16-virus release the degree row beside it is no longer the
problem.** It used to be undetermined as well ([7, 8] from a tie block), so
the virus read BASELINE NOT DETERMINED and no ratio against it could be
quoted at all. It is now determined at 6 hits, lift 10.02, p=3.1e-05, and
SARS-CoV reads CONFOUNDED. The reason M1 cannot be reported is unchanged
and sits entirely on M1's own side: a hit count of `[2, 23]` out of 23 is
not a measurement whatever the baseline beside it is doing. A row can stop
being unreadable for one reason while staying unreadable for another, and
the verdict column shows only the baseline's half of that.

## The channel partly IS the baseline, and that is a defect in the scorer

Everything above compares a channel against a degree baseline on the
assumption that they are different quantities. They are not entirely.

`walk_metapath` multiplies by each ENTERED node's `degree ** -w`. Two
consequences follow, neither of which is a judgement call:

**The start node's degree is never divided out.** Himmelstein's path-degree
product divides by the degree of every node on the path, source included, and
does it once per incident metaedge -- the reference implementation builds an
eight-element degree list for a four-hop metapath, at the same `w = 0.4` this
file uses. Here three of those six terms are applied and the compound's own
is not among them. So a compound with more edges accumulates more score. On a
controlled graph, ten targets against one scores **11.58x** under this
implementation and **5.13x** under the published one.

**`g.degree` counts every predicate.** The published form uses the degree
within the specific metaedge being traversed. Using the global degree means
an edge of a type the metapath never walks still moves the score: adding one
`CHEMICALLY_SIMILAR_TO` edge to a traversed protein moved a probe's M1 score
0.644394 -> 0.574349, a **10.9%** shift. Under the published form that probe
is invariant to the digit: 0.757858 before and after. Mechanically, then,
shipping the fold layer moved M1 scores without anything in M1 changing.

`docs/prereg_m1_potency.md` records the related half of this -- "nothing
drug-specific ever enters", hence 12 distinct M1 values across 3,020 drugs.
The two above are not recorded anywhere and are not what that prereg
proposes to fix.

### Measured on v0.1, by `scripts/diagnose_scorer.py`

**M4 is the degree channel, and the correction does not rescue it.** Spearman
between a channel's score and the compound's global degree:

| virus | M4 legacy | M4 pdp | M5 legacy | M5 pdp |
|---|---|---|---|---|
| HCoV-229E | +0.685 | +0.487 | +0.539 | +0.258 |
| MERS-CoV | +0.675 | +0.456 | +0.399 | +0.076 |
| SARS-CoV-2 | +0.685 | +0.492 | +0.496 | +0.289 |
| HCoV-NL63 | +0.620 | +0.284 | +0.257 | −0.095 |
| HCoV-OC43 | +0.669 | +0.470 | +0.444 | +0.138 |

Five independent viruses and M4's legacy figure sits in +0.62 to +0.69 on
every one. That consistency is structural, not biological. M4 is also the
channel the table above credited with its largest lifts before correction --
MERS-CoV 11.86, HCoV-229E 8.48, HCoV-OC43 7.40 -- all of which are
`trustworthy: False` on their own expectation.

**Two channels do not rank at all.** M8 returns ONE distinct value across
every drug it reaches (11 for SARS-CoV-2, 29 for Zika, 11 for SARS-CoV), as
does M1 on HCoV-229E across 17. An AUC computed for those is separating
reached from unreached, not ordering anything within reach.

**The scorings disagree about the ranking, sometimes badly.** Spearman
between the two scorings' own scores falls to +0.171 for M2 on HCoV-NL63 and
+0.310 for M3 on HCoV-HKU1. This is not a cosmetic rescaling.

**Four channels flip the sign, and it is not an artifact.** M2, M3, M6 and
M7 are positive under the legacy form and negative under pdp -- M2 on
HCoV-NL63 goes +0.256 to −0.827. The candidate explanation was that this
came from correlating against the GLOBAL degree while pdp divides by the
degree in the layer the channel walks, since a compound can be globally busy
and locally sparse. It does not: across the 36 channel-virus cells the two
correlations **agree in sign 36 times out of 36**, and for M2, M3 and M5
they agree to within 0.013 to 0.023 in median absolute difference. The flip
is a property of the scorer, not of which degree it is measured against.

**So neither form is degree-neutral, in opposite directions.** The legacy
form correlates +0.3 to +0.7 with degree; pdp correlates −0.04 to −0.87 on
25 of 36 cells. `w` is the knob between them -- the legacy form applies no
damping to the source and pdp applies `w = 0.4` -- and some intermediate
value would cross zero. Tuning `w` until the correlation vanishes would be
fitting the scorer to the confound, so it is noted here and not done.

**M4 is worse against the degree it actually uses.** Measured against the
`TARGETS` degree rather than the global one, M4's legacy correlation is
+0.831 / +0.814 / +0.830 / +0.700 / +0.812 across the five viruses it
reaches -- higher than the global figures above -- and it remains +0.317 to
+0.593 under pdp.

**Four cells are a perfect monotone function of the compound's degree.** M8
on SARS-CoV-2, Zika and SARS-CoV, and M1 on HCoV-229E, all read exactly
**−1.000** against the layer degree under pdp. Not "correlated with": rank
for rank, those channels ARE the compound's own edge count, inverted. They
carry no information beyond it, in either direction. M1 on MERS-CoV (−0.966)
and M7 on Zika (−0.955) are within noise of the same statement.

**What follows, and is registered in `docs/prereg_partial_degree.md`.**
"Does this channel beat degree" cannot be settled by choosing a scoring,
because both scorings are degree with a sign. It needs degree as an explicit
covariate in the same model as the channel -- which is what Rephetio did,
carrying node degrees for 14 metaedges as candidate features beside the
DWPCs.

### That measurement has now run. Verdict in the prereg; three results here

**No channel clears the bar.** Across 78 cells, 16 tested under measured
negatives and 26 under the looser reading, not one is significant after
correction under BOTH scorings with the same sign. The registered criterion
for proceeding to a full regularized model was not met.

**The closest cell is SARS-CoV-2 M4**, measured negatives: `b_c` +0.401
legacy and +0.403 pdp, raw p 0.0422 and 0.0093, corrected 0.6758 and
0.1488. Same sign, near-identical magnitude, and the only cell in the table
with that property -- and it fails Bonferroni, so it is not a result. It is
the one cell worth more labels.

**How a "negative" is defined decides what the model measures.** `b_d` on
SARS-CoV-2 is negative on 5 of 7 channels when negatives are MEASURED
inactive, and positive on 6 of 7 -- reaching +2.23 on SARS-CoV -- when
unlabelled compounds are counted as inactive. Counting the unlabelled as
inactive turns the model into a predictor of which compounds were ever
measured, which is the 92.0% shared-PMID confound wearing another costume.
That is why the measured-negative reading is primary.

The verdict also records a defect in the registered statistic: the
penalised likelihood ratio can come out negative, because the two models
carry the same lambda over different numbers of coefficients. It is
conservative -- those cells get p = 1 -- so the null survives it, but the
null is weaker than the corrected p-values suggest.

Nothing in `data/results` is affected yet: `baselines.SCORING_LEGACY` remains
the default and produced every number in this document. `SCORING_PDP` exists
beside it so the two can be measured against each other. Switching the
default is a separate decision that needs that measurement first, and a
registration.

**Neither form breaks the M1 tie block, and that is now evidence.** Two
compounds with one edge each to the same protein have equal degrees on both
sides of every metaedge, so they score identically under both. Measured:
M1's largest tie on SARS-CoV-2 goes 2,442 -> 2,380 and on Zika 195 -> 194.
Surviving the correction means those compounds have IDENTICAL `INHIBITS`
degree, which corroborates the prereg's "every drug whose only route is
nsp5" from a direction that prereg did not have. The degeneracy is a
different defect and potency weighting remains its remedy.

**Other channels' ties do collapse**, and that was not predicted. M2 on
HCoV-OC43 goes from a largest tie of 137 to 20, M3 on SARS-CoV from 114 to
16, M3 on SARS-CoV-2 from 81 to 13. The mechanism is the same one: pdp
divides by the compound's degree in the layer being walked, so compounds
differing in how many OTHER edges of that predicate they hold stop being
tied. M1's compounds hold one edge each and M2's do not. The expectation
written into `diagnose_scorer.py` before the first run -- that this column
would not move, and that if it did something was wrong -- generalised from
the single-edge case and was wrong; the docstring now records that.

## The temporal protocol: the first signal above the baseline

Pre-2021 graph, post-2021 measurements. **This protocol is
publication-disjoint by construction** -- evidence is pre-cutoff and the label
post-cutoff, so the 92.0% shared-PMID confound that took cross-sectional M1 to
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

**64 of 1,580 post-2021 compounds exist in the pre-2021 graph: a 4.1%
ceiling on recall.** `p` is the exact hypergeometric tail, P(hits >=
observed) under random ranking -- the lift equivalent of the AUC intervals
above.

That ceiling read **14.7% (231 of 1,575)** until 2026-11-09, and the
difference is the similarity leak: 167 of those 231 compounds "existed" in
the pre-2021 graph only through a `CHEMICALLY_SIMILAR_TO` edge computed
over the whole release and stamped with the placeholder date. A prospective
protocol cannot see 14.7% of its test set; it can see 4.1%. Every figure
below is on the 4.1%.

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
test inactives, and passed only `test_pos` to `lift_row`. Calling it on the
release of 2026-11-08:

<!-- TABLE: temporal_include_auc — asserted cell-by-cell by
     tests/test_temporal_table.py against
     data/results/temporal_2021_include.json, SARS-CoV-2. Do not
     hand-edit; re-derive by running the Reproducing block.

     This table had NO test until 2026-11-08 and every cell in it was
     stale: it read M1 2.49/0.0033, M4 4.16/0.081, M7 0.85, COMBINED
     3.36/0.0154 with AUC 0.512, and degree 0.97/0.5909 with AUC 0.504 --
     figures from a release two rebuilds earlier. The claim block at the
     top of this document quoted that COMBINED/degree pair as the summary
     of the whole temporal protocol. The marker column is derived: `~ties`
     when hits_determined is False, `(!)` when trustworthy is False. -->

| scorer | pool | pos | lift | p | AUC [95% CI] | |
|---|---|---|---|---|---|---|
| M1 | 284 | 8 | 0.36 | 0.9706 | 0.501 [0.478, 0.524] | ~ties |
| M2 | 1,200 | 17 | 1.41 | 0.4209 | 0.505 [0.482, 0.528] |  |
| M3 | 1,146 | 13 | 2.64 | 0.0970 | 0.504 [0.481, 0.527] |  |
| M4 | 3,870 | 27 | 2.87 | 0.1532 | 0.509 [0.486, 0.532] | (!) |
| M6 | 255 | 0 | 0.00 | — | 0.500 [0.477, 0.523] | (!) |
| M7 | 66 | 0 | 0.00 | — | 0.500 [0.477, 0.523] | (!) |
| M8 | 24 | 0 | 0.00 | — | 0.500 [0.477, 0.523] | (!) |
| **COMBINED** | 4,452 | 32 | 2.78 | 0.1608 | 0.509 [0.486, 0.532] | (!) |
| degree | 5,922 | 64 | 0.93 | 0.6657 | 0.516 [0.493, 0.539] | ~ties |

**Every figure in that table moved, and the two that matter moved in
opposite directions.** COMBINED's lift fell from 3.36 (p=0.0154,
significant) to 2.39 (p=0.0574, not) while its AUC rose from 0.512 to
0.526 -- and at 96 reached positives rather than 48, the larger sample
tightened the interval enough that **0.526 [0.503, 0.549] now excludes
0.5**. Degree's lift fell from 0.97 to 0.44 and its hit count became an
interval, while its AUC rose from 0.504 to 0.555, also excluding 0.5 and
by a wider margin than COMBINED. So the ensemble gained an AUC result and
lost its lift result in the same regeneration, and the baseline it is
compared against gained a larger AUC result than it did.

**THE PRE-REGISTERED BAR WAS COMPUTED FOR THE WRONG n, AND THAT IS WHY
COMBINED NOW CLEARS IT.** `docs/prereg_temporal_auc.md` worked out, before
the run, that a Hanley-McNeil interval excludes 0.5 only from **0.586**
upward "at COMBINED's 48 reached positives" against 3,630 negatives. The
AUC is not computed over the reached positives. The results file records
`auc_n_pos` and `auc_n_neg` directly, and on this split they are
**1,571 and 979** for every scorer: `evaluate_against_negatives` ranks the
whole labelled test set and imputes the floor to the compounds a scorer
cannot reach. At n=1,571/979 the interval half-width is 0.023, so the bar
-- the smallest AUC whose interval excludes 0.5 -- is **0.523**, not 0.586.
COMBINED's 0.526 clears it with 0.003 to spare. The prereg's arithmetic was done carefully and shown in
full; it was done on the reach count, which is the denominator of the LIFT.

**AND IT NO LONGER CLEARS THE BAR AT ALL, ON A GRAPH WITH THE LEAK
REMOVED.** COMBINED now reaches 32 of 1,571 test positives, 2.0%, and 3 of
979 test negatives, 0.3%, and reads **0.509 [0.486, 0.532]** -- below the
0.523 bar and containing chance. The arithmetic of its floor is exact: if
reached compounds all outrank unreached ones and are internally unordered,
the AUC is `0.5 + (coverage_pos - coverage_neg) / 2`.

    coverage gap   2.0% - 0.3%           =  1.72 points
    floor          0.5 + 0.0086          =  0.5086
    measured AUC                            0.5086
    ceiling, ranking the reached set
      PERFECTLY                             0.5086

**Floor, measurement and ceiling agree to four decimal places.** There is
nothing between what COMBINED's reach gives it for free and the best it
could possibly do, so the figure contains no ranking at all.

**On the leaking graph it DID clear the bar, and on exactly the same
arithmetic.** It read 0.526 [0.503, 0.549] against a floor of 0.5263 and a
ceiling of 0.5266: an interval excluding chance, a bar cleared by 0.003,
and not one thousandth of it from the ranking. **That is the clearest case
this project has produced of a protocol returning a significant result with
no information in it**, and it is recorded here rather than deleted because
the decomposition caught it one day before the leak fix did, by a different
route, and either alone would have been enough.

**Degree reads 0.516 [0.493, 0.539] on the same split, and has needed three
corrections.** The original text said "degree reads 0.504 on the same
split, which is the point -- the temporal protocol removes the confound it
was built to remove, degree is at chance here", on a figure from a release
two rebuilds earlier. On 2026-11-08 the re-derived figure was 0.555
[0.533, 0.578], which excludes chance, and the replacement text argued that
the protocol removes degree's LIFT advantage and not its AUC. It also
claimed degree "reaches every compound so none of its figure is imputed" --
true cross-sectionally, where degree's floor is exactly 0.5, and false
temporally, where degree reaches only what the pre-2021 graph holds. With
the leak fixed degree reads 0.516 against a floor of 0.5162, containing
chance, so the original CONCLUSION was right and every figure and every
mechanism offered for it was wrong.

What the protocol does remove, and this part has held throughout, is
degree's lift advantage: 0.93 at p=0.6628 here against 10.02 at p<.0001
cross-sectionally for SARS-CoV.

The pre-registered decision rule survives all of this intact, which is the
argument for having written it down: *"Claim a result for COMBINED only if
its interval clears 0.5 AND the degree baseline on the same split does
not."* On the leaking graph the first condition was met and the second was
not; on the fixed graph neither is. **No result is claimed**, by the rule
as written, under either graph, with no amendment and no judgement call
after seeing the number.

WHY M1 AND NOT THE OTHERS. On the full release M1 resolves **12 distinct
scores across 3,020 drugs**, with 2,441 of them (80.8%) sharing one value,
because `walk_metapath` multiplies node degrees and nothing drug-specific
ever enters: M1's score is a property of the TARGET, so every drug whose
only route is nsp5 earns the identical number. M4 and M5 resolve 1,835 and
1,202 distinct values over comparable pools. Letting measured potency break
those ties was pre-registered, implemented, measured and NOT ADOPTED --
`docs/prereg_m1_potency.md` carries the predictions and the verdict.

The same table's MERS-CoV rows do NOT qualify: M4 reaches p=0.0013 there, but
degree reaches p=0.0015 on the same labels at a nearly identical lift --
12.79 against M4's 12.83 -- so that virus is still inside the confound even
temporally, and the channel and the null hypothesis are not separable on it
at all.

**What that does and does not establish.** **23** scorer/virus rows were
evaluated in this run -- down from 51, because twelve of the fifteen viruses
now fall below the reachability floor -- so a 0.05 threshold expects about
one row at p<0.05 by chance alone. **Four** are below it: MERS-CoV M4,
COMBINED and degree, and SARS-CoV degree. Bonferroni over 23 is
alpha=0.00218, and **all four survive it**: p = 0.001276, 0.001131,
0.001452 and 0.000683.

**And not one of the four is readable.** Every one carries
`trustworthy: False`. Two of the four are the degree baseline, which is
the confound rather than a result, and in both viruses carrying the other
two degree is itself significant, so `confound_warning` returns CONFOUNDED
for MERS-CoV and SARS-CoV under all three undated policies. **Four rows
surviving a Bonferroni correction and none of them meaning anything is a
sharper statement than four rows failing it**, because it is the
readability gates rather than the multiplicity that refuse them, and those
gates were fixed before the p-values were in hand.

The protocol removes the confound on SARS-CoV-2, where degree drops to
lift 0.93 at p=0.6628 and an AUC of 0.516 containing chance, and does not
remove it on MERS-CoV and SARS-CoV, where degree still separates.

**These figures replace two earlier sets, and the shrinking family is the
story.** Before the 2026-11-07 rebuild: 28 rows, M1 at p=0.020, MERS-CoV M4
at p=0.00048 on 3 hits from 5 positives, lift 16.09. On 2026-11-08: 51 rows,
eight significant, only Zika COMBINED at p=1.95e-05 surviving correction.
Today: 23 rows, four significant, all four surviving and all four
unreadable. **Zika is not in this run at all** -- with the similarity leak
removed it reaches 7 pre-2021 compounds of 580 test positives, below the
floor of 10 -- so the row that survived correction yesterday was computed
on a virus this protocol cannot evaluate. All three sets are written down
because the paragraph's conclusion -- suggestive, unreplicated, resting on
a handful of compounds -- is the one thing that has not moved.

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

## The pathway landscape: nothing is enriched, and half of it was not asked

Everything above scores *compounds*. This section scores *pathways*, and it is
the one measurement in this document that was designed to come back negative
and did.

The question is whether the host proteins a virus depends on cluster in
particular Reactome pathways — because if they did, a pathway would be a
repurposing target in its own right, and a compound hitting any member of it
would inherit the hypothesis. Two channels ask it two ways. The **taxa**
channel takes the dependency factors asserted directly against the organism
(`HOST_FACTOR_FOR`). The **interactome** channel takes the host partners of
viral proteins (`PHYSICALLY_INTERACTS_WITH`), which is a larger and much
softer set.

**Why the null is degree-matched, and why that mattered more than the test.**
A protein in many pathways is also a protein in many screens, many
interaction datasets and many papers. Compare observed overlap against a
*uniform* random draw and well-studied proteins carry the signal; every
answer is then a restatement of how much attention a protein has already had.
So the null permutes within connectivity bins. An earlier run of this binned
against the whole graph, which is circular — the background included the
proteins under test — and it produced a uniform, confident-looking depletion
across every virus. That finding was an artifact and is retracted. The run
below uses `degree_basis = study`, binning within the studied set.

**The correction.** Each row is tested in both directions, enrichment and
depletion, so `n_tests = 44` over 22 rows and the Bonferroni threshold is
`alpha = 0.00114`. 2,000 permutations, seed 0.

<!-- TABLE: pathway_landscape — asserted cell-by-cell by
     tests/test_landscape_table.py against data/results/pathway_landscape.json.
     Do not hand-edit; re-derive with `make pathway-landscape`. The verdict
     column is ENRICHED when p enriched < alpha, DEPLETED when p depleted <
     alpha, and an em-dash otherwise, where alpha is 0.05/n_tests read from
     the file rather than restated here. -->

| channel | virus | dependency factors | observed | expected | p enriched | p depleted | verdict |
|---|---|---|---|---|---|---|---|
| taxa | SARS-CoV-2 | 3,460 | 384 | 359.19 | 0.0710 | 0.9390 | — |
| taxa | DENV-2 | 50 | 2 | 6.40 | 0.9890 | 0.0365 | — |
| taxa | West Nile | 18 | 0 | 2.55 | 1.0000 | 0.0610 | — |
| taxa | HCoV-NL63 | 57 | 10 | 7.34 | 0.1869 | 0.8911 | — |
| taxa | Zika | 815 | 76 | 80.37 | 0.7241 | 0.3168 | — |
| taxa | HCoV-OC43 | 657 | 62 | 65.29 | 0.6887 | 0.3568 | — |
| taxa | MERS-CoV | 678 | 60 | 66.39 | 0.8256 | 0.2164 | — |
| taxa | HCoV-229E | 1,367 | 128 | 134.91 | 0.7596 | 0.2694 | — |
| interactome | West Nile | 6 | 4 | 0.87 | 0.0055 | 1.0000 | — |
| interactome | HCoV-229E | 397 | 34 | 49.94 | 0.9990 | 0.0040 | — |
| interactome | HCoV-HKU1 | 279 | 25 | 36.71 | 0.9910 | 0.0165 | — |
| interactome | MERS-CoV | 622 | 64 | 79.69 | 0.9775 | 0.0310 | — |
| interactome | HCoV-OC43 | 289 | 28 | 38.25 | 0.9765 | 0.0365 | — |
| interactome | SARS-CoV | 1,238 | 167 | 149.51 | 0.0605 | 0.9500 | — |
| interactome | DENV-1 | 53 | 11 | 7.15 | 0.0860 | 0.9555 | — |
| interactome | HCoV-NL63 | 311 | 32 | 40.73 | 0.9440 | 0.0815 | — |
| interactome | SARS-CoV-2 | 3,494 | 438 | 426.23 | 0.2354 | 0.7796 | — |
| interactome | DENV-2 | 271 | 40 | 35.27 | 0.2154 | 0.8361 | — |
| interactome | TBEV | 35 | 3 | 4.59 | 0.8641 | 0.3063 | — |
| interactome | DENV-3 | 33 | 6 | 4.41 | 0.2814 | 0.8701 | — |
| interactome | Zika | 353 | 47 | 45.09 | 0.3938 | 0.6672 | — |
| interactome | JEV | 7 | 1 | 0.83 | 0.5937 | 0.8186 | — |

**0 of 22** rows clear the threshold, in either direction.

**The uncorrected view is where the discoveries live.** DENV-1's interactome
overlap is `0.0102 naive` and `0.0860 matched`; DENV-2's is `0.0036 naive`
and `0.2154 matched`. Neither gap opens because the data changed. The
naive p-value asks whether these proteins are over-represented in pathways;
the matched one asks whether they are over-represented *compared with
proteins as well-connected as they are*. A third-of-a-percent result becomes
nothing, and so does a one-percent one. That gap is the finding of this
section, and it is the same confound the rest of this document measures on
compounds.

**Every figure in this table moved on 2026-11-02, and the reason is worth
stating.** `n_observed` counts how many of a virus's dependency factors are
*drugged*, so the Guide to PHARMACOLOGY layer's 3,601 new `TARGETS` edges
enlarged the drugged proteome and every observed count rose with it --
SARS-CoV-2's taxa row from 333 to 384. The degree-matched expectation rose
too, from 308.24 to 358.98, because the same edges raised the degrees the
null draws against. Neither verdict changed, and nothing crossed the
threshold in either direction: a finding that survives a 7% increase in the
drugged proteome is a finding about the method, not about the sample.

**The two closest rows, and why neither is a result.** West Nile's
interactome is the smallest p in the table at 0.0065 — nearly six times the
threshold, and it rests on **six** interaction partners, four of which land
in a pathway where 0.88 were expected. Six. The agent attaches a thin-route
caveat to exactly this virus for exactly this reason. On the other side,
HCoV-229E reads depleted at 0.0075, more than six times the threshold.

**A pattern worth naming and not claiming.** Seven coronaviruses have an
interactome row, and the five with the smaller interactomes all lean
depleted: HCoV-229E 0.0075, HCoV-HKU1 0.0200, MERS-CoV 0.0275, HCoV-OC43
0.0320, HCoV-NL63 0.0695. The two SARS viruses lean the other way and
nowhere near significance (SARS-CoV-2 depleted at 0.8006, SARS-CoV at
0.9550), which is worth saying out loud: an earlier version of this
paragraph called these "the five coronaviruses with data" and left the
reader to assume there were only five. Same direction, five related viruses,
and four of the five would be significant uncorrected.
Not one clears correction. Five correlated tests pointing the same way is the
situation multiple-testing correction exists for, not evidence that survives
it. It is recorded here as an observation to re-examine when the interactome
grows, and nothing more.

**Half the register was never asked.** The taxa channel tested **8 of 16**
viruses; the interactome channel **14 of 16**. A negative result over half a
register is a weaker statement than a negative result over a register, and
the table above only lists what was testable. The two viruses missing from
*both* channels are **DENV-4** and **Yellow fever** — which is the same pair
`guards.check` refuses to rank, arrived at independently. The statistics and
the agent agree on which viruses this graph cannot reach.

**One figure in the table should not be read as it looks.** DENV-2's
interactome row counts 271 partners across 276 edges, of which **zero are
binary**. Every one is spoke-expanded from a complex, so "DENV-2 interacts
with 271 host proteins" means "DENV-2 was found in complexes containing 271
host proteins". It is the second-largest interactome row in the table and the
weakest-evidenced, and a reader comparing it against SARS-CoV-2's 3,494 is
comparing two different kinds of claim.

**What this section does not claim.** Not that viral dependency factors are
unstructured — a test that cannot reach half the register and rests on six
partners for one of its rows is not positioned to say so. Not that Reactome
is the wrong decomposition. Not that the depletion lean is real. The claim is
narrower: on this release, with this null, no pathway-level enrichment
survives correction, and the uncorrected version of the same test would have
reported four.

## Seven measured causes

| Cause | Measurement |
|---|---|
| Evidence co-reported with labels | 92.0% shared PMIDs (1,020/1,109); a further 50 compounds have no `INHIBITS` publication and are outside the control |
| M1 was never purely direct-acting | 1,884 of 7,677 protein-target rows (24.5%) described a CELLULAR readout and were filed against a viral protein; now reattributed to the organism |
| `INHIBITS` is effectively two proteins | nsp5 and nsp3 dominate; `INHIBITS` reaches 22 distinct viral proteins in all |
| The graph is two compound populations | 44 of 18,075 compounds have both `INHIBITS` and `TARGETS` |
| M8 walks, and reaches almost nothing | FIXED 2026-10-16: the layer assembles (104 `MEMBER_OF_CLASS`, 6 `FOLD_SIMILAR_TO`, 2 carrying `catalytic_type_match`). M8 now reaches 11 drugs for SARS-CoV-2, 29 for Zika, 11 for SARS-CoV and `no paths` for the other four coronaviruses; cross-sectional AUC 0.501 and 0.508 at 0.1% and 1.5% coverage, every row flagged. 34 of the 36 allowed transfers are blocked by comparator chains the spine does not build |
| M3/M5 untestable temporally | both PPI layers wholly post-cutoff (STRING 2023-08-28, VirHostNet 2024-01-01); the protocol itself was blocked by L5 until 2026-10-11 |
| The chemistry source is a medicinal-chemistry corpus | 5 journals supply 87.1% of ChEMBL 37's 93,489 documents; 19 named repurposing venues supply 2.5% of them and **7 documents in 2020-2021**. A congeneric series on one target is the shape it curates; a compound measured on two viral proteases is not |

Note which are *evaluation* faults rather than graph faults. The graph
validates, provenance is intact, and 500 host proteins that drugs target are
dependency host factors by CRISPR screen -- a real repurposing surface that no
protocol here has scored against, because almost no labelled compound reaches
it.

## What is not claimed

- **Not** that knowledge-graph repurposing fails. Four of eight channels
  were never evaluated; the four that were are thin, reaching 31, 26, 16 and
  11 of 89 actives.
- **Not** that the host-directed hypothesis is refuted. M4 and M5 are the
  best-reached channels, at 31 and 26 compounds, and both at chance on those
  numbers is weak evidence either way.
- **Not** that M1 works. Its temporal lift reads 0.36 at p=0.9706, and the
  hit count behind it is NOT DETERMINED -- 0 to 8 of 8 positives, all of
  them inside one tie block at the cut -- so there is no figure for
  multiple-comparison correction to be applied to. The same cell read 2.49
  at p=0.0033 before the 2026-11-07 rebuild and 2.22 at p=0.020 before that,
  and it did not survive correction at either value.
- **Not** that the significant cross-sectional rows are results. Degree is
  significant in one of the five viruses carrying them -- MERS-CoV -- and
  every one of the other rows fails `trustworthy` or `hits_determined`
  instead, which is a stronger objection than the confound. This bullet said
  "degree is significant in every one of them" until 2026-11-08, which was
  the D2 miscount.
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

**One field in that stamp was reporting on its own output.** `code.dirty`
is meant to say whether the CODE differed from the commit beside it, and it
was computed from a bare `git status --porcelain` -- but `data/results/*.json`
are tracked, so the first results file a run writes dirties the tree and
every stamp after it reads `True`. In the 2026-10-16 regeneration 13 of 14
files say `dirty: True` and the single `False` is simply the first written
(`temporal_2021_exclude.json`, 05:20:02; the next is 05:20:45). That clean
first stamp is positive evidence the run started from a clean checkout at
`e7d1451`, so those numbers are citable -- but the flag could not have told
us otherwise. The status check now excludes the generated-results directory.
**The 14 committed files keep the old flag**; the next regeneration records
it correctly.

It had already happened. These figures were recorded against an earlier
release and were wrong by the time anyone re-ran the diagnostic.

**Read the second column as a date, not as the present.** It records what
each figure read when the drift was found, and four of them have moved
again since: the shared-PMID share to 92.0% (1,020/1,109), the proteins
`INHIBITS` reaches to 22 on 2026-11-02, the BASELINE AT CHANCE row to two
temporal pairs on the same day and then to **nowhere at all** on
2026-11-08, and the six cross-sectional degree lifts a third time -- to
10.02 / 13.36 / 4.20 / 0.00 / 0.63 / 1.55 -- when the baseline stopped
counting held-out label edges. A history table whose "is" column is read
as current is the defect it is recording, one layer up. The
live values are the ones asserted against `data/results/*.json`, which is
what the rest of this section is about; these are kept because a project
that deletes its wrong figures cannot be checked on how it found them.

| was | became, when the drift was found | why it moved |
|---|---|---|
| 94.9% shared PMIDs (1,015/1,070) | 92.2% (933/1,012) | the assay-readout fix removed 540 `INHIBITS` edges |
| 45 of 12,821 with both edges | 44 of 18,075 | same, plus a different denominator |
| `INHIBITS` reaches 19 proteins | 18 | same |
| cross-sectional degree lifts 9.51 / 10.87 / 4.48 / 36.15 / 0.45 / 2.20 | 11.40 / 13.03 / 5.37 / 43.33 / 0.54 / 2.64 | every one of the six drifted across releases while the table sat unchanged; the 2026-10-16 `NSP_PATTERNS` fix moved only SARS-CoV, so the rest had been wrong for longer |
| "best channel" named for 6 of 7 viruses | readable for 3 | the rows were printed without checking `trustworthy` and `hits_determined`, the two fields this document requires of a reader |
| Zika absent from the table | present | it was never added; the guard test now asserts the table's virus set equals the JSON's |
| BASELINE AT CHANCE "no longer holds anywhere" | holds for SARS-CoV-2 under `computed_only` | not stale -- never true |

**The first table is no longer prose.** Every cell of the cross-sectional
lift table is now asserted against `data/results/metapath_comparison.json`
by `tests/test_evaluation_table.py`, including that the set of viruses in
the table equals the set in the JSON. The rule below -- name the script that
re-derives a figure -- leaves the reader to re-run it; for that one table the
suite re-runs it, so the next release that moves a degree lift fails CI
instead of quietly rotting. The rows above are what happens without that.

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
| the cross-sectional lift table above | `compare_metapaths.py`, then asserted cell-by-cell by `tests/test_evaluation_table.py` |
| flavivirus volumes and substrate concentration | `diagnose_flavivirus_join.py`, `diagnose_flavivirus_targets.py` |
| fold layer counts, allowed transfers, blocking chains, M8's reach | `ingest_folds.py`, `diagnose_comparator_cost.py` |
| journal profile, repurposing-venue coverage, rows by year per protease role | `diagnose_source_scope.py` |
| how much of a channel's score is the compound's own degree | `diagnose_scorer.py` |
| whether the channel beats degree once both are in one model | `partial_degree.py` |
| labels required for 80% power, against what the reach can supply | `label_headroom.py` |
| whether a reach ceiling is missing data or a retrieval limit | `diagnose_reach_ceiling.py` |

## Reproducing

**Order matters, for two reasons.** `temporal_split.py` cuts the training
graph that `hard_negatives.py` reads, so the split comes first — and the
three `--undated` protocols write to *one* output directory, so the last one
run is the split that stays on disk. Every `hard_negatives.py` protocol
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

# 2. per-virus baselines -- one file each, every taxon in
#    config/viruses.yaml. Six of the eight flaviviruses promoted on
#    2026-10-23 have no HOST_FACTOR_FOR edge, so their host-mediated
#    channels are empty by construction; the file is still written, and an
#    empty channel recorded is what lets a later ingest show up as a change.
#
#    ONE INVOCATION, not sixteen. dwpc_scores computes
#    {metapath: {drug: {virus}}} for EVERY virus at once and the script then
#    used one slice, so a loop recomputed the same traversal over 25,383
#    drugs sixteen times -- about 24 of this block's 60 minutes. --all-viruses
#    reads the register, so a virus promoted later is evaluated without
#    anyone editing this list, which went stale once already when it was
#    eight taxa long.
python scripts/run_baselines.py --all-viruses

# 3. metapath comparison, and the hard-negative protocols
python scripts/compare_metapaths.py

# rows 1-5, IN ORDER, in ONE process.
#
# Each --protocol string is parsed by hard_negatives.py's own parser, so
# nothing here could not be run on its own: copy any one of them onto the
# command line and it behaves exactly as it did when these were five
# commands. The results file is rewritten after each protocol, so the LAST
# one is still the row that stays on disk -- row 5, the protocol the rest
# of this document quotes. It used to run first, which meant the no-op
# control overwrote it and hard_negatives.json held row 2 while being read
# as row 5.
#
# ONE PROCESS because dwpc_scores depends on the GRAPH alone. These five
# differ only in their LABELS -- which compounds count as active, which as
# negative, which screen they came from -- so run separately they recompute
# the same traversal: five protocols over the full graph and the training
# graph is ten DWPC computations of which four are distinct. One invocation
# measured 2m41s on v0.1, so the five were about thirteen and a half
# minutes of which roughly six were duplicated work.
#
# The no-op control is NOT in this list. It used to be a sixth run, and it
# verified nothing: every protocol writes to one hard_negatives.json, so
# row 2's file had been overwritten two steps before the control ran and
# the control's own output was overwritten one step after. It cost ~2.5
# minutes a rebuild and nobody could read the comparison it was supposed to
# make. It runs in tests/test_hard_negatives.py instead, on every commit,
# with a second arm proving the fixture could have shown the two settings
# disagreeing.
#
# EACH ROW WRITES ITS OWN FILE. They all wrote hard_negatives.json, so only
# row 5 survived and rows 1-4 left no record at all -- which is how the
# table in "## The sequence" came to report four rows the graph had stopped
# producing, with check-results certifying the one file it could see. Row 5
# keeps the canonical name, because calibrate.py and the agent read it.
python scripts/hard_negatives.py \
    --protocol "--selectivity selective-only --results-name hard_negatives_row1.json" \
    --protocol "--selectivity all --results-name hard_negatives_row2.json" \
    --protocol "--selectivity verified-only --selectivity-applies both --results-name hard_negatives_row3.json" \
    --protocol "--selectivity verified-only --selectivity-applies both --publication-disjoint --results-name hard_negatives_row4.json" \
    --protocol "--selectivity all --publication-disjoint --label-source infores:ncats-opendata"
python scripts/calibrate.py          # which channels were evaluable, and why not

# 4. the pathway landscape. It depends on nothing above -- it reads the
#    release directly and scores no compounds -- but check-results globs
#    data/results and audits it like everything else, so leaving it out
#    meant a reader could run this block top to bottom and still be told
#    a file was stale. It was missing here from the day the file was added
#    until 2026-10-28; it looked fine only because it was being rebuilt by
#    hand. tests/test_reproduce_block.py now fails on a results file that
#    no step here regenerates.
make pathway-landscape

#    the shared-PMID confound, for the same reason. The 92.0% figure is why
#    this document calls its own M1 result partly self-scoring, and until
#    2026-11-02 it existed nowhere but that script's console output while
#    being quoted in six places here and one in kg_audit_report.txt.
#    Nothing could check any of them. It is a stamped results file now, so
#    it is audited below and tests/test_label_overlap_figures.py holds both
#    documents to it.
make label-overlap

# 5. every results file, AND the derived release they depend on, must
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
