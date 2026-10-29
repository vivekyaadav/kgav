# Evaluation: what the metapath scores actually separate

## The claim

**Of seven reasoning channels, two were ever evaluated. Both are
indistinguishable from chance. The other five reach too few measured actives
for their AUCs to describe anything but the imputed floor.**

    evaluated   M4  0.507 [0.447, 0.568]   reaches 15 of 89 actives
                M5  0.516 [0.455, 0.577]   reaches 14 of 89

    not         M3  reaches  9    M2 reaches 5    M1 reaches 3
    evaluated   M6  reaches  1    M7 reaches 0    M8 emits no paths at all

    baseline    degree 0.584 [0.522, 0.647] — the only scorer clear of 0.5,
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
| SARS-CoV-2 | M4 | 1.49 | 0.3266 | M7 0.0004 | 0.63 | 0.9077 | BASELINE AT CHANCE |
| DENV-2 | M4 | 3.75 | 0.0190 | M4 0.0190 | 1.21 | 0.3953 | BASELINE AT CHANCE |
| Zika | — | — | — | — | 0.99 | 0.5909 | BASELINE AT CHANCE |
| West Nile | — | — | — | — | 0.97 | 0.6464 | — |
| DENV-3 | — | — | — | — | 2.38 | 0.2050 | — |
| DENV-1 | — | — | — | — | 2.60 | 0.1793 | — |
| HCoV-OC43 | — | — | — | — | 3.10 | 0.1367 | — |
| DENV-4 | — | — | — | — | 3.32 | 0.1222 | — |
| SARS-CoV | M3 | 2.17 | 0.0974 | — | 11.69 | <.0001 | CONFOUNDED |
| Yellow fever | — | — | — | — | 1.88 | 0.4139 | — |
| MERS-CoV | M2 | 1.13 | 0.5661 | — | 15.27 | <.0001 | CONFOUNDED |
| HCoV-229E | — | — | — | — | 6.29 | 0.0122 | CONFOUNDED |
| TBEV | — | — | — | — | 3.53 | 0.2477 | — |
| JEV | — | — | — | — | 0.00 | 1.0000 | — |
| HCoV-NL63 | — | — | — | — | 50.77 | 0.0195 | CONFOUNDED |

**Eleven of fifteen viruses have no readable channel at all**, and that is
the first thing this table should say. It was four of seven when this table
covered seven viruses; adding the eight remaining flaviviruses on 2026-10-23
did not dilute the finding, it concentrated it. Only four rows name a
channel, none is significant, and their lifts run 1.13 to 2.17.

A row is unreadable for one of two recorded reasons, and on this release
both are severe.

`hits_determined: False` means the published hit count is one draw from a
tie block straddling k, so the lift computed from it is one draw too. The
tie blocks here are not marginal. **DENV-2 M1** has a nominal lift of 2.04
at p=3.6e-11 -- the smallest p-value anywhere in this project -- and its
Hits@100 range is `[0, 44]` out of 44 positives, with 166 of a 236-compound
pool tied at the cut and every one of the 44 positives inside that block.
The hit count could as truthfully have been zero. **SARS-CoV M1** reads lift
2.08 at p=3.9e-06 with range `[2, 23]` of 23, 219 of a 239 pool tied, 21
positives in the block. **Zika M1** is `[7, 29]`, **West Nile M1** is
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

### The two readable cells that are nominally significant, and why only one survives correction

Across all fifteen viruses and eight metapaths, **thirteen cells are both
`trustworthy` and `hits_determined`** -- readable, in the sense the rest of
this document means. Two of them are nominally significant at p<0.05:

**SARS-CoV-2 M7, lift 1.34, p=0.000379**, reaching 234 positives in a pool
of 515, against a degree baseline at 0.63 (p=0.9077) that is itself
indistinguishable from random. It is the only cell in this table that
clears any threshold while being readable, and the only reason the
cross-sectional verdict for SARS-CoV-2 is BASELINE AT CHANCE.

**DENV-2 M4, lift 3.75, p=0.019033**, which did not exist as a result
until 2026-10-30. Before that it read lift 1.88 at p=0.2892, computed from
a walk that stopped at an undocumented 20,000-path cap -- and which of
M4's paths it had seen by then depended on the order edges.jsonl happened
to be written in. 114 of M4's 2,983 reaching drugs truncated that way. The
cap was not only inflating figures: here it suppressed one, and the
complete traversal nearly doubled the lift while dropping p by a factor of
fifteen. It clears NONE of the correction families below -- the loosest is
0.00385 -- so it is recorded as a lead to examine, not a result.

It survives correction under most reasonable families and not the
strictest. This document has never declared one for this table, which is
itself the defect; the sensitivity is given rather than a choice made after
seeing the number:

| family | n | alpha | M7 at p=0.000379 |
|---|---|---|---|
| readable cells only | 13 | 0.00385 | clears, by a factor of 10 |
| every cell computed | 91 | 0.00055 | clears |
| every cell planned (15 viruses x 8 metapaths) | 120 | 0.00042 | clears narrowly |
| planned, plus COMBINED and degree | 150 | 0.00033 | **fails** |

Correcting over cells that could be read is the defensible choice -- an
unreadable cell yields no inference, so it is not a test -- but the family
belongs in the prereg, before the next run, not here after this one.

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
highest lift among readable channels, so it prints M4 at 1.99 (p=0.1401)
while the verdict beside it is about M7 at 1.34 (p=0.0004). A reader sees a
channel named, a verdict given, and no way to tell they are about different
things. Highest lift was the wrong selector on a page whose entire argument
is that lift without a determined hit count and a credible expectation is
not a result; the column is kept, because removing it would hide the
comparison, and the row now carries the significant channel beside it.

**Degree's own p-value is the diagnostic.** Where it is significant, the
labels track how much a compound has been studied, and a channel beating
random on those labels says nothing about mechanism -- only beating *degree*
does. `confound_warning()` prints the verdict under every table, both ways:

    CONFOUNDED            degree is itself significant -- SARS-CoV,
                          MERS-CoV, HCoV-229E and HCoV-NL63
                          cross-sectionally, and MERS-CoV, SARS-CoV and
                          Zika temporally (both the computed_only and
                          include policies)
    BASELINE AT CHANCE    degree is indistinguishable from random AND a
                          TRUSTWORTHY, DETERMINED channel is significant --
                          four (run, virus) pairs, where before the
                          complete M4 traversal of 2026-10-30 it was one:
                          SARS-CoV-2 cross-sectionally, degree lift 0.63
                          (p=0.9077) against M7's 1.34 (p=0.0004);
                          DENV-2 cross-sectionally, degree lift 1.21
                          (p=0.3953) against M4's 3.75 (p=0.0190);
                          Zika cross-sectionally, degree lift 0.99
                          (p=0.5909) against COMBINED's 5.98 (p=0.0000);
                          and SARS-CoV-2 under temporal computed_only,
                          degree lift 1.08 (p=0.4693) against COMBINED's
                          2.87 (p=0.0170)
    BASELINE NOT          the degree row's own Hits@k is an interval, so no
    DETERMINED            comparison against it can be read -- including
                          every ratio in the vs-deg column, which divides
                          by it. SARS-CoV-2 under temporal include. It
                          held under computed_only too until 2026-10-30;
                          degree's hit count there is now determined and
                          that pair moved to BASELINE AT CHANCE
    (nothing)             neither, so there is nothing to report

**An earlier version of this block said BASELINE AT CHANCE "no longer holds
anywhere".** It does hold, for SARS-CoV-2 -- and the artifact diff proved
the claim was already false when it was written rather than overtaken by a
later change: regenerating across the 2026-10-16 `NSP_PATTERNS` fix moved
SARS-CoV and nothing else, leaving that entire block byte-identical. The
sentence was not stale -- it was never right.

**Where it holds has since moved, which is the better argument for printing
it.** When this was written the condition held for SARS-CoV-2 under the
temporal `computed_only` policy -- degree lift 1.15 (p=0.4262) against
COMBINED's 3.31 (p=0.0162) -- and nowhere cross-sectionally. On the
16-virus release it is the reverse. `computed_only` now reads NOT
DETERMINED, because degree's Hits@100 is `[7, 8]`: a 7-compound tie block
straddles the cut with one positive inside it, so its lift 1.08 and
p=0.4693 are not measurements and nothing can be compared against them.
Cross-sectionally the condition now holds, carried by M7. The verdict did
not strengthen or weaken; it changed which run it describes, which no
reader could have inferred from a table that prints the verdict and not
the channel it refers to.

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
2.37x a degree lift of **0.63**, which is M4 lift 1.49 at p=0.3266, and the
cell carries `vs_degree_readable: False` rather than the bare ratio. A ratio
against a noisy denominator describes; it does not test.

Note SARS-CoV especially. M1 there reads p<0.0001 and looks like the
strongest result in the project -- and it is not readable at all:
`hits_determined` is `False`, its Hits@100 is one draw from [2, 23] of 23
positives, so its lift is one draw from [0.21, 2.39], an interval that
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
from [0.21, 2.39]; quoting either one, and then saying in the next clause
that the interval contains 1.00, is the document disagreeing with itself
in two sentences.

**On the 16-virus release the degree row beside it is no longer the
problem.** It used to be undetermined as well ([7, 8] from a tie block), so
the virus read BASELINE NOT DETERMINED and no ratio against it could be
quoted at all. It is now determined at 7 hits, lift 11.69, p=2.4e-06, and
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
measured, which is the 92.2% shared-PMID confound wearing another costume.
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
here, where cross-sectionally it was 0.584 and significant in four of four
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
| taxa | SARS-CoV-2 | 3,460 | 333 | 308.24 | 0.0560 | 0.9500 | — |
| taxa | DENV-2 | 50 | 2 | 5.66 | 0.9825 | 0.0655 | — |
| taxa | West Nile | 18 | 0 | 2.31 | 1.0000 | 0.0940 | — |
| taxa | HCoV-NL63 | 57 | 8 | 6.68 | 0.3473 | 0.7801 | — |
| taxa | Zika | 815 | 65 | 67.85 | 0.6697 | 0.3798 | — |
| taxa | HCoV-OC43 | 657 | 53 | 55.17 | 0.6442 | 0.4118 | — |
| taxa | MERS-CoV | 678 | 54 | 56.08 | 0.6357 | 0.4173 | — |
| taxa | HCoV-229E | 1,367 | 113 | 114.80 | 0.5897 | 0.4428 | — |
| interactome | West Nile | 6 | 4 | 0.86 | 0.0035 | 0.9995 | — |
| interactome | HCoV-229E | 397 | 30 | 45.74 | 0.9965 | 0.0050 | — |
| interactome | HCoV-HKU1 | 279 | 22 | 33.73 | 0.9940 | 0.0115 | — |
| interactome | MERS-CoV | 622 | 57 | 72.48 | 0.9820 | 0.0255 | — |
| interactome | HCoV-OC43 | 289 | 25 | 35.19 | 0.9715 | 0.0440 | — |
| interactome | SARS-CoV | 1,238 | 152 | 133.91 | 0.0480 | 0.9645 | — |
| interactome | DENV-1 | 53 | 11 | 6.57 | 0.0600 | 0.9715 | — |
| interactome | HCoV-NL63 | 311 | 29 | 37.21 | 0.9455 | 0.0800 | — |
| interactome | SARS-CoV-2 | 3,494 | 401 | 381.98 | 0.1114 | 0.9030 | — |
| interactome | DENV-2 | 271 | 38 | 32.19 | 0.1499 | 0.8836 | — |
| interactome | TBEV | 35 | 2 | 4.22 | 0.9330 | 0.1919 | — |
| interactome | DENV-3 | 33 | 6 | 4.05 | 0.2139 | 0.9015 | — |
| interactome | Zika | 353 | 45 | 40.88 | 0.2574 | 0.7886 | — |
| interactome | JEV | 7 | 1 | 0.76 | 0.5562 | 0.8366 | — |

**0 of 22** rows clear the threshold, in either direction.

**The uncorrected view is where the discoveries live.** DENV-1's interactome
overlap is `0.0029 naive` and `0.0600 matched`; DENV-2's is `0.0005 naive`
and `0.1499 matched`. Neither number moved because the data changed. The
naive p-value asks whether these proteins are over-represented in pathways;
the matched one asks whether they are over-represented *compared with
proteins as well-connected as they are*. A twentieth-of-a-percent result
becomes nothing, twice. That gap is the finding of this section, and it is
the same confound the rest of this document measures on compounds.

**The two closest rows, and why neither is a result.** West Nile's
interactome is the smallest p in the table at 0.0035 — three times the
threshold, and it rests on **six** interaction partners, four of which land
in a pathway where 0.86 were expected. Six. The agent attaches a thin-route
caveat to exactly this virus for exactly this reason. On the other side,
HCoV-229E reads depleted at 0.0050, four times the threshold.

**A pattern worth naming and not claiming.** On the interactome side the five
coronaviruses with data all lean depleted: HCoV-229E 0.0050, HCoV-HKU1
0.0115, MERS-CoV 0.0255, HCoV-OC43 0.0440, HCoV-NL63 0.0800. Same direction,
five related viruses, and four of the five would be significant uncorrected.
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
| Evidence co-reported with labels | 92.2% shared PMIDs (933/1,012); a further 50 compounds have no `INHIBITS` publication and are outside the control |
| M1 was never purely direct-acting | 1,884 of 7,677 protein-target rows (24.5%) described a CELLULAR readout and were filed against a viral protein; now reattributed to the organism |
| `INHIBITS` is effectively two proteins | nsp5 and nsp3 dominate; `INHIBITS` reaches 18 distinct viral proteins in all |
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
release and were wrong by the time anyone re-ran the diagnostic:

| was | is | why it moved |
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

# 2. per-virus baselines -- one file each, all SIXTEEN taxa in
#    config/viruses.yaml. It was eight until 2026-10-23, when the remaining
#    flaviviruses were promoted. Six of the eight added here have no
#    HOST_FACTOR_FOR edge, so their host-mediated channels are empty by
#    construction; the file is still written, and an empty channel recorded
#    is what lets a later ingest show up as a change.
for t in 2697049 694009 1335626 11137 277944 31631 290028 64320 \
         11060 11053 11069 11070 11082 11072 11089 11084; do
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

# 4. the pathway landscape. It depends on nothing above -- it reads the
#    release directly and scores no compounds -- but check-results globs
#    data/results and audits it like everything else, so leaving it out
#    meant a reader could run this block top to bottom and still be told
#    a file was stale. It was missing here from the day the file was added
#    until 2026-10-28; it looked fine only because it was being rebuilt by
#    hand. tests/test_reproduce_block.py now fails on a results file that
#    no step here regenerates.
make pathway-landscape

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
