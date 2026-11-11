# Pre-registration: a structure-derived scorer, and the leak it has to survive

**SIGNED 2026-11-11, BEFORE THE CODE AND BEFORE THE RUN.** No score has been
computed, no AUC exists, and `scripts/diagnose_structure_scorer_design.py`
was written so that the facts below could be gathered without producing one.

**What is distinctive here is narrow, and the first draft of this paragraph
overclaimed it.** It said this was the first registration in this directory
written before its measurement; it is not. `prereg_m1_potency.md`,
`prereg_partial_degree.md`, `prereg_phase4_hcv.md` and
`prereg_temporal_auc.md` all open by saying they were written before the
code and committed before the measurement, and nothing here contradicts
them. The one that was signed second is `prereg_pathway_overlap.md`, which
records that at length. What is true of this document and of no other is
smaller: **the scorer it registers does not exist in any form**, so there is
no exploratory version of this number anywhere in the repository, and the
diagnostic that produced the design facts below was built with a boundary it
cannot cross -- `tests/test_prereg_structure_scorer.py` fails if it ever
computes a similarity.

Date: 2026-11-11. Release: v0.1 (78,789 nodes / 392,375 edges, `nodes_sha`
c800b9005c8b074c, `edges_sha` ad7121bc2195a93a). Scripts to be written:
`scripts/structure_scorer.py`. Design facts below: `make
diagnose-structure-design`.

## When the code and the run happened

Recorded here so the ordering above can be checked rather than taken, and
appended rather than replacing anything: a registration edited to match what
happened later is worth nothing.

| | when |
|---|---|
| this document signed | 2026-11-11, commit `9a7b4ee` |
| `scripts/structure_scorer.py` was written on | 2026-11-11, after the signature |
| first run | recorded below, under *Outcome* |

The design facts in this document predate the scorer by construction:
`scripts/diagnose_structure_scorer_design.py` went in with the signature and
a test fails if it ever computes a similarity, so nothing resembling the
outcome was available when the parameters were fixed.

## Why this measurement exists

Every evaluated channel in this project has almost no room to rank. Across
14 viruses and 56 evaluated channel cells the **largest margin above the
imputed floor is +0.0050 and the median is +0.0000**; metapath spans
(`auc_ceiling` - `auc_floor`) run 0.0001 to 0.0377 against degree's 0.5000.
That is not a ranking failure, it is a coverage failure: **74.8% of labelled
compounds carry no traversable edge at all**, so an imputed AUC is almost
entirely `auc_floor`, the free credit a coverage gap pays, and the ranking
contributes a rounding error.

A structure-derived scorer is the only candidate measured to escape that.
**Every labelled compound carries a SMILES: 13,417 of 13,417 under the
ChEMBL-only population registered below, and 18,998 of 18,998 pooled.**
Both are 100.0%, and the registered figure is the first -- quoting the
pooled one for a ChEMBL-only study is the cross-population slip this
document made in its own first draft. So the coverage gap is zero, the floor
is exactly 0.5000, the ceiling exactly 1.0000, and the span is 0.5000: the
same span as degree, for the first time in this project. Whatever it reads
will be a statement about ranking.

**AND THAT IS PRECISELY WHY IT NEEDS REGISTERING RATHER THAN TRYING.** A
scorer with a real span can produce a real-looking number for reasons that
have nothing to do with antiviral activity, and this label set offers it
three. They are named and controlled below, and the controls are fixed here
because each one lowers the headline figure.

## The measurement

For each (virus, compound): the **maximum Tanimoto similarity from the
compound to that virus's known actives**, under a reference set from which
the compound itself and two classes of near-duplicate have been removed.
Candidates are ranked by that score and scored against the measured
negatives.

### Fixed in advance

- **The fingerprint: ECFP4, Morgan radius 2, 2048 bits.** Not chosen here --
  imported from `kgav.similarity` (`FP_RADIUS`, `FP_BITS`), which is what
  the `CHEMICALLY_SIMILAR_TO` layer already uses. Registering a second
  fingerprint would make the scorer and the graph's own similarity edges
  answer the same question two ways, which this project has twice found is
  not corroboration.

- **The statistic: the MAXIMUM over the reference set**, i.e. nearest
  neighbour. Not mean, not mean-of-top-k, not a sum. Max is the standard
  ligand-based virtual-screening baseline and the only one with no free
  parameter; every alternative adds a k or a weighting that would then need
  its own registration. **A single alternative is permitted and is
  SECONDARY: the mean of the top 3.** It is registered now so that it cannot
  be introduced later as the version that worked.

- **The reference set, and the three exclusions.** For candidate *c* and
  virus *v*, the reference set is the positives of *v* minus:

      1. c itself                      -- a positive otherwise scores 1.0
      2. every id sharing c's `inchikey_skel`
                                       -- its own stereoisomers and
                                          tautomers, which the identity
                                          policy groups but does not merge
      3. every id citing a publication c also cites
                                       -- the congeneric-series leak

  Exclusion 3 is the one that matters and the mechanism already exists:
  `baselines.label_publication_index` feeds
  `Graph.load(drop_pubs_shared_with=...)`, and its docstring records that
  **1,015 of 1,070 compounds -- 94.9% -- share at least one PMID between
  their label and the edge scored against it.** The same hazard applies with
  more force here, because a max-similarity scorer's nearest neighbour is
  very often the analog published in the same table.

- **A candidate with no PMID excludes nothing under (3), is KEPT, and is
  COUNTED.** The alternative -- dropping compounds whose reference set
  cannot be made publication-disjoint -- is a restriction on an annotation
  that correlates with the label, which is exactly what fix-list B3 measured
  and recommended declining: PMID coverage is 184 of 213 positives against
  18 of 95 negatives for DENV-3 (86% against 19%) and 1,423 of 1,786 against
  2,278 of 5,362 for SARS-CoV-2 (80% against 42%), so the size and the sign
  of the induced bias depend on the virus. **So the population stays fixed
  and the control's incompleteness becomes a reported per-cell number**
  (`n_candidates_without_publication`), not a filter. This follows E6 and
  L10: record the comparison, keep the channel.

- **The population: ChEMBL-only, per virus, never pooled.** This is the
  hazard that would otherwise dominate. `labels.build_labels`' own docstring
  states it: pooling ChEMBL with the NCATS panel mixes *antiviral research
  compounds* against *library compounds*, "and a scorer can separate the two
  classes by recognising which collection a compound came from". A
  structure-derived scorer is the scorer that does exactly that, so the
  pooled task is **not registered as any part of this analysis** -- not as
  primary, not as sensitivity. `build_labels(sources={"infores:chembl"})`
  keeps both classes within one knowledge source and leaves **14 of 15 cells
  readable**.

- **And the NCATS panel alone is the SECONDARY population, because it is
  matched by construction rather than by a filter.** One screen, one
  protocol, one concentration range, both classes off the same plates:
  `build_labels(sources={"infores:ncats-opendata"})` gives **one cell,
  SARS-CoV-2, 89 positives against 7,343 negatives**. It is the strongest
  control available and the weakest in power, which is why it is secondary
  and reported whatever it says.

- **The metric: the tie-safe rank-sum AUC already in
  `labels.evaluate_against_negatives`**, ties counted as half. Registered
  explicitly because max-Tanimoto ties heavily -- many candidates share an
  identical nearest-neighbour value -- and `prereg_temporal_auc.md` exists
  because a tie block straddling *k* decided a figure by dict insertion
  order. Reported with the Hanley & McNeil interval from
  `labels.auc_interval`.

- **The decomposition is reported for every cell**, `auc_floor` /
  `auc_ceiling` / `auc_skill` from `kgav.labels`, even though the floor is
  predicted to be exactly 0.5000 everywhere. A floor that is not 0.5000 is a
  coverage gap the premise of this study says cannot exist, so printing it is
  the check on the premise.

- **Correction family: the 14 readable ChEMBL cells, two-sided.** alpha =
  0.05/14 = 0.003571, **z = 2.9137**. The family is the 14 cells of the
  primary analysis; the secondary NCATS cell and the top-3-mean variant are
  reported at the same corrected level but are NOT added to the family,
  because adding a planned secondary to the denominator penalises the
  primary for a measurement made to be conservative. Two-sided because a
  below-chance AUC is a finding here, not a null -- it would mean known
  actives are structurally *unlike* the actives.

- **The readability gate: `calibration.MIN_REACHED_POSITIVES` = 10 of each
  class,** unchanged, and every cell is reported beside the smallest AUC its
  class sizes could distinguish from chance. That figure, uncorrected at
  z=1.96 from this project's own interval, is:

      SARS-CoV-2     (2697049)    1,786 / 5,362    0.5156
      DENV-2         (11060)        996 / 3,043    0.5208
      Zika           (64320)        748 / 2,577    0.5237
      West Nile      (11082)        261 /   311    0.5475
      Yellow fever   (11089)        135 /   882    0.5533
      DENV-1         (11053)        195 /   123    0.5640
      DENV-3         (11069)        213 /    95    0.5680
      DENV-4         (11070)        153 /   109    0.5697
      SARS-CoV       (694009)       143 /    96    0.5732
      HCoV-229E      (11137)        120 /    88    0.5778
      HCoV-OC43      (31631)        164 /    70    0.5781
      TBEV           (11084)         72 /    65    0.5948
      MERS-CoV       (1335626)      133 /    41    0.5957
      JEV            (11072)         19 /    27    0.6632

  **HCoV-NL63 (277944) is NOT readable and is not in the family** (5
  positives against 7; it needs AUC >= 0.7842 to clear chance uncorrected).
  Reported as excluded, with its counts, because a cell dropped silently is
  the "powerless row" this project has already had to go back and label
  twice.

  **THE NAMES IN THAT TABLE ARE OFF THE REGISTER AND THE FIRST DRAFT'S WERE
  NOT.** Written by hand it called 11082 a flavivirus (West Nile), 11072 an
  enterovirus (JEV), 11084 a morbillivirus (TBEV), and attributed a PMID
  figure to HCoV-NL63 that belongs to West Nile. The diagnostic now reads
  `virus_register` and prints the name beside every taxon id, so a bare id
  cannot be mislabelled into this document again.

- **The comparator: degree, recomputed on the identical population in the
  same run.** The existing degree figures -- 0.5767 and 0.6585 on SARS-CoV-2,
  0.7397 on HCoV-229E -- come from other protocols with other populations,
  and quoting one as this study's baseline would be the cross-protocol
  comparison `vs_degree_readable` exists to prevent. Degree gets the same
  labels, the same tie-safe AUC and the same interval. **It is the only
  comparison with equal span**, 0.5000 against 0.5000, which is the whole
  reason this scorer is interesting.

## The hazards, quantified before signing

Measured by `make diagnose-structure-design` on `edges_sha
ad7121bc2195a93a`, with no score computed:

- **Near-duplicate chemistry inside the label sets is LOW by exact skeleton
  identity.** On the registered ChEMBL-only population, positives run 0% to
  7% of members sitting in a shared-`inchikey_skel` group and negatives 0% to
  14%, with no group larger than 3. The extremes are HCoV-229E negatives (12
  of 88, 14%), SARS-CoV-2 (118 of 1,786 positives and 385 of 5,362 negatives,
  7% each) and SARS-CoV positives (10 of 143, 7%); HCoV-OC43, Yellow fever
  and JEV have none at all. So exclusion (2) will remove little, and its
  value is that it cannot be blamed afterwards.

- **AND THAT MEASUREMENT UNDERSTATES THE HAZARD ON PURPOSE.** The
  informative version -- the distribution of maximum Tanimoto from each
  class to the actives -- is the leave-one-out score of the positives and
  the score of the negatives, and the separation between the two
  distributions *is* this study's result to four decimal places. It was not
  taken. Exact skeleton identity is a cruder proxy that cannot be turned
  into a ranking, and understating a hazard before signing is the correct
  direction to be wrong in.

- **A publication-disjoint split is buildable in every readable cell.**
  PMIDs reach 1,423 of 1,786 SARS-CoV-2 positives across 354 papers, and 164
  of 164 HCoV-OC43 positives across 182. The weakest are **West Nile (104 of
  261 positives)** and **Zika (331 of 748)**, and under the rule above those
  cells report a larger `n_candidates_without_publication` rather than a
  smaller population. Negative-side coverage is the more uneven half -- 18 of
  95 for DENV-3 against 88 of 88 for HCoV-229E -- which is why the rule
  applies per candidate rather than per cell.

## Registered predictions

**P1. The publication-disjoint scorer reads LOWER than the same scorer with
exclusion (3) removed, in at least 8 of the 14 cells.** Both variants are
computed in one run; the no-exclusion variant exists only to measure the
leak and is reported as a leak estimate, never as a result. If P1 fails --
if removing same-paper actives does not lower the score -- then the
congeneric-series leak is not what dominates this scorer and the rest of the
reading changes.

**P2. The primary scorer clears the corrected alpha in at least one cell.**
Zero metapath channels do, so any cell at all is a new kind of result in
this project. Stated as the weakest interesting claim rather than a
confident one, because nothing here has been run.

**P3. `auc_floor` is exactly 0.5000 and `auc_ceiling` exactly 1.0000 in
every one of the 14 cells.** A violation means a labelled compound failed to
fingerprint, which the premise of this study -- 100.0% SMILES coverage --
says cannot happen but RDKit can still cause. This is the prediction most
likely to fail and the cheapest to check.

**P4. SARS-CoV-2 and DENV-2, the two best-powered cells, clear the corrected
alpha if any cell does; JEV (19 positives against 27) does not clear it
under any variant.** It needs 0.6632 uncorrected and more corrected.

**P5. The NCATS single-screen cell reads LOWER than the ChEMBL SARS-CoV-2
cell.** This is the load-bearing prediction. The panel is matched by
construction, so the collection confound that the ChEMBL cell can only
partly control is absent there. If the matched cell reads as high as or
higher than the unmatched one, the confound is not driving the result and
the finding is stronger than this registration expects.

**P6. Degree does not clear the corrected alpha in a cell where the
structure scorer does not.** Registered because the inverse has never been
checked: degree's published figures were computed on other populations, and
this is the first time it is measured on these.

## What would change the conclusion

A cell clearing the corrected alpha under the publication-disjoint primary
**and** retaining a positive margin in the NCATS matched cell, with
`auc_skill` above degree's on the same population. That would be the first
scorer in this project with both a span and a margin, and it would justify
the held-out temporal evaluation this registration does not attempt.

## What this cannot do

- **It is not a mechanism and cannot become one.** A compound scoring well
  because it resembles a known active says nothing about what it binds. The
  graph channels at least name a path; this names a neighbour.
- **It does not remove the collection confound, only the worst of it.**
  ChEMBL-only matches the source, not the assay, the cell line or the
  concentration range. Only the NCATS cell matches those, and it is one cell
  with 89 positives.
- **It does not test the temporal question.** No held-out-by-year evaluation
  is registered here, because the reference set under a temporal split is
  the pre-2021 actives and that is a different and smaller study. It is the
  obvious follow-up and must be registered separately, before it runs.
- **It does not settle `MIN_REACHED_POSITIVES` or `MIN_CLASS_SIZE`.** Both
  are 10 and the pair is circular, which `prereg_pathway_overlap.md` already
  names. This analysis uses the gate as it stands and does not relitigate it.
- **It cannot rescue the metapath channels.** If this scorer works, the
  graph channels are still coverage-limited; a structural scorer that ranks
  well is a reason to fix coverage, not a substitute for it.

## Implementation constraints

No new dependency: RDKit is already required by `kgav.similarity`, and
`DataStructs.BulkTanimotoSimilarity` is already the call used. The largest
cell is 13,417 candidates against 1,786 references -- 24.0M similarities,
about 96 MB as float32 -- on a machine with 242 GB available and 32 cores,
and the per-virus work is independent. So the full matrix is computed once
per virus and masked per candidate rather than recomputed, and the
exclusions are set operations on indices.

**The run must print the fingerprint parse rate**, because P3 depends on it
and a SMILES is not a fingerprint. A compound whose SMILES RDKit rejects is
unreachable for this scorer and must appear in `auc_floor` rather than
disappear.

## What is still the owner's to change -- but only now

Every value above is fixed and this document is signed. **Changing any of
them is legitimate until the first score exists and indefensible
afterwards**, which is the asymmetry `prereg_pathway_overlap.md` lost by
being signed second. The three most contestable, with the reasoning already
given: max rather than mean-of-top-3 as primary; ChEMBL-only rather than
NCATS-only as primary; and keeping no-PMID candidates rather than dropping
them. A change to any of these before the run is a re-registration and
should be recorded as one.

---

# Verdict — the scorer has a span, and the control that gives it one also decides its sign

Run 2026-11-11, `make structure-scorer`, written to
`data/results/structure_scorer.json`. Release v0.1, `edges_sha`
ad7121bc2195a93a. The registration above is unchanged; this section is
appended beneath it.

| | prediction | outcome |
|---|---|---|
| **P1** | publication-disjoint reads lower in >= 8 of 14 cells | **HELD, 13 of 14** |
| **P2** | at least one cell clears the corrected alpha | **HELD, 5 cells above chance** |
| **P3** | floor exactly 0.5000, ceiling exactly 1.0000, every cell | **HELD exactly, all 14** |
| **P4** | SARS-CoV-2 and DENV-2 clear if any does; JEV clears under no variant | **HALF-FAILED** — both clear, and JEV clears too |
| **P5** | the matched NCATS cell reads lower than the ChEMBL cell | **CANNOT BE EVALUATED** |
| **P6** | degree does not clear where the structure scorer does not | **FAILED, 4 cells** |

## The numbers

    virus                 pos     neg      AUC   corrected 95%      skill    leak   degree
    SARS-CoV-2          1,786   5,362   0.7693^  [0.7488, 0.7897]   0.5385  +0.1003  0.5131
    MERS-CoV              133      41   0.6677^  [0.5364, 0.7990]   0.3354  -0.0001  0.7085
    Zika                  748   2,577   0.6316^  [0.5965, 0.6666]   0.2631  +0.1722  0.5534
    West Nile             261     311   0.6308^  [0.5626, 0.6990]   0.2616  +0.1123  0.5657
    HCoV-229E             120      88   0.5951   [0.4806, 0.7097]   0.1902  +0.0386  0.5908
    TBEV                   72      65   0.5764   [0.4344, 0.7184]   0.1528  +0.1840  0.4203
    DENV-2                996   3,043   0.5735^  [0.5425, 0.6046]   0.1471  +0.2899  0.4016
    Yellow fever          135     882   0.5545   [0.4752, 0.6337]   0.1089  +0.2598  0.7264
    HCoV-OC43             164      70   0.5396   [0.4210, 0.6582]   0.0792  +0.0992  0.5345
    SARS-CoV              143      96   0.5186   [0.4077, 0.6294]   0.0372  +0.0998  0.6812
    DENV-4                153     109   0.4644   [0.3587, 0.5701]  -0.0712  +0.3036  0.2302
    DENV-1                195     123   0.4489   [0.3518, 0.5461]  -0.1021  +0.2619  0.3144
    DENV-3                213      95   0.3861v  [0.2825, 0.4896]  -0.2279  +0.3848  0.1715
    JEV                    19      27   0.1452v  [0.0000, 0.3061]  -0.7096  +0.6832  0.6374

    ^ clears the corrected alpha above chance   v clears it below

**THE SPAN IS REAL AND P3 IS THE REASON TO BELIEVE ANY OF THIS.** 19,085 of
19,085 labelled compounds fingerprinted, zero unparseable, zero without
SMILES; `coverage_pos` and `coverage_neg` are 1.0000 in all 14 cells, so
`auc_floor` is exactly 0.5000 and `auc_ceiling` exactly 1.0000 everywhere.
Every thousandth of SARS-CoV-2's 0.7693 is ranking. For comparison, the
best metapath channel in this project has a span of 0.0377 and a margin
above its floor of +0.0050.

## P5 could not be evaluated, and the reason is a defect in the registration

The matched NCATS cell **refuses rather than reporting**. All 7,432
candidates lost their entire reference set: the panel's 89 actives and
7,343 inactives all cite publications from one screen, so exclusion (3)
removed every admissible reference for every candidate.
`excluded_publication` fired 661,358 times, `empty_reference` is 7,432, and
`coverage_pos` is 0.0.

**The first run of this script reported that cell as `0.5000 [0.4103,
0.5897]`** -- an AUC with an interval, reading exactly like a scorer that
found nothing, when in fact nothing had been ranked and a rank-sum AUC over
a single tie block is 0.5 by arithmetic. `auc_skill` was None and the
printed line showed `0.0000`, because the formatter wrote `or 0.0`. Both are
fixed: a cell with zero coverage now carries `measurable: false` and a
refusal string, and the formatter prints `--` for a quantity that does not
exist.

**The registration's two strongest controls are incompatible by
construction, and nothing in the design facts would have shown it.** A
single-screen population is matched because one paper reports it; a
publication-disjoint reference set requires the actives to come from
different papers. `make diagnose-structure-design` reported PMID coverage of
89/89 and 7,343/7,343 across 232 papers for that cell, which looked like
ample material -- the question it does not answer is whether the ACTIVES'
papers are disjoint from each other, and they are not. That question is now
the one worth asking before any future registration pairs the two.

## P6 failed, and P4 failed for a reason I registered and then did not think through

**Degree clears the corrected alpha in 4 cells where the structure scorer
does not** -- DENV-1, DENV-4, Yellow fever and SARS-CoV -- so P6 is false.
Two of those four clear by ranking actives BELOW inactives (DENV-1 0.3144,
DENV-4 0.2302), which is the same two-sided reading that broke P4.

**P4 said JEV "does not clear it under any variant", reasoning that it needs
0.6632 uncorrected.** It clears at **0.1452**, by anti-ranking, with the
interval [0.0000, 0.3061] excluding 0.5 from below. The registration fixed a
two-sided test on purpose -- *"a below-chance AUC is a finding here, not a
null"* -- and then the prediction was written as though only the upper tail
existed. The test was right and my reasoning about it was not.

## What the leak measurement says, which is the most useful thing here

**P1 held in 13 of 14 cells and the magnitudes are large:** removing
same-paper actives costs +0.6832 AUC on JEV, +0.3848 on DENV-3, +0.3036 on
DENV-4, +0.2899 on DENV-2. Only MERS-CoV is unaffected (-0.0001).

**And the leak predicts the result.** Across the 14 cells, Spearman
correlation between the final AUC and the leak magnitude is **-0.684**,
which clears the ~0.544 critical value at n=14: the cells that lose the most
to the control are the cells that end up lowest. Read plainly, most of what
a naive run of this scorer would have reported is same-paper retrieval.

**A SECOND MECHANISM IS SUGGESTIVE AND DOES NOT CLEAR THE BAR, AND IT IS A
FLAW IN THE CONTROL ITSELF.** Exclusion (3) is not symmetric between the
classes: a positive's nearest known actives are its own paper-mates, which
are removed, while a negative sitting in a paper containing no active loses
nothing. References excluded per candidate, positives against negatives:

    JEV            11.16 / 2.00   = 5.6x      AUC 0.1452
    DENV-3         23.66 / 1.94   = 12.2x     AUC 0.3861
    DENV-1         13.68 / 2.75   = 5.0x      AUC 0.4489
    DENV-4         11.93 / 1.98   = 6.0x      AUC 0.4644
    ...
    MERS-CoV       10.36 / 12.61  = 0.8x      AUC 0.6677
    SARS-CoV-2     20.63 / 32.23  = 0.6x      AUC 0.7693

The two cells where the control handicaps NEGATIVES more are the two
highest AUCs, and the four where it handicaps positives by 5x or more are
the four lowest. Spearman against log ratio is **-0.490**, which does NOT
clear 0.544 at n=14, so this is a hypothesis and is reported as one. It is
also the obvious thing to fix: a control that removes a candidate's closest
analogs should remove them for both classes on the same terms, and this one
removes them in proportion to how much company each compound keeps.

## What this does and does not license

- **It licenses the span claim and nothing about mechanism.** A scorer that
  reaches every labelled compound and separates the classes at 0.7693 on
  SARS-CoV-2 is the first thing in this project with room to rank and a
  margin to show. It still only names a neighbour.
- **It does not license the SARS-CoV-2 figure as an estimate of
  prospective performance.** 3,447 of its 7,148 candidates -- 48% -- carry
  no PMID, so the leak control did not apply to half of them; the registered
  rule kept them and counted them, and the count is the caveat.
- **It does not establish that the collection confound was controlled.** The
  one cell that could have tested it is the one that refuses. ChEMBL-only
  matches the source, not the assay.
- **The four below-chance cells are not a finding about chemistry** until
  the asymmetry above is ruled out. They may be the control's artifact.

## What would change the conclusion, restated for the next registration

A control that handicaps both classes equally -- for example, excluding each
candidate's k nearest references regardless of publication, with k fixed in
advance -- run on the same 14 cells. If SARS-CoV-2 survives that at a
comparable margin, the result is about structure. If it collapses towards
the dengue serotypes, it was about who publishes with whom. **That must be
registered before it runs, and this verdict is exactly the kind of result
that makes the next registration hard to write honestly.**
