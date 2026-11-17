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

Date: 2026-11-11. Release: v0.1 (78,649 nodes / 395,646 edges, `nodes_sha`
dc6d76a2b3cfb2a1, `edges_sha` c0c4daccd48e9979). Scripts to be written:
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
**Every labelled compound carries a SMILES: 13,345 of 13,345 under the
ChEMBL-only population registered below, and 18,900 of 18,900 pooled.**
Both are 100.0%, and the registered figure is the first -- quoting the
pooled one for a ChEMBL-only study is the cross-population slip this
document made in its own first draft. So the coverage gap is zero, the floor
is exactly 0.5000, the ceiling exactly 1.0000, and the span is 0.5000: the
same span as degree, for the first time in this project. Whatever it reads
will be a statement about ranking.

**[The pooled figure above is the pre-entry-assay release's. Superseded for
the current release by the amendment of 2026-11-17 at the end of this
document; the ChEMBL-only figure is unmoved.]**

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

      SARS-CoV-2    (2697049)    1,761 / 5,320    0.5157
      DENV-2        (11060)        995 / 3,043    0.5208
      Zika          (64320)        748 / 2,577    0.5237
      West Nile     (11082)        261 /   311    0.5475
      Yellow fever  (11089)        135 /   882    0.5533
      DENV-1        (11053)        195 /   123    0.5640
      DENV-3        (11069)        213 /    95    0.5680
      DENV-4        (11070)        153 /   109    0.5697
      SARS-CoV      (694009)       142 /    96    0.5733
      HCoV-229E     (11137)        119 /    88    0.5780
      HCoV-OC43     (31631)        162 /    70    0.5783
      TBEV          (11084)         72 /    65    0.5948
      MERS-CoV      (1335626)      131 /    41    0.5959
      JEV           (11072)         19 /    27    0.6632

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
545fb7d87b0aec08`, with no score computed:

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
545fb7d87b0aec08. The registration above is unchanged; this section is
appended beneath it.

| | prediction | outcome |
|---|---|---|
| **P1** | publication-disjoint reads lower in >= 8 of 14 cells | **HELD, 13 of 14** |
| **P2** | at least one cell clears the corrected alpha | **HELD, 5 cells above chance** |
| **P3** | floor exactly 0.5000, ceiling exactly 1.0000, every cell | **HELD exactly, all 14** |
| **P4** | SARS-CoV-2 and DENV-2 clear if any does; JEV clears under no variant | **HALF-FAILED** — both clear, and JEV clears too |
| **P5** | the matched NCATS cell reads lower than the ChEMBL cell | **CANNOT BE EVALUATED** |
| **P6** | degree does not clear where the structure scorer does not | **FAILED, 3 cells** (4 before the comparator was fixed; see the amendment) |

## The numbers

    virus                 pos     neg      AUC   corrected 95%      skill    leak   degree
    SARS-CoV-2          1,761   5,320   0.7718^  [0.7512, 0.7924]   0.5436  +0.1021  0.5151
    MERS-CoV              131      41   0.6851^  [0.5566, 0.8135]   0.3701  -0.0067  0.7042
    Zika                  748   2,577   0.6316^  [0.5965, 0.6666]   0.2631  +0.1722  0.5534
    West Nile             261     311   0.6308^  [0.5626, 0.6990]   0.2616  +0.1123  0.5662
    HCoV-229E             119      88   0.5917   [0.4767, 0.7067]   0.1834  +0.0345  0.5879
    TBEV                   72      65   0.5764   [0.4344, 0.7184]   0.1528  +0.1840  0.4203
    DENV-2                995   3,043   0.5742^  [0.5431, 0.6053]   0.1484  +0.2901  0.4013
    Yellow fever          135     882   0.5545   [0.4752, 0.6337]   0.1089  +0.2598  0.7264
    HCoV-OC43             162      70   0.5257   [0.4062, 0.6453]   0.0515  +0.0998  0.5363
    SARS-CoV              142      96   0.5098   [0.3986, 0.6210]   0.0197  +0.1020  0.6826
    DENV-4                153     109   0.4644   [0.3587, 0.5701]  -0.0712  +0.3036  0.2302
    DENV-1                195     123   0.4489   [0.3518, 0.5461]  -0.1021  +0.2619  0.3144
    DENV-3                213      95   0.3861v  [0.2825, 0.4896]  -0.2279  +0.3848  0.1715
    JEV                    19      27   0.1452v  [0.0000, 0.3061]  -0.7096  +0.6832  0.6374

    ^ clears the corrected alpha above chance   v clears it below

**THE SPAN IS REAL AND P3 IS THE REASON TO BELIEVE ANY OF THIS.** 19,168 of
19,168 labelled compounds fingerprinted, zero unparseable, zero without
SMILES; `coverage_pos` and `coverage_neg` are 1.0000 in all 14 cells, so
`auc_floor` is exactly 0.5000 and `auc_ceiling` exactly 1.0000 everywhere.
Every thousandth of SARS-CoV-2's 0.7718 is ranking. For comparison, the
best metapath channel in this project has a span of 0.0377 and a margin
above its floor of +0.0050.

## P5 could not be evaluated, and the reason is a defect in the registration

The matched NCATS cell **refuses rather than reporting**. All 7,427
candidates lost their entire reference set: the panel's 89 actives and
7,338 inactives all cite publications from one screen, so exclusion (3)
removed every admissible reference for every candidate.
`excluded_publication` fired 660,913 times, `empty_reference` is 7,427, and
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

**Degree clears the corrected alpha in 3 cells where the structure scorer
does not** -- HCoV-229E, HCoV-OC43 and SARS-CoV -- so P6 is false.

**THIS PARAGRAPH SAID "4 cells -- DENV-1, DENV-4, Yellow fever and
SARS-CoV", AND THREE OF THOSE FOUR WERE AN ARTEFACT OF A DEFECTIVE
COMPARATOR.** It went on: *"Two of those four clear by ranking actives
BELOW inactives (DENV-1 0.3144, DENV-4 0.2302), which is the same
two-sided reading that broke P4."* Those two anti-rankings were the
degree baseline counting the label edges it was being scored against --
a ChEMBL negative carries more label edges than a positive, so the null
ranked negatives first. With the held-out predicates excluded, DENV-1
reads 0.5251 and DENV-4 0.5245: at chance, clearing nothing. Yellow
fever went 0.7262 -> 0.5079 the other way. Two cells that could not
clear before now do -- HCoV-229E 0.6434 and HCoV-OC43 0.6222 -- so the
count is 3 and only SARS-CoV survives from the original four. The
defect, and every moved cell, is recorded in the amendment at the foot
of this document. P6 is still false; it is false about different cells
and for a reason that is about study volume rather than about an
arithmetic artefact.

**P4 said JEV "does not clear it under any variant", reasoning that it needs
0.6632 uncorrected.** It clears at **0.1452**, by anti-ranking, with the
interval [0.0000, 0.3061] excluding 0.5 from below. The registration fixed a
two-sided test on purpose -- *"a below-chance AUC is a finding here, not a
null"* -- and then the prediction was written as though only the upper tail
existed. The test was right and my reasoning about it was not.

## What the leak measurement says, which is the most useful thing here

**P1 held in 13 of 14 cells and the magnitudes are large:** removing
same-paper actives costs +0.6832 AUC on JEV, +0.3848 on DENV-3, +0.3036 on
DENV-4, +0.2901 on DENV-2. Only MERS-CoV is unaffected (-0.0067).

**And the leak predicts the result.** Across the 14 cells, Spearman
correlation between the final AUC and the leak magnitude is **-0.723**,
which clears the ~0.544 critical value at n=14: the cells that lose the most
to the control are the cells that end up lowest. Read plainly, most of what
a naive run of this scorer would have reported is same-paper retrieval.

**A SECOND MECHANISM IS SUGGESTIVE AND DOES NOT CLEAR THE BAR, AND IT IS A
FLAW IN THE CONTROL ITSELF.** Exclusion (3) is not symmetric between the
classes: a positive's nearest known actives are its own paper-mates, which
are removed, while a negative sitting in a paper containing no active loses
nothing. References excluded per candidate, positives against negatives:

    JEV            11.16 / 2.00   =   5.6x      AUC 0.1452
    DENV-3         23.66 / 1.94   =  12.2x      AUC 0.3861
    DENV-1         13.68 / 2.75   =   5.0x      AUC 0.4489
    DENV-4         11.93 / 1.98   =   6.0x      AUC 0.4644
    SARS-CoV        9.23 / 6.71   =   1.4x      AUC 0.5098
    HCoV-OC43      25.36 / 26.99  =   0.9x      AUC 0.5257   <- counterexample
    Yellow fever   12.13 / 2.99   =   4.1x      AUC 0.5545
    DENV-2         26.99 / 2.13   =  12.7x      AUC 0.5742   <- counterexample
    TBEV           10.31 / 5.91   =   1.7x      AUC 0.5764
    HCoV-229E      10.97 / 3.49   =   3.1x      AUC 0.5917
    West Nile       3.16 / 2.18   =   1.5x      AUC 0.6308
    Zika            6.89 / 0.84   =   8.2x      AUC 0.6316   <- counterexample
    MERS-CoV       10.95 / 13.02  =   0.8x      AUC 0.6851
    SARS-CoV-2     22.16 / 37.45  =   0.6x      AUC 0.7718

    <- counterexample: a cell the asymmetry hypothesis does not
       explain -- handicapped 5x or more against POSITIVES yet above
       chance, or handicapped against NEGATIVES yet not among the two
       highest cells.

**ALL FOURTEEN ROWS ARE SHOWN AND THE COUNTEREXAMPLES ARE MARKED.** Until
2026-11-11 this table showed four rows, an elision, and two more -- the four
most positive-handicapped low cells and the two negative-handicapped high
ones -- and the paragraph beneath it read *"the two cells where the control
handicaps NEGATIVES more are the two highest AUCs, and the four where it
handicaps positives by 5x or more are the four lowest"*. **Both halves were
false when written, and the elision is what hid it:**

- **Three cells handicap negatives more, not two.** The third is HCoV-OC43
  at 0.94x, AUC 0.5257, which is sixth of fourteen from the bottom.
- **Five cells are handicapped 5x or more against positives**, not four, and
  two of them -- DENV-2 at 12.67x and Zika at 8.20x -- read ABOVE chance.
  Those are the same two counterexamples the symmetric registration's copy
  of this table was corrected for on the same day. This copy was not,
  because no test read it.
- **And the four lowest cells are not the four at 5x or more.** DENV-1 is
  third from the bottom at 4.98x, below the bar; DENV-4 at 6.02x is fourth.
  The sentence named a set of four that does not exist.

**WHAT THE TABLE ACTUALLY SUPPORTS IS WEAKER AND IS STILL WORTH STATING.**
At the extremes the association holds: the two highest cells (SARS-CoV-2 at
0.62x, MERS-CoV at 0.84x) are both handicapped more against negatives, and
the two lowest (JEV at 5.58x, DENV-3 at 12.22x) are both handicapped 5x or
more against positives. In the middle it breaks, in both directions.
Spearman against log ratio is **-0.490**, which does NOT clear 0.544 at
n=14, so this is a hypothesis and is reported as one. It is
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

**IT IS NOW REGISTERED:** `docs/prereg_structure_scorer_symmetric.md`,
signed 2026-11-11 before its code, with k = 10 and a curve over
k ∈ {1, 3, 10}. Writing it immediately found two things this verdict got
wrong. **The k = 30 suggested above is infeasible** -- JEV holds 18
references after self-exclusion, so k = 30 empties the cell and reproduces
the P5 failure recorded here; the class sizes said so without a run. And
**the asymmetry table in this verdict is incomplete in a way that favours
its own hypothesis**: it shows the four highest ratios and the two lowest,
omitting DENV-2 at 12.7x and Zika at 8.2x, which are heavily asymmetric and
ABOVE chance. Those two are counterexamples. The new registration carries
all fourteen rows with the counterexamples marked, and the -0.490
correlation reported here already included them -- the table did not.

---

# Amendment, 2026-11-17: the secondary cell needs a second coordinate to name it

**No AUC has been computed on this release.** Appended; nothing above is
rewritten.

## The expression that defined the secondary cell no longer defines one cell

Signed above: *"`build_labels(sources={"infores:ncats-opendata"})` gives one
cell, SARS-CoV-2, 89 positives against 7,343 negatives."* On the release of
2026-11-17 the same expression gives **three** cells, and the SARS-CoV-2 one
reads **149 / 7,474**.

| virus | by source alone |
|---|---:|
| SARS-CoV-2 | 149 / 7,474 |
| SARS-CoV | 95 / 1,524 |
| MERS-CoV | 71 / 1,951 |

`infores:ncats-opendata` now publishes four screens rather than one
(`docs/prereg_ncats_entry_assays.md`), and 149 is a **pooled** figure: a
cytopathic-effect screen and a pseudotyped-particle entry screen, two
protocols over two concentration ranges, in one cell.

**This is the registration's own hazard, one level down.** The population
decision above refuses to pool ChEMBL with the panel because it "mixes
antiviral research compounds against library compounds". Pooling CPE with
entry mixes *protect-the-cell* against *block-the-spike* -- and it does it
inside the cell that exists specifically to be **matched by construction**.
The sentence "one screen, one protocol, one concentration range, both
classes off the same plates" stopped being true of that expression, while
every count it produces still looks reasonable. Nothing in the old code
could have told the difference.

## The decision — **DECISION, 2026-11-17**

**The secondary cell is `sources={"infores:ncats-opendata"},
assay_types={"cell_based_antiviral"}`: 89 positives against 7,338
negatives.** Source answers whose measurement; `assay_type` answers
measuring what; a cell needs both pinned before "matched by construction" is
true of it.

- **The positives are the registered ones, unchanged and byte-identical** --
  `sha256 f0dfc331cc232729a0100323dd09cf8269b44baf027cda67e440db0140d5da99`,
  the same digest the pre-entry snapshot recorded.
- **The negatives are 7,338, not the registered 7,343.** Minus five,
  accounted for compound by compound in
  `docs/prereg_ncats_entry_assays.md`'s second amendment: four weak-hit
  refusals and one edge lost to a merge. The direction is conservative --
  it removes measured inhibitors from the class of non-inhibitors.
- **The two pseudotyped-entry cells are NOT added to this registration.**
  They are registered in `prereg_ncats_entry_assays.md` with their own
  family of three, and this document's family stays the 14 readable ChEMBL
  cells at alpha = 0.05/14. A second registration's cells do not join a
  first registration's family after the fact.

## P5 still cannot be evaluated, and that was re-measured rather than assumed

The verdict above records P5 refusing because every candidate lost its
entire reference set: one screen, one paper, so publication-disjointness
removes every admissible reference. The merge has since changed the
publications on that cell's edges, so the refusal was re-derived:

| edges | publications |
|---|---|
| 4,805 CPE-typed | `PMID:33708112` only |
| 2,622 CPE-typed | `PMID:33708112` **and** `PMID:35086793` |

**All 7,427 still carry `PMID:33708112`**, and a reference is excluded on
*any* shared PMID, so every reference is still excluded for every candidate
and `coverage_pos` is still 0. **P5's refusal stands, for the reason
originally recorded.** The second PMID does not loosen the control; it was
worth checking because it could only have tightened or loosened it, and a
verdict that assumed either would be guessing.

## Registered in advance for the two new plate cells

Measured before any run: the SARS-CoV entry cell cites exactly one
publication and the MERS-CoV entry cell cites exactly one, and **it is the
same paper, `PMID:33330839`.** Both cells therefore refuse under the
publication control for precisely the reason P5 refused.

**Registered here so it is not reported as a discovery:** all three plate
cells are unreadable under the publication-disjoint control. The only
registered control that can read them is the symmetric k-nearest-reference
exclusion of `docs/prereg_structure_scorer_symmetric.md`, which handicaps
both classes regardless of who published what. A plate AUC quoted under the
publication control is a refusal, not a number.

## The premise figures, re-derived

The premise above states SMILES coverage over two populations. On this
release:

- **ChEMBL-only: 13,345 of 13,345 — unmoved.** The rebuild did not touch the
  registered population, which is the same finding the entry-assay amendment
  records for the 15 ChEMBL label cells.
- **19,091 of 19,091 pooled**, up from the signed 18,900 by the 191
  compounds the entry assays added. Still 100.0%, so the premise holds and
  the floor is still exactly 0.5000.

The signed sentence's pooled figure is **superseded for this release** and
left in place, because a registration edited to match the release it is
being read against is not a registration.

## Filed, not fixed: the citation on a merged edge

The 2,622 CPE-typed edges citing `PMID:35086793` are the visible shadow of
the identity rule. `assemble` unions publications by design -- "one fact
with two witnesses" -- and that is correct *given* that both sides assert
the same fact. They do not: `edge_key` excludes `assay_type`, so a CPE
measurement and an entry measurement of one compound against one virus are
one key, the CPE qualifiers win, and the entry paper is appended as a
witness to a measurement it does not report.

**It is a useful diagnostic and a wrong citation at the same time.** Useful:
2,622 + 196 = 2,818, exactly the SARS-CoV-2 entry measurements the ingest
produced, so the count of two-PMID edges is the count of measurements the
merge discarded. Wrong: a reader following `PMID:35086793` from an edge
labelled `cell_based_antiviral` finds a different assay, and cannot tell
from the edge which paper reports the EC50 it carries.

**Not fixed here.** The fix is to make `assay_type` an identity qualifier,
which is fix-list L2 and splits every merged edge in the graph. No claim in
this document rests on the citation list of an NCATS activity edge, and the
one control that reads publications refuses this cell for an unrelated
reason. Stated so that it is on the record before anything is reported, not
after someone follows a link.

---

# Amendment, 2026-11-17: the comparator was reading the answer sheet

**APPENDED, NOT EDITED.** Nothing above is changed. This records a defect
in the COMPARATOR this registration defined, found after the verdict was
written, and the verdict figures it moves.

Release this amendment was measured on: v0.1, **78,649 nodes 396,223
edges**, `nodes_sha` dc6d76a2b3cfb2a1, `edges_sha` 873623c9a9df2049 —
the release above plus the 577 clinical-label edges of
`docs/prereg_clinical_labels.md`, which is how the defect surfaced. The
line the registration was signed on, 395,646 edges, is kept above as
signed. **No structure figure moves:
161 non-degree channel cells across the three studies are identical
before and after, to machine precision. Every moved number is a `degree`
number.**

## The defect

This registration specified *"the comparator: degree, recomputed on the
identical population in the same run"*. `baselines.degree_ranking` has
stated since fix-list L7 that it ranks over the **traversable** graph,
not the whole one, because *"a null hypothesis that reads the held-out
answers is not a null hypothesis"*.

It delivers that only if its caller asks. `Graph.load` takes
`skip_predicates` and **defaults it to the empty set**, in which case
`degree_unlabelled` is `degree`. **Every other `Graph.load` caller in
`scripts/` that touches the held-out set passes
`schema.held_out_predicates()`** — twelve of them on the day this was
found. `scripts/structure_scorer.py` — the only caller that computes the
comparator for these three registrations — passed nothing, on line 565,
from the day it was written, which `git show` of its first commit
confirms.

The claim above is written as an invariant rather than a count because a
tally of callers goes stale the next time a script is added, and one was
added the same day: `scripts/clinical_labels.py`, which passes it.

So every `degree` figure in `structure_scorer.json`,
`structure_scorer_symmetric.json` and
`structure_scorer_strength_matched.json`, and every verdict sentence
comparing a structure channel against "its own degree baseline", was
measured against a null that counted the `HAS_ANTIVIRAL_ACTIVITY_AGAINST`
and `MEASURED_INACTIVE_AGAINST` edges the studies are scored on.

**IT WAS FOUND BY COINCIDENCE AND THAT IS THE WORST PART.** The
`indications` layer added 577 `TREATS`/`IN_TRIAL_FOR` edges on
2026-11-17. Those predicates are `held_out: true`, so the degree
baseline was documented as unable to see them — and the plate cells'
degree AUC moved anyway, 0.5513 -> 0.5524. A number that could not move
moved, which is the only reason anyone looked. **Twelve hundred tests
passed over this for six days**, because the docstring described
behaviour the default did not provide and no test compared the two.

`tests/test_degree_baseline_skips_labels.py` now holds the general form:
`degree_ranking` REFUSES a graph loaded with nothing skipped, and a
commit-time check fails any script that ranks by degree over such a
graph. It fails against the line as it stood.

## What it cost

All seventeen readable cells moved. **Eleven changed significance
status.** Largest moves, publication-control run (the figures are
identical in all three studies, because `degree` depends on the graph
and not on the exclusion protocol):

    virus           before    after     delta   alpha verdict
    DENV-3          0.1736   0.4970   +0.3234   clears -> does not
    DENV-4          0.2332   0.5245   +0.2913   clears -> does not
    Yellow fever    0.7262   0.5079   -0.2184   clears -> does not
    DENV-1          0.3164   0.5251   +0.2087   clears -> does not
    DENV-2          0.4012   0.5638   +0.1626
    JEV             0.6374   0.5263   -0.1111
    HCoV-OC43       0.5382   0.6222   +0.0840   does not -> clears
    SARS-CoV-2      0.4874   0.5535   +0.0660   does not -> clears
    HCoV-229E       0.5879   0.6434   +0.0556   does not -> clears
    TBEV            0.4203   0.4753   +0.0550
    MERS-CoV        0.7001   0.6533   -0.0467
    West Nile       0.5655   0.5319   -0.0337
    Zika            0.5540   0.5342   -0.0198   clears -> does not
    SARS-CoV        0.6743   0.6908   +0.0165
      plate cells
    SARS-CoV-2 CPE  0.5524   0.5700   +0.0176
    SARS-CoV entry  0.5536   0.6070   +0.0534   does not -> clears
    MERS-CoV entry  0.5839   0.6326   +0.0488   does not -> clears

**THE DIRECTION IS EXPLAINED BY THE DEFECT, which is what says this is a
fix and not a tuning.** Four flavivirus cells had degree AUCs between
0.17 and 0.33 — far below chance, and "clearing" their alpha in the
wrong direction. A ChEMBL negative is a compound measured and found
inactive, often in a large panel, so negatives carry more label edges
than positives; a null counting label edges therefore ranked negatives
above positives and produced a strongly anti-predictive baseline out of
nothing but assay volume. With the label edges out, those four sit
between 0.4970 and 0.5638 — at chance, which is what a null hypothesis
on a structureless feature should look like.

## What it costs THIS registration's verdict

P6 read *"degree does not clear where the structure scorer does not"* and
was recorded **FAILED, 4 cells**. On the corrected comparator the
publication-control run has **degree at or above the structure statistic
in 7 of 17 cells**, and clearing its alpha in three of them — HCoV-229E
0.6434 against 0.5917, HCoV-OC43 0.6222 against 0.5257, SARS-CoV 0.6908
against 0.5098. P6 does not become passable; it fails differently, and
the cells it fails on are not the same cells.

The headline SARS-CoV-2 ChEMBL result is **unmoved and the gap widens in
its favour**: structure 0.7718 publication-disjoint and 0.8304 symmetric,
against a degree baseline that went 0.4874 -> 0.5535. A corrected null
nearer chance makes that comparison cleaner, not weaker.

**AND ONE PLATE VERDICT CHANGES.** Under the symmetric control the three
plate cells were reported as all three clearing with "degree clearing in
none". With the comparator fixed, degree clears in **two of three** and
in the MERS-CoV entry cell it is **above** the structure statistic:
0.6326 against 0.5873. That cell can no longer be read as evidence about
chemistry.

**AND THE STATEMENT IT FALSIFIES IS NOT IN A DOCUMENT, WHICH IS ITS OWN
FINDING.** "Degree clears in none of the three plate cells" was reported
when the plate AUCs were first computed, and neither
`docs/prereg_structure_scorer_symmetric.md` nor
`docs/prereg_structure_scorer_strength_matched.md` contains the word
`degree` at all — nor does the plate family's own registration,
`docs/prereg_ncats_entry_assays.md`, beyond one aside that degree
changes when compounds enter. So the comparator that decides whether any
plate result is about chemistry was quoted in conversation and in
`data/results`, and written into none of the three registrations that
depend on it. `docs/scorer_aucs.txt` is generated from the results files
and prints the degree column beside every structure column for exactly
this reason; it is regenerated by `make scorer-table` as the last step of
`make reproduce`.
