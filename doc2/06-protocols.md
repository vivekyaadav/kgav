# 06 — The evaluation protocols

Every protocol here was registered before its p-values existed, writes a
stamped artifact, and fixes its correction family in advance. Where a protocol
refused to run, **the refusal is the artifact** —
1 of 36 files in `data/results/` is a
recorded refusal:

- `temporal_2021_include_overlap-keep.json`: L5 label overlap: 12 positive compounds for NCBITaxon:2697049 appear in BOTH train and test labels; L5 label overlap: 38 negative compounds for NCBITaxon:2697049 appear in BOTH train and test labels; L5 label overlap: 1 negative compounds for NCBITaxon:694009 appear in BOTH train and test labels

## 1. The three non-learned baselines

`src/kgav/baselines.py`, `scripts/run_baselines.py`. *A learned scorer that
does not clearly beat all three is not a result.* Most published
knowledge-graph repurposing work never runs this comparison, which is why so
much of it is unreliable.

1. **degree ranking** — the null hypothesis.
2. **DWPC over the metapaths** — degree-weighted path count, Hetionet-style.
   Each path contributes the product of its nodes' degrees raised to −0.4, so
   a route through a promiscuous protein counts for less than one through a
   specific target. Without it, path counts simply measure how well studied
   the intermediates are.
3. **network proximity** — drug targets to host factors through the
   interactome.

**Held-out predicates are never traversed.** `HAS_ANTIVIRAL_ACTIVITY_AGAINST`, `IN_TRIAL_FOR`, `MEASURED_INACTIVE_AGAINST`, `TREATS` are the
only edges connecting most compounds directly to a virus, so excluding them is
what makes this a prediction rather than a lookup.

`make degree-sensitivity` answers the obvious follow-up — does the DWPC step
weight itself read the answer sheet? Three arms over 61 cells:
as published, with held-out edges out of the step weight only, and out of the
hub cut as well. **4 cells move and 0 cross alpha.**
The five channels that never move are the control, and they are what licenses
the claim that the finding is small. Running only the channels that *can* move
was considered and rejected: it would save five minutes and the ability to say
so.

## 2. Channel calibration

`src/kgav/calibration.py`. `baselines.combine` ranks a drug by its **best**
percentile across metapaths, which was right when the alternative was summing
incomparable scales and is wrong once some channels carry no signal: max over
*k* channels is biased upward by *k*, so a drug at the top of a pure-noise pool
scores 1.0 and outranks one with a real route. The module's own docstring
records the case that prompted it — a protocol on which the ensemble scored
below its own best single channel, so **pooling made the ranking worse**. On
the artifact as it stands the ensemble reads 0.5666 against
M4's 0.5400, and the reason the ensemble is
higher is its coverage: its floor alone is 0.5578.

What a channel must earn is not a point AUC above 0.5 — half of pure
noise clears that — but a 95% interval that **excludes** it.
0 of the 8 channels does
(none), and `degree` and the ensemble are excluded from
calibration **by role** — the first is the null hypothesis a ranking has to
beat, the second is the thing being built — so a calibrated combine over this
release returns **nothing**. That is the correct output: an agent whose
channels are all at chance should say it cannot rank, not emit a confident
order derived from noise.

Weights are **loaded, not hardcoded**. A weight written into the source is a
fitted parameter with no record of what fitted it. These come from a results
file carrying its own protocol, label source, n, and whether the
publication-disjoint control was applied — so a weight can always be traced to
the run that earned it, and a weight measured under a confounded protocol can
be **refused** by name. The module records four different AUCs produced for
one channel in one evening by four protocols, as the controls were added; a
number with no protocol attached is not a measurement.

## 3. The hard-negative evaluation

`scripts/hard_negatives.py`. The registered headline protocol:

| | |
|:---|:---|
| label source | `infores:ncats-opendata` |
| assay type | `cell_based_antiviral` |
| inactive floor | 10,000 nM |
| selectivity | all / positives |
| publication-disjoint | True |
| temporal section | skipped — the section is empty in the artifact |

Restricting to one source **and one protocol** is what makes a cell
plate-matched. One source can publish several assays against one virus, and
grouped by virus alone those merge into a single cell — which destroys the only
thing such a cell is worth, because its value is that **one** protocol over
**one** concentration range produced both classes. The two filters are
independent on purpose: `sources` answers *whose measurement*, `assay_types`
answers *measuring what*, and a cell needs both pinned before "matched by
construction" is true of it.

| virus | channel | + | − | AUC | 95% interval | positives reached | floor | skill |
|:---|:---|---:|---:|:---|:---|---:|---:|---:|
| SARS-CoV-2 | M1 | 89 | 7,338 | 0.5129 | [0.4520, 0.5737] | 3.4% | 0.5128 | 0.220 |
| SARS-CoV-2 | M2 | 89 | 7,338 | 0.5063 | [0.4458, 0.5669] | 12.4% | 0.5045 | 0.260 |
| SARS-CoV-2 | M3 | 89 | 7,338 | 0.5339 | [0.4724, 0.5954] | 16.9% | 0.5328 | 0.129 |
| SARS-CoV-2 | M4 | 89 | 7,338 | 0.5400 | [0.4783, 0.6016] | 33.7% | 0.5292 | 0.229 |
| SARS-CoV-2 | M5 | 89 | 7,338 | 0.5391 | [0.4774, 0.6007] | 28.1% | 0.5350 | 0.136 |
| SARS-CoV-2 | M6 | 89 | 7,338 | 0.5155 | [0.4546, 0.5765] | 5.6% | 0.5155 | 0.080 |
| SARS-CoV-2 | M7 | 89 | 7,338 | 0.4984 | [0.4381, 0.5586] | 0.0% | 0.4984 | — |
| SARS-CoV-2 | M8 | 89 | 7,338 | 0.5056 | [0.4451, 0.5662] | 1.1% | 0.5056 | — |
| SARS-CoV-2 | COMBINED | 89 | 7,338 | 0.5666 | [0.5044, 0.6289] | 42.7% | 0.5578 | 0.133 |
| SARS-CoV-2 | degree | 89 | 7,338 | 0.5700 | [0.5077, 0.6323] | 100.0% | 0.5000 | 0.140 |

## 4. The structure-derived scorer, under three controls

`scripts/structure_scorer.py`. The statistic is **max Tanimoto to the known
actives**, ECFP4 at radius 2, folded to 2,048 bits. The three
runs are three *controls* on the same statistic, each with its own
registration, each writing its own file — because the publication-disjoint run
is the comparator for the others and overwriting it would destroy the only
thing they are being compared against.

| run | control | family | corrected z | registration |
|:---|:---|---:|---:|:---|
| `structure_scorer.json` | `publication` — withhold references sharing a publication with the candidate | 14 | 2.9137 | `docs/prereg_structure_scorer.md` |
| `structure_scorer_symmetric.json` | `knn` — withhold each candidate's 10 nearest references, **symmetrically** | 14 | 2.9137 | `docs/prereg_structure_scorer_symmetric.md` |
| `structure_scorer_strength_matched.json` | `strength` — withhold each cell's own measured asymmetry, `k_cell` | 14 | 2.9137 | `docs/prereg_structure_scorer_strength_matched.md` |

The strength-matched run **derives every `k` from the publication-disjoint
run** and refuses if that file is absent, because there is then no strength to
match:

| virus | k_cell |
|:---|---:|
| SARS-CoV-2 | 38 |
| DENV-2 | 27 |
| HCoV-OC43 | 27 |
| DENV-3 | 24 |
| DENV-1 | 14 |
| MERS-CoV | 14 |
| Yellow fever | 13 |
| DENV-4 | 12 |
| JEV | 12 |
| TBEV | 11 |
| HCoV-229E | 11 |
| SARS-CoV | 10 |
| Zika | 7 |
| West Nile | 4 |

The plate-matched secondary populations carry their own smaller family
(3 tests, z 2.394), because they are a different question
asked of fewer cells.

## 5. The prospective (temporal) protocol

`src/kgav/temporal.py` cuts the split; `scripts/temporal_structure_scorer.py`
scores it. Registered in advance by
`docs/prereg_temporal_structure_scorer.md`, signed before any AUC against a
post-cutoff label existed.

| | |
|:---|:---|
| split | `v0.1-train2021` |
| cutoff | 2021 |
| undated policy | `include` |
| overlap policy | `drop-from-test` |
| family | 10 cells, alpha 0.005, z 2.807 |
| evaluable | 10 of 15 candidate cells |

**The cutoff is chosen from the label distribution, not from convention.**
5,181 activity edges are dated 2007–2021 and 2,027 are dated 2022–2025, so
both sides are large enough to measure. The 2021 spike of 4,541 labels is the
COVID screening wave, which makes the test side the steadier post-pandemic
tail — arguably the right shape, since *what gets screened next* is the real
question.

**The trap this module exists to avoid.** Filtering only the **labels** leaks
the future: if post-cutoff `INHIBITS` and `TARGETS` edges stay in the graph, a
compound whose nsp5 activity was published in 2023 gets an M1 path built from
2023 evidence to predict a 2023 label. The training graph is filtered on the
same date as the labels, and the audits check that it was.

**Undated edges are a declared policy, not a default.** Computed edges
(`CHEMICALLY_SIMILAR_TO`) carry the `1970` placeholder deliberately, so they
appear in every training graph. Three policies exist — `include`, `exclude`,
`computed_only` — and all three are run and written, so the sensitivity is on
disk rather than asserted.

**The overlap policy is why one artifact is a refusal.** Under
`overlap_policy: keep`, 12 positive and 38 negative compounds for SARS-CoV-2
appear on **both** sides of the cut. The protocol refuses rather than
reporting, and `temporal_2021_include_overlap-keep.json` carries the reasons
and **no results block** — because a refusal that still carried numbers would
be read as numbers.

**The degree comparator is required to be usable.** A comparator reaching
under 20% of the compounds being ranked cannot
be beaten informatively, and on this release
**0** cells have a usable one — the training
graph reaches 0.3%–6.5% of the post-cutoff compounds. So the registered
decision rule's second clause ("and beats a usable comparator") is satisfied
by **nothing**, and the result is reported as *measurable, not comparable*.

## 6. The clinical-label axis

`scripts/clinical_labels.py`, registered by
`docs/prereg_clinical_labels.md`. One cell, family of
1, alpha 0.05. It is the **one independent
label axis**: did a compound reach a trial or an approval for the disease this
virus causes, rather than did it score under a threshold in a cell assay.

| population | compounds |
|:---|---:|
| labelled | 12,816 |
| invitro active | 1,768 |
| invitro inactive | 11,048 |
| positives | 221 |
| negatives | 12,595 |
| positives invitro active | 21 |
| positives invitro inactive | 200 |

## 7. The pathway-overlap landscape

`src/kgav/landscape.py`, `scripts/pathway_landscape.py`, registered by
`docs/prereg_pathway_overlap.md`. The overlap between the drugged human
proteome (1,944 proteins reached by at least one `TARGETS`
edge) and each virus's CRISPR dependency factors, measured twice: naive, and
against a background sampled to **match the observed degree distribution** —
because well-studied proteins are both likelier to be drugged and likelier to
be called dependency factors.

100,000 permutations, 44 tests, degree basis
`study`, seeds `per-row: Random(f'{seed}|{channel}|{virus}')`.

**The test is one-sided and says so.** `empirical_p` measures the
**enrichment** tail only. A ratio below 1 is therefore a description, not a
result: nothing here has tested depletion, and printing "ratios 0.00–0.94"
beside "no enrichment" invites the reader to hear a finding that was never
measured. The depletion column is reported separately and labelled.

## 8. Correction families, fixed in advance

| family | tests | corrected alpha | survivors |
|:---|---:|:---|:---|
| readable_cells | 13 | 0.0038462 | SARS-CoV-2 M7 |
| cells_computed | 61 | 0.0008197 | none |
| cells_planned | 120 | 0.0004167 | none |
| planned_plus_baselines | 150 | 0.0003333 | none |

> the family is reported, not chosen: this project fixes a correction family in a pre-registration before the p-values are in hand, and these were not

## 9. The post-hoc measurements, labelled as such

Three measurements were computed **after** the registered AUCs were read, on
the same held-out labels. They qualify a claim rather than being one, they are
labelled post-hoc, and they are deliberately **not** written to
`data/results/` where they would read as claim-carrying results:

- `make hit-rate` — precision@k, at three exclusion strengths. AUC is a
  statement about a whole ordering and a screen only ever sees the top of a
  list.
- `make protocol-disagreement` — why one cell can anti-rank cross-sectionally
  and clear prospectively.
- `make prospective-analogy` — the registered prospective AUC recomputed with
  the *k* nearest pre-cutoff actives withheld, k = 0, 5, 10, 20, 40. Two rungs
  reproduce published columns, so the ladder is checkable against the run it
  extends.
