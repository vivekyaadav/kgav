<!-- GENERATED FILE. Do not edit.
       make results-doc         regenerate
       make results-doc-check   fail if stale, write nothing

     Every figure is read from data/results or from the release.
     scripts/write_results_doc.py computes no AUC and applies no
     threshold of its own. -->

# kgav — what was measured

Release **v0.1** · `nodes_sha dc6d76a2b3cfb2a1` · `edges_sha 873623c9a9df2049`  
Generated 2026-11-21 from the results files on disk. `make check-results` says whether they describe this graph.

This document states what the project measured. `docs/evaluation.md` records how each result was arrived at, including the versions that were wrong; this one is the summary, with every qualification attached to the claim it qualifies.

## The four findings, in one table

| finding | direction | what it rests on |
| --- | --- | --- |
| **Chemical similarity predicts antiviral activity** | **positive** | 8 of 17 registered cells clear a corrected alpha above a degree baseline that does not |
| **It predicts FORWARD, out of time** | **positive** | 6 of 10 prospective cells clear, trained on evidence published at or before 2021 |
| The knowledge graph adds nothing over a degree baseline | negative | 0 of 61 metapath channel-cells both beat degree and are readable against it |
| An AUC above chance is not an expected hit rate | qualification | on the one assay-matched cell the top of the ranking holds no actives once near neighbours are withheld |

---

## 1. What was built

| node class | count |
| --- | ---: |
| `SmallMolecule` | 25,243 |
| `Protein` | 20,888 |
| `Gene` | 19,190 |
| `Pathway` | 13,236 |
| `TargetClass` | 66 |
| `OrganismTaxon` | 16 |
| `Disease` | 10 |
| **total** | **78,649** |

396,223 edge records over 8 declared reasoning metapaths. The full block-matrix view is `docs/figures/kgav_metagraph.png`, generated from the same release.

## 2. What counts as active

A compound is a **positive** for a virus when a measured potency is at or below **10 µM**, and a **negative** when a measurement exists and is above **10 µM**. A compound with no measurement is neither: it is a candidate, not a negative. Dose-series policy `include`.

**The labels are held out of every reasoning path.** `HAS_ANTIVIRAL_ACTIVITY_AGAINST`, `MEASURED_INACTIVE_AGAINST`, `TREATS` and `IN_TRIAL_FOR` are excluded when the graph is loaded for scoring, so no channel can read the answer sheet it is being scored against.

## 3. The reasoning channels

Eight metapaths, declared in `schema/kg_schema.yaml` and walked by degree-weighted path count (DWPC).

| channel | path |
| --- | --- |
| **M1** | `SmallMolecule → INHIBITS → Protein → ENCODED_BY → Gene → BELONGS_TO → OrganismTaxon` |
| **M2** | `SmallMolecule → TARGETS → Protein → HOST_FACTOR_FOR → OrganismTaxon` |
| **M3** | `SmallMolecule → TARGETS → Protein → PHYSICALLY_INTERACTS_WITH → Protein → ENCODED_BY → Gene → BELONGS_TO → OrganismTaxon` |
| **M4** | `SmallMolecule → TARGETS → Protein → PARTICIPATES_IN → Pathway → <PARTICIPATES_IN → Protein → HOST_FACTOR_FOR → OrganismTaxon` |
| **M5** | `SmallMolecule → TARGETS → Protein → PHYSICALLY_INTERACTS_WITH → Protein → HOST_FACTOR_FOR → OrganismTaxon` |
| **M6** | `SmallMolecule → CHEMICALLY_SIMILAR_TO → SmallMolecule → TARGETS → Protein → HOST_FACTOR_FOR → OrganismTaxon` |
| **M7** | `SmallMolecule → CHEMICALLY_SIMILAR_TO → SmallMolecule → INHIBITS → Protein → ENCODED_BY → Gene → BELONGS_TO → OrganismTaxon` |
| **M8** | `SmallMolecule → INHIBITS → Protein → FOLD_SIMILAR_TO → Protein → ENCODED_BY → Gene → BELONGS_TO → OrganismTaxon` |

## 4. The leakage controls

A structure scorer ranks a candidate by its similarity to known actives. Left uncontrolled it can score highly by finding a candidate's own near-twin — which is retrieval, not prediction. Three controls were registered, each removing a different route for that to happen.

| control | what it removes |
| --- | --- |
| **publication-disjoint** | every known active citing a publication the candidate also cites is dropped from that candidate's reference set |
| **symmetric k = 10** | the ten nearest known actives are dropped from EVERY candidate alike, active and inactive |
| **strength-matched** | k_cell nearest are dropped, where k_cell is that cell's own measured asymmetry under the publication control |

The publication control is **asymmetric between the classes** — a positive's nearest actives are often its own paper-mates, a negative's are not — which is why the other two were registered. Every result below is reported under all three.

## 5. Multiple comparisons

Every family is declared before its cells are read, and the corrected z follows from the family size rather than being chosen.

| family | size | alpha | where |
| --- | ---: | ---: | --- |
| ChEMBL structure cells | 14 | 0.003571 | `prereg_structure_scorer.md` |
| NCATS plate cells | 3 | 0.016667 | `prereg_ncats_entry_assays.md` |
| prospective cells | 10 | 0.005000 | `prereg_temporal_structure_scorer.md` |
| clinical labels | 1 | 0.050000 | `prereg_clinical_labels.md` |
| metapath channel-cells | 61 | 0.000820 | `evaluation.md` |

---

## Result 1 — chemical similarity separates actives from inactives

Registered by `docs/prereg_structure_scorer.md` (publication-disjoint), `docs/prereg_structure_scorer_symmetric.md` (symmetric k10), `docs/prereg_structure_scorer_strength_matched.md` (strength-matched) and `docs/prereg_ncats_entry_assays.md` (the plate cells).

**Positive.** Each cell is one virus in one assay population. `clears` means the corrected interval excludes chance ABOVE it. `degree` is the baseline: a compound's edge count, no chemistry in it. A structure column that does not beat its own degree column is not evidence about chemistry.

| cell | n+ | n− | pub-disj | sym k10 | strength | degree | verdict |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Yellow fever *(ChEMBL)* | 135 | 882 | 0.5545 | 0.8388\* | 0.8117\* | 0.5079 | **clears, degree does not** |
| JEV *(ChEMBL)* | 19 | 27 | 0.1452↓ | 0.8002\* | 0.7885\* | 0.5263 | **clears, degree does not** |
| SARS-CoV-2 *(ChEMBL)* | 1,761 | 5,320 | 0.7718\* | 0.8304\* | 0.7803\* | 0.5535\* | clears, so does degree |
| TBEV *(ChEMBL)* | 72 | 65 | 0.5764 | 0.7464\* | 0.7389\* | 0.4753 | **clears, degree does not** |
| West Nile *(ChEMBL)* | 261 | 311 | 0.6308\* | 0.6600\* | 0.7079\* | 0.5319 | **clears, degree does not** |
| Zika *(ChEMBL)* | 748 | 2,577 | 0.6316\* | 0.6812\* | 0.7032\* | 0.5342 | **clears, degree does not** |
| DENV-2 *(ChEMBL)* | 995 | 3,043 | 0.5742\* | 0.7888\* | 0.6921\* | 0.5638\* | clears, so does degree |
| HCoV-229E *(ChEMBL)* | 119 | 88 | 0.5917 | 0.7100\* | 0.6758\* | 0.6434\* | clears, so does degree |
| DENV-3 *(ChEMBL)* | 213 | 95 | 0.3861↓ | 0.6496\* | 0.6653\* | 0.4970 | **clears, degree does not** |
| SARS-CoV *(plate entry)* | 95 | 1,524 | refused | 0.6224\* | 0.6224\* | 0.6070\* | clears, so does degree |
| DENV-1 *(ChEMBL)* | 195 | 123 | 0.4489 | 0.6066\* | 0.5967\* | 0.5251 | **clears, degree does not** |
| HCoV-OC43 *(ChEMBL)* | 162 | 70 | 0.5257 | 0.5735 | 0.5903 | 0.6222\* | does not clear |
| SARS-CoV-2 *(plate CPE)* | 89 | 7,338 | refused | 0.6062\* | 0.5834\* | 0.5700 | **clears, degree does not** |
| SARS-CoV *(ChEMBL)* | 142 | 96 | 0.5098 | 0.5763 | 0.5763 | 0.6908\* | does not clear |
| MERS-CoV *(plate entry)* | 71 | 1,951 | refused | 0.5873\* | 0.5655 | 0.6326\* | does not clear |
| MERS-CoV *(ChEMBL)* | 131 | 41 | 0.6851\* | 0.4845 | 0.5045 | 0.6533\* | does not clear |
| DENV-4 *(ChEMBL)* | 153 | 109 | 0.4644 | 0.5868 | 0.4792 | 0.5245 | does not clear |

`*` = clears above chance · `↓` = clears BELOW chance · **8 of 17** cells clear above a degree baseline that does not, and that membership is **identical under both registered controls**.

---

## Result 2 — it predicts forward, out of time

**Positive, and the only out-of-time measurement here.** Everything above is cross-sectional: labels split by publication or by assay, never by time, so a scorer could be ranking chemistry that was already known when the test compound was published. This trains on activity published at or before **2021** and tests on activity published after it — the arrangement that answers *would this have worked*.

Registered by `docs/prereg_temporal_structure_scorer.md`. The screen this result justifies is registered separately, in `docs/prereg_denv2_prospective_screen.md`.

Split `v0.1-train2021` · undated `include` · overlap `drop-from-test` · family of 10, alpha 0.005000, z 2.807

| virus | references | test + | test − | AUC | corrected 95% CI | degree | verdict |
| --- | ---: | ---: | ---: | ---: | --- | ---: | --- |
| DENV-1 | 124 | 71 | 82 | 0.8025 | [0.7004, 0.9046] | — | **clears above chance** |
| DENV-4 | 86 | 67 | 86 | 0.7896 | [0.6830, 0.8961] | — | **clears above chance** |
| DENV-2 | 809 | 186 | 2,669 | 0.7518 | [0.6926, 0.8109] | — | **clears above chance** |
| DENV-3 | 143 | 70 | 79 | 0.7300 | [0.6132, 0.8468] | — | **clears above chance** |
| HCoV-229E | 43 | 76 | 32 | 0.7136 | [0.5713, 0.8560] | — | **clears above chance** |
| SARS-CoV-2 | 359 | 1,432 | 1,093 | 0.6209 | [0.5899, 0.6519] | — | **clears above chance** |
| HCoV-OC43 | 15 | 147 | 70 | 0.4348 | [0.3163, 0.5533] | — | does not clear |
| Zika | 181 | 567 | 2,495 | 0.4218 | [0.3859, 0.4577] | — | clears BELOW chance |
| Yellow fever | 105 | 30 | 51 | 0.3722 | [0.1962, 0.5482] | — | does not clear |
| West Nile | 80 | 181 | 189 | 0.3532 | [0.2735, 0.4329] | — | clears BELOW chance |

**6 of 10** clear above chance. Two clear BELOW it — they rank post-cutoff actives reliably last, which is a finding and not noise, and nothing here establishes that inverting them would hold forward.

**0 of 10 cells have a usable degree comparator.** The pre-cutoff graph reaches too few of the compounds being ranked to rank them at all, so the baseline column is empty by measurement rather than by omission — and these clearing cells are therefore **not shown to beat a baseline**, because there is none to beat.

### Three of those six are analogy, and three are not

The six clearing cells are **not one result**. The registered AUC, recomputed with the k nearest pre-cutoff actives withheld from every candidate alike — `k=0` is the registration's primary, `k=10` its secondary:

| virus | k=0 | k=5 | k=10 | k=20 | k=40 | clears up to |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| DENV-2 | 0.752\* | 0.719\* | 0.708\* | 0.687\* | 0.655\* | k≤40 |
| HCoV-229E | 0.714\* | 0.737\* | 0.711\* | 0.678\* | 0.662\* | k≤40 |
| SARS-CoV-2 | 0.621\* | 0.635\* | 0.636\* | 0.623\* | 0.614\* | k≤40 |
| DENV-4 | 0.790\* | 0.727\* | 0.687\* | 0.516 | 0.436 | k≤10 |
| DENV-1 | 0.802\* | 0.690\* | 0.618 | 0.462 | 0.357 | k≤5 |
| DENV-3 | 0.730\* | 0.654\* | 0.590 | 0.452 | 0.368 | k≤5 |

**3 survive to k = 40** — DENV-2, HCoV-229E, SARS-CoV-2 — with 40 of their closest known actives removed and the order still holding. That is a chemotype signal.

**2 stop clearing before k = 10** — DENV-1, DENV-3. Their signal was carried by a handful of near neighbours: retrieval of a known series, which is useful and is a narrower claim than the AUC sounds.

**The headline AUC is close to anti-correlated with robustness:**

```
by the published AUC (k=0) : DENV-1 > DENV-4 > DENV-2 > DENV-3 > HCoV-229E > SARS-CoV-2
by robustness      (k=40)  : HCoV-229E > DENV-2 > SARS-CoV-2 > DENV-4 > DENV-3 > DENV-1
```

DENV-1 has the highest primary AUC of the six and loses the most (0.802 to 0.357). SARS-CoV-2 has the lowest and barely moves (0.621 to 0.614). **A cell's headline number does not say which kind of result it is.**

Registered as H3 in `docs/prereg_temporal_structure_scorer.md`, left open there as "a pointer, not a finding", and closed by the amendment of 2026-11-21. `make prospective-analogy`.

---

## Result 3 — the knowledge graph adds nothing over a degree baseline

**Negative, and it is the finding this project was built to test.** The graph is the thing under test; `degree` — a compound's edge count, with no biology in it — is what any reasoning channel has to beat.

| metapath channel-cells | count |
| --- | ---: |
| measured | 61 |
| beat the degree baseline on lift | 16 |
| readable against degree at all | 4 |
| **both beat degree and are readable against it** | **0** |

On **SARS-CoV-2**, the assay-matched plate cell under the publication-disjoint control:

| scorer | AUC | 95% CI |
| --- | ---: | --- |
| strongest single channel (M4) | 0.5400 | [0.478, 0.602] |
| all eight channels COMBINED | 0.5666 | [0.504, 0.629] |
| **degree baseline** | 0.5700 | [0.508, 0.632] |

The ensemble of all eight reasoning channels scores **-0.0034** against the baseline that uses no biology at all.

## Result 4 — the pathway overlap is null

Registered by `docs/prereg_pathway_overlap.md`.

**Negative.** Whether a drug's targets sit in the same pathways as a virus's host-dependency factors, against a degree-matched permutation null.

| channel | rows measured | clearing the corrected alpha |
| --- | ---: | ---: |
| dependency | 8 | 0 |
| interactome | 14 | 0 |
| **total** | **22** | **0** |

Permutations: 100,000, seed `0`, degree basis `study`.

## Result 5 — the one independent label axis does not corroborate it

Registered by `docs/prereg_clinical_labels.md`.

**Negative, and it points the other way.** Every AUC above rests on `HAS_ANTIVIRAL_ACTIVITY_AGAINST` — one predicate, in-vitro concentration-response at one threshold. The clinical axis is independent of it: compounds in a trial or approved for COVID-19 (`MONDO:0100096`).

| scorer | AUC | 95% CI | verdict |
| --- | ---: | --- | --- |
| structure (symmetric k10) | 0.4087 | [0.3734, 0.4440] | **below chance** |
| degree baseline | 0.7751 | [0.7385, 0.8116] | well above chance |

221 positives against 12,595 negatives, of which 200 of the positives were measured INACTIVE in vitro. Reaching a trial is predicted by how much a compound has been studied, not by its resemblance to the compounds active in vitro — so this axis measures something else rather than confirming the one above, and cannot be used to validate a structure scorer.

---

## The qualification that applies to every AUC above

**An AUC is not a hit rate.** It is the probability that a random active outranks a random inactive — a statement about the WHOLE ordering. A screen only ever orders the top of a list, and the two can disagree completely.

On the one population matched by assay rather than by source — SARS-CoV-2, pool 7,427, base rate 1.20%:

| control | precision@10 | enrichment over base |
| --- | ---: | ---: |
| no exclusion | 50.0% | 41.72× |
| symmetric k = 10 | 0.0% | 0.00× |
| strength-matched (k_cell 38) | 0.0% | 0.00× |

That cell's strength-matched AUC is **0.5834 and it CLEARS** its registered alpha. The AUC cannot show this, because it averages over a ranking nobody screens. A precision that vanishes when a candidate's nearest known actives are withheld was analogue retrieval, not prediction.

Prospectively, at k = 10 under the registered secondary control:

| virus | base rate | precision@10 | enrichment |
| --- | ---: | ---: | ---: |
| DENV-2 | 6.5% | 70% | **10.74×** |
| DENV-4 | 43.8% | 80% | 1.83× |
| SARS-CoV-2 | 56.7% | 90% | 1.59× |
| HCoV-229E | 70.4% | 100% | 1.42× |
| DENV-1 | 46.4% | 60% | 1.29× |
| DENV-3 | 47.0% | 60% | 1.28× |
| HCoV-OC43 | 67.7% | 55% | 0.81× |
| Zika | 18.5% | 10% | 0.54× |
| West Nile | 48.9% | 20% | 0.41× |
| Yellow fever | 37.0% | 10% | 0.27× |

**Enrichment, not precision, is the column that matters.** A precision of 90% in a pool that is already 70% active is a ranking that added almost nothing. Only DENV-2 shows a lift that would survive a feasible screen; `make screen-power` sizes one for every cell, and six of the ten cannot be powered at any arm size up to 120.

`make hit-rate` reproduces this. It is **post-hoc**: precision@k was computed after the registered AUCs were read, so it qualifies a claim rather than being one.

---

## What is claimed, and what is not

**Claimed.** On this release, ranking compounds by maximum ECFP4 Tanimoto to a virus's known actives separates measured actives from measured inactives in 8 of 17 registered cells, above a degree baseline that does not clear, with the same membership under both registered controls — and does so **prospectively** in 6 of 10 cells trained on pre-2021 evidence.

**Not claimed.**

- *That the knowledge graph contributes to this.* It does not. The scorer uses no graph edge, and `make auc-contribution` measures the graph's contribution to every scorer here as zero or negative against its own degree baseline.
- *That a high AUC implies a usable hit rate.* On the one assay-matched cell it does not, measured above.
- *That a clearing cell beats a baseline prospectively.* None of them is shown to: no prospective cell has a usable degree comparator.
- *Anything clinical.* The endpoint is an in-vitro EC50 in one assay, often in Vero E6, which lacks meaningful TMPRSS2 — and the one independent label axis anti-correlates with this statistic.
- *That a suggestion is a prediction of activity.* It is the claim that a compound RESEMBLES a known active, with the known active named.

## Reproducing all of it

```
make rebuild                 # the graph, ~27 min, byte-reproducible
make check-results           # which results describe this graph
make scorer-table            # Result 1, every control side by side
make temporal-structure-scorer   # Result 2
make auc-contribution        # Result 3, per-approach shares
make pathway-landscape       # Result 4
make clinical-labels         # Result 5
make hit-rate                # the qualification
make screen-power            # what a wet-lab screen would need
```

## Every registration

Each was signed before the result it governs existed. `reported here` is derived from whether this document names it above.

| registration | subject | reported here |
| --- | --- | --- |
| `prereg_clinical_labels.md` | a label that is not an in-vitro assay | yes |
| `prereg_denv2_prospective_screen.md` | a prospective wet-lab screen, DENV-2 | yes |
| `prereg_m1_potency.md` | letting potency into the direct-acting channel | no — see `evaluation.md` |
| `prereg_ncats_entry_assays.md` | the NCATS pseudotyped-entry assays as plate-matched labels | yes |
| `prereg_partial_degree.md` | does any channel carry signal beyond degree | no — see `evaluation.md` |
| `prereg_pathway_overlap.md` | the pathway and interactome overlap, re-run | yes |
| `prereg_phase4_hcv.md` | the HCV retrospective, as the evidence permits it | no — see `evaluation.md` |
| `prereg_structure_scorer.md` | a structure-derived scorer, and the leak it has to survive | yes |
| `prereg_structure_scorer_strength_matched.md` | the same strength, symmetrically | yes |
| `prereg_structure_scorer_symmetric.md` | a control that handicaps both classes equally | yes |
| `prereg_temporal_auc.md` | a tie-safe metric for the temporal protocol | no — see `evaluation.md` |
| `prereg_temporal_structure_scorer.md` | does the structure scorer predict forward, or only separate what is already known? | yes |

