# 07 — Every measured result

Each table is read out of the artifact named above it. A verdict is the run's
own `clears_corrected_alpha` **and** its `direction`, never a comparison a
reader is asked to make: an interval that excludes chance **below** it is a
measured anti-ranking, which is a different finding from a null.

---

## The structure scorer, cross-sectionally

The same statistic — max Tanimoto to the known actives, ECFP4 at
radius 2 folded to 2,048 bits — under three controls of increasing severity. Read
them left to right: a cell that survives the third has had its own chemical
series taken away.

### Control 1 — publication-disjoint (`structure_scorer.json`)

19,168 fingerprinted compounds, family 14, corrected
z 2.9137. **5 of 14 clear above chance**
(DENV-2, MERS-CoV, SARS-CoV-2, West Nile, Zika); **2 clear below** (DENV-3, JEV).

| virus | active | inactive | AUC | corrected 99.643% interval | verdict |
|:---|---:|---:|:---|:---|:---:|
| SARS-CoV-2 | 1,761 | 5,320 | 0.7718 | [0.7513, 0.7924] | **clears** |
| MERS-CoV | 131 | 41 | 0.6851 | [0.5566, 0.8135] | **clears** |
| Zika | 748 | 2,577 | 0.6316 | [0.5965, 0.6666] | **clears** |
| West Nile | 261 | 311 | 0.6308 | [0.5626, 0.6990] | **clears** |
| HCoV-229E | 119 | 88 | 0.5917 | [0.4767, 0.7067] | no verdict |
| TBEV | 72 | 65 | 0.5764 | [0.4344, 0.7184] | no verdict |
| DENV-2 | 995 | 3,043 | 0.5742 | [0.5431, 0.6053] | **clears** |
| Yellow fever | 135 | 882 | 0.5545 | [0.4752, 0.6337] | no verdict |
| HCoV-OC43 | 162 | 70 | 0.5257 | [0.4062, 0.6453] | no verdict |
| SARS-CoV | 142 | 96 | 0.5098 | [0.3986, 0.6210] | no verdict |
| DENV-4 | 153 | 109 | 0.4644 | [0.3587, 0.5701] | no verdict |
| DENV-1 | 195 | 123 | 0.4489 | [0.3518, 0.5461] | no verdict |
| DENV-3 | 213 | 95 | 0.3861 | [0.2825, 0.4896] | **clears BELOW chance** |
| JEV | 19 | 27 | 0.1452 | [0.0000, 0.3061] | **clears BELOW chance** |

### Control 2 — symmetric k=10 (`structure_scorer_symmetric.json`)

**10 of 14 clear above chance**
(DENV-1, DENV-2, DENV-3, HCoV-229E, JEV, SARS-CoV-2, TBEV, West Nile, Yellow fever, Zika); 0 below (none).

| virus | active | inactive | AUC | corrected 99.643% interval | verdict |
|:---|---:|---:|:---|:---|:---:|
| Yellow fever | 135 | 882 | 0.8388 | [0.7751, 0.9026] | **clears** |
| SARS-CoV-2 | 1,761 | 5,320 | 0.8304 | [0.8119, 0.8488] | **clears** |
| JEV | 19 | 27 | 0.8002 | [0.5973, 1.0000] | **clears** |
| DENV-2 | 995 | 3,043 | 0.7888 | [0.7621, 0.8154] | **clears** |
| TBEV | 72 | 65 | 0.7464 | [0.6253, 0.8675] | **clears** |
| HCoV-229E | 119 | 88 | 0.7100 | [0.6069, 0.8131] | **clears** |
| Zika | 748 | 2,577 | 0.6812 | [0.6470, 0.7154] | **clears** |
| West Nile | 261 | 311 | 0.6600 | [0.5932, 0.7268] | **clears** |
| DENV-3 | 213 | 95 | 0.6496 | [0.5554, 0.7439] | **clears** |
| DENV-1 | 195 | 123 | 0.6066 | [0.5139, 0.6994] | **clears** |
| DENV-4 | 153 | 109 | 0.5868 | [0.4842, 0.6895] | no verdict |
| SARS-CoV | 142 | 96 | 0.5763 | [0.4676, 0.6850] | no verdict |
| HCoV-OC43 | 162 | 70 | 0.5735 | [0.4568, 0.6901] | no verdict |
| MERS-CoV | 131 | 41 | 0.4845 | [0.3328, 0.6361] | no verdict |

### Control 3 — strength-matched (`structure_scorer_strength_matched.json`)

The sternest run and the one every downstream decision quotes: each
candidate's `k_cell` nearest known actives withheld, with `k_cell` derived
from control 1. **10 of 14 clear above chance**
(DENV-1, DENV-2, DENV-3, HCoV-229E, JEV, SARS-CoV-2, TBEV, West Nile, Yellow fever, Zika); **0 below** (none).

| virus | active | inactive | AUC | corrected 99.643% interval | verdict |
|:---|---:|---:|:---|:---|:---:|
| Yellow fever | 135 | 882 | 0.8117 | [0.7443, 0.8791] | **clears** |
| JEV | 19 | 27 | 0.7885 | [0.5810, 0.9960] | **clears** |
| SARS-CoV-2 | 1,761 | 5,320 | 0.7803 | [0.7600, 0.8005] | **clears** |
| TBEV | 72 | 65 | 0.7389 | [0.6164, 0.8614] | **clears** |
| West Nile | 261 | 311 | 0.7079 | [0.6440, 0.7717] | **clears** |
| Zika | 748 | 2,577 | 0.7032 | [0.6695, 0.7369] | **clears** |
| DENV-2 | 995 | 3,043 | 0.6921 | [0.6624, 0.7219] | **clears** |
| HCoV-229E | 119 | 88 | 0.6758 | [0.5684, 0.7832] | **clears** |
| DENV-3 | 213 | 95 | 0.6653 | [0.5728, 0.7579] | **clears** |
| DENV-1 | 195 | 123 | 0.5967 | [0.5034, 0.6901] | **clears** |
| HCoV-OC43 | 162 | 70 | 0.5903 | [0.4750, 0.7056] | no verdict |
| SARS-CoV | 142 | 96 | 0.5763 | [0.4676, 0.6850] | no verdict |
| MERS-CoV | 131 | 41 | 0.5045 | [0.3538, 0.6552] | no verdict |
| DENV-4 | 153 | 109 | 0.4792 | [0.3735, 0.5850] | no verdict |

#### Two definitions of "below", and the read-out gates on the stricter one

The verdict column above is the **run's**: the corrected interval excludes
chance on the wrong side. `scripts/repurposing_candidates.py` gates a
shortlist on a stricter tier — `auc < 0.5 OR direction == below` — because a
ranking of unlabelled compounds out of a cell whose point estimate sits under
chance is not something to print with a caveat attached. The two disagree on
exactly one cell, and conflating them is how a document comes to contradict
the tool it describes.

| virus | AUC | interval verdict | read-out tier |
|:---|:---|:---|:---|
| Yellow fever | 0.8117 | clears | `clears` |
| JEV | 0.7885 | clears | `clears` |
| SARS-CoV-2 | 0.7803 | clears | `clears` |
| TBEV | 0.7389 | clears | `clears` |
| West Nile | 0.7079 | clears | `clears` |
| Zika | 0.7032 | clears | `clears` |
| DENV-2 | 0.6921 | clears | `clears` |
| HCoV-229E | 0.6758 | clears | `clears` |
| DENV-3 | 0.6653 | clears | `clears` |
| DENV-1 | 0.5967 | clears | `clears` |
| HCoV-OC43 | 0.5903 | no verdict | `weak` |
| SARS-CoV | 0.5763 | no verdict | `weak` |
| MERS-CoV | 0.5045 | no verdict | `weak` |
| DENV-4 | 0.4792 | no verdict | `below` |

By that tier: **10 clears**, 3 weak,
**1 below** (DENV-4). The read-out prints no
list for the last group.

### The plate-matched populations

Matched by **construction** rather than by a filter: one source, one protocol,
one concentration range produced both classes. Family 3,
z 2.394.

| control | population | virus | active | inactive | AUC | verdict |
|:---|:---|:---|---:|---:|:---|:---:|
| symmetric k=10 | NCATS CPE | SARS-CoV-2 | 89 | 7,338 | 0.6062 | **clears** |
| symmetric k=10 | NCATS entry | MERS-CoV | 71 | 1,951 | 0.5873 | **clears** |
| symmetric k=10 | NCATS entry | SARS-CoV | 95 | 1,524 | 0.6224 | **clears** |
| strength-matched | NCATS CPE | SARS-CoV-2 | 89 | 7,338 | 0.5834 | **clears** |
| strength-matched | NCATS entry | MERS-CoV | 71 | 1,951 | 0.5655 | no verdict |
| strength-matched | NCATS entry | SARS-CoV | 95 | 1,524 | 0.6224 | **clears** |

---

## The structure scorer, prospectively

`temporal_structure_scorer.json`. Trained on activity published at or before
2021, scored on compounds measured after. This is the strongest
positive result in the project and the one with the most ways to be
misread.

**6 of 10 cells clear above chance**:
DENV-1, DENV-2, DENV-3, DENV-4, HCoV-229E, SARS-CoV-2. **2 clear below chance**:
West Nile, Zika — the scorer, trained on what was known and asked to rank
what came next, put the actives **last**.

| virus | references | test + | test − | AUC | corrected interval | k=10 | degree | verdict |
|:---|---:|---:|---:|:---|:---|:---|:---|:---:|
| DENV-1 | 124 | 71 | 82 | 0.8025 | [0.7004, 0.9046] | 0.6181 | — | **clears** |
| DENV-2 | 809 | 186 | 2,669 | 0.7518 | [0.6926, 0.8109] | 0.7076 | 0.5127 | **clears** |
| DENV-3 | 143 | 70 | 79 | 0.7300 | [0.6132, 0.8468] | 0.5901 | 0.5020 | **clears** |
| DENV-4 | 86 | 67 | 86 | 0.7896 | [0.6830, 0.8961] | 0.6868 | — | **clears** |
| HCoV-229E | 43 | 76 | 32 | 0.7136 | [0.5713, 0.8560] | 0.7105 | — | **clears** |
| HCoV-OC43 | 15 | 147 | 70 | 0.4348 | [0.3163, 0.5533] | 0.4838 | 0.5031 | no verdict |
| SARS-CoV-2 | 359 | 1,432 | 1,093 | 0.6209 | [0.5899, 0.6519] | 0.6360 | 0.5071 | **clears** |
| West Nile | 80 | 181 | 189 | 0.3532 | [0.2735, 0.4329] | 0.4521 | 0.4923 | **clears BELOW chance** |
| Yellow fever | 105 | 30 | 51 | 0.3722 | [0.1962, 0.5482] | 0.3863 | 0.5317 | no verdict |
| Zika | 181 | 567 | 2,495 | 0.4218 | [0.3859, 0.4577] | 0.4572 | 0.5044 | **clears BELOW chance** |

Three things this table is not:

- **Not comparable to a baseline.** 0 cells have
  a usable degree comparator, because the pre-cutoff training graph reaches
  0.3%–6.5% of the compounds being ranked. The registered rule's second clause
  is satisfied by nothing, so the result is *measurable, not comparable*.
- **Not uniform in kind.** `make prospective-analogy` recomputes each clearing
  cell's AUC with the *k* nearest pre-cutoff actives withheld (k = 0, 5, 10,
  20, 40; k=0 and k=10 reproduce the published columns). Three of the six
  survive to k=40 — DENV-2, HCoV-229E, SARS-CoV-2 — two stop clearing before
  k=10, and one clears at k=10 but not k=20. **The headline AUC is close to
  anti-correlated with robustness**: the cell with the highest primary loses
  the most.
- **Not a contradiction where it disagrees with the cross-sectional run.**
  DENV-4 anti-ranks cross-sectionally (0.4792, with
  `k_cell` = 12) and clears prospectively
  (0.7896). `make protocol-disagreement` shows the two
  protocols are nested
  exactly in their compounds and differ only in their **reference set**: for a
  post-cutoff active, the cross-sectional control removes its own chemical
  series while the prospective one cannot, because that series is not in the
  pre-cutoff reference set at all. Both results stand; they answer different
  questions.

---

## The graph's own reasoning channels

`hard_negatives.json`, under the registered publication-disjoint protocol on
the plate-matched population.

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

Read the `floor` and `skill` columns before the AUC.
0 of the 8 metapath channels has an
interval excluding 0.5 (none). The two scorers that
do are the ensemble (0.5666 [0.5044, 0.6289], skill
13.3%) and the **degree baseline** (0.5700
[0.5077, 0.6323], skill 14.0%) — and the difference between them
is the floor, not the AUC: degree reaches 100.0% of the positives so
its floor is exactly 0.5 and all of its margin is ranking, while the
ensemble reaches 42.7% and its coverage asymmetry **alone** scores
0.5578.

`metapath_comparison.json` covers 15 viruses and the channels
M1, M2, M3, M4, M5, M6, M7, M8, COMBINED, degree. With the correction family fixed in advance, the survivors
are:

| family | tests | corrected alpha | survivors |
|:---|---:|:---|:---|
| readable_cells | 13 | 0.0038462 | SARS-CoV-2 M7 |
| cells_computed | 61 | 0.0008197 | none |
| cells_planned | 120 | 0.0004167 | none |
| planned_plus_baselines | 150 | 0.0003333 | none |

---

## The clinical-label axis

`clinical_labels.json`. SARS-CoV-2 / MONDO:0100096,
271 `IN_TRIAL_FOR` and 1 `TREATS` edges,
family 1, alpha 0.05.

| scorer | AUC | corrected interval | verdict |
|:---|---:|:---|:---|
| max Tanimoto, symmetric k=10 | 0.4087 | [0.3734, 0.4440] | clears **below** chance |
| degree | 0.7751 | [0.7385, 0.8116] | clears above chance |

Degree at or above structure: **True**. What predicts
whether a compound reached a trial is how well studied it is. This is the one
independent label axis in the repository and it does not support the
structure scorer.

| population | compounds |
|:---|---:|
| labelled | 12,816 |
| invitro active | 1,768 |
| invitro inactive | 11,048 |
| positives | 221 |
| negatives | 12,595 |
| positives invitro active | 21 |
| positives invitro inactive | 200 |

---

## The pathway-overlap landscape

`pathway_landscape.json`. 100,000 degree-matched permutations,
44 tests, 1,944 proteins in the drugged
proteome.

| virus | dependency factors | observed overlap | expected | naive p | degree-matched p | p (depletion) |
|:---|---:|---:|---:|---:|---:|---:|
| DENV-2 | 50 | 2 | 6.28 | 0.9581 | 0.9907 | 0.0394 |
| HCoV-229E | 1,367 | 128 | 135.13 | 0.5963 | 0.7643 | 0.2665 |
| HCoV-NL63 | 57 | 10 | 7.39 | 0.0412 | 0.1963 | 0.8881 |
| HCoV-OC43 | 657 | 63 | 65.16 | 0.4964 | 0.6336 | 0.4176 |
| MERS-CoV | 678 | 59 | 66.41 | 0.7883 | 0.8530 | 0.1809 |
| SARS-CoV-2 | 3,460 | 383 | 358.31 | 0.0005 | 0.0665 | 0.9412 |
| West Nile | 18 | 0 | 2.51 | 1.0000 | 1.0000 | 0.0663 |
| Zika | 815 | 76 | 80.01 | 0.5956 | 0.7046 | 0.3396 |

**8 viruses have no cell at all**: DENV-1, DENV-3, DENV-4, HCoV-HKU1, JEV, SARS-CoV, TBEV, Yellow fever. A null
result is only a finding if the reader knows what would have survived, so
`docs/scope_and_limits.md` carries the minimum detectable effect per cell —
failing to detect a 2× excess is not evidence against a 2× excess when the
test needed 7×.

`landscape_count_move.json` records why two cells moved between releases, and
that the mechanism accounts for **every** move: a protein already inside the
factor set joined the drugged proteome, with no edge in its own channel
changing.

---

## The circularity confound, measured

`label_overlap.json`. Of 943 compounds carrying both an `INHIBITS`
edge and an antiviral label, **779 cite the same publication on
both** — 88.5% of those where both sides have one — because one
paper reports the target IC50 and the cell-based EC50 together.
58 pairs cannot be withheld at all, because one side
carries no PMID. 22 viral proteins are reached.

This is why `--publication-disjoint` exists, why the NCATS panel was ingested,
and why an M1 figure measured without the control (0.660) and with it (0.499)
are not two estimates of one quantity.

---

## The read-out: what the graph actually suggests

`repurposing_candidates.json`, `make repurposing-candidates`. **Not a
protocol** — it points the scorer that survived the registered controls at the
question the graph exists for and prints compound names, with the nearest known
active named beside each one, because a similarity score with no neighbour
named is unauditable.

| virus | verdict | known actives | candidates | resemblance floor | failures scored | rows the floor cost | rows shown |
|:---|:---|---:|---:|:---|---:|---:|---:|
| DENV-1 | clears | 195 | 24,883 | 0.5246 | 123 | 5 | 40 |
| DENV-2 | clears | 995 | 21,163 | 0.4267 | 3,043 | 0 | 45 |
| DENV-3 | clears | 213 | 24,893 | 0.5156 | 95 | 15 | 30 |
| DENV-4 | below | 153 | 24,939 | 0.4143 | 109 | 0 | 0 |
| HCoV-229E | clears | 119 | 24,994 | 0.6664 | 88 | 20 | 25 |
| HCoV-OC43 | weak | 162 | 24,969 | 0.7264 | 70 | 25 | 20 |
| JEV | clears | 19 | 25,155 | 0.2222 | 27 | 0 | 45 |
| MERS-CoV | weak | 199 | 23,044 | 0.2115 | 1,989 | 0 | 45 |
| SARS-CoV | weak | 232 | 23,386 | 0.2316 | 1,610 | 0 | 45 |
| SARS-CoV-2 | clears | 1,768 | 12,420 | 0.2667 | 11,048 | 0 | 45 |
| TBEV | clears | 72 | 25,064 | 0.5968 | 65 | 25 | 20 |
| West Nile | clears | 261 | 24,629 | 0.5634 | 311 | 12 | 33 |
| Yellow fever | clears | 135 | 24,184 | 0.3008 | 882 | 0 | 45 |
| Zika | clears | 748 | 21,876 | 0.5195 | 2,577 | 0 | 45 |

483 rows in total, top 15 per phase tier per cell. A
cell whose verdict is `below` or `unknown` prints **no list**, because a caveat
beside fifteen named drugs loses to the drugs.

The chemistry filter, over 25,243 molecules:
23,005 annotated (3,121 of them by InChIKey,
carrying no ChEMBL id in the release), **2,238 unchecked** —
unchecked is not clean. 42 compounds are diagnostic, contrast
or radiotracer agents and are **dropped and named**
(17 of them resolvable only by structure, which is
where the dyes are). 205 antineoplastics and
219 withdrawn compounds are **labelled and kept**. The
non-therapeutic ATC classes are `S01J` (ophthalmological diagnostic agents), `V04` (diagnostic agents), `V07` (all other non-therapeutic products), `V08` (contrast media), `V09` (diagnostic radiopharmaceuticals), and the rule is **every** class,
not any: TOLBUTAMIDE is a diagnostic *and* a sulfonylurea.

Each surviving row also carries where it sits in that cell's own
**failed-resemblance distribution** — the max-Tanimoto its screened-and-
inactive compounds reach — because a bare Tanimoto cannot be compared across
cells. The floor is the 50% line on it, and the floors in the table above span
a factor of three.

---

## What is licensed, derived from both runs at once

| virus | cross-sectional | AUC | prospective | AUC | shortlist licensed | contradicted |
|:---|:---|:---|:---|:---|:---:|:---:|
| DENV-1 | clears | 0.5967 | clears | 0.8025 | **yes** | no |
| DENV-2 | clears | 0.6921 | clears | 0.7518 | **yes** | no |
| DENV-3 | clears | 0.6653 | clears | 0.7300 | **yes** | no |
| DENV-4 | below | 0.4792 | clears | 0.7896 | no | no |
| HCoV-229E | clears | 0.6758 | clears | 0.7136 | **yes** | no |
| HCoV-OC43 | weak | 0.5903 | no_verdict | 0.4348 | no | no |
| JEV | clears | 0.7885 | not_evaluable | — | **yes** | no |
| MERS-CoV | weak | 0.5045 | not_evaluable | — | no | no |
| SARS-CoV | weak | 0.5763 | not_evaluable | — | no | no |
| SARS-CoV-2 | clears | 0.7803 | clears | 0.6209 | **yes** | no |
| TBEV | clears | 0.7389 | not_evaluable | — | **yes** | no |
| West Nile | clears | 0.7079 | below | 0.3532 | no | **yes** |
| Yellow fever | clears | 0.8117 | no_verdict | 0.3722 | **yes** | no |
| Zika | clears | 0.7032 | below | 0.4218 | no | **yes** |

**8 viruses are licensed for a shortlist**:
DENV-1, DENV-2, DENV-3, HCoV-229E, JEV, SARS-CoV-2, TBEV, Yellow fever. **West Nile, Zika** clear cross-sectionally and
anti-rank prospectively, and are refused by name — a shortlist ranks compounds
nobody has screened, which is the prospective act those cells have been
measured getting backwards.
