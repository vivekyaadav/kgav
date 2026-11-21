# 08 — Every measured result

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

#### Every registered cell beside its own degree comparator

The headline claim is not "the scorer clears" but "the scorer clears **where
degree does not**", and those are different counts. All
17 readable cells, in three groups, so that a table of the
first group alone — the failure this project is built against — is not what
a reader gets:

<!-- TABLE: scorer_vs_degree — generated by scripts/write_doc2.py and
     re-derived independently by tests/test_doc2_claims.py. Do not hand-edit:
     the marker is how that test finds the table, and the rows are a
     function of data/results/. -->
| virus | population | structure AUC | clears | degree AUC | clears | group |
|:---|:---|:---|:---:|:---|:---|:---:|
| Yellow fever | ChEMBL | 0.8117 | **yes** | 0.5079 | no | clears, degree does not |
| JEV | ChEMBL | 0.7885 | **yes** | 0.5263 | no | clears, degree does not |
| TBEV | ChEMBL | 0.7389 | **yes** | 0.4753 | no | clears, degree does not |
| West Nile | ChEMBL | 0.7079 | **yes** | 0.5319 | no | clears, degree does not |
| Zika | ChEMBL | 0.7032 | **yes** | 0.5342 | no | clears, degree does not |
| DENV-3 | ChEMBL | 0.6653 | **yes** | 0.4970 | no | clears, degree does not |
| DENV-1 | ChEMBL | 0.5967 | **yes** | 0.5251 | no | clears, degree does not |
| SARS-CoV-2 | plate CPE | 0.5834 | **yes** | 0.5700 | no | clears, degree does not |
| SARS-CoV-2 | ChEMBL | 0.7803 | **yes** | 0.5535 | yes | clears, degree clears too |
| DENV-2 | ChEMBL | 0.6921 | **yes** | 0.5638 | yes | clears, degree clears too |
| HCoV-229E | ChEMBL | 0.6758 | **yes** | 0.6434 | yes | clears, degree clears too |
| SARS-CoV | plate entry | 0.6224 | **yes** | 0.6070 | yes | clears, degree clears too |
| HCoV-OC43 | ChEMBL | 0.5903 | no | 0.6222 | yes | does not clear |
| SARS-CoV | ChEMBL | 0.5763 | no | 0.6908 | yes | does not clear |
| MERS-CoV | plate entry | 0.5655 | no | 0.6326 | yes | does not clear |
| MERS-CoV | ChEMBL | 0.5045 | no | 0.6533 | yes | does not clear |
| DENV-4 | ChEMBL | 0.4792 | no | 0.5245 | no | does not clear |

The two registered near-neighbour controls — symmetric k = 10 and strength-matched — agree on every cell the claim rests on, and disagree on 1: MERS-CoV (plate entry). In that cell the degree baseline clears and sits above the structure statistic under both readings, so it is not evidence about chemistry either way and is outside the claim in both. The membership of the 8 cells that clear above a non-clearing baseline is **identical under both controls**.

**8 clear above a degree baseline that does not clear.**
4 clear while their baseline clears too, and are therefore
not evidence about chemistry. 5 do not clear, and in
4 of those the degree baseline is both above the
structure statistic and clearing.

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
  cell's AUC with the *k* nearest pre-cutoff actives withheld
  (k = 0, 5, 10, 20, 40; k=0 and k=10 reproduce the published columns).
  **Three of those six are analogy, and three are not** — and a results document that reported
  "6 of 10 clear" without this
  would be true and misleading.

| virus | k=0 | k=5 | k=10 | k=20 | k=40 | clears up to |
|:---|---:|---:|---:|---:|---:|:---|
| DENV-1 | 0.8025\* | 0.6900\* | 0.6181 | 0.4617 | 0.3574 | k≤5 |
| DENV-4 | 0.7896\* | 0.7274\* | 0.6868\* | 0.5158 | 0.4362 | k≤10 |
| DENV-2 | 0.7518\* | 0.7195\* | 0.7076\* | 0.6869\* | 0.6555\* | k≤40 |
| DENV-3 | 0.7300\* | 0.6543\* | 0.5901 | 0.4517 | 0.3680 | k≤5 |
| HCoV-229E | 0.7136\* | 0.7366\* | 0.7105\* | 0.6776\* | 0.6616\* | k≤40 |
| SARS-CoV-2 | 0.6209\* | 0.6346\* | 0.6360\* | 0.6233\* | 0.6143\* | k≤40 |

  `\*` marks a cell clearing the corrected alpha at that exclusion strength.
  **3 of the 6 survive to
  k=40** (DENV-2, HCoV-229E, SARS-CoV-2) — forty of their closest known
  actives removed and the order still holding, which is a chemotype signal.
  **2 stop clearing before k=10** (DENV-1, DENV-3):
  their signal was carried by a handful of near neighbours, which is retrieval
  of a known series — useful, and a narrower claim than the AUC sounds.

  **The headline AUC is close to anti-correlated with robustness:**

      by the published AUC (k=0): DENV-1 > DENV-4 > DENV-2 > DENV-3 > HCoV-229E > SARS-CoV-2
      by robustness (k=40): HCoV-229E > DENV-2 > SARS-CoV-2 > DENV-4 > DENV-3 > DENV-1

  DENV-1 has the highest primary AUC of the
  6 and loses the most (0.8025 to 0.3574). **A cell's
  headline number does not say which kind of result it is.** Registered as H3
  in `registrations/prereg_temporal_structure_scorer.md`, left open there as
  "a pointer, not a finding", and closed by the amendment of 2026-11-21.
  `make prospective-analogy`; it is a post-hoc diagnostic on a registered run
  and is labelled as one.
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

### The uncontrolled view, which is where the discoveries live

`metapath_comparison.json` covers 15 viruses and the channels
M1, M2, M3, M4, M5, M6, M7, M8, COMBINED, degree. **It applies none of the controls above**: it runs on
unfiltered `build_labels` output, and `lift_pvalue` reports the exact
hypergeometric tail — the probability of that many hits under a *random
ranking*. That is all it reports. It cannot see a confound that inflates the
scorer and the labels together, which is the confound this whole project is
about. The table is printed because it is what a reader would have produced
without the rest of this repository, and the reasons to distrust it are
visible in its own columns.

<!-- TABLE: metapath_comparison — generated by scripts/write_doc2.py and
     re-derived independently by tests/test_evaluation_table.py. Do not hand-edit:
     the marker is how that test finds the table, and the rows are a
     function of data/results/. -->
| virus | best readable channel | lift | p | significant | degree | p | baseline verdict |
|:---|:---|---:|---:|:---|---:|---:|:---|
| DENV-1 | — | — | — | — | 2.59 | 0.1809 | — |
| DENV-2 | — | — | — | — | 0.76 | 0.7597 | BASELINE NOT DETERMINED |
| DENV-3 | — | — | — | — | 2.37 | 0.2067 | — |
| DENV-4 | — | — | — | — | 3.30 | 0.1233 | — |
| HCoV-229E | — | — | — | — | 4.24 | 0.0810 | — |
| HCoV-NL63 | — | — | — | — | 0.00 | 1.0000 | — |
| HCoV-OC43 | — | — | — | — | 1.56 | 0.4754 | — |
| JEV | — | — | — | — | 0.00 | 1.0000 | — |
| MERS-CoV | M4 | 3.58 | 0.0122 | M4 0.0122 | 10.15 | <.0001 | BASELINE NOT DETERMINED |
| SARS-CoV | M3 | 2.30 | 0.0107 | M3 0.0107 | 11.97 | <.0001 | CONFOUNDED |
| SARS-CoV-2 | M4 | 2.30 | 0.0449 | M7 0.0008 | 0.86 | 0.7096 | BASELINE AT CHANCE |
| TBEV | — | — | — | — | 0.00 | 1.0000 | — |
| West Nile | — | — | — | — | 0.97 | 0.6470 | — |
| Yellow fever | — | — | — | — | 0.00 | 1.0000 | — |
| Zika | — | — | — | — | 0.67 | 0.8004 | — |

**Twelve of fifteen viruses have no readable channel at all** (DENV-1, DENV-2, DENV-3, DENV-4, HCoV-229E, HCoV-NL63, HCoV-OC43, JEV, TBEV, West Nile, Yellow fever, Zika), and that is the first
thing this table says. Of the whole family,
thirteen cells are both `trustworthy` and `hits_determined`, and **Four of them are nominally significant** at an
uncorrected 0.05: SARS-CoV-2 M7, lift 1.30, p=0.000841, SARS-CoV M3, lift 2.30, p=0.010747, MERS-CoV M4, lift 3.58, p=0.012164, SARS-CoV-2 M4, lift 2.30, p=0.044856.

Two prose figures beside the table, because they are the ones a reader
quotes. SARS-CoV-2 M4's lift is 2.69x a degree baseline that is itself at chance here — that cell's
**degree lift of **0.86**** is indistinguishable from random, which
is the only condition under which a channel's p means what it appears to
mean. On SARS-CoV-2, degree lift 0.86 (p=0.7096) is indistinguishable from random, so M7's 1.30 (p=0.0008) is not explained by study volume. That cell's **degree lift of **0.86**** is
the reason its channel row can be read at all. A row
is unreadable for one of two recorded reasons, and on this release both are
severe.

**`hits_determined: False`** means the published hit count is one draw from a
tie block straddling the cut, so the lift computed from it is one draw too.
8 channel–virus cells are in that state; the ones with the
smallest p-values:

| virus | channel | nominal lift | p | hits range | positives |
|:---|:---|---:|---:|:---|---:|
| DENV-2 | M5 | 2.04 | 0.030564 | [6, 8] | 8 |
| Zika | M1 | 1.37 | 0.085960 | [7, 18] | 18 |
| DENV-2 | M1 | 1.13 | 0.238741 | [0, 32] | 32 |
| DENV-2 | M3 | 1.51 | 0.437086 | [1, 2] | 2 |
| SARS-CoV | M1 | 1.12 | 0.446930 | [2, 15] | 15 |
| HCoV-229E | M3 | 0.97 | 0.766519 | [0, 2] | 2 |

A hit count whose range runs from zero to every positive in the cell is not a
measurement.

Note SARS-CoV especially. SARS-CoV M1's Hits@100 is [2, 15] of 15 positives and its lift is [0.32, 2.39], which contains 1.00 — no enrichment at all is inside the range.

**NO POINT ESTIMATE IS PUBLISHED FOR THAT CELL**, and the reason is worth
stating: writing the release in a canonical order changed nothing about the
graph and moved the figure from 2.08 to 0.83, because `rank()` sorts on
score with no tiebreak and Python's sort is stable — so which of the tied
compounds fall inside Hits@100 follows the order `edges.jsonl` happened to
be written in. A document quoting one draw from that range, and saying two
clauses later that the range contains 1.00, is publishing an artifact of
file ordering. Several of the largest-looking results in this project are a tie
block being reported as one.

**`trustworthy: False`** means the expectation at the cut is a fraction of a
hit, so a single hit produces an enormous lift. 22 cells are in
that state, and every eye-catching double-digit lift here is of that shape:

| virus | channel | lift | hits expected at the cut |
|:---|:---|---:|---:|
| West Nile | M4 | 10.05 | 0.30 |
| HCoV-OC43 | M5 | 6.14 | 0.16 |
| HCoV-OC43 | M4 | 5.60 | 0.18 |
| Zika | M5 | 4.63 | 0.65 |
| HCoV-229E | M4 | 4.38 | 0.23 |
| DENV-2 | M4 | 4.15 | 0.96 |

**Seven of fifteen viruses have no positive reachable by any metapath** (`COMBINED.pos == 0`):
DENV-1, DENV-3, DENV-4, HCoV-NL63, JEV, TBEV, Yellow fever. That is not a weak result, it is
the absence of a measurement, and it comes in two shapes a reader must not
conflate: **2 reach no compound at all**
(DENV-4 and Yellow fever), so nothing was ranked; the rest rank compounds and none
of them carries a label, so a ranking exists, is complete, and cannot be
scored. A ranking over those is not a weak prediction — it is a prediction
nothing in this graph can mark right or wrong.

That distinction is load-bearing outside this document: the agent's no-route
guard refuses exactly the viruses with an empty pool, and answers for the
others. Edge counts are what that guard reads, and edge counts cannot see
this.

### What the baseline says, per (run, virus) pair

The verdict column above is `confound_warning`'s, and it is a statement
about the **baseline** rather than about any channel. Across every run in
the project:

    CONFOUNDED
        Two (run, virus) pairs — the degree baseline is itself significant, so a channel beating random says nothing about mechanism.
        cross-sectional / SARS-CoV, temporal include / SARS-CoV.
    BASELINE AT CHANCE
        One (run, virus) pair — degree is indistinguishable from random, which is the only condition under which a channel's p means what it appears to mean.
        cross-sectional / SARS-CoV-2.
    BASELINE NOT DETERMINED
        Three (run, virus) pairs — the degree row's own hit count is one draw from a tie block, so no comparison against it can be read.
        cross-sectional / DENV-2, cross-sectional / MERS-CoV, temporal include / SARS-CoV-2.
    (nothing)
        every other pair: the comparison is readable and no verdict applies.

### Every cell in the project whose interval excludes chance

The one-line version of every table here:

> **No individual reasoning channel, on any protocol, on any virus, has an
> AUC whose interval excludes chance.** Every scorer that does is either the
> degree baseline or the COMBINED ensemble of all eight channels. There are
> **2** such cells on this release and they are all of
> those two scorers, on SARS-CoV-2:
>
    cross-sectional          SARS-CoV-2   COMBINED 0.567 [0.504, 0.629]  floor 0.5578  87% of the margin is reach
    cross-sectional          SARS-CoV-2   degree 0.570 [0.508, 0.632]  floor 0.5000  0% of the margin is reach
>
> The last column is how much of each margin above chance is the imputation
> floor rather than the ranking. Degree has no floor at all, because it is
> the only scorer here that reaches every compound in both classes.

Five figures for M1 — produced by five protocols, each of which failed a
control the next one applied — are reported in the sequence above rather
than withdrawn, because the sequence is the result.

### The correction families, fixed in advance

The smallest p among the readable cells is **at p=0.000841**, and
the table below judges that one — the **smallest**, not the only one, since
Four of them are nominally significant uncorrected. Four families were declared before
the p-values existed, and the family is reported rather than chosen:

| family | tests | corrected alpha | the smallest readable cell, at p=0.000841 | how the size is arrived at | survivors |
|:---|---:|:---|:---|:---|:---|
| readable cells only | 13 | 0.00385 | clears, by a factor of 5 | cells that are both `trustworthy` and `hits_determined` | SARS-CoV-2 M7 |
| every cell computed | 61 | 0.00082 | does not clear | every (virus, channel) the run produced | none |
| every cell planned | 120 | 0.00042 | does not clear | (15 viruses x 8 metapaths) | none |
| planned, plus COMBINED and degree | 150 | 0.00033 | does not clear | (15 viruses x 10 scorers) | none |

the family is reported, not chosen: this project fixes a correction family in a pre-registration before the p-values are in hand, and these were not

**The family is reported, not chosen.** Printing four and then reading the
smallest as the answer is exactly the ordering this project's rule forbids:
which family applies **belongs in the prereg**, before the p-values exist,
and the table above is a sensitivity on a decision already made — not a menu.

### The metapath channels under the same split, which is a different study

`temporal_2021_include.json`. The structure scorer above and the channels
below split the labels by time in the same way and differ in what they rank,
so neither bears on the other — which is why both are printed rather than
one being reconciled into the other.

<!-- TABLE: temporal_include_auc — generated by scripts/write_doc2.py and
     re-derived independently by tests/test_temporal_table.py. `~ties`
     marks a row whose hit count is one draw from a tie block at the cut;
     `(!)` marks one whose expectation at the cut is a fraction of a hit.
     Do not hand-edit: the marker is how that test
     finds the table, and the rows are a function of data/results/. -->
| scorer | pool | pos | lift | p | AUC [95% CI] |  |
|:---|---:|---:|---:|---:|:---|:---|
| M1 | 280 | 8 | 0.35 | 0.9724 | 0.501 [0.479, 0.524] | ~ties |
| M2 | 1,182 | 24 | 2.46 | 0.0456 | 0.506 [0.483, 0.528] |  |
| M3 | 1,144 | 21 | 2.18 | 0.1036 | 0.505 [0.482, 0.528] |  |
| M4 | 3,772 | 47 | 3.21 | 0.0347 | 0.505 [0.482, 0.527] |  |
| M6 | 234 | 2 | 1.17 | 0.6731 | 0.499 [0.477, 0.522] | (!) |
| M7 | 61 | 0 | 0.00 | — | 0.500 [0.477, 0.523] | (!) |
| M8 | 24 | 0 | 0.00 | — | 0.500 [0.477, 0.523] | (!) |
| COMBINED | 4,323 | 54 | 3.20 | 0.0353 | 0.504 [0.482, 0.527] |  |
| degree | 5,756 | 80 | 2.16 | 0.1617 | 0.507 [0.484, 0.530] | ~ties |

**0 of these scorers has an AUC interval
excluding chance** (none).

**The AUC is computed over the whole labelled test set, not over what a
scorer reaches** — `evaluate_against_negatives` imputes the floor to
everything a scorer cannot reach, so the *n* is the same for every row and
is **1,432 and 1,093**. That matters because the bar follows from it: the
smallest AUC whose 95% interval excludes chance at that sample size is
**0.523**, and no row above reaches it.

**And the decomposition leaves nothing behind.** `temporal_split.py` records
the coverages — as `auc_coverage_pos` and `auc_coverage_neg`, which three
documents once said did not exist — so the imputation floor costs nothing
extra to compute. **7 of these cells equal their own
floor to four decimal places** (M1, M2, M3, M4, M6, M7, M8), and across all
41 temporal AUC cells in every policy,
**31 do**, with the largest deviation from a
floor anywhere being 0.0191. A scorer on this split reaches a
few per cent of the test positives and a smaller fraction of the test
negatives; the gap between those two numbers is the whole figure, and the
ranking inside the reached set has almost no room to show itself and
demonstrably uses none of it.

That is a stronger negative than "no channel separates". It says the
tie-safe AUC **cannot** separate on this protocol: the statistic is
dominated by imputation before any scorer is applied, so a future channel
could rank its reached set perfectly and still read near chance.

---

### The one readable cell that is nominally significant, under both weightings

`degree_weighting_sensitivity.json`. The step weight in `walk_metapath` uses
each entered node's **global** degree, which includes the held-out label
edges of a labelled compound. Whether that inflates the one readable
significant cell is a question about the scorer, so it is measured rather
than argued:

| arm | what it changes |
|:---|:---|
| `A_as_published` | g.degree for the step weight; hub set from it. What data/results holds. |
| `B_weight_only` | degree_unlabelled for the step weight, hub set HELD AT A's. The isolated change. |
| `C_weight_and_hubs` | degree_unlabelled for both. |

**SARS-CoV-2 M7 reads lift 1.30 at p=0.000841 as published
(arm A, which is what `data/results` holds), and lift 1.28 at
p=0.001788 with the weighting changed (arm B).** Arm B still clears the
readable-cells alpha, by a factor of 2. The hub set is
identical across the arms — 62 nodes — which is what makes arm B
an isolated change rather than two at once: every held-degree protein is
viral and the hub cut is class-scoped to host proteins.

**4 of the 61 cells move at all and
0 cross their alpha**, so no cell crosses the threshold in
either direction. That denominator is what keeps the finding small: a
document stating the claim cell's own move without it would read as a
general problem with the scorer.

This artifact **registers nothing**. Switching the weighting is a decision
that needs a registration of its own, and this is a sensitivity on a figure
that already exists.

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

`pathway_landscape.json`. 100,000 permutations,
`degree_basis = study`, 1,944 proteins in the
drugged proteome. Each row is tested in both directions, so
`n_tests = 44` over 22 rows and the Bonferroni
threshold is `alpha = 0.00114`.

**0 of 22** rows clear it, in either direction.

<!-- TABLE: pathway_landscape — generated by scripts/write_doc2.py and
     re-derived independently by tests/test_landscape_table.py. Do not hand-edit:
     the marker is how that test finds the table, and the rows are a
     function of data/results/. -->
| channel | virus | dependency factors | observed overlap | expected | naive p | degree-matched p | p (depletion) | verdict |
|:---|:---|---:|---:|---:|---:|---:|---:|:---|
| taxa | DENV-2 | 50 | 2 | 6.28 | 0.9581 | 0.9907 | 0.0394 | — |
| taxa | HCoV-229E | 1,367 | 128 | 135.13 | 0.5963 | 0.7643 | 0.2665 | — |
| taxa | HCoV-NL63 | 57 | 10 | 7.39 | 0.0412 | 0.1963 | 0.8881 | — |
| taxa | HCoV-OC43 | 657 | 63 | 65.16 | 0.4964 | 0.6336 | 0.4176 | — |
| taxa | MERS-CoV | 678 | 59 | 66.41 | 0.7883 | 0.8530 | 0.1809 | — |
| taxa | SARS-CoV-2 | 3,460 | 383 | 358.31 | 0.0005 | 0.0665 | 0.9412 | — |
| taxa | West Nile | 18 | 0 | 2.51 | 1.0000 | 1.0000 | 0.0663 | — |
| taxa | Zika | 815 | 76 | 80.01 | 0.5956 | 0.7046 | 0.3396 | — |
| interactome | DENV-1 | 53 | 11 | 7.13 | 0.0102 | 0.0907 | 0.9547 | — |
| interactome | DENV-2 | 271 | 40 | 35.45 | 0.0035 | 0.2259 | 0.8240 | — |
| interactome | DENV-3 | 33 | 6 | 4.43 | 0.0884 | 0.2775 | 0.8556 | — |
| interactome | HCoV-229E | 397 | 34 | 50.24 | 0.7682 | 0.9969 | 0.0052 | — |
| interactome | HCoV-HKU1 | 279 | 25 | 36.72 | 0.6561 | 0.9896 | 0.0174 | — |
| interactome | HCoV-NL63 | 311 | 32 | 40.69 | 0.3486 | 0.9466 | 0.0769 | — |
| interactome | HCoV-OC43 | 289 | 28 | 38.25 | 0.4909 | 0.9754 | 0.0377 | — |
| interactome | JEV | 7 | 1 | 0.83 | 0.5037 | 0.5880 | 0.8045 | — |
| interactome | MERS-CoV | 622 | 64 | 79.86 | 0.2728 | 0.9817 | 0.0253 | — |
| interactome | SARS-CoV | 1,238 | 167 | 149.42 | 0.0000 | 0.0588 | 0.9509 | — |
| interactome | SARS-CoV-2 | 3,494 | 438 | 426.36 | 0.0000 | 0.2451 | 0.7734 | — |
| interactome | TBEV | 35 | 4 | 4.62 | 0.4305 | 0.6989 | 0.4989 | — |
| interactome | West Nile | 6 | 4 | 0.91 | 0.0011 | 0.0058 | 0.9996 | — |
| interactome | Zika | 353 | 47 | 45.00 | 0.0117 | 0.3968 | 0.6646 | — |

**8 viruses have no cell in the host-factor channel**: DENV-1, DENV-3, DENV-4, HCoV-HKU1, JEV, SARS-CoV, TBEV, Yellow fever, because a virus with no `HOST_FACTOR_FOR` edge has no dependency-factor set to overlap. the interactome channel tested 14 of 16 and the host-factor channel 8 of 16; 2 have no interactome cell (DENV-4, Yellow fever), and **2 are missing from both** (DENV-4, Yellow fever) — the same pair the agent refuses to rank, arrived at independently. A null
result is only a finding if the reader knows what would have survived, so
[09-limits.md](09-limits.md) carries the minimum detectable effect per cell —
failing to detect a 2× excess is not evidence against a 2× excess when the
test needed 7×.

**The uncorrected view is where the discoveries live, and the gap is the
finding.** 7 rows are significant at an uncorrected
0.05 against a *uniform* null and are nothing against a degree-matched one:
SARS-CoV-2's interactome overlap is 0.0000 naive and 0.2451 matched; SARS-CoV's interactome overlap is 0.0000 naive and 0.0588 matched; SARS-CoV-2's taxa overlap is 0.0005 naive and 0.0665 matched; DENV-2's interactome overlap is 0.0035 naive and 0.2259 matched; DENV-1's interactome overlap is 0.0102 naive and 0.0907 matched; Zika's interactome overlap is 0.0117 naive and 0.3968 matched; HCoV-NL63's taxa overlap is 0.0412 naive and 0.1963 matched. Neither gap opens because the data changed. The naive
p asks whether these proteins are over-represented in pathways; the matched
one asks whether they are over-represented *compared with proteins as
well-connected as they are*. A third-of-a-percent result becomes nothing,
and so does a one-percent one — the same confound the rest of this document
measures on compounds, measured here on proteins.

`landscape_count_move.json` records why two cells moved between releases, and
that the mechanism accounts for **every** move: a protein already inside the
factor set joined the drugged proteome, with no edge in its own channel
changing.

---

## The circularity confound, measured

`label_overlap.json`. A paper reporting an Mpro IC50 reports a cell-based
EC50 in the same table, and ChEMBL files them as two records: the channel
traverses a measurement taken alongside the answer.

Of 943 compounds carrying both an `INHIBITS` edge and an antiviral
label, 880 have a publication on both sides and
**779 -- 88.5% -- share a PMID**
(88.5% shared PMIDs (779/880)). `INHIBITS`
reaches 22 distinct viral proteins in all, and two of them
dominate it.

**THE FIGURE IS CONDITIONAL, and the condition is the hole.** It is the
share among compounds where *both* sides cite a publication, not among all
labelled compounds — quoting it without that qualifier is recorded in the
script's own comment as a past mistake.

**58 compounds (6.2%) carry no
publication on the `INHIBITS` side at all**, and those are outside the
control rather than cleared by it: `Graph.load(drop_pubs_shared_with=...)`
withholds an edge that *shares* a PMID with its label, so an edge with no
PMID is never a candidate for withholding. At 6.2% that
does not move row 4 of the sequence, which is why the figure is stated
rather than used to qualify it — but a virus can look clean under this
control by having nothing to be disjoint from.

One line: evidence co-reported with labels is 88.5% shared PMIDs (779/880), and
a further 58 compounds have no `INHIBITS` publication
and are outside the control rather than cleared by it.

An earlier version of this section quoted **94.9% shared PMIDs**, over a
different denominator. The figure is retired and kept here so the
correction is readable; everything above is the current run's.

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

## M8 is motif-only on this release, which is weaker than the channel sounds

<!-- FIGURE: m8_fold_strength -->
`folds.gate()` makes two kinds of comparison and labels what each
establishes. NUCLEOPHILE mode requires an identical catalytic residue and is
marked `strong`; MOTIF mode has no nucleophile at either end, rests on shared
catalytic-domain entries, and the gate marks it `weak` and refuses to promote
it — *polymerases share motif C across all three families, which does not
discriminate the way Cys against Ser does*. Of the **36** allowed
cross-family pairs, **20** are `nucleophile/strong` and
**16** are `motif/weak` — but the strong ones are protease pairs
whose comparator chains the spine does not build, so **0**
of them reach the release. Every one of the **4**
`FOLD_SIMILAR_TO` edges M8 can walk is `motif/weak`.

M8 deliberately does **not** constrain on `fold_strength`. Doing so would not
tighten the channel, it would remove it — which is a claim about whether a
shared polymerase motif is evidence at all, not a schema correction. The
channel's own AUC row above already reports it as ranking nothing.

---

## What crosses the cutoff, and what the overlap policy drops

<!-- FIGURE: overlap_census -->
A compound measured on both sides of 2021 is the thing that
makes a temporal protocol not temporal. Counted on this release, the overlap
is **76** pairs; of those, **51** agree on polarity
across the cutoff and **25 disagree** — the same compound reading
active on one side and inactive on the other, which is a real measurement
and not a bookkeeping error. The policy `drop-from-test` therefore
removes them from the TEST side, recorded as 12 positive and
39 negative labels dropped. (A pair holding both polarities
on both sides is dropped once on each side, so the two totals need not
agree with the pair count and must not be conflated with it.)

The withdrawn version of this paragraph said the overlap was evenly split,
"19/19 by source", and it was neither 19 nor even. The figure is retired and
quoted here so the correction is readable.

---

## What the dose-series policy costs, both ways

<!-- FIGURE: dose_series_cost -->
ChEMBL's curators flag a potency read outside the concentrations actually
tested, and this release carries that flag as a qualifier rather than
dropping the edge. **4,839** label edges carry it,
**23** of them positives. Excluding them moves
positives from **5,276** to **5,253** and
negatives from **22,044** to **17,353** — a tenth
of the negatives, which is why the choice is recorded in every results file
rather than left implicit.

The policy is `include` and the alternative is one argument
away (``include`, `exclude``). The reasoning is the same as for every
other weak-evidence decision here: **excluding the evidence deletes it,
while recording the flag makes the weakness legible** and leaves the
exclusion to whoever wants it. What is lost either way: a curator's flag is
a fact about one source row and an edge collapses rows, so a flagged edge
means *at least one* contributing row was extrapolated, not that all were.

---

## Seven measured causes

Why the graph's own channels read the way they do, each with the script
that re-derives it. **Treat a figure with no named source as unverified**
— that rule is this section's, and a cause added without registering its
script is exactly the figure it forbids.

| Cause | Measurement | Re-derived by |
|---|---|---|
| Evidence co-reported with labels | 88.5% shared PMIDs (779/880); a further 58 compounds have no `INHIBITS` publication and are outside the control | `diagnose_label_overlap.py` |
| M1 was never purely direct-acting | 1,884 of 7,677 protein-target rows (24.5%) described a CELLULAR readout and were filed against a viral protein; now reattributed to the organism | `diagnose_assay_type.py` |
| `INHIBITS` is effectively two proteins | nsp5 and nsp3 dominate it; `INHIBITS` reaches 22 distinct viral proteins in all | `diagnose_label_overlap.py` |
| The graph is two compound populations | 44 of 18,075 compounds carry both `INHIBITS` and `TARGETS`, so almost nothing is reachable by a direct-acting and a host-directed route at once | `diagnose_hop_attrition.py` |
| M8 walks and reaches almost nothing | the layer assembles — 104 `MEMBER_OF_CLASS` and 6 `FOLD_SIMILAR_TO` edges, 2 carrying `catalytic_type_match` — and 34 of the 36 allowed transfers are blocked by comparator chains the spine does not build | `ingest_folds.py` |
| M3 and M5 are untestable prospectively | both PPI layers are wholly post-cutoff by their own dates, so the pre-2021 graph holds no host–host interaction at all | `diagnose_temporal_leak.py` |
| The chemistry source is a medicinal-chemistry corpus | 5 journals supply 87.1% of ChEMBL 37's 93,489 documents, and 19 named repurposing venues supply 2.5% of them — 7 documents in 2020–2021. A congeneric series on one target is the shape it curates; a compound measured on two viral proteases is not | `diagnose_source_scope.py` |

Note which of these are **evaluation** faults rather than graph faults. The
graph validates, provenance is intact, and the host proteins that drugs
target and CRISPR screens call dependencies are a real repurposing surface
that no protocol here has scored against — because almost no labelled
compound reaches it. `diagnose_comparator_cost.py` measures what that costs.

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
