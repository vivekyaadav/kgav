# 01 — Overview

## What this is

A knowledge graph of 78,649 nodes and 396,223 edges over
16 viruses from Coronaviridae and Flaviviridae, built to answer one
question: **can a drug already in the clinic be repurposed against one of
these viruses, and does the graph provide any evidence for the answer?**

Three things are in the repository, and they are different kinds of object:

1. **The graph.** Assembled from 12 of the
   13 declared ingest layers, validated against a frozen schema,
   and stamped with content hashes so a result can never be attributed to a
   graph that has since moved. Absent from this release:
   `curated` — nothing matches data/raw/curated_screens/*.tsv.
2. **The evaluation.** Pre-registered protocols, each with its correction
   family fixed before the p-values existed, each writing a stamped artifact.
   Where a protocol refused to run, the refusal is the artifact.
3. **The agent.** A question-answering layer whose refusals are
   keyword-deterministic rather than model-decided, because the one thing the
   graph must not do is answer the question it cannot answer fluently.

## The findings, in one paragraph each

**Chemistry separates actives from measured-inactives, cross-sectionally,
in most cells.** Under the sternest registered control — each candidate's own
`k_cell` nearest known actives withheld — the max-Tanimoto scorer clears its
corrected threshold in **10 of 14** cells
(DENV-1, DENV-2, DENV-3, HCoV-229E, JEV, SARS-CoV-2, TBEV, West Nile, Yellow fever, Zika). Cells whose interval excludes chance on the **wrong** side:
none. This is the project's most robust positive result and it is a
statement about chemistry, not about the graph: the graph's own reasoning
channels do not reach it.

**It also transfers forward, in six cells of ten.** Trained only on activity
published at or before 2021 and scored on compounds measured
after, the same statistic clears in **6 of
10** cells (DENV-1, DENV-2, DENV-3, DENV-4, HCoV-229E, SARS-CoV-2). In
**2** it clears *below* chance — West Nile, Zika — which
is a measured anti-ranking and not a null result.

**The graph's reasoning channels are at chance.** Under the registered
publication-disjoint protocol on the one plate-matched population,
0 of the 8 metapath channels has an
interval excluding 0.5 (none). Two scorers do clear
— the ensemble at 0.5666 [0.5044, 0.6289] and the **degree
baseline** at 0.5700 [0.5077, 0.6323] — and the ensemble's margin is
not ranking: its coverage asymmetry alone scores 0.5578 of that
0.5666, reaching 42.7% of the positives, against
degree's floor of exactly 0.5 at 100.0% coverage. The
correction family was fixed in advance; across 61 cells the only
survivor at the readable-cell family is SARS-CoV-2 M7.

**Clinical labels go the other way.** Against 271
`IN_TRIAL_FOR` and 1 `TREATS` edges for
SARS-CoV-2/MONDO:0100096, the structure scorer reads
0.4087 [0.3734, 0.4440] — clearing
**below** chance — while degree reads
0.7751 [0.7385, 0.8116]. What predicts whether a
compound reached a trial is how well studied it is.

**The host-factor overlap is a null with its power stated.** Over
44 tests at 100,000 degree-matched permutations,
nothing survives correction, and 8 viruses have no cell at
all (DENV-1, DENV-3, DENV-4, HCoV-HKU1, JEV, SARS-CoV, TBEV, Yellow fever).

## What follows from that

A shortlist is a **prospective** act dressed as a cross-sectional one, so it
is licensed only where the cross-sectional control clears **and** the
prospective control does not contradict it. On this release that is
**8 viruses**: DENV-1, DENV-2, DENV-3, HCoV-229E, JEV, SARS-CoV-2, TBEV, Yellow fever. Two clear
cross-sectionally and anti-rank prospectively — West Nile, Zika — and
the read-out refuses them by name.

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

## What this is not

- **Not a predictor of activity.** Every row of every shortlist is the claim
  "this resembles that known active, at this Tanimoto", with the known active
  named so a reader can judge.
- **Not a hit rate.** AUC is a statement about a whole ordering; a screen only
  ever sees the top of a list. `make hit-rate` measures the top, and on the
  one assay-matched cell the precision of the head falls to zero once each
  candidate's nearest known actives are withheld.
- **Not a mechanism.** Where the graph reaches a mechanism it says so and
  attaches the cell-context caveat that made chloroquine look effective in
  2020.
