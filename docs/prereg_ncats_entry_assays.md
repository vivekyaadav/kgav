# Pre-registration: the NCATS pseudotyped-entry assays as plate-matched labels

**SIGNED 2026-11-17, BEFORE THE INGEST CODE AND BEFORE ANY AUC.** No score
exists on these labels. `src/kgav/ncats.py` cannot read a loss-of-signal
assay today, so there is no exploratory version of this measurement
anywhere in the repository and none can be produced until the code below is
written.

**What HAS been looked at, stated plainly, because a registration that hides
its reconnaissance is worth nothing.** Six files were downloaded and
counted: row totals, curve-class distributions, efficacy ranges, AC50
medians, and the class sizes each admission rule would produce. Those are
SAMPLE SIZES, not outcomes -- no similarity was computed, no compound was
ranked, and the scorer was not run. The power figures below were derived
from the class sizes alone, which is the one thing a registration is
supposed to do with them.

Date: 2026-11-17. Release at signature: v0.1, 78,500 nodes / 391,727 edges,
`nodes_sha` b0d6fc69d63942c2, `edges_sha` 545fb7d87b0aec08. Code to be
written: per-assay polarity in `src/kgav/ncats.py`, new inputs in
`scripts/ingest_ncats.py`.

## Why this exists

`docs/prereg_structure_scorer.md` registers the NCATS panel as the SECONDARY
population, "the strongest control available and the weakest in power", and
fixes it at **one cell, SARS-CoV-2, 89 positives against 7,343 negatives**.
It reads AUC 0.5870 with a corrected interval of [0.4941, 0.6799] -- it
cannot reject chance, and could not have at any AUC below 0.5769.

That cell is the only place in this project where both label classes come
off the same plates under one protocol. Every other evaluated cell uses
ChEMBL labels, where a positive is a compound somebody published as an
antiviral and a negative is a compound somebody published as inactive, so
the populations differ by more than their activity. The ChEMBL cells read
0.7822 and 0.8304; the one plate-matched cell reads 0.5870 with an interval
crossing chance. **Whether the first number means what it looks like is the
open scientific question in this repository, and 89 positives is why it
stays open.**

The NCATS OpenData portal carries 18 assay files. Three are cell-based
antiviral assays this project does not read.

## Fixed in advance

### Which assays are admitted — **DECISION**

Admitted under `HAS_ANTIVIRAL_ACTIVITY_AGAINST`, because each measures
inhibition of a virus in cells:

| file | virus | taxon |
|---|---|---|
| `SARS-CoV-2_Pseudotyped_Particle_Entry.csv` | SARS-CoV-2 | NCBITaxon:2697049 |
| `SARS-CoV_Pseudotyped_Particle_Entry.csv` | SARS-CoV | NCBITaxon:694009 |
| `MERS_Pseudotyped_Particle_Entry.csv` | MERS-CoV | NCBITaxon:1335626 |

**REFUSED, and the refusal is the load-bearing decision here:** the
enzymatic and binding assays -- `3CL_Enzymatic_Activity`,
`RdRp_Enzymatic_Activity`, `TMPRSS2_Enzymatic_Activity`,
`ACE2_Enzymatic_Activity`, the three `Spike-ACE2_Protein-Protein_
Interaction` files. They measure target engagement. A compound that
inhibits 3CL protease in a tube has not been shown to inhibit a virus, and
labelling it an antiviral positive would inflate the class with compounds
that have never met a virus. They would roughly double the positive count,
which is exactly why the rule is fixed here rather than after seeing how
thin the cells are.

`HEK293_Cell_Line_Toxicity` and `Human_Fibroblast_Toxicity` are refused as
antiviral labels and are not used at all; each assay's OWN counterscreen is
used for cytotoxicity, below. `All_SARS_CoV2_Datasets.csv` is refused as a
pooled re-export of files already admitted individually.

### A pseudotype is a surrogate, and is labelled as one — **DECISION**

A pseudotyped particle carries the spike of the named virus on a different
viral core. Inhibiting its entry is evidence about entry, not about
replication of the authentic virus, and it is weaker evidence than a
cytopathic-effect assay on live virus. These edges therefore carry
`assay_type: cell_based_pseudotype_entry`, DISTINCT from the CPE layer's
`cell_based_antiviral`, so any consumer can include or exclude them and
neither this document nor a later one has to relitigate the question.

### No pooling across assays within a virus — **DECISION**

SARS-CoV-2 will have two plate cells, CPE and entry. They are **not
merged**. The entire value of a plate-matched cell is that both classes come
off the same plates under one protocol over one concentration range, and
pooling two protocols destroys exactly that property -- it is the confound
`build_labels`' docstring warns about, one level down. **The registered
SARS-CoV-2 CPE cell is unchanged and remains the one quoted by
`prereg_structure_scorer.md`.** The entry cell is a new, separate cell.

### Polarity is DERIVED per assay, never assumed — **DECISION**

The CPE assay is gain-of-signal: a compound that protects cells raises the
readout, so activity is POSITIVE efficacy, and `src/kgav/ncats.py` encodes
that. **The three entry assays are loss-of-signal**: inhibiting entry lowers
the reporter, so inhibition appears as NEGATIVE efficacy and a NEGATIVE
curve class. Measured in the staged files before signing: SARS-CoV-2 entry
efficacy runs -130.2 to +101.3 with **1,451 rows at or below -50 against 60
at or above +50**, and its curve classes are -3, -2.4, -2.3, -1.3 where CPE
carries 3, 2.4, 2.2, 1.4.

**Applying the CPE rule to these files yields ZERO actives in all three** --
measured, not supposed -- so a layer that assumed one convention would have
emitted three all-negative cells that validate, count plausibly and encode
the opposite conclusion. That is the failure `scripts/ingest_ncats.py`'s own
header warns about in those words. The polarity is therefore read from the
sign of the fitted curve classes in each file and **asserted per assay in
code, with a test that fails if a file's observed sign disagrees with its
declared one** -- not inferred at runtime, because a file that arrives
empty or truncated would then silently redefine activity.

### The activity call — **DECISION**, and it is the CPE rule with the sign mapped

- **ACTIVE** requires all of: a fitted curve class, `|curve_class2|` in
  {1.1, 1.2, 1.3, 1.4, 2.1, 2.2, 2.3, 2.4}; efficacy at or beyond ±50% in
  the assay's active direction; `ac50` below **10 µM**
  (`labels.DEFAULT_ACTIVE_BELOW_NM` = 10,000 nM, unchanged); a parseable
  SMILES; and the cytotoxicity condition below.
- **INACTIVE** requires curve class exactly 4 and that the assay reached a
  top concentration at or above 10 µM (`ncats.DEFAULT_INACTIVE_AT_NM`,
  unchanged), so an untested compound is never called inactive.
- **REFUSED, counted, and neither class:** `|curve_class2|` of 3 or 5 (an
  incomplete fit -- NCATS itself treats these as inconclusive, and the CPE
  layer already refuses them, 638 rows there); efficacy short of the
  threshold; no SMILES; no concentration.
- **CYTOTOXICITY.** A compound whose counterscreen AC50 is at or below its
  antiviral AC50 is killing cells at the concentration where it "works", and
  is REFUSED rather than called active. Each assay uses its own
  counterscreen: Tox for SARS-CoV-2 entry, Vero E6 for SARS-CoV, Huh7 for
  MERS. Measured cost before signing: it removes 73 of 408, 10 of 107 and
  9 of 82 candidate actives respectively.

### The resulting cells — class sizes fixed before any score

| cell | positives | negatives |
|---|---:|---:|
| SARS-CoV-2, CPE (existing, unchanged) | 89 | 7,343 |
| SARS-CoV-2, pseudotype entry | 335 | 2,857 |
| SARS-CoV, pseudotype entry | 97 | 1,614 |
| MERS-CoV, pseudotype entry | 73 | 2,059 |

These are predictions about the ingest, not about the scorer: if the built
layer does not produce these four cells at these sizes, the ingest is wrong
and this document says so before it is written.

### Correction family and the readability gate — **DECISION**

- **The plate-matched analysis is its own family of four cells, two-sided.**
  alpha = 0.05/4 = **0.012500**, **z = 2.4977**. Two-sided because a
  below-chance AUC on plate labels is a finding and not a null.
- **They are NOT added to the 14-cell ChEMBL family.** That convention is
  already registered in `prereg_structure_scorer.md` -- a planned secondary
  does not penalise the primary -- and this document does not change it.
- `calibration.MIN_REACHED_POSITIVES` = 10 of each class, unchanged. All
  four cells clear it.
- The metric is unchanged: the tie-safe rank-sum AUC with ties at half, the
  Hanley & McNeil interval, and `auc_floor` / `auc_ceiling` / `auc_skill`
  reported for every cell.

## The hazards, quantified before signing

**H1. Three of the four cells are still underpowered, and adding them does
not fix that.** The smallest AUC each cell could distinguish from chance at
alpha = 0.05/4:

| cell | pos/neg | smallest detectable AUC | half-width |
|---|---|---:|---:|
| SARS-CoV-2 CPE | 89 / 7,343 | 0.5769 | 0.0769 |
| **SARS-CoV-2 entry** | 335 / 2,857 | **0.5416** | **0.0416** |
| SARS-CoV entry | 97 / 1,614 | 0.5754 | 0.0754 |
| MERS-CoV entry | 73 / 2,059 | 0.5859 | 0.0859 |

**Only the SARS-CoV-2 entry cell materially improves on what exists.** The
other two buy coverage -- a plate cell for two viruses that had none -- and
not power. This is registered because the obvious way to misreport this work
is to say "the plate-matched evidence tripled" when what tripled is the
number of underpowered cells.

**H2. The ingest moves `nodes_sha` and `edges_sha`, so the ChEMBL analysis
moves too, by a route that is not the labels.** The primary cells are
`build_labels(sources={"infores:chembl"})` and are unaffected by NCATS
labels by construction. They are NOT unaffected by the graph: new compounds
enter, degree changes, and the structure scorer's candidate pool and
publication control both read the graph. **Registered in advance: every
claim-carrying ChEMBL AUC must be shown moved or unmoved by `make
auc-delta` against the pre-ingest release, and a cell that stops carrying
its claim is a result to report, not a regression to tune away.**

**H3. An entry inhibitor is not an antiviral.** Blocking spike-mediated
entry of a pseudotyped particle is one mechanism against one step. A cell
built from it answers "does structural similarity to known entry inhibitors
rank entry inhibitors", which is a narrower question than the CPE cell's.
The `assay_type` qualifier keeps the two separable; no claim may pool them.

**H4. Compound overlap between the assays is unmeasured at signature.** The
three entry files draw on the same three libraries as CPE
(Annotated/Bioactive, Approved Drugs NPC, Anti-infective), so the cells will
share compounds. They are separate cells and no cell pools them, so this
does not confound any single cell -- but **the four AUCs are not
independent, and the family correction above treats them as if they were,
which is conservative for each cell and wrong for any statement about all
four together.** No such statement is registered.

## Registered predictions

Falsifiable, and recorded before the code exists.

- **P1. The ingest produces exactly the four cells at the sizes tabled
  above**, ±0 compounds. A different size means the admission rule as built
  is not the one registered here.
- **P2. The SARS-CoV-2 entry cell's corrected interval is narrower than the
  CPE cell's** -- half-width below 0.05 against 0.0769. This follows from
  the class sizes and is a check on the arithmetic, not on the scorer.
- **P3. Every entry cell reads BELOW its virus's ChEMBL cell.** SARS-CoV-2's
  ChEMBL cells read 0.7822 and 0.8304; the prediction is that the plate
  cells read lower, because the ChEMBL figure is inflated by the literature
  selecting which compounds get published as actives and the plate labels
  cannot be. **If an entry cell matches the ChEMBL figure, the
  publication-selection explanation of the 0.78 is wrong** and this project
  has been wrong about its own main caveat.
- **P4. The SARS-CoV-2 entry cell clears alpha = 0.05/4**, i.e. AUC above
  0.5416. This is the prediction the whole exercise exists to test and the
  one most likely to fail: the CPE cell at 0.5870 did not clear its own
  looser bar.
- **P5. The SARS-CoV and MERS-CoV entry cells do NOT clear.** Their bars are
  0.5754 and 0.5859 and their classes are thin. Registered so that failing
  to clear is recorded as predicted rather than as disappointment.

## What would change the conclusion

- P4 clearing with P3 holding is the result this measurement is for:
  structural similarity ranks plate-generated actives above plate-generated
  inactives, at a lower level than the literature-labelled cells.
- P4 failing while the ChEMBL cells hold at 0.78 is **also a result**, and
  the more consequential one: it would say the ChEMBL figure does not
  survive contact with independently generated labels, and every repurposing
  list in this repository would need that stated beside it.
- P3 failing -- a plate cell matching 0.78 -- would mean the
  publication-selection caveat is wrong and should be withdrawn.

## What this cannot do

It cannot make a claim about any virus other than the three named. It cannot
say a ranked compound works: an entry assay at 10 µM in a pseudotype system
is a screening hit, not a drug. And it cannot repair the ChEMBL cells'
exposure to publication selection -- it can only measure whether a different
kind of label agrees.
