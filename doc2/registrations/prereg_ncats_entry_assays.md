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

`doc2/registrations/prereg_structure_scorer.md` registers the NCATS panel as the SECONDARY
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

---

# Amendment, 2026-11-17: P1 failed, and the reason was worth the registration

**STILL BEFORE ANY AUC.** No score has been computed on any of these cells.
This section is appended and replaces nothing above, because a registration
edited to match what happened later is worth nothing.

## What P1 said and what the built code gives

P1 predicted the ingest would reproduce the four cells at the tabled sizes,
"±0 compounds", and said that a different size means the ingest is wrong.
It was half right and the half it got wrong was mine, not the code's.

| cell | registered | built | positives | negatives |
|---|---|---|---:|---:|
| SARS-CoV-2 entry | 335 / 2,857 | 336 / 3,326 | **+1** | **+469** |
| SARS-CoV entry | 97 / 1,614 | 97 / 1,778 | **0** | **+164** |
| MERS-CoV entry | 73 / 2,059 | 73 / 2,098 | **0** | **+39** |

**The positives hold.** Two cells match exactly. The +1 on SARS-CoV-2 is a
duplicate `sample_id`: 44 compounds appear twice in that file, one pair at
3.16 µM and 0.001 µM, and the reconnaissance filter and the layer resolved
the pair differently. Replicate handling is the existing layer's
(`replicate_disagreement`), not a new decision.

**The negatives were wrong in my table, for two reasons, both of them mine.**

1. The prose above registers that an INACTIVE "requires that the assay
   reached a top concentration at or above 10 µM". **My table did not apply
   that rule when it was computed** -- it counted every curve-class-4 row
   with a structure. The layer applies it and refuses 52, 38 and 48 rows
   respectively. The rule was registered correctly and the figure beside it
   was not derived from it, which is this repository's most frequent defect
   appearing inside its own pre-registration.
2. **The larger half, and a real finding.** My table excluded compounds the
   assay FITTED as inhibitors whose EC50 lands past the 10 µM threshold. The
   pipeline did not exclude them -- it counted them as NEGATIVES.

## The finding: a confirmed inhibitor was being counted as a non-inhibitor

`ncats.classify_row` returned "active" for any fitted hit in the right
direction. `labels.classify` then read the project's threshold off the
emitted edge and returned "inactive". So a compound with a real
concentration-response curve at 12-25 µM was a NEGATIVE in every cell this
layer feeds.

**It also broke the contract `classify_row`'s own docstring states** -- "the
qualifiers are written so that labels.classify() reaches the SAME verdict
from the emitted edge alone ... that symmetry is the contract" -- and
nothing checked it. Measured on the existing CPE layer: **17 edges, EC50
11,220 to 25,119 nM, 10 distinct compounds after normalisation.** On the
three entry assays it would have been 596, 209 and 96 edges.

**This is pre-existing and was not introduced by the entry assays.** It was
found while adding them because they make it twenty times larger.

## The decision taken — **DECISION, 2026-11-17**

**A fitted hit past the active threshold is refused from BOTH classes,
counted and named** (`weak_hit_above_active_threshold`). Not active: its
potency is past the threshold the project applies everywhere. Not inactive:
the assay fitted a curve, so calling it inactive asserts the opposite of
what was measured. The refusal still records the EC50 it refused on, so the
count can be audited rather than taken.

**CONSEQUENCE FOR AN ALREADY-REGISTERED CELL, stated rather than buried.**
`prereg_structure_scorer.md` fixes the secondary NCATS cell at **89
positives against 7,343 negatives**. After this amendment it is **89 against
7,333**: ten compounds leave the negative class and **no positive moves**.
The paragraph above saying "the registered SARS-CoV-2 CPE cell is unchanged"
was written before this was found and is superseded by this one. The
direction is the conservative one -- it removes known inhibitors from the
class of non-inhibitors -- and it slightly *reduces* that cell's power,
which was already too low to reject.

**SCOPE, FIXED HERE SO IT CANNOT GROW LATER.** This is the NCATS layer only.
The identical question exists for ChEMBL -- a published "= 25 µM" with a
stated relation becomes an inactive label there too -- over a population
orders of magnitude larger. It is **not changed, not measured, and filed**
rather than folded into this amendment, because changing `labels.classify`
moves every cell in the project and that is a separate registration.

## P1, restated

- **P1a. The three entry cells' POSITIVE counts are 336, 97 and 73.**
  Verified against the built layer.
- **P1b. Their negative counts are not predicted here**, because the
  label-level figure depends on structure normalisation and dedup that only
  a rebuild produces. The edge-level counts are 2,805, 1,576 and 2,011 and
  the refused weak hits are 596, 209 and 96.
- **P1c. The CPE cell becomes 89 / 7,333.** Measured on the current release
  by removing the ten compounds rather than predicted.

Everything else above -- P2 through P5, the admission rule, the refusals,
the correction family of four at alpha = 0.0125 -- is unchanged and still
precedes any AUC.

## What is pinned

`tests/test_ncats.py::test_every_emitted_edge_classifies_the_way_the_layer_
said_it_would` asserts the contract over **every row of every declared
assay**, not a fixture -- the disagreement needs a real fitted curve whose
potency lands past the threshold, and which rows do that is a property of
the screen. It reports 17 disagreements against the module committed in
`a9c05f9` and 0 against the module as amended.

---

# Second amendment, 2026-11-17: one of the four cells cannot exist, and the schema had already said why

**STILL BEFORE ANY AUC.** No score has been computed on any of these cells.
Appended, replacing nothing above.

The rebuild converged (two passes, `edges_sha 545fb7d87b0aec08 ->
c0c4daccd48e9979`, second pass byte-identical; 78,649 nodes, 395,646 edges,
schema 0.16.0). It produced three usable plate cells, not four.

## The finding: edge identity does not include `assay_type`

`assemble.edge_key` is `(subject, predicate, object)` plus
`primary_knowledge_source` plus the declared `identity_qualifiers`, and
`assay_type` is **not** one of them. `schema/kg_schema.yaml` says so in
terms, and says what making it one would cost: *"NOT an identity qualifier.
Making it one is exactly L2's decision and would split every merged edge."*

So a CPE measurement and an entry measurement of the same compound against
SARS-CoV-2 are **one fact** to this graph, and they merge. Measured in the
built release:

| virus | `cell_based_antiviral` | `cell_based_pseudotype_entry` |
|---|---:|---:|
| SARS-CoV-2 | 7,427 | **196** |
| SARS-CoV | 0 | 1,619 |
| MERS-CoV | 0 | 2,022 |

`(subject, object)` pairs carrying both assay types: **0**. Not "few" --
none, and none is the signature of a merge rather than of a thin screen.
The ingest offered 2,818 SARS-CoV-2 entry compounds (321 active, 2,497
inactive after normalisation); 196 reached the release.

**Why the other two cells are intact.** The object differs, so the edge key
differs, so there is nothing to collide with -- even though the compounds
are largely the same ones: 1,554 of the 1,619 SARS-CoV entry compounds and
1,936 of the 2,022 MERS ones are also in the SARS-CoV-2 CPE compound set.
The collision is specific to the one virus screened twice.

## The decision — **DECISION, 2026-11-17: the SARS-CoV-2 entry cell is REFUSED**

It survives as 60 positives against 136 negatives. It is refused **not for
its size but for its membership**. Measured: the overlap between the
surviving entry compound set and the CPE compound set is **0 of 196**. 195
of the 196 were never in the pre-entry CPE screen at all, and the one that
was (`INCHIKEY:VYQNWZOUAUKGHI-UHFFFAOYSA-N`) lost its CPE edge to the entry
edge in the merge.

That is the definition of a complement: the cell contains exactly the entry
compounds CPE did not screen. It is neither a random sample nor a complete
one of the entry assay, so an AUC computed on it measures the structural
character of "what the CPE library missed" and not the structural character
of entry inhibition. **It must not be reported as a plate-matched cell**, and
the refusal is pinned in code rather than left to discipline.

The registration above says SARS-CoV-2 "will have two plate cells, CPE and
entry... they are not merged." That paragraph was written without checking
whether the graph could express the distinction. It cannot, short of L2.
**It is superseded by this one.** Nothing was tuned to produce this outcome:
the merge is a property of the identity rule, which predates the
registration by months.

## Correction family — **DECISION, superseding the family of four**

- **Three cells, two-sided.** alpha = 0.05/3 = **0.016667**, **z = 2.3940**.
- The smallest AUC each surviving cell could distinguish from chance:

| cell | positives | negatives | AUC must exceed | by `detectable_auc` |
|---|---:|---:|---:|---:|
| SARS-CoV-2, CPE | 89 | 7,338 | 0.5737 | 0.5762 |
| SARS-CoV, entry | 95 | 1,524 | 0.5731 | 0.5753 |
| MERS-CoV, entry | 71 | 1,951 | 0.5835 | 0.5864 |
| ~~SARS-CoV-2, entry~~ | ~~60~~ | ~~136~~ | ~~0.6074~~ | **refused** |

**Two columns because the two methods disagree, and H1 above used the
looser one.** Column four is the interval half-width evaluated at AUC =
0.5, which is what H1's table was computed with and is therefore the
apples-to-apples comparison against its 0.5769. Column five is the
project's own `diagnose_structure_scorer_design.detectable_auc`, which
re-evaluates the Hanley-McNeil variance at each candidate AUC and is the
figure the evaluation would actually print -- sterner by about 0.0025.
**No statement in this amendment turns on that gap**, and column five is
the one to quote outside this document.

- `MIN_REACHED_POSITIVES` = 10, unchanged; all three clear it. The metric is
  unchanged.

## What this exercise bought, stated against what it was for

**P4 is no longer testable.** It named the SARS-CoV-2 entry cell -- 335
positives, half-width below 0.05, bar 0.5416 -- as "the prediction the whole
exercise exists to test." That cell does not exist. No surviving cell has a
bar below 0.573, which is where the 89-positive CPE cell already sat
(0.5769 at the family of four). **The entry assays did not buy detection
power for SARS-CoV-2.**

What they did buy is narrower and real: **SARS-CoV and MERS-CoV now have
complete plate-matched cells where they had none.** Those are the two
viruses whose ChEMBL cells read 0.5763 (weak) and 0.4926 (below chance).
Independently generated labels for them are new evidence about those two
figures. They are not more power for the one cell this document was written
to strengthen.

P2 falls with P4 (it was arithmetic about that cell). P3 and P5 stand, and
P5 now covers two cells rather than two of four.

## The label route of H2 is clean — measured

H2 asked whether the rebuild moves the ChEMBL cells. **It does not move them
at all.** All **15** ChEMBL cells are byte-identical to the pre-entry
snapshot (`doc2/records/release_operands_2026-11-17_pre_entry_assays.json`),
positive and negative sets both, 5,111 positives and 12,854 negatives.

This closes only the label route. The graph-structure route -- new compound
nodes changing the similarity neighbourhoods that the scorer ranks with --
is untouched by this and is what `make auc-delta` is for. It has not run.

## P1c failed too, and the accounting is worth keeping

The first amendment predicted the CPE cell becomes **89 / 7,333**: 89
positives, ten weak-hit compounds removed from 7,343 negatives. The built
release gives **89 / 7,338**.

The positives are byte-identical (`pos_sha256
f0dfc331cc232729a0100323dd09cf8269b44baf027cda67e440db0140d5da99`, the
snapshot's). The negatives lost 5 and gained 0, traced compound by compound
against the archived pre-entry release:

- **4** of the ten weak-hit compounds left both classes, as predicted.
- **6** remain negatives, correctly: a second row for the same normalised
  structure carries a genuine curve-class-4 measurement, and all six now
  hold a `> 20,000 nM` inactive edge. The weak-hit refusal removed the
  fitted row; the inactive row was always there and still is.
- **1** further compound (`VYQNWZOUAUKGHI`) left because its CPE edge lost
  the merge to an entry edge -- the collision above, costing the CPE cell
  one negative.

−4 − 1 = −5. **The prediction's error was treating a merged edge as a row.**
Ten edges carried a weak hit; six of those ten structures had another row
saying inactive, and refusing the first does not delete the second. That is
the layer behaving correctly and the prediction being written one level too
coarse.

## What is pinned

- `labels.build_labels(assay_types=...)` -- a cell is single-protocol by
  construction or it is not single-protocol. Without it this release's
  SARS-CoV-2 ncats cell is **149 / 7,474**, the pooled figure the "no
  pooling" decision forbids, and nothing in the old code could have told
  the difference.
- `scripts/structure_scorer.py` `POPULATIONS` -- the protocols are separate
  populations, and `secondary_ncats` keeps its name so the cell quoted by
  `prereg_structure_scorer.md` is still the cell being compared.
- The refusal of the 60/136 remnant, so it cannot be reported by running the
  scorer with different arguments.
