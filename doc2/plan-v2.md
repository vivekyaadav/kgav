# The v2 plan

**Not generated.** This is what v2 is for: ten fixes to the inference in
v1's four measurements, each with the measurement that shows it is needed,
hand-maintained because what it records is a decision and a priority, which
no artifact can derive. Every figure in it is held to the artifacts by
`tests/test_v2_plan_figures.py`, which re-derives all of them and fails if
one has moved.

Written 2026-11-22, against release `v0.1` and the 36 artifacts in
`data/results/` as published in v1.

---

## What v2 is

v1 published four controlled measurements and the apparatus that keeps them
honest. **v2 does not add a measurement. It repairs the inference drawn from
the four already made**, because a review of the published paper against its
own artifacts found that the strongest claim in it rests on a comparison
that was never tested, and the strongest-sounding negative rests on a
population with no power to detect the effect it reports as absent.

**None of the ten items needs new data.** Every one runs against release
`v0.1` and the artifacts already published. Nine of them are hours to a day
of work; one, C9, is three to four days, because it is the only item that
builds a scorer that does not exist yet.

The ordering is by how much a published claim moves if the item is wrong,
not by cost.

---

## A. The inference is the wrong shape

### A1 — There is no paired test anywhere in this repository

**WHAT IS WRONG.** Every scorer-against-baseline comparison in the paper is
two *marginal* confidence intervals, compared by whether each separately
excludes 0.5. That is the "a difference in significance is not a significant
difference" fallacy. It is load-bearing for §3.1, for the whole of §3.3, and
therefore for which viruses get a shortlist.

**THE MEASUREMENT.** `grep -rn "delong\|DeLong\|paired bootstrap" src/
scripts/` returns nothing. On the assay-matched population the paper reports
COMBINED at 0.5666 and degree at 0.5700 and calls them
"three thousandths apart". Nothing in the repository has ever computed an
interval on that difference, and the two scorers are evaluated on **the same
89 actives and the same 7,338 inactives**, so a paired
statistic has far more power than the comparison actually used.

**WHAT V2 DOES.** `scripts/paired_auc.py`: a stratified bootstrap over
compounds, resampling positives and negatives independently, reporting
`AUC_a − AUC_b` with a percentile interval, for every (scorer, degree) pair
in every published cell. DeLong as a closed-form cross-check on the cells
where ties permit it. The statistic goes in `src/kgav/` beside
`auc_interval` so no caller can reimplement it, and
`evaluate_against_negatives` gains a paired sibling rather than being
changed.

**COST.** A day. The bootstrap is embarrassingly parallel and the label sets
are already in memory.

**WHAT WOULD CHANGE.** Every "clears above a degree baseline that does not
clear" sentence becomes a signed difference with an interval. Some of the
eight cells in §3.3 will not survive it — see B5, where one of them clears
chance by 0.0034.

**THIS IS PRE-REGISTERED BEFORE THE RUN.** `prereg_paired_auc.md` fixes the
statistic, the resampling unit, the number of resamples, the correction
family and the decision rule before any difference is computed. The reason
is item B6: the published cells already disagree across controls by up to
0.6550 AUC, so a protocol chosen after seeing a difference could
pick any answer it liked.

### A2 — The negative result has no power statement, and the abstract outruns it

**WHAT IS WRONG.** The paper's §3.1 concludes that 0 of 8 channels clears
chance and the abstract reads this as the graph's contribution being "zero
or negative". On that population the study could not have detected a real
effect of the size three of the channels report.

**THE MEASUREMENT.** Hanley–McNeil at AUC 0.5 on
89 positives and 7,338 negatives gives
SE 0.0308, so the minimum detectable AUC at
uncorrected 95% is **0.5603**.

**All eight channels read below that floor.** The highest of them is
0.5400; the lowest is 0.4984. Not one of them
could have cleared chance on this population at any ranking quality short
of a much larger effect, so **"0 of 8 channels clears" was fixed by the
size of the population before a single path was counted.** It is a true
sentence that carries no information about the channels.

Under the correction family the same paper declares for the metapath
comparison the bar rises to 0.5897, which neither COMBINED
(0.5666) nor the degree baseline (0.5700)
would clear either.

**WHAT V2 DOES.** Every AUC cell gains an `auc_mde` field — the minimum
AUC its own n could have detected, at both the uncorrected and the corrected
alpha — computed in `src/kgav/labels.py` beside `auc_floor`, `auc_ceiling`
and `auc_skill`, so it is impossible to publish a cell without it. The
pathway landscape already does this (`landscape.py` carries the minimum
detectable effect per row); the AUC protocols do not, and the asymmetry was
never a decision. A test fails on any results file holding an AUC with no
`auc_mde` beside it.

The paper's §3.1 verdict is then rewritten from "is zero or negative" to
"is not detectably positive at n=89, where the detection
floor is 0.5603", and the channels inside the blind spot
are reported as *unresolved* rather than as *at chance* — the same
distinction the project already draws with *not evaluable*.

**COST.** Half a day, most of it rewriting prose and re-deriving.

**WHAT WOULD CHANGE.** The headline negative becomes narrower and true. What
survives untouched is the §3.1 result that does not depend on power at all:
the floor/ceiling/skill decomposition, which shows 87% of the
ensemble's margin is coverage asymmetry rather than ranking.

### A3 — One paper, two significance standards, unstated

**WHAT IS WRONG.** §3.1 reports uncorrected 95% intervals and says COMBINED
and degree "clear". §3.3 reports Bonferroni-corrected intervals at
z 2.9137. The word *clears* means two different things in two
sections and the paper does not say so.

**THE MEASUREMENT.** `hard_negatives.json` and its four `_row*` companions
carry no `clears_corrected_alpha` field on any cell;
`structure_scorer_strength_matched.json` carries one on every cell.

**WHAT V2 DOES.** The metapath runs gain the corrected interval, from the
correction family already registered for them, and §3.1's verdict column is
read off it. Where the two standards give different answers both are
printed, because the gap between them is the same kind of fact as the gap
between a naive and a degree-matched p-value, which §3.2 already treats as
its main finding.

**COST.** Hours. `corrected()` already exists and is already imported by the
structure scorer.

**WHAT WOULD CHANGE.** §3.1's "two scorers do clear" becomes "neither clears
the corrected bar". This moves the paper's conclusion further in the
direction it already argues, which is why it is a correction and not a
result.

---

## B. The positive result is narrower than it reads

### B4 — The licensed cells are not independent, and nothing measures it

**WHAT IS WRONG.** §3.3 reports "8 of 17 cells" and the discussion licenses
8 of 16 viruses. Both read as eight independent confirmations. The label
sets overlap heavily and no artifact has ever measured it.

**THE MEASUREMENT**, over the eight licensed cells' positive sets, derived
from `build_labels` on release `v0.1`:

| pair | Jaccard | share of the smaller set |
|:---|---:|---:|
| DENV-1 / DENV-3 | 0.693 | **85.6%** |
| DENV-2 / DENV-3 | 0.163 | 79.3% |
| DENV-1 / DENV-2 | 0.144 | 76.9% |
| JEV / Yellow fever | 0.085 | **63.2%** |
| HCoV-229E / SARS-CoV-2 | 0.023 | 36.1% |

`{DENV-1, DENV-2, DENV-3}` is close to one experiment measured three times;
12 of JEV's 19 positives are also Yellow
fever positives. The effective number of independent clearing cells is
nearer four or five than eight.

**WHAT V2 DOES.** `scripts/cell_dependence.py` writes
`data/results/cell_dependence.json`: the full pairwise positive-set and
negative-set overlap, a clustering at a registered Jaccard cut, and an
effective-cell count. Every "N of M cells" sentence in the documents is then
derived with the effective count beside the raw one.
`data/results/label_overlap.json` measures evidence-against-label
publication overlap and is a different quantity; this is the one nothing
covers.

The Bonferroni correction is **not** changed. Over dependent tests it is
conservative on alpha, so the existing verdicts stand; what changes is the
count a reader is given.

**COST.** Half a day. `build_labels` runs in seconds.

**WHAT WOULD CHANGE.** "8 of 16 viruses" acquires a second number. The
shortlist licence itself is unaffected — it is granted per virus, and a
virus whose actives resemble its neighbour's is still a virus somebody wants
a shortlist for.

### B5 — Two licensed cells rest on almost nothing, and the paper does not say which

**WHAT IS WRONG.** §3.3's table marks eight cells **yes** with no indication
that their margins differ by two orders of magnitude.

**THE MEASUREMENT**, corrected lower bound minus 0.5, from
`structure_scorer_strength_matched.json`:

| virus | n+ | n− | AUC | corrected lower bound | margin |
|:---|---:|---:|---:|---:|---:|
| DENV-1 | 195 | 123 | 0.5967 | 0.5034 | **+0.0034** |
| JEV | 19 | 27 | 0.7885 | 0.5810 | +0.0810 |
| TBEV | 72 | 65 | 0.7389 | 0.6164 | +0.1164 |

DENV-1's licence rests on **three thousandths of AUC**. JEV's rests on
19 positives, of which k_cell 12 are withheld
from every candidate as references — leaving
7. Its minimum detectable AUC at the corrected
alpha is 0.7546, so JEV could only ever have produced a very
large number or none at all.

**WHAT V2 DOES.** The margin column is published in §3.3 beside the verdict,
and a registered *minimum evidential n* gate — the same shape as the
existing `MIN_REACHED_POSITIVES = 10` — is applied to the
shortlist licence rather than to the AUC. A cell may clear and still be too
thin to license an act; those are two decisions and v1 ran them together.

**COST.** Hours for the column, a registration and a day for the gate.

**WHAT WOULD CHANGE.** DENV-1 and JEV are the two most likely to lose the
licence. That is the point: the read-out currently offers compounds for JEV
on the strength of seven admissible reference compounds.

### B6 — The paper drops the control that disagrees most

**WHAT IS WRONG.** `doc2/08-results.md` publishes all three controls. The
paper publishes one — strength-matched — and §3.3 refers to "the two
registered near-neighbour controls", silently excluding the first
registration's. A reader of `paper.md` alone cannot learn that JEV's first
registered reading was a strong anti-ranking.

**THE MEASUREMENT.** The licensed cells under all three controls:

| virus | publication-disjoint | symmetric k=10 | strength-matched | spread |
|:---|---:|---:|---:|---:|
| JEV | **0.1452** | 0.8002 | 0.7885 | **0.6550** |
| Yellow fever | 0.5545 | 0.8388 | 0.8117 | 0.2844 |
| DENV-3 | **0.3861** | 0.6496 | 0.6653 | 0.2793 |
| DENV-2 | 0.5742 | 0.7888 | 0.6921 | 0.2146 |
| TBEV | 0.5764 | 0.7464 | 0.7389 | 0.1700 |
| DENV-1 | 0.4489 | 0.6066 | 0.5967 | 0.1577 |
| HCoV-229E | 0.5917 | 0.7100 | 0.6758 | 0.1183 |
| SARS-CoV-2 | 0.7718 | 0.8304 | 0.7803 | 0.0585 |

**THIS IS NOT CONTROL-SHOPPING AND V2 DOES NOT TREAT IT AS SUCH.** The
replacement is sound and was signed before the run:
`prereg_structure_scorer_symmetric.md` measured the publication control
stripping 12.7x more references from positives than from negatives
on DENV-2, tested the asymmetry hypothesis directly, and named its two
counterexamples. The defect is one of **reporting**, not of protocol choice.

**WHAT V2 DOES.** The three-control table goes in the paper, with the
asymmetry argument in two sentences and a pointer to the registration. By
this project's own standard — that a disagreement between protocols is a
result, which is why `make protocol-disagreement` exists — a
0.6550 swing on a licensed cell is a finding and belongs in
§3.3 rather than only in `08-results.md`.

**COST.** Hours. Every figure already exists in three artifacts.

### B7 — The exclusion ladder runs only on the temporal cells

**WHAT IS WRONG.** `diagnose_prospective_analogy.py` recomputes each
temporal cell's AUC at k = 0, 5, 10, 20, 40 and it is the best diagnostic in
the paper: it showed DENV-1 falling 0.8025 to 0.3574 and the
robustness ordering nearly inverting the headline ordering, which is how the
paper can say a cell's headline AUC does not tell you what kind of result it
is. **The cross-sectional cells of §3.3 get exactly one exclusion strength
each** — their own `k_cell` — and the cells most at risk from series
retrieval are precisely the ones being licensed.

**WHAT V2 DOES.** The same ladder on the cross-sectional cells, as
`scripts/diagnose_cross_sectional_analogy.py`, reusing the k-exclusion code
in `structure_scorer.py` rather than reimplementing it, with the registered
`k_cell` and k = 10 as two rungs that must reproduce published figures — the
same self-check the temporal ladder uses.

**COST.** A day. The scoring is already written; this varies one parameter.

**WHAT WOULD CHANGE.** On the temporal cells 2 of 6 stopped
clearing before k=10. If the cross-sectional cells behave the same way, two
or three of the eight licensed viruses are retrieving a known chemical
series rather than generalising, and should be described that way.

---

## C. Missing comparators and untested mechanisms

### C8 — 60.2% of positives are threshold-contingent and no AUC is recomputed

**WHAT IS WRONG.** `09-limits.md` states the label set's sensitivity to the
cut and stops there. `active_below_nm` is never varied from its default in
any script in this repository, so every published AUC rests on one
threshold and no reader can tell which results are an artifact of it.
10 µM in a cell-based CPE assay is also where cytotoxicity begins doing the
work.

**THE MEASUREMENT.** Re-running `build_labels` at a 1 µM cut:
5,276 compound–virus pairs are positive at 10 µM and
2,102 of them are still positive at 1 µM —
3,174 flip, which is **60.2%**. By unique
compound the figure is 59.3%.

**THE DOCUMENT'S OWN FIGURE FOR THIS WAS STALE, AND THAT HALF IS ALREADY
FIXED.** `09-limits.md` read *"56.5% of positives at 10 µM are not positive
at 1 µM"*, and neither of the two readings above is that. The
56.5% was typed in a docstring in `src/kgav/labels.py` and
quoted into the template as literal prose rather than as a placeholder the
generator fills — so `write_doc2.py` had nothing to compare it against and
`--check` passed on it for as long as it was there. **The one figure in doc2
about how much rests on the label threshold was itself the defect doc2 was
built to prevent.**

It is now derived at render time, the docstring's number is deleted rather
than corrected, and the instance is in `14-defect-ledger.md` under *a figure
with no path back to the data*.

**AND THE SWEEP THAT SHOULD HAVE CAUGHT IT ONLY LOOKS AT FOUR-DECIMAL
FIGURES.** `test_doc2_does_not_claim_a_figure_it_did_not_derive` matches
`0.dddd`, because those are AUCs, so **every integer and percentage typed
into a template is unguarded.** The same sentence of `09-limits.md` held
four more of them. Two — 4,839 dose-series edges and
23 of them positives — were typed copies of
the `ds_flagged` and `ds_flagged_positive` placeholders, which the generator already
computes and which `08-results.md` already prints from the derived values:
one quantity, two surfaces, agreeing by luck. A third,
169 contested-potency edges, is in
`build_labels().stats` and had no placeholder. All three are now derived.
The fourth, *157 verdicts rest on a relation nobody stated*, needs a
per-edge sweep and is still typed; it is the only typed figure left on that
line and it is named in the document as such.

**WHAT V2 DOES.** `make label-threshold-sensitivity`: §3.3 and §3.4
recomputed at 1 µM and at 10 µM, published side by side. A cell whose
verdict changes with the cut is reported as threshold-contingent. And the
four-decimal sweep is widened to thousands-separated integers and
percentages, which is what would have caught all five of the figures above
without anybody reading the line.

**COST.** A day. `build_labels` takes seconds per threshold and the
scorers take minutes.

**WHAT WOULD CHANGE.** Possibly a great deal, and this is the cheapest
unknown in the project. More than half the positive class is at stake and
nobody has looked.

### C9 — No embedding or GNN comparator exists

**WHAT IS WRONG.** The paper tests eight degree-weighted path-count
channels and the title generalises to what a knowledge graph contributes.
The field's actual practice is knowledge-graph embedding and GNN link
prediction. A reviewer will say so in the first paragraph, and the reply
cannot be that metapaths were chosen for interpretability — true as that is,
it is not a measurement.

**WHAT V2 DOES.** One comparator, not a survey: DistMult on the training
graph with the held-out predicates stripped, scored on the same cells,
through the **same** floor/ceiling/skill decomposition and the same degree
baseline. The decomposition is the point. The argument v2 should be able to
make is not "embeddings are worse" but *"an embedding's margin above chance
decomposes into coverage and ranking the same way, and the ranking fraction
is the same size"* — which would make §3.1's result about the graph rather
than about one scorer family, and is the honest version of the claim the
title currently makes.

The title changes either way, to name metapath channels rather than
knowledge graphs, because the discussion already scopes the claim correctly
and the title does not.

**COST.** Three to four days, the only item here that is not trivial.

**WHAT WOULD CHANGE.** The scope of the paper's central negative. If
DistMult clears where the metapaths do not, §3.1 is a statement about DWPC
and the paper must say so.

### C10 — The two anti-rankings have a mechanism nobody has tested

**WHAT IS WRONG.** West Nile reads 0.3532 and Zika
0.4218 prospectively, clearing **below** chance, and the paper
says only that "something systematic changed between the training period and
the test period". That is true and it is not an explanation.

**THE HYPOTHESIS.** Post-cutoff negatives are disproportionately *analogues
of pre-cutoff actives*: somebody screened the follow-up series and it
failed. Max-Tanimoto-to-a-known-active then ranks negatives above positives
by construction, and the anti-ranking is a property of how medicinal
chemistry proceeds rather than of the virus.

**WHAT V2 DOES.** One distribution comparison per anti-ranked cell: the
max-Tanimoto-to-pre-cutoff-active distribution for post-cutoff actives
against post-cutoff inactives, with the gap and its interval. Registered
before it is run, because the direction is predicted.

**COST.** Half a day. Every fingerprint and both label sets already exist.

**WHAT WOULD CHANGE.** If it holds, this is the most useful finding in the
project and it is currently filed as a limitation: **following up a known
active by chemical similarity is measurably worse than chance on these two
viruses.** That is an instruction a medicinal chemist can act on, and it
also explains §3.3's vanishing top-of-list precision by the same mechanism.
If it does not hold, the anti-rankings are unexplained and should be said to
be.

---

## The ten items, by what they cost and what they can overturn

| # | item | cost | new data | can overturn a published claim |
|:---|:---|:---|:---:|:---:|
| A1 | paired AUC test | 1 day | no | **yes** — all of §3.3 |
| A2 | minimum detectable effect on every AUC | ½ day | no | **yes** — §3.1's verdict |
| A3 | one significance standard | hours | no | **yes** — §3.1's two clearing scorers |
| B4 | cell dependence measured | ½ day | no | no — changes a count |
| B5 | margins published, licence gated on n | 1 day | no | **yes** — DENV-1, JEV |
| B6 | three-control table in the paper | hours | no | no — changes what is shown |
| B7 | exclusion ladder on the cross-sectional cells | 1 day | no | **yes** — 2–3 licences |
| C8 | 1 µM threshold sensitivity | 1 day | no | **yes** — unknown, largest |
| C9 | one embedding comparator | 3–4 days | no | **yes** — the paper's scope |
| C10 | the anti-ranking mechanism | ½ day | no | no — adds a finding |

Eight to eleven working days, and **no item requires a new screen, a new
ingest or a rebuild.** Every one of them runs against release `v0.1` and the
artifacts already published, which is what makes them worth doing before
anything is added to the graph: v1's constraint was data, and v2's
constraint is that the data already collected has not been fully read.

---

## What v2 is not

- **Not a new release of the graph.** No ingest changes, no new source, no
  rebuild. If `v0.1` moves, every figure above is invalidated by
  `make check-results` and the items are re-derived, not re-argued.
- **Not a retraction.** Three of v1's four measurements stand as published.
  What moves is §3.1's wording, §3.3's width, and the shortlist licence for
  at most three viruses.
- **Not tuning.** No item adjusts a scorer, a damping term or a threshold to
  improve a number. C8 varies the label threshold to *measure* the
  dependence, and publishes both cuts; a cut will not be chosen because its
  AUC is better.
- **Not more channels.** The metapath set is frozen. §3.1's problem is
  power and inference, and adding a ninth channel would enlarge the
  correction family without adding reach.

## Where the new work is registered

Each item that computes a comparison gets a registration signed before its
code exists, in `doc2/registrations/`, on the existing convention —
amended by appending, never rewritten, with a test asserting the code obeys
its alpha, its family size and its exclusion rule:

| item | registration |
|:---|:---|
| A1 | `prereg_paired_auc.md` |
| B5 | `prereg_minimum_evidential_n.md` |
| B7 | `prereg_cross_sectional_analogy.md` |
| C8 | `prereg_label_threshold_sensitivity.md` |
| C9 | `prereg_embedding_comparator.md` |
| C10 | `prereg_analogue_followup.md` |

A2, A3, B4 and B6 compute no new comparison — they report a quantity beside
an existing one, or print an artifact the repository already holds — so they
need a test, not a registration.
