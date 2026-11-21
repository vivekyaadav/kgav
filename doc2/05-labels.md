# 05 — Labels, thresholds, and the AUC helpers

`src/kgav/labels.py`. Everything the evaluation measures is measured against
these sets, so this is the one module where a wrong default moves every figure
in the repository.

## Negatives are MEASURED inactive, not untested

Every evaluation before this treated the positives as compounds measured
active and **everything else in the pool** as negative. Most of that pool was
never tested against the virus, so "lift 1.40" meant *ranks actives above
untested compounds* — which a graph can achieve by separating well-studied
chemical matter from obscure chemical matter, with no antiviral insight at all.

ChEMBL's censored values fix it. A row reading `EC50 > 10000 nM` is a compound
that **was** assayed and found inactive. They are the best negatives available
anywhere in this project — better than mining failed trials, because they are
measurements rather than inferences about why a programme stopped.

## The two thresholds, and which of them is derived

| constant | value | status |
|:---|---:|:---|
| `inactive_above_nm` | 10,000 nM | **measured** from the corpus |
| `active_below_nm` | 10,000 nM | a decision about the task, **not** derived |

**The inactive floor is measured, not conventional.** The original
justification was that 10 µM "is the concentration above which a cell-based
antiviral hit is conventionally not pursued" — an appeal to convention with no
citation, and the weakest sentence in the module. It is not needed, because a
censored row states the answer: `EC50 > 20000` is a screen reporting **where it
stopped**. The distribution of those bounds is therefore a direct measurement
of the assay ceilings in this corpus, and `censored_bound_census()` computes
it. Over the 21,656 censored rows of the label predicate on this
release, across 194 distinct bounds:

- **98.5% sit at or above 10,000 nM** —
  the evidence the floor admits;
- 1.5% sit below it — exactly the weak bounds the rule
  exists to refuse;
- the commonest bounds are 20 µM (9,901 rows, 45.7%), 57.1 µM (3,312 rows, 15.3%), 99.5 µM (2,549 rows, 11.8%), 100 µM (2,276 rows, 10.5%).

The corpus would support anywhere in 10–20 µM and 10 is the conservative end
of that range.

**The two are tied, and that tie is a consistency requirement rather than a
convention.** A censored row `> X` proves inactivity only if
`X >= active_below_nm`; below that the compound's true value could lie between
`X` and the active cut, which would make it **active**. Violate it and one
compound can be both:

```
classify({">": 2000}, inactive_above=1000, active_below=10000) -> "inactive"
classify({"=": 5000},  inactive_above=1000, active_below=10000) -> "active"
```

**What has no derivation is the value of `active_below_nm`.** "Above what
bound can an inactive call be believed" is a fact about assay ceilings and is
measured; "what counts as a hit" is a decision about the task, and no
distribution in this corpus answers it. It is also only **half** of the
convention it invokes: published antiviral screening criteria pair the potency
cut with a selectivity requirement (SI ≥ 5–10), and this module applies the
potency clause alone. Of 5,465 positive (compound, virus) pairs, 2,167 (39.7%)
pass SI ≥ 5, 392 of the 2,559 that have an index **fail** it, 65 of those have
SI < 1 — cytotoxic below the effective concentration — and 2,906 (53.2%) have
no cytotoxicity data at all.

Two consequences are left as declared open choices, because they would move
every published figure: there is **no undecided band** (9,999 nM is active and
10,001 nM is inactive), and **56.5% of positives at 10 µM are not positive at
1 µM**, so the majority of the positive set is weak-potency chemical matter
grouped with sub-nanomolar compounds under one label.

## Four ways a measurement decides nothing

`classify()` returns `None` — and the reasons are different questions with
different remedies:

1. **The relation and value do not settle it.** `> 100 nM` is a weak lower
   bound, not an inactive compound.
2. **The edge's potency fields disagree.** `merge_edges` joins the records of
   one pair into a single edge, so a compound with an EC50 of 400 nM and an
   IC50 of 40,000 nM arrives as **one** edge carrying both. `_potency`
   returned the first field present in a module-level tuple, so that edge
   classified **active** — and had the tuple been written the other way round
   it would have classified **inactive**. The verdict came from the
   declaration order of a constant. Measured: 169 of 1,478 multi-field edges
   disagree, 110 of them becoming positives. The protein layer is unaffected,
   because `chembl.py` calls `classify()` on single rows *before* the merge —
   the contradiction is **created** by the merge.
3. **The relation was never stated.** `chembl.py` reads a NULL
   `standard_relation` as `=` and marks the edge `relation_imputed`. For most
   rows that costs nothing; for a value below the floor it decides everything
   — `= 400` is active and `> 400` proves nothing. Of 1,069 imputed rows, 912
   are inactive either way and **157 verdicts rest on the imputation**. Those
   157 are refused.
4. **No potency value at all.**

The rule is derived rather than thresholded: it asks whether the verdict
*changes* when the relation is read as a lower bound, because the two cuts are
independent parameters and a future run with a gap between them needs the
comparison rather than the shortcut.

## Ambiguous compounds are excluded, not assigned

A compound with both an active and an inactive measurement against the same
virus is real — different cell lines, strains, labs — and forcing it to one
side would manufacture a label the data does not support. It is dropped and
counted. `contested_potency` applies the same rule one level down, to a single
edge whose merged fields disagree, which is that disagreement surviving a
merge.

## The dose-series policy

`dose_series_policy` is ``include`, `exclude``, default
**`include`**. A potency ChEMBL flagged as
`outside_dose_series` was not read off the concentrations actually tested; it
is an extrapolation past them. On this release 4,839 label edges carry the
flag — 4,691 inactive, 125 undecidable and **23 active** — so the published
figures rest partly on extrapolated potencies, 23 positives' worth.

The decision was to keep them and make the alternative runnable, for the same
reason the fold-strength question was decided the same way: excluding them
deletes evidence, while recording the flag makes the weakness legible. *"A
consumer can filter"* is only true if something can, so the policy is a named
parameter with a counter on **both** branches. An unknown value is **refused**
rather than treated as the default — a misspelled policy silently reading as
`include` would publish figures under a protocol nobody chose.

## The AUC helpers, and why three of them exist

`evaluate_against_negatives` ranks only measured compounds and returns AUC
with a Hanley–McNeil interval. Three further quantities travel with it, and
they are the difference between a number and a result:

**`auc_ceiling`** — the highest AUC a scorer with this coverage could attain,
ranking everything it reaches perfectly. Unreached compounds take the imputed
floor, so a scorer that reaches few positives cannot separate the classes
however good its ranking is. On the plate-matched cell,
**3 of the 8 channels sit on their own
ceiling** (M1, M7, M8): the thinnest, M7, reaches
0 of 89 positives (0.0%) and
measured 0.4984 against a ceiling of 0.4984. It
achieved the best score it was capable of, and that score is
0.4984.

Channels whose **ceiling** is below chance: M7. A
channel reaching negatives and no positives has those negatives necessarily
outranking the floor-tied positives, so "slightly below chance" is arithmetic
there rather than a result.

**`auc_floor`** — what the same coverage produces with **no ranking skill at
all**, and the more important of the two. The ensemble reads
0.5666 and the degree baseline 0.5700, which reads as
"the graph ties with counting edges". It is not a comparison: degree reaches
100.0% of the positives so its floor is exactly 0.5 and all
of its margin is ranking, while the ensemble reaches 42.7% and
its coverage asymmetry **alone** scores 0.5578 —
**87% of its margin above chance is the imputation**,
not the paths.

**`auc_skill`** — `(auc − floor) / (ceiling − floor)`, the only part of an
imputed AUC that reflects the ranking. Those two cells read
13.3% and 14.0%, which is the defensible version
of "indistinguishable". It returns `None` when nothing was available, rather
than a 0.0 that would rank a channel reaching one compound alongside a genuine
failure.

`MIN_REACHED_POSITIVES = 10` is the floor below which an
AUC describes the imputation rather than the channel. The count alone does not
determine the ceiling — both coverages do — which is why the ceiling is
reported beside it rather than replacing it.

## Imputation is not optional

A metapath scorer only reaches compounds that have paths, and
measured-inactive compounds frequently have none. Evaluating each scorer on
just the compounds it happens to reach lets it **choose its own test set**,
which is how a scorer with 3% coverage posts a perfect AUC. A compound with no
path has no evidence, which *is* a prediction of "not active", so it takes the
floor score.
