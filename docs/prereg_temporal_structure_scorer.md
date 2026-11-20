# Pre-registration: does the structure scorer predict forward, or only separate what is already known?

**Signed 2026-11-19, before any AUC against a post-2021 label exists.**
Every operand below was counted on `data/releases/v0.1-train2021` on that
date and is reproduced by `tests/test_prereg_temporal_structure_scorer.py`.
No ranking of a post-2021 compound by this statistic has been computed.

## Why this exists

`docs/prereg_structure_scorer.md` declines the temporal question in
writing, under *What this cannot do*:

> **It does not test the temporal question.** No held-out-by-year
> evaluation is registered here, because the reference set under a
> temporal split is the pre-2021 actives and that is a different and
> smaller study. It is the obvious follow-up and must be registered
> separately, before it runs.

This is that registration. It matters because **every claim in this
repository is retrospective.** The structure scorer separates known
actives from known inactives measured at the same time, from largely the
same papers. Nothing yet shows it would have called a compound active
*before* anyone measured it, and the one independent axis that could
have corroborated it — clinical labels — came back at 0.4087, below
chance, with degree at 0.7751 (`docs/prereg_clinical_labels.md`).

It also matters because it is the **only** question left where this
project's own measurements say the answer is not already known. The
graph channels cannot answer it: of 2,945 post-2021 test compounds,
**2,803 (95.2%) have no pre-2021 evidence in the graph at all**, so the
metapath channels can rank at most 4.8% of them. That ceiling is
recorded in `data/results/temporal_2021_include.json`'s audits and is a
property of coverage, not of any algorithm. A fingerprint needs only the
molecule, so the structure scorer reaches essentially all of them. **This
study is the one place the two are not comparable, and that asymmetry is
the finding it is built to expose rather than to exploit** — see H1.

## What is held fixed from the cross-sectional study

Changed: the reference set and the evaluated population. **Nothing else.**

- The statistic is `max Tanimoto to known actives`, ECFP4, 2048 bits,
  radius 2, exactly as `scripts/structure_scorer.py` computes it. It is
  imported, not reimplemented.
- Skeleton exclusion (`inchikey_skel`) stays on.
- A compound with no fingerprint is counted, not dropped, and the count
  is reported.
- The AUC is imputed the same way, and `auc_skill`, `auc_floor` and
  reach are reported beside it, because an imputed AUC is not comparable
  across cells with different reach (`docs/evaluation.md`).

## The split, and why this one

`data/releases/v0.1-train2021`, built by `scripts/temporal_split.py
--cutoff 2021 --undated include`, derived from `v0.1`
(`nodes_sha dc6d76a2b3cfb2a1`, `edges_sha 873623c9a9df2049`).

`--undated include` is chosen **in advance** and for a stated reason: it
is the policy the rest of this document set reads, so the split is the
one `docs/evaluation.md` already describes. The two alternatives
(`exclude`, `computed_only`) exist and are **not** run here. Running
three and reporting the best would be three studies reported as one.

`--overlap-policy drop-from-test` is the default and stays. The
`keep` variant refuses on this release by design (12 positive and 38
negative SARS-CoV-2 compounds, 1 negative SARS-CoV, appearing in both
sides), and that refusal is a current record in `data/results`.

## The cells, counted before anything is scored

A cell is evaluable when it has **≥10 pre-2021 reference actives, ≥10
post-2021 test positives and ≥10 post-2021 test negatives** — the same
`MIN_CLASS = 10` the cross-sectional study uses, taken rather than
re-chosen.

    virus          train_act  test_pos  test_neg   cell
    ---------------------------------------------------
    DENV-1               124        71        82   YES
    DENV-2               809       186     2,669   YES
    DENV-3               143        70        79   YES
    DENV-4                86        67        86   YES
    HCoV-229E             43        76        32   YES
    HCoV-OC43             15       147        70   YES
    SARS-CoV-2           359     1,432     1,093   YES
    West Nile             80       181       189   YES
    Yellow fever         105        30        51   YES
    Zika                 181       567     2,495   YES
    ---------------------------------------------------
    HCoV-NL63              0         5         1   no
    JEV                    3        16         4   no
    MERS-CoV             152        48         6   no
    SARS-CoV             195        38         3   no
    TBEV                  61        11         9   no

**THE FAMILY IS TEN.** alpha = 0.05/10 = **0.005000**, two-sided
z = **2.8070**, by the same `NormalDist().inv_cdf(1 - alpha/2)` that
reproduces the registered 2.9137 for the family of 14 and 2.3940 for the
family of 3. The number is written here so a run that finds a different
count of evaluable cells **says so and refuses**, rather than dividing by
whatever it found.

Five cells drop out, and four of the five are coronaviruses or TBEV
losing their negatives, not their positives. MERS-CoV has 152 references
and 48 post-2021 actives and is excluded on **6** negatives. That is a
loss worth naming: the cross-sectional family's two most degree-driven
cells are both among the five, so this study cannot check them.

## The comparator, and the hazard that it may not exist

`degree_ranking` on the **training graph**, loaded with
`skip_predicates=schema.held_out_predicates()`. This is the control that
has decided every claim in this project so far — it beat the structure
statistic in 5 of 17 cross-sectional cells — and it is the reason any of
this is readable.

**H1 — THE DEGREE BASELINE IS EXPECTED TO BE NEAR-DEGENERATE HERE, AND
THAT IS NOT A WIN FOR STRUCTURE.** 95.2% of post-2021 test compounds have
no pre-2021 graph evidence, so most will have training degree 0 and tie
at the floor. A structure AUC that beats a degree baseline made of ties
has beaten nothing. **Registered in advance:** the fraction of each
cell's test compounds with non-zero training degree is reported beside
every AUC, and **a cell where under 20% of test compounds have any
training degree is reported as HAVING NO USABLE COMPARATOR** rather than
as a cell the scorer won. If that applies to most cells, the honest
output of this study is "prospectively measurable, not prospectively
comparable", and that is an acceptable result.

## Controls

**Primary: none beyond the split.** Under a temporal split the candidate
is not in its own reference set by construction, and a near neighbour in
the reference set is *earlier* than the candidate — which is exactly the
prospective situation a medicinal chemist is in. Removing it would
remove the signal rather than a leak.

**Secondary, reported beside it: symmetric k=10 exclusion**, the same
`K_PRIMARY` the cross-sectional symmetric study registered. It answers a
different and real question: how much of any result is the nearest
analogue alone. It does **not** gate the primary.

There is no publication-disjoint arm. 47% of candidates carry no PMID in
the cross-sectional study, the control is asymmetric between classes,
and under a temporal split the date already does most of its work.

## Predictions, recorded before the run

- **P1.** At least one cell's primary AUC interval excludes chance at
  z = 2.8070. *Rationale: the cross-sectional cells run 0.48–0.81 and
  seven clear a sterner alpha.*
- **P2.** The median primary AUC across the ten cells is **lower** than
  the median cross-sectional strength-matched AUC (0.6756 over the same
  viruses). *Rationale: the reference sets are 3–5x smaller and the
  prospective task is harder.*
- **P3.** Fewer than half the cells will have a usable degree comparator
  by the 20% rule in H1.
- **P4.** SARS-CoV-2, the largest cell, will NOT be the strongest. *Its
  1,432 test positives are overwhelmingly post-2021 COVID-era chemistry
  with only 359 pre-2021 references to resemble.*

P1 and P2 can both be true. Recording them separately is the point: a
registration that predicts only success cannot be surprised.

## Hazards

- **H2 — a deposition date is not a discovery date.** ChEMBL records when
  a measurement was published. A 2022 paper may report 2019 work. This
  biases toward leakage, not against it, so a null is safe and a positive
  result is weaker than it looks.
- **H3 — series continuity is a legitimate hit and the easiest one.** A
  post-2021 active that is the next analogue in a pre-2021 series is a
  real prospective success and is also what any chemist would have
  predicted. **The Tanimoto distribution of the true positives is
  reported**, so "it found the next analogue" and "it found a new
  chemotype" cannot be confused. The secondary k=10 arm bounds this.
- **H4 — the undated compounds.** 10,574 undated labels are kept by
  `include`, and `undated_kept_computed` is 298. A compound with no date
  sits in training. If the undated are disproportionately old, this is
  conservative; if disproportionately recent, it leaks. **Not known, and
  not resolvable from this data** — which is why only one policy is run
  and why no amount of this study licenses a claim about the other two.
- **H5 — HCoV-OC43 has 15 references and 147 test positives.** The most
  inverted cell in the family. It is kept because it clears the
  registered floor, and it is flagged here so a strong or weak result
  there is read against 15 reference compounds.
- **H6 — this is one cutoff.** 2021 is taken from the existing split, not
  chosen for this study, but it is still one choice. No sensitivity
  across cutoffs is registered, and a result here does not survive being
  re-run at 2019 without its own registration.

## Refusal conditions

The protocol refuses, and writes a refusal rather than numbers, when:

- the evaluable-cell count is not **10**, for any reason, including a
  rebuild changing a label;
- any cell's reference set is empty after skeleton exclusion;
- the split's `DERIVED_FROM.json` does not name
  `edges_sha 873623c9a9df2049`, i.e. the split is stale against the
  release being read;
- fingerprint coverage of either class in a cell falls below 95%.

## What this cannot do, whatever the numbers

- **It cannot make the graph channels prospective.** Their 4.8% ceiling
  is coverage and is untouched by anything here.
- **It cannot license a clinical claim.** The labels are in-vitro
  concentration-response at one threshold, the same single predicate
  everything else rests on, and the clinical axis already anti-correlates.
- **It cannot establish mechanism.** It names a neighbour.
- **A null here does not retract the cross-sectional cells.** They
  measure separation among known compounds, which is a real and
  different question; it would mean that separation does not transfer
  forward, which is worth as much.

## Decision rule, fixed now

**If P1 holds and the cells that clear also have a usable comparator
they beat**, the structure scorer has a prospective result and the
shortlist in `data/results/repurposing_candidates.json` acquires its
first forward-looking justification.

**If P1 holds but no cell has a usable comparator**, the result is
reported as "measurable, not comparable" and the shortlist stays
retrospectively justified. That is the outcome H1 says to expect, and
reporting it as a success would be the failure this document exists to
prevent.

**If P1 fails**, the scorer does not predict forward on this data, and
`docs/evaluation.md`'s claim section must say so beside the 8 cells it
currently claims.

In all three cases the result is written to
`data/results/temporal_structure_scorer.json` and the verdict is appended
to this document. Nothing above is edited.


---

# Verdict, 2026-11-19: it predicts forward, and there is nothing to compare it against

Run once, on the split and release this document names, by
`make temporal-structure-scorer`. Figures in
`data/results/temporal_structure_scorer.json`, derived by
`tests/test_temporal_structure_scorer_verdict.py`. **Nothing above is
edited.**

    virus          refs    pos    neg      AUC   95% CI (corrected)   clears    k10   degree  reach
    -----------------------------------------------------------------------------------------------
    DENV-1          124     71     82   0.8025   [0.7004, 0.9046]      YES    0.6181     --    1.3%
    DENV-4           86     67     86   0.7896   [0.6830, 0.8961]      YES    0.6868     --    1.3%
    DENV-2          809    186  2,669   0.7518   [0.6926, 0.8109]      YES    0.7076  0.5127   0.3%
    DENV-3          143     70     79   0.7300   [0.6132, 0.8468]      YES    0.5901  0.5020   2.7%
    HCoV-229E        43     76     32   0.7136   [0.5713, 0.8560]      YES    0.7105     --    6.5%
    SARS-CoV-2      359  1,432  1,093   0.6209   [0.5899, 0.6519]      YES    0.6360  0.5071   5.0%
    HCoV-OC43        15    147     70   0.4348   [0.3163, 0.5533]       -     0.4838  0.5031   1.8%
    Zika            181    567  2,495   0.4218   [0.3859, 0.4577]      YES    0.4572  0.5044   0.5%
    Yellow fever    105     30     51   0.3722   [0.1962, 0.5482]       -     0.3863  0.5317   6.2%
    West Nile        80    181    189   0.3532   [0.2735, 0.4329]      YES    0.4521  0.4923   1.9%

`--` is a degree column the evaluator could not produce at all, not a
zero. alpha = 0.005000, z = 2.8070, family of 10 as registered.

## All four predictions hold, and one of them is the result

- **P1 holds. Six cells clear alpha above chance**, and four of them —
  the dengue serotypes — sit between 0.73 and 0.80 on reference sets of
  86 to 809 pre-2021 actives. **The structure scorer does rank post-2021
  actives above post-2021 inactives it has never seen.** That is the
  first forward-looking result in this repository.
- **P2 holds.** Median primary AUC **0.6673** against the
  cross-sectional strength-matched median of **0.6840** over the same
  ten viruses. *The figure written above the run, 0.6756, was wrong; the
  derived value is 0.6840. It was typed rather than computed, and the
  prediction survives either way.*
- **P3 holds, and maximally. ZERO of ten cells have a usable degree
  comparator.** Training-graph reach runs **0.3% to 6.5%**, against a
  registered floor of 20%. The prediction was "fewer than half".
- **P4 holds.** SARS-CoV-2 ranks **6th of 10** at 0.6209, as the largest
  cell with the thinnest reference base.

## The registered decision rule returns MEASURABLE, NOT COMPARABLE

H1 said to expect this and said in advance that reporting it as a
success would be the failure this registration exists to prevent. So:

**The scorer has a prospective signal. It does NOT have a prospective
result, because nothing in this graph can serve as its null.** A degree
baseline that sees 0.3% of DENV-2's test compounds is not a baseline;
its 0.5127 is 2,849 ties and nine ranked compounds. The comparison that
decided every cross-sectional claim — structure against degree — cannot
be made here at all.

**This is the clearest statement of the knowledge graph's contribution
that the repository has produced.** Of 2,945 post-2021 test compounds,
2,803 have no pre-2021 evidence of any kind in the graph. The scorer
reaches 7,556 of 7,556 (100.0%) because a fingerprint needs only the
molecule. The graph is not losing a comparison; it is absent from it.

## Two cells rank significantly BACKWARDS, and that is not noise

**West Nile 0.3532 [0.2735, 0.4329] and Zika 0.4218 [0.3859, 0.4577]**
both clear the corrected alpha *below* chance: the pre-2021 actives'
nearest post-2021 neighbours are systematically the compounds that
turned out INACTIVE. Zika's interval is tight on 567 positives and 2,495
negatives, so this is not a small-sample artifact.

Cross-sectionally both are clean structure cells — West Nile 0.7079 and
Zika 0.7032 strength-matched, each above a non-clearing degree baseline.
**The same two cells invert under a temporal split.** Nothing in this
study says why. The obvious reading — that post-2021 Zika and West Nile
screening deliberately explored away from the known actives, so
resembling them predicts inactivity — is a hypothesis, not a measurement,
and it is not tested here.

## What this changes in the rest of the repository

- `docs/evaluation.md`'s claim section claims 8 cross-sectional cells.
  **Those claims stand**; this measures a different and harder question.
  What it adds is that the separation transfers forward for the dengue
  serotypes, HCoV-229E and SARS-CoV-2, and **inverts** for West Nile and
  Zika.
- `data/results/repurposing_candidates.json` stays **retrospectively**
  justified, as the decision rule requires. Its DENV-1/2/3 lists now
  have a forward-looking argument behind them that its West Nile and
  Zika lists do not — and this study says the Zika list may be
  anti-predictive.
- The 4.8% ceiling is no longer a prediction about the graph channels.
  It is a measured property of every cell in this table.

## What is still not licensed

Everything the registration said, unchanged: no mechanism, no clinical
claim, one cutoff, and a deposition date is not a discovery date (H2).
H3's analogue question is now answerable from the stored scores and is
**not** answered here; the k=10 secondary column drops the dengue cells
by 0.06 to 0.18 and leaves HCoV-229E almost untouched, which is a
pointer, not a finding.

---

## Amendment, 2026-11-21 — DENV-4 clears here and anti-ranks cross-sectionally

**DENV-4 is the one cell where the two protocols disagree, and the
disagreement is not a contradiction.** It is recorded because the pair
reads as one, and because a reader who sees only this registration would
take the prospective pass to license a shortlist. It does not.

Measured by `scripts/diagnose_protocol_disagreement.py`
(`make protocol-disagreement`), which computes no new claim and amends
no verdict.

### The two results

| protocol | AUC | verdict |
|---|---:|---|
| cross-sectional, strength-matched (k_cell 12) | 0.4792 | below chance — the shortlist tool REFUSES this cell |
| prospective, no exclusion | 0.7896 | clears |
| prospective, symmetric k10 | 0.6868 | clears |

### It is not the populations

They are nested exactly. Every prospective test positive is a
cross-sectional positive (67 of
67), and so is every reference
(86 of 86).
The cross-sectional cell's 153 positives are
the 86 pre-2021 actives — which the prospective
split uses as REFERENCES — plus the 67 post-2021
actives it tests. Train and test overlap in
0 compounds.

### It is not era pooling

A pooled AUC can sit below both of its parts, so the published cell was
split by era under its own control:

| scored | AUC |
|---|---:|
| all vs all — *the published cell* | 0.4792 |
| post-2021 vs post-2021 | 0.4842 |
| pre-2021 vs pre-2021 | 0.6105 |
| pre-2021 positives vs post-2021 negatives | 0.3791 |

The post-cutoff half is **also below chance**. Pooling is not the
explanation.

### It is the reference set

Those post-cutoff compounds are exactly the prospective test set, so the
two numbers are directly comparable — identical compounds, identical
labels:

    post-cutoff compounds, references = all 153 actives : 0.4842
    post-cutoff compounds, references = 86 pre-2021 actives : 0.7896

**The mechanism, measured.** For a post-2021 DENV-4 active, of the 12
nearest references the cross-sectional control removes, the share that
are its own post-2021 era-mates is **median 100%, mean
83%** over 67 compounds.

The DENV-4 actives are tight chemical series. Cross-sectionally, a
post-cutoff active's nearest references ARE its own series, and the
strength-matched control removes them by construction — that is what it
is for. Prospectively, the series is not in the pre-cutoff reference set
at all, so the control cannot remove it. **The cross-sectional control is
strictly stronger for exactly these compounds.**

### Both results stand, and they answer different questions

- **cross-sectional, strength-matched** — *setting aside a compound's own
  analogues, does chemistry separate actives from inactives?* For DENV-4,
  no: it anti-separates. Strip the series and an active resembles the
  remaining actives LESS than an inactive does. That is a real statement
  about how clustered this virus's actives are.
- **prospective** — *does pre-cutoff chemistry rank post-cutoff
  compounds?* Yes, 0.7896.

### What this amendment forecloses

**A prospective pass does not license a shortlist.** The shortlist ranks
never-screened compounds against ALL known actives, which is the
cross-sectional setting — the one where DENV-4 is below chance and
`scripts/repurposing_candidates.py` refuses it. That refusal is correct
and stays.

What the prospective pass licenses is narrower and stated here in full:
that someone holding the pre-2021 DENV-4 series would have ranked what
came next, by analogy to that series. It rests on cross-era analogy,
which is how medicinal chemistry proceeds and is precisely what the
strength-matched control exists to remove.

**This also sharpens H3, which the verdict above left open.** The
analogue question is not a pointer any more for this cell: the
prospective signal here IS analogy to a known series, and the measurement
is the era-mate share above.

---

## Amendment, 2026-11-21 — H3 answered: three of the six are analogy

The verdict above left H3 open in as many words: *"H3's analogue
question is now answerable from the stored scores and is **not**
answered here; the k=10 secondary column drops the dengue cells by 0.06
to 0.18 and leaves HCoV-229E almost untouched, which is a pointer, not a
finding."*

It is answered now, by `scripts/diagnose_prospective_analogy.py`
(`make prospective-analogy`). The pointer was right and understated.

### The ladder

The registered prospective AUC, recomputed with the **k nearest
pre-cutoff actives withheld from every candidate alike**. `k=0` is the
registration's primary and `k=10` its secondary, so two rungs reproduce
published figures; the rest extend them. `\*` = still clears the
registered corrected alpha.

| virus | k=0 | k=5 | k=10 | k=20 | k=40 | clears up to |
|---|---:|---:|---:|---:|---:|---|
| DENV-2 | 0.752\* | 0.719\* | 0.708\* | 0.687\* | 0.655\* | k≤40 |
| HCoV-229E | 0.714\* | 0.737\* | 0.711\* | 0.678\* | 0.662\* | k≤40 |
| SARS-CoV-2 | 0.621\* | 0.635\* | 0.636\* | 0.623\* | 0.614\* | k≤40 |
| DENV-4 | 0.790\* | 0.727\* | 0.687\* | 0.516 | 0.436 | k≤10 |
| DENV-1 | 0.802\* | 0.690\* | 0.618 | 0.462 | 0.357 | k≤5 |
| DENV-3 | 0.730\* | 0.654\* | 0.590 | 0.452 | 0.368 | k≤5 |

**3 of 6 survive to k = 40**: DENV-2, HCoV-229E, SARS-CoV-2.
Forty of their closest known actives removed, and the order holds.

**2 stop clearing before k = 10**: DENV-1, DENV-3. Their
signal was carried by a handful of near neighbours.
**DENV-4** clears at k = 10 and not at k = 20.

### The headline AUC is close to anti-correlated with robustness

    ranked by the published AUC (k=0) : DENV-1 > DENV-4 > DENV-2 > DENV-3 > HCoV-229E > SARS-CoV-2
    ranked by robustness     (k=40)   : HCoV-229E > DENV-2 > SARS-CoV-2 > DENV-4 > DENV-3 > DENV-1

DENV-1 has the **highest** primary AUC of the six and loses the most
(0.802 to
0.357). SARS-CoV-2 has the
**lowest** and barely moves
(0.621 to
0.614, a drop of
-0.007).

**A cell's headline number does not say which kind of result it is.**
That is the finding, and it applies to how every figure in this
registration should be read.

### How close the nearest known active actually is

Median max-Tanimoto to the pre-cutoff reference set, no exclusion:

| virus | positives | negatives | gap |
|---|---:|---:|---:|
| DENV-2 | 0.369 | 0.271 | +0.098 |
| DENV-4 | 0.290 | 0.209 | +0.081 |
| DENV-1 | 0.271 | 0.208 | +0.063 |
| SARS-CoV-2 | 0.263 | 0.229 | +0.034 |
| DENV-3 | 0.251 | 0.210 | +0.042 |
| HCoV-229E | 0.208 | 0.170 | +0.039 |

All six sit well below identity, so this is not near-twin matching in
the median. Combined with the ladder, the reading is that a **minority**
of positives have a very close pre-cutoff neighbour and those few carry
the AUC in the cells that collapse.

### What this changes

- **Nothing is retracted.** Every rung is the registered statistic on
  the registered split at a different exclusion strength.
- **The six clearing cells are not one result.** Three are a chemotype
  signal; three are retrieval of a series. A summary that reports "6 of
  10 clear" without this is true and misleading.
- **It strengthens the DENV-2 choice.** DENV-2 is the only cell that is
  simultaneously the most enriched at the top of its ranking (10.74x),
  powerable at a feasible screen size, and robust to k = 40. The screen
  registered in `docs/prereg_denv2_prospective_screen.md` was chosen on
  the first two; the third was not known when it was signed.
- **H3 is closed.**
