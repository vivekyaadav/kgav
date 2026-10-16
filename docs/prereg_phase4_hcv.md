# Pre-registration: the HCV retrospective, as the evidence permits it

Written **before** the code, and committed before the measurement, because a
prediction made afterwards is a description. If the result contradicts what
is written here, this file stays as it is and the verdict goes beneath it.

Date: 2026-10-16. Release: v0.1 (74,935 nodes / 374,795 edges, schema
v0.13.0). Source for every figure below: ChEMBL 37, read by
`scripts/phase4_survey_chemistry.py`, which writes nothing and recommends
nothing.

**The outcome has not been computed.** Every number in this file describes
how much evidence exists and how it is distributed. Not one of them is a
result. That line is the only thing that makes the prediction below a
prediction.

## What phase 4 was supposed to test

The folds layer ships a catalytic-type gate whose whole job is one refusal.
Coronavirus 3CLpro and HCV NS3 share the trypsin-like protease fold — phase
0 measured it, and the shipped layer records 54 refused-despite-shared-fold
pairs — and differ in the catalytic residue, cysteine against serine.
`schema/kg_schema.yaml` states the consequence of removing the gate: M8
"traverses from SARS-CoV-2 Mpro to HCV NS3 — same superfamily, opposite
nucleophile — and recommends nirmatrelvir for hepatitis C".

`phase0_verify_folds.py` and `ingest_folds.py` both say phases 5–10 must not
be built if the conservation premise is refuted. Both measure ANNOTATION —
which InterPro entries two chains share. Annotation cannot say whether a
compound that inhibits one protease inhibits the other. Phase 4 is where the
chemistry was to be asked, as a 2×2: hold the fold constant and vary the
nucleophile (A against B), then hold the nucleophile constant and vary the
fold (A against C1).

## What the evidence permits, measured

Compounds with a protein-level, exact-relation, unit-convertible potency at
BOTH ends of each pair. `units` counts INDEPENDENT observations: one per
contributing publication, plus one for each pair whose two ends come from
different papers.

    shared fold  same nucleophile  pair                      compounds  units
    yes          yes               A   3CLpro + picorna 3C          12      4
    yes          yes               A2  HCV NS3 + flavi NS3           0      0
    yes          no                B   3CLpro + HCV NS3              4      4
    yes          no                B'  3CLpro + flavi NS3            3      3
    no           yes               C1  3CLpro + CoV PLpro           49     13
    no           no                C2  PLpro + HCV NS3               0      0

**The paired design is refused.** No arm reaches ten independent
observations, the floor `calibration.MIN_REACHED_POSITIVES` already sets for
a class worth evaluating. Three facts, in order of how much they cost:

1. **Arm A is one paper.** Twelve compounds, ten of them with both ends in a
   single ChEMBL document, nine from one publication. Four independent
   observations. That is a medicinal-chemistry series, almost certainly
   designed as a dual inhibitor, and a correlation over it would measure the
   series.
2. **Arm B is four.** The gate's central refusal — the pair the schema names
   — has four paired compounds across ends holding 3,040 and 2,125.
3. **The only well-populated cell is the one where the fold is absent.** C1,
   two cysteine proteases of the same virus sharing no fold, carries 49
   compounds over 13 independent observations: four times arm A's pairs.
   Overlap tracks who screened what together, not what the enzymes have in
   common. The COVID-era papers that ran 3CLpro ran PLpro beside it; HCV
   protease chemistry and flavivirus NS3 chemistry never met at all, which
   is why A2 and C2 are zero across thousands of compounds.

## What is therefore NOT tested, and will not be

**The catalytic-type gate is untestable on ChEMBL 37.** Arm B has four
independent observations and arm B′ has three. M8's safety rule — the one
thing standing between the graph and recommending nirmatrelvir for hepatitis
C — has never been checked against chemistry and cannot be with this
database. No result in this file bears on it.

That narrows `ingest_folds`'s own gate on phases 5–10. The conservation
premise has two halves, fold and nucleophile, and only one of them was ever
testable here. Anywhere the project relies on the catalytic-type gate it
should say that the gate is an untested safety measure, not a validated one.

## The test that remains

The broad set: compounds with a protein-level HCV NS3 potency AND at least
one coronavirus measurement of any kind in ChEMBL 37 — enzyme or antiviral,
any `target_type`. This is the retrospective as the literature actually ran
it in 2020: the HCV protease inhibitors that got tested against SARS-CoV-2.

    identity                    compounds   units
    exact InChIKey (PRIMARY)           12      12
    InChIKey skeleton (sensitivity)    17      17

Compounds equal units, which is the one piece of luck in this whole survey:
the HCV enzyme literature (2005–2015, pharma) and the coronavirus literature
(2020 on, repurposing) are disjoint, so **every pairing is a cross-paper
observation.** That is exactly the independence arm A lacks.

It answers a WEAKER question than arm B. Arm B asks whether an HCV NS3
inhibitor inhibits the coronavirus protease. The broad set asks whether it
stops the coronavirus — by any mechanism, including ones that have nothing
to do with Mpro. A reader will take "the HCV retrospective" to mean the
enzyme test unless told otherwise, and this paragraph is the telling.

### Definitions, fixed here

- **Compound identity**: `INCHIKEY:<standard_inchi_key>`, the node id
  `ingest_activities` writes. Rows without one are dropped, as the ingest
  drops them. The 14-character skeleton — `normalize_chem.same_as_rule`'s
  declared `shared_inchikey_skeleton:salt_or_stereo_or_tautomer` merge — is
  the pre-declared sensitivity, not the primary.
- **HCV NS3 potency**: median of log10 nM over protein-level rows
  (`chembl.assay_readout` == "protein"), within one assay type chosen
  IC50 > Ki > Kd > EC50, relation `=` only.
- **Coronavirus activity**: `labels.classify` on each coronavirus row, with
  its declared defaults — active at or below 10,000 nM, inactive above. The
  implementation CALLS that function; it does not restate the rule.
  A compound is active if classify returns "active" for at least one
  coronavirus row.
- **Background**: every other compound with a coronavirus measurement in the
  same extract, classified the same way. Numerator and denominator come from
  one query, so no second extraction can disagree with the first.

## The prediction

**P1.** HCV NS3 protease inhibitors are NOT enriched for coronavirus
activity relative to background. One-sided Fisher exact test, α = 0.05,
alternative = enrichment.

This is the gate's prediction, and it is a null. That asymmetry is the whole
character of this test and must not be lost:

> **This test can refute the gate. It cannot confirm it.**

A significant enrichment says the gate is deleting real signal. A
non-significant result says nothing was detected at n = 12 — not that
nothing is there.

**P2.** A significant DEPLETION is not evidence for the gate either, and
will not be reported as such. HCV NS3 inhibitors are a narrow chemotype of
large macrocyclic peptidomimetics; being less active than a repurposing
library is what that chemotype predicts regardless of any fold or
nucleophile.

## The arithmetic, done before the data

Minimum actives out of n for a one-sided binomial test to reach p < 0.05,
against a background rate the measurement has not yet produced:

    background   n=12: k   as rate   n=17: k   as rate
            5%         3       25%         4       24%
           10%         4       33%         5       29%
           15%         5       42%         6       35%
           20%         6       50%         7       41%
           25%         7       58%         8       47%
           30%         7       58%         9       53%

**The test detects roughly a three-fold enrichment and nothing finer.** At a
15% background it needs 5 of 12. A true enrichment to 25% would be missed
more often than not. Written here so that a null cannot later be read as a
tight bound.

## Decision rule

- **Enrichment significant (P1 rejected)** → the conservation premise's
  nucleophile rule is refuted at the level this evidence reaches. The
  catalytic-type gate is removing signal, and phases 5–10 must not be built
  on it as it stands. This is the outcome `ingest_folds` names.
- **Not significant** → nothing is confirmed. Phase 4 returns NO RESULT, the
  gate stays an untested safety measure, and that wording goes wherever the
  project describes it.
- **Depletion significant** → reported as a fact about the chemotype, per
  P2, and claimed as evidence for nothing.

In every branch the paired arms stay refused and the file above stands.

## What would make the real test possible

Not pursued here; recorded so the next person does not rediscover the dead
end.

1. **A source with cross-enzyme panels.** The missing data is one compound
   measured on two viral proteases in one protocol. ChEMBL 37 has four such
   compounds for the pair that matters.
2. **A set-level design** — compare the chemical space of 3,040 3CLpro
   inhibitors against 2,125 HCV NS3 inhibitors, requiring no paired
   compounds at all. It tests a necessary condition rather than transfer
   itself, and it inherits the same co-screening confound in a different
   shape, so it needs its own control arm measured before it is registered.
3. **A prospective assay.** Twelve compounds on two enzymes in one protocol
   would settle arm A better than the whole of ChEMBL does.

---

# Amendment 1 — the arms are partly a property of the source

2026-10-16, measured by `scripts/diagnose_source_scope.py` on chembl_37.
The prereg body above is unchanged. This amendment corrects one inference
in it and records what the correction rests on.

## What I got wrong

Point 3 of "What the evidence permits" explains control C1 — 49 paired
compounds across two proteases that share no fold — with "The COVID-era
papers that ran 3CLpro ran PLpro beside it". **That attribution is wrong.**

Per-year, C1's two ends are not COVID-era at all. `plpro` carries 8 rows
from 2020–2021 out of 809, and 493 from 2024. `3clpro`'s mass is 2022
(1,325), 2023 (932) and 2024 (1,403). C1 is **post-pandemic medicinal
chemistry running two protease assays inside one programme**, which is a
different thing from a 2020 repurposing screen and explains the same 47-of-49
single-document structure for a different reason.

I also predicted, from a top-fifteen table, that the coronavirus protease
data post-dates the repurposing window entirely. **It does not**: 3CLpro
carries 489 rows from 2020–2021, 10.3% of its total. The complete histogram
contradicted the truncated one, which is why the script prints every year
and a test asserts it.

## What the source measurement found

ChEMBL 37 holds 93,489 documents under 304 distinct journal titles. Five
titles supply **87.1%** of them; fifteen supply 98.2%; the remaining 289
titles supply 1.8%.

    J Med Chem             26,839     28.7%
    Bioorg Med Chem Lett   24,365     54.8%
    Eur J Med Chem         11,799     67.4%
    Bioorg Med Chem         9,438     77.5%
    J Nat Prod              8,945     87.1%

Against a lookup of nineteen named venues where repurposing screens are
published — PLoS ONE, Sci Rep, Nat Commun, eLife, Antiviral Res, bioRxiv and
others — the whole set contributes 2,353 documents, 2.5% of the database,
and **seven documents in 2020–2021**. PLoS ONE contributes 13 documents in
total and none from 2021, which is why the screen that is arm B's experiment
is absent.

Every one of this project's protease roles draws from the same five or six
medicinal-chemistry journals, in every year it has data for.

## What that changes, and what it does not

**It does not weaken any figure in the body.** The arms, the units and the
publication structure are unchanged, and the test the prereg registers is
unaffected.

**It narrows one claim.** The body says the catalytic-type gate is
"untestable on ChEMBL 37", which was correctly scoped and is now explained:
a medicinal-chemistry paper is a congeneric series — one scaffold, one
target — and a cross-target observation is not a thing that shape of paper
produces. ChEMBL has a great deal of coronavirus protease chemistry. What it
does not have is the same compound measured on two viral proteases, because
the experiments that do that are published in venues it does not curate.

**It does not say ChEMBL is deficient.** It answers "how does this series
vary on this target" and answers it extremely well. This project asks
"would this existing drug work on that virus". Those are different
questions, and the second one's evidence is largely outside this source.

**It raises the priority of the backlog's second screening panel.** The
NCATS CPE layer is already the only source here generated independently of
ChEMBL annotation, and `docs/evaluation.md` already asks for a second. The
measurement above is the reason: such a panel is not more data of the kind
we have, it is the only data of a kind this source structurally lacks.

## Arm A, in the same light

`3cpro` carries 157 rows, **none** from 2020–2021, spread over 2009–2018
with 63 rows each in 2016 and 2018. So arm A is not a COVID artifact either.
It pairs a pre-2019 picornavirus 3C series with coronavirus 3CLpro data —
a broad-spectrum protease programme, which is exactly the kind of paper that
would measure both, and exactly why nine of its twelve compounds are one
document.

---

# Verdict — NO RESULT, and two defects in the registration

2026-10-16, measured by `scripts/phase4_hcv_retrospective.py` on chembl_37,
after this file was committed. The decision rule registered above requires
P1 to be rejected at the PRIMARY identity. It is not.

    identity            n    active   rate    background   p (enrichment)
    exact InChIKey     11     7        63.6%   44.4%        0.1630   PRIMARY
    InChIKey skeleton  14    10        71.4%   44.0%        0.0364   sensitivity

**P1 is not rejected. Phase 4 returns NO RESULT.** Nothing is confirmed, the
catalytic-type gate remains an untested safety measure, and the paired arms
stay refused.

## The sensitivity crossed alpha, and it is not the result

The skeleton identity reads p = 0.0364. It was pre-declared as the
sensitivity and not the primary, in this file, before the numbers existed —
and the reason it is not promoted now is the reason it was written down
then: an identity chosen after seeing which one clears 0.05 is not an
identity, it is a result looking for a protocol.

It is also one compound wide:

    n=14, 10 active   p = 0.0364   <- observed
    n=14,  9 active   p = 0.1049
    n=11,  8 active   p = 0.0560
    n=11,  7 active   p = 0.1630   <- observed

A single compound moving from active to inactive takes the significant
reading to 0.105. Nothing should be built on a figure with that spread.

## Two defects in the registration, both mine

**The registered n was wrong.** This file's table says the broad set is 12
compounds at the exact identity and 17 at the skeleton. The test ran on 11
and 14. The definitions section defines the HCV end by the potency rules —
protein-level, relation `=` only, median within one assay type — and the
counts in the table were quoted from the survey's MEMBERSHIP line, which
applies none of them. The definition and the count in one document
disagreed, which is the same defect as the two floors one paragraph apart
that the survey's own verdict had.

It does not change the decision: both n clear the floor, the code
implemented the definition rather than the quoted count, and the primary is
not significant at either. But the registered sample size was wrong before
anything was measured.

**The power table did not span the real background.** It tabulated
backgrounds of 5% to 30% and the measured background is 44.4% — compounds
with a coronavirus measurement in ChEMBL are largely compounds somebody made
FOR a coronavirus, so about half of them work. At the real background the
test needs 9 of 11 (1.84x) or 10 of 14 (1.62x), not the "roughly three-fold"
this file registered. The claim was conservative in the wrong direction and
the table should have been computed against a background estimated first.

## What the direction means, and what it does not

The tilt is towards REFUTING the gate: 63.6% against 44.4%, in the direction
that says HCV NS3 inhibitors do carry coronavirus activity and the
catalytic-type refusal removes real signal. That is the uncomfortable
direction and it is recorded as such. It is not significant at the registered
rule and it is not a finding.

**An unregistered confound, noticed after the fact and labelled as such.**
An HCV NS3 inhibitor only enters this set if somebody tested it against a
coronavirus, and compounds get tested because they look promising. Selection
into the set correlates with the outcome being measured, which would inflate
the enrichment. This was not in the registration and cannot be corrected
retrospectively; a properly powered test would need a panel where the set
was chosen before the outcome was known.

## Decision

Per the rule above: **NO RESULT.** No wording in the project changes except
that the gate is described as untested, which this file already required.
Phases 5-10 gain no licence from phase 4 and lose none.

What would settle it is the same thing `docs/evaluation.md` already asks
for and `diagnose_source_scope.py` now explains: a screening panel that
measures a library against both enzymes, chosen before the outcome, and
published somewhere a database like ChEMBL would read.
