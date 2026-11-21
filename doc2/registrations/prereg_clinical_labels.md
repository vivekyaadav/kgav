# Pre-registration: a label that is not an in-vitro assay

**SIGNED 2026-11-17, BEFORE THE INGEST AND BEFORE ANY AUC.** No score
under this protocol has been computed. `scripts/ingest_indications.py`
does not exist; `TREATS` and `IN_TRIAL_FOR` carry **0 edges** on the
assembled release, which `doc2/kg_release_census.txt` prints and
`tests/test_release_census.py` holds it to. Every operand counted below
was counted from `data/raw` and the release and is reproducible from
them; **no ranking, no interval and no verdict exists**, and
`tests/test_prereg_clinical_labels.py` fails if this document ever
claims one without the results file to carry it.

Date: 2026-11-17. Release: v0.1 (78,649 nodes / 395,646 edges,
`nodes_sha` dc6d76a2b3cfb2a1, `edges_sha` c0c4daccd48e9979).

## The defect this registers a fix for

`doc2/kg_release_census.txt` ends its predicate census with its own
verdict on this project's evaluations:

> 2 of 3 evaluation label predicate(s) carry no edges, so every
> evaluation rests on HAS_ANTIVIRAL_ACTIVITY_AGAINST alone.

`docs/scope_and_limits.md` states the consequence: *"There is no
independent label against which to ask whether a result is about the
biology or about that one source."* Every AUC this repository reports —
145 claim-carrying cells across the baselines, the hard-negative
protocols, the temporal split and the three scorer studies — is one
predicate, derived from in-vitro concentration-response measurements, at
one threshold. A scorer that recovered the idiosyncrasies of
concentration-response screening and nothing else would pass all of
them.

**AND ONE FIGURE IN THAT LIMITS SECTION WAS ALREADY STALE WHEN READ.**
It says *"13 of 16 register entries carry `mondo: null`, which is what
leaves the two labels above empty."* All 16 now carry a MONDO id, and
the census shows **10** `Disease` nodes, not three. The mapping is done;
what is missing is the drug-to-disease edges.

**AND THIS WAS ALREADY FOUND, TEN DAYS BEFORE THIS REGISTRATION.**
`doc2/records/data_audit_2026-11-07.md` records it in a table row: *"all 16
curated; 10 Disease nodes, 16 CAUSES. The stated reason
`TREATS`/`IN_TRIAL_FOR` are empty no longer holds — `drug_indication` is
simply never queried."* The limits section kept the superseded cause
anyway, and `docs/v1_release_plan.md` kept the item filed as *"blocked
behind the Disease axis"* on the strength of it. What let that survive
was the test beneath the sentence: it asserted both predicate names
appeared and that the document said every evaluation "rests on" one
label — both of which stayed true while the explanation went wrong.
**A presence check cannot fail on a wrong explanation.** Both documents
are corrected in the same change as this registration, and
`tests/test_scope_table.py` now derives the Disease count from the
release and the MONDO count from the register instead.

## The source, and why it costs nothing to fetch

ChEMBL 37's `drug_indication` table, in the SQLite distribution **already
on disk** at
`data/raw/chembl/chembl_37/chembl_37_sqlite/chembl_37.db`, stamped
`sha256 33c20374...c290d281` by `scripts/record_raw_provenance.py` on
2026-11-10 and already read by `scripts/ingest_chembl.py` for the
activity layer. 60,055 rows, each a (molecule, indication, phase) triple
with MeSH and EFO identifiers.

This is **not a new dependency and not a new fetch.** The provenance
record for that file is `reconstructed` — URL and retrieval date null,
checksum recorded after the fact — and that limitation is inherited here
unchanged rather than restated as if it were new.

## The join, measured before any score

Register MONDO id → `data/raw/mondo/mondo.json` xrefs → MeSH / EFO /
DOID → `drug_indication` → `molregno` →
`compound_structures.standard_inchi_key` → the release's
`INCHIKEY:` node identity.

**EVERY DISEASE IS SHOWN, INCLUDING THE FIVE THAT YIELD NOTHING.** A
table of the three that work would be the failure this project is built
against.

    MONDO          disease                      mols  w/structure  in release
    MONDO:0100096  COVID-19                      551          360         273
    MONDO:0005091  severe acute respiratory sy.  537          364         277
    MONDO:0005709  common cold                    43           39          21
    MONDO:0005502  dengue disease                 14            7           4
    MONDO:0019209  Japanese encephalitis           5            2           1
    MONDO:0002282  West Nile fever                 3            1           0
    MONDO:0018661  Zika virus infectious dis.      3            2           0
    MONDO:0100116  Middle East respiratory sy.     3            2           2
    MONDO:0017572  tick-borne encephalitis         0            0           0
    MONDO:0020502  yellow fever                    0            0           0

1,358 indication rows over the ten diseases, 612 distinct molecules, 405
with a structure in ChEMBL, **294 already `SmallMolecule` nodes** in the
release.

Edges that would land, `max_phase_for_ind` 4 → `TREATS`, 0.5–3 →
`IN_TRIAL_FOR`:

    IN_TRIAL_FOR   569
    TREATS           9
    total          578

## The family is ONE cell, and that is a measurement not a concession

**SARS AND COVID-19 ARE NOT TWO CELLS.** Their indicated-compound sets
share **354 of 364** — ChEMBL files a large share of COVID-era trials
under the SARS MeSH term `D045169`. Treating them as independent would
double-count one cohort and halve the alpha for nothing. The other eight
diseases yield 0 to 21 compounds in the release, none of them a testable
cell.

So: **one evaluable cell, COVID-19 against SARS-CoV-2. Family size 1,
alpha = 0.05, two-sided, z = 1.9600.** No Bonferroni correction is
applied because there is one test, and that is registered here so that a
second cell found later cannot be added without restating the alpha.

## The population, and the number that decides the design

Candidates are the compounds carrying a SARS-CoV-2 in-vitro label from
`kgav.labels.build_labels` at the registered 10,000 nM threshold — 1,768
positive, 11,048 negative, 12,816 in all. A compound in that population
was **screened**, so the absence of a trial is at least a measured
absence rather than an absence of attention.

Of the 360 COVID-19-indicated compounds with structures:

    in-vitro ACTIVE   vs SARS-CoV-2       22
    in-vitro INACTIVE vs SARS-CoV-2      200
    not labelled at all                  138

**200 OF THE 222 LABELLED IN-TRIAL COMPOUNDS WERE MEASURED INACTIVE.**
This is counted before any ranking exists and it drives two design
decisions:

1. The primary cell is the **whole labelled population**: 222 positives
   (in trial) against 12,594 negatives (screened, not in trial). The
   alternative — "among in-vitro actives, which reached a trial" — is
   22 positives, which no interval would resolve, and choosing it after
   seeing it is 22 would be choosing a population by its power.
2. The two axes are close to orthogonal in this graph. That is the
   independence the census said was missing, and it is **not flattering
   to the in-vitro axis**. It is registered as the expectation, so a
   null result below is a null and not a surprise to be explained away.

## Hypotheses

**P1 (primary).** The registered structure statistic — max Tanimoto to
admissible known actives, symmetric k = 10 exclusion, as registered by
`doc2/registrations/prereg_structure_scorer_symmetric.md` — separates in-trial from
screened-not-in-trial compounds: the 95% Hanley-McNeil interval on AUC
excludes 0.5.

**P2 (the fame control, and it is load bearing).** The `degree` channel
is computed on the same cell and reported beside P1. A compound in
COVID-19 trials is a compound that was noticed, published and edged.
**If `degree` matches or beats the structure statistic, P1 is reported as
a fame result and not a chemistry one** — the same rule the scorer
registrations apply, and here it is the likelier outcome rather than a
formality.

**P3 (independence of the label).** Reported on the 138 compounds
carrying no SARS-CoV-2 in-vitro label at all. These cannot enter P1's
population, so this is a separate statistic over a separate population
and is reported as an observation, **not as a test** — there is no
screened negative class for them.

**P4 (the mechanism stratum is NOT testable, registered as such).** Of
the 222 labelled in-trial compounds, **13** carry an `INHIBITS` or
`TARGETS` edge into any of the 63 SARS-CoV-2 proteins. A stratum of 13
is not a cell. It is counted here so that it cannot be introduced later
as though it had been planned, and so the mechanism confound below is
quantified rather than asserted.

## Hazards

**H1 — the indication label is mostly not an antiviral mechanism.** 13 of
222. COVID-19 trial registries contain immunomodulators, anticoagulants
and supportive care. A structure scorer trained toward direct-acting
antivirals has no reason to rank dexamethasone highly, and should not.
This bounds what P1 can mean in either direction and is the reason P1 is
a ranking claim about one cell rather than a claim about repurposing.

**H2 — reverse causation, and it cannot be controlled here.** A compound
entered COVID-19 trials in 2020 *because* of an in-vitro result, and that
same in-vitro result is in this graph as the label the scorer's
references are drawn from. `drug_indication` carries no per-row date this
repository can use, so **no temporal split is available** and no cutoff
can separate "the scorer predicts trials" from "trials followed the
assays". This is registered as an unresolved confound, not as a
limitation to be noted in passing: **a positive P1 does not license a
predictive claim**, and the verdict must say so in those words.

**H3 — the reference set contains the cohort.** 245 of the 370
COVID∪SARS indicated compounds also carry a SARS-CoV-2 activity edge, so
a similarity statistic would be scoring indicated compounds against
indicated compounds. The registered self-exclusion and the symmetric
k = 10 exclusion both carry over unchanged; the publication exclusion
does not apply, because the label does not come from a publication.

**H4 — presence-only labels.** A compound with no `IN_TRIAL_FOR` edge may
have been in a trial ChEMBL has not curated. The negative class is "no
curated indication", and the verdict must use that phrase and not "did
not reach a trial".

**H5 — one threshold, inherited.** The 10,000 nM activity cut and the
dose-series policy come from the existing label layer unchanged. Changing
either to improve P1 is out of scope for this registration.

## What would make this REFUSE rather than report

- Fewer than 150 in-trial compounds in the labelled population, or fewer
  than 1,000 negatives, on the release the run reads. The counts above
  are from v0.1 and a rebuild can move them; the run must recount and
  refuse rather than report a cell it did not register.
- A `TREATS`/`IN_TRIAL_FOR` edge whose object is not one of the ten
  register MONDO ids.
- Any compound entering the positive class that is not already a
  `SmallMolecule` node. The ingest **adds no compounds**; it only adds
  edges between nodes that exist. A new node would change every other
  population in the repository and would make `make auc-delta`
  uninterpretable for this release.
- The release's `nodes_sha` changing between the ingest and the run.

## What counts as a result

A null is a result. The expectation registered above is that `degree`
performs comparably to the structure statistic and that both sit near
chance, because H1 says most of the positive class is not an antiviral
cohort. **If that is what is measured, it is reported as the finding** —
that this repository's one independent label axis does not corroborate
its in-vitro rankings — and `docs/evaluation.md` says so in its claim
block.

What is not permitted: narrowing the population, moving the threshold,
splitting SARS from COVID-19, or promoting P3 or P4 to a test, after any
AUC exists.

## Enforcement

`tests/test_prereg_clinical_labels.py` asserts, from the release and
`data/raw` rather than from this document:

1. the ten MONDO ids are exactly the register's distinct non-null ids;
2. the SARS/COVID overlap is at least 90% of either set, so the
   family-of-one argument above is still true;
3. the in-trial/labelled counts in this document equal
   `build_labels` crossed with the ingest, each with its disease;
4. 13 of 222 for P4, so the untestable stratum cannot silently become
   testable;
5. this document claims no AUC while `data/results` holds no clinical
   label results file.

---

# First amendment, 2026-11-17: the ingest ran, and it placed 577 of 578

**APPENDED, NOT EDITED.** Everything above is what was signed before
`scripts/ingest_indications.py` existed. One registered operand came out
different and this records which, why, and what it cost.

Release after the ingest: v0.1, `nodes_sha` **dc6d76a2b3cfb2a1 —
unchanged**, `edges_sha` c0c4daccd48e9979 -> 873623c9a9df2049.

## The registered count was 578 and the layer placed 577

    registered    placed    predicate
           569       568    IN_TRIAL_FOR
             9         9    TREATS
           578       577

**ONE EDGE WAS REFUSED BY THE SCHEMA, NOT BY THE DATA.**
`IN_TRIAL_FOR.required_qualifiers` is `[nct_id, status]`, and one
COVID-19 (compound, disease) pair carries no `ClinicalTrials` row in
`indication_refs` — ChEMBL cites a DailyMed label or an ATC code for it
instead. The operand count above was taken from `drug_indication` before
the schema's requirement was checked against it, which is the defect:
**a count of what a source holds is not a count of what a schema
admits.** The edge is dropped and counted as `in_trial_without_nct`,
because the alternatives were to invent an identifier or to relax the
requirement, and either would make an unverifiable trial edge
indistinguishable from a verifiable one. 568 of 569 is the measured
yield of this join, and it is the number every figure below uses.

The other registered operands came out exactly as registered: 9 `TREATS`,
276 SARS, 4 dengue, 14 common cold, 1 Japanese encephalitis, 2 MERS, and
`compound_not_in_graph` dropping every indicated compound the graph does
not already hold.

## Two things the registration did not say, both measured now

**THE XREF MAP IS DERIVED FROM THE ONTOLOGY, NOT TYPED.** The
registration's table was produced from a hand-written
MONDO-to-MeSH/EFO map, and that map was wrong once in a way worth
recording: it put the COVID-19 MeSH term `D000086382` under
`MONDO:0100116` as well, because MERS-CoV's term carries no MeSH xref
and the nearest-looking one was reached for. That attributed **66
COVID-19 trials to MERS**. `kgav.indications.disease_xrefs` reads the
xrefs out of `mondo.json` instead and **refuses outright** if two
registered diseases ever claim one identifier. The derived map is wider
than the hand-written one — it picks up `Orphanet:` and extra `DOID:`
terms — and produces identical counts, which is what says the
substitution is sound rather than merely safer.

**`status` IS `unknown` ON ALL 568, AND THAT IS THE ENUM'S PURPOSE.**
ChEMBL records that an indication was studied at a phase, not whether
the trial recruited, completed or was terminated. `trial_status`
includes `unknown`; every edge carries it, so a consumer filtering for
`completed` gets nothing rather than something invented. `phase` carries
ChEMBL's `max_phase_for_ind` as a string, with 0.5 kept as `0.5`.

**46 pairs cite more than one trial.** `nct_id` is not an identity
qualifier, so one (compound, disease) pair is ONE fact; the lowest NCT
id is the representative and `support_count` carries the number of
registrations, so the second one is counted rather than lost. 10 pairs
are filed at several phases and the furthest is the edge — the others
are not contradicted, they are earlier.

## The release moved by exactly the new edges, and that was verified

The claim that this ingest adds no compounds is not an argument from the
code. The assembled release was built to a scratch directory first and
compared against the one on disk:

    nodes.jsonl   sorted, identical            78,649 both
    edges.jsonl   sorted, identical once the two new
                  predicates are removed      395,646 both
    nodes_sha     dc6d76a2b3cfb2a1  unchanged

So `make auc-delta` across this ingest reads a graph that differs from
the pre-ingest one by 577 label edges and nothing else, which is what
the registration's refusal conditions were for.

**A FULL REBUILD WAS DECLINED, AND THE REASON IS THE SAME CLAUSE.**
`rebuild_release.py --plan` reports all eleven existing layers as
"written by different code", because registering a layer edits
`pipeline.py` and `pipeline.py` is in every layer's source closure.
Re-running them would re-fetch BioGRID ORCS over the network, and a
newer screen index would change the graph in ways unrelated to this
ingest — which would confound the one comparison this label exists to
make. The scratch assembly was promoted instead, after the byte
comparison above. The layer stamps are therefore older than the code;
`indications` is unstamped for the same reason.

## One effect on existing numbers, named rather than discovered later

`TREATS` and `IN_TRIAL_FOR` are `held_out: true`, so neither is
traversable and `baselines.degree_unlabelled` — what the degree
baseline ranks on — skips both. The degree BASELINE is therefore
unmoved.

`baselines.Graph.degree` is a different counter and does include
held-out predicates, by a decision recorded in `baselines.py`. So the
294 compounds carrying an indication gain degree there, and DWPC
weighting reads it. The effect is 577 edges against 395,646 and is
expected to be invisible, but **it is not zero and it is not asserted to
be**: `make auc-delta` against the pre-ingest results is what measures
it, and whatever it reports is reported.

## Still no AUC

No ranking has been computed against these labels. The analysis
registered above — one cell, COVID-19 against SARS-CoV-2, 222 positives
against 12,594 negatives, `degree` as the load-bearing control — is
unrun, and `tests/test_prereg_clinical_labels.py` still fails if this
document acquires a score while no results file carries one.

---

# Verdict, 2026-11-17: the independent label does not corroborate the
# chemistry, and it is very well predicted by fame

**RUN ONCE, BY `scripts/clinical_labels.py`, ON THE PROTOCOL ABOVE.**
Release v0.1, `nodes_sha` dc6d76a2b3cfb2a1, `edges_sha`
873623c9a9df2049. Results: `data/results/clinical_labels.json`.

    channel         n+       n-      AUC     95% interval   clears alpha
    structure      221   12,595   0.4087   [0.3734, 0.4440]     YES
    degree         221   12,595   0.7751   [0.7385, 0.8116]     YES

| | prediction | outcome |
|---|---|---|
| **P1** | the structure statistic separates in-trial from screened-not-in-trial | **TRUE, AND IN THE WRONG DIRECTION** |
| **P2** | degree is the control; if it matches or beats P1, P1 is a fame result | **DEGREE WINS BY 0.37 AUC** |
| **P3** | the clinical-label-only compounds, reported not tested | **no signal: median 0.2059 against the negatives' 0.2069** |
| **P4** | the mechanism stratum is not a cell | **confirmed, 13 of 221** |

## P1 is true and useless, which the registration made sayable in advance

The interval excludes 0.5, so the registered prediction is literally
satisfied — **by anti-ranking.** Chemical similarity to the compounds
measured active in vitro against SARS-CoV-2 puts the compounds that
reached a COVID-19 trial **below** the compounds that were screened and
did not, 0.4087 with the whole interval under chance.

That is not a paradox and it is not a bug. **200 of the 221 positives
were measured INACTIVE in vitro**, which this document counted before
any ranking existed, and 13 of them carry an `INHIBITS`/`TARGETS` edge
into any SARS-CoV-2 protein. The in-trial class is largely approved
drugs repurposed on clinical reasoning — immunomodulators,
anticoagulants, supportive care — and a scorer built to resemble
direct-acting antivirals has no reason to rank dexamethasone highly and
does not. The registration said this in H1 and registered the
expectation; the measurement is the expectation, with a sign.

## P2 decides the reading, and it is not close

**Degree reads 0.7751 [0.7385, 0.8116]** on the identical population —
0.37 AUC above the structure statistic and the strongest baseline figure
anywhere in this repository. Whether a compound reached a COVID-19 trial
is very well predicted by how many edges it has, which is to say by how
much it has been studied.

**AND THIS NUMBER WOULD HAVE BEEN WRONG BY CONSTRUCTION A DAY AGO.**
`baselines.degree_ranking` ranks over the traversable graph only if its
caller passes `skip_predicates`, and the three registered scorer studies
never did — found the same day, recorded in the amendment to
`doc2/registrations/prereg_structure_scorer.md`. This run passes
`schema.held_out_predicates()`, which includes `TREATS` and
`IN_TRIAL_FOR`, so the fame control cannot read its own answer. A degree
baseline counting the clinical label edges would have been scoring the
positives on the fact that they are positives.

## So what this says about the project

The census line this registration was written against —
*"every evaluation rests on `HAS_ANTIVIRAL_ACTIVITY_AGAINST` alone ...
there is no independent label against which to ask whether a result is
about the biology or about that one source"* — now has an answer, and it
is a null for the chemistry:

- **The one independent label does not corroborate the in-vitro
  rankings.** It anti-correlates with them, for a reason that is about
  what clinical trials select for rather than about the scorer.
- **It is not independent of study volume.** Degree explains it far
  better than chemistry does, so this axis cannot be used to validate a
  structure scorer — it would reward the well-studied.
- **The in-vitro axis is not impeached by this.** The two axes measure
  different things, which is what near-orthogonality means, and the
  SARS-CoV-2 in-vitro result (0.8304 symmetric against a 0.5535 degree
  baseline) is untouched by it.

**WHAT IT DOES NOT LICENSE, in the registration's own words.**
`drug_indication` carries no usable per-row date, there is no temporal
split, and nothing here separates "the scorer predicts trials" from
"trials followed the assays the references are drawn from". **A positive
P1 would not have licensed a predictive claim, and an anti-predictive
one does not license its negation either.**

## The controls that were supposed to hold, and did

- **The symmetric exclusion is exactly symmetric**: `excluded_knn` =
  128,160 = 10 x 12,816 candidates, `knn_short` = 0. Every candidate of
  both classes lost exactly ten references.
- **Nothing was imputed**: 12,867 of 12,867 compounds fingerprinted, 0
  unparseable, 0 without SMILES, `empty_reference` = 0, so the floor
  that distorts the metapath channels does not arise here.
- **Self-exclusion fired on every reference**: `excluded_self` = 1,768.
- **The publication exclusion was off** (hazard H3) and is recorded as
  off: `excluded_publication` = 0.

## One registered operand came out one lower

**221 positives against 12,595 negatives, where this document
registered 222 and 12,594.** The missing compound is the single
`IN_TRIAL_FOR` edge the schema refused for carrying no ClinicalTrials
reference, recorded in the first amendment — so the one edge dropped at
ingest was one of this cell's positives, and it moved to the negative
class because it is screened and now carries no curated indication. The
refusal thresholds this document registered (150 positives, 1,000
negatives) are met with room, and `clinical_labels.py` recounts and
refuses rather than reporting a cell it did not register.


---

# Amendment, 2026-11-18: the reason this document gave for declining the rebuild was false

**NOT AN AMENDMENT TO THE PROTOCOL.** No registered decision changes and
no figure moves. This corrects a factual claim made above, in the
section *A FULL REBUILD WAS DECLINED, AND THE REASON IS THE SAME
CLAUSE*. A registration is appended to and never rewritten, so the
sentence stands where it was written and is wrong there.

## The sentence, and what is wrong with it

> Re-running them would re-fetch BioGRID ORCS over the network, and a
> newer screen index would change the graph in ways unrelated to this
> ingest — which would confound the one comparison this label exists to
> make.

`scripts/ingest_orcs.py`'s `fetch_screen` is **cache-first**: it reads
`data/raw/orcs/screens/screen_<sid>.json` when that file exists and only
calls the API for a screen it does not have. All **81** screens are
cached. Re-running the `orcs` layer alone reproduced its `edges.jsonl`
hash `68306c710ed7eff0` **byte-identical in 6.07 seconds**, with no
network call. Nothing would have been re-fetched, and no newer screen
index could have entered the graph.

**THE DECISION WAS STILL RIGHT. THE REASON WAS NEVER CHECKED.** What
actually made the rebuild a bad idea at that moment was that
`make auc-delta` compares the release on disk against `HEAD~1`, and
`HEAD~1` was a release built by older code — so the gate could not
distinguish a scorer regressing from a comparison spanning two label
populations. Promoting the scratch assembly after the byte comparison
preserved the one comparison this label was registered to make. That
argument needs none of the ORCS claim, which is why the false sentence
survived being read twice: it was decoration on a conclusion that held.

## What it cost, and the rule it earns

Nothing, in figures — and that is the point. An unchecked reason
attached to a correct decision is harder to catch than a wrong decision,
because the outcome vindicates it. It then becomes a cited fact: this
sentence was repeated in `docs/v1_release_plan.md`'s open-work table as
grounds for deferring the rebuild, so one unverified clause propagated
into a second document and sat in the way of a 25-minute job for two
sessions.

**A reason given for a decision is a claim, and carries the same
obligation as a figure.** If it is about what the code does, it is
checkable by running the code.

## The rebuild it was declining, now measured

Run on 2026-11-18. **25m29s.** The release came back byte-identical —
`nodes_sha dc6d76a2b3cfb2a1`, `edges_sha 873623c9a9df2049`, both
unchanged — and so did all 16 layer directories. `rebuild_release.py
--plan` now reports no layer "written by different code" and nothing
UNSTAMPED, so the staleness detector reads as a detector again.
`make auc-delta` PASSES: 157 claim-carrying cells, 0 stopped carrying a
claim, 0 started, largest fall +0.0000 against a floor of -0.0100. The
only diff in the tree was `MANIFEST.json`'s layer ordering.

The 577 edges this document registered therefore cost nothing in the
existing numbers, as the section *One effect on existing numbers* said
they were expected to and declined to assert.
