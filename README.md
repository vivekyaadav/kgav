# kgav — antiviral drug repurposing knowledge graph

A knowledge graph over the 16 viruses in `config/viruses.yaml` — seven
coronaviruses and nine flaviviruses — built to test one claim: that
assembling published facts into chains — *this compound inhibits that
protein; that protein is required by this virus* — will rank candidate
compounds for repurposing.

**It was tested against the controls that claim requires, and it did not
survive them.** That is the result this repository exists to report.

## What was measured, and what came out

Four pre-registered evaluations. Two are negative, and both publish the
effect size they had the power to detect.

| | finding |
|:---|:---|
| **the graph's own reasoning** | none of the declared metapath channels separates actives from measured-inactives on the assay-matched population once evidence sharing a publication with the label is withheld. The ensemble that does clear is matched by a baseline that counts a compound's edges, and most of its margin above chance is the asymmetry of what it reaches rather than how it ranks |
| **the pathway route** | no enrichment of a virus's dependency factors in Reactome pathways survives a degree-matched null, in either direction. Against a *uniform* null, several rows would have looked significant |
| **chemistry** | a scorer using no graph structure at all — the maximum ECFP4 Tanimoto a candidate reaches against that virus's known actives — does separate, in most registered cells, above a degree baseline that does not. It is a **resemblance** claim and not a hit rate: on the one assay-matched population the head of its ranking holds no actives once each candidate's nearest known actives are withheld |
| **prospectively** | trained before 2021 and tested after it, the chemical scorer clears for most evaluable viruses, **clears below chance for two**, and has no usable degree comparator anywhere — so none of it is shown to beat a baseline |

The one independent label axis — whether a compound reached a trial — is
**anti**-predicted by in-vitro resemblance and well predicted by how much
the compound has been studied. It measures something else, so it cannot
validate the chemistry.

Every figure behind those sentences, with its interval and its verdict, is
in [doc2/08-results.md](doc2/08-results.md); the argument is
[doc2/paper.md](doc2/paper.md); the weaknesses still in the *inference* are
[doc2/plan-v2.md](doc2/plan-v2.md).

## The register

16 viruses across two families, and the list below is checked against
`config/viruses.yaml` on every run — nine of them were once missing from
the front page of their own repository.

- **Coronaviridae (7)** — SARS-CoV-2, SARS-CoV, MERS-CoV, HCoV-229E,
  HCoV-OC43, HCoV-NL63, HCoV-HKU1
- **Flaviviridae (9)** — Zika, DENV-1, DENV-2, DENV-3, DENV-4, West Nile,
  JEV, Yellow fever, TBEV

## Layout

| | |
|:---|:---|
| `schema/kg_schema.yaml` | the frozen schema: 11 node classes, 24 edge classes, 8 metapaths |
| `config/viruses.yaml` | the virus register — taxa, aliases, proteomes, MONDO terms |
| `src/kgav/` | the validator, alias resolution, label construction, per-source ingest logic |
| `scripts/ingest_*.py` | one driver per layer |
| `data/releases/` | `nodes.jsonl` and `edges.jsonl` per layer, with a stamped MANIFEST |
| `data/results/` | every measured figure, tracked on purpose |
| `doc2/` | the documentation, generated |

## Run

```bash
pip install -e ".[dev]"
pytest -q
make gate     # planted-error release; the validator MUST reject it
```

## Why there is no count in this file

**A count in a README is pinned to one release and goes stale silently**,
which is the defect this project has spent most of its history finding in
its own documents — thirteen stale figures in one hand-drawn figure, a scope
statement wrong twice, a claim enumerating seven cells that described no run
that had ever happened. So the descriptive half of this file carries no
release figure at all, and `tests/test_readme.py` fails if one appears.

The live figures are one command away:

```bash
make audit           # per-layer, per-virus, what is built and what is missing
make kg-report       # regenerates the release census in doc2/kg_release_census.txt
make check-results   # every results file against the graph on disk
```

## The documentation

**`doc2/` is the documentation, and it is generated.** The prose lives in
`doc2/_templates/*.md.in`; every figure in it is derived at render time from
the schema, the register, the release MANIFEST, the artifacts in
`data/results/` and the constants imported from `src/kgav/`.

```bash
make doc2          # regenerate
make doc2-check    # fail if any document no longer matches the code
```

Where to start:

| | |
|:---|:---|
| [doc2/00-README.md](doc2/00-README.md) | the index |
| [doc2/06-methodology.md](doc2/06-methodology.md) | **if you are not a programmer, start here.** Every idea the results rest on, from first principles |
| [doc2/11-user-guide.md](doc2/11-user-guide.md) | how to drive the question-answering agent |
| [doc2/08-results.md](doc2/08-results.md) | every measured number, with its interval and its verdict |
| [doc2/09-limits.md](doc2/09-limits.md) | what this graph cannot do, each limit with the figure behind it |
| [doc2/paper.md](doc2/paper.md) | the manuscript: four controlled measurements |
| [doc2/plan-v2.md](doc2/plan-v2.md) | **what v2 is for.** Ten fixes to the inference drawn from those four measurements, each with the measurement that shows it is needed and the published claim it can overturn |
| [doc2/14-defect-ledger.md](doc2/14-defect-ledger.md) | the recurring defect shapes, named, with the instances that named them — the most reusable thing here |
| [doc2/registrations/](doc2/registrations) | the signed pre-registrations, amended by appending and never rewritten |
| [doc2/figures/](doc2/figures) | the architecture and the metagraph, each with an A4 PDF that prints on one page |

`doc2/kg_audit_report.txt` is the full specification: schema, architecture,
every data source and how it was used, the build pipeline, the metapaths,
how the graph supports pathway analysis, and the per-virus state. Its
release census is generated from the assembled graph and stamped with its
content hashes, so a stale block is detectable rather than merely
regrettable.

**`docs/` was retired on 2026-11-21** and everything in it moved into
`doc2/`: the registrations and the dated audit records verbatim, the
generated text files unchanged, the figures regenerated, and the five
hand-written documents superseded by generated ones.
`doc2/registrations/README.md` carries the full mapping, because a path
written inside a signed registration is not edited to follow a move.

## What the build actually taught us

**A DATED RECORD, NOT A DESCRIPTION OF THE RELEASE.** Every figure below is
as measured when the defect it describes was found, and is deliberately not
re-derived — the same convention `doc2/records/` follows, and the reason the
no-counts rule above stops at this heading. These are the five lessons that
generalise beyond this project.

**The dominant bug class is a discarded identifier.** Four times: pp1a chain
ids, single-chain ORF chain ids, MERS filed under an isolate taxon, STRING's
Ensembl keys. Every instance failed SILENTLY — no error, no missing-data
signal, just edges that never joined. Schema validation cannot see this: a
graph missing half its edges is structurally perfect. The rule now is that any
identifier the graph does not keep as a node must still RESOLVE, via SAME_AS.

**A layer must not assert a property only the merged graph can know.** This is
not the drift class above, where a stale file described a graph that had moved
on; every layer here was internally correct and current. `unquantified` meant
"this activity has no measured value". ChEMBL's coronavirus data holds 2 CC50
records in 21,900, so the chembl layer wrote true on every edge — a true
statement about that layer. The selectivity layer then supplied same-document
CC50 for 1,568 of those pairs, and assembly resolves qualifier conflicts by
layer precedence, which put chembl first. The layer HOLDING THE EVIDENCE lost
the merge: 12,756 edges claimed no measured value existed, not one of them
actually lacked one, and 1,264 carried a selectivity index alongside the claim
that none was computable. Nothing was stale, nothing disagreed with its
source, and no count looked wrong — the detector above cannot see this one.
The rule now is that such a property is DERIVED after merging, declared in the
schema as `derived_from_absence_of`, and a stated value that contradicts its
own edge is a validation error.

**Cross-entity count comparison is the primary detector.** MERS showing
mature=2 against siblings at 18-23. 229E and NL63 missing nsp5 while having
nsp12. Seventeen CRISPR screens tripping a sentinel at once. Every real bug in
this build was found by comparing counts across things that should be similar,
not by a test and not by the validator.

**BioGRID ORCS has no global sign convention.** 21 coronavirus screens score
positive = dependency, 11 score positive = restriction, and no metadata field
distinguishes them. A global rule inverts the biology for a third of them
silently, and a mis-signed restriction factor becomes a recommendation to help
the virus. Each screen is calibrated against anchor genes whose polarity is not
in doubt. This is a general obstacle for anyone aggregating ORCS, not a
project-specific quirk.

**A document is a figure holder, and figures drift.** Six degree lifts moved
across releases while the table quoting them sat unchanged; a summary column
named four rows the code refuses to read; two reports quoted a figure from a
graph that no longer existed. Every number in `doc2/` that a rebuild could
move is now re-derived from the results files by a test, and the release
census is generated rather than typed.

## Known limitations

These are the structural ones, stated without figures for the reason given
above. The measured ones — which viruses can be scored, which channels clear
chance, what the pathway test found, and what each null could have detected
— are in [doc2/09-limits.md](doc2/09-limits.md) and section 9 of the audit
report.

- **All three evaluation label predicates carry edges, and the newest of
  them contradicts the chemistry.** This file previously said "`TREATS` and
  `IN_TRIAL_FOR` are declared and empty, so every evaluation rests on
  `HAS_ANTIVIRAL_ACTIVITY_AGAINST` alone". The `indications` layer filled
  them from ChEMBL's `drug_indication`, registered in advance by
  `doc2/registrations/prereg_clinical_labels.md`. Measured on the registered
  COVID-19 cell, the structure statistic reads **below** chance while the
  degree control reads well above it: reaching a trial is predicted by how
  much a compound has been studied and anti-predicted by its resemblance to
  the compounds active in vitro. The axis is independent and it measures
  something else, so it cannot validate a structure scorer.
- **One of thirteen declared layers is not built.** `curated` needs
  published-screen TSVs whose provenance has to be established first, and it
  is what would give DENV-1, DENV-3, DENV-4, JEV and TBEV their first
  dependency factors. `gtop` WAS in this sentence and is built;
  `doc2/kg_release_census.txt` lists it among the assembled layers. The prose
  said "not built" for months after the edges landed, which is the drift this
  release's figures were re-derived to find.
- **Four of eleven node classes are empty** — `CellLine`, `Publication`,
  `Variant`, `VirusStrain` — so resistance and variant modelling have
  nothing to attach to, and assay context and publications exist only as
  edge qualifiers.
- **Half the register has no dependency factor at all**, so every
  host-directed channel is empty for those viruses by construction, and the
  pathway null is a negative result over part of a register rather than over
  a register.
- **Restriction-factor findings are NOT supported.** Screen-level
  calibration cannot separate the two tails within a screen, and
  interferon-pathway genes stay contested. The dependency side is validated
  against non-anchor controls; the per-gene tallies are in the audit report.
- **The virus-host bridge is mostly proximity and co-purification
  evidence**, not binary assays. Paths through it mean "same complex
  neighbourhood as a viral protein", not "binds", and every edge carries its
  detection method so the distinction survives into an answer.
- **Coverage is heavily skewed by study bias.** SARS-CoV-2 holds an order of
  magnitude more of everything than any other virus in the register, and two
  viruses have no curated virus-host interaction in any of the ten sources
  surveyed. On this graph, a compound's degree *is* study volume, which is
  why degree is the null hypothesis on every protocol and not a rival method.
- **Most drug–target edges carry no direction**, so "this drug targets a
  dependency factor" cannot always become "inhibiting it should block
  infection".
- **Labels are one measurement at one threshold.** A single cut on in-vitro
  concentration–response decides active from measured-inactive, and the
  majority of positives at that cut are not positive at a cut ten times
  lower. No result in this repository has been recomputed at the lower cut —
  [doc2/plan-v2.md](doc2/plan-v2.md) item C8.
- **The four endemic human coronaviruses share one MONDO term**, as do the
  four dengue serotypes. That is correct — the disease is one entity and the
  viruses are four — but it means disease-level joins cannot separate them.

## Data licences

Code is MIT. Data is not — see DATA_LICENSES.md. DrugBank is excluded by
default (CC BY-NC). The chemistry is ChEMBL; the full source list, with what
each contributed, is section 5 of `doc2/kg_audit_report.txt`.
