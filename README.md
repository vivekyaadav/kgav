# kgav — antiviral drug repurposing knowledge graph

16 viruses across two families.

- **Coronaviridae (7)** — SARS-CoV-2, SARS-CoV, MERS-CoV, HCoV-229E,
  HCoV-OC43, HCoV-NL63, HCoV-HKU1
- **Flaviviridae (9)** — Zika, DENV-1, DENV-2, DENV-3, DENV-4, West Nile,
  JEV, Yellow fever, TBEV

Layout: `schema/kg_schema.yaml` holds the frozen schema (11 node classes, 24 edge
classes, 8 metapaths); `config/viruses.yaml` is the virus register;
`src/kgav/` holds the validator, alias resolution and per-source ingest
logic; `scripts/ingest_*.py` is one driver per layer; `data/releases/` holds
nodes.jsonl and edges.jsonl per layer; `data/results/` holds every measured
figure, tracked on purpose.

## Run

```bash
pip install -e ".[dev]"
pytest -q
make gate     # planted-error release; the validator MUST reject it
```

## What the graph holds

**Not restated here.** A count in a README is pinned to one release and goes
stale silently, which is the defect this project has spent most of its
history finding in its own documents. The live figures live in two places:

```bash
make audit           # per-layer, per-virus, what is built and what is missing
make kg-report       # regenerates the release census in docs/kg_audit_report.txt
make check-results   # every results file against the graph on disk
```

`docs/kg_audit_report.txt` is the full specification: schema, architecture,
every data source and how it was used, the build pipeline, the metapaths,
how the graph supports pathway analysis, and the per-virus state. Its
release census is generated from the assembled graph and stamped with its
content hashes, so a stale block is detectable rather than merely
regrettable.

`docs/evaluation.md` is what the metapath scores actually separate, and is
mostly a record of things that did not hold up.

## What the build actually taught us

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
graph that no longer existed. Every number in `docs/` that a rebuild could
move is now re-derived from the results files by a test, and the release
census is generated rather than typed.

## Known limitations

These are the structural ones. The measured ones — which viruses can be
scored, which channels clear chance, what the pathway test found — are in
`docs/evaluation.md` and section 9 of the audit report.

- **Two of three evaluation label predicates carry no edges.** `TREATS` and
  `IN_TRIAL_FOR` are declared and empty, so every evaluation rests on
  `HAS_ANTIVIRAL_ACTIVITY_AGAINST` alone.
- **Two of twelve layers are not built.** `gtop` needs a Guide to
  PHARMACOLOGY download that requires registration; `curated` needs
  published-screen TSVs whose provenance has to be established first. The
  second is what would give DENV-1, DENV-3, DENV-4, JEV and TBEV their first
  dependency factors.
- **Four of eleven node classes are empty** — `CellLine`, `Publication`,
  `Variant`, `VirusStrain` — so resistance and variant modelling have
  nothing to attach to, and assay context and publications exist only as
  edge qualifiers.
- **Restriction-factor findings are NOT supported.** Screen-level
  calibration cannot separate the two tails within a screen;
  interferon-pathway genes stay contested (IFNAR1 2:3, IRF9 2:2). The
  dependency side is validated against non-anchor controls: LY6E 1:7
  restriction, NPC1 8:0, TMEM106B 7:0, PIK3C3 4:0, SCAP 6:1, SREBF2 3:0
  dependency.
- **The virus-host bridge is mostly proximity and co-purification
  evidence**, not binary assays. Paths through it mean "same complex
  neighbourhood as a viral protein", not "binds", and every edge carries its
  detection method so the distinction survives into an answer.
- **Coverage is heavily skewed by study bias.** SARS-CoV-2 holds an order of
  magnitude more of everything than any other virus in the register, and
  DENV-4 and Yellow fever have no curated virus-host interaction in any of
  the ten sources surveyed.
- **The four endemic human coronaviruses share one MONDO term**
  (MONDO:0005709), as do the four dengue serotypes (MONDO:0005502). That is
  correct — the disease is one entity and the viruses are four — but it
  means disease-level joins cannot separate them.

## Data licences

Code is MIT. Data is not — see DATA_LICENSES.md. DrugBank is excluded by
default (CC BY-NC). The chemistry is ChEMBL; the full source list, with what
each contributed, is section 5 of `docs/kg_audit_report.txt`.
