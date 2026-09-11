# kgav — antiviral drug repurposing knowledge graph

Coronaviridae v1: SARS-CoV-2, SARS-CoV, MERS-CoV, HCoV-229E, HCoV-OC43,
HCoV-NL63, HCoV-HKU1.

Layout: schema/kg_schema.yaml holds the frozen schema (10 node classes, 21 edge
classes, 6 metapaths); src/kgav/ holds the validator, alias resolution and
per-source ingest logic; scripts/ingest_*.py is one driver per layer;
data/releases/ holds nodes.jsonl and edges.jsonl per layer.

## Run

```bash
pip install -e ".[dev]"
pytest -q
make gate     # planted-error release; the validator MUST reject it
```

## Progress

| Day | Layer | Nodes | Edges |
|---|---|---|---|
| 3 | Spine: 7 coronaviruses, proteomes, mature peptides | 423 | 411 |
| 4 | Host: proteins, genes, pathways, STRING interactome | 52,763 | 317,382 |
| 5 | Virus-host bridge (VirHostNet) | — | 10,889 |
| 6 | CRISPR host factors (BioGRID ORCS) | — | 10,180 |

## What the build actually taught us

**The dominant bug class is a discarded identifier.** Four times: pp1a chain
ids, single-chain ORF chain ids, MERS filed under an isolate taxon, STRING's
Ensembl keys. Every instance failed SILENTLY — no error, no missing-data
signal, just edges that never joined. Schema validation cannot see this: a
graph missing half its edges is structurally perfect. The rule now is that any
identifier the graph does not keep as a node must still RESOLVE, via SAME_AS.

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

## Known limitations

- Restriction-factor findings are NOT supported in v1. Screen-level
  calibration cannot separate the two tails within a screen; interferon-pathway
  genes stay contested (IFNAR1 2:3, IRF9 2:2). The dependency side is validated
  against non-anchor controls: LY6E 1:7 restriction, NPC1 8:0, TMEM106B 7:0,
  PIK3C3 4:0, SCAP 6:1, SREBF2 3:0 dependency.
- The virus-host bridge is 98.7% proximity and co-purification evidence; only
  147 of 10,889 edges come from a binary assay. Paths through it mean "same
  complex neighbourhood as a viral protein", not "binds".
- Coverage is heavily skewed by study bias: SARS-CoV-2 has 7,815 bridge edges,
  the seasonal coronaviruses 270-405 each.
- The four seasonal coronaviruses have no MONDO disease term and are left
  unmapped rather than collapsed onto a generic "common cold" node.

## Data licences

Code is MIT. Data is not — see DATA_LICENSES.md. DrugBank is excluded by
default (CC BY-NC); the graph is built on DrugCentral and ChEMBL.
