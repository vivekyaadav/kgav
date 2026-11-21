# 02 — The graph model

`schema/kg_schema.yaml`, version **0.16.0**, loaded by
`src/kgav/schema.py`. The schema's single job: *make it impossible for a
triple that is not declared to enter the graph, and impossible for an edge to
enter without provenance and a date.*

17 CURIE prefixes, 21 controlled vocabularies,
11 node classes, 24 declared
predicates, 8 whitelisted metapaths. Nothing validates by
convention: every `validate_*` function returns a list of `Violation` and
raises on nothing, so a caller decides whether to drop the row or fail the
batch.

## Node classes

7 of 11 declared classes carry
nodes on this release.

| node class | nodes | id prefixes |
|:---|---:|:---|
| SmallMolecule | 25,243 | INCHIKEY |
| Protein | 20,888 | UniProtKB |
| Gene | 19,190 | NCBIGene, KGAV |
| Pathway | 13,236 | REACT, GO |
| TargetClass | 66 | SCOP, CATH, PFAM, INTERPRO, SSF, KGAV |
| OrganismTaxon | 16 | NCBITaxon |
| Disease | 10 | MONDO |

A `SmallMolecule` id is an InChIKey, and three required properties record
**which key the id is** rather than anything about the compound:
`salt_collapsed`, `stereo_collapsed`, `tautomer_collapsed`. They are constant
across the graph, every layer must declare the same values, and
`tests/test_identity_policy.py` derives them from what `normalize_chemical`
actually does — because a policy stated in a property and a policy implemented
in a function are two things that can disagree.

`smiles_as_supplied` keeps the source's own structure string **only where
normalisation changed the substance** — a counterion stripped, a charge
neutralised, a proton moved — and not where it merely changed the spelling.
The first version of that rule compared against RDKit's canonical form and so
set the property on 98% of a layer, which says nothing.

## Edge classes

16 of 24 declared predicates carry
edges.

| predicate | edges | subject → object | flags |
|:---|---:|:---|:---|
| `PARTICIPATES_IN` | 157,211 | Protein → Pathway | — |
| `PHYSICALLY_INTERACTS_WITH` | 153,533 | Protein → Protein | symmetric |
| `HAS_ANTIVIRAL_ACTIVITY_AGAINST` | 29,781 | SmallMolecule → OrganismTaxon | held out, label |
| `ENCODED_BY` | 19,498 | Protein → Gene | — |
| `TARGETS` | 17,635 | SmallMolecule → Protein | — |
| `HOST_FACTOR_FOR` | 8,263 | Protein → OrganismTaxon | — |
| `INHIBITS` | 3,874 | SmallMolecule → Protein | — |
| `CHEMICALLY_SIMILAR_TO` | 3,035 | SmallMolecule → SmallMolecule | symmetric |
| `MEASURED_INACTIVE_AGAINST` | 2,420 | SmallMolecule → Protein | held out |
| `IN_TRIAL_FOR` | 568 | SmallMolecule → Disease | held out, label |
| `MEMBER_OF_CLASS` | 158 | Protein → TargetClass | — |
| `SAME_AS` | 119 | any → any | symmetric |
| `BELONGS_TO` | 79 | Gene → OrganismTaxon | — |
| `FOLD_SIMILAR_TO` | 24 | Protein → Protein | symmetric |
| `CAUSES` | 16 | OrganismTaxon → Disease | — |
| `TREATS` | 9 | SmallMolecule → Disease | held out, label |

### Two flags that are not the same flag

- `held_out` is about **traversal**: the edge must be stripped from the
  training graph. `HAS_ANTIVIRAL_ACTIVITY_AGAINST`, `IN_TRIAL_FOR`, `MEASURED_INACTIVE_AGAINST`, `TREATS`
- `evaluation_label` is about **measurement**: the edge *is* a label, a
  measurement about a (compound, virus) or (compound, disease) pair.
  `HAS_ANTIVIRAL_ACTIVITY_AGAINST`, `IN_TRIAL_FOR`, `TREATS`

Reading one as the other was a live bug. `MEASURED_INACTIVE_AGAINST` is held
out and is **not** a label: it terminates on a protein. A split that read
`held_out` as "this is a label" bucketed it by `e["object"]` and put viral
protein nodes where viruses belong — the reported negative counts rose
5,164 → 6,789 while positives did not move. So a label now declares itself,
`Schema.validate_labels()` requires every label to terminate on an
`OrganismTaxon` or a `Disease` **and** to be held out, and no predicate is
promoted to a label by being held out.

### The predicate whose name is wrong on purpose

`HAS_ANTIVIRAL_ACTIVITY_AGAINST` means *"this compound was measured against
this virus"*, active **or** inactive. Counted over the
29,781 edges carrying it on this release:

| verdict | edges | share |
|:---|---:|---:|
| inactive | 23,803 | 79.9% |
| active | 5,426 | 18.2% |
| undecidable | 552 | 1.9% |

So **79.9% of this predicate is evidence of inactivity**
and any consumer reading the name as a claim of activity is reading the
opposite of the measurement. The name is kept and the schema says so in
its own `doc:` field, because the alternative — routing inactives to a
different predicate at ingest — would bake an evaluation threshold into the
graph's topology, and the predicate would stop matching the verdict the moment
the threshold moved. `kgav.labels.classify` decides; the predicate name never
does.

### Identity: when two edges are one fact

`assemble.edge_key` keys a fact on `(subject, predicate, object)` plus:

- **identity qualifiers.** `HOST_FACTOR_FOR.direction` is identity-defining: a
  dependency edge and a restriction edge for one (protein, virus) are
  different facts, not two witnesses to one. 320 pairs are contested across
  screens — ACE2 is 34 dependency to 5 restriction — so the arbitrary winner
  was not a rare edge case.
- **the source**, where `identity_includes_source` is set. An NCATS panel
  inactive and a ChEMBL active for one compound must not collapse into one
  edge, or one real measurement is discarded before `build_labels` can see
  that the two disagree and drop the compound as ambiguous. Keying on
  `cell_line` would have split the merge the selectivity layer depends on;
  keying on the source keeps that merge and separates only the panel.

### Derived qualifiers

A qualifier declared `derived_from_absence_of: [...]` is **recomputed at
assembly and may never be ingested**. `unquantified` once had three
definitions in three modules; assembly merges by layer precedence, and
`chembl`'s "no CC50 exists" overwrote the selectivity layer's CC50 on 1,264
edges — the losing layer was the one holding the evidence. Measured on v0.1 at
the time: 12,756 edges claimed `unquantified: true` and **not one** of them
lacked a measured value. `validate_edge` now rejects an edge whose derived
qualifier contradicts its own evidence.

### Symmetry is declared, never inferred

`symmetric` is required on every class a metapath traverses, and
`validate_metapaths()` rejects a hop over a predicate that left it unset.
The walker used to guess: follow outgoing edges, fall back to incoming ones
only when a node had none. That works for a node with edges in one direction
and silently truncates every node with both — 9,016 host proteins carry
`PHYSICALLY_INTERACTS_WITH` in both directions and lost their inbound half,
costing M5 113 drugs and a third of its path mass.

A reverse marker (`<PREDICATE`) on a **symmetric** predicate is itself a
schema error: symmetric predicates are walked both ways by definition, so the
marker signals a misunderstanding rather than an intent.

## Metapaths

| code | what it is | path | constraints | hub exclusion |
|:---|:---|:---|:---|:---|
| **M1** | direct-acting | `SmallMolecule → INHIBITS → Protein → ENCODED_BY → Gene → BELONGS_TO → OrganismTaxon` | `INHIBITS.target_resolution` = `chain_from_assay_description` | no |
| **M2** | host-directed, single hop | `SmallMolecule → TARGETS → Protein → HOST_FACTOR_FOR → OrganismTaxon` | `HOST_FACTOR_FOR.direction` = `dependency` | no |
| **M3** | host-directed via virus-host PPI | `SmallMolecule → TARGETS → Protein → PHYSICALLY_INTERACTS_WITH → Protein → ENCODED_BY → Gene → BELONGS_TO → OrganismTaxon` | — | no |
| **M4** | pathway-mediated | `SmallMolecule → TARGETS → Protein → PARTICIPATES_IN → Pathway → <PARTICIPATES_IN → Protein → HOST_FACTOR_FOR → OrganismTaxon` | `HOST_FACTOR_FOR.direction` = `dependency` | yes |
| **M5** | host PPI-mediated | `SmallMolecule → TARGETS → Protein → PHYSICALLY_INTERACTS_WITH → Protein → HOST_FACTOR_FOR → OrganismTaxon` | `HOST_FACTOR_FOR.direction` = `dependency` | yes |
| **M6** | cold start via chemical similarity, host-directed | `SmallMolecule → CHEMICALLY_SIMILAR_TO → SmallMolecule → TARGETS → Protein → HOST_FACTOR_FOR → OrganismTaxon` | `HOST_FACTOR_FOR.direction` = `dependency` | no |
| **M7** | cold start via chemical similarity, direct-acting | `SmallMolecule → CHEMICALLY_SIMILAR_TO → SmallMolecule → INHIBITS → Protein → ENCODED_BY → Gene → BELONGS_TO → OrganismTaxon` | `INHIBITS.target_resolution` = `chain_from_assay_description` | no |
| **M8** | cross-family direct-acting via conserved fold, catalytic type gated | `SmallMolecule → INHIBITS → Protein → FOLD_SIMILAR_TO → Protein → ENCODED_BY → Gene → BELONGS_TO → OrganismTaxon` | `FOLD_SIMILAR_TO.catalytic_type_match` = `True`; `INHIBITS.target_resolution` = `chain_from_assay_description` | yes |

Two constraints are load-bearing rather than cosmetic:

- **`INHIBITS.target_resolution = chain_from_assay_description`** on M1, M7
  and M8. ChEMBL has no chain-level coronavirus targets: every activity maps
  to the whole replicase polyprotein, and `chembl.py` recovers the specific
  nsp from the assay *description* where it can. Without the constraint, DWPC
  did not merely tolerate unresolved edges, it **preferred** them: rank is
  inversely monotone in target degree, and a parent accession's degree is low
  *because* resolution succeeded for everything else. Every improvement to
  resolution lowered the parent's degree and raised the score of whatever
  still failed — a denominator that rewards its own gaps.
- **`FOLD_SIMILAR_TO.catalytic_type_match = true`** on M8. Coronavirus,
  picornavirus and flavivirus proteases share SSF50494 and
  G3DSA:2.40.10.10 — by fold they are indistinguishable, and only the
  nucleophile separates them. Without the gate, M8 walks from SARS-CoV-2 Mpro
  to HCV NS3 and recommends nirmatrelvir for hepatitis C.

`fold_strength: strong` is deliberately **not** a constraint: of 36 allowed
cross-family pairs, the 20 strong ones are protease pairs whose chains the
spine does not build, so all four hops M8 can actually walk are motif/weak.
Adding the line would not tighten M8, it would delete it.

## Retrieval limits

| limit | value | why |
|:---|---:|:---|
| `hub_exclude_classes` | Protein | compounds top out at degree 9 and need no cut; viral proteins are destinations, never hubs |
| `hub_degree_cut` | 460 | the `Protein` p99.9 on this release |
| `max_pathway_members` | 200 | a GO term with 1,389 members generates ~1.9M protein pairs in M4 |
| `max_neighbors_per_node` | 50 | an unbounded expansion from a well-connected protein fills a context window with ubiquitin |
| `max_subgraph_edges` | 100 | the same, for a retrieved subgraph |

Degree on this release: median 2, p99 127,
p99.9 400, max 19,498. A **global** degree cut would have
excluded the virus nodes and the viral replicase, which is where M1–M6 *end*.
The real connectors are the ubiquitin/ribosomal cluster.

## The viruses

| name | taxon | family | polyprotein | chains | chains from a strain | disease |
|:---|:---|:---|---:|---:|:---|:---|
| SARS-CoV-2 | NCBITaxon:2697049 | Coronaviridae | P0DTD1 | 16 | no | MONDO:0100096 |
| SARS-CoV | NCBITaxon:694009 | Coronaviridae | P0C6X7 | 15 | no | MONDO:0005091 |
| MERS-CoV | NCBITaxon:1335626 | Coronaviridae | K9N7C7 | 15 | yes | MONDO:0100116 |
| HCoV-229E | NCBITaxon:11137 | Coronaviridae | P0C6X1 | 15 | no | MONDO:0005709 |
| HCoV-NL63 | NCBITaxon:277944 | Coronaviridae | P0C6X5 | 15 | no | MONDO:0005709 |
| HCoV-OC43 | NCBITaxon:31631 | Coronaviridae | P0C6X6 | 15 | no | MONDO:0005709 |
| HCoV-HKU1 | NCBITaxon:290028 | Coronaviridae | P0C6X2 | 15 | yes | MONDO:0005709 |
| Zika | NCBITaxon:64320 | Flaviviridae | Q32ZE1 | 13 | no | MONDO:0018661 |
| DENV-2 | NCBITaxon:11060 | Flaviviridae | P29990 | 13 | yes | MONDO:0005502 |
| DENV-1 | NCBITaxon:11053 | Flaviviridae | P27909 | 14 | yes | MONDO:0005502 |
| DENV-3 | NCBITaxon:11069 | Flaviviridae | P27915 | 13 | yes | MONDO:0005502 |
| DENV-4 | NCBITaxon:11070 | Flaviviridae | Q2YHF2 | 13 | yes | MONDO:0005502 |
| West Nile | NCBITaxon:11082 | Flaviviridae | P06935 | 13 | no | MONDO:0002282 |
| JEV | NCBITaxon:11072 | Flaviviridae | P27395 | 13 | yes | MONDO:0019209 |
| Yellow fever | NCBITaxon:11089 | Flaviviridae | P03314 | 13 | yes | MONDO:0020502 |
| TBEV | NCBITaxon:11084 | Flaviviridae | Q01299 | 13 | yes | MONDO:0017572 |

The register (`config/viruses.yaml`, validated by
`src/kgav/virus_register.py`) currently reports
**0 violations**. Its central rule: *a virus whose
mature-peptide chains live under a descendant taxon must declare that taxon.*
UniProt's `taxonomy_id:` search includes descendants, so querying a species
silently returns whatever strain happens to be curated underneath it — MERS
lost Mpro and RdRp that way, HKU1 the same. Seven of nine surveyed
flaviviruses are curated under a strain taxon, so a lineage-walking lookup
would have "worked" for all of them and left no record of which strain it
picked.
