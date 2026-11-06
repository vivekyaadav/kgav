# Data audit — every source, 2026-11-07

What each external database contributes to `data/releases/v0.1`, how much of it
survives each filter, what provenance is recorded, and what is not.

Two kinds of evidence, never mixed:

- **measured** — counted directly from `data/releases/v0.1/{nodes,edges}.jsonl`
  and `data/raw/`. Reproduce with the scripts named in each section.
- **read** — taken from the ingest code, cited as `file:line`. Nothing was
  re-ingested; a figure a script would print at the next run is named, not
  guessed.

Companion to `docs/audit_2026-11-03_fixes.md` (the fix list). New items found
here are numbered to extend its groups and are collected in section 9.

---

## 0. The release on disk is not the release the code would build

`data/releases/v0.1/MANIFEST.json` records `code.commit 2589e1d0…`, which is
**three commits behind `main` (ed1e9da)**. The fix in `1b24254` ("a
contradiction must survive the next witness") is not in the shipped graph.

Consequence, measured: the MANIFEST counts
`qualifier_irreconcilable:direction: 6`, and `grep -c contested edges.jsonl`
returns **0**. Those six `TARGETS` edges — two witnesses asserting opposite
directions — ship as `direction: "unknown"`, indistinguishable from the 10,389
edges that never had a direction at all. The contradiction exists only as a
counter in the manifest.

**Every number below describes the graph on disk.** A rebuild will move some of
them.

---

## 1. One row per database

Edge counts are per *component* source, so a merged edge carrying two sources is
counted under each; they therefore sum to more than 392,907.

| Database | Version recorded? | Raw bytes | Edges | Dating | % with a PMID |
|---|---|---|---|---|---|
| UniProt (human + 16 viral proteomes) | **no** | 315 MB | 176,907 | per-edge, entry date | 0% |
| STRING v12.0 | filename only | 139 MB | 141,019 | **one date** 2023-08-28 | 0% |
| ChEMBL 37 | path string only | 34 GB | 39,589 | per-paper year, `-01-01` | 33% |
| VirHostNet / IntAct / MINT / DIP / HPIDB | **no** | 1.5 GB | 12,514 | **one date** 2024-01-01 | 100% |
| BioGRID ORCS | **no** | 318 MB | 8,263 | per-screen year | 98% |
| NCATS OpenData | **no** | 13 MB | 7,432 | **one date** 2020-08-18 | 100% |
| Guide to PHARMACOLOGY | **parsed, then discarded** | 14 MB | 6,765 | **all 1970** | 86% |
| InterPro + UniProt (folds) | **no** | 404 KB | 182 | **all 1970** | 0% |
| MONDO | free text in YAML | 215 MB (two copies, **unread**) | 16 | all 1970 | 0% |
| RDKit (computed similarity) | **no** | — | 3,384 | **all 1970** | 0% |
| curated screens | n/a | **empty directory** | **0** | — | — |

**Measured:** 346,792 of 392,907 edges (**88.3%**) carry no publication
identifier. Everything traceable to a paper comes from ORCS, NCATS, the PPI
databases, GtoP, or the third of ChEMBL that has a `docs` row.

Nodes, by the namespace they were minted in:

| Class | Namespace | Count |
|---|---|---|
| SmallMolecule | INCHIKEY | 25,383 |
| Protein | UniProtKB | 20,888 |
| Gene | NCBIGene | 19,111 |
| Pathway | GO | 10,931 |
| Pathway | REACT | 2,305 |
| Gene | **KGAV** | **79** |
| TargetClass | **KGAV** | **39** |
| TargetClass | INTERPRO | 27 |
| OrganismTaxon | NCBITaxon | 16 |
| Disease | MONDO | 10 |

118 identifiers are locally minted and resolve to nothing outside this repo.

---

## 2. Version provenance, source by source

| Source | Where the version lives | Recoverable by a reader? |
|---|---|---|
| ChEMBL 37 | default path string, `scripts/ingest_chembl.py:35`, `pipeline.py:43` | only by reading the path |
| STRING v12.0 | hardcoded filenames, `scripts/ingest_host.py:51,57` | only by reading the code |
| GtoPdb | **parsed from the file's line 1** (`src/kgav/gtop.py:264-291`, "GtoPdb Version: 2026.3.1 — published 2026-10-12") and only `print`ed (`scripts/ingest_gtop.py:97-98`) | **no — thrown away** |
| UniProt | nowhere. `human_proteome.json.gz` has no version in its name and `version_in_filename: null` in its sidecar | **no** |
| BioGRID ORCS | nowhere; `screens_all.json` has 38 fields, none a release | **no** |
| NCATS | nowhere; only a hardcoded publication constant `PMID:33708112` (`src/kgav/ncats.py:67`) | **no** |
| InterPro | nowhere; two live REST endpoints (`src/kgav/folds.py:86-87`), cache has no date field | **no** |
| MONDO | free text inside `mondo_source` strings in `config/viruses.yaml` | partly |
| RDKit | nowhere | **no** |

`string_release_date = "2023-08-28"` (`scripts/ingest_host.py:22-28`) is a
**file mtime**, by its own comment — not a publisher release date.

Sidecars (`url`, `retrieved`, `sha256`) exist for `host`, `proteomes` and
`vhppi` only. For `host` and `proteomes` they are marked `RECONSTRUCTED`:
checksums added after the fact, with `url: null` and `retrieved: null` because
neither is recoverable. **Spine and host read no sidecar at build time**; only
`vhppi` does, and only for `source` and `taxon`. Nothing from any sidecar
reaches a node, an edge, or the MANIFEST.

`make check-raw` passes, and that is what passing means today.

---

## 3. Row accounting

### 3.1 Gates that drop rows with no counter

Each of these discards input and records nothing, so the loss cannot be
reported or regression-tested.

| Where | What is dropped |
|---|---|
| `src/kgav/host.py:119-120` | every GO molecular-function and cellular-component annotation (`P:` aspect only) |
| `src/kgav/host.py:129` | proteins with no GeneID xref or gene name — no Gene node, no `ENCODED_BY` |
| `src/kgav/host.py:154-158` | STRING alias rows that are not `UniProt_AC`, or whose accession is outside the reviewed proteome |
| `src/kgav/host.py:186` vs `:188` | rows below `experimental ≥ 400` — derivable only as `rows − passed_filter` |
| `scripts/filter_intact.py:113-115,146-147` | IntAct rows with <15 columns, and rows naming no declared taxon |
| `scripts/ingest_spine.py:310-312,327` | proteome entries discarded when curated accessions substitute for the file |
| `scripts/ingest_spine.py:383-384` | single-chain features with no `PRO_` id (the multi-chain branch counts the same condition at `:401`) |
| `src/kgav/chembl.py:390-391` (SQL) | every `standard_type` outside IC50/EC50/Ki/Kd, and every NULL `standard_value` |
| `src/kgav/host_targets.py:106-110` (SQL) | non-human, non-single-protein, pChEMBL < 6, `max_phase` < 2 — four gates, no counts |
| `src/kgav/selectivity.py:77-78,96-97` (SQL) | every censored EC50, and **every row not reported literally in `nM`** (no unit conversion) |
| `src/kgav/similarity.py:111-112` | above-threshold neighbours beyond `top_k=10` |
| `src/kgav/gtop.py:434-445` | the affinity value when its unit is unmapped — **the edge still ships**, with a direction and no number |
| `src/kgav/folds.py:504` | InterPro entries beyond 200 per chain (no pagination) |
| `src/kgav/orcs.py:310` | screens matching no configured virus — counted as `no_virus_match` but absent from the driver's printed list (`scripts/ingest_orcs.py:211-216`) |

Only the ChEMBL layer asserts that its accounting closes:
`emitted + dropped == rows`, exit 1 otherwise (`scripts/ingest_chembl.py:127-138`).
`hosttargets`, `ncats` and `selectivity` have no such check.

### 3.2 BioGRID ORCS: the screen funnel (measured over `screens_all.json`)

| Stage | Screens |
|---|---|
| Human screens in the index | 1,952 |
| Pass `is_virus_resistance_screen` (substring `"resistance to virus"` in the rationale, `src/kgav/orcs.py:175-183`) | **144** |
| Match a configured virus (11 ordered regexes, `:98-146`) | 100 |
| Pass both | **86** |
| …of which polarity resolves | 46 screen-level, 35 signed, **5 unresolvable → dropped** |

Per taxon, matched → usable: SARS-CoV-2 68→58, OC43 7→7, 229E 7→7, Zika 5→2,
MERS 4→4, DENV-2 3→3, NL63 3→3, West Nile 2→2, **yellow fever 1→0**.

**There is no score cutoff.** Hit-calling is `HIT == "YES"`
(`src/kgav/orcs.py:334`). `Score.1` is read only for its sign; ORCS's own
`SIGNIFICANCE_CRITERIA` and `SIGNIFICANCE_INDICATOR` fields are never read.

`config/orcs_polarity_overrides.json` contains `{}` — zero entries.

The module docstring's "41 of 96 coronavirus screens" (`src/kgav/orcs.py:21`) is
stale against the index; the "35 signed" figure still holds.

### 3.3 Dedup keys, and what they throw away

| Layer | Key | What is lost |
|---|---|---|
| `host_targets.py:180` | `(drug, target, standard_type)` | all assays and all papers for a pair collapse to one edge; the surviving potency is whichever row SQL returned first |
| `chembl.py:598,671,700` | includes `assay_chembl_id` | correct — but the assay id is then **not written to the edge** |
| `ncats.py:402-416` | normalized structure only | a later active plate **overwrites an already-emitted inactive edge in place**; no trace of the inactive measurement |
| `orcs.py:380-385` | `(protein, taxon, screen, direction)` | nothing at this stage; lost later in the merge (§6) |
| `similarity.py:84` | exact node-id identity only | nothing — but see §5 on same-molecule pairs |

---

## 4. Labels: what "has antiviral activity" actually means

**Measured** over the 26,276 `HAS_ANTIVIRAL_ACTIVITY_AGAINST` edges
(compound → virus taxon), applying the project's own 10 µM rule from
`src/kgav/labels.py:53-54`:

| Measurement | Edges | Share |
|---|---|---|
| Potent, below 10 µM | **6,126** | **23.3%** |
| `>` a value ≥ 10 µM — measured inactive | 17,568 | 66.9% |
| `=` a value ≥ 10 µM — inactive by the same rule | 2,547 | 9.7% |
| other censored relations ≥ 10 µM | 35 | 0.1% |

**76.7% of the edges named "has antiviral activity" are not active.**

Why, in code: the organism-level branch emits the predicate
**unconditionally and never calls `labels.classify`**
(`src/kgav/chembl.py:587-606`; NCATS the same by construction,
`src/kgav/ncats.py:429-431`). A verdict is consulted on exactly one branch —
compound → *viral protein* — where `active → INHIBITS` and
`inactive → MEASURED_INACTIVE_AGAINST` (`src/kgav/chembl.py:688-707`).

So `MEASURED_INACTIVE_AGAINST` (2,421 edges) is the inactive counterpart of
`INHIBITS`, not of `HAS_ANTIVIRAL_ACTIVITY_AGAINST`. **NCATS contributes 0 of
them**, and at the virus level — the level the evaluation uses — no inactive
result has been separated out at all. The verdict is recoverable only by
re-running `classify` on the qualifiers, which the evaluation does and the
agent and every export do not.

Per virus:

| Virus | Label edges | Potent < 10 µM | Potent share |
|---|---|---|---|
| SARS-CoV-2 | 14,996 | 2,306 | 15.4% |
| DENV-2 | 4,173 | 1,131 | 27.1% |
| Zika | 3,414 | 926 | 27.1% |
| Yellow fever | 1,030 | 143 | 13.9% |
| West Nile | 594 | 280 | 47.1% |
| DENV-1 | 349 | 225 | 64.5% |
| DENV-4 | 312 | 200 | 64.1% |
| DENV-3 | 311 | 216 | 69.5% |
| HCoV-OC43 | 257 | 187 | 72.8% |
| SARS-CoV | 250 | 150 | 60.0% |
| HCoV-229E | 212 | 124 | 58.5% |
| MERS-CoV | 183 | 142 | 77.6% |
| TBEV | 137 | 71 | 51.8% |
| JEV | 46 | 20 | 43.5% |
| HCoV-NL63 | 12 | 5 | 41.7% |
| HCoV-HKU1 | 0 | 0 | — |

The potent share rises as the dataset shrinks: the large sets are screening
campaigns reporting everything, the small ones are curated actives. Any metric
computed across viruses mixes two different sampling regimes.

**Reachability, measured:** 15,905 of 21,589 labelled compounds (**73.7%**) have
no edge besides their label.

### 4.1 Relation handling

`relation` is coerced, not validated:

```python
relation = r["standard_relation"] if r["standard_relation"] in VALID_RELATIONS else "="
```

`src/kgav/chembl.py:564` and `src/kgav/host_targets.py:155`. A NULL — or any
unrecognised string — becomes `"="`: a bound read as an exact value. There is
**no counter for the coercion**, and `stats["censored"]` is incremented *after*
it, so a coerced NULL is tallied as exact. Nothing on the edge records that the
relation was imputed, and the repo holds no measurement of how many ChEMBL rows
have a NULL relation.

### 4.2 No ChEMBL quality column is ever consulted

`confidence_score`, `data_validity_comment`, `potential_duplicate`,
`standard_flag` and `activity_comment` appear **nowhere in the repository**.
There is also no `assay_type` (B/F/A/T/P) filter. A flagged, duplicate or
low-confidence activity enters the label set exactly like a curated one, and no
assay-confidence value travels with the edge.

### 4.3 Selectivity pairs two unrelated experiments

`src/kgav/selectivity.py:93-99`: the cytotoxicity CTE is restricted by compound
and document only — **no cell line, no assay, not even the same organism**. An
antiviral EC50 from Vero cells can be divided by a CC50 from an unrelated human
line reported in the same paper. Three further defects on top:

- a censored (`>`) CC50 survives whenever its bound-derived ratio clears 10, and
  is then written as a point `cc50_nm` with `relation: "="` and
  `selectivity_verified: true` (`:273-294`), with nothing recording that it was
  a bound. All 3,159 edges carry both fields;
- the `cell_line` qualifier comes from `MIN(a.description)` (`:71`), an
  aggregate independent of `MIN(standard_value)` — the alphabetically first
  description in the group, which need not belong to the assay that produced
  the retained EC50;
- in the exact-value path (`:243-253`) the first exact measurement after a
  censored one is dropped from the median set, so `selectivity_n_studies`
  undercounts.

`src/kgav/ncats.py:332-337` is the one place in the project where a censored
cytotoxicity bound is correctly refused rather than imputed
(`cytotoxicity_not_detected: true`, no finite index written).

---

## 5. Dates, and the temporal split

### 5.1 Where each date comes from

| Source | Value | Code |
|---|---|---|
| UniProt-derived | `entryAudit.firstPublicDate` — the **protein entry** date, not the annotation date | `src/kgav/host.py:74`, `scripts/ingest_spine.py:171-174` |
| STRING | `2023-08-28` on every edge, from a download mtime | `src/kgav/host.py:211`, `scripts/ingest_host.py:22-28` |
| vhppi | `2024-01-01` on every edge | `src/kgav/vhppi.py:424` |
| ChEMBL / selectivity | `docs.year`, with month and day fabricated as `-01-01`; `1970-01-01` when the `LEFT JOIN docs` misses | `src/kgav/chembl.py:577`, `host_targets.py:188`, `selectivity.py:275` |
| ORCS | publication year parsed out of the `AUTHOR` string | `src/kgav/orcs.py:245-257` |
| NCATS | `2020-08-18`, the bioRxiv first-assertion date (deliberate) | `src/kgav/ncats.py:77` |
| GtoP | `1970-01-01`, unconditional | `src/kgav/gtop.py:375,486` |
| folds | `1970-01-01` | `src/kgav/folds.py:712,766`, `scripts/ingest_folds.py:926` |
| similarity | `1970-01-01` (deliberate) | `src/kgav/similarity.py:142,161` |
| spine `CAUSES` | `1970-01-01`, hardcoded | `scripts/ingest_spine.py:452` |

**The vhppi `--date` is not a fallback.** `scripts/ingest_vhppi.py:156` defaults
it to `2024-01-01` and its help text says "fallback first_asserted_date when a
row has no PMID year", but no branch anywhere derives a date from a PMID. It is
the value for every row. And the PMID is already in hand: it is passed as
`pmids=` on the *same* `em.edge()` call that sets the date
(`src/kgav/vhppi.py:424-425`). `2024-01-01` corresponds to no VirHostNet, IntAct
or IMEx release named anywhere in the repo.

**Measured:** 9,901 edges (2.52%) carry `1970-01-01`. 170,865 carry a single
release-date stamp. Only 222,042 (56.5%) have a per-edge date, and of those the
157,211 `PARTICIPATES_IN` edges are dated by their protein's entry date.

### 5.2 What survives a 2021 cut (measured)

| Source | Total | ≤ 2021 | After 2021 | Undated |
|---|---|---|---|---|
| UniProt | 176,907 | **176,669 (99.9%)** | 238 | 0 |
| STRING | 141,019 | **0** | 141,019 | 0 |
| ChEMBL | 36,425 | 13,191 (36%) | 20,657 | 2,577 |
| vhppi (all services) | 12,514 | **0** | 12,514 | 0 |
| BioGRID ORCS | 8,263 | 4,566 (55%) | 3,697 | 0 |
| NCATS | 7,432 | 7,432 (100%) | 0 | 0 |
| GtoP | 3,601 | 0 | 0 | **3,601** |
| similarity | 3,384 | 0 | 0 | **3,384** |
| InterPro | 182 | 0 | 0 | **182** |

The 2021 training graph is therefore: all the pathway annotations, all of NCATS,
a third of ChEMBL, and — depending on the `--undated` policy — the computed and
GtoP edges. **Zero protein–protein and zero virus–host interaction edges.** The
mechanistic backbone the metapaths M3 and M5 need is absent by construction, not
by evidence.

`src/kgav/temporal.py:49`: `COMPUTED_SOURCES = {"infores:kgav-computed"}` — one
member. All 182 folds edges carry `infores:interpro` and are therefore dropped
under `computed_only`, while the same module's docstrings (`:24-26`, `:56-57`)
name fold similarity as the paradigm computed edge.

---

## 6. The merge

Order (`src/kgav/assemble.py:73-75`), earlier wins a qualifier conflict:

`spine → host → vhppi → orcs → curated → chembl → selectivity → gtop → hosttargets → similarity → ncats → folds`

Run order is different and lives in `src/kgav/pipeline.py:64-107`: phase 1
builds the raw layers, phase 2 assembles, **phase 3 (selectivity, similarity,
gtop) reads the assembled release**, phase 4 assembles again. So gtop and
similarity are computed against a v0.1 that already contains a gtop and a
similarity layer, and `pipeline.py:151-162` records the related hazard: a
skipped phase-3 layer whose directory still exists is silently **reused**.

**Merge key** (`assemble.py:84-109`): `(subject, predicate, object)`, plus
`primary_knowledge_source` where the predicate sets `identity_includes_source`
(only `HAS_ANTIVIRAL_ACTIVITY_AGAINST`), plus the predicate's
`identity_qualifiers` (`HOST_FACTOR_FOR: [direction]`). This is what keeps
genuine contradictions alive as separate facts.

**`edge_merged: 28,490`** counts merge *events* — duplicate rows folded into a
fact already seen. A fact with five witnesses contributes four. 392,907 + 28,490
≈ the 421,397 rows read from the layer files.

**Qualifier conflicts outside `direction` are resolved by position and never
counted.** `resolve_qualifier` (`:220-261`) has a ranked ladder for exactly one
qualifier: `QUALIFIER_LADDER = {"direction": ("unknown", "modulator")}`
(`:176-178`). Every other qualifier takes the branch at `:247-248`, which
returns `hi` whenever it is not None — and `hi` is whichever edge arrived first:
across layers that is the precedence order, but **within one layer file it is
literally the earlier line of `edges.jsonl`**. So two rows disagreeing about
`ic50_nm`, `effect_size`, `cell_line`, `start`/`end` or `tanimoto_ecfp4` are
settled by row order, the loser is discarded, and no counter fires. This is the
quietest drop in the pipeline.

**Correction to the fix list (L2).** `src/kgav/labels.py:23` is cited there as a
claim that `assemble.py` contradicts. It is not about `assemble.py`: it
describes label construction, and it is accurate for what it describes —
`labels.py:168-177` computes `both = active & inactive`, counts
`ambiguous_dropped`, and removes the compound from both sets. `assemble.py`
drops nothing. The real defect is narrower and still real: a contradictory
*qualifier value* on a single fact has no safeguard outside the `direction`
ladder.

**ORCS witnesses collapse.** The layer writes 11,112 edges; the release holds
8,263. Because `direction` is identity-defining but the screen id is not, every
screen asserting the same `(protein, virus, direction)` becomes one edge, and
one surviving `cell_line` and one surviving `effect_size` stand for all
witnesses. Only `support_count` records how many there were. **The ORCS screen
id is never written to the edge** (`src/kgav/orcs.py:386-389`), so an edge
cannot be traced to the screen that produced it.

More generally: **no edge in any layer carries an assay or document
identifier.** ChEMBL's assay id is used in dedup keys and discarded;
selectivity's `doc_id` — the entire basis of its EC50/CC50 pairing — is never
emitted; NCATS's `sample_id` is resolved and used only for the tox join. No edge
in the graph can be traced back to the row that produced it.

### 6.1 ORCS direction vs effect_size

Confirmed and quantified. `direction` is computed through the calibration
multiplier (`src/kgav/orcs.py:366` → `:191-198`), while `effect_size` is the raw
`Score.1` with no multiplier (`:375-378`); both are written into the same
qualifier dict (`:386-389`). **Measured: 1,247 of 8,263 `HOST_FACTOR_FOR` edges
(15.1%) carry a direction that contradicts the sign of their own
`effect_size`** — 720 `restriction` with a positive score, 527 `dependency`
with a negative one. Any consumer reading the score's sign gets the biology
backwards on one edge in seven. (Fix list P4.)

---

## 7. Unused, duplicated and dead data

| Item | Size | Status |
|---|---|---|
| `data/raw/chembl/chembl_37_sqlite.tar.gz` + extracted `.db` | 5.7 GB + 30.4 GB | the same database twice; **92% of the raw footprint** |
| `data/raw/mondo/mondo.json` | 103 MB | read only by `scripts/verify_mondo.py`, which emits nothing |
| `data/raw/refs/mondo.json` | 103 MB | **no code reference anywhere**; a stale predecessor of the above |
| `data/raw/vhppi/intact.zip` | 1.35 GB | filtered once to a 139 MB working file, which is also kept |
| `data/raw/refs/sec_ac.txt`, `delac_sp.txt` | 49 MB | read by `tests/test_normalize.py` only |
| `data/raw/refs/taxdump.tar.gz` | 79 MB | `TaxonResolver` has no production caller |
| `data/raw/ncats/cpe.tsv` | 4 MB | deliberately never read — pools three protocols (`src/kgav/ncats.py:24-32`) |
| `data/raw/curated_screens/` | empty | **the layer is dead**: no `v0.1-curated`, absent from the MANIFEST, yet it holds a precedence slot and is documented as the source for DENV-1/3/4, JEV and TBEV host factors |
| `web/` | — | unused Manus template (fix list P11) |

**No layer resolves secondary or obsolete UniProt accessions.**
`normalize_bio.UniProtResolver` has **no caller outside the test suite** — the
fix list's P3 ("only `vhppi` resolves secondary accessions") overstates the
state. `vhppi` resolves *aliases* through `AliasIndex` (`SAME_AS` chains, 8-hop
cap); the host layer resolves nothing and silently discards any STRING alias
whose accession is outside the reviewed proteome.
`taxonomy_source: "infores:ncbi-taxonomy"` is declared at
`config/viruses.yaml:785` and carried by no edge.

**Same-molecule similarity edges: 434 of 3,384 confirmed.** `similarity.py:84`
tests node-id identity only. There is no skeleton check — while `gtop.py:182,222`
uses the same 14-char InChIKey skeleton as *proof two keys are one compound*.
The two layers disagree about what a shared skeleton means.

---

## 8. Documentation that the release contradicts

Beyond group D of the fix list:

| Claim | Measured truth |
|---|---|
| `scope_and_limits.md:130` 6,287 undated edges (1.6%) | **9,901** (2.52%) |
| `v1_release_plan.md:147` 150,432 `PHYSICALLY_INTERACTS_WITH` | **153,533** |
| `v1_release_plan.md:16` 386,192 edges | **392,907** |
| `scope_and_limits.md:210-212` "three Disease nodes; 13 of 16 `mondo: null`" | **all 16 curated; 10 Disease nodes, 16 CAUSES.** The stated reason `TREATS`/`IN_TRIAL_FOR` are empty no longer holds — `drug_indication` is simply never queried |
| `v1_release_plan.md:154` "`ingest_mitab` … discards" the PMID | it doesn't; all 12,514 edges carry one. Only the *date* ignores it |
| `v1_release_plan.md:64-65` host layer has "no recorded checksum" | sidecars with sha256 exist since 2026-10-25; URL and date are still null |
| `src/kgav/host.py:7-24` four row counts (474k at combined ≥ 700, 284k, 120k, ~309k) justifying the `experimental ≥ 400` cutoff | **no code in the repo computes any of them.** `ingest_string` never evaluates a combined-score threshold. The only checkable figure is the 141,019 edges emitted |
| `src/kgav/orcs.py:21` "41 of 96 coronavirus screens" | stale against `screens_all.json` (144 pass the rationale gate, 86 usable) |
| `scripts/filter_intact.py:14-20` "20,895 edges" | the layer contributes 12,514 |
| `config/viruses.yaml:386-770` "No MONDO term pinned" | nine entries whose comment sits directly above a pinned `mondo:` value |

Also confirmed, with specifics (fix list P5): `MONDO:0005709` is labelled
**"HCoV-HKU1 disease"** and shared by 229E, NL63, OC43 and HKU1;
`MONDO:0005502` is labelled **"DENV-4 disease"** and shared by all four dengue
serotypes. Six of sixteen viruses point at a disease node named after a
different virus.

---

## 9. New items for the fix list

Numbered to extend `docs/audit_2026-11-03_fixes.md`.

### Group V — version provenance (new group)

| # | Problem | Fix |
|---|---|---|
| V1 | GtoPdb's release is parsed from the file and only printed | Persist `file_version` into the layer and the MANIFEST (`scripts/ingest_gtop.py:97`) |
| V2 | No UniProt release recorded for any proteome | Record the release in each sidecar at fetch time; stamp it into the MANIFEST |
| V3 | No ORCS, NCATS, InterPro or RDKit version recorded | Write a sidecar at fetch time; for ORCS record the API query date |
| V4 | `string_release_date` is a download mtime | Use STRING's published release date, or rename the field to `downloaded` |
| V5 | MANIFEST records only layer SHAs and the code commit | Add a `sources` block: one entry per database with version, retrieval date and checksum |
| V6 | No edge carries an assay, document or screen identifier | Emit `assay_chembl_id`, `doc_id`, `screen_id` as qualifiers; they are already in hand at emit time |

### Extending group C — dates

| # | Problem | Fix |
|---|---|---|
| C9 | `PARTICIPATES_IN` is dated by the protein's UniProt entry date, so 99.9% of UniProt edges land in every training graph | Date the annotation, or exclude pathway edges from the temporal protocol and say so |
| C10 | ChEMBL and `TARGETS` dates fabricate `-01-01`, indistinguishable from a real 1 January; 3,989 carry the 1970 placeholder | Carry year-only dates as a year, and count the placeholder |

### Extending group L — labels

| # | Problem | Fix |
|---|---|---|
| L8 | No ChEMBL quality column is consulted anywhere (`confidence_score`, `data_validity_comment`, `potential_duplicate`, `standard_flag`, `activity_comment`) | Filter on `confidence_score`, drop flagged rows, carry the score onto the edge |
| L9 | The NULL → `"="` coercion has no counter and is tallied as exact | Count it; mark the edge `relation_imputed: true` |
| L10 | Selectivity's tox CTE is restricted by compound and document only | Restrict to the same cell line and assay; carry the CC50's `all_exact` flag |
| L11 | `selectivity.py:71` `MIN(a.description)` is independent of `MIN(standard_value)`, so `cell_line` may come from another assay | Select the description of the retained row |
| L12 | `host_targets` writes `assay_type: "binding"` on all 18,954 edges, including 1,183 EC50 rows | Derive it from `standard_type` |
| L13 | No `TARGETS` edge has a publication — `docs.pubmed_id` is not selected | Add the column (`host_targets.py:91-111`) |
| L14 | NCATS dedup overwrites an emitted inactive edge in place | Keep both measurements, or record the overwrite |

### Extending group E — statistics and reporting

| # | Problem | Fix |
|---|---|---|
| E6 | Non-`direction` qualifier conflicts are settled by row order with no counter | Add a ladder or a counter for measurement qualifiers |
| E7 | 14 gates drop rows with no counter (§3.1) | Add a counter per gate; extend the `emitted + dropped == rows` assertion to all layers |
| E8 | `similarity`'s `top_k=10` cut is unmeasured | Count the above-threshold pairs discarded |

### Extending group P — process

| # | Problem | Fix |
|---|---|---|
| P13 | The shipped v0.1 was built by `2589e1d`, three commits behind `main`; six contested directions ship as `unknown` | Rebuild, or record in the release notes which commit built it and what changed since |
| P14 | ChEMBL is stored twice (5.7 GB tarball + 30.4 GB db); two MONDO copies total 215 MB and one is referenced nowhere | Delete the tarball and `data/raw/refs/mondo.json` after checksumming |
| P15 | The curated-screen layer is empty but holds a precedence slot and is the documented source for five viruses' host factors | Transcribe the screens, or state in `scope_and_limits.md` that the layer is unbuilt |
| P16 | `orcs_polarity_overrides.json` is `{}`, so five usable screens are dropped | Resolve the five by hand, or record that they are unresolvable |
| P17 | GtoP and similarity read the assembled release, which already contains their own previous output | Document the fixed point, or build them from the phase-1 layers |
| P18 | Evidence tier is a per-layer constant; a weak ORCS hit and ACE2 carry the same tier | Derive the tier from the row, or document it as a layer label |

---

## Reproducing

```
make data-timeline          # acquisition dates, evidence dates, per-source dating mode
make data-audit             # schema declares vs release holds
make check-raw              # sidecars and checksums
python scripts/data_timeline.py | sed -n '/by source/,$p'
```

The two measurement scripts written for this audit are in the session
scratchpad; they read only `data/releases/v0.1` and take about 2 seconds each.
Their outputs are the tables in §1, §4 and §5.2.
