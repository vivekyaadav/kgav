# Scope and limits of release v0.1

What this graph covers, what it does not, and why — for the methods
section of the pathway-analysis paper and for anyone reading a number out
of this repository.

Every figure below names the command that produces it. The table is
checked cell by cell against the release by `tests/test_scope_table.py`,
which hard-codes no figure: a count that moves fails the suite until this
document is corrected, and a virus that appears in the release fails it
until this document mentions one.

The scope is most of the contribution. A negative result's value is
almost entirely in knowing what it was measured over.

## The sixteen viruses

Seven coronaviruses and nine flaviviruses, built from
`config/viruses.yaml` by `scripts/ingest_spine.py`. Counts from
`data/releases/v0.1`, reproduced by `make data-audit`.

<!-- TABLE: scope_by_virus — generated, checked by tests/test_scope_table.py -->

| virus        | taxon | family | dependency | restriction | antiviral | chain source | in overlap table |
|--------------|------:|--------|-----------:|------------:|----------:|--------------|------------------|
| SARS-CoV-2   | 2697049 | Coronaviridae | 3,460 | 1,122 | 14,996 | species | yes |
| HCoV-229E    | 11137 | Coronaviridae | 1,367 | 9 | 212 | species | yes |
| Zika         | 64320 | Flaviviridae | 815 | 0 | 3,414 | species | yes |
| MERS-CoV     | 1335626 | Coronaviridae | 678 | 0 | 183 | strain 1263720 | yes |
| HCoV-OC43    | 31631 | Coronaviridae | 657 | 6 | 257 | species | yes |
| HCoV-NL63    | 277944 | Coronaviridae | 57 | 24 | 12 | species | yes |
| DENV-2       | 11060 | Flaviviridae | 50 | 0 | 4,173 | strain 31634 | yes |
| West Nile    | 11082 | Flaviviridae | 18 | 0 | 594 | species | yes |
| Yellow fever | 11089 | Flaviviridae | 0 | 0 | 1,030 | strain 11090 | **no** |
| DENV-1       | 11053 | Flaviviridae | 0 | 0 | 349 | strain 408685 | **no** |
| DENV-4       | 11070 | Flaviviridae | 0 | 0 | 312 | strain 408689 | **no** |
| DENV-3       | 11069 | Flaviviridae | 0 | 0 | 311 | strain 408870 | **no** |
| SARS-CoV     | 694009 | Coronaviridae | 0 | 0 | 250 | species | **no** |
| TBEV         | 11084 | Flaviviridae | 0 | 0 | 137 | strain 70733 | **no** |
| JEV          | 11072 | Flaviviridae | 0 | 0 | 46 | strain 11073 | **no** |
| HCoV-HKU1    | 290028 | Coronaviridae | 0 | 0 | 0 | strain 443239 | **no** |

<!-- END TABLE: scope_by_virus -->

7,102 dependency and 1,161 restriction factors, 26,276 antiviral activity
edges. **8 of the 16 viruses have no dependency factor at all** and are
therefore absent from the overlap analysis entirely.

### Why each of the eight is absent

Not one reason but three, and the distinction matters:

- **No screen in this ORCS release.** DENV-1, DENV-3, DENV-4 and yellow
  fever. `src/kgav/orcs.py` matches dengue serotype 2, Zika, West Nile
  and yellow fever by pattern; of 1,952 human screens, 11 name a
  flavivirus and 7 pass `is_virus_resistance_screen`. A yellow fever
  pattern exists and no yellow fever screen survives the gate.
- **No genome-wide screen in the literature at all.** JEV and TBEV. This
  is a limit of the field. No ingest and no curation reaches it.
- **Coronaviruses with no screen reaching them.** SARS-CoV and HCoV-HKU1.

HCoV-HKU1 is the extreme case and is worth naming rather than leaving a
reader to find it in a node count: it carries **eight edges in the whole
graph** — its `BELONGS_TO` relations — and no host factor, no antiviral
activity, no PPI. It is in the register and answers no question.

## Limits that apply to every virus

### One strain stands for each flavivirus

Where the table says `strain N`, that virus's proteins come from one
strain's polyprotein, reassigned to the species taxon. "DENV-1 NS3" in
this graph is strain 408685's NS3.

This is what `chain_source` has always meant and it is not a defect: for
seven of the nine flaviviruses, UniProt's species proteome is a
TrEMBL-only placeholder with one protein and zero `Chain` features, and
those placeholders were retired on or before 2026-10-23. The mature
peptides are curated under a strain. But a chain-level claim about a
species is a claim about one isolate of it.

### Most drug–target edges do not say which direction

12,155 of 14,381 `TARGETS` edges state `direction: unknown` — **84.5%**.
ChEMBL's bulk activity table records potency; only `drug_mechanism`
carries an action type, for a minority of approved drugs.

So "this drug targets a dependency factor the virus needs" does **not**
become "inhibiting it should block infection" for five edges in six. The
landscape's own output says this, and it is the single largest
qualification on any mechanistic reading of a path.

`make data-audit` reports it; `src/kgav/gtop.py` exists to reduce it and
its input is not downloaded.

### The analysable set is 11.8% of the molecules

Of 25,383 small molecules, **3,002** carry a human `TARGETS` edge and can
be placed in a pathway. The rest have no human target in this graph:
633 approved, 1,046 clinical, 58 phase 1, 108 preclinical and 20,536
unknown.

A curated drug-target source covers the approved and clinical end, so
**1,679** is the ceiling on what one could add — not 22,381. The
remainder are bulk screening compounds no curated source describes.

### Two rows in the overlap table carry no statistical weight

DENV-2 (50 dependency factors, 2 drugged, 7.5 expected) and West Nile
(18, 0, 3.1) are reported for scope, not evidence. A ratio computed on
an expectation of three is one Poisson draw. If the true ratio were 0.90
like the six adequately-powered taxa, observing 0 of an expected 2.8 has
probability about 0.06, and ≤2 of 6.75 about 0.08 — both unremarkable.

**They must not be read as depletion.** What they establish is that two
more viruses, in a second family, were examined.

### What the evidence dates mean differs by source

`make data-timeline` separates three kinds:

- **per-edge** — ChEMBL (publication year), UniProt (entry first-public),
  ORCS (screen date)
- **one date on every edge** — STRING `2023-08-28`, VirHostNet
  `2024-01-01`, NCATS `2020-08-18`. These are *file release* dates, not
  when any fact was established.
- **undated** — 6,287 edges (1.6%) carry `1970-01-01`, the placeholder
  every ingest writes when it has none.

The 1986–2008 shape of the histogram is UniProt's curation history —
176,907 of 386,192 edges — not a record of when anything about these
drugs was learned.

## What the graph cannot be asked

`make data-audit` reports the declared schema against the built graph.
Four declared node classes have no nodes and ten predicates have no
edges. Two of those matter for interpretation:

- **`TREATS` and `IN_TRIAL_FOR` are empty**, and the schema marks both
  `evaluation_label`. Every evaluation in this project therefore rests on
  `HAS_ANTIVIRAL_ACTIVITY_AGAINST` alone. There is no independent label
  against which to ask whether a result is about the biology or about
  that one source.
- **Three `Disease` nodes for sixteen viruses.** 13 of 16 register
  entries carry `mondo: null`, which is what leaves the two labels above
  empty.

The resistance axis — `VirusStrain`, `Variant`, `CONFERS_RESISTANCE_TO`,
`OCCURS_IN`, `HAS_VARIANT`, `HAS_STRAIN` — is entirely absent. No claim
of the form "this mutation confers resistance to this drug" can be made
or checked here.

## Provenance

`make check-raw` verifies every raw input. As of 2026-10-25 the host layer
— the human proteome and STRING v12.0, 329 MB, feeding over
three-quarters of the edges — carries a reconstructed record: a checksum,
recorded after the fact, with the URL and retrieval date null because the
files were placed by hand and neither is recoverable. Eight coronavirus
proteomes are in the same state and not yet stamped.

A reconstructed record makes drift detectable from now on. It does not
make the build reproducible from a clean machine; that needs a fetcher
for those inputs, which does not exist.
