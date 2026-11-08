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
| SARS-CoV-2   | 2697049 | Coronaviridae | 3,460 | 1,122 | 14,772 | species | yes |
| HCoV-229E    | 11137 | Coronaviridae | 1,367 | 9 | 211 | species | yes |
| Zika         | 64320 | Flaviviridae | 815 | 0 | 3,401 | species | yes |
| MERS-CoV     | 1335626 | Coronaviridae | 678 | 0 | 183 | strain 1263720 | yes |
| HCoV-OC43    | 31631 | Coronaviridae | 657 | 6 | 257 | species | yes |
| HCoV-NL63    | 277944 | Coronaviridae | 57 | 24 | 12 | species | yes |
| DENV-2       | 11060 | Flaviviridae | 50 | 0 | 4,163 | strain 31634 | yes |
| West Nile    | 11082 | Flaviviridae | 18 | 0 | 591 | species | yes |
| Yellow fever | 11089 | Flaviviridae | 0 | 0 | 1,030 | strain 11090 | **no** |
| DENV-1       | 11053 | Flaviviridae | 0 | 0 | 349 | strain 408685 | **no** |
| DENV-4       | 11070 | Flaviviridae | 0 | 0 | 312 | strain 408689 | **no** |
| DENV-3       | 11069 | Flaviviridae | 0 | 0 | 311 | strain 408870 | **no** |
| SARS-CoV     | 694009 | Coronaviridae | 0 | 0 | 242 | species | **no** |
| TBEV         | 11084 | Flaviviridae | 0 | 0 | 137 | strain 70733 | **no** |
| JEV          | 11072 | Flaviviridae | 0 | 0 | 46 | strain 11073 | **no** |
| HCoV-HKU1    | 290028 | Coronaviridae | 0 | 0 | 0 | strain 443239 | **no** |

<!-- END TABLE: scope_by_virus -->

7,102 dependency and 1,161 restriction factors, 26,017 antiviral activity
edges. **8 of the 16 viruses have no dependency factor at all** and are
therefore absent from the overlap analysis entirely.

**"Antiviral activity edges" is the predicate's name, not a description of
its contents.** `HAS_ANTIVIRAL_ACTIVITY_AGAINST` means "this compound was
measured against this virus", and `kgav.labels.classify` over those 26,017
edges reads **20,239 inactive (77.8%)**, 5,232 active (20.1%) and 546 that
decide nothing. Three quarters of this predicate is evidence that a compound
does NOT work -- which is the most valuable part of it, because measured
inactives are the negatives every evaluation here depends on, and they are
also what makes the name misleading to anyone reading it from outside.

The inactives are not routed to `MEASURED_INACTIVE_AGAINST` because the cut
is a protocol parameter: `hard_negatives.py --inactive-above-nm` varies it
and records the value used. Routing by verdict at ingest would put an
evaluation threshold into the graph's topology. The schema's `doc` for the
predicate says all of this, and `tests/test_agent_activity_sentence.py`
re-derives the share from the release so the figure cannot rot.

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

10,376 of 17,982 `TARGETS` edges state `direction: unknown` — **57.7%**
on the assembled v0.1. It was 12,155 of 14,381 — **84.5%** — before the
gtop layer landed, and that pre-gtop figure was quoted here, and in six
other places, as though it described v0.1.
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
- **undated** — 14,838 edges (3.8%) carry `1970-01-01`, the placeholder
  every ingest writes when it has none. By source: VirHostNet 5,222, Guide
  to PHARMACOLOGY 3,601, `kgav-computed` 3,178, ChEMBL 2,522, InterPro 158,
  MONDO 16.

  **That count was 6,287 (1.6%) and the rise is a correction, not a
  regression.** Every VirHostNet edge used to carry a blanket
  `2024-01-01` applied to all 12,514 rows regardless of evidence
  (fix-list C1); dating them by the publication year in their own
  reference field left 5,222 with no parseable year, which are now
  honestly undated instead of confidently wrong. The same holds for the
  Guide to PHARMACOLOGY `TARGETS` edges (C3): they are undated because
  GtoP publishes no per-interaction year, and nothing in this graph can
  supply one offline — a PMID-to-year index built from the release is
  unsound, since 42.9% of PMIDs here get different years from different
  dated edges.

The 1986–2008 shape of the histogram is UniProt's curation history —
176,907 of 392,624 edges — not a record of when anything about these
drugs was learned.

## What the pathway null could have detected

**A null result is only a finding if the reader knows what would have
survived.** This table did not exist anywhere in the project until
2026-11-02, and "nothing survived correction" was being reported without
it. Failing to detect a 2x excess is not evidence against a 2x excess
when the test needed 7x.

Twenty-two tests in two channels, each also tested downward, so the
family is forty-four and the Bonferroni threshold is 0.05/44 =
**0.001136**. Minimum detectable effect is the smallest observed count
whose Poisson upper tail falls below that threshold, as a multiple of the
degree-matched expectation.

| channel | virus | expected | observed | min. detectable | depletion |
|---|---|---:|---:|---:|---|
| interactome | JEV | 0.83 | 1 | 6 hits = **7.26x** | **no rejection region** |
| interactome | West Nile | 0.87 | 4 | 6 hits = **6.91x** | **no rejection region** |
| pathway | West Nile | 2.55 | 0 | 10 hits = **3.92x** | **no rejection region** |
| interactome | DENV-3 | 4.41 | 6 | 13 hits = **2.95x** | **no rejection region** |
| interactome | TBEV | 4.59 | 3 | 13 hits = **2.83x** | **no rejection region** |
| pathway | DENV-2 | 6.40 | 2 | 16 hits = **2.50x** | **no rejection region** |
| pathway | HCoV-NL63 | 7.34 | 10 | 18 hits = **2.45x** | testable |
| interactome | DENV-1 | 7.15 | 11 | 18 hits = **2.52x** | testable |
| interactome | DENV-2 | 35.27 | 40 | 56 hits = **1.59x** | testable |
| interactome | HCoV-HKU1 | 36.71 | 25 | 58 hits = **1.58x** | testable |
| interactome | HCoV-OC43 | 38.25 | 28 | 59 hits = **1.54x** | testable |
| interactome | HCoV-NL63 | 40.73 | 32 | 63 hits = **1.55x** | testable |
| interactome | Zika | 45.09 | 47 | 68 hits = **1.51x** | testable |
| interactome | HCoV-229E | 49.94 | 34 | 74 hits = **1.48x** | testable |
| pathway | HCoV-OC43 | 65.29 | 62 | 92 hits = **1.41x** | testable |
| pathway | MERS-CoV | 66.39 | 60 | 94 hits = **1.42x** | testable |
| pathway | Zika | 80.37 | 76 | 110 hits = **1.37x** | testable |
| interactome | MERS-CoV | 79.69 | 64 | 109 hits = **1.37x** | testable |
| pathway | HCoV-229E | 134.91 | 128 | 173 hits = **1.28x** | testable |
| interactome | SARS-CoV | 149.51 | 167 | 189 hits = **1.26x** | testable |
| pathway | SARS-CoV-2 | 359.19 | 384 | 419 hits = **1.17x** | testable |
| interactome | SARS-CoV-2 | 426.23 | 438 | 492 hits = **1.15x** | testable |

**Six of the twenty-two downward tests have no rejection region at all.**
Under a Poisson background with mean *m* the most extreme downward
observation is zero, with probability exp(−*m*); where that exceeds
0.001136 — i.e. *m* < 6.78 — no observation whatsoever can clear the bar.
Those six tests are in the correction family and cannot contribute a
finding to it, which makes the sixteen that can harder to pass for no
return. Dropping them is defensible, because emptiness depends on the
matched expectation and not on any outcome, so it is not choosing tests
after seeing them — but the family belongs in a preregistration that does
not exist yet, and it is left alone until it does.

**Seven of the twenty-two need a 2.4x excess or more**, and two need
roughly 7x. Those are the rows where a real but moderate effect would be
invisible, and they are all flaviviruses or the small coronaviruses — the
taxa this release covers thinnest. The large-*m* rows at the bottom of
the table are sensitive to a 15–30% excess and found none, which is the
part of the null that carries weight.

A separate limit compounds this: at 2,000 permutations only **two**
attainable p-values lie below 0.001136, so the test can distinguish
"significant" from "not" and little else. `src/kgav/landscape.py` reports
this and states that 8,799 permutations would give ten.

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
