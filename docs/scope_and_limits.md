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

### The analysable set is 16.3% of the molecules

<!-- FIGURE: analysable_set -->
Of 25,383 small molecules, **4,130** carry a human `TARGETS` edge and can
be placed in a pathway. The rest have no human target in this graph:
**520** approved, **973** clinical, **35** phase 1, **95** preclinical and
**19,630** unknown.

A curated drug-target source covers the approved and clinical end, so
**1,493** is the ceiling on what one could add — not 21,253. The
remainder are bulk screening compounds no curated source describes.

**Every figure in this paragraph was wrong until 2026-11-10, and not
because anyone mismeasured.** It read 3,002 and 11.8%, with the partition
633 / 1,046 / 58 / 108 / 20,536 and a ceiling of 1,679 — all correct when
written against a graph whose only drug–target source was ChEMBL. The
GtoPdb layer then began contributing `TARGETS` edges (3,601 of its own and
3,164 shared with ChEMBL), 1,128 more molecules became analysable, and
nothing re-derived the paragraph. This is the rot that
`evaluation.md`'s *"The figures in this document are not under
`make check-results`"* section describes, caught in the act. All eight
figures now come from the release via a test, so the next layer that moves
them fails the suite instead of ageing them.

### Two rows in the overlap table carry no statistical weight

DENV-2 (50 dependency factors, 2 drugged, 7.5 expected) and West Nile
(18, 0, 3.1) are reported for scope, not evidence. A ratio computed on
an expectation of three is one Poisson draw. If the true ratio were 0.90
like the six adequately-powered taxa, observing 0 of an expected 2.8 has
probability about 0.06, and ≤2 of 6.75 about 0.08 — both unremarkable.

**They must not be read as depletion.** What they establish is that two
more viruses, in a second family, were examined.

### Nine nodes hold two compounds' screening results

<!-- FIGURE: tautomer_conflation -->
The node id is an InChIKey, and standard InChI's mobile-H layer treats
annular tautomers as one compound — `normalize_chemical` also canonicalises
the tautomer before keying. Every compound node therefore declares
`tautomer_collapsed: true`, and on **9 nodes that policy pools two distinct
NCATS library members**: 18 substances on 9 ids. The clearest pair differs
only in whether a guanidine is written `CCCN=C(N)N` or `CCCNC(N)=N`
(`NCGC00390264` and `NCGC00386242`).

**They are counted by the build, not estimated here.** `ingest_ncats`
reports `tautomer_conflated_nodes` and `tautomer_conflated_members` on every
run, and `make diagnose-tautomer-conflation` lists them.

**It is not fixable by changing the identifier, which is why it is stated
rather than repaired.** Canonical-tautomer SMILES would separate these nine
and split compounds every source treats as one; tautomers are frequently not
separable experimentally either. The alternative is a different identity for
the whole graph to resolve nine nodes out of 25,383, and the honest trade is
to declare the policy on every node and count the cases.

**What it costs a reader, measured rather than assumed:** on those nine, an
activity or inactivity call pools two molecules — and they are not off to one
side of the evaluation. **8 of the 9 carry a `HAS_ANTIVIRAL_ACTIVITY_AGAINST`
label and 5 carry a traversable edge**, so a channel can reach them;
`MQQNFDZXWVTQEH` holds 16 `TARGETS` edges and a similarity edge. Against
21,371 labelled compounds that is 8 labels, so no figure in
`docs/evaluation.md` moves measurably — but the right statement is the one
with the count in it, not the one that puts them outside the analysis.

**THE FIRST DRAFT OF THAT PARAGRAPH SAID THE OPPOSITE.** It read "None is
among the compounds any evaluated channel reaches, so no figure in
docs/evaluation.md depends on them", written from the assumption that nine
tautomer oddities would be peripheral. They are ordinary screened drugs.
Checked before the sentence shipped, and recorded here because the asymmetry
matters: the wrong version was the comfortable one.

**THE NUMBER TOOK FOUR READINGS TO GET RIGHT** (8, 11, 12, then 9), and the
record is fix-list L19. Two defects of method, not of chemistry: a classifier
that discarded any key containing a salt spelling, which made a subset report
more conflations than its superset; and counting NCGC sample *lots* as
distinct compounds, where 815 of 8,596 base ids appear under more than one
lot. Both were caught by two measurements disagreeing, neither by review.

### 4,839 label edges carry an extrapolated potency

ChEMBL's curators flag a value read outside the dose range they tested
("Falls Outside The Dose Series"). On release v0.1, 4,839 edges carrying a
label do: 4,691 inactive, 125 undecidable and **23 positives**. The number
on such an edge is an extrapolation past the concentrations measured, not a
reading from them.

They are kept, deliberately (fix-list L10, option (b)): the flag is a
declared qualifier, `outside_dose_series`, so a consumer can exclude them,
and `build_labels(..., dose_series_policy="exclude")` does it in one
argument -- costing 23 positives and 4,691 negatives. Excluding them by
default would delete evidence rather than describe it, which is the same
choice made for M8's weak fold mode. The cost of the choice is that every
AUC in `docs/evaluation.md` rests partly on extrapolated potencies.

Until 2026-11-09 the field was read by nothing (fix-list L9), so this limit
existed and was invisible.

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
| interactome | West Nile | 0.91 | 4 | 6 hits = **6.62x** | **no rejection region** |
| pathway | West Nile | 2.50 | 0 | 10 hits = **3.99x** | **no rejection region** |
| interactome | DENV-3 | 4.42 | 6 | 13 hits = **2.94x** | **no rejection region** |
| interactome | TBEV | 4.63 | 3 | 13 hits = **2.81x** | **no rejection region** |
| pathway | DENV-2 | 6.27 | 2 | 16 hits = **2.55x** | **no rejection region** |
| interactome | DENV-1 | 7.11 | 11 | 17 hits = **2.39x** | testable |
| pathway | HCoV-NL63 | 7.40 | 10 | 18 hits = **2.43x** | testable |
| interactome | DENV-2 | 35.38 | 40 | 56 hits = **1.58x** | testable |
| interactome | HCoV-HKU1 | 36.62 | 25 | 57 hits = **1.56x** | testable |
| interactome | HCoV-OC43 | 38.17 | 28 | 59 hits = **1.55x** | testable |
| interactome | HCoV-NL63 | 40.66 | 32 | 62 hits = **1.52x** | testable |
| interactome | Zika | 44.98 | 47 | 68 hits = **1.51x** | testable |
| interactome | HCoV-229E | 50.16 | 34 | 74 hits = **1.48x** | testable |
| pathway | HCoV-OC43 | 65.20 | 62 | 92 hits = **1.41x** | testable |
| pathway | MERS-CoV | 66.44 | 60 | 94 hits = **1.41x** | testable |
| interactome | MERS-CoV | 79.77 | 64 | 109 hits = **1.37x** | testable |
| pathway | Zika | 80.08 | 76 | 110 hits = **1.37x** | testable |
| pathway | HCoV-229E | 135.28 | 128 | 173 hits = **1.28x** | testable |
| interactome | SARS-CoV | 149.22 | 167 | 189 hits = **1.27x** | testable |
| pathway | SARS-CoV-2 | 358.43 | 384 | 419 hits = **1.17x** | testable |
| interactome | SARS-CoV-2 | 426.06 | 438 | 491 hits = **1.15x** | testable |

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
roughly 7x. (This count has now moved three times and the figure has never
been wrong by more than one row. It read *seven* until 2026-11-10, when a
rebuild made it eight and exposed the test beside it asserting the word
`Seven` — a check written with the same literal it was checking, so the
pair agreed with each other and not with `pathway_landscape.json`. The
2026-11-11 rebuild moved it back to seven, and **the row that crossed is
`interactome` DENV-1, from 2.52x to 2.39x** -- not the one the earlier note
pointed at. That note named `pathway` HCoV-NL63 as the eighth row at 2.45,
said the next one down was 1.59, and concluded "nothing about this sits on
a rounding boundary". Two of those three statements were wrong: DENV-1 at
2.52 was the next one down, and it was 0.12 from the boundary. HCoV-NL63
has since moved to 2.43 and is still above. **Two rows now sit within 0.05
of the line**, so this count should be expected to move whenever the graph
does. The count is derived from the file on every run, which is the only
reason any of this is visible.) Those are the rows where a real but moderate effect would be
invisible, and they are all flaviviruses or the small coronaviruses — the
taxa this release covers thinnest. The large-*m* rows at the bottom of
the table are sensitive to a 15–30% excess and found none, which is the
part of the null that carries weight.

That limit was compounded by resolution until 2026-11-11: at 2,000
permutations only **two** attainable p-values lay below 0.001136, so the
test could distinguish "significant" from "not" and little else. The run
now on disk uses the **100,000** permutations
`docs/prereg_pathway_overlap.md` fixes in advance, giving **113** attainable
values below alpha. No verdict changed and no p-value moved by more than
0.015, so what the extra resolution bought is not a different answer but
the right to state this one: the test was powered to reject and did not.

**The registration is signed retrospectively and the weight it carries is
correspondingly lower.** The permutation count, the correction family and
the drugged-proteome operand were all fixed before the run -- 9, 9 and 1
day(s) before it -- but the signature was added after it, against the
draft's own instruction not to run first. `docs/prereg_pathway_overlap.md`
opens with the ordering, commit by commit. It is not what a prospective
registration would be worth, and the null above is reported as a null
either way.

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
