# 09 — The agent

`src/kgav/agent/`, `scripts/agent_server.py`, `scripts/kgav_ask.py`,
`web/chat.html`. Five modules, one HTML file, no build step.

```
make agent-ui                                  # http://127.0.0.1:8765
make agent-ask Q="what is known about DENV-2?"
python scripts/kgav_ask.py --no-model "..."    # facts only, no synthesis
```

## Why it exists in this shape

The evaluation produced a specific, measured account of what this graph can
and cannot support, and **an agent that ignored it would be worse than no
agent**. The failures are not random: they are systematic, plausible-looking,
and in one case actively harmful. So the agent is five layers and the model
touches only the last one.

```
question
  │
  ├─ guards.classify_intent      keyword-deterministic routing
  ├─ orchestrator.build_brief    entity resolution, the refusal, the carry
  ├─ tools.*                     fixed functions, bounded output, verbalised
  ├─ guards.check                the verdict and the mandatory caveats
  └─ verify.verify_response      every sentence checked against the facts
```

## 1. Intents are classified by keyword, not by a model

`explain`, `triage`, `evidence`, `profile`, `pathway`, `prospect`, `capability`, `method`, `unknown`

*Deliberately keyword-based rather than model-based.* A language model asked
to classify intent will occasionally route a prospective question into a
permitted category, and the whole point of this layer is that **the refusal
must not be probabilistic**.

The ordering inside the classifier is load-bearing. The candidate-seeking
backstop runs **before** the permitted branches, because it used to run last
and appending any permitted keyword to a prospective question defeated it
entirely:

```
which drugs could work against Zika and what is the mechanism?   -> EXPLAIN
suggest compounds for JEV and rank them                          -> TRIAGE
which inhibitors could be repurposed for yellow fever, which process -> PATHWAY
```

That is the exact failure the backstop exists to prevent, reached by one extra
word, and the answer came back under full caveats that make it read as vetted.

`explain`, `pathway`, `triage` are the kinds that need a route to the virus
at all — a question of one of those kinds about a virus with no route is
refused rather than answered with an empty ranking.

## 2. Refusals are scoped, not categorical

The prospective refusal says the graph *cannot rank **novel** chemical
matter*, and names the viruses where a ranking **is** licensed. It used to say
"this graph cannot support prospective candidate generation", full stop —
while the same agent answered "name a drug for yellow fever" with twelve
ranked approved and clinical compounds. Both turns were defensible and
together they read as an arbitrary gate.

A virus with **no route at all** gets a different refusal, and the difference
matters: an empty ranking reads as *no compound looks promising*, when the
truth is that no compound can be reached. The refusal names the counts —
antiviral measurements held out as labels, and the absence of any
protein-level inhibitor, CRISPR host factor or virus-host interaction.

A refusal also carries its caveats **in its text**, because a refusal reaches
the user as `brief.refusal` and nothing else does.

## 3. The caveats are labelled, and verification enforces the label

20 caveats. The **label** is what verification checks, so the prose
below it can be rewritten freely without changing the requirement:

| key | label a synthesised answer must carry |
|:---|:---|
| `cell_context` | **CELL CONTEXT** |
| `host_directed_route` | **HOST-DIRECTED ROUTE** |
| `cell_line` | **CELL LINE** |
| `selectivity_unknown` | **SELECTIVITY UNKNOWN** |
| `selectivity_cytotoxic` | **CYTOTOXIC** |
| `selectivity_selective` | **SELECTIVE** |
| `proximity_evidence` | **PROXIMITY EVIDENCE** |
| `evidence_tier` | **EVIDENCE** |
| `virus_not_validated` | **VIRUS NOT VALIDATED** |
| `virus_calibrated` | **CALIBRATION** |
| `virus_assumed` | **VIRUS ASSUMED** |
| `contradicted` | **CONTRADICTED** |
| `ranking` | **RANKING** |
| `pathway_route` | **PATHWAY ROUTE** |
| `target_direction` | **TARGET DIRECTION** |
| `thin_route` | **THIN ROUTE** |
| `labels_unreachable` | **LABELS UNREACHABLE** |
| `held_out_label` | **ANSWER SHEET** |
| `entities_omitted` | **NAMES OMITTED** |
| `carried_forward` | **CARRIED FORWARD** |

The trap this closes: verification used to derive its key from the caveat's
own text — the part before the first colon, or the first 40 characters when
there was no colon. So editing any caveat's wording silently changed what was
required, without touching the verifier. Editing a label is now a deliberate
act; editing prose is free.

`held_out_label` is deliberately **not** the same string as the per-line
marker the tools put on each held-out fact. Verification enforces a caveat by
searching the answer for its label, so if the two matched, a model that merely
copied a marked fact line would satisfy the caveat without ever stating what
it means — enforcement that can never bite.

## 4. The biology claims carry sources; the measured ones carry none

Every performance figure in the agent is **rendered from a calibration** at
answer time, because a hardcoded figure goes stale: `TEMPORAL_AUC = 0.500` and
`ABSENT_FRACTION = 0.77` were both deleted from the source, and the second was
worse than stale — it was the wrong **shape**. The post-cutoff absent fraction
is 85.3% for SARS-CoV-2 and 99.0% for Zika, so there is no single number for a
module-level constant to hold.

The biology statements are different in kind — published virology that does
not change when the release is rebuilt — so they stay in prose, and they carry
references:

| claim | reference |
|:---|:---|
| `vero_ifn_deletion` | Osada N, Kohara A, Yamaji T, Hirayama N, Kasai F, Sekizuka T, Kuroda M, Hanada K (2014). The genome landscape of the African green monkey kidney-derived Vero cell line. DNA Research 21(6):673-683. doi:10.1093/dnares/dsu029. PMID:25267831. Reports a homozygous ~9 Mb deletion on chromosome 12 removing the type I interferon gene cluster, along with CDKN2A and CDKN2B. |
| `chloroquine_not_in_airway` | Hoffmann M, Mosbauer K, Hofmann-Winkler H, Kaul A, Kleine-Weber H, Kruger N, et al. (2020). Chloroquine does not inhibit infection of human lung cells with SARS-CoV-2. Nature 585(7826):588-590. doi:10.1038/s41586-020-2575-3. Chloroquine blocks entry in Vero and not in TMPRSS2-expressing Calu-3 lung cells. |
| `flavivirus_entry` | van der Schaar HM, Rust MJ, Chen C, van der Ende-Metselaar H, Wilschut J, Zhuang X, Smit JM (2008). Dissecting the cell entry pathway of dengue virus by single-particle tracking in living cells. PLoS Pathogens 4(12):e1000244. doi:10.1371/journal.ppat.1000244. Dengue enters by clathrin-mediated endocytosis and fuses from within acidified endosomes. |
| `flavivirus_ifn_restriction` | Lazear HM, Govero J, Smith AM, Platt DJ, Fernandez E, Miner JJ, Diamond MS (2016). A mouse model of Zika virus pathogenesis. Cell Host & Microbe 19(5):720-730. doi:10.1016/j.chom.2016.03.010. Zika produces lethal neurological disease in mice lacking the type I interferon receptor (Ifnar1-/-), which is why an interferon-deficient background is needed to model it. The claim supported is that restriction, not a wild-type comparison the abstract does not report. |

`cite()` raises on an unknown key, because a claim whose source has been
removed must break the render rather than quietly lose its citation.

## 5. The tools are fixed functions

A language model asked to generate traversal code against a
24-predicate schema fails **quietly**: a malformed
traversal returns an empty result, and an empty result is indistinguishable
from "no evidence exists". The model then reports the latter.

| tool | what it answers |
|:---|:---|
| `calibration_for` | The calibration measured ON this virus, or None. |
| `evidence` | Every edge supporting a specific claim, with provenance. |
| `explain` | Mechanistic paths from a compound to a virus, with provenance. |
| `family_of` | The viral family of a taxon, from the graph's own OrganismTaxon. |
| `neighbours` | Adjacent nodes, capped. |
| `pathways` | Which pathways this compound's targets sit in, and whether the |
| `profile` | What the graph holds about one entity. |
| `resolve` | Free text to node identifiers. |
| `routes_for` | Which feature channels this release holds for one virus. |
| `temporal_for` | The PROSPECTIVE calibration measured on this virus, or None. |
| `triage` | Rank compounds you already have. |
| `virus_pathways` | Which pathways this virus depends on, and which of them are |

Every tool returns a structured record **plus a verbalisation**, and the
verbalisation is rendered here rather than by the model, because leaving it to
the model means it invents phrasing for facts it half understands. The facts
are finished English written in code from the edge tuples — direction markers
were tried first and taught the model nothing about arrows: it read
`<- BELONGS_TO` and wrote "Zika belongs to the gene". *A text marker is not an
enforcement mechanism*, and that took two predicates to establish.

Bounded by construction: 50 neighbours,
20 paths, 100 subgraph edges. A route
reaching fewer than 25 edges is called thin and says so.
Selectivity is reported against SI 10.0.

## 6. Verification withholds, it does not annotate

Every sentence the model writes is checked against the retrieved facts. If a
sentence is not backed by one, the paragraph is **withheld** —
`ANSWER WITHHELD — FAILED VERIFICATION` — and the facts stay. Withholding is
right, but withholding on the first attempt made it the normal outcome, so the
model gets **one** retry with the omitted caveat label named. Nothing else is
retried: an invented citation or a fabricated number is not an oversight to
nudge, and re-rolling until a check passes is how a verifier becomes a
sampler.

The machine surface withholds on the same condition as the terminal one. It
did not: `render` withheld the answer while the `--json` branch printed the
payload whole, so `jq -r .answer` returned text that had failed verification
with the failure in a sibling key.

## 7. Follow-ups: three conditions, and the third is not optional

Every question is resolved on its own, so a follow-up resolves nothing: *why
that one?*, *and for DENV-3?*, *tell me more* each used to return "No entity
in this question matches a node in the graph". 53 reference
phrases are recognised, and a carry happens only when **all three** hold:

1. the question contains an explicit reference back — a pronoun, a
   demonstrative, an ordinal, or a continuation like *what about*;
2. it names nothing of its own for the slot being filled, so
   *and what about DENV-3?* is never answered about DENV-2;
3. the reply **says what it took**, as a `CARRIED FORWARD` caveat — which is
   the channel verification enforces, so the synthesised prose cannot drop it
   either.

Only node ids travel: a virus, the compounds the previous **answer** named,
and the previous intent. **No free text**, deliberately — a context carrying a
sentence would be a way to put words into the brief, and the brief is the one
thing the model is given. Everything is revalidated on arrival, because the
context is echoed back by a web page.

An ambiguous reference **asks**: *why that one?* after a list of twelve prints
the twelve and waits. Picking the first would answer about a compound nobody
asked about, under a notice saying it had been carried, which is worse than
refusing.

The page keeps the previous turn, not the server. The server is one process on
a loopback port and may be read by more than one person, so "the last
question" held there would make one person's follow-up resolve against
another's topic.

## 8. The frontend

`web/chat.html` — one HTML file, stdlib `http.server` behind it, no npm
install step to forget. It binds to `127.0.0.1` by default because the
payloads name compounds, targets and unpublished AUCs; `BIND=0.0.0.0` is an
explicit choice. To reach it from another machine, forward the port rather
than binding wider:

```
ssh -N -L 8765:127.0.0.1:8765 <host>
```

The reply is rendered in the same order as the terminal: facts, then caveats,
then the written answer. That order is deliberate — **a warning printed under
a fluent paragraph is a warning nobody reads**.
