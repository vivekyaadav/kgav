# 11 — How to use it

This page is for everybody on the team, including people who do not usually
work in a terminal. You type a question in plain English; the system looks it
up in the graph; it prints the facts it found, each with the paper it came
from, and the warnings that say how far to trust them.

- It runs **entirely on the workstation**. Nothing is sent to the internet,
  including to the language model, which is local.
- It listens on **this machine only**. Nobody else on the network can open it.
- It will **refuse** questions it cannot answer honestly, and say why. A
  refusal is the system working.

[10-agent.md](10-agent.md) is the same system described for whoever maintains
it. This page is how to drive it.

---

## 1. What is in the graph

| node class | nodes | id prefixes |
|:---|---:|:---|
| SmallMolecule | 25,243 | INCHIKEY |
| Protein | 20,888 | UniProtKB |
| Gene | 19,190 | NCBIGene, KGAV |
| Pathway | 13,236 | REACT, GO |
| TargetClass | 66 | SCOP, CATH, PFAM, INTERPRO, SSF, KGAV |
| OrganismTaxon | 16 | NCBITaxon |
| Disease | 10 | MONDO |

Those counts are read out of release `v0.1` as this page was
generated; the predecessor of this page typed them in and was wrong by 140
compounds for a fortnight. The live figures are also on the web page's
release chip — hover it.

The whole picture on one printable page: `figures/kgav_metagraph_a4.pdf`.
How the pieces fit together: `figures/kgav_architecture_a4.pdf`.

**`HAS_ANTIVIRAL_ACTIVITY_AGAINST`, `IN_TRIAL_FOR`, `MEASURED_INACTIVE_AGAINST`, `TREATS` are held out.** They are the laboratory results the
system is graded against, and they are deliberately kept out of the
reasoning. When an answer quotes one, it is labelled **ANSWER SHEET** —
meaning *this is real data, but it is the answer, not a prediction*.

---

## 2. Starting it (step by step)

You only need to copy and paste.

**Step 1 — open a terminal.** On a Mac: `Cmd + Space`, type `Terminal`,
press Enter.

**Step 2 — connect to the workstation.** Replace `YOUR-USERNAME` with the
account you were given:

```
ssh YOUR-USERNAME@10.10.48.150
```

Nothing appears while you type the password. That is normal.

**Step 3 — go to the project and switch on its Python:**

```
cd /home/vikas/projects/kgav
source /home/vikas/anaconda3/bin/activate kgav
```

(If you are logged in as `vikas`, `conda activate kgav` does the same.)

**Step 4 — ask something:**

```
python scripts/kgav_ask.py --no-model "explain how nirmatrelvir works against SARS-CoV-2"
```

To ask another, press the up-arrow, edit the text inside the quotes, Enter.

**Step 5 — when you are finished,** type `exit`.

---

## 3. Three ways to ask

| how | command | what you get |
|:---|:---|:---|
| **the web page** (easiest) | `make agent-ui`, then <http://127.0.0.1:8765> | a chat window; the "evidence only" tickbox is the second row |
| **facts only** | `python scripts/kgav_ask.py --no-model "…"` | the facts and the caveats, identical every time |
| **facts + written answer** | `python scripts/kgav_ask.py "…"` | also a paragraph from a local model, every sentence of which is checked against a retrieved fact before you see it |

All three run the same agent and apply the same guards. The web page is
faster per question for a structural reason rather than a clever one: the
server loads the graph once and keeps it, while each command-line run loads
it again.

### The web page

```
make agent-ui          # defaults: port 8765, bound to 127.0.0.1
```

Then open <http://127.0.0.1:8765> in a browser **on the workstation**.

- **Nobody else can open it.** It binds loopback. Somebody would have to
  pass `BIND=0.0.0.0` on purpose, and then anyone on the network can read
  whatever you ask it.
- **To use it from your laptop, forward the port — do not bind wider.**
  With the workstation reachable as `z8` in your `~/.ssh/config`:

  ```
  ssh -N -L 8765:127.0.0.1:8765 z8
  ```

  Leave that running and open <http://localhost:8765> on the laptop. The
  server still listens only on the workstation's loopback; the tunnel
  carries it, so nothing is exposed to either network. `-N` means "no
  shell, just the forward". Without an ssh-config entry, spelled out:

  ```
  ssh -N -J <user>@<jump-host> -p <port> <user>@localhost -L 8765:127.0.0.1:8765
  ```

- **The order of the reply is deliberate**: facts, then caveats, then the
  written paragraph. A warning printed underneath a fluent paragraph is a
  warning nobody reads.
- **The two chips at the top** show which release is loaded and whether the
  local model is running. Hover them. A grey model chip means you still get
  facts, without the paragraph.
- If you edit the page itself (`web/chat.html`), just reload the browser.
  There is no build step and no restart.

---

## 4. Follow-up questions

Ask *"why that one?"*, *"and for DENV-3?"*, *"the third"*, *"tell me more"*
and the question is resolved against the turn before it. Three things have
to be true before anything is inherited, and the third is the point:

1. the question has to **refer back** — a pronoun, a demonstrative, an
   ordinal, or a continuation like "what about". There are
   53 recognised phrasings;
2. it has to name **nothing of its own** for the slot being filled, so
   "what is known about DENV-3?" is never answered about DENV-2 merely
   because the previous question was;
3. the reply **says what it took**, as a **CARRIED FORWARD** caveat.

Nothing else travels. The previous turn contributes a virus, the compounds
its answer named, and its intent — all as node identifiers, re-checked
against the graph on arrival. If a reference is ambiguous you are **asked**
rather than guessed at: "why that one?" after a list of twelve prints the
twelve and waits.

**The page remembers the last turn; the server does not.** The server is one
process on a loopback port and may be read by more than one person, so
holding "the last question" there would make your follow-up resolve against
somebody else's topic. **Clear** forgets it.

On the command line each run is its own process, so a follow-up is given the
previous reply's `context` explicitly:

```
python scripts/kgav_ask.py --no-model --json "what is known about DENV-2?" | jq -c .context
python scripts/kgav_ask.py --no-model --context '{"virus":"NCBITaxon:11060","intent":"profile"}' "tell me more"
```

---

## 5. Naming a virus

The system recognises **16** viruses and
78 spellings of their names. **If it does not recognise the
virus it answers about SARS-CoV-2 instead** — it warns you, under
**VIRUS ASSUMED**, so always check the caveats for that line.

| virus | write it as (any of these) | deliberately not accepted |
|:---|:---|:---|
| SARS-CoV-2 | `sarscov2`, `sars-cov2`, `sars cov 2`, `sars-cov-2`, `severe acute respiratory syndrome coronavirus 2` | — |
| SARS-CoV | `sars-cov`, `sars-cov-1`, `sars coronavirus`, `severe acute respiratory syndrome coronavirus` | — |
| MERS-CoV | `mers-cov`, `mers coronavirus`, `middle east respiratory syndrome coronavirus` | — |
| HCoV-229E | `hcov 229e`, `hcov-229e`, `coronavirus 229e`, `human coronavirus 229e` | — |
| HCoV-NL63 | `hcov nl63`, `hcov-nl63`, `coronavirus nl63`, `human coronavirus nl63` | — |
| HCoV-OC43 | `hcov oc43`, `hcov-oc43`, `coronavirus oc43`, `human coronavirus oc43` | — |
| HCoV-HKU1 | `hcov hku1`, `hcov-hku1`, `coronavirus hku1`, `human coronavirus hku1` | — |
| Zika | `zika`, `zikv`, `zika virus` | — |
| DENV-2 | `denv2`, `denv 2`, `denv-2`, `dengue 2`, `dengue type 2`, `dengue virus 2`, `dengue virus-2`, `dengue serotype 2`, `dengue virus type 2` | — |
| DENV-1 | `denv1`, `denv 1`, `denv-1`, `dengue 1`, `dengue type 1`, `dengue virus 1`, `dengue virus-1`, `dengue serotype 1`, `dengue virus type 1` | — |
| DENV-3 | `denv3`, `denv 3`, `denv-3`, `dengue 3`, `dengue type 3`, `dengue virus 3`, `dengue virus-3`, `dengue serotype 3`, `dengue virus type 3` | — |
| DENV-4 | `denv4`, `denv 4`, `denv-4`, `dengue 4`, `dengue type 4`, `dengue virus 4`, `dengue virus-4`, `dengue serotype 4`, `dengue virus type 4` | — |
| West Nile | `west nile`, `west nile virus` | — |
| JEV | `japanese encephalitis`, `japanese encephalitis virus` | `JEV` |
| Yellow fever | `yellow fever`, `yellow fever virus` | — |
| TBEV | `tbev`, `tick borne encephalitis`, `tick-borne encephalitis`, `tick borne encephalitis virus`, `tick-borne encephalitis virus` | — |

**1 short name is refused on purpose**, in the last column. A name shorter than 4 characters, or one that is a substring of ordinary English, matches text that is not a virus — and a wrong match is worse than no match, because it answers confidently about the wrong organism.

Drugs are found by their common names — "nirmatrelvir", "remdesivir",
"molnupiravir", "baricitinib".

---

## 6. Questions it answers

| you want | start the question with |
|:---|:---|
| how a drug might act on a virus | `explain how DRUG works against VIRUS` |
| the raw evidence behind a claim | `what evidence links DRUG to VIRUS` |
| everything about one drug | `profile of DRUG` |
| the pathways involved | `which pathways does DRUG touch` |
| compare drugs you already have | `rank DRUG1, DRUG2 and DRUG3 for VIRUS` |
| what the system itself can do | `what can you do for Zika?` |
| how a result was arrived at | `how was the prospective test done?` |

Internally these become `explain`, `triage`, `evidence`, `profile`, `pathway`, `prospect`, `capability`, `method`, `unknown`. Try:

```
python scripts/kgav_ask.py --no-model "explain how nirmatrelvir works against SARS-CoV-2"
python scripts/kgav_ask.py --no-model "what evidence links remdesivir to SARS-CoV-2"
python scripts/kgav_ask.py --no-model "profile of molnupiravir"
python scripts/kgav_ask.py --no-model "which pathways does remdesivir touch"
python scripts/kgav_ask.py --no-model "rank remdesivir, favipiravir and ribavirin for DENV-2"
python scripts/kgav_ask.py --no-model "which compounds should be considered for Yellow fever"
```

The last one is the system **refusing on purpose**: Yellow fever has no route
in this release, and an empty ranking reads as "nothing looks promising"
when the truth is "nothing could be looked at".

---

## 7. Shortlists, and which viruses are allowed one

A shortlist is the only output here that proposes something you have not
already named, so it is the most tightly gated. It is produced by
`scripts/repurposing_candidates.py` and offered for a virus only where the
cross-sectional control clears **and** the prospective control does not
contradict it — **8 of 16 viruses**:
DENV-1, DENV-2, DENV-3, HCoV-229E, JEV, SARS-CoV-2, TBEV, Yellow fever.

**2 viruses clear cross-sectionally and anti-rank prospectively — West Nile, Zika — and are refused by name.** A cross-sectional verdict alone does not license a shortlist, because a shortlist is a prospective act.

**What a shortlist row claims, exactly.** That the compound *resembles* a
named known active for that virus, at a stated Tanimoto, and where that
Tanimoto sits among the compounds that virus has already **screened and
found inactive**. It is not an expected hit rate; see
[06-methodology.md](06-methodology.md) §11.

Compounds that are not therapies are dropped, and compounds whose class is a
real confound are kept and labelled:

| virus | verdict | known actives | candidates | resemblance floor | failures scored | rows the floor cost | rows shown |
|:---|:---|---:|---:|:---|---:|---:|---:|
| DENV-1 | clears | 195 | 24,883 | 0.5246 | 123 | 5 | 40 |
| DENV-2 | clears | 995 | 21,163 | 0.4267 | 3,043 | 0 | 45 |
| DENV-3 | clears | 213 | 24,893 | 0.5156 | 95 | 15 | 30 |
| DENV-4 | below | 153 | 24,939 | 0.4143 | 109 | 0 | 0 |
| HCoV-229E | clears | 119 | 24,994 | 0.6664 | 88 | 20 | 25 |
| HCoV-OC43 | weak | 162 | 24,969 | 0.7264 | 70 | 25 | 20 |
| JEV | clears | 19 | 25,155 | 0.2222 | 27 | 0 | 45 |
| MERS-CoV | weak | 199 | 23,044 | 0.2115 | 1,989 | 0 | 45 |
| SARS-CoV | weak | 232 | 23,386 | 0.2316 | 1,610 | 0 | 45 |
| SARS-CoV-2 | clears | 1,768 | 12,420 | 0.2667 | 11,048 | 0 | 45 |
| TBEV | clears | 72 | 25,064 | 0.5968 | 65 | 25 | 20 |
| West Nile | clears | 261 | 24,629 | 0.5634 | 311 | 12 | 33 |
| Yellow fever | clears | 135 | 24,184 | 0.3008 | 882 | 0 | 45 |
| Zika | clears | 748 | 21,876 | 0.5195 | 2,577 | 0 | 45 |

---

## 8. Reading the output

The answer arrives in labelled blocks.

**EVIDENCE** — the facts. A line looks like:

```
nirmatrelvir inhibits the viral protein nsp5 (IC50 75 nM, potent)
  [chembl, PMID:34798775 …, 2021; E83e08310cff11af7]
```

- `IC50 75 nM` — how much drug is needed. **Smaller is stronger.** The cut
  here is 10,000 nM, which is 10 µM: at or below it a measured
  potency counts as **active**, above it as **measured-inactive**. 75 nM is
  more than a hundred times below the cut.
- `chembl` — which database the fact came from.
- `PMID:34798775` — the paper. Search it at pubmed.ncbi.nlm.nih.gov.
- `2021` — when the fact was first reported.
- `E83e…` — the edge's identifier in the graph, for the technical team.

A **Route** is a chain of facts: *drug → blocks a viral protein → that
protein belongs to the virus*. [06-methodology.md](06-methodology.md) §5
explains what each chain does and does not establish.

**CAVEATS** — warnings. They are part of the answer, not decoration. There
are 20:

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

**REFUSED** — it declined, and says why.

**ANSWER WITHHELD — FAILED VERIFICATION** — the model wrote a paragraph that
contained a sentence no retrieved fact supported, so the paragraph was
dropped. **The facts above it are still correct.** This happens fairly often
and it is the check working, not a crash.

---

## 9. What it can and cannot do

**Can:**

- show, with sources, how a known drug is linked to a virus;
- gather everything the graph holds about one compound;
- rank compounds *you* name, for a virus, and say what the ranking rests on;
- propose a shortlist for the 8 viruses where a control
  licenses one, as a resemblance claim with the resembled compound named;
- say what it cannot do for a given virus, and why.

**Cannot:**

- rank **novel** chemical matter. Compounds enter screening because they are
  novel, and novelty is what a prior-evidence graph cannot represent;
- answer for a virus with no route at all — it refuses rather than returning
  an empty list;
- tell you a drug works in patients. Everything here is a laboratory
  measurement, and [09-limits.md](09-limits.md) is the full list of what
  that excludes.

---

## 10. When something goes wrong

| you see | do this |
|:---|:---|
| `conda: command not found` / `python: command not found` | run `source /home/vikas/anaconda3/bin/activate kgav` again |
| `No such file or directory` | wrong folder: `cd /home/vikas/projects/kgav` |
| a note about Ollama not being reachable | use `--no-model`, or ask for `sudo systemctl start ollama` |
| **VIRUS ASSUMED** when you did name a virus | use a spelling from §5 |
| it sits for more than two minutes | `Ctrl + C`, retry with `--no-model` |
| `note: … is absent` | a message for the technical team; ignore it |
| the port is already in use | somebody is already running the server; just open the page |

---

## 11. For whoever maintains it

- Entry point `scripts/kgav_ask.py`; logic in `src/kgav/agent/`. `--json`
  gives machine-readable output, `--show-prompt` prints what is sent to the
  model, `--context` supplies the previous turn.
- Server `scripts/agent_server.py`; `make agent-ui` passes `PORT` and
  `BIND`, which default to `8765` and `127.0.0.1`.
- Model: Ollama on `localhost:11434`, ``qwen2.5:7b-instruct-q5_K_M``, temperature 0.
  Change it with `--model` / `--host`. The model never chooses a refusal —
  every guard is keyword-deterministic and runs before the model is called.
- Data: release `v0.1` (78,649 nodes, 396,223 edges),
  plus the artifacts in `data/results/` that supply the calibration caveats.
- The tools the orchestrator can call, and what each answers:

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

- Retrieval is capped: 20 paths, 50
  neighbours, 100 subgraph edges. A cap reached is
  reported, never silently truncated.
