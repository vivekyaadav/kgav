# 13 — The tests

3,058 test functions in 138 files,
57,099 lines — against 21,673 lines of library and
31,730 lines of scripts. 26 are marked `slow` and
deselected by `make test`.

The ratio is the point. Almost none of these tests check that a function
returns the right value for a toy input. They check **properties of the
measurement**: that a figure in a document is still the one the run produced,
that a check can fail, that a refusal is reachable, that two implementations
of one calculation do not exist.

## The kinds of property

### 1. A published figure is re-derived from its artifact

Not *is the number plausible* but *is this number the one the run wrote*. The
pattern is a test that reads the results file, recomputes the figure, and
asserts the document contains exactly it. It exists because a hand-typed AUC
sat in a banner beside a cell that had since moved by two thousandths, and
because "none below chance" sat in the same banner for a family with two cells
below chance.

Stronger than presence: several of these assert that **every** four-decimal
figure in a section is one the run produced, because checking only presence
passes while a second occurrence says something else — one figure appeared
twice in a single amendment, and perturbing one occurrence left the other to
satisfy the test. `tests/test_doc2.py` applies the same rule to this folder.

### 2. A document's prose counts are checked, not just its numbers

"3 of 6 survive" is the sentence a reader takes away, and a mutation of it
passed every figure check in the same file. So the counts in prose are
asserted against the derivation too, and every named cell has to be named.

### 3. A check is proven able to fail

`make gate` is the clearest instance: the validator is run against a release
carrying four planted violations and the **codes and counts** are asserted. It
used to be *"if the validator fails, OK"*, which accepted an `ImportError` as
proof. `tests/test_negative_control.py`, `tests/test_auc_delta.py` and
`tests/test_compare_results.py` do the same for three other checks.

### 4. One calculation has one implementation

`tests/test_no_duplicate_logic.py`. Unit conversion existed three times and
the three did not agree — two of them differed in the last bit for 13.75% of
values, and the units they supported had already diverged. The threshold logic
that decides a label exists once, in `labels.classify`, and the ingest that
routes an edge by polarity calls **that** function, so the features and the
labels cannot disagree about what "acts on" means.

### 5. A registration is what the code actually does

`tests/test_prereg_is_obeyed.py` and the seven `test_prereg_*` files. A
registration is binding only if something reads it: these assert that the
operands a registration states are real, that the command it records is the
command that runs, and that a verdict quoted in it is derivable from the file
it registered.

### 6. A refusal is reachable

A unit test asks the question its author thought to ask.
`make tier-check` prints what the agent says for one virus from **each** tier,
and the first run of it found that the refusal written for Yellow fever was
never reached — which every unit test had missed.

### 7. The Makefile's own targets resolve

`tests/test_makefile_targets_resolve.py`. Every variable a recipe expands must
carry a default, be wrapped in `$(or)`, or be guarded by a `test -n` in its
own recipe. `make release-validate` expanded to `data/releases/` and exited 2
on **every invocation since the target was written** — a target that cannot
run is the same defect class as a check that cannot fail.

### 8. No script crashes while printing its own help

`tests/test_every_script_help_runs.py` greps for a traceback rather than
checking the exit code, because four scripts exit non-zero on `--help` **on
purpose** and each of those is a deliberate message. A list of four names
would go stale the next time a script was added.

`hard_negatives.py --help` once raised `TypeError: %o format: an integer is
required, not dict`, because a help string contained "94.9% of compounds" and
argparse %-formats help text. The flag it documents defines the headline
protocol, and its help had never once printed.

## The largest files, and what they hold

| test file | test functions | what it holds |
|:---|---:|:---|
| `tests/test_baselines.py` | 112 | Gate for the baselines: held-out edges are never traversed, DWPC |
| `tests/test_landscape.py` | 100 | The agent may not imply its ranking is predictive. |
| `tests/test_folds.py` | 91 | Gate for phase 2: the two-mode comparison, and the rule that a shared fold |
| `tests/test_gtop.py` | 81 | Gate for the Guide to PHARMACOLOGY layer. |
| `tests/test_agent_guards.py` | 78 | Gate for the agent's guardrails: the refusal is deterministic, and every |
| `tests/test_pathway_landscape.py` | 69 | Gate for the pathway landscape. |
| `tests/test_assemble.py` | 68 | Day 8 gate: cross-layer problems are caught, not merged away. |
| `tests/test_agent_pathways.py` | 67 | A5: which pathways a compound's targets sit in, and what that is not. |
| `tests/test_ncats.py` | 65 | The panel layer's two contracts. |
| `tests/test_vhppi.py` | 65 | Day 5 gate: the bridge layer joins, and alias resolution is what makes it. |
| `tests/test_agent_tools.py` | 62 | Gate for the tool layer: bounded output, provenance in every sentence, and |
| `tests/test_agent_orchestrator.py` | 61 | Gate for the orchestrator: the pipeline is deterministic, the refusal comes |
| `tests/test_agent_verify.py` | 61 | Gate for verification: hallucination is detected deterministically, not by |
| `tests/test_orcs.py` | 61 | Day 6 gate: polarity is never guessed, and a flipped sign is rejected. |
| `tests/test_phase4_survey.py` | 59 | Phase 4's feasibility survey: the classifiers, the rejections, the verdict. |
| `tests/test_chembl.py` | 55 | Day 7 gate: nsp resolution from assay text, censored values preserved, |
| `tests/test_host_targets.py` | 52 | Gate for the drug -> human protein layer: direction is never invented. |
| `tests/test_kg_audit_report.py` | 50 | doc2/kg_audit_report.txt, checked against the code it describes. |
| `tests/test_labels.py` | 50 | Gate for measured negatives: a lower bound is only inactivity if it is high |
| `tests/test_fetch_vhppi.py` | 47 | Gate for the virus-host PPI fetcher. |

## Mutation testing, as a working habit

A test that cannot distinguish the fix from the bug is the defect this
repository finds most often, and the counter-measure is to break the code on
purpose and watch the test fail. Instances recorded in the commit log include:
a test that grepped a source file and passed with the logic replaced by
`if False: pass`; a guard clause that **survived** every mutation and turned
out to be dead code reading as load-bearing; and a condition-2 test that
stayed green when the guard it was written for was deleted outright, because
its example could not reach that branch.

## What the suite does not do

It does not rebuild the graph — `make rebuild` is ~25 minutes and reads 30 GB
of ChEMBL — and it does not run the long diagnostics. Those are
`make test-all`, the `slow` mark, and the gates in
[12-reproducing.md](12-reproducing.md). What the suite guarantees is that
**every figure in every document still matches the artifact it came from**,
and that every check in the repository can still fail.
