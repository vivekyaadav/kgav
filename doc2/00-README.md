# kgav, described from the code

An antiviral drug-repurposing knowledge graph, its evaluation, and an agent
that answers questions about it without overstating what it measured.

**Release `v0.1`** — 78,649 nodes, 396,223 edges,
schema `0.16.0`, built at commit `e2bc1305fb19`
(nodes `dc6d76a2b3cfb2a1`, edges `873623c9a9df2049`).

---

## This folder is generated

Every figure below is computed from the repository at the moment of writing:
from `schema/kg_schema.yaml`, `config/viruses.yaml`,
`data/releases/v0.1/MANIFEST.json`, the 36 artifacts in
`data/results/`, and the constants imported from `src/kgav/`. The prose lives
in `doc2/_templates/*.md.in`; `scripts/write_doc2.py` renders it.

```
make doc2          # regenerate
make doc2-check    # fail if any file no longer matches the code
```

Three rules the generator enforces, each of which this repository learned by
shipping its opposite:

- **a placeholder with no value is an error**, not an empty string — a figure
  that renders blank reads like prose that was always meant to read that way;
- **a derived value no template uses is reported** — a derivation nothing
  renders is dead code that looks like coverage;
- **`--check` diffs and names the first line that moved**, so a stale document
  fails a gate instead of being read.

**`docs/prereg_*.md` is not reproduced here.** Those are signed
pre-registrations, amended by appending and never rewritten. doc2 describes
them, quotes their verdicts out of the results files they registered, and
links to them. Rewriting a signed registration from scratch would destroy the
only property that makes it one.

---

## The documents

| file | what it covers |
|:---|:---|
| [01-overview.md](01-overview.md) | what the project claims, what it refuses to claim, and the one-paragraph version of every result |
| [02-graph-model.md](02-graph-model.md) | the schema: 11 node classes, 24 predicates, qualifiers, the two held-out flags, 8 metapaths |
| [03-data-sources.md](03-data-sources.md) | every source, its version, what it contributed, and what it cost |
| [04-build-pipeline.md](04-build-pipeline.md) | 13 ingest layers in dependency order, assembly, and the three provenance stamps |
| [05-labels.md](05-labels.md) | what counts as active, what counts as measured-inactive, and where both thresholds come from |
| [06-protocols.md](06-protocols.md) | every evaluation protocol: baselines, the three structure-scorer controls, the temporal split, clinical labels, the pathway landscape |
| [07-results.md](07-results.md) | every measured number, per virus, with its interval and its verdict |
| [08-limits.md](08-limits.md) | what this graph cannot do, each limit with the figure behind it |
| [09-agent.md](09-agent.md) | the agent: guards, tools, verification, the frontend, follow-ups |
| [10-reproducing.md](10-reproducing.md) | 60 make targets, what each writes, and what checks it |
| [11-tests.md](11-tests.md) | 3,061 tests in 138 files, and the classes of property they enforce |
| [12-defect-ledger.md](12-defect-ledger.md) | the recurring defect patterns, named, with the instances that named them |

## The repository at a glance

| | |
|:---|---:|
| library modules (`src/kgav/`) | 36 files, 21,669 lines |
| scripts (`scripts/`) | 87 files, 30,762 lines |
| tests (`tests/`) | 138 files, 56,973 lines, 3,061 test functions |
| viruses built | 16 (Coronaviridae and Flaviviridae) |
| viruses surveyed, not built | 0 |
| registered results | 36 artifacts, 1 of them a recorded refusal |
