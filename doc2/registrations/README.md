# The pre-registrations

Each document here was signed **before** the result it governs existed, and
each is amended by **appending**. A registration edited to match the result
it registered is not a registration, which is why none of them is generated
and none is rewritten — including for the move described below.

`scripts/write_doc2.py` builds the index in
[../15-appendix.md](../15-appendix.md) from this directory, so a
registration added here appears there without anyone editing a table.
`tests/test_prereg_is_obeyed.py` and the per-registration verdict tests read
these files and assert the running code obeys them: the permutation count,
the alpha, the family size, the exclusion rule.

## Paths inside these documents, as of 2026-11-21

On 2026-11-21 `docs/` was retired and its contents moved into `doc2/`. These
documents were **not** edited, so a path written inside one still points at
where that file lived when the registration was signed. The mapping:

| written in a registration | where it is now |
|:---|:---|
| `docs/prereg_*.md` | `doc2/registrations/prereg_*.md` |
| `docs/evaluation.md` | superseded by [../06-methodology.md](../06-methodology.md), [../07-protocols.md](../07-protocols.md) and [../08-results.md](../08-results.md); its Reproducing block is in [../12-reproducing.md](../12-reproducing.md) |
| `docs/results.md` | superseded by [../08-results.md](../08-results.md) |
| `docs/scope_and_limits.md` | superseded by [../09-limits.md](../09-limits.md) |
| `docs/agent_guide.md` | superseded by [../11-user-guide.md](../11-user-guide.md) |
| `docs/v1_release_plan.md` | [../plan.md](../plan.md) |
| `docs/kg_audit_report.txt`, `docs/kg_release_census.txt`, `docs/scorer_aucs.txt` | `doc2/` |
| `docs/figures/` | [../figures/](../figures) |
| `docs/audit_2026-11-03_fixes.md`, `docs/data_audit_2026-11-07.md` | `doc2/records/` |

A redirect table is the right repair here and an edit is not: the value of
these files is that their bytes have not changed since they were signed.
