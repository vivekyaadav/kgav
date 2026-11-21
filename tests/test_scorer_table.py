"""The generated AUC table says what the results files say.

A table of AUCs is the one artifact in this repository that a reader
will quote without opening anything else, so the ways it can be wrong
matter more than usual: a stale figure, a refusal formatted as a low
AUC, a cell quietly dropped because its channel key moved, or a
significance mark carried over from the wrong alpha. Each of those has
a test here, and each is checked against data/results rather than
against the file's own text.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import write_scorer_table as W  # noqa: E402

TABLE = ROOT / "doc2" / "scorer_aucs.txt"


@pytest.fixture(scope="module")
def table() -> str:
    assert TABLE.exists(), "doc2/scorer_aucs.txt is missing; run `make scorer-table`"
    return TABLE.read_text()


def _docs() -> dict[str, dict]:
    return {name: W._load(name) for _, name, _ in W.COLUMNS}


def test_the_file_on_disk_is_what_the_results_files_produce(table):
    """The whole point of generating it. `make scorer-table-check` is this
    assertion as a gate; this is it on every commit."""
    assert table == W.build(), (
        "doc2/scorer_aucs.txt is stale. Run `make scorer-table`.")


def test_every_evaluable_cell_in_every_results_file_appears(table):
    """A cell is dropped by a channel key moving, which looks like nothing.

    Builds the expected set from the results files and requires each
    AUC to be present in the table at four decimals, so a renamed
    channel removes a row and fails here rather than shrinking the
    table silently.
    """
    docs = _docs()
    missing = []
    for pop, _, _, _, _ in W.POPULATIONS:
        for _, fname, channel in W.COLUMNS:
            block = docs[fname].get(pop) or {}
            for virus, cell in block.items():
                if not isinstance(cell, dict):
                    continue
                m = W._cell(docs[fname], pop, virus, channel)
                if m is None:
                    continue
                if f"{m['auc']:.4f}" not in table:
                    missing.append((pop, virus, channel, m["auc"]))
    assert not missing, f"AUCs absent from the table: {missing}"


def test_a_refused_cell_is_printed_as_refused_and_not_as_an_auc(table):
    """The failure mode that would misread as evidence.

    The publication-disjoint protocol refuses all three plate cells, as
    its registration requires. A refusal carries no "auc" key, so a
    formatter that reached for it would print 0.0000 -- a number a
    reader would take as "the scorer failed" rather than "the protocol
    declined".
    """
    docs = _docs()
    pub = docs["structure_scorer.json"]
    refused = [
        (pop, virus)
        for pop, _, _, _, _ in W.POPULATIONS
        for virus, cell in (pub.get(pop) or {}).items()
        if isinstance(cell, dict)
        and W._cell(pub, pop, virus, "structure_max") is None
    ]
    assert refused, (
        "no cell refuses under the publication control, so this test is "
        "vacuous -- the three plate cells are registered to refuse")
    assert "refused" in table
    assert "0.0000" not in table, (
        "a 0.0000 in the table is a refusal formatted as an AUC")


def test_the_significance_mark_is_the_alpha_of_the_block_it_sits_in(table):
    """Two families, two z values, and the plate cells belong to both.

    MERS-CoV's entry cell is the one that distinguishes them: it clears
    at the plate alpha (z = 2.3940) and does not at the ChEMBL one
    (z = 2.9137). A table printing one mark would have to choose, and
    whichever it chose would be wrong for some reader.
    """
    docs = _docs()
    sym = docs["structure_scorer_symmetric.json"]
    m = W._cell(sym, "secondary_ncats_entry", "NCBITaxon:1335626",
                "structure_max_k10")
    assert m is not None
    assert m["clears_corrected_alpha"] is True
    assert m["clears_alternate_alpha"] is False, (
        "MERS-CoV's entry cell no longer separates the two alphas, so this "
        "test no longer proves the table prints both")

    block = table.split("NCATS PLATE, PSEUDOTYPED_ENTRY")[1]
    primary, alternate = block.split("AND AT THE 14-CELL ALPHA")
    mers_primary = [ln for ln in primary.splitlines() if "MERS-CoV" in ln][0]
    mers_alt = [ln for ln in alternate.splitlines() if "MERS-CoV" in ln][0]
    assert "0.5873*" in mers_primary, mers_primary
    assert "0.5873" in mers_alt and "0.5873*" not in mers_alt, mers_alt


def test_the_table_quotes_the_release_each_file_was_measured_on(table):
    """A table assembled from two releases is the failure the stamp guards.

    Not a comparison that passes silently: the hashes are printed, so a
    mixture is readable in the artifact itself and not only in a test.
    """
    docs = _docs()
    shas = {
        ((d.get("provenance") or {}).get("release") or {}).get("nodes_sha")
        for d in docs.values()
    }
    for sha in shas:
        assert sha and sha in table, f"{sha} not stamped in the table"


def test_degree_is_one_baseline_across_the_three_runs():
    """The column reads one file, so the agreement has to be checked.

    If a future change made a study's degree baseline depend on its own
    protocol, a single-source column would show a number that still
    looked like a baseline. The script prints a disagreement block; this
    asserts there is nothing to print.
    """
    docs = _docs()
    for pop, _, _, _, _ in W.POPULATIONS:
        viruses = sorted({v for d in docs.values() for v in (d.get(pop) or {})
                          if isinstance((d.get(pop) or {}).get(v), dict)})
        out = W._degree_agreement(docs, pop, viruses, {})
        assert not out, "\n".join(out)


def test_check_mode_fails_on_a_perturbed_file(tmp_path):
    """Proof the gate can fail, which is the only way to know it is a gate."""
    original = TABLE.read_text()
    try:
        TABLE.write_text(original.replace("0.8304", "0.9999"))
        r = subprocess.run(
            [sys.executable, "scripts/write_scorer_table.py", "--check"],
            cwd=ROOT, capture_output=True, text=True)
        assert r.returncode == 1, r.stdout + r.stderr
        assert "STALE" in r.stdout
    finally:
        TABLE.write_text(original)
