"""Gate for results provenance.

The failure this closes: metapath_comparison.json was generated the day
before the selectivity layer added 1,568 activity edges, and was still read
as a baseline afterwards. A positive count that had fallen 1,011 -> 889 read
as a rise from 736, because nothing in the file said which graph it described.
"""
import json
import re
from pathlib import Path

from kgav.provenance import (
    PROVENANCE_KEY,
    check_dir,
    check_file,
    release_fingerprint,
    write_results,
)

ROOT = Path(__file__).resolve().parents[1]


def _release(tmp_path, edges="a\n", nodes="n\n"):
    d = tmp_path / "rel"
    d.mkdir(exist_ok=True)
    (d / "nodes.jsonl").write_text(nodes)
    (d / "edges.jsonl").write_text(edges)
    return d


def test_a_results_file_records_the_release_it_came_from(tmp_path):
    rel = _release(tmp_path)
    out = tmp_path / "r" / "x.json"
    write_results(out, {"auc": 0.826}, rel, "0.9.0")
    doc = json.loads(out.read_text())
    assert doc["auc"] == 0.826
    prov = doc[PROVENANCE_KEY]
    assert prov["release"] == release_fingerprint(rel)
    assert prov["schema_version"] == "0.9.0"
    assert "generated_utc" in prov and "code" in prov


def test_a_matching_release_passes(tmp_path):
    rel = _release(tmp_path)
    out = tmp_path / "r" / "x.json"
    write_results(out, {"auc": 0.826}, rel)
    assert check_file(out, rel)[0] == "ok"


def test_a_release_rebuilt_in_place_is_caught(tmp_path):
    """THE ACTUAL BUG. The directory keeps its name; the graph is different.
    Fingerprinting content rather than the name is what catches it."""
    rel = _release(tmp_path)
    out = tmp_path / "r" / "x.json"
    write_results(out, {"auc": 0.826}, rel)
    (rel / "edges.jsonl").write_text("a\nb\n")      # selectivity layer lands
    status, msg = check_file(out, rel)
    assert status == "stale"
    assert "edges" in msg


def test_a_file_with_no_provenance_is_reported_not_assumed_fine(tmp_path):
    rel = _release(tmp_path)
    p = tmp_path / "old.json"
    p.write_text(json.dumps({"auc": 0.826}))
    status, msg = check_file(p, rel)
    assert status == "unstamped"
    assert "cannot tell which release" in msg


def test_check_dir_reports_every_file(tmp_path):
    rel = _release(tmp_path)
    d = tmp_path / "res"
    d.mkdir()
    write_results(d / "good.json", {}, rel)
    (d / "bare.json").write_text("{}")
    rows = {name: status for name, status, _ in check_dir(d, rel)}
    assert rows == {"good.json": "ok", "bare.json": "unstamped"}


def test_nodes_and_edges_are_both_fingerprinted(tmp_path):
    rel = _release(tmp_path)
    out = tmp_path / "r" / "x.json"
    write_results(out, {}, rel)
    (rel / "nodes.jsonl").write_text("n\nn2\n")     # edges untouched
    assert check_file(out, rel)[0] == "stale"


def test_an_unreadable_file_is_stale_not_silently_skipped(tmp_path):
    rel = _release(tmp_path)
    p = tmp_path / "broken.json"
    p.write_text("{not json")
    assert check_file(p, rel)[0] == "stale"


def test_a_missing_release_is_not_reported_as_current(tmp_path):
    """Two Nones compared equal, so a typo'd release path reported every
    results file as current -- a false "current" from the tool that exists to
    catch false "currents"."""
    from kgav.provenance import check_file, write_results

    release = tmp_path / "rel"
    release.mkdir()
    (release / "nodes.jsonl").write_text('{"id":"A","class":"Gene"}')
    (release / "edges.jsonl").write_text("")
    out = tmp_path / "r.json"
    write_results(out, {"n": 1}, release, "0.12.0")
    assert check_file(out, release)[0] == "ok"

    status, msg = check_file(out, tmp_path / "does-not-exist")
    assert status == "stale"
    assert "cannot verify" in msg


# ------------------- the audit must check the directory in actual use
def test_the_results_audit_checks_the_directory_the_scripts_write_to():
    """make check-results pointed at data/results-corrected, a snapshot from
    2026-09-17 that nothing writes to, while every script defaults to
    data/results. It reported 11 stale files in an unused directory and
    nothing about the one in use -- which held eight files with no provenance
    block at all, the exact case provenance.py was written for.

    Compares the Makefile target against the scripts' own defaults, because
    two places that must agree and are edited separately will drift.
    """
    mk = (ROOT / "Makefile").read_text()
    m = re.search(r"check-results:\s*\n\s*python scripts/check_results\.py "
                  r"\$\(or \$\(RESULTS\),([^)]+)\)", mk)
    assert m, "could not find the check-results target"
    audited = m.group(1).strip()

    written = set()
    for script in ("hard_negatives", "run_baselines", "temporal_split",
                   "compare_metapaths"):
        src = (ROOT / "scripts" / f"{script}.py").read_text()
        for d in re.findall(r'default=root / "(data/results[^"]*)"', src):
            if not d.endswith(".json"):      # a file, not an output directory
                written.add(d)

    assert written, "no script declares a results output directory"
    assert written == {audited}, (
        f"make check-results audits {audited!r} but the scripts write to "
        f"{sorted(written)}. An integrity check aimed at the wrong directory "
        f"reports a verdict about neither.")
# ------------------------------------------- a refusal is a recordable state
def test_a_refusal_replaces_the_numbers_it_declined_to_produce(tmp_path):
    """temporal_split.py aborts on a fatal leakage audit without writing,
    which left the PREVIOUS run's numbers on disk -- from a graph three
    revisions old -- with nothing saying the protocol now refuses to run.
    check_results then said "STALE, regenerate it", advice that cannot be
    followed because the refusal is why it cannot be regenerated.
    """
    from kgav.provenance import REFUSAL_KEY, write_refusal

    rel = _release(tmp_path)
    out = tmp_path / "temporal_2021_include.json"

    # a previous, successful run
    write_results(out, {"cutoff": 2021, "results": {"SARS-CoV-2": {"M1": 1}}},
                  rel)
    assert check_file(out, rel)[0] == "ok"
    assert "results" in json.loads(out.read_text())

    # the same protocol, now refusing
    write_refusal(out, {"cutoff": 2021, "results": {"SARS-CoV-2": {"M1": 1}}},
                  ["L5 label overlap: 5 positive compounds"], rel)
    doc = json.loads(out.read_text())

    # the dead numbers are gone, not kept alongside the refusal
    assert "results" not in doc, "a refusal that still carries numbers reads as numbers"
    assert doc[REFUSAL_KEY] == ["L5 label overlap: 5 positive compounds"]
    assert doc["cutoff"] == 2021, "the split's own parameters are still recorded"

    status, msg = check_file(out, rel)
    assert status == "refused", f"expected refused, got {status}"
    assert "L5 label overlap" in msg, "the reason must reach the operator"


def test_a_refusal_from_an_older_graph_still_reads_stale(tmp_path):
    """Staleness is checked before refusal: a refusal recorded against a graph
    that has since moved describes a refusal that may no longer happen, so it
    is not a current statement about this release.
    """
    from kgav.provenance import write_refusal

    rel = _release(tmp_path, edges="a\n")
    out = tmp_path / "temporal_2021_include.json"
    write_refusal(out, {"cutoff": 2021}, ["L5 label overlap: 5"], rel)
    assert check_file(out, rel)[0] == "refused"

    (rel / "edges.jsonl").write_text("a\nb\n")      # release moves
    status, msg = check_file(out, rel)
    assert status == "stale", f"expected stale, got {status}"
    assert "edges" in msg, "the stale message must name what moved"


# ------------------------------------------------- derived releases
# A RESULT'S PROVENANCE CHAIN IS ONLY AS LONG AS ITS WEAKEST LINK, and until
# these landed it stopped one step short. temporal_split.py cuts
# v0.1-train2021 from v0.1; hard_negatives.py reads the cut and writes a
# result; write_results stamps that result against v0.1. So a v0.1 rebuild
# left a training graph from the previous graph, and the number computed from
# it was stamped CURRENT -- a false "current" from the one tool that exists to
# prevent false "current" readings. Same class as
# test_a_missing_release_is_not_reported_as_current, one link further out.

def test_a_derived_release_records_the_release_it_was_cut_from(tmp_path):
    from kgav.provenance import DERIVED_STAMP_NAME, check_derived, write_derived

    rel = _release(tmp_path)
    train = tmp_path / "train"
    train.mkdir()
    write_derived(train, rel, "0.12.0", cutoff=2021, undated="include")

    assert (train / DERIVED_STAMP_NAME).exists()
    assert check_derived(train, rel)[0] == "ok"


def test_a_derived_release_goes_stale_when_its_parent_is_rebuilt(tmp_path):
    """THE ACTUAL BUG. v0.1 is rebuilt in place; v0.1-train2021 still sits on
    disk and still loads. Before the stamp, `exists()` was the whole check."""
    from kgav.provenance import check_derived, write_derived

    rel = _release(tmp_path)
    train = tmp_path / "train"
    train.mkdir()
    write_derived(train, rel, "0.12.0", cutoff=2021, undated="include")

    (rel / "edges.jsonl").write_text("a\nb\n")      # the Zika gene merge
    status, msg = check_derived(train, rel)
    assert status == "stale", f"expected stale, got {status}"
    assert "edges" in msg, "the message must name what moved"


def test_an_unstamped_derived_release_is_not_assumed_fine(tmp_path):
    """Every v0.1-train2021 built before the stamp existed lands here. It
    cannot be attributed to any graph, so it must not be used as if it could."""
    from kgav.provenance import check_derived

    rel = _release(tmp_path)
    train = tmp_path / "train"
    train.mkdir()
    (train / "edges.jsonl").write_text("a\n")      # a real split, no stamp
    assert check_derived(train, rel)[0] == "unstamped"


def test_a_derived_release_that_was_never_built_is_absent_not_stale(tmp_path):
    """Different situations need different instructions: 'run temporal_split'
    against 're-run temporal_split'. A checker that returns one status for
    both makes the caller print the wrong one."""
    from kgav.provenance import check_derived

    rel = _release(tmp_path)
    assert check_derived(tmp_path / "never-built", rel)[0] == "absent"


def test_the_derived_stamp_records_which_protocol_cut_the_split(tmp_path):
    """The three --undated protocols write to ONE output directory, so the
    stamp is the only record of which split is currently on disk."""
    from kgav.provenance import derived_policy, write_derived

    rel = _release(tmp_path)
    train = tmp_path / "train"
    train.mkdir()
    write_derived(train, rel, "0.12.0", cutoff=2021, undated="computed_only",
                  overlap_policy="drop-from-test")

    pol = derived_policy(train)
    assert pol["undated"] == "computed_only"
    assert pol["cutoff"] == 2021
    assert PROVENANCE_KEY not in pol, "the policy view must not leak the stamp"


def test_write_split_stamps_the_training_graph(tmp_path):
    """The stamp is written by write_split itself, not by its caller: a second
    consumer of build_split must not be able to produce an unstamped release."""
    import json as _json

    from kgav.provenance import check_derived
    from kgav.temporal import build_split, write_split

    rel = tmp_path / "rel"
    rel.mkdir()
    edges = [{"subject": "DRUG:a", "object": "NCBITaxon:1",
              "predicate": "HAS_ANTIVIRAL_ACTIVITY_AGAINST",
              "publications": [], "sources": ["infores:x"],
              "properties": {"publication_date": "2019-01-01",
                             "ec50_nm": 100.0, "relation": "="}}]
    (rel / "edges.jsonl").write_text("\n".join(_json.dumps(e) for e in edges))
    (rel / "nodes.jsonl").write_text(_json.dumps(
        {"id": "DRUG:a", "category": "SmallMolecule", "properties": {}}) + "\n"
        + _json.dumps({"id": "NCBITaxon:1", "category": "OrganismTaxon",
                       "properties": {}}) + "\n")

    s = build_split(rel, 2021, {"HAS_ANTIVIRAL_ACTIVITY_AGAINST"}, "include")
    out = tmp_path / "train"
    write_split(s, rel, out, schema_version="0.12.0",
                policy={"undated": "include", "overlap_policy": "drop-from-test"})
    assert check_derived(out, rel)[0] == "ok"


def test_one_derived_release_path_across_the_three_scripts_that_touch_it():
    """temporal_split WRITES it, hard_negatives READS it, check_results
    AUDITS it -- three defaults in three files that must name one directory.
    Exactly the drift that pointed make check-results at data/results-corrected
    while every script wrote to data/results: an audit aimed one directory
    off returns a verdict about neither.
    """
    paths = {}
    for script, flag in (("temporal_split", "--out"),
                         ("hard_negatives", "--train")):
        src = (ROOT / "scripts" / f"{script}.py").read_text()
        m = re.search(rf'"{re.escape(flag)}".*?default=root / "([^"]+)"',
                      src, re.S)
        assert m, f"{script}.py declares no {flag} default"
        paths[script] = m.group(1)

    src = (ROOT / "scripts" / "check_results.py").read_text()
    m = re.search(r'args\.derived is not None\s*\\?\s*\n?\s*else '
                  r'\[root / "([^"]+)"\]', src)
    assert m, "check_results.py declares no default derived release"
    paths["check_results"] = m.group(1)

    assert len(set(paths.values())) == 1, (
        f"the derived release is named differently by each script: {paths}. "
        "The one that disagrees is either writing where nobody reads or "
        "auditing a directory nothing produces.")


def test_the_reproducing_block_regenerates_every_audited_results_file():
    """docs/evaluation.md's Reproducing block named four of the fourteen
    files make check-results audits. The other ten -- eight per-virus
    baselines, metapath_comparison.json and the --undated exclude protocol --
    had no command anywhere, so following the documented procedure after a
    rebuild left ten files describing a graph that no longer existed while
    the operator believed they had regenerated everything.

    Checked by script and by taxon rather than by filename, because the
    filenames are built at runtime from arguments.
    """
    from kgav import virus_register

    block = (ROOT / "docs" / "evaluation.md").read_text()
    m = re.search(r"## Reproducing\n(.*?)\n```\n", block, re.S)
    assert m, "could not find the Reproducing block"
    repro = m.group(1)

    for script in ("hard_negatives", "run_baselines", "temporal_split",
                   "compare_metapaths"):
        assert f"scripts/{script}.py" in repro, (
            f"{script}.py writes into data/results but the Reproducing block "
            "never invokes it, so its output cannot be regenerated by "
            "following the documented procedure")

    # run_baselines writes one file per virus and takes one virus per run.
    built = [str(v["taxon"]) for v in virus_register.load()["viruses"]]
    missing = [t for t in built if t not in repro]
    assert not missing, (
        f"taxa {missing} are built by the spine but never appear in the "
        "Reproducing block, so their baselines_NCBITaxon_*.json files go "
        "stale with no command to regenerate them")


def test_the_reproducing_block_cuts_the_split_before_reading_it():
    """hard_negatives.py reads the training graph temporal_split.py writes.
    The block used to invoke hard_negatives first, so on a fresh rebuild its
    temporal column came from the PREVIOUS release's split -- and the result
    was then stamped against the current one.
    """
    block = (ROOT / "docs" / "evaluation.md").read_text()
    m = re.search(r"## Reproducing\n(.*?)\n```\n", block, re.S)
    assert m
    repro = m.group(1)
    assert repro.index("scripts/temporal_split.py") < \
        repro.index("scripts/hard_negatives.py"), (
        "the training graph must be cut before the script that reads it runs")
