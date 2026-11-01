"""Gate for the hard-negative driver: the two label sides are filtered by one
function, the metapath reach survives into the results, and the no-op control
is actually compared against something.
"""
import importlib.util
import json
import shlex
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _mod():
    spec = importlib.util.spec_from_file_location(
        "hard_negatives", ROOT / "scripts" / "hard_negatives.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


class _Labels:
    def __init__(self, pos, neg):
        self.positives, self.negatives = pos, neg


V = "NCBITaxon:2697049"


def _labels():
    return _Labels({V: {"A", "B", "C", "D"}}, {V: {"W", "X", "Y", "Z"}})


# si: A and W are selective, B and X have an index but are cytotoxic,
# C, D, Y, Z have no paired cytotoxicity data at all.
SI = {("A", V): 50.0, ("B", V): 2.0, ("W", V): 80.0, ("X", V): 1.5}


def test_default_filters_positives_only_which_is_not_like_for_like():
    """The published 0.826 is computed this way, and the asymmetry is real.

    Requiring a verified index of a POSITIVE requires paired same-document
    CC50/EC50, which selects for well-characterised compounds -- the ones most
    likely to carry the target annotations M1 reads. Measured on v0.1 the
    positive set halves 1,726 -> 808 and M1 rises 0.729 -> 0.823, while
    cov_pos goes 51.5% -> 70.3% and cov_neg stays at 5.3%.
    """
    m = _mod()
    lab = _labels()
    args = SimpleNamespace(selectivity="selective-only", si_threshold=10.0,
                           selectivity_applies="positives")
    m._filter_labels(lab, SI, args)
    assert lab.positives[V] == {"A"}          # only the selective active
    assert lab.negatives[V] == {"W", "X", "Y", "Z"}   # untouched


def test_both_applies_the_same_requirement_to_negatives():
    """The control. A compound must clear the same evidence bar to be a
    negative as to be a positive."""
    m = _mod()
    lab = _labels()
    args = SimpleNamespace(selectivity="selective-only", si_threshold=10.0,
                           selectivity_applies="both")
    m._filter_labels(lab, SI, args)
    assert lab.positives[V] == {"A"}
    assert lab.negatives[V] == {"W"}


def test_verified_only_keeps_cytotoxic_compounds_on_both_sides():
    """verified-only isolates 'has matched data' from 'passed the threshold',
    so a cytotoxic compound with an index stays in."""
    m = _mod()
    lab = _labels()
    args = SimpleNamespace(selectivity="verified-only", si_threshold=10.0,
                           selectivity_applies="both")
    m._filter_labels(lab, SI, args)
    assert lab.positives[V] == {"A", "B"}
    assert lab.negatives[V] == {"W", "X"}


def test_one_helper_serves_both_the_cross_sectional_and_temporal_sides():
    """The temporal block held a second hand-written copy of the loop, which
    is how the two sides came to differ without anyone deciding they should."""
    src = (ROOT / "scripts" / "hard_negatives.py").read_text()
    assert src.count("_filter_labels(") >= 3      # definition + both call sites
    assert "si_t.get(" not in src                 # the duplicated loop is gone


def test_cross_sectional_reach_is_not_discarded_before_it_is_written():
    """Batch 2a reset cross_reach to {} one line AFTER run() returned it, so
    the cross-sectional metapath reach never reached the results file -- the
    exact thing 2a was written to record."""
    src = (ROOT / "scripts" / "hard_negatives.py").read_text()
    assigned = src.index("cross, cross_reach = run(")
    reinit = src.find("cross_reach: dict = {}")
    assert reinit == -1 or reinit < assigned, \
        "cross_reach is re-initialised after run() set it"


# --------------------------------------------------------- the no-op control
#
# docs/evaluation.md states that `--selectivity all --selectivity-applies
# both` "reproduces the row-2 figure, so the flag does nothing it should
# not." The Reproducing block ran it and NOTHING COMPARED THE RESULT. Every
# hard_negatives.py variant writes to one hard_negatives.json, so by the
# time the control ran, row 2's file had been overwritten two steps earlier
# -- and the control's own output was overwritten one step later by row 5.
# About two and a half minutes per rebuild producing a file nobody read,
# reported in the document in the past tense as a check that had passed.
#
# It runs here instead, on every commit rather than once per rebuild, on a
# fixture small enough to be cheap and large enough to answer the question.
#
# TWO ARMS, and the second one is the point. That two documents agree proves
# nothing unless the fixture could have shown them disagreeing: the smoke
# release yields `cross_sectional: []` and would pass the agreement arm by
# comparing two empty documents, which is the fixture-too-small failure this
# suite keeps finding. So the same pair of settings must DIVERGE under
# `--selectivity selective-only`, where the filter actually runs.

VIRUS = "NCBITaxon:2697049"
VIRUS_LABEL = "Test virus"
PROTEINS = [f"UniProtKB:P{i}" for i in range(1, 5)]

#: 12 selective actives and 4 cytotoxic ones; 12 selective inactives and 6
#: with no verified index. Every arm leaves both classes at or above
#: MIN_CLASS_SIZE, so no arm is skipped for want of labels -- a skip would
#: make the two documents agree by both being empty.
COMPOUNDS: list[tuple[str, float | None]] = (
    [(f"INCHIKEY:A{i:02d}", 50.0) for i in range(1, 13)]
    + [(f"INCHIKEY:A{i:02d}", 2.0) for i in range(13, 17)]
    + [(f"INCHIKEY:N{i:02d}", 80.0) for i in range(1, 13)]
    + [(f"INCHIKEY:N{i:02d}", None) for i in range(13, 19)]
)


#: Structure is pre-cutoff and every label is post-cutoff, so a 2021 split
#: leaves the traversal intact and puts the whole test set on the far side --
#: which is what gives the temporal column something to evaluate.
STRUCTURE_DATE = "2020-01-01"
LABEL_DATE = "2022-06-01"


def _edge(subject, predicate, obj, qualifiers, date=STRUCTURE_DATE):
    return {"subject": subject, "predicate": predicate, "object": obj,
            "qualifiers": qualifiers, "primary_knowledge_source": "infores:test",
            "evidence_tier": 1, "first_asserted_date": date}


def _build_release(out: Path) -> Path:
    nodes = [{"id": VIRUS, "class": "OrganismTaxon",
              "properties": {"label": VIRUS_LABEL}}]
    nodes += [{"id": p, "class": "Protein",
               "properties": {"taxon_id": "NCBITaxon:9606", "is_viral": False}}
              for p in PROTEINS]
    nodes += [{"id": d, "class": "SmallMolecule", "properties": {}}
              for d, _ in COMPOUNDS]

    # Two of the four proteins are host dependencies, so M2 reaches some
    # compounds and not others. A channel that reaches everything or nothing
    # would make coverage constant and the comparison less sensitive.
    edges = [_edge(p, "HOST_FACTOR_FOR", VIRUS, {"direction": "dependency"})
             for p in PROTEINS[:2]]
    for k, (drug, si) in enumerate(COMPOUNDS):
        targets = [PROTEINS[k % 4]]
        if k % 3 == 0:                      # varied degree: the degree
            targets.append(PROTEINS[(k + 1) % 4])   # baseline is not constant
        # Every third compound cites one paper on BOTH its evidence and its
        # label, so --publication-disjoint has something to withhold and the
        # filtered graph is really a different graph. Without this the two
        # would be identical and a cache that ignored the filter would go
        # unnoticed.
        pub = [f"PMID:{k}"] if k % 3 == 0 else []
        edges += [dict(_edge(drug, "TARGETS", p, {}), publications=pub)
                  for p in targets]
        active = drug.startswith("INCHIKEY:A")
        q = {"assay_type": "cell_based_antiviral",
             "ec50_nm": 200.0 if active else 50_000.0,
             "relation": "=" if active else ">"}
        if si is not None:
            q["selectivity_index"] = si
        edges.append(dict(_edge(drug, "HAS_ANTIVIRAL_ACTIVITY_AGAINST", VIRUS,
                                q, date=LABEL_DATE), publications=pub))

    _write(out, nodes, edges)
    return out


def _write(out: Path, nodes: list[dict], edges: list[dict]) -> Path:
    out.mkdir(parents=True, exist_ok=True)
    (out / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n, sort_keys=True) for n in nodes))
    (out / "edges.jsonl").write_text(
        "\n".join(json.dumps(e, sort_keys=True) for e in edges))
    return out


def _build_train(release: Path, out: Path) -> Path:
    """A stamped pre-cutoff split, so the temporal column actually runs.

    DELIBERATELY NOT A COPY. Four compounds lose their TARGETS edges here, so
    the training graph's DWPC differs from the full graph's -- which is what
    makes it possible to notice a cache handing one graph's traversal to the
    other. Against an identical copy that mistake would be invisible.
    """
    from kgav.provenance import write_derived
    from kgav.schema import load_schema

    nodes = [json.loads(ln) for ln
             in (release / "nodes.jsonl").read_text().splitlines() if ln.strip()]
    severed = {d for d, _ in COMPOUNDS[-4:]}
    edges = [e for e in (json.loads(ln) for ln
                         in (release / "edges.jsonl").read_text().splitlines()
                         if ln.strip())
             # A 2021 cut drops every label (all post-cutoff) and keeps the
             # structure (all pre-cutoff).
             if e["predicate"] != "HAS_ANTIVIRAL_ACTIVITY_AGAINST"
             and not (e["predicate"] == "TARGETS" and e["subject"] in severed)]
    _write(out, nodes, edges)
    write_derived(out, release, load_schema().version,
                  cutoff=2021, undated="include")
    return out


def _run(release: Path, results: Path, selectivity: str, applies: str) -> dict:
    """The real entry point, as the block invokes it -- not main() imported.

    main() calls parse_args() with no argv, so importing it would test a
    different surface from the one the document runs.
    """
    proc = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "hard_negatives.py"),
         "--release", str(release),
         # Absent on purpose: the temporal column needs a stamped training
         # graph, both arms are denied one identically, and the skip and its
         # reason are recorded in the document being compared.
         "--train", str(release.parent / "no-training-graph"),
         "--results", str(results),
         "--selectivity", selectivity,
         "--selectivity-applies", applies],
        cwd=ROOT, capture_output=True, text=True, check=False)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    doc = json.loads((results / "hard_negatives.json").read_text())
    # The one field that moves between two runs of the SAME configuration.
    # Dropped by name rather than by dropping `provenance`, so a difference
    # in release fingerprint or code commit would still show up.
    doc["provenance"].pop("generated_utc", None)
    return doc


@pytest.fixture(scope="module")
def arms(tmp_path_factory) -> dict[tuple[str, str], dict]:
    release = _build_release(tmp_path_factory.mktemp("hn_release") / "release")
    out = tmp_path_factory.mktemp("hn_results")
    return {(f, a): _run(release, out / f"{f}_{a}", f, a)
            for f in ("all", "selective-only")
            for a in ("positives", "both")}


def test_the_fixture_actually_evaluates_something(arms):
    """GUARDS EVERY COMPARISON BELOW. If the virus were skipped for want of
    labels, or no channel reached a compound, the documents would agree by
    recording nothing -- a pass that means the opposite of what it looks
    like."""
    per_virus = arms[("all", "positives")]["cross_sectional"]
    assert list(per_virus) == [VIRUS_LABEL], per_virus
    per = per_virus[VIRUS_LABEL]
    assert {"M2", "COMBINED", "degree"} <= set(per), sorted(per)
    floor = _mod().MIN_CLASS_SIZE
    assert per["degree"]["n_pos"] >= floor
    assert per["degree"]["n_neg"] >= floor
    assert 0 < per["M2"]["coverage_pos"] < 1, "M2 reaches all or none"


def test_the_fixture_can_tell_the_two_settings_apart(arms):
    """THE CONTROL FOR THE CONTROL. Where the filter actually runs,
    --selectivity-applies both must change the evaluation. Without this arm
    the agreement below would be satisfied by any fixture too small to
    disagree."""
    a = arms[("selective-only", "positives")]
    b = arms[("selective-only", "both")]
    assert a["cross_sectional"] != b["cross_sectional"]
    pa = a["cross_sectional"][VIRUS_LABEL]["degree"]
    pb = b["cross_sectional"][VIRUS_LABEL]["degree"]
    assert pa["n_pos"] == pb["n_pos"], "the positives side is not the one 'both' adds"
    assert pb["n_neg"] < pa["n_neg"], "the negatives were never filtered"


def test_selectivity_applies_is_a_no_op_when_the_filter_is_all(arms):
    """The claim docs/evaluation.md makes. `all` short-circuits the filter
    block entirely, so the flag is recorded in the output and otherwise
    unread."""
    a, b = arms[("all", "positives")], arms[("all", "both")]
    differing = {k for k in set(a) | set(b) if a.get(k) != b.get(k)}
    # EXACTLY the recorded flag, not 'no evaluation difference'. Requiring
    # that one field to differ is what proves the two subprocesses really ran
    # different configurations; without it this assertion would be satisfied
    # just as well by comparing a document with itself.
    assert differing == {"selectivity_applies"}, sorted(differing)
    assert (a["selectivity_applies"], b["selectivity_applies"]) \
        == ("positives", "both")


# --------------------------------------------------------- the protocol batch
#
# dwpc_scores depends on the GRAPH alone. Every protocol in the Reproducing
# block varies only the labels -- which compounds count as active, which as
# negative, which screen they come from -- so five protocols over two graphs
# each is ten traversals of which four are distinct. --protocol runs them in
# one process and computes each distinct one once.
#
# THE ONLY THING THAT MATTERS IS THAT IT CHANGES NOTHING. A batch must write
# exactly what separate invocations write, after EVERY protocol and not only
# the last: the block relies on last-write-wins, so an interrupted batch must
# leave a document some protocol actually produced.

#: The block's five rows, in the block's order, with the fixture's own screen
#: standing in for infores:ncats-opendata.
BLOCK_PROTOCOLS = [
    "--selectivity selective-only",
    "--selectivity all",
    "--selectivity verified-only --selectivity-applies both",
    # Parenthesised, not a bare implicit concatenation: a stray comma here
    # would silently make this two protocols, and ruff's ISC004 flags the
    # shape for exactly that reason.
    ("--selectivity verified-only --selectivity-applies both "
     "--publication-disjoint"),
    "--selectivity all --publication-disjoint --label-source infores:test",
]


def _invoke(release, train, results, *flags):
    return subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "hard_negatives.py"),
         "--release", str(release), "--train", str(train),
         "--results", str(results), *flags],
        cwd=ROOT, capture_output=True, text=True, check=False)


def _doc(results: Path) -> dict:
    doc = json.loads((results / "hard_negatives.json").read_text())
    doc["provenance"].pop("generated_utc", None)
    return doc


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    """A release WITH a stamped training graph, so the temporal column runs.

    The arms fixture above deliberately has none; here the training graph is
    the point, because it is the second of the two graphs whose traversal the
    batch is meant to share.
    """
    base = tmp_path_factory.mktemp("hn_world")
    release = _build_release(base / "release")
    return release, _build_train(release, base / "train")


@pytest.fixture(scope="module")
def separate(world, tmp_path_factory):
    """Each protocol run on its own, snapshotting the file after every one."""
    release, train = world
    out = tmp_path_factory.mktemp("hn_separate") / "results"
    snapshots = []
    for spec in BLOCK_PROTOCOLS:
        proc = _invoke(release, train, out, *shlex.split(spec))
        assert proc.returncode == 0, proc.stdout + proc.stderr
        snapshots.append(_doc(out))
    return snapshots


def test_the_training_graph_is_not_the_full_graph(separate):
    """GUARDS THE CACHE KEY. The two graphs share a release and differ only
    in their path, so a key that collided would hand the temporal column the
    cross-sectional traversal. Against an identical copy that would be
    invisible; here the training graph is four compounds short of its
    TARGETS edges and the numbers must differ."""
    doc = separate[1]
    assert doc["temporal_skipped"] is None, doc["temporal_skipped"]
    assert doc["temporal"], "the temporal column did not run"
    assert (doc["metapath_reach"]["cross_sectional"]
            != doc["metapath_reach"]["temporal"])
    assert (doc["cross_sectional"][VIRUS_LABEL]["M2"]["auc"]
            != doc["temporal"][VIRUS_LABEL]["M2"]["auc"])


@pytest.mark.parametrize("n", range(1, len(BLOCK_PROTOCOLS) + 1))
def test_a_batch_writes_what_separate_invocations_write(
        n, world, separate, tmp_path_factory):
    """Truncated at every length, so the intermediate documents are compared
    and not only the last. Getting the final file right while clobbering the
    ones before it would still break the block, which depends on each run
    leaving a complete document behind."""
    release, train = world
    out = tmp_path_factory.mktemp(f"hn_batch{n}") / "results"
    flags = [f for spec in BLOCK_PROTOCOLS[:n] for f in ("--protocol", spec)]
    proc = _invoke(release, train, out, *flags)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert _doc(out) == separate[n - 1]


def test_the_batch_shares_the_traversals_it_is_there_to_share(
        world, tmp_path_factory):
    """The whole point, asserted as a count rather than as a stopwatch.

    Five protocols touch two graphs each: ten traversals, four distinct --
    the full and training graphs, with and without the publication filter.
    An assertion of 'fewer than ten' would pass on a cache that worked once,
    so this pins the number."""
    release, train = world
    out = tmp_path_factory.mktemp("hn_count") / "results"
    flags = [f for spec in BLOCK_PROTOCOLS for f in ("--protocol", spec)]
    proc = _invoke(release, train, out, *flags)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    computed = proc.stdout.count("computing the graph and DWPC")
    reused = proc.stdout.count("reusing the graph and DWPC")
    # Four distinct: {full, training} x {filtered, not}. The total is derived
    # rather than written down, so adding a protocol changes one number here
    # and not two.
    assert computed == 4, proc.stdout
    assert computed + reused == 2 * len(BLOCK_PROTOCOLS), proc.stdout


def test_a_shared_flag_reaches_every_protocol(world, tmp_path_factory):
    """--release and --results are given once. If they did not reach each
    protocol, every protocol would silently read data/releases/v0.1 -- the
    real release, not the fixture -- and this suite would depend on one
    being present."""
    release, train = world
    out = tmp_path_factory.mktemp("hn_shared") / "results"
    proc = _invoke(release, train, out, "--protocol", "--selectivity all")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert _doc(out)["provenance"]["release"]["id"] == release.name


def test_a_protocol_overrides_a_shared_flag(world, tmp_path_factory):
    """Each spec is parsed AFTER the shared flags, so a protocol can set its
    own. Without this the batch could express less than five invocations
    can, which would make it a different protocol set rather than the same
    one run faster."""
    release, train = world
    out = tmp_path_factory.mktemp("hn_override") / "results"
    proc = _invoke(release, train, out, "--inactive-above-nm", "10000",
                   "--protocol", "--selectivity all --inactive-above-nm 99999")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert _doc(out)["inactive_above_nm"] == 99999.0


def test_a_protocol_cannot_nest(world, tmp_path_factory):
    release, train = world
    out = tmp_path_factory.mktemp("hn_nest") / "results"
    proc = _invoke(release, train, out,
                   "--protocol", "--selectivity all --protocol '--selectivity all'")
    assert proc.returncode != 0
    assert "cannot nest" in proc.stderr


def test_a_failing_protocol_does_not_stop_the_ones_after_it(
        world, tmp_path_factory):
    """Run separately these are five commands in a block with no `set -e`, so
    one exiting 1 does not stop the next and the LAST protocol is still the
    one left on disk. A batch that aborted early would change which numbers
    survive a rebuild, which is the one thing it must not do."""
    release, train = world
    out = tmp_path_factory.mktemp("hn_rc") / "results"
    proc = _invoke(
        release, train, out,
        "--protocol", "--selectivity all --label-source infores:nonesuch",
        "--protocol", "--selectivity verified-only --selectivity-applies both")
    assert proc.returncode == 1, "a protocol that wrote nothing was not reported"
    assert "NO POSITIVES" in proc.stdout
    # The protocol AFTER the failure still ran, and its document is the one
    # on disk.
    doc = _doc(out)
    assert doc["selectivity_filter"] == "verified-only"
    assert doc["label_source"] is None


def test_shared_argv_drops_only_the_protocol_tokens():
    """The reconstruction that makes shared flags work. Getting it wrong by
    one token would feed a protocol string into the shared list, where it
    would be parsed as flags rather than as a protocol."""
    m = _mod()
    assert m.shared_argv(
        ["--release", "r", "--protocol", "--selectivity all",
         "--results", "o", "--protocol", "--selectivity none"]
    ) == ["--release", "r", "--results", "o"]
    assert m.shared_argv(["--protocol=--selectivity all", "--cutoff", "2021"]) \
        == ["--cutoff", "2021"]
    assert m.shared_argv(["--release", "r"]) == ["--release", "r"]


def test_run_uses_the_dwpc_it_is_given_instead_of_recomputing(world):
    """The saving itself, asserted directly.

    Everything above would still pass if run() accepted a DWPC and then
    quietly recomputed its own: the numbers would be identical and the
    cache's own stdout lines would still say 'reusing'. This fails instead
    of merely being slow."""
    from kgav.baselines import Graph, dwpc_scores, metapath_reach
    from kgav.labels import build_labels
    from kgav.schema import load_schema

    release, _ = world
    schema = load_schema()
    g = Graph.load(release, skip_predicates=schema.held_out_predicates(),
                   symmetry=schema.symmetry())
    dwpc = dwpc_scores(g, schema.metapaths, schema.retrieval_limits)
    reach = metapath_reach(g, schema.metapaths, dwpc)

    m = _mod()

    def no_dwpc(*a, **k):
        raise AssertionError("run() recomputed a DWPC it had been handed")

    def no_reach(*a, **k):
        raise AssertionError("run() recomputed a reach it had been handed")

    m.dwpc_scores, m.metapath_reach = no_dwpc, no_reach
    out, out_reach = m.run(g, schema, build_labels(release), "t",
                           dwpc=dwpc, reach=reach)
    assert out, "the fixture evaluated nothing, so nothing was proved"
    assert out_reach == reach


def test_run_still_computes_its_own_when_given_none(world):
    """The other half: every existing caller passes neither, and must keep
    working. A signature change that made them mandatory would be caught
    here rather than by a script failing at the end of a rebuild."""
    from kgav.baselines import Graph
    from kgav.labels import build_labels
    from kgav.schema import load_schema

    release, _ = world
    schema = load_schema()
    g = Graph.load(release, skip_predicates=schema.held_out_predicates(),
                   symmetry=schema.symmetry())
    out, reach = _mod().run(g, schema, build_labels(release), "t")
    assert out and reach


# ------------------------------------------------ why a column scored nothing
#
# TWO SILENCES, FOUND ON 2026-11-02 and neither distinguishable from a
# healthy run by reading the results file.
#
#   A VIRUS BELOW MIN_CLASS_SIZE printed to stdout and was recorded nowhere.
#   Rows 1-4 each drop viruses this way: row 2 scored 10 of 15, and JEV,
#   TBEV, MERS-CoV, HCoV-NL63 and SARS-CoV left no trace at all. A reader
#   counting entries sees ten viruses and no reason to think there were
#   fifteen.
#
#   AN EMPTY LABEL SET printed NOTHING -- the per-virus loop never ran. Row
#   5's temporal column is this: NCATS-only labels carry no post-cutoff
#   dates, so `temporal: {}` beside `temporal_skipped: null` read as a
#   column that computed fine, while docs/evaluation.md quoted temporal
#   AUCs for that protocol.
def test_a_virus_below_the_class_floor_is_recorded_with_its_reason(world):
    """The reason names WHICH class was short and by how much. "needs 10
    of each" sends a reader to look for both when one is the problem --
    and on the real release four of the five skips are inactives only."""
    from kgav.baselines import Graph
    from kgav.labels import build_labels
    from kgav.schema import load_schema

    release, _ = world
    schema = load_schema()
    g = Graph.load(release, skip_predicates=schema.held_out_predicates(),
                   symmetry=schema.symmetry())
    notes: dict = {}
    # A floor above the fixture's class sizes, so every virus is skipped.
    out, _ = _mod().run(g, schema, build_labels(release), "t",
                        k_min=10_000, notes=notes)
    assert not out, "nothing should score under an impossible floor"
    assert notes["unscored"], "every virus was skipped and none was recorded"
    assert notes["empty_label_set"] is False, (
        "the label set is not empty -- the viruses were unscorable, which "
        "is a different statement and needs a different fix")
    # THE REASON IS DERIVED, NOT PATTERN-MATCHED. `"active" in why or
    # "inactive" in why` survived a mutant that replaced the ACTIVE branch
    # with "short by N": every virus here is short on BOTH classes, so the
    # surviving inactive branch satisfied the `or` and the broken half went
    # unexamined. Each class is now checked on its own, against the counts.
    labels = build_labels(release)
    for virus, why in notes["unscored"].items():
        pos = next(len(p) for v, p in labels.positives.items()
                   if (g.nodes.get(v, {}).get("properties", {})
                       or {}).get("label", v) == virus)
        neg = next((len(n) for v, n in labels.negatives.items()
                    if (g.nodes.get(v, {}).get("properties", {})
                        or {}).get("label", v) == virus), 0)
        assert f"active {pos} < 10000" in why, f"{virus}: {why!r}"
        assert f"inactive {neg} < 10000" in why, f"{virus}: {why!r}"
    assert notes["labelled_viruses"] == len(notes["unscored"])


def test_an_empty_label_set_is_recorded_as_such_not_as_an_empty_result(world):
    """ROW 5'S TEMPORAL COLUMN. `empty_label_set` is what separates "no
    virus could be scored" from "no virus was labelled at all": the first
    is a floor to maybe lower, the second is a fact about the label
    source that no rerun will change."""
    from kgav.baselines import Graph
    from kgav.labels import build_labels
    from kgav.schema import load_schema

    release, _ = world
    schema = load_schema()
    g = Graph.load(release, skip_predicates=schema.held_out_predicates(),
                   symmetry=schema.symmetry())
    empty = build_labels(release, sources={"infores:not-a-real-source"})
    assert not empty.positives, "the fixture's labels survived the filter"
    notes: dict = {}
    out, _ = _mod().run(g, schema, empty, "t", notes=notes)
    assert not out
    assert notes["empty_label_set"] is True
    assert notes["labelled_viruses"] == 0
    assert notes["unscored"] == {}, (
        "an empty label set has no virus to report unscored; filling this "
        "would make the two silences indistinguishable again")
    # AGAINST THE LABEL OBJECT, not just "> 0". A mutant hard-coding
    # active to 0 survived a bare positivity check on one other field.
    assert notes["labels"] == {
        "active": empty.stats.get("active", 0),
        "inactive": empty.stats.get("inactive", 0),
        "other_source": empty.stats.get("other_source", 0)}, (
        "the counts that explain the emptiness are not the label set's own, "
        "so the file says it is empty and not why")
    assert notes["labels"]["other_source"] > 0, (
        "nothing was dropped as another source, so this fixture is not "
        "exercising the case row 5's temporal column is in")


def test_the_notes_say_nothing_is_wrong_when_nothing_is(world):
    """A field that is always populated is a field nobody reads."""
    from kgav.baselines import Graph
    from kgav.labels import build_labels
    from kgav.schema import load_schema

    release, _ = world
    schema = load_schema()
    g = Graph.load(release, skip_predicates=schema.held_out_predicates(),
                   symmetry=schema.symmetry())
    labels = build_labels(release)
    notes: dict = {}
    out, _ = _mod().run(g, schema, labels, "t", notes=notes)
    assert out, "the fixture scored nothing, so nothing is proved"
    assert notes["empty_label_set"] is False
    assert notes["unscored"] == {}
    # ON A NON-EMPTY SET. The equality check in the empty-set test above
    # cannot see a mutant that hard-codes `active` to 0, because the
    # correct answer there IS 0. This fixture has actives, so it can.
    assert labels.stats.get("active", 0) > 0, "fixture has no actives"
    assert notes["labels"] == {
        "active": labels.stats.get("active", 0),
        "inactive": labels.stats.get("inactive", 0),
        "other_source": labels.stats.get("other_source", 0)}
    assert notes["labelled_viruses"] == len(labels.positives)


def test_the_notes_default_is_none_and_not_a_shared_dict():
    """A MUTABLE DEFAULT WOULD BE SHARED BY EVERY CALL. `notes: dict = {}`
    passed every test here -- nothing observes the default object -- while
    being a real defect: one dict would accumulate across both columns of
    every protocol in a --protocol run, so row 5's notes would carry row
    1's unscored viruses. The signature is checked directly, because the
    bug is invisible from the outside until it has already corrupted a
    results file."""
    import inspect

    sig = inspect.signature(_mod().run)
    default = sig.parameters["notes"].default
    assert default is None, (
        f"run()'s `notes` defaults to {default!r}; a mutable default is "
        f"shared by every call that omits it")


def test_run_without_notes_is_unchanged(world):
    """Every existing caller passes none, exactly as with dwpc and reach.
    The idiom is Assembly.stats: fill in place, never change the arity."""
    from kgav.baselines import Graph
    from kgav.labels import build_labels
    from kgav.schema import load_schema

    release, _ = world
    schema = load_schema()
    g = Graph.load(release, skip_predicates=schema.held_out_predicates(),
                   symmetry=schema.symmetry())
    out, reach = _mod().run(g, schema, build_labels(release), "t")
    assert out and reach


def test_both_columns_notes_reach_the_results_file():
    """THE POINT OF ALL OF IT. A note computed and not written is the
    silence it was meant to end, one layer in -- which is how
    `temporal_skipped` came to be the only one of the three states the
    file could express."""
    src = Path(_mod().__file__).read_text()
    for key in ('"cross_sectional_notes": cross_notes',
                '"temporal_notes": temporal_notes'):
        assert key in src, f"{key} is not written into the results file"
    # And distinct from the pre-existing third thing.
    assert '"temporal_skipped": temporal_skipped' in src, (
        "temporal_skipped means the column had NO TRAINING GRAPH to vouch "
        "for, which is neither of the two states above")
    assert src.index("temporal_notes: dict = {}") < src.index(
        "notes=temporal_notes"), "temporal_notes is used before it exists"
