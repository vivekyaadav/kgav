"""Gate for the pathway landscape.

This is the analysis the paper rests on, so two things must not be wrong.
The census must account for every molecule, because a landscape over a
minority that called itself "all drugs" is false in its first line. And the
matched test must actually control for degree, because the confound that
killed the scoring is the same one that inflates a naive enrichment.
"""
from __future__ import annotations

import json
import os
import random
import sys
from importlib import import_module
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

L = import_module("pathway_landscape")


def _release(tmp_path: Path, n_drugged=6, n_undrugged=10, n_deps=4) -> Path:
    """Proteins, pathways, a virus, and molecules of three kinds: targeting,
    INHIBITS-only, and neither."""
    d = tmp_path / "rel"
    d.mkdir(parents=True, exist_ok=True)
    nodes = [{"id": "NCBITaxon:2697049", "class": "OrganismTaxon",
              "properties": {}},
             {"id": "REACT:P1", "class": "Pathway", "properties": {}},
             {"id": "UniProtKB:VIRAL", "class": "Protein",
              "properties": {"taxon_id": "NCBITaxon:2697049",
                             "is_viral": True}}]
    edges = []
    for i in range(n_drugged + n_undrugged):
        p = f"UniProtKB:H{i}"
        nodes.append({"id": p, "class": "Protein",
                      "properties": {"taxon_id": "NCBITaxon:9606",
                                     "is_viral": False}})
        edges.append((p, "PARTICIPATES_IN", "REACT:P1", {}))
        if i < n_deps:
            edges.append((p, "HOST_FACTOR_FOR", "NCBITaxon:2697049",
                          {"direction": "dependency"}))
    # a restriction factor, which must NOT be counted as a dependency
    edges.append(("UniProtKB:H0", "HOST_FACTOR_FOR", "NCBITaxon:2697049",
                  {"direction": "restriction"}))
    for i in range(n_drugged):
        c = f"INCHIKEY:T{i}"
        nodes.append({"id": c, "class": "SmallMolecule", "properties": {}})
        edges.append((c, "TARGETS", f"UniProtKB:H{i}",
                      {"direction": "inhibitor"} if i == 0 else {}))
    for i in range(3):
        c = f"INCHIKEY:I{i}"
        nodes.append({"id": c, "class": "SmallMolecule", "properties": {}})
        edges.append((c, "INHIBITS", "UniProtKB:VIRAL", {}))
    # A COMPOUND WITH BOTH. The first fixture had none, so the census test
    # asserted that four counts summed to the population and passed; on
    # v0.1, where 47 compounds carry both, they summed to 21,714 against
    # 21,667. The fixture is what let that through.
    edges.append(("INCHIKEY:T0", "INHIBITS", "UniProtKB:VIRAL", {}))
    for i in range(5):
        nodes.append({"id": f"INCHIKEY:N{i}", "class": "SmallMolecule",
                      "properties": {}})
    (d / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in nodes))
    (d / "edges.jsonl").write_text("\n".join(json.dumps(
        {"subject": s, "predicate": p, "object": o, "qualifiers": q,
         "primary_knowledge_source": "infores:chembl", "evidence_tier": 1,
         "first_asserted_date": "2021-01-01"}) for s, p, o, q in edges))
    return d


def _graph(rel: Path):
    from kgav.baselines import Graph
    from kgav.schema import load_schema
    return Graph.load(rel, symmetry=load_schema().symmetry())


# ----------------------------------------------------------------- census
def test_the_census_partitions_the_population(tmp_path):
    """THE FIRST LINE OF THE PAPER, and it was wrong. `both` is a SUBSET of
    the analysable set, not a fourth bucket: the three parts partition, and
    adding `both` over-counts. On v0.1 that printed 21,714 against 21,667
    molecules."""
    g = _graph(_release(tmp_path))
    c = L.census(g)
    assert c["molecules"] == 6 + 3 + 5
    assert c["analysable"] == 6          # includes the one with both
    assert c["inhibits_only"] == 3
    assert c["neither"] == 5
    assert c["both"] == 1
    assert c["parts_sum"] == c["molecules"]
    assert c["analysable"] + c["inhibits_only"] + c["neither"] + c["both"] \
        != c["molecules"], "the four must NOT sum; that was the bug"


def test_both_is_inside_the_analysable_set(tmp_path):
    g = _graph(_release(tmp_path))
    c = L.census(g)
    assert c["both"] <= c["analysable"]


def test_the_census_counts_edges_that_state_a_direction(tmp_path):
    """Potency is not agonist versus antagonist. Targeting a dependency
    factor is only a claim about helping if the direction is known, and the
    paper has to show how often it is not."""
    g = _graph(_release(tmp_path))
    c = L.census(g)
    assert c["target_edges"] == 6
    assert c["directed_edges"] == 1


# ------------------------------------------------------------- the surface
def test_a_restriction_factor_is_not_a_dependency_factor(tmp_path):
    """Opposite claims. Pooling them would make every overlap count
    meaningless, and the ORCS layer's own module says a sign error here
    produces a confidently wrong recommendation."""
    g = _graph(_release(tmp_path, n_deps=4))
    deps = L.dependency_factors(g, "NCBITaxon:2697049")
    assert len(deps) == 4
    assert "UniProtKB:H0" in deps          # it is BOTH in this fixture
    g2 = _graph(_release(tmp_path / "b", n_deps=0))
    assert L.dependency_factors(g2, "NCBITaxon:2697049") == set()


def test_oversized_pathways_are_dropped_at_the_schemas_cap(tmp_path):
    g = _graph(_release(tmp_path, n_drugged=6, n_undrugged=10))
    assert len(L.pathway_members(g, cap=None)) == 1
    assert L.pathway_members(g, cap=4) == {}
    assert len(L.pathway_members(g, cap=100)["REACT:P1"]) == 16


def test_drug_targets_are_proteins_not_molecules(tmp_path):
    g = _graph(_release(tmp_path))
    got = L.drug_targets(g)
    assert all(t.startswith("UniProtKB:") for t in got)
    assert len(got) == 6


# -------------------------------------------------------- the matched test
def test_degree_bins_spread_a_heavy_tail(tmp_path):
    """Equal-width bins on a heavy-tailed degree would put nearly every
    protein in the first bin and match nothing."""
    g = _graph(_release(tmp_path))
    universe = {n for n in g.nodes if g.klass(n) == "Protein"}
    bins = L.degree_bins(g, universe, bins=4)
    assert set(bins) == universe
    assert len(set(bins.values())) > 1


def test_the_matched_null_preserves_the_degree_profile(tmp_path):
    """THE POINT OF THE WHOLE TEST. If the null set did not carry the same
    degree profile as the observed set, the matching would be decorative
    and the enrichment would still be measuring study volume."""
    g = _graph(_release(tmp_path, n_drugged=6, n_undrugged=10))
    universe = {n for n in g.nodes if g.klass(n) == "Protein"}
    bins = L.degree_bins(g, universe, bins=4)
    observed = L.drug_targets(g)
    from collections import Counter
    want = Counter(bins[n] for n in observed)
    rng = random.Random(1)
    by_bin = {}
    for n in universe:
        by_bin.setdefault(bins[n], []).append(n)
    for b, k in want.items():
        assert k <= len(by_bin[b]), "a bin cannot supply its matched draw"
    null = L.matched_null(observed, universe, bins, set(), 50, rng)
    assert len(null) == 50
    assert set(null) == {0}, "no hits means every draw must score zero"


def test_the_matched_null_finds_hits_when_they_exist(tmp_path):
    g = _graph(_release(tmp_path))
    universe = {n for n in g.nodes if g.klass(n) == "Protein"}
    bins = L.degree_bins(g, universe, bins=4)
    null = L.matched_null(L.drug_targets(g), universe, bins, universe, 20,
                          random.Random(2))
    assert all(v == len(L.drug_targets(g)) for v in null), \
        "every drawn protein is a hit, so every draw must be full"


@pytest.mark.parametrize("observed,null,lo,hi", [
    (10, [0] * 99, 0.0, 0.011),       # never beaten: 1/100
    (0, [5] * 99, 0.99, 1.0),         # always beaten
    (5, [5] * 99, 0.99, 1.0),         # ties count against
])
def test_the_empirical_p_uses_the_plus_one_correction(observed, null, lo, hi):
    """Without it a p of exactly 0 is reportable, which claims more
    precision than the permutation count can support."""
    got = L.empirical_p(observed, null)
    assert lo <= got <= hi


def test_the_exact_tail_is_imported_not_reimplemented():
    """Phase 4 validated fisher_one_sided against scipy over 400 random
    tables. A second copy here would drift from the validated one."""
    src = (ROOT / "scripts" / "pathway_landscape.py").read_text()
    assert "from phase4_hcv_retrospective import fisher_one_sided" in src
    assert "def fisher_one_sided" not in src


# ------------------------------------------------------------- end to end
def test_it_runs_and_reports_both_tests(tmp_path, capsys):
    rel = _release(tmp_path)
    assert L.main(["--release", str(rel), "--permutations", "200"]) == 0
    out = capsys.readouterr().out
    assert "CENSUS" in out and "DRUGGED PROTEOME" in out
    assert "p naive" in out and "p matched" in out
    assert "NCBITaxon:2697049" in out


def test_it_never_reads_a_label_predicate():
    """The claim that the shared-PMID confound cannot apply rests on this.
    A label predicate appearing here would quietly reintroduce it."""
    src = (ROOT / "scripts" / "pathway_landscape.py").read_text()
    body = src.split('"""', 2)[2]
    for pred in ("HAS_ANTIVIRAL_ACTIVITY_AGAINST", "MEASURED_INACTIVE_AGAINST",
                 "IN_TRIAL_FOR", "TREATS", "build_labels"):
        assert pred not in body, f"{pred} reaches the analysis"


def test_a_release_with_no_dependency_factors_says_so(tmp_path, capsys):
    """An empty table is not a null result, and the two must not look
    alike."""
    rel = _release(tmp_path, n_deps=0)
    assert L.main(["--release", str(rel), "--permutations", "50"]) == 1
    out = capsys.readouterr().out
    assert "NO VIRUS HAS DEPENDENCY FACTORS" in out
    assert "not an" in out


def test_it_states_what_it_cannot_say(tmp_path, capsys):
    L.main(["--release", str(_release(tmp_path)), "--permutations", "50"])
    out = capsys.readouterr().out
    assert "NOT a claim that the drug would help" in out
    assert "ranks a drug" in out
    src = (ROOT / "scripts" / "pathway_landscape.py").read_text()
    assert "Read-only. Writes nothing, ranks no drug" in src


def test_a_missing_release_exits_two(tmp_path, capsys):
    assert L.main(["--release", str(tmp_path / "nope")]) == 2
    assert "no release at" in capsys.readouterr().out


# --------------------------------- defects found on the first v0.1 run
def test_the_null_universe_excludes_viral_proteins(tmp_path):
    """A viral protein can never be a human drug target nor a host
    dependency factor, so every one left in a degree bin is a guaranteed
    miss that lowers the matched expectation and flatters the observed
    overlap. The first version sampled from every Protein node."""
    g = _graph(_release(tmp_path))
    universe = L.human_proteins(g)
    assert "UniProtKB:VIRAL" not in universe
    assert all(not g.is_viral(n) for n in universe)
    every_protein = {n for n in g.nodes if g.klass(n) == "Protein"}
    assert universe < every_protein, "the fixture must contain a viral one"


def test_a_taxon_with_no_dependency_factors_is_named_not_dropped(
        tmp_path, capsys):
    """Three of eight taxa dropped out of the v0.1 table in silence, so a
    reader counted five rows against eight and had nothing to tell them
    why. The CRISPR layer's coverage is a scope limit for the paper."""
    rel = _release(tmp_path)
    nodes = [json.loads(l) for l in
             (rel / "nodes.jsonl").read_text().splitlines()]
    nodes.append({"id": "NCBITaxon:64320", "class": "OrganismTaxon",
                  "properties": {}})
    (rel / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in nodes))
    L.main(["--release", str(rel), "--permutations", "50"])
    out = capsys.readouterr().out
    assert "NCBITaxon:64320" in out
    assert "NO dependency" in out
    assert "scope" in out


def test_the_census_assertion_fires_if_the_parts_stop_partitioning(tmp_path):
    """The printed census carries its own check. If a future edge kind
    breaks the partition, the run must stop rather than print numbers that
    do not add up."""
    src = (ROOT / "scripts" / "pathway_landscape.py").read_text()
    assert 'assert c["parts_sum"] == c["molecules"]' in src


def test_main_actually_uses_the_human_only_universe():
    """A test on the helper alone is decorative: reverting the call site
    back to every Protein node left all 21 tests green. The caller is what
    decides which universe the null samples from."""
    src = (ROOT / "scripts" / "pathway_landscape.py").read_text()
    body = src.split('def main(', 1)[1]
    assert "universe = human_proteins(g)" in body
    assert 'g.klass(n) == "Protein"' not in body, \
        "main is building its own universe again, viral proteins included"


def test_coverage_separates_no_screen_from_restriction_only(tmp_path):
    """THE DISTINCTION THE v0.1 RUN COULD NOT MAKE. Three taxa dropped out
    of the overlap table and the report said only that they had no
    dependency factor. A taxon nothing screened and a taxon whose every hit
    was a restriction factor look identical from there, and they call for
    completely different work."""
    g = _graph(_release(tmp_path, n_deps=2))
    cov = L.host_factor_coverage(g)
    c = cov["NCBITaxon:2697049"]
    assert c["dependency"] == 2
    assert c["restriction"] == 1
    assert "NCBITaxon:64320" not in cov, "a taxon with no screen is absent"


def test_an_unstated_direction_is_counted_not_dropped(tmp_path):
    """A HOST_FACTOR_FOR edge with no direction is neither a dependency nor
    a restriction claim, and silently discarding it would make the totals
    disagree with the graph."""
    rel = _release(tmp_path)
    edges = [json.loads(l) for l in
             (rel / "edges.jsonl").read_text().splitlines()]
    edges.append({"subject": "UniProtKB:H9", "predicate": "HOST_FACTOR_FOR",
                  "object": "NCBITaxon:2697049", "qualifiers": {},
                  "primary_knowledge_source": "infores:biogrid-orcs",
                  "evidence_tier": 1, "first_asserted_date": "2021-01-01"})
    (rel / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in edges))
    cov = L.host_factor_coverage(_graph(rel))
    assert cov["NCBITaxon:2697049"]["unstated"] == 1


def test_the_coverage_table_flags_a_taxon_nothing_screened(tmp_path, capsys):
    rel = _release(tmp_path)
    nodes = [json.loads(l) for l in
             (rel / "nodes.jsonl").read_text().splitlines()]
    nodes.append({"id": "NCBITaxon:11060", "class": "OrganismTaxon",
                  "properties": {}})
    (rel / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in nodes))
    L.main(["--release", str(rel), "--permutations", "50"])
    out = capsys.readouterr().out
    assert "HOST-FACTOR COVERAGE" in out
    assert "NO screen reached this taxon" in out
    assert "restriction hits only has data" in out


# ------------------------------------------------- reproducibility
def test_the_same_seed_gives_the_same_answer_across_processes():
    """THE DEFECT THREE v0.1 RUNS EXPOSED. One release, one seed, and
    SARS-CoV-2's matched expectation came back 353.3, then 352.6, then
    352.4. Degree is heavily tied, so sorting on degree alone left ties in
    SET order, and the pools rng.sample draws from were built by iterating
    a set -- both vary between processes through string hash
    randomisation.

    A figure in a paper has to come back the same on a rerun, so this runs
    the computation under three different PYTHONHASHSEEDs and requires one
    answer. An in-process test cannot catch this: set order is stable
    within a process, which is exactly why it survived twenty-two tests.
    """
    import subprocess
    prog = (
        "import sys; sys.path.insert(0, 'scripts');"
        "import random;"
        "from pathway_landscape import degree_bins, matched_null;"
        "deg = {'P' + str(i): (i % 5) for i in range(60)};"
        "G = type('G', (), {'degree': deg});"
        "u = {'P' + str(i) for i in range(60)};"
        "b = degree_bins(G(), u, bins=4);"
        "print(matched_null({'P1','P7','P13','P22'}, u, b,"
        " {'P2','P9','P21','P33'}, 8, random.Random(0)))"
    )
    out = set()
    for hashseed in ("0", "1", "2", "12345"):
        env = {**os.environ, "PYTHONHASHSEED": hashseed}
        r = subprocess.run([sys.executable, "-c", prog], cwd=str(ROOT),
                           capture_output=True, text=True, env=env)
        assert r.returncode == 0, r.stderr
        out.add(r.stdout.strip())
    assert len(out) == 1, f"seed 0 gave {len(out)} different answers: {out}"


def test_degree_bins_breaks_ties_on_the_node_id():
    """Thousands of proteins share a low degree. Without the id in the
    key, which bin a tied protein lands in is set order."""
    src = (ROOT / "scripts" / "pathway_landscape.py").read_text()
    assert "key=lambda n: (g.degree.get(n, 0), n)" in src
    assert "for node in sorted(universe)" in src
    assert "wanted = sorted(want.items())" in src, \
        "a Counter iterates in insertion order; it must be sorted too"
