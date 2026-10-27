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


@pytest.mark.parametrize("observed,null,lo,hi", [
    (0, [5] * 99, 0.0, 0.011),        # never undercut: 1/100
    (10, [0] * 99, 0.99, 1.0),        # always undercut
    (5, [5] * 99, 0.99, 1.0),         # ties count against, same as greater
])
def test_the_lower_tail_counts_the_other_side(observed, null, lo, hi):
    """greater=False must count v <= observed, not v >= observed.

    Flipping the comparison left every test passing: nothing exercised the
    tail this flag exists to compute, so the agent would have reported the
    enrichment p twice under two different names.
    """
    assert lo <= L.empirical_p(observed, null, greater=False) <= hi


def test_the_two_tails_are_not_the_same_number():
    """The cheap check that would have caught the flipped comparison."""
    null = list(range(10))
    up = L.empirical_p(2, null)
    down = L.empirical_p(2, null, greater=False)
    assert up != down
    assert up == pytest.approx(9 / 11) and down == pytest.approx(4 / 11)


def test_the_default_is_still_the_enrichment_tail():
    """Every existing caller was written against it."""
    null = list(range(10))
    assert L.empirical_p(2, null) == L.empirical_p(2, null, greater=True)


def test_the_written_file_carries_both_tails(tmp_path):
    """And the loader counts two directions because of it.

    Dropping p_matched_less from the row left 175 tests passing: the loader
    had a forward-compatible path for it and nothing checked the writer ever
    produced one.
    """
    from kgav import landscape as LS
    rel = _release(tmp_path)
    out = tmp_path / "landscape.json"
    L.main(["--release", str(rel), "--permutations", "300", "--out", str(out)])
    doc = json.loads(out.read_text())
    assert doc["taxa"], "fixture must measure at least one taxon"
    for row in doc["taxa"].values():
        assert row["p_matched_less"] is not None
        assert row["p_matched"] is not None
    land = LS.load(out, release=rel.name)
    assert land.depletion_tested()
    assert land.n_tests == 2 * len(doc["taxa"])


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
    # Both tails, not just the enrichment one. The table reported "p matched"
    # alone for a long time, and every sub-1 ratio in it was unreadable as a
    # result because nothing had measured the other direction.
    assert "p naive" in out and "p up" in out and "p down" in out
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
    assert "ranks no drug, recommends no compound" in src


def test_it_writes_nothing_unless_asked(tmp_path):
    """The script gained --out; without it the old guarantee still holds.

    Asserted on the directory rather than on the docstring. The sentence this
    replaces ("Read-only. Writes nothing") was true when written and silently
    became false the moment --out existed, which is a doc checking itself
    against its own prose rather than against the behaviour.
    """
    rel = _release(tmp_path)
    before = sorted(p.name for p in tmp_path.iterdir())
    L.main(["--release", str(rel), "--permutations", "50"])
    assert sorted(p.name for p in tmp_path.iterdir()) == before


def test_out_writes_a_file_the_loader_accepts(tmp_path):
    from kgav import landscape as LS
    rel = _release(tmp_path)
    out = tmp_path / "landscape.json"
    L.main(["--release", str(rel), "--permutations", "300", "--out", str(out)])
    assert out.exists()
    land = LS.load(out, release=rel.name)
    assert land.usable, land.refused_reason
    assert land.permutations == 300


def test_a_single_virus_run_marks_its_file_partial(tmp_path):
    """And the loader refuses it, so one taxon cannot speak for the release."""
    from kgav import landscape as LS
    rel = _release(tmp_path)
    out = tmp_path / "landscape.json"
    L.main(["--release", str(rel), "--permutations", "300",
            "--virus", "NCBITaxon:2697049", "--out", str(out)])
    doc = json.loads(out.read_text())
    assert doc["partial"] is True
    assert "NCBITaxon:2697049" in doc["partial_reason"]
    assert not LS.load(out, release=rel.name).usable


def test_the_written_file_is_strict_json(tmp_path):
    """An infinite ratio would make it Python-only. None is stored.

    json.loads with parse_constant rejects Infinity/NaN, which json.loads
    otherwise accepts silently -- so a file that round-trips in Python but is
    not JSON would pass an ordinary read.
    """
    rel = _release(tmp_path)
    out = tmp_path / "landscape.json"
    L.main(["--release", str(rel), "--permutations", "300", "--out", str(out)])

    def reject(const):
        raise ValueError(f"not strict JSON: {const}")

    json.loads(out.read_text(), parse_constant=reject)


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


# ------------------------------- what a new source could actually reach
def test_the_unreachable_set_is_broken_down_by_clinical_phase(tmp_path):
    """"Add another source" is only worth doing if that source covers the
    molecules that are missing. DrugCentral, IUPHAR and DrugBank curate
    approved and clinical-stage drugs, not bulk screening compounds, so
    the phase mix of the untargeted set is the ceiling on what they could
    add."""
    rel = _release(tmp_path)
    nodes = [json.loads(l) for l in
             (rel / "nodes.jsonl").read_text().splitlines()]
    for n in nodes:
        if n["id"] == "INCHIKEY:I0":
            n["properties"]["max_phase"] = 4
        elif n["id"] == "INCHIKEY:I1":
            n["properties"]["max_phase"] = 2
        elif n["id"] == "INCHIKEY:N0":
            n["properties"]["max_phase"] = 0
    (rel / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in nodes))
    c = L.census(_graph(rel))
    ph = c["untargeted_by_phase"]
    assert ph["approved"] == 1
    assert ph["clinical"] == 1
    assert ph["preclinical"] == 1
    assert sum(ph.values()) == c["molecules"] - c["analysable"]


def test_the_ceiling_is_approved_plus_clinical_not_the_whole_gap(
        tmp_path, capsys):
    """Reporting the whole untargeted count as addressable would overstate
    what a curated source can do by an order of magnitude."""
    rel = _release(tmp_path)
    nodes = [json.loads(l) for l in
             (rel / "nodes.jsonl").read_text().splitlines()]
    for n in nodes:
        if n["id"].startswith("INCHIKEY:I"):
            n["properties"]["max_phase"] = 4
    (rel / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in nodes))
    L.main(["--release", str(rel), "--permutations", "50"])
    out = capsys.readouterr().out
    assert "is the CEILING on what one could add here" in out
    assert "bulk screening compounds no such source curates" in out


def test_a_molecule_with_no_phase_is_unknown_not_preclinical(tmp_path):
    """Absent is not zero. ChEMBL leaves max_phase null for most of its
    compounds, and filing those as preclinical would assert something the
    source does not say."""
    c = L.census(_graph(_release(tmp_path)))
    ph = c["untargeted_by_phase"]
    assert ph["unknown"] > 0
    assert ph["preclinical"] == 0


# ------------------------------------------------- the interactome channel
def _ppi_release(tmp_path: Path) -> Path:
    """A virus with physical interactors and NO CRISPR screen.

    The shape the channel exists for: six of nine flaviviruses have no
    screen in this ORCS release, so the dependency channel cannot speak
    about them however good the PPI data is.
    """
    d = tmp_path / "ppi"
    d.mkdir(parents=True, exist_ok=True)
    nodes = [{"id": "NCBITaxon:11053", "class": "OrganismTaxon",
              "properties": {}},
             {"id": "REACT:P1", "class": "Pathway", "properties": {}},
             {"id": "UniProtKB:POLY", "class": "Protein",
              "properties": {"taxon_id": "NCBITaxon:11053",
                             "is_viral": True}},
             # a viral protein of ANOTHER virus, to prove the channel
             # does not pool taxa
             {"id": "UniProtKB:OTHER", "class": "Protein",
              "properties": {"taxon_id": "NCBITaxon:64320",
                             "is_viral": True}},
             {"id": "NCBITaxon:64320", "class": "OrganismTaxon",
              "properties": {}}]
    edges = []
    for i in range(8):
        p = f"UniProtKB:H{i}"
        nodes.append({"id": p, "class": "Protein",
                      "properties": {"taxon_id": "NCBITaxon:9606",
                                     "is_viral": False}})
        edges.append((p, "PARTICIPATES_IN", "REACT:P1", {}))
    # H0..H3 touch the DENV-1 polyprotein; H4 touches the other virus only
    for i in range(4):
        edges.append((f"UniProtKB:H{i}", "PHYSICALLY_INTERACTS_WITH",
                      "UniProtKB:POLY",
                      {"detection_method": "MI:0018" if i == 0 else "MI:1314"}))
    edges.append(("UniProtKB:H4", "PHYSICALLY_INTERACTS_WITH",
                  "UniProtKB:OTHER", {"detection_method": "MI:1314"}))
    # a host-host edge carrying the SAME predicate, from the STRING layer
    edges.append(("UniProtKB:H5", "PHYSICALLY_INTERACTS_WITH",
                  "UniProtKB:H6", {"detection_method": "MI:1314"}))
    # a viral-viral edge
    edges.append(("UniProtKB:POLY", "PHYSICALLY_INTERACTS_WITH",
                  "UniProtKB:OTHER", {"detection_method": "MI:1314"}))
    for i in range(3):
        c = f"INCHIKEY:T{i}"
        nodes.append({"id": c, "class": "SmallMolecule", "properties": {}})
        edges.append((c, "TARGETS", f"UniProtKB:H{i}", {}))
    (d / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in nodes))
    (d / "edges.jsonl").write_text("\n".join(json.dumps(
        {"subject": s, "predicate": p, "object": o, "qualifiers": q,
         "primary_knowledge_source": "infores:virhostnet", "evidence_tier": 1,
         "first_asserted_date": "2021-01-01"}) for s, p, o, q in edges))
    return d


def test_partners_are_the_host_side_of_a_viral_edge(tmp_path):
    g = _graph(_ppi_release(tmp_path))
    partners, hubs = L.interactome_partners(
        g, L.viral_proteins(g)["NCBITaxon:11053"])
    assert partners == {f"UniProtKB:H{i}" for i in range(4)}
    assert hubs == set()


def test_a_host_host_edge_is_not_an_interactome_partner(tmp_path):
    """STRING writes PHYSICALLY_INTERACTS_WITH between two human proteins.
    The predicate is shared across layers, so the viral side has to be
    established from the NODE, never from the edge."""
    g = _graph(_ppi_release(tmp_path))
    partners, _ = L.interactome_partners(
        g, L.viral_proteins(g)["NCBITaxon:11053"])
    assert "UniProtKB:H5" not in partners and "UniProtKB:H6" not in partners


def test_a_viral_viral_edge_contributes_no_partner(tmp_path):
    g = _graph(_ppi_release(tmp_path))
    partners, _ = L.interactome_partners(
        g, L.viral_proteins(g)["NCBITaxon:11053"])
    assert not any(g.is_viral(p) for p in partners)


def test_another_virus_partners_are_not_pooled(tmp_path):
    """Each taxon's interactome is its own. H4 touches Zika's protein and
    must not appear under DENV-1."""
    g = _graph(_ppi_release(tmp_path))
    vp = L.viral_proteins(g)
    denv, _ = L.interactome_partners(g, vp["NCBITaxon:11053"])
    zika, _ = L.interactome_partners(g, vp["NCBITaxon:64320"])
    assert "UniProtKB:H4" in zika and "UniProtKB:H4" not in denv


def test_hubs_are_returned_rather_than_silently_dropped(tmp_path):
    """Ubiquitin and ribosomal proteins are the sticky partners of
    pull-down data. How much of an interactome they account for is a
    property of the assay and belongs in the report."""
    g = _graph(_ppi_release(tmp_path))
    partners, hubs = L.interactome_partners(
        g, L.viral_proteins(g)["NCBITaxon:11053"], hub_cut=1)
    assert hubs, "every partner here has degree > 1"
    assert not (partners & hubs)


def test_the_binary_share_counts_only_viral_host_edges(tmp_path):
    """One of the four DENV-1 edges is two-hybrid; the rest are proximity.
    The host-host and viral-viral edges are not this virus's evidence."""
    g = _graph(_ppi_release(tmp_path))
    binary, total = L.binary_share(
        g, L.viral_proteins(g)["NCBITaxon:11053"])
    assert (binary, total) == (1, 4)


def test_a_virus_with_no_ppi_has_no_partners(tmp_path):
    g = _graph(_release(tmp_path))
    partners, _ = L.interactome_partners(
        g, L.viral_proteins(g).get("NCBITaxon:2697049", set()))
    assert partners == set()


def test_the_two_channels_stay_separate_in_the_output(tmp_path):
    """A loader that knows only the dependency channel must not pick up
    interactome numbers under the keys it already reads. The two answer
    different questions -- what a virus NEEDS versus what it TOUCHES --
    and a mixed table means nothing.
    """
    out = tmp_path / "r.json"
    rel = _ppi_release(tmp_path)
    rc = L.main(["--release", str(rel), "--out", str(out),
                 "--permutations", "50"])
    assert rc in (0, 1)
    doc = json.loads(out.read_text())
    assert "interactome" in doc
    # DENV-1 has partners but no dependency factor, so it appears in one
    # table and not the other -- the whole reason for the second channel.
    assert "NCBITaxon:11053" in doc["interactome"]
    assert "NCBITaxon:11053" not in (doc.get("taxa") or {})
    assert "NCBITaxon:11053" in doc["uncovered"]


def test_the_overlap_statistic_has_one_implementation(tmp_path):
    """Both channels call overlap_row. Two copies of a statistic drift,
    and this project has paid that twice: two trustworthy rules, five
    taxon maps."""
    import ast
    tree = ast.parse((ROOT / "scripts" / "pathway_landscape.py").read_text())
    calls = [n.func.id for n in ast.walk(tree)
             if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)]
    # CALLS, not substrings: a definition, an import and a comment all
    # contain the name, and counting text instead of calls is how the
    # provenance check came to be satisfied by its own comment.
    for fn in ("matched_null", "fisher_one_sided", "empirical_p"):
        assert calls.count(fn) == (2 if fn == "empirical_p" else 1), (
            f"{fn} is called {calls.count(fn)} times; both channels "
            f"should reach it through overlap_row")
