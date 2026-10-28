"""Gate for the build audit.

The facts it joins live in five places -- pipeline.LAYERS, data/raw,
the per-layer releases, the assembled release and the register -- and
answering "can the agent handle DENV-1 yet" meant holding that join in
your head. That is how a virus with a taxon in the map and no file on
disk went unnoticed for weeks.
"""
from __future__ import annotations

import json
import sys
from importlib import import_module
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

A = import_module("audit_build")

CFG = {"viruses": [
    {"taxon": "11053", "name": "DENV-1", "family": "Flaviviridae",
     "mondo": None,
     "chain_source": {"accession": "P27909", "taxon": "11053", "chains": 14,
                      "verified": "2026-09-17"}},
    # NO per-entry family, exactly like the real register: the
    # coronaviruses inherit it from defaults below. The fixture used to
    # declare it, so the audit's family column never took the defaults
    # path and a mutant reverting it survived.
    # YELLOW FEVER'S SHAPE: held-out labels and no feature route at
    # all. The fixture had DENV-1 (labels + PPI) and SARS-CoV-2
    # (routes, no labels) and nothing in between, so the branch that
    # names the unreachable case was never taken and a mutant removing
    # it survived.
    {"taxon": "11089", "name": "Yellow fever", "family": "Flaviviridae",
     "mondo": None,
     "chain_source": {"accession": "P03314", "taxon": "11089", "chains": 13,
                      "verified": "2026-09-17"}},
    {"taxon": "2697049", "name": "SARS-CoV-2",
     "mondo": None,
     "chain_source": {"accession": "P0DTD1", "taxon": "2697049",
                      "chains": 16, "verified": "2026-09-17"}},
], "defaults": {"family": "Coronaviridae"}}


def _release(tmp_path: Path) -> Path:
    d = tmp_path / "v0.1"
    d.mkdir(parents=True)
    nodes = [
        {"id": "NCBITaxon:11053", "class": "OrganismTaxon", "properties": {}},
        {"id": "NCBITaxon:2697049", "class": "OrganismTaxon",
         "properties": {}},
        {"id": "NCBITaxon:11089", "class": "OrganismTaxon",
         "properties": {}},
        {"id": "MONDO:1", "class": "Disease", "properties": {}},
        {"id": "UniProtKB:POLY", "class": "Protein",
         "properties": {"is_viral": True, "taxon_id": "NCBITaxon:11053"}},
        {"id": "UniProtKB:PRO_1", "class": "Protein",
         "properties": {"is_viral": True, "taxon_id": "NCBITaxon:11053"}},
        {"id": "UniProtKB:NSP5", "class": "Protein",
         "properties": {"is_viral": True, "taxon_id": "NCBITaxon:2697049"}},
        {"id": "UniProtKB:H1", "class": "Protein",
         "properties": {"is_viral": False, "taxon_id": "NCBITaxon:9606"}},
        {"id": "UniProtKB:H2", "class": "Protein",
         "properties": {"is_viral": False, "taxon_id": "NCBITaxon:9606"}},
        {"id": "INCHIKEY:D1", "class": "SmallMolecule", "properties": {}},
        {"id": "INCHIKEY:D2", "class": "SmallMolecule", "properties": {}},
    ]
    edges = [
        # DENV-1: two PPI edges, two distinct partners, no CRISPR, no drug
        ("UniProtKB:H1", "PHYSICALLY_INTERACTS_WITH", "UniProtKB:POLY", {}),
        ("UniProtKB:H2", "PHYSICALLY_INTERACTS_WITH", "UniProtKB:PRO_1", {}),
        # a host-host edge, which is NOT virus-host: excluded because
        # neither endpoint resolves to a viral taxon
        ("UniProtKB:H1", "PHYSICALLY_INTERACTS_WITH", "UniProtKB:H2", {}),
        # a VIRAL-VIRAL edge, which is excluded by a DIFFERENT check --
        # both endpoints resolve to a viral taxon, so without the guard
        # each virus would count the other's protein as a host partner.
        # The fixture had none, and a mutant removing that guard
        # survived while the host-host case went on passing.
        ("UniProtKB:POLY", "PHYSICALLY_INTERACTS_WITH", "UniProtKB:NSP5",
         {}),
        # SARS-CoV-2: CRISPR both directions, a drug, a disease
        ("UniProtKB:H1", "HOST_FACTOR_FOR", "NCBITaxon:2697049",
         {"direction": "dependency"}),
        ("UniProtKB:H2", "HOST_FACTOR_FOR", "NCBITaxon:2697049",
         {"direction": "restriction"}),
        ("INCHIKEY:D1", "INHIBITS", "UniProtKB:NSP5", {}),
        # DENV-1 has NO protein-level inhibitor and two organism-level
        # antivirals, which is the real shape of Yellow fever: 1,070
        # ChEMBL compounds and zero INHIBITS, because a cell-based
        # assay has no protein to point at.
        ("INCHIKEY:D1", "HAS_ANTIVIRAL_ACTIVITY_AGAINST", "NCBITaxon:11053",
         {}),
        ("INCHIKEY:D2", "HAS_ANTIVIRAL_ACTIVITY_AGAINST", "NCBITaxon:11053",
         {}),
        ("INCHIKEY:D2", "MEASURED_INACTIVE_AGAINST", "UniProtKB:NSP5", {}),
        # labels only, nothing to reach them with
        ("INCHIKEY:D1", "HAS_ANTIVIRAL_ACTIVITY_AGAINST", "NCBITaxon:11089",
         {}),
        ("NCBITaxon:2697049", "CAUSES", "MONDO:1", {}),
    ]
    (d / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in nodes))
    (d / "edges.jsonl").write_text("\n".join(json.dumps(
        {"subject": s, "predicate": p, "object": o, "qualifiers": q,
         "primary_knowledge_source": "infores:test", "evidence_tier": 1,
         "first_asserted_date": "2021-01-01"}) for s, p, o, q in edges))
    return d


def _run(tmp_path, capsys):
    cfg = tmp_path / "v.yaml"
    import yaml
    cfg.write_text(yaml.safe_dump(CFG))
    rc = A.main(["--config", str(cfg), "--release", str(_release(tmp_path)),
                 "--root", str(tmp_path)])
    return rc, capsys.readouterr().out


def test_it_counts_virus_host_ppi_and_nothing_else(tmp_path, capsys):
    """The predicate is shared with the STRING layer, so the fixture
    carries all three shapes: virus-host, host-host and viral-viral.
    Only the first is this layer's evidence, and the other two are
    excluded by different checks."""
    _rc, out = _run(tmp_path, capsys)
    denv = [ln for ln in out.splitlines() if "DENV-1" in ln and "Flav" in ln]
    assert denv, out
    # proteins 2, chains 1, vhPPI 2, partners 2
    assert denv[0].split()[2:6] == ["2", "1", "2", "2"]


def test_a_virus_with_no_data_on_an_axis_is_named(tmp_path, capsys):
    """A zero that appears nowhere is the failure this whole audit is
    for: DENV-1 has PPI and no CRISPR, SARS-CoV-2 the reverse."""
    _rc, out = _run(tmp_path, capsys)
    assert "no CRISPR dependency factor" in out
    assert "no virus-host PPI" in out
    gap = out.split("5. WHAT IS MISSING")[1]
    assert "DENV-1" in gap and "SARS-CoV-2" in gap


def test_dependency_and_restriction_are_not_pooled(tmp_path, capsys):
    """Opposite claims. One of each here, and they must not add to two
    in a single column."""
    _rc, out = _run(tmp_path, capsys)
    row = [ln for ln in out.splitlines() if "SARS-CoV-2" in ln and "Coro" in ln]
    assert row[0].split()[6:8] == ["1", "1"]


def test_an_unassembled_release_says_so_and_stops(tmp_path, capsys):
    """Rather than printing an empty per-virus table that reads as 'no
    data for any virus'."""
    import yaml
    cfg = tmp_path / "v.yaml"
    cfg.write_text(yaml.safe_dump(CFG))
    rc = A.main(["--config", str(cfg), "--release", str(tmp_path / "nope"),
                 "--root", str(tmp_path)])
    out = capsys.readouterr().out
    assert rc == 1 and "NOT ASSEMBLED" in out
    assert "PER-VIRUS" not in out


def test_a_layer_that_built_nothing_is_distinguished_from_one_not_built(
        tmp_path, capsys):
    """Identical in a directory listing, and they need different work."""
    d = tmp_path / "data" / "releases" / "v0.1-vhppi"
    d.mkdir(parents=True)
    (d / "nodes.jsonl").write_text("")
    (d / "edges.jsonl").write_text("")
    _rc, out = _run(tmp_path, capsys)
    assert "ran and joined nothing" in out
    assert "NOT BUILT" in out


def test_a_blocked_layer_names_what_is_absent(tmp_path, capsys):
    """"skipped" sends the reader to look for a bug in a layer that
    simply had no input."""
    _rc, out = _run(tmp_path, capsys)
    assert "BLOCKED:" in out
    assert "nothing matches" in out or "no data/raw" in out


def test_a_viral_viral_edge_is_not_a_host_partner(tmp_path, capsys):
    """POLY touches NSP5 in the fixture. Counting it would give DENV-1
    a host partner that is a SARS-CoV-2 protein."""
    _rc, out = _run(tmp_path, capsys)
    for name in ("DENV-1", "SARS-CoV-2"):
        row = next(ln for ln in out.splitlines()
                   if name in ln and ("Flav" in ln or "Coro" in ln))
        cols = row.split()
        vhppi, partners = int(cols[4]), int(cols[5])
        assert vhppi == (2 if name == "DENV-1" else 0), row
        assert partners == (2 if name == "DENV-1" else 0), row


def test_the_family_column_applies_the_registers_defaults(tmp_path, capsys):
    """Coronaviruses carry no per-entry family. Reading the key alone
    printed "?" for all seven of them."""
    _rc, out = _run(tmp_path, capsys)
    row = next(ln for ln in out.splitlines() if "SARS-CoV-2" in ln
               and "prot" not in ln)
    assert "Coro" in row, row
    assert "?" not in row.split()[1]


def test_organism_level_antivirals_are_counted(tmp_path, capsys):
    """A cell-based assay measuring "compound inhibits yellow fever
    virus" has no protein to point at, so chembl.py lands it on the
    TAXON. Counting only INHIBITS reported Yellow fever as having no
    compound while ChEMBL holds 1,070 of them."""
    _rc, out = _run(tmp_path, capsys)
    row = next(ln for ln in out.splitlines() if "DENV-1" in ln and "Flav" in ln)
    cols = row.split()
    assert cols[8] == "0", f"no protein-level inhibitor: {row}"
    assert cols[9] == "2", f"two organism-level antivirals: {row}"


def test_the_two_kinds_of_absence_are_not_the_same_gap(tmp_path, capsys):
    """"No compound of any kind" and "organism-level only" need
    different work: the first is a source limit, the second is a
    resolution limit on data that is present."""
    _rc, out = _run(tmp_path, capsys)
    gaps = out.split("5. WHAT IS MISSING")[1]
    assert "no protein-level inhibitor (organism-level only)" in gaps
    assert "DENV-1" in gaps.split("organism-level only")[1].split("\n")[1]


def test_an_inactive_measurement_is_reported_separately(tmp_path, capsys):
    """Measured and found inactive is evidence, not absence, and must
    not read as either a hit or a gap."""
    _rc, out = _run(tmp_path, capsys)
    row = next(ln for ln in out.splitlines()
               if "SARS-CoV-2" in ln and "Coro" in ln)
    assert row.split()[10] == "1", row


def test_labels_without_a_feature_route_are_named(tmp_path, capsys):
    """Yellow fever holds 1,030 held-out antiviral labels and has no
    PPI, no CRISPR and no INHIBITS. Nothing can reach a compound for
    it, so every prediction is a miss by construction -- and in the
    per-virus table that row reads like its best-covered column.
    """
    _rc, out = _run(tmp_path, capsys)
    # DENV-1 in the fixture: two antiviral labels, and PPI as a route
    assert "LABELS vs FEATURES" in out
    assert "evaluable" in out


def test_a_virus_with_features_and_no_labels_is_distinguished(tmp_path,
                                                               capsys):
    """SARS-CoV-2 in the fixture has CRISPR and an inhibitor but no
    antiviral label: it cannot be scored, which is the opposite
    problem and needs the opposite work."""
    _rc, out = _run(tmp_path, capsys)
    block = out.split("LABELS vs FEATURES")[1]
    row = next(ln for ln in block.splitlines() if "SARS-CoV-2" in ln)
    assert "features, nothing to score" in row, row


def test_the_held_out_set_comes_from_the_schema(tmp_path, capsys):
    """Not a list in this script. The schema already marks which
    predicates are labels, and a second copy would drift from it."""
    _rc, out = _run(tmp_path, capsys)
    head = out.split("LABELS vs FEATURES")[1].split("`labels` counts")[0]
    for pred in ("HAS_ANTIVIRAL_ACTIVITY_AGAINST", "TREATS",
                 "MEASURED_INACTIVE_AGAINST"):
        assert pred in head, head


def test_the_negative_label_is_not_presented_as_evidence():
    """MEASURED_INACTIVE_AGAINST is held out too, so the `inact`
    column is answer sheet rather than feature. The legend called it
    "evidence, not absence", which is true and reads as a feature."""
    src = (ROOT / "scripts" / "audit_build.py").read_text()
    assert "negative label" in src


def test_a_feature_route_is_shown_with_its_size(tmp_path, capsys):
    """JEV reads `evaluable` on seven PPI edges, which is true and
    reads identically to eight thousand. No threshold is invented:
    the size is shown and the reader judges."""
    _rc, out = _run(tmp_path, capsys)
    block = out.split("LABELS vs FEATURES")[1]
    row = next(ln for ln in block.splitlines() if "DENV-1" in ln)
    assert "ppi(2)" in row, row


def test_a_virus_with_labels_and_no_route_is_called_out(tmp_path, capsys):
    """Yellow fever: 1,030 held-out antiviral labels, no PPI, no
    CRISPR, no INHIBITS. No metapath terminates there, so every
    prediction is a miss by construction -- and in the per-virus table
    that row reads like its best-covered column."""
    _rc, out = _run(tmp_path, capsys)
    block = out.split("LABELS vs FEATURES")[1]
    row = next(ln for ln in block.splitlines() if "Yellow fever" in ln)
    assert "LABELS, NO FEATURE ROUTE" in row, row
    assert "miss by" in block
    assert "measures the graph's reach, not the model" in block


def test_the_legend_is_not_split_by_the_verdict_section(tmp_path, capsys):
    """Every legend line belongs together. The verdict block was
    inserted between `inact` and `dis`, so the last legend entry
    printed below a whole section and read as part of it."""
    _rc, out = _run(tmp_path, capsys)
    legend_end = out.index("dis         CAUSES edges to a Disease node")
    section = out.index("LABELS vs FEATURES")
    assert legend_end < section, "the legend must finish before the section"


# ------------- a route that exists and reaches none of the held-out labels
def _comparison(dirpath, release, pos_by_virus, stale=False):
    """A metapath_comparison.json stamped against `release`.

    Written through provenance.write_results so the stamp is the one
    check_file reads, rather than a hand-made dict that happens to have
    the right keys today.
    """
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from kgav.provenance import write_results
    dirpath.mkdir(parents=True, exist_ok=True)
    per = {v: {"COMBINED": {"pos": n, "pool": max(n, 1)}}
           for v, n in pos_by_virus.items()}
    write_results(dirpath / "metapath_comparison.json", {"per_virus": per},
                  release=release)
    if stale:
        doc = json.loads((dirpath / "metapath_comparison.json").read_text())
        doc["provenance"]["release"]["edges_sha"] = "0" * 16
        (dirpath / "metapath_comparison.json").write_text(json.dumps(doc))
    return dirpath


def _run_with_reach(tmp_path, capsys, pos_by_virus, stale=False):
    import yaml
    cfg = tmp_path / "v.yaml"
    cfg.write_text(yaml.safe_dump(CFG))
    rel = _release(tmp_path)
    res = _comparison(tmp_path / "results", rel, pos_by_virus, stale=stale)
    rc = A.main(["--config", str(cfg), "--release", str(rel),
                 "--root", str(tmp_path), "--results", str(res)])
    return rc, capsys.readouterr().out


def test_a_route_reaching_no_label_is_not_called_evaluable(tmp_path, capsys):
    """DENV-1 on the real release: 349 held-out antiviral labels, 53
    virus-host interactions, and the scorer reaches NONE of the 349.
    This column read `evaluable`, because edge counts cannot see
    reachability. Seven of fifteen viruses are in that state.
    """
    _rc, out = _run_with_reach(tmp_path, capsys, {"DENV-1": 0})
    block = out.split("LABELS vs FEATURES")[1]
    row = next(ln for ln in block.splitlines() if "DENV-1" in ln)
    assert "ROUTES REACH NO LABEL" in row, row
    assert "evaluable" not in row, row
    # Whitespace-normalised: the explanation is a wrapped indented
    # paragraph, so the phrase spans a newline and a raw search misses it.
    assert "nothing in this graph can mark it right or wrong" in \
        " ".join(block.split())


def test_a_route_reaching_labels_is_still_evaluable(tmp_path, capsys):
    """The other side, so the verdict is not simply renamed."""
    _rc, out = _run_with_reach(tmp_path, capsys, {"DENV-1": 11})
    row = next(ln for ln in out.split("LABELS vs FEATURES")[1].splitlines()
               if "DENV-1" in ln)
    assert "evaluable" in row and "NO LABEL" not in row, row


def test_an_absent_comparison_is_unmeasured_not_zero(tmp_path, capsys):
    """"Cannot be known" and "reaches nothing" are opposite claims. With
    no comparison file the audit must say the first; returning 0 would
    print the second about every virus in the release.

    --results POINTS AT AN EMPTY DIRECTORY, deliberately. This used
    _run(), which leaves --results at its default -- the repository's
    real data/results -- so the file existed, the absent-file branch
    never ran, and the test passed on a stale comparison while its name
    claimed an absent one. Mutating that branch to return {} survived.
    """
    import yaml
    cfg = tmp_path / "v.yaml"
    cfg.write_text(yaml.safe_dump(CFG))
    empty = tmp_path / "no-results"
    empty.mkdir()
    rc = A.main(["--config", str(cfg), "--release", str(_release(tmp_path)),
                 "--root", str(tmp_path), "--results", str(empty)])
    out = capsys.readouterr().out
    assert rc is not None
    block = out.split("LABELS vs FEATURES")[1]
    row = next(ln for ln in block.splitlines() if "DENV-1" in ln)
    assert "reach unmeasured" in row, row
    assert "ROUTES REACH NO LABEL" not in block
    assert "absent, stale or" in " ".join(block.split())


def test_a_stale_comparison_is_not_quoted(tmp_path, capsys):
    """It is the same check `make check-results` runs. A comparison
    measured on a graph that has since moved describes that graph, and
    this script has twice shipped a figure from a previous run."""
    _rc, out = _run_with_reach(tmp_path, capsys, {"DENV-1": 0}, stale=True)
    block = out.split("LABELS vs FEATURES")[1]
    assert "ROUTES REACH NO LABEL" not in block, (
        "a stale comparison was quoted as a measurement")
    assert "reach unmeasured" in block


def _verdict_sections(block: str) -> dict[str, str]:
    """Each verdict heading mapped to the text until the next heading.

    The point is attachment. Asserting a paragraph is somewhere in the
    section passes while it sits under the wrong list, which is exactly
    what happened: "No metapath terminates at these" is true of DENV-4
    and Yellow fever and printed under `features, nothing to score`,
    then under ROUTES REACH NO LABEL where it is false.
    """
    heads = [v for v, _ in A.VERDICT_NOTES]
    out, current = {}, None
    for line in block.splitlines():
        hit = next((h for h in heads if line.strip().startswith(h + ":")),
                   None)
        if hit:
            current = hit
            out[hit] = ""
        elif current:
            out[current] += line + "\n"
    return out


def test_each_explanation_sits_under_its_own_verdict(tmp_path, capsys):
    """A route that reaches no label is NOT a virus no metapath reaches,
    and the two paragraphs say opposite things about whether a path
    exists. Printed after the loop they stacked under whichever list
    came last."""
    _rc, out = _run_with_reach(tmp_path, capsys, {"DENV-1": 0})
    secs = _verdict_sections(out.split("LABELS vs FEATURES")[1])
    assert "ROUTES REACH NO LABEL" in secs, secs.keys()
    mine = " ".join(secs["ROUTES REACH NO LABEL"].split())
    assert "A metapath does terminate at these" in mine, mine
    assert "No metapath terminates at these" not in mine, (
        "the no-route explanation is attached to the reaches-no-label "
        "list, where it is false")
    if "LABELS, NO FEATURE ROUTE" in secs:
        theirs = " ".join(secs["LABELS, NO FEATURE ROUTE"].split())
        assert "No metapath terminates at these" in theirs, theirs
        assert "reaches none of the held-out" not in theirs


def test_every_called_out_verdict_has_a_paired_note():
    """VERDICT_NOTES is the only list of these headings. A tier added to
    the verdict logic and not to it prints no list at all, which is a
    silent omission rather than a visible one."""
    src = (ROOT / "scripts" / "audit_build.py").read_text()
    for verdict, _note in A.VERDICT_NOTES:
        assert f'verdict = "{verdict}"' in src, (
            f"{verdict!r} is in VERDICT_NOTES and no branch assigns it")
    for line in src.splitlines():
        if line.strip().startswith('verdict = "'):
            name = line.split('"')[1]
            if name in ("evaluable", "neither",
                        "evaluable (reach unmeasured)"):
                continue
            assert any(name == v for v, _ in A.VERDICT_NOTES), (
                f"{name!r} is assigned as a verdict and is not in "
                f"VERDICT_NOTES, so its list is never printed")
