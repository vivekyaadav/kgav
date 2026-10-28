"""Day 5 gate: the bridge layer joins, and alias resolution is what makes it."""
import json
from pathlib import Path

import pytest
import yaml

from kgav.aliases import AliasIndex
from kgav.emit import Emit
from kgav.schema import load_schema
from kgav.vhppi import classify_row, ingest_mitab, parse_identifier, parse_row

ROOT = Path(__file__).resolve().parents[1]
VHPPI = ROOT / "data" / "raw" / "vhppi"
SPINE = ROOT / "data" / "releases" / "v0.1-spine"
needs_data = pytest.mark.skipif(not (VHPPI.exists() and SPINE.exists()),
                                reason="vhppi sources or spine release missing")

NSP5 = "UniProtKB:PRO_0000449623"
NSP5_PP1A = "UniProtKB:PRO_0000449639"      # pp1a id for the same protein
ORF7A = "UniProtKB:P0DTC7"
ORF7A_CHAIN = "UniProtKB:PRO_0000449654"
ACE2, TMPRSS2 = "UniProtKB:Q9BYF1", "UniProtKB:O15393"


def _row(a, b, method, pmid, ta, tb, itype, score):
    f = ["-"] * 15
    f[0], f[1] = a, b
    f[6] = f'psi-mi:"{method}"(x)'
    f[8] = f"pubmed:{pmid}"
    f[9], f[10] = f"taxid:{ta}", f"taxid:{tb}"
    f[11] = f'psi-mi:"{itype}"(y)'
    f[14] = f"virhostnet-miscore:{score}"
    return "\t".join(f)


@pytest.fixture
def world(tmp_path):
    def prot(pid, viral, **extra):
        p = {"taxon_id": "NCBITaxon:2697049" if viral else "NCBITaxon:9606",
             "is_viral": viral, "sequence_hash": pid[-4:], "reviewed": True}
        p.update(extra)
        return {"id": pid, "class": "Protein", "properties": p}

    spine_nodes = [
        {"id": "NCBITaxon:2697049", "class": "OrganismTaxon",
         "properties": {"label": "SARS-CoV-2", "family": "Coronaviridae",
                        "baltimore_class": "IV", "is_enveloped": True}},
        prot(NSP5, True, mature_peptide="3C-like proteinase", protein_family="nsp5"),
        prot(NSP5_PP1A, True), prot(ORF7A, True), prot(ORF7A_CHAIN, True),
    ]
    same_as = lambda s, o, rule: {
        "subject": s, "predicate": "SAME_AS", "object": o,
        "qualifiers": {"merge_rule": rule},
        "primary_knowledge_source": "infores:uniprot",
        "evidence_tier": 1, "first_asserted_date": "2020-04-22"}
    spine_edges = [same_as(NSP5_PP1A, NSP5, "shared_chain_between_polyproteins"),
                   same_as(ORF7A_CHAIN, ORF7A, "single_chain_equals_parent_protein")]
    host_nodes = [prot(ACE2, False), prot(TMPRSS2, False)]

    for name, nodes, edges in [("spine", spine_nodes, spine_edges),
                               ("host", host_nodes, [])]:
        d = tmp_path / name
        d.mkdir()
        (d / "nodes.jsonl").write_text("\n".join(json.dumps(x) for x in nodes))
        (d / "edges.jsonl").write_text("\n".join(json.dumps(x) for x in edges))

    rows = [
        _row(f"uniprotkb:P0DTD1-{NSP5_PP1A.split(':')[1]}", f"uniprotkb:{ACE2.split(':')[1]}",
             "MI:1314", "32838362", 2697049, 9606, "MI:0915", "0.51"),
        _row(f"uniprotkb:P0DTC7-{ORF7A_CHAIN.split(':')[1]}", f"uniprotkb:{TMPRSS2.split(':')[1]}",
             "MI:0018", "33110215", 2697049, 9606, "MI:0407", "0.72"),
        _row(f"uniprotkb:P0DTD1-{NSP5_PP1A.split(':')[1]}", f"uniprotkb:{ACE2.split(':')[1]}",
             "MI:1314", "32838362", 2697049, 9606, "MI:0915", "0.51"),      # duplicate
        _row(f"uniprotkb:{ACE2.split(':')[1]}", f"uniprotkb:{TMPRSS2.split(':')[1]}",
             "MI:0400", "1234", 9606, 9606, "MI:0915", "0.3"),              # host-host
        _row(f"uniprotkb:{NSP5.split(':')[1]}", "uniprotkb:P0DTC7",
             "MI:0400", "1235", 2697049, 2697049, "MI:0915", "0.3"),        # viral-viral
        _row("uniprotkb:NOPE9", f"uniprotkb:{ACE2.split(':')[1]}",
             "MI:0400", "1236", 2697049, 9606, "MI:0915", "0.3"),           # unresolved
        _row(f"uniprotkb:P0DTD1-{NSP5_PP1A.split(':')[1]}", f"uniprotkb:{ACE2.split(':')[1]}",
             "MI:0400", "1237", 99999, 9606, "MI:0915", "0.4"),             # unmapped taxon
    ]
    (tmp_path / "vhn_test.tab27").write_text("\n".join(rows))
    (tmp_path / "viruses.yaml").write_text(
        yaml.safe_dump({"viruses": [{"taxon": "2697049", "name": "SARS-CoV-2"}]}))
    return tmp_path


@pytest.fixture
def built(world):
    index = AliasIndex.from_releases(world / "spine", world / "host")
    em = Emit()
    stats = ingest_mitab(em, [world / "vhn_test.tab27"], index,
                         {"2697049": "2697049"}, "infores:virhostnet", "2024-01-01")
    return em, stats, index


# ------------------------------------------------------------------- parsing
def test_chain_reference_resolves_to_the_chain_not_the_parent():
    """The source is saying WHICH mature peptide it measured. Collapsing to the
    polyprotein throws that away and loses all target specificity in M3."""
    assert parse_identifier("uniprotkb:P0DTD1-PRO_0000449624") == "UniProtKB:PRO_0000449624"
    assert parse_identifier("uniprotkb:P0DTD1") == "UniProtKB:P0DTD1"
    assert parse_identifier("-") is None


def test_row_parsing_extracts_psi_mi_codes():
    r = parse_row(_row("uniprotkb:A", "uniprotkb:B", "MI:0018", "999",
                       2697049, 9606, "MI:0407", "0.9"))
    assert r["detection_method"] == "MI:0018"
    assert r["interaction_type"] == "MI:0407"
    assert r["pmid"] == "999" and r["score"] == "0.9"


# ------------------------------------------------------------------- aliases
def test_alias_index_follows_same_as(world):
    ix = AliasIndex.from_releases(world / "spine", world / "host")
    assert ix.resolve(NSP5_PP1A) == NSP5
    assert ix.resolve(ORF7A_CHAIN) == ORF7A
    assert ix.resolve(NSP5) == NSP5
    assert ix.resolve("UniProtKB:NOPE") is None


def test_alias_cycle_returns_none(tmp_path):
    """A cycle is a data bug. Looping forever would hide it."""
    d = tmp_path / "r"
    d.mkdir()
    nodes = [{"id": f"X:{i}", "class": "Protein",
              "properties": {"taxon_id": "NCBITaxon:9606", "is_viral": False,
                             "sequence_hash": "z", "reviewed": True}} for i in (1, 2)]
    edges = [{"subject": "X:1", "predicate": "SAME_AS", "object": "X:2",
              "qualifiers": {"merge_rule": "r"}, "primary_knowledge_source": "s",
              "evidence_tier": 1, "first_asserted_date": "2020-01-01"},
             {"subject": "X:2", "predicate": "SAME_AS", "object": "X:1",
              "qualifiers": {"merge_rule": "r"}, "primary_knowledge_source": "s",
              "evidence_tier": 1, "first_asserted_date": "2020-01-01"}]
    (d / "nodes.jsonl").write_text("\n".join(json.dumps(x) for x in nodes))
    (d / "edges.jsonl").write_text("\n".join(json.dumps(x) for x in edges))
    assert AliasIndex.from_releases(d).resolve("X:1") is None


# -------------------------------------------------------------------- bridge
def test_aliased_endpoints_join(built):
    """The whole point: without alias resolution both of these rows drop."""
    em, stats, _ = built
    assert stats["edges"] == 2
    objects = {e["object"] for e in em.edges}
    assert objects == {NSP5, ORF7A}


def test_edges_are_host_to_viral(built):
    em, _, index = built
    for e in em.edges:
        assert index.is_viral(e["subject"]) is False
        assert index.is_viral(e["object"]) is True


def test_same_species_pairs_are_excluded(built):
    """host-host belongs to the STRING pass; viral-viral is out of scope."""
    _, stats, _ = built
    assert stats["host_host"] == 1 and stats["viral_viral"] == 1


def test_duplicates_unresolved_and_unmapped_taxa_are_counted(built):
    _, stats, _ = built
    assert stats["duplicate"] == 1
    assert stats["unresolved_id"] == 1
    assert stats["unmapped_viral_taxon"] == 1


def test_assay_distinction_is_preserved_not_averaged(built):
    """Proximity and binary evidence must stay distinguishable at query time."""
    em, stats, _ = built
    assert stats["direct"] == 1 and stats["co_complex"] == 1
    # interaction_type is deliberately absent: VirHostNet reports MI:0915 for
    # every row regardless of assay, so storing it would fabricate a
    # distinction the source does not make.
    assert all("interaction_type" not in e["qualifiers"] for e in em.edges)
    methods = {e["qualifiers"]["detection_method"] for e in em.edges}
    assert methods == {"MI:1314", "MI:0018"}


def test_all_bridge_edges_are_tier_one(built):
    """Tier means how the fact was established, not how good the assay was.
    BioID is curated experimental; demoting it would corrupt the Day 12
    tier ablation."""
    em, _, _ = built
    assert {e["evidence_tier"] for e in em.edges} == {1}


def test_publications_and_score_are_carried(built):
    em, _, _ = built
    for e in em.edges:
        assert e["publications"] and e["publications"][0].startswith("PMID:")
        assert isinstance(e["source_score"], float)


def test_output_validates_against_schema(built):
    em, _, index = built
    violations = load_schema().validate_batch(index.nodes.values(), em.edges)
    assert violations == [], [str(v) for v in violations[:8]]


# ----------------------------------------------------------------- real data
@needs_data
def test_real_bridge_resolves_nearly_everything(tmp_path):
    index = AliasIndex.from_releases(SPINE, ROOT / "data" / "releases" / "v0.1-host")
    cfg = yaml.safe_load((ROOT / "config" / "viruses.yaml").read_text())
    tmap = {}
    for v in cfg["viruses"]:
        tmap[str(v["taxon"])] = str(v["taxon"])
        for iso in v.get("isolate_taxa") or []:
            tmap[str(iso)] = str(v["taxon"])

    em = Emit()
    stats = ingest_mitab(em, sorted(VHPPI.glob("vhn_*.tab27")), index, tmap,
                         "infores:virhostnet", "2024-01-01")
    assert stats["edges"] > 3_000, f"only {stats['edges']} bridge edges"
    # Unresolved ids are the silent-failure mode this layer exists to avoid.
    assert stats["unresolved_id"] / max(stats["rows"], 1) < 0.25


# --------------------------------------------- one gate order, one implementation
def test_the_classifier_and_the_ingest_count_the_same_gates(world):
    """A ROW IS COUNTED AT THE FIRST GATE IT FAILS, which is why `which gate`
    cannot be re-derived by re-testing one condition: a row that would also
    fail `unresolved_id` may have died at `unparseable` first and been
    counted there. diagnose_flavivirus_join.py carries a paragraph about
    this, written after an earlier version asserted rows were lost at a gate
    they never reached.

    So the order lives in classify_row and ingest_mitab consumes it. This
    pins the two together on a fixture that exercises every gate.
    """
    from collections import Counter

    from kgav.vhppi import classify_row

    index = AliasIndex.from_releases(world / "spine", world / "host")
    taxon_map = {"2697049": "2697049"}
    lines = [x for x in (world / "vhn_test.tab27").read_text().splitlines()
             if x.strip()]

    em = Emit()
    stats = ingest_mitab(em, [world / "vhn_test.tab27"], index, taxon_map,
                         "infores:virhostnet", "2024-01-01")

    by_gate = Counter(classify_row(x, index, taxon_map).gate for x in lines)
    assert by_gate["unparseable"] == 0, "fixture sanity: none are unparseable"
    for gate in ("no_detection_method", "unresolved_id", "self_interaction",
                 "not_a_protein", "host_host", "viral_viral",
                 "unmapped_viral_taxon"):
        assert stats[gate] == by_gate[gate], (
            f"{gate}: ingest counted {stats[gate]}, classifier {by_gate[gate]}")

    # Deduplication is deliberately NOT a gate -- it needs state across rows
    # -- so survivors exceed edges by exactly the duplicates.
    assert by_gate[None] == stats["edges"] + stats["duplicate"]


def test_an_unresolved_row_records_which_side_failed_and_its_taxid(world):
    """The counter alone cannot answer the question the counter raises. An
    unresolved HOST id is interactome this build is missing; an unresolved
    VIRAL id is a virus it does not carry. 513 of either reads identically
    as `unresolved_id 513`.
    """
    from kgav.vhppi import classify_row

    index = AliasIndex.from_releases(world / "spine", world / "host")
    line = _row("uniprotkb:NOPE9", f"uniprotkb:{ACE2.split(':')[1]}",
                "MI:0400", "1236", 2697049, 9606, "MI:0915", "0.3")

    v = classify_row(line, index, {"2697049": "2697049"})
    assert v.gate == "unresolved_id"
    assert v.unresolved == (("UniProtKB:NOPE9", "2697049", "a"),), \
        "only the failing side, with its taxid column verbatim"


# ------------------------------------- a virus with no edges is a finding
import sys as _sys

_sys.path.insert(0, str(ROOT / "scripts"))
from ingest_vhppi import coverage_lines

_NODES = {"NCBITaxon:2697049": {"properties": {"label": "SARS-CoV-2"}},
          "NCBITaxon:64320": {"properties": {"label": "Zika"}},
          "NCBITaxon:11060": {"properties": {"label": "DENV-2"}}}
_MAP = {"2697049": "2697049", "64320": "64320", "11060": "11060",
        "999999": "999999"}


def test_a_virus_with_no_edges_is_named_not_omitted():
    """The report was per_virus.most_common(), so a zero did not appear.

    Before 2026-10-27 no flavivirus file had been fetched for VirHostNet
    -- the vhn_*.tab27 files were one coronavirus each -- and Zika sat in
    the taxon map producing nothing while the build said so nowhere.
    """
    from collections import Counter
    lines = coverage_lines(Counter({"NCBITaxon:2697049": 12}), _MAP, _NODES)
    text = "\n".join(lines)
    assert "Zika" in text and "DENV-2" in text
    assert "have NO virus-host PPI edge" in text


def test_it_counts_how_many_viruses_have_nothing():
    from collections import Counter
    text = "\n".join(coverage_lines(
        Counter({"NCBITaxon:2697049": 12}), _MAP, _NODES))
    assert "2 of 3 built viruses have NO virus-host PPI edge" in text


def test_a_zero_with_no_rows_fetched_is_called_a_fetch_gap():
    """A zero reads as "the source does not cover this virus" when the
    truth is that no file for it was ever downloaded. Those need different
    work, and the report must not let them look the same."""
    from collections import Counter
    text = "\n".join(coverage_lines(Counter(), _MAP, _NODES, fetched={}))
    assert "nothing fetched (3)" in text
    assert "FETCH gap" in text
    assert "JOIN gap" not in text


def test_a_zero_with_rows_fetched_is_called_a_join_gap_and_names_the_gate():
    """The failure this was written for. 56 DENV-1 rows were fetched and
    produced no edge; the report still said the cause was that no file had
    been fetched, because that sentence was hardcoded. A virus whose rows
    are on disk has a different problem and a different fix.
    """
    from collections import Counter
    fetched = {"11060": Counter({"rows": 56, "unresolved_id": 56})}
    text = "\n".join(coverage_lines(Counter(), _MAP, _NODES, fetched))
    assert "no edge from 56 fetched row(s)" in text
    assert "unresolved_id 56" in text
    assert "JOIN gap" in text
    assert "Fetching again changes nothing" in text


def test_the_two_kinds_of_zero_are_reported_separately():
    """Both causes in one run, each counted against its own viruses."""
    from collections import Counter
    fetched = {"11060": Counter({"rows": 9, "viral_viral": 9})}
    text = "\n".join(coverage_lines(Counter(), _MAP, _NODES, fetched))
    assert "nothing fetched (2)" in text          # SARS-CoV-2, Zika
    assert "fetched but nothing survived (1)" in text
    assert "3 of 3 built viruses" in text


def test_no_cause_is_claimed_when_there_is_no_tally_to_derive_it_from():
    """Older callers pass no per-file tally. The report must then say it
    does not know, rather than fall back on the sentence that was wrong."""
    from collections import Counter
    text = "\n".join(coverage_lines(Counter(), _MAP, _NODES))
    assert "cannot say whether they were fetched" in text
    assert "FETCH gap" not in text and "JOIN gap" not in text


def test_the_old_hardcoded_cause_is_gone():
    """It survived into a run where every one of the five viruses it
    blamed had a file on disk."""
    from collections import Counter
    for f in (None, {}, {"11060": Counter({"rows": 1, "duplicate": 1})}):
        text = "\n".join(coverage_lines(Counter(), _MAP, _NODES, f))
        assert "one coronavirus each" not in text


def test_a_declared_but_unbuilt_taxon_is_not_counted_as_missing_data():
    """999999 is in the register and not in this release. Counting it
    among the viruses with no PPI would report a scope decision as a data
    gap."""
    from collections import Counter
    text = "\n".join(coverage_lines(
        Counter({"NCBITaxon:2697049": 1}), _MAP, _NODES))
    assert "2 of 3 built viruses" in text
    assert "1 further declared taxa are not in this release" in text


def test_full_coverage_says_nothing_about_gaps():
    """The branch must stay quiet when there is nothing to report."""
    from collections import Counter
    full = Counter({"NCBITaxon:2697049": 1, "NCBITaxon:64320": 2,
                    "NCBITaxon:11060": 3})
    text = "\n".join(coverage_lines(full, _MAP, _NODES))
    assert "NO virus-host PPI" not in text
    assert "have NO virus-host PPI edge" not in text


# ------------------------------------- attributing a file back to a virus
from ingest_vhppi import rows_per_taxon

_VIRUSES = [{"name": "DENV-1", "taxon": "11053"},
            {"name": "SARS-CoV-2", "taxon": "2697049"},
            {"name": "Yellow fever", "taxon": "11089"}]


def test_the_sidecar_says_which_virus_a_file_belongs_to(tmp_path):
    """fetch_vhppi records the virus it asked for. That is better evidence
    than the filename because it survives a rename."""
    from collections import Counter
    f = tmp_path / "whatever.tab27"
    f.write_text("")
    (tmp_path / "whatever.tab27.source.json").write_text(
        '{"virus": "DENV-1", "taxon": "11053"}')
    got = rows_per_taxon({f: Counter({"rows": 56})}, _VIRUSES, {})
    assert got["11053"]["rows"] == 56


def test_the_filename_is_the_fallback_for_a_file_with_no_sidecar(tmp_path):
    """The coronavirus files were placed by hand before the fetcher existed
    and carry no sidecar. They still have to be attributable."""
    from collections import Counter
    f = tmp_path / "vhn_sarscov2.tab27"
    f.write_text("")
    got = rows_per_taxon({f: Counter({"rows": 7})}, _VIRUSES, {})
    assert got["2697049"]["rows"] == 7


def test_a_stem_with_a_space_in_the_register_name_still_matches(tmp_path):
    """vhn_yellowfever.tab27 <- "Yellow fever". The stem rule lives in one
    place precisely so this cannot drift."""
    from collections import Counter
    f = tmp_path / "vhn_yellowfever.tab27"
    f.write_text("")
    got = rows_per_taxon({f: Counter({"rows": 4})}, _VIRUSES, {})
    assert got["11089"]["rows"] == 4


def test_an_isolate_taxon_in_a_sidecar_folds_to_the_species(tmp_path):
    """The map is what the rest of the report counts in."""
    from collections import Counter
    f = tmp_path / "x.tab27"
    f.write_text("")
    (tmp_path / "x.tab27.source.json").write_text('{"taxon": "31634"}')
    got = rows_per_taxon({f: Counter({"rows": 3})}, _VIRUSES,
                         {"31634": "11060"})
    assert got["11060"]["rows"] == 3


def test_a_multi_virus_file_is_unattributed_rather_than_guessed(tmp_path):
    """intact_*.tab27 holds many viruses. Spreading its rows over them, or
    pinning them on one, would invent the evidence the report then uses to
    name a cause."""
    from collections import Counter
    f = tmp_path / "intact_viral.tab27"
    f.write_text("")
    got = rows_per_taxon({f: Counter({"rows": 900})}, _VIRUSES, {})
    assert list(got) == ["?"] and got["?"]["rows"] == 900


def test_a_corrupt_sidecar_falls_back_rather_than_raising(tmp_path):
    from collections import Counter
    f = tmp_path / "vhn_denv1.tab27"
    f.write_text("")
    (tmp_path / "vhn_denv1.tab27.source.json").write_text("{not json")
    got = rows_per_taxon({f: Counter({"rows": 2})}, _VIRUSES, {})
    assert got["11053"]["rows"] == 2


# --------------------------------------- the per-file tally is really filled
def test_the_per_file_tally_matches_the_totals(world):
    """If this stayed empty every zero would be reported as "nothing
    fetched" -- the same confident wrong answer in a new shape, and one
    that reads as correct because the old answer said the same thing.
    """
    index = AliasIndex.from_releases(world / "spine", world / "host")
    per_file: dict = {}
    stats = ingest_mitab(Emit(), [world / "vhn_test.tab27"], index,
                         {"2697049": "2697049"}, "infores:virhostnet",
                         "2024-01-01", per_file=per_file)
    assert list(per_file) == [world / "vhn_test.tab27"]
    c = per_file[world / "vhn_test.tab27"]
    assert c["rows"] == stats["rows"] > 0
    assert c["edges"] == stats["edges"] > 0
    for gate in ("duplicate", "host_host", "viral_viral", "unresolved_id",
                 "unmapped_viral_taxon"):
        assert c[gate] == stats[gate], gate


def test_rows_split_across_two_files_are_counted_against_each(world):
    """One Counter for the lot would attribute every drop to whichever
    virus happened to be read first."""
    second = world / "vhn_other.tab27"
    second.write_text((world / "vhn_test.tab27").read_text())
    index = AliasIndex.from_releases(world / "spine", world / "host")
    per_file: dict = {}
    stats = ingest_mitab(Emit(), [world / "vhn_test.tab27", second], index,
                         {"2697049": "2697049"}, "infores:virhostnet",
                         "2024-01-01", per_file=per_file)
    assert len(per_file) == 2
    assert sum(c["rows"] for c in per_file.values()) == stats["rows"]
    # The second file is a copy, so every one of its rows is a duplicate of
    # one already seen -- and the duplicate must land on the file that holds
    # the later copy, not be smeared over both.
    assert per_file[second]["edges"] == 0
    assert per_file[world / "vhn_test.tab27"]["edges"] > 0


def test_omitting_per_file_still_works(world):
    """Every other caller passes nothing and must be unaffected."""
    index = AliasIndex.from_releases(world / "spine", world / "host")
    stats = ingest_mitab(Emit(), [world / "vhn_test.tab27"], index,
                         {"2697049": "2697049"}, "infores:virhostnet",
                         "2024-01-01")
    assert stats["edges"] > 0


# ------------------- attribution to a reference polyprotein, not identity
def test_a_mapped_accession_resolves_where_it_did_not_before(world):
    """VirHostNet curates some flavivirus rows against an unreviewed,
    strain-specific polyprotein. Without the map the row dies unresolved;
    DENV-1 and TBEV lost 89 real host-viral rows that way."""
    index = AliasIndex.from_releases(world / "spine", world / "host")
    line = _row("uniprotkb:STRAIN9", f"uniprotkb:{ACE2.split(':')[1]}",
                "MI:0018", "99", 2697049, 9606, "MI:0407", "0.5")
    assert classify_row(line, index, {"2697049": "2697049"}).gate == \
        "unresolved_id"
    v = classify_row(line, index, {"2697049": "2697049"},
                     {"STRAIN9": "P0DTC7"})
    assert v.gate is None
    assert v.viral == ORF7A


def test_the_measured_accession_is_remembered_not_discarded(world):
    """An attribution that left no trace is indistinguishable from a
    measurement on the reference entry -- the silent identity claim this
    design exists to avoid."""
    index = AliasIndex.from_releases(world / "spine", world / "host")
    line = _row("uniprotkb:STRAIN9", f"uniprotkb:{ACE2.split(':')[1]}",
                "MI:0018", "99", 2697049, 9606, "MI:0407", "0.5")
    v = classify_row(line, index, {"2697049": "2697049"},
                     {"STRAIN9": "P0DTC7"})
    assert v.measured == "STRAIN9"


def test_an_unmapped_row_carries_no_measured_accession(world):
    """Only a rewritten row gets the qualifier, so its presence means
    something."""
    index = AliasIndex.from_releases(world / "spine", world / "host")
    line = _row(f"uniprotkb:P0DTD1-{NSP5_PP1A.split(':')[1]}",
                f"uniprotkb:{ACE2.split(':')[1]}",
                "MI:0018", "99", 2697049, 9606, "MI:0407", "0.5")
    assert classify_row(line, index, {"2697049": "2697049"},
                        {"STRAIN9": "P0DTC7"}).measured is None


def test_the_qualifier_reaches_the_edge(world):
    index = AliasIndex.from_releases(world / "spine", world / "host")
    f = world / "vhn_mapped.tab27"
    f.write_text(_row("uniprotkb:STRAIN9", f"uniprotkb:{ACE2.split(':')[1]}",
                      "MI:0018", "99", 2697049, 9606, "MI:0407", "0.5"))
    em = Emit()
    stats = ingest_mitab(em, [f], index, {"2697049": "2697049"},
                         "infores:virhostnet", "2024-01-01",
                         accession_map={"STRAIN9": "P0DTC7"})
    assert stats["edges"] == 1 and stats["attributed"] == 1
    assert em.edges[0]["qualifiers"]["measured_accession"] == "STRAIN9"


def test_no_same_as_edge_is_emitted_for_an_attribution(world):
    """THE POINT OF THE WHOLE DESIGN. A SAME_AS would make AliasIndex treat
    two different UniProt entries, from different isolates with different
    sequences, as one protein. This project refuses that claim elsewhere
    and must not make it here by a side door."""
    index = AliasIndex.from_releases(world / "spine", world / "host")
    f = world / "vhn_mapped2.tab27"
    f.write_text(_row("uniprotkb:STRAIN9", f"uniprotkb:{ACE2.split(':')[1]}",
                      "MI:0018", "99", 2697049, 9606, "MI:0407", "0.5"))
    em = Emit()
    ingest_mitab(em, [f], index, {"2697049": "2697049"},
                 "infores:virhostnet", "2024-01-01",
                 accession_map={"STRAIN9": "P0DTC7"})
    assert [e for e in em.edges if e["predicate"] == "SAME_AS"] == []
    assert all(e["predicate"] == "PHYSICALLY_INTERACTS_WITH"
               for e in em.edges)


def test_the_attribution_is_not_a_fallback_for_a_resolvable_id(world):
    """Applying it only when resolution fails would make the graph's shape
    depend on which proteomes a build happened to load: the same row would
    attribute to the reference in one build and to a strain entry in
    another, with nothing saying which."""
    index = AliasIndex.from_releases(world / "spine", world / "host")
    line = _row(f"uniprotkb:{ORF7A.split(':')[1]}",
                f"uniprotkb:{ACE2.split(':')[1]}",
                "MI:0018", "99", 2697049, 9606, "MI:0407", "0.5")
    # ORF7A resolves perfectly well on its own; the map still fires and
    # sends it somewhere visibly different.
    v = classify_row(line, index, {"2697049": "2697049"},
                     {"P0DTC7": NSP5.split(":")[1]})
    assert v.viral == NSP5
    assert v.measured == "P0DTC7"


def test_the_map_applies_to_either_column(world):
    """A MITAB row puts the viral protein in whichever column the curator
    chose, and the real data uses both: the DENV-1 rows include partners
    recorded as 9606 x 11053, with the virus in column B. A map that only
    rewrote column A would lose exactly those rows while appearing to
    work on the others.
    """
    index = AliasIndex.from_releases(world / "spine", world / "host")
    amap = {"STRAIN9": "P0DTC7"}
    host = f"uniprotkb:{ACE2.split(':')[1]}"
    in_a = _row("uniprotkb:STRAIN9", host, "MI:0018", "99",
                2697049, 9606, "MI:0407", "0.5")
    in_b = _row(host, "uniprotkb:STRAIN9", "MI:0018", "98",
                9606, 2697049, "MI:0407", "0.5")
    for line in (in_a, in_b):
        v = classify_row(line, index, {"2697049": "2697049"}, amap)
        assert v.gate is None, line
        assert v.viral == ORF7A and v.host == ACE2
        assert v.measured == "STRAIN9"


def test_a_host_side_rewrite_is_not_reported_as_the_viral_measurement(world):
    """`measured` describes the VIRAL endpoint. If a map ever rewrote the
    host side too, attributing that to the viral protein would put a human
    accession in a qualifier that says which virus entry was assayed."""
    index = AliasIndex.from_releases(world / "spine", world / "host")
    line = _row(f"uniprotkb:{ORF7A.split(':')[1]}", "uniprotkb:HOSTX",
                "MI:0018", "97", 2697049, 9606, "MI:0407", "0.5")
    v = classify_row(line, index, {"2697049": "2697049"},
                     {"HOSTX": ACE2.split(":")[1]})
    assert v.gate is None and v.host == ACE2 and v.viral == ORF7A
    assert v.measured is None


# ------------------------------- which service a raw file actually came from
def test_the_sidecar_names_the_service(tmp_path):
    from kgav.vhppi import file_source
    f = tmp_path / "anything.tab27"
    f.write_text("")
    (tmp_path / "anything.tab27.source.json").write_text('{"source": "infores:mint"}')
    assert file_source(f) == "infores:mint"


def test_the_sidecar_beats_a_misleading_filename(tmp_path):
    """THE WHOLE POINT. fetch_vhppi takes --base and will query any
    PSICQUIC service. A MINT file landing under VirHostNet's naming
    scheme must not be recorded as VirHostNet: the tier ablation reads
    these labels as fact."""
    from kgav.vhppi import file_source
    f = tmp_path / "vhn_denv1.tab27"
    f.write_text("")
    (tmp_path / "vhn_denv1.tab27.source.json").write_text('{"source": "infores:mint"}')
    assert file_source(f) == "infores:mint"


@pytest.mark.parametrize("name,src", [
    ("vhn_denv1.tab27", "infores:virhostnet"),
    ("intact_viral.tab27", "infores:intact"),
])
def test_a_file_predating_sidecars_falls_back_to_its_prefix(tmp_path, name, src):
    """The coronavirus files were placed by hand and their names are the
    only record of where they came from."""
    from kgav.vhppi import file_source
    f = tmp_path / name
    f.write_text("")
    assert file_source(f) == src


def test_an_unknown_file_has_no_source_rather_than_a_guessed_one(tmp_path):
    """None is the signal the ingest refuses on. Inventing a label here
    is the false provenance claim this function exists to stop."""
    from kgav.vhppi import file_source
    f = tmp_path / "mint_denv4.tab27"
    f.write_text("")
    assert file_source(f) is None


def test_a_corrupt_sidecar_falls_back_to_the_prefix(tmp_path):
    from kgav.vhppi import file_source
    f = tmp_path / "vhn_x.tab27"
    f.write_text("")
    (tmp_path / "vhn_x.tab27.source.json").write_text("{not json")
    assert file_source(f) == "infores:virhostnet"


def test_a_sidecar_with_no_source_field_falls_back(tmp_path):
    """Sidecars written before --source existed have every other field."""
    from kgav.vhppi import file_source
    f = tmp_path / "vhn_y.tab27"
    f.write_text("")
    (tmp_path / "vhn_y.tab27.source.json").write_text('{"url": "http://x"}')
    assert file_source(f) == "infores:virhostnet"


def test_the_ingest_refuses_a_file_whose_source_nothing_records(tmp_path,
                                                                world):
    """Rather than labelling it with whichever group it was globbed into.
    Guessing here is the false provenance claim the grouping exists to
    prevent, and the tier ablation reads these labels as fact.
    """
    import subprocess
    import sys
    raw = tmp_path / "raw"
    raw.mkdir()
    (raw / "mint_denv4.tab27").write_text("")
    out = subprocess.run(
        [sys.executable, str(ROOT / "scripts/ingest_vhppi.py"),
         "--vhppi-dir", str(raw), "--spine", str(world / "spine"),
         "--host", str(world / "host"), "--out", str(tmp_path / "o")],
        capture_output=True, text=True, timeout=120, check=False)
    assert out.returncode != 0
    assert "no source recorded for: mint_denv4.tab27" in (out.stdout + out.stderr)
    assert "Refusing rather than guessing" in (out.stdout + out.stderr)


# ------------------------- who curated a row, versus who served the file
def _row_src(line: str, src: str) -> str:
    """Set MITAB column 12 on a row built by _row."""
    f = line.split("\t")
    f[12] = src
    return "\t".join(f)


def test_the_curator_comes_from_the_row_not_the_file(world):
    """THE LIVE DEFECT. IntAct's release carries IMEx partner records:
    of 12,496 viral rows, 1,059 are MINT, 110 DIP, and two others. All
    were filed as infores:intact, because the label came from the glob
    the file matched. The tier ablation reads these labels.
    """
    index = AliasIndex.from_releases(world / "spine", world / "host")
    f = world / "intact_mixed.tab27"
    f.write_text("\n".join(
        _row_src(_row(f"uniprotkb:P0DTD1-{NSP5_PP1A.split(':')[1]}",
                      f"uniprotkb:{ACE2.split(':')[1]}", "MI:0018", str(p),
                      2697049, 9606, "MI:0407", "0.5"), src)
        for p, src in ((1, 'psi-mi:"MI:0469"(IntAct)'),
                       (2, 'psi-mi:"MI:0471"(MINT)'),
                       (3, 'psi-mi:"MI:0465"(DIP)'))))
    em = Emit()
    ingest_mitab(em, [f], index, {"2697049": "2697049"},
                 "infores:intact", "2024-01-01")
    got = sorted(e["primary_knowledge_source"] for e in em.edges)
    assert got == ["infores:dip", "infores:intact", "infores:mint"]


def test_a_row_with_no_source_column_falls_back_and_is_counted(world):
    """A silent fallback is how 1,171 rows came to be filed wrongly.
    The fallback is correct -- the serving service is the best
    available answer -- but it has to be visible."""
    index = AliasIndex.from_releases(world / "spine", world / "host")
    stats = ingest_mitab(Emit(), [world / "vhn_test.tab27"], index,
                         {"2697049": "2697049"}, "infores:virhostnet",
                         "2024-01-01")
    assert stats["source_db_absent"] == stats["edges"] > 0
    assert stats["src:infores:virhostnet"] == stats["edges"]


def test_a_self_curated_file_keeps_its_label(world):
    """VirHostNet's own rows say virhostnet in column 12, so the fix
    must not move the 12,638 edges that were already right."""
    index = AliasIndex.from_releases(world / "spine", world / "host")
    f = world / "vhn_self.tab27"
    f.write_text(_row_src(
        _row(f"uniprotkb:P0DTD1-{NSP5_PP1A.split(':')[1]}",
             f"uniprotkb:{ACE2.split(':')[1]}", "MI:0018", "9",
             2697049, 9606, "MI:0407", "0.5"),
        'psi-mi:"MI:1114"(virhostnet)'))
    em = Emit()
    ingest_mitab(em, [f], index, {"2697049": "2697049"},
                 "infores:virhostnet", "2024-01-01")
    assert em.edges[0]["primary_knowledge_source"] == "infores:virhostnet"


@pytest.mark.parametrize("name,infores", [
    ("virhostnet", "infores:virhostnet"), ("IntAct", "infores:intact"),
    ("MINT", "infores:mint"), ("DIP", "infores:dip"),
    ("HPIDb", "infores:hpidb"), ("matrixdb", "infores:matrixdb"),
])
def test_every_curator_in_this_corpus_slugs_cleanly(name, infores):
    """All six names present across both files, derived rather than
    mapped: a MI-code table is a second place for the truth to drift."""
    from kgav.vhppi import source_infores
    assert source_infores(name) == infores


def test_an_unreadable_source_column_yields_no_curator():
    """None is the signal to fall back, and it must not be an empty
    infores id that joins to nothing."""
    from kgav.vhppi import source_infores
    assert source_infores(None) is None
    assert source_infores("") is None
    assert source_infores("!!!") is None


# ---------------------------------- an aggregator serves the same fact twice
def _dup_pair(world, src_a: str, src_b: str):
    """The same interaction in two files, with the given curators."""
    line = _row(f"uniprotkb:P0DTD1-{NSP5_PP1A.split(':')[1]}",
                f"uniprotkb:{ACE2.split(':')[1]}", "MI:0018", "77",
                2697049, 9606, "MI:0407", "0.5")
    a = world / "intact_one.tab27"
    b = world / "imex_one.tab27"
    a.write_text(_row_src(line, src_a))
    b.write_text(_row_src(line, src_b))
    return a, b


def test_an_aggregator_reserving_the_same_record_does_not_double_it(world):
    """IMEx serves IntAct's records: its SARS-CoV-2 count is within six
    of IntAct's own. With the label now read from the row, both copies
    say infores:intact -- one fact twice, not two sources. A per-call
    dedup set cannot see across services."""
    index = AliasIndex.from_releases(world / "spine", world / "host")
    a, b = _dup_pair(world, 'psi-mi:"MI:0469"(IntAct)',
                     'psi-mi:"MI:0469"(IntAct)')
    em, shared = Emit(), set()
    for f, svc in ((a, "infores:intact"), (b, "infores:imex")):
        ingest_mitab(em, [f], index, {"2697049": "2697049"}, svc,
                     "2024-01-01", seen=shared)
    assert len(em.edges) == 1
    assert em.edges[0]["primary_knowledge_source"] == "infores:intact"


def test_two_real_curators_of_one_fact_are_both_kept(world):
    """The opposite case, and it must not regress: the same interaction
    curated by VirHostNet AND by IntAct is two provenance records, and
    assemble_release merges them into one edge carrying both."""
    index = AliasIndex.from_releases(world / "spine", world / "host")
    a, b = _dup_pair(world, 'psi-mi:"MI:0469"(IntAct)',
                     'psi-mi:"MI:1114"(virhostnet)')
    em, shared = Emit(), set()
    for f, svc in ((a, "infores:intact"), (b, "infores:virhostnet")):
        ingest_mitab(em, [f], index, {"2697049": "2697049"}, svc,
                     "2024-01-01", seen=shared)
    assert len(em.edges) == 2
    assert {e["primary_knowledge_source"] for e in em.edges} == {
        "infores:intact", "infores:virhostnet"}


def test_without_a_shared_set_each_call_still_dedups_itself(world):
    """Callers that pass nothing keep the old behaviour."""
    index = AliasIndex.from_releases(world / "spine", world / "host")
    a, b = _dup_pair(world, 'psi-mi:"MI:0469"(IntAct)',
                     'psi-mi:"MI:0469"(IntAct)')
    em = Emit()
    for f in (a, b):
        ingest_mitab(em, [f], index, {"2697049": "2697049"},
                     "infores:intact", "2024-01-01")
    assert len(em.edges) == 2, "no shared set means no cross-call dedup"


def test_a_source_covered_by_an_aggregator_reports_its_real_content(world):
    """`edges` is order-dependent: groups run in sorted order and the
    first to carry a fact claims it. IMEx preceded IntAct, the same
    12,496 records read as "12,496 rows -> 10 edges", and that looks
    like a source contributing nothing rather than one fully covered.
    Renaming a file would have changed the figure.
    """
    index = AliasIndex.from_releases(world / "spine", world / "host")
    a, b = _dup_pair(world, 'psi-mi:"MI:0469"(IntAct)',
                     'psi-mi:"MI:0469"(IntAct)')
    em, shared = Emit(), set()
    first = ingest_mitab(em, [b], index, {"2697049": "2697049"},
                         "infores:imex", "2024-01-01", seen=shared)
    second = ingest_mitab(em, [a], index, {"2697049": "2697049"},
                          "infores:intact", "2024-01-01", seen=shared)
    assert first["edges"] == 1 and second["edges"] == 0
    # the order-free number is the same for both
    assert first["eligible"] == second["eligible"] == 1
    assert second["duplicate"] == 1


def test_eligible_counts_only_rows_that_passed_every_gate(world):
    """Not rows read. A row dropped as unparseable or viral-viral was
    never going to be an edge and must not inflate the figure."""
    index = AliasIndex.from_releases(world / "spine", world / "host")
    stats = ingest_mitab(Emit(), [world / "vhn_test.tab27"], index,
                         {"2697049": "2697049"}, "infores:virhostnet",
                         "2024-01-01")
    assert stats["eligible"] == stats["edges"] + stats["duplicate"]
    assert stats["eligible"] < stats["rows"], "the fixture has dropped rows"
