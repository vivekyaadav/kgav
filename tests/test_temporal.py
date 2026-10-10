"""Gate for the temporal split: the future is removed from the training graph,
and the audits fail loudly when it is not.
"""
import json

import pytest

from kgav.schema import load_schema
from kgav.temporal import audit, build_split, edge_year, write_split

V = "NCBITaxon:2697049"
OLD = "INCHIKEY:OLD"
NEW = "INCHIKEY:NEW"


def _e(s, p, o, date, q=None, src="s", tier=1):
    return {"subject": s, "predicate": p, "object": o, "qualifiers": q or {},
            "primary_knowledge_source": src, "evidence_tier": tier,
            "first_asserted_date": date}


def _label(drug, date, relation="="):
    return _e(drug, "HAS_ANTIVIRAL_ACTIVITY_AGAINST", V, date,
              {"assay_type": "cell_based_antiviral", "ec50_nm": 100.0,
               "relation": relation})


@pytest.fixture
def held():
    return load_schema().held_out_predicates()


@pytest.fixture
def release(tmp_path):
    def prot(pid, viral, fam=None):
        p = {"taxon_id": V if viral else "NCBITaxon:9606", "is_viral": viral,
             "sequence_hash": "h", "reviewed": True}
        if fam:
            p["protein_family"] = fam
        return {"id": pid, "class": "Protein", "properties": p}

    def drug(n):
        return {"id": n, "class": "SmallMolecule",
                "properties": {"smiles": "C", "inchikey_skel": n[-2:],
                               "is_approved": True, "salt_collapsed": False,
                               "stereo_collapsed": False}}

    nodes = [{"id": V, "class": "OrganismTaxon",
              "properties": {"label": "SARS-CoV-2", "family": "Coronaviridae",
                             "baltimore_class": "IV", "is_enveloped": True}},
             prot("UniProtKB:NSP5", True, "nsp5"),
             {"id": "KGAV:G", "class": "Gene",
              "properties": {"symbol": "rep", "taxon_id": V, "is_viral": True}},
             drug(OLD), drug(NEW)]
    edges = [
        _e(OLD, "INHIBITS", "UniProtKB:NSP5", "2020-01-01",
           {"assay_type": "biochemical", "ic50_nm": 25.0}),
        _e(NEW, "INHIBITS", "UniProtKB:NSP5", "2023-01-01",
           {"assay_type": "biochemical", "ic50_nm": 10.0}),
        _e("UniProtKB:NSP5", "ENCODED_BY", "KGAV:G", "2020-01-01"),
        _e("KGAV:G", "BELONGS_TO", V, "2020-01-01"),
        _e(OLD, "CHEMICALLY_SIMILAR_TO", NEW, "1970-01-01",
           {"tanimoto_ecfp4": 0.9}, "infores:kgav-computed", 3),
        _e(NEW, "TARGETS", "UniProtKB:NSP5", "1970-01-01",
           {"direction": "unknown", "assay_type": "binding", "ic50_nm": 5.0}),
        _label(OLD, "2020-06-01"),
        _label(NEW, "2023-06-01"),
        _label(NEW, "2024-01-01", relation=">"),
    ]
    (tmp_path / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in edges))
    return tmp_path


# ------------------------------------------------------------------- parsing
def test_placeholder_date_is_not_a_year():
    assert edge_year({"first_asserted_date": "2021-03-04"}) == 2021
    assert edge_year({"first_asserted_date": "1970-01-01"}) is None
    assert edge_year({}) is None


# --------------------------------------------------------------------- split
def test_future_evidence_is_removed_from_the_training_graph(release, held):
    """THE trap this module exists for. Filtering only the labels leaves a
    compound whose nsp5 activity was published in 2023 with an M1 path built
    from 2023 evidence predicting a 2023 label."""
    s = build_split(release, 2021, held)
    assert s.stats["edge_dropped_future"] == 1
    for e in s.train_edges:
        y = edge_year(e)
        assert y is None or y <= 2021


def test_labels_are_partitioned_by_year(release, held):
    s = build_split(release, 2021, held)
    assert s.train_labels[V] == {OLD}
    assert s.test_labels[V] == {NEW}


def test_censored_labels_are_never_positives(release, held):
    """'EC50 > 9999 nM' is a measurement of INACTIVITY; using it as a positive
    would invert the target.

    It is no longer DISCARDED either. Dropping every censored row threw away
    4,685 measurements -- the best negatives the project has -- and left the
    temporal split with no negatives of its own. Whether such a row becomes a
    negative or stays undecidable is labels.classify's judgement: this one is
    bounded below 10 uM, so it decides nothing.
    """
    s = build_split(release, 2021, held)
    censored = "INCHIKEY:CENSORED"
    assert all(censored not in d for d in s.test_labels.values())
    assert all(censored not in d for d in s.train_labels.values())
    assert s.stats["label_undecidable"] == 1


def test_held_out_predicates_never_enter_the_training_graph(release, held):
    s = build_split(release, 2021, held)
    assert not [e for e in s.train_edges if e["predicate"] in held]


# ------------------------------------------------------------ undated policy
def test_undated_include_keeps_computed_and_real(release, held):
    s = build_split(release, 2021, held, "include")
    assert s.stats["undated_kept_computed"] == 1   # similarity
    assert s.stats["undated_kept"] == 1            # the dateless TARGETS edge


def test_undated_exclude_drops_everything_dateless(release, held):
    """Conservative: an unknown date may be post-cutoff."""
    s = build_split(release, 2021, held, "exclude")
    assert s.stats["undated_dropped"] == 2
    assert all(edge_year(e) is not None for e in s.train_edges)


def test_undated_computed_only_keeps_similarity_but_not_unknown_dates(release, held):
    s = build_split(release, 2021, held, "computed_only")
    assert s.stats["undated_kept_computed"] == 1
    assert s.stats["undated_dropped"] == 1


def test_unknown_policy_is_rejected(release, held):
    with pytest.raises(ValueError):
        build_split(release, 2021, held, "nonsense")


# -------------------------------------------------------------------- audits
def test_clean_split_passes_the_audits(release, held):
    s = build_split(release, 2021, held)
    problems = audit(s, release, held)
    assert not [p for p in problems if p.startswith(("L3", "L5", "label leak"))]


def test_l3_detects_planted_future_bleed(release, held):
    s = build_split(release, 2021, held)
    s.train_edges.append(_e("X", "INHIBITS", "Y", "2024-01-01"))
    assert any(p.startswith("L3") for p in audit(s, release, held))


def test_label_leak_is_detected(release, held):
    s = build_split(release, 2021, held)
    s.train_edges.append(_label(NEW, "2020-01-01"))
    assert any(p.startswith("label leak") for p in audit(s, release, held))


def test_l5_detects_a_compound_in_both_label_sets(release, held):
    s = build_split(release, 2021, held)
    s.train_labels[V].add(NEW)
    assert any(p.startswith("L5") for p in audit(s, release, held))


def test_recall_ceiling_is_reported_not_hidden(release, held):
    """A test compound with no pre-cutoff evidence is unpredictable by
    construction. Not a leak, but it caps achievable recall."""
    s = build_split(release, 2021, held)
    s.test_labels[V].add("INCHIKEY:GHOST")
    assert any(p.startswith("ceiling") for p in audit(s, release, held))


# --------------------------------------------------------------------- output
def test_written_graph_keeps_only_reachable_nodes(release, held, tmp_path):
    s = build_split(release, 2021, held)
    out = tmp_path / "train"
    write_split(s, release, out)
    ids = {json.loads(x)["id"] for x in (out / "nodes.jsonl").read_text().splitlines()
           if x.strip()}
    edges = [json.loads(x) for x in (out / "edges.jsonl").read_text().splitlines()
             if x.strip()]
    for e in edges:
        assert e["subject"] in ids and e["object"] in ids
    labels = json.loads((out / "TEST_LABELS.json").read_text())
    assert labels["cutoff"] == 2021
    assert labels["test"][V] == [NEW]


# -------------------------------------------- label polarity via classify (M2)
def _activity(subject, obj, value, relation="=", date="2023-01-01"):
    return {"subject": subject, "predicate": "HAS_ANTIVIRAL_ACTIVITY_AGAINST",
            "object": obj, "qualifiers": {"ec50_nm": value, "relation": relation,
                                          "assay_type": "cell_based_antiviral"},
            "primary_knowledge_source": "s", "evidence_tier": 1,
            "first_asserted_date": date}


def _release(tmp_path, edges):
    (tmp_path / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in edges))
    (tmp_path / "nodes.jsonl").write_text("")
    return tmp_path


def test_an_exact_measurement_above_the_threshold_is_a_negative(tmp_path):
    """578 labels in v0.1 read "= 50000 nM" and were counted as POSITIVES
    because the relation was exact. They are measurements of inactivity."""
    r = _release(tmp_path, [_activity("INCHIKEY:A", V, 50_000.0)])
    s = build_split(r, 2021, {"HAS_ANTIVIRAL_ACTIVITY_AGAINST"})
    assert s.test_labels.get(V, set()) == set()
    assert s.test_negatives[V] == {"INCHIKEY:A"}


def test_a_censored_measurement_is_kept_as_a_negative(tmp_path):
    """4,685 censored rows were discarded as "label_censored". A compound
    assayed and found inactive is the best negative this project has."""
    r = _release(tmp_path, [_activity("INCHIKEY:B", V, 20_000.0, relation=">")])
    s = build_split(r, 2021, {"HAS_ANTIVIRAL_ACTIVITY_AGAINST"})
    assert s.test_negatives[V] == {"INCHIKEY:B"}
    assert not s.stats["label_undecidable"]


def test_a_weak_lower_bound_stays_undecidable(tmp_path):
    """"> 100 nM" does not prove inactivity: the bound sits below the point
    where a hit would stop being pursued."""
    r = _release(tmp_path, [_activity("INCHIKEY:C", V, 100.0, relation=">")])
    s = build_split(r, 2021, {"HAS_ANTIVIRAL_ACTIVITY_AGAINST"})
    assert s.stats["label_undecidable"] == 1
    assert not s.test_labels.get(V) and not s.test_negatives.get(V)


def test_a_compound_measured_both_ways_is_dropped_from_both_sides(tmp_path):
    """labels.py's rule, applied within each side of the split."""
    r = _release(tmp_path, [_activity("INCHIKEY:D", V, 50.0),
                            _activity("INCHIKEY:D", V, 50_000.0)])
    s = build_split(r, 2021, {"HAS_ANTIVIRAL_ACTIVITY_AGAINST"})
    assert "INCHIKEY:D" not in s.test_labels.get(V, set())
    assert "INCHIKEY:D" not in s.test_negatives.get(V, set())
    assert s.stats["test_ambiguous_dropped"] == 1


def test_negatives_are_written_to_test_labels_json(tmp_path):
    (tmp_path / "rel").mkdir()
    r = _release(tmp_path / "rel", [_activity("INCHIKEY:A", V, 50_000.0),
                                    _activity("INCHIKEY:E", V, 5.0)])
    s = build_split(r, 2021, {"HAS_ANTIVIRAL_ACTIVITY_AGAINST"})
    out = tmp_path / "split"
    write_split(s, r, out)
    payload = json.loads((out / "TEST_LABELS.json").read_text())
    assert payload["test"][V] == ["INCHIKEY:E"]
    assert payload["test_negatives"][V] == ["INCHIKEY:A"]
    assert "train_negatives" in payload


def test_split_label_sets_are_ordered_deterministically(tmp_path, held):
    """Same defect as labels.build_labels: the per-virus order came from set
    iteration, and it reaches the audit list and TEST_LABELS.json. Two runs
    over one release must produce the same file, or the provenance stamps
    cannot tell drift from interpreter noise."""
    viruses = ["NCBITaxon:2697049", "NCBITaxon:694009", "NCBITaxon:11137",
               "NCBITaxon:31631", "NCBITaxon:277944", "NCBITaxon:1335626"]
    nodes, edges = [], []
    for v in viruses:
        nodes.append({"id": v, "class": "OrganismTaxon",
                      "properties": {"label": v, "family": "Coronaviridae",
                                     "baltimore_class": "IV", "is_enveloped": True}})
        for drug, rel, date in ((OLD, "=", "2019-01-01"), (NEW, "=", "2023-01-01"),
                                ("INCHIKEY:NEG", ">", "2023-01-01")):
            e = _label(drug, date, rel)
            e["object"] = v
            if rel == ">":
                e["qualifiers"]["ec50_nm"] = 50_000.0
            edges.append(e)
    (tmp_path / "nodes.jsonl").write_text("\n".join(json.dumps(n) for n in nodes))
    (tmp_path / "edges.jsonl").write_text("\n".join(json.dumps(e) for e in edges))
    s = build_split(tmp_path, cutoff=2021, held_out=held)
    for side in (s.test_labels, s.train_labels, s.test_negatives, s.train_negatives):
        assert list(side) == sorted(side)


# --------------------------------------------------------------------------
# held_out is about TRAVERSAL; evaluation_label is about LABELS. Reading one
# as the other put viral protein nodes into the negative sets.
# --------------------------------------------------------------------------
NSP5 = "UniProtKB:NSP5"


@pytest.fixture
def labels_only():
    return load_schema().evaluation_label_predicates()


def _inactive_against_protein(drug, date):
    """A MEASURED_INACTIVE_AGAINST edge: held out, and NOT a label.

    Its object is a viral protein, so a split that treats every held-out
    predicate as a label keys its negative sets on a protein.
    """
    return _e(drug, "MEASURED_INACTIVE_AGAINST", NSP5, date,
              {"assay_type": "biochemical", "ic50_nm": 80000.0, "relation": "="})


@pytest.fixture
def release_with_protein_measurements(tmp_path, release):
    extra = [_inactive_against_protein(NEW, "2023-01-01"),
             _inactive_against_protein(OLD, "2019-01-01")]
    path = tmp_path / "rel2"
    path.mkdir()
    (path / "nodes.jsonl").write_text((release / "nodes.jsonl").read_text())
    (path / "edges.jsonl").write_text(
        (release / "edges.jsonl").read_text().rstrip("\n") + "\n"
        + "\n".join(json.dumps(e) for e in extra))
    return path


def test_label_sets_are_keyed_only_on_viruses(
        release_with_protein_measurements, held, labels_only):
    """Every key of every label set must be a virus, never a protein.

    Before this, MEASURED_INACTIVE_AGAINST was bucketed by e["object"], so
    'UniProtKB:NSP5' appeared alongside the taxa in train_negatives and
    test_negatives. On the real release that moved the reported negative count
    from 5,164 to 6,789 while positives did not move -- the signature of a
    predicate that only ever classifies inactive.
    """
    s = build_split(release_with_protein_measurements, 2021, held,
                    label_predicates=labels_only)
    for name, sets in (("test_negatives", s.test_negatives),
                       ("train_negatives", s.train_negatives),
                       ("test_labels", s.test_labels),
                       ("train_labels", s.train_labels)):
        assert all(k.startswith("NCBITaxon:") for k in sets), \
            f"{name} is keyed on something that is not a virus: {sorted(sets)}"


def test_held_out_non_label_is_stripped_but_never_counted_as_a_label(
        release_with_protein_measurements, held, labels_only):
    """It must leave the training graph AND leave the label counts alone."""
    s = build_split(release_with_protein_measurements, 2021, held,
                    label_predicates=labels_only)
    assert s.stats["held_out_not_a_label"] == 2
    assert not [e for e in s.train_edges
                if e["predicate"] == "MEASURED_INACTIVE_AGAINST"]

    baseline = build_split(release_with_protein_measurements, 2021, held,
                           label_predicates={"HAS_ANTIVIRAL_ACTIVITY_AGAINST"})
    for key in ("label_test_active", "label_test_inactive",
                "label_train_active", "label_train_inactive"):
        assert s.stats[key] == baseline.stats[key]


def test_a_label_that_is_not_held_out_is_refused(release, labels_only):
    """A label left traversable predicts itself. Caught before the training
    graph is built, not by an audit that runs after it."""
    with pytest.raises(ValueError, match="not in held_out"):
        build_split(release, 2021, held_out=set(),
                    label_predicates=labels_only)


def test_schema_rejects_a_label_on_the_wrong_object_class():
    """The flag cannot be put somewhere that reintroduces the bug."""
    from kgav.schema import Schema
    doc = dict(load_schema().raw)
    doc["edge_classes"] = [
        dict(ec, evaluation_label=True) if ec["predicate"] == "INHIBITS" else ec
        for ec in doc["edge_classes"]]
    codes = {v.code for v in Schema(doc).validate_labels()}
    assert "LABEL_BAD_OBJECT" in codes
    assert "LABEL_NOT_HELD_OUT" in codes


def test_the_shipped_schema_declares_consistent_labels():
    assert load_schema().validate_labels() == []


# ----------- a compound measured on both sides is not a prospective case
HELD = {"HAS_ANTIVIRAL_ACTIVITY_AGAINST"}


def _both_sides(tmp_path, value=100.0):
    """One compound measured pre- and post-cutoff, as the NCATS screen and a
    later ChEMBL paper did for 19 compounds in v0.1."""
    return _release(tmp_path, [
        _activity("INCHIKEY:BOTH", V, value, date="2020-01-01"),   # train
        _activity("INCHIKEY:BOTH", V, value, date="2024-01-01"),   # test
        _activity("INCHIKEY:ONLYTEST", V, value, date="2024-06-01"),
    ])


def test_a_compound_known_before_the_cutoff_leaves_the_test_set(tmp_path):
    """audit()'s L5 check already said a compound measured on both sides is
    already known, then refused -- the module stated the rule and declined to
    apply it. In v0.1 this was 19 compounds: the NCATS CPE screen is one
    publication dated 2020, so every label it contributes routes to train,
    and ChEMBL re-measured 19 of them in 2022-2024.
    """
    s = build_split(_both_sides(tmp_path), 2021, HELD)

    assert s.train_labels[V] == {"INCHIKEY:BOTH"}, "train keeps the evidence"
    assert s.test_labels[V] == {"INCHIKEY:ONLYTEST"}, \
        "the already-known compound must leave TEST, not train"
    assert s.stats["test_positive_already_known_dropped"] == 1
    assert any("already measured pre-2021" in v for v in s.violations), \
        "the drop must be recorded, not silent"

    # and the audit that used to abort the run now finds nothing
    assert not [p for p in audit(s, tmp_path, HELD)
                if p.startswith("L5")], "L5 should have nothing left to report"


def test_the_direction_of_the_drop_is_test_not_train(tmp_path):
    """Removing the pre-cutoff measurement instead would delete evidence that
    genuinely existed at the time, which is a different (and wrong) protocol:
    the question is what the pre-cutoff graph could have predicted.
    """
    s = build_split(_both_sides(tmp_path), 2021, HELD)
    assert "INCHIKEY:BOTH" in s.train_labels[V]
    assert "INCHIKEY:BOTH" not in s.test_labels[V]


def test_negatives_are_dropped_on_the_same_rule(tmp_path):
    """A compound already known INACTIVE is no more a prospective test case
    than one already known active. 14 of the 19 v0.1 overlaps were negatives:
    both sides carried HAS_ANTIVIRAL_ACTIVITY_AGAINST and both classified as
    inactive, because polarity comes from classify() and not the predicate.
    """
    s = build_split(_both_sides(tmp_path, value=50_000.0), 2021, HELD)
    assert s.train_negatives[V] == {"INCHIKEY:BOTH"}
    assert s.test_negatives[V] == {"INCHIKEY:ONLYTEST"}
    assert s.stats["test_negative_already_known_dropped"] == 1
    assert not s.stats["test_positive_already_known_dropped"]


def test_keep_reproduces_the_refusal(tmp_path):
    """The previous behaviour stays reachable, so the refusal this replaced
    can be reproduced rather than only described.
    """
    r = _both_sides(tmp_path)
    s = build_split(r, 2021, HELD, overlap_policy="keep")
    assert s.test_labels[V] == {"INCHIKEY:BOTH", "INCHIKEY:ONLYTEST"}
    assert not s.stats["test_positive_already_known_dropped"]
    assert [p for p in audit(s, r, HELD) if p.startswith("L5")], \
        "keep must still trip L5, or the refusal cannot be reproduced"


def test_an_unknown_overlap_policy_is_refused(tmp_path):
    import pytest
    with pytest.raises(ValueError, match="overlap_policy"):
        build_split(_both_sides(tmp_path), 2021, HELD, overlap_policy="nope")


def test_the_temporal_driver_defines_no_scoring_rules_of_its_own():
    """Three figures in this project have had two implementations, and in
    every case only one got fixed: AUC 0.826 in guards.py and tools.py; the
    coverage floor; and lift trustworthiness, where temporal_split.py kept a
    copy of lift_row from before `and pool > k` was added, so SARS-CoV-2 M7
    (pool 91), SARS-CoV M7 (pool 34) and MERS-CoV M7 (pool 14) printed
    unflagged while compare_metapaths.py flagged the same channel.

    The first version of this test ALLOWLISTED the duplicated
    MIN_TRUSTWORTHY_EXPECTATION on the grounds that it was not a coverage
    floor. That was true and beside the point: it was a duplicated constant
    backing a divergent copy of lift_row, which is what the test was for.
    """
    import ast
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "scripts/temporal_split.py"
    tree = ast.parse(src.read_text())

    for module, wanted in (("kgav.calibration", "MIN_REACHED_POSITIVES"),
                           ("kgav.baselines", "lift_row")):
        assert any(isinstance(n, ast.ImportFrom) and n.module == module
                   and any(a.name == wanted for a in n.names)
                   for n in ast.walk(tree)), \
            f"temporal_split must import {wanted} from {module}"

    # No local threshold or floor. DEFAULT_OVERLAP_POLICY is a policy name,
    # not a numeric rule, so the check is on numeric assignments only.
    own = [t.id for n in ast.walk(tree) if isinstance(n, ast.Assign)
           for t in n.targets if isinstance(t, ast.Name) and t.id.isupper()
           and isinstance(n.value, ast.Constant)
           and isinstance(n.value.value, (int, float))]
    assert not own, f"numeric scoring rule defined locally: {own}"


def test_a_pool_no_larger_than_k_is_flagged_in_the_temporal_protocol():
    """The behavioural half: lift = pool/k exactly when the pool fits inside
    k, so the figure is a function of pool size and nothing else. This is the
    case the temporal driver printed clean for every M7 row.
    """
    from kgav.baselines import lift_row

    scores = {f"D{i}": 1.0 / (i + 1) for i in range(30)}
    pos = {f"D{i}" for i in range(8)}

    small = lift_row(scores, pos, k=100)          # pool 30 <= k
    assert small["pool"] <= 100
    assert small["lift_is_pool_artifact"] is True
    assert small["trustworthy"] is False, \
        "a lift that is pool/k by construction must never read as trustworthy"
    assert small["expected"] >= 1.0, \
        "and expectation alone would have passed it, which was the bug"

    big = lift_row(scores, pos, k=5)              # pool 30 > k
    assert big["lift_is_pool_artifact"] is False
    assert big["trustworthy"] is True


def test_a_control_run_does_not_overwrite_the_result_it_controls():
    """--overlap-policy keep wrote to the same path as the real run, so
    reproducing the refusal destroyed the evaluation it was a control for.
    """
    import argparse
    import importlib.util
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "scripts/temporal_split.py"
    spec = importlib.util.spec_from_file_location("_ts", src)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    default = argparse.Namespace(cutoff=2021, undated="include",
                                 overlap_policy=mod.DEFAULT_OVERLAP_POLICY)
    control = argparse.Namespace(cutoff=2021, undated="include",
                                 overlap_policy="keep")

    assert mod._results_name(default) == "temporal_2021_include.json", \
        "the default must keep its existing filename"
    assert mod._results_name(control) != mod._results_name(default), \
        "the control must not land on the real run's path"
    assert "keep" in mod._results_name(control)
