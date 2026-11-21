"""B3's measurement, and the figures section B states in prose.

WHY THIS FILE EXISTS. Section B of the fix list opened with "13,916 of
19,506 labelled compounds (71%) have no edge besides their label" and
nothing re-derived it. On 2026-11-10 the release said 15,979 of 21,371
(74.8%) -- both operands wrong, in the direction that UNDERSTATED the
problem, which is why no reader would have questioned it. That is the third
figure set this repository has found rotting in prose (L15 found the other
two), so this one is derived.

The pure functions are tested on a fixture release rather than only on
data/releases/v0.1, because a test that needs the 30 GB sqlite or a built
graph is a test that gets skipped and then gets believed.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import diagnose_label_reachability as D  # noqa: E402

DOC = ROOT / "doc2/records/audit_2026-11-03_fixes.md"
FIGURE = "<!-- FIGURE: label_reachability -->"
RELEASE = ROOT / "data/releases/v0.1"


# ------------------------------------------------------- the document figures
@pytest.mark.skipif(not (RELEASE / "edges.jsonl").exists(),
                    reason="no release on disk")
def test_section_Bs_header_is_the_releases():
    """BOTH OPERANDS AND THE PERCENTAGE, from the release.

    Proven to fail against the old sentence: this file's own history is the
    13,916 / 19,506 / 71% it replaced.
    """
    import re

    from kgav.baselines import Graph
    from kgav.schema import load_schema

    schema = load_schema()
    g = Graph.load(RELEASE, skip_predicates=schema.held_out_predicates(),
                   symmetry=schema.symmetry())
    traversable = {n for n, d in g.degree_unlabelled.items() if d > 0}
    fig = D.header_figures(RELEASE, traversable)

    text = DOC.read_text()
    assert FIGURE in text, "section B's header lost its marker"
    para = " ".join(text.split(FIGURE, 1)[1].split("\n\n")[0].split())
    m = re.search(r"\*\*([\d,]+) of ([\d,]+) labelled compounds "
                  r"\(([\d.]+)%\) have no edge besides their label\*\*", para)
    assert m, ("the header no longer states its three figures in the form this "
               "test reads; update both together")
    bare, labelled, share = m.group(1), m.group(2), m.group(3)
    assert int(bare.replace(",", "")) == fig["bare"], (
        f"header says {bare} bare, release says {fig['bare']:,}")
    assert int(labelled.replace(",", "")) == fig["labelled"], (
        f"header says {labelled} labelled, release says {fig['labelled']:,}")
    assert share == f"{fig['share_bare'] * 100:.1f}", (
        f"header says {share}%, release says {fig['share_bare']:.1%}")


@pytest.mark.skipif(not (RELEASE / "edges.jsonl").exists(),
                    reason="no release on disk")
def test_the_stale_figures_are_not_still_quoted_as_current():
    """The old numbers may appear as HISTORY and must not appear as fact.
    The header records them with the date they stopped being true; a bare
    "13,916 of 19,506" anywhere in the row would be the rot returning."""
    para = " ".join(DOC.read_text().split(FIGURE, 1)[1].split("\n\n")[0].split())
    assert "13,916" not in para and "19,506" not in para, (
        "the superseded figures are back in the header paragraph")


# -------------------------------------------------------- the pure functions
def _release(tmp_path, nodes, edges):
    d = tmp_path / "rel"
    d.mkdir()
    (d / "nodes.jsonl").write_text(
        "".join(json.dumps(n) + "\n" for n in nodes))
    (d / "edges.jsonl").write_text(
        "".join(json.dumps(e) + "\n" for e in edges))
    return d


def test_header_figures_counts_label_subjects_and_not_edges(tmp_path):
    """A compound with THREE label edges is one labelled compound. The
    figure is a compound count and the file is an edge list, which is the
    confusion that made L9 read 217 where the release said 23."""
    rel = _release(tmp_path, [], [
        {"subject": "A", "predicate": "HAS_ANTIVIRAL_ACTIVITY_AGAINST",
         "object": "V1"},
        {"subject": "A", "predicate": "HAS_ANTIVIRAL_ACTIVITY_AGAINST",
         "object": "V2"},
        {"subject": "A", "predicate": "MEASURED_INACTIVE_AGAINST", "object": "V3"},
        {"subject": "B", "predicate": "MEASURED_INACTIVE_AGAINST", "object": "V1"},
    ])
    fig = D.header_figures(rel, traversable=set())
    assert fig["labelled"] == 2
    assert fig["bare"] == 2
    assert fig["share_bare"] == 1.0


def test_header_figures_calls_a_compound_reachable_only_via_TRAVERSABLE(tmp_path):
    """`bare` is "no edge besides its label", so the caller passes the nodes
    with a traversable degree -- computed with the label predicates held
    out. Counting a compound as reachable because of its OWN label edge
    would make the figure 0 by construction and the section meaningless."""
    rel = _release(tmp_path, [], [
        {"subject": "A", "predicate": "HAS_ANTIVIRAL_ACTIVITY_AGAINST",
         "object": "V1"},
        {"subject": "B", "predicate": "HAS_ANTIVIRAL_ACTIVITY_AGAINST",
         "object": "V1"},
    ])
    assert D.header_figures(rel, traversable={"A"})["bare"] == 1
    assert D.header_figures(rel, traversable={"A", "B"})["bare"] == 0


def test_the_label_predicates_are_the_schemas_held_out_ones():
    """THE TWO LISTS MUST AGREE OR THE FIGURE COUNTS THE WRONG POPULATION.
    This script names the label predicates to find labelled compounds; the
    graph holds predicates out of traversal by the schema's declaration. If
    the schema ever holds out a predicate this does not count, a compound
    labelled only by that predicate would be invisible here while still
    unreachable in the graph."""
    from kgav.schema import load_schema
    held = load_schema().held_out_predicates()
    extra = set(D.LABEL_PREDICATES) - held
    assert not extra, (
        f"{sorted(extra)} are counted as label predicates but are TRAVERSABLE; "
        f"a path across them is retrieval, not prediction")
    unused = held - set(D.LABEL_PREDICATES)
    assert unused == {"IN_TRIAL_FOR", "TREATS"}, (
        f"the schema holds out {sorted(unused)}, which this script does not "
        f"count as labels. That was true of IN_TRIAL_FOR and TREATS because "
        f"v0.1 has 0 edges of either; a new one needs a decision, not a "
        f"silent omission")


@pytest.mark.skipif(not (RELEASE / "edges.jsonl").exists(),
                    reason="no release on disk")
def test_the_clinical_labels_are_excluded_by_decision_and_the_cost_is_stated():
    """THIS TEST USED TO ASSERT BOTH PREDICATES WERE EMPTY, AND IT FIRED.

    Its note said: "if either predicate gains edges, the header's
    population is wrong and the assertion above is the one that should
    fail first". The indications layer gave them 577 edges on
    2026-11-17, so the question it was holding open had to be answered
    rather than re-asserted.

    THE ANSWER: this diagnostic's population stays the in-vitro label,
    because its figures describe the compounds the baselines score and
    the clinical labels are a separate population with their own
    registration, cell and control. The decision is recorded in the
    script beside the constant, together with what it hides -- the
    compounds carrying ONLY a clinical label, which this script cannot
    see. This test holds that number to the release so the cost cannot
    drift away from the comment that states it.
    """
    invitro: set[str] = set()
    clinical: set[str] = set()
    for line in (RELEASE / "edges.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        if e["predicate"] in D.LABEL_PREDICATES:
            invitro.add(e["subject"])
        elif e["predicate"] in D.CLINICAL_LABEL_PREDICATES:
            clinical.add(e["subject"])
    hidden = len(clinical - invitro)
    src = Path(D.__file__).read_text()
    assert f"**{hidden} carry only a clinical label**" in src, (
        f"{hidden} compounds carry only a clinical label and are invisible "
        "to this diagnostic; the constant's comment states a different "
        "number, so the recorded cost of the decision is stale")
    assert f"{len(clinical)} compounds" in src, (
        f"{len(clinical)} compounds carry a clinical label; the comment "
        "states a different total")


def test_mechanism_nodes_joins_on_the_chembl_id_the_ingest_used(tmp_path, monkeypatch):
    """The graph is keyed on InChIKey and drug_mechanism on chembl_id, so
    the join goes through the node's own `chembl_id` property. A node
    without one cannot be matched and must not be guessed at."""
    monkeypatch.setattr(D, "__name__", D.__name__)  # no-op, keeps import local
    import kgav.host_targets as H
    monkeypatch.setattr(H, "load_mechanisms",
                        lambda db: {("CHEMBL25", "P00001"): "INHIBITOR"})
    got = D.mechanism_nodes(Path("/nowhere"), {
        "CHEMBL25": "INCHIKEY:AAA-BBB-N",
        "CHEMBL999": "INCHIKEY:CCC-DDD-N",
    })
    assert got == {"INCHIKEY:AAA-BBB-N"}


def test_the_definitions_are_ordered_loosest_last():
    """The report is read as a narrowing sequence and the ratio column is
    compared across it, so the order is part of the output, not cosmetic."""
    assert D.DEFINITIONS == ("mechanism", "targets", "mech_or_tgt", "any_edge")


def test_the_driver_refuses_a_missing_db_rather_than_skipping_it(tmp_path):
    """A --db that does not exist must fail, not fall through to the
    no-mechanism path: the caller asked for those columns and silence would
    read as "no compound has a mechanism"."""
    ap = D.build_parser(ROOT)
    args = ap.parse_args(["--db", str(tmp_path / "absent.db")])
    assert args.db is not None and not args.db.exists()
    assert D.main(["--db", str(tmp_path / "absent.db")]) == 1
