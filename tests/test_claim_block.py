"""The claim block at the top of docs/evaluation.md, re-derived.

WHAT THIS IS FOR. Those three indented lines are the first thing a reader
sees, and nothing had ever checked them. On 2026-11-01 they were found to
carry three wrong figures:

    M6 reaches 1        -> reaches 7
    M8 emits no paths   -> reaches 82 drugs and 1 labelled active
    degree 0.584        -> 0.606 [0.544, 0.669]

The M8 line was the costly one. It is why the sentence above the block said
SEVEN reasoning channels: M8 was excluded from the count for emitting
nothing. It emits, so the count is eight and the "other five" are six. One
figure going stale had silently changed a number in the headline claim.

M4 and M5 -- the two the block is actually about -- were correct to the last
digit, which is what makes the others drift rather than a different protocol
being quoted.

WHY IT DRIFTED. Same reason as the sequence table: check-results compares a
results FILE against the graph, and the file was current. The document is
not the file. So this parses the block and re-derives every figure in it,
hard-coding none.

The block is left as indented text rather than converted to a table. The
parse is regular enough to check, and the document's opening should read the
way its author wrote it.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "evaluation.md"
RESULTS = ROOT / "data" / "results" / "hard_negatives.json"
sys.path.insert(0, str(ROOT / "src"))

from kgav.calibration import MIN_REACHED_POSITIVES
from kgav.labels import auc_interval

MARKER = "<!-- BLOCK: claim_channels"
VIRUS = "SARS-CoV-2"
CHANNELS = [f"M{i}" for i in range(1, 9)]
#: The two scorers that are NOT reasoning channels but DO have a line in
#: the block. COMBINED's was never checked by this file, so it sat at
#: 0.575 [0.512, 0.637] after the release moved it to 0.573 [0.510, 0.635]
#: -- and the census block ten lines further down carried the right one,
#: so the document disagreed with itself. `scripts/diagnose_doc_figures.py`
#: found it once that script stopped matching figures by value (D5); no
#: test could, because this list stopped at M8.
NON_CHANNELS = ["COMBINED", "degree"]

#: Only as far as the channel count can reach. A missing word must fail
#: loudly rather than fall back to a digit the document does not use.
WORDS = {0: "no", 1: "one", 2: "two", 3: "three", 4: "four", 5: "five",
         6: "six", 7: "seven", 8: "eight", 9: "nine", 10: "ten"}


# ---------------------------------------------------------------- derivation
def channels() -> dict[str, dict]:
    """Every channel's auc, interval and reached-positive count."""
    doc = json.loads(RESULTS.read_text())
    per = doc["cross_sectional"][VIRUS]
    out = {}
    # COMBINED TOO. The fixture built every M-channel plus degree, so
    # `chans["COMBINED"]` raised KeyError and the block's COMBINED line
    # was unreachable by any assertion here -- which is why it went
    # stale while the census block below it did not.
    for name in CHANNELS + NON_CHANNELS:
        m = per[name]
        lo, hi = auc_interval(m["auc"], m["n_pos"], m["n_neg"])
        out[name] = {"auc": m["auc"], "lo": lo, "hi": hi,
                     "n_pos": m["n_pos"],
                     # THE FREE COMPONENT, carried so a channel clearing
                     # chance can be attributed rather than announced: an
                     # imputed AUC rises on coverage asymmetry alone, and
                     # auc_floor is exactly that part. Added 2026-11-10,
                     # when three channels cleared 0.5 and nothing here
                     # could say whether any of it was ranking.
                     "floor": m.get("auc_floor"),
                     "reached": round(m["coverage_pos"] * m["n_pos"])}
    return out


def evaluated(chans: dict) -> list[str]:
    return [n for n in CHANNELS
            if chans[n]["reached"] >= MIN_REACHED_POSITIVES]


def fmt_auc(c: dict) -> str:
    return f"{c['auc']:.3f} [{c['lo']:.3f}, {c['hi']:.3f}]"


# -------------------------------------------------------------------- parsing
def block() -> str:
    text = DOC.read_text()
    assert MARKER in text, "the claim block lost its machine-readable marker"
    after = text.split(MARKER, 1)[1].split("-->", 1)[1]
    # Up to the first unindented prose line after the block.
    lines = []
    started = False
    for line in after.splitlines():
        if line.startswith("    "):
            started, _ = True, lines.append(line)
        elif started and line.strip():
            break
    assert lines, "no indented block follows the marker"
    return "\n".join(lines)


def parsed_evaluated(body: str) -> dict[str, str]:
    """channel -> the rendered 'auc [lo, hi]' the block prints for it."""
    return {m[1]: f"{m[2]} [{m[3]}, {m[4]}]" for m in re.finditer(
        r"(M\d)\s+(\d\.\d{3}) \[(\d\.\d{3}), (\d\.\d{3})\]\s+reaches", body)}


def parsed_reach(body: str) -> dict[str, int]:
    """channel -> the reached count the block prints, from either shape."""
    out = {n: int(k) for n, k in
           re.findall(r"(M\d)\s+reaches\s+(\d+)", body)}
    out.update({n: int(k) for n, k in
                re.findall(r"(M\d)\s+\d\.\d{3} \[[^\]]+\]\s+reaches\s+(\d+)",
                           body)})
    return out


def parsed_degree(body: str) -> str:
    m = re.search(r"degree\s+(\d\.\d{3}) \[(\d\.\d{3}), (\d\.\d{3})\]", body)
    assert m, f"no degree line in the block:\n{body}"
    return f"{m[1]} [{m[2]}, {m[3]}]"


# -------------------------------------------------------------- the checks
@pytest.fixture(scope="module")
def chans() -> dict:
    if not RESULTS.exists():
        pytest.skip("hard_negatives.json is not built in this checkout")
    return channels()


def test_the_block_names_every_channel(chans):
    """GUARDS EVERY CHECK BELOW. A channel missing from the block is a
    channel whose figure cannot be wrong, which is how M8 sat there
    claiming to emit nothing."""
    seen = parsed_reach(block())
    assert set(seen) == set(CHANNELS), (
        f"block names {sorted(seen)}; the results file has {CHANNELS}")


def test_the_evaluated_channels_are_the_ones_over_the_floor(chans):
    """Which channels get an AUC is a derived decision, not an editorial
    one: reached >= MIN_REACHED_POSITIVES. Printing an AUC for a channel
    below it is quoting the imputed floor."""
    assert sorted(parsed_evaluated(block())) == evaluated(chans)


@pytest.mark.parametrize("channel", CHANNELS)
def test_every_reach_count_matches_the_results_file(channel, chans):
    got = parsed_reach(block())[channel]
    want = chans[channel]["reached"]
    assert got == want, (
        f"the block says {channel} reaches {got}; the results file gives "
        f"{want}. Re-derive the block, do not edit the file to match it.")


def test_every_printed_interval_is_the_one_the_codebase_computes(chans):
    """Not a second Hanley-McNeil implementation living in a document."""
    printed = parsed_evaluated(block())
    assert printed, "no channel in the block carries an AUC"
    for name, cell in printed.items():
        assert cell == fmt_auc(chans[name]), f"{name}: {cell}"


def test_the_degree_line_matches_the_results_file(chans):
    assert parsed_degree(block()) == fmt_auc(chans["degree"])


def test_a_channel_clearing_chance_clears_it_only_through_COVERAGE(chans):
    """THE GUARD THIS REPLACES SAID "degree is the only scorer clear of
    chance", and its docstring said that if a reasoning channel ever
    cleared, "that failure is the result, not a test to relax". On
    2026-11-10 M3, M4 and M5 cleared, because fix-list B1 added 6,608
    TARGETS edges and no labels, so every host-directed channel reached
    further. The failure was the result and it is recorded in the claim
    above, not relaxed away.

    So this is not a weaker test, it is the sharp version of the same
    question. An imputed AUC rises for free when a scorer reaches more
    positives than negatives, and `auc_floor` is exactly that free
    component. A channel clearing 0.5 whose FLOOR also clears 0.5 has
    demonstrated reach. A channel clearing 0.5 whose floor does NOT --
    i.e. one that beats chance without a coverage asymmetry to explain it
    -- would be the first ranking result this project has ever produced,
    and it must not land quietly. That is what this fails on now.
    """
    clear = [n for n in CHANNELS
             if chans[n]["lo"] > 0.5 or chans[n]["hi"] < 0.5]
    text = DOC.read_text()
    for n in clear:
        floor = chans[n].get("floor")
        assert floor is not None, f"{n} records no auc_floor to attribute to"
        assert floor > 0.5, (
            f"{n} reads {fmt_auc(chans[n])} with a floor of {floor:.4f}, "
            f"which does NOT exclude chance -- so its separation is not "
            f"explained by coverage. This is a ranking result; write it up "
            f"rather than making this test pass")
        assert f"floor {floor:.4f}" in text, (
            f"{n} clears chance and the block does not print its floor; a "
            f"reader cannot tell reach from ranking without it")
    d = chans["degree"]
    assert d["lo"] > 0.5, (
        f"degree reads {fmt_auc(d)}, which contains chance -- the block "
        f"calls it clear of 0.5")


def test_the_counts_in_the_sentence_are_the_channel_totals(chans):
    """THE COUNT THAT WENT STALE WITHOUT BEING EDITED. The sentence said
    seven channels because M8 was excluded for emitting nothing; M8 emits,
    so it is eight and the remainder is six. A figure in the block had
    silently changed a number in the headline above it."""
    text = DOC.read_text()
    n_eval = len(evaluated(chans))
    sentence = (f"Of {WORDS[len(CHANNELS)]} reasoning channels, "
                f"{WORDS[n_eval]} were ever evaluated")
    assert sentence in text, f"expected: {sentence!r}"
    assert f"The other {WORDS[len(CHANNELS) - n_eval]} reach too few" in text


def test_the_block_is_marked_as_machine_checked():
    text = DOC.read_text()
    assert text.count(MARKER) == 1
    head = text.split(MARKER, 1)[1].split("-->", 1)[0]
    # Whitespace-normalised: the marker is a wrapped comment, and "Do not\n
    # hand-edit" is the phrase being present. A raw substring search calls
    # it missing -- which it did, on the first run of this test.
    flat = " ".join(head.split())
    assert "tests/test_claim_block.py" in flat
    assert "Do not hand-edit" in flat


# ------------------------------------------- the parser, on synthetic input
#
# The real-release checks above skip without a results file. These do not:
# the parsing is what has to be right for any of it to mean anything, and a
# parser that silently matched nothing would make every check vacuous.
SYNTHETIC = """
    evaluated   M4  0.507 [0.447, 0.568]   reaches 15 of 89 actives
                M5  0.516 [0.455, 0.577]   reaches 14 of 89

    not         M3  reaches  9    M6 reaches 7    M2 reaches 5
    evaluated   M1  reaches  3    M8 reaches 1    M7 reaches 0

    baseline    degree 0.606 [0.544, 0.669] — the only scorer clear of 0.5,
                and it is the null hypothesis, not a reasoning channel
"""


def test_the_parser_finds_every_channel_in_both_shapes():
    """Two shapes in one block: a channel with an AUC and a bare reach
    count. A parser that read only the second would never check an
    interval, and one that read only the first would miss six channels."""
    reach = parsed_reach(SYNTHETIC)
    assert reach == {"M1": 3, "M2": 5, "M3": 9, "M4": 15,
                     "M5": 14, "M6": 7, "M7": 0, "M8": 1}
    assert sorted(parsed_evaluated(SYNTHETIC)) == ["M4", "M5"]
    assert parsed_degree(SYNTHETIC) == "0.606 [0.544, 0.669]"


def test_the_degree_line_is_not_read_as_a_channel():
    """`degree` carries an interval in the same shape. Picking it up as a
    channel would put a ninth entry in a set compared for equality."""
    assert "degree" not in parsed_reach(SYNTHETIC)
    assert not any(k.startswith("deg") for k in parsed_evaluated(SYNTHETIC))


@pytest.mark.parametrize("bad", [
    "    evaluated   M4  0.507 [0.447, 0.568]   reaches 15 of 89 actives",
    "    not         M3  reaches  9",
])
def test_a_partial_block_does_not_silently_parse_as_complete(bad):
    """Half a block must not look like a whole one: the channel-set check
    is what turns a dropped line into a failure."""
    assert set(parsed_reach(bad)) != set(CHANNELS)


def test_the_parser_raises_rather_than_returning_an_empty_degree():
    with pytest.raises(AssertionError):
        parsed_degree("    not  M3 reaches 9\n")


@pytest.mark.parametrize("scorer", NON_CHANNELS)
def test_the_block_prints_the_baseline_and_ensemble_intervals(scorer, chans):
    """COMBINED and degree, asserted like every M-channel.

    `parsed_degree` covered `degree` alone. COMBINED had a line in the
    block, a figure in it, and nothing reading either.
    """
    cell = fmt_auc(chans[scorer])
    body = block()
    assert f"{scorer} {cell}" in " ".join(body.split()), (
        f"the block does not print {scorer} as {cell}")


def test_the_gap_between_them_is_the_one_the_results_give(chans):
    """The block states how far COMBINED sits below degree, and that
    difference is a third number which can disagree with the two it comes
    from -- it read "0.002" while the two figures beside it implied
    0.004."""
    gap = chans["degree"]["auc"] - chans["COMBINED"]["auc"]
    flat = " ".join(DOC.read_text().split())
    assert f"{gap:.3f} below" in flat or f"{gap:.3f} apart" in flat, (
        f"degree leads COMBINED by {gap:.3f}; the document states a "
        f"different gap")


# --------------------------------------------------------------------------
# THE SECOND CLAIM: a structure-derived scorer, added to this section on
# 2026-11-17 by the owner's decision on fix-list A1. Everything in it is
# derived here from the three results files, for the reason the channel
# block above is: a claim block is the one part of this document a reader
# quotes without opening anything else.
# --------------------------------------------------------------------------

RESULTS_DIR = RESULTS.parent
SCORER_SYM = RESULTS_DIR / "structure_scorer_symmetric.json"
SCORER_STR = RESULTS_DIR / "structure_scorer_strength_matched.json"
CLINICAL = RESULTS_DIR / "clinical_labels.json"

_SCORER_POPS = (("primary_chembl", "ChEMBL"),
                ("secondary_ncats", "plate CPE"),
                ("secondary_ncats_entry", "plate entry"))


def _scorer_cells(path, channel: str) -> dict:
    """{(population label, virus name): (structure, degree)} for one run."""
    import json as _j
    from kgav import virus_register as _vr

    names = {"NCBITaxon:" + v["taxon"]: v["name"]
             for v in _vr.load()["viruses"]}
    doc = _j.loads(path.read_text())
    out = {}
    for pop, label in _SCORER_POPS:
        for virus, cell in (doc.get(pop) or {}).items():
            if not isinstance(cell, dict):
                continue
            s, g = cell.get(channel), cell.get("degree")
            if not (isinstance(s, dict) and "auc" in s
                    and isinstance(g, dict) and "auc" in g):
                continue
            out[(label, names.get(virus, virus))] = (s, g)
    return out


def _group(sym: dict, key) -> int:
    s, g = sym[key]
    if not s.get("clears_corrected_alpha"):
        return 3
    if g["auc"] >= s["auc"]:
        return 2
    return 0 if not g.get("clears_corrected_alpha") else 1


@pytest.fixture(scope="module")
def scorer():
    for p in (SCORER_SYM, SCORER_STR, CLINICAL):
        if not p.exists():
            pytest.skip(f"{p.name} absent")
    return (_scorer_cells(SCORER_SYM, "structure_max_k10"),
            _scorer_cells(SCORER_STR, "structure_max"))


@pytest.fixture(scope="module")
def scorer_block() -> list[str]:
    text = DOC.read_text()
    assert "BLOCK: claim_scorer" in text, "the scorer claim block is gone"
    body = text.split("either protocol. -->")[1]
    body = body.split("**1 cell")[0] if "**1 cell" in body else body
    return [ln for ln in body.splitlines() if ln.startswith("    ")]


def test_the_scorer_block_holds_every_cell_of_the_family(scorer, scorer_block):
    """NOT ALMOST EVERY CELL. The first draft printed three of the four
    groups -- 13 of 17 -- while its own comment claimed the table was
    complete, which is the selection failure this project is built
    against wearing a completeness notice."""
    sym, _ = scorer
    assert len(scorer_block) == len(sym), (
        f"the block prints {len(scorer_block)} cells and the runs hold "
        f"{len(sym)}")
    for _label, virus in sym:
        assert any(virus in ln for ln in scorer_block), (
            f"{virus} is absent from the claim block")


def _parse_scorer_rows(lines: list[str]) -> dict:
    """{(population, virus): (mark, three figures)} by fixed columns.

    MATCHING ON THE VIRUS NAME ALONE IS WRONG and the first version of
    this did it: SARS-CoV-2, SARS-CoV and MERS-CoV each appear in two
    populations, so three cells were compared against another cell's
    row and two tests passed or failed for the wrong reason. The writer
    emits fixed-width fields, so the reader uses them.
    """
    out = {}
    for ln in lines:
        mark, virus, pop = ln[4:18].strip(), ln[18:32].strip(), ln[32:44].strip()
        out[(pop, virus)] = (mark, tuple(ln[44:].split()))
    return out


def test_every_number_in_the_scorer_block_is_the_runs(scorer, scorer_block):
    sym, strength = scorer
    rows = _parse_scorer_rows(scorer_block)
    assert set(rows) == set(sym), (
        f"block keys {sorted(set(rows) - set(sym))} / "
        f"run keys {sorted(set(sym) - set(rows))}")
    for key, (s, g) in sym.items():
        want = (f"{s['auc']:.4f}", f"{strength[key][0]['auc']:.4f}",
                f"{g['auc']:.4f}")
        assert rows[key][1] == want, (key, rows[key][1], want)


def test_each_row_carries_the_verdict_its_numbers_imply(scorer, scorer_block):
    """The leading words are a verdict, so they are derived from the
    intervals rather than typed beside them."""
    sym, _ = scorer
    marks = {0: "clears, and", 1: "degree also", 2: "degree wins",
             3: "no clear"}
    rows = _parse_scorer_rows(scorer_block)
    for key in sym:
        assert rows[key][0] == marks[_group(sym, key)], (
            key, marks[_group(sym, key)], rows[key][0])


def test_the_counts_in_the_scorer_sentences_are_the_groups(scorer):
    sym, strength = scorer
    text = " ".join(DOC.read_text().split())
    counts = {g: sum(1 for k in sym if _group(sym, k) == g)
              for g in (0, 1, 2, 3)}
    assert sum(counts.values()) == len(sym)
    assert f"In **{counts[0]} of the {len(sym)}**" in text, (
        f"{counts[0]} of {len(sym)} cells clear above a non-clearing "
        f"baseline and the claim states something else")
    assert f"**{counts[2]} cell carries the opposite finding" in text


def test_the_two_controls_agree_on_the_membership_the_claim_rests_on(scorer):
    """The claim says "identical under both registered controls". That is
    a fact about the data, not a property of the protocols."""
    sym, strength = scorer
    a = {k for k in sym if _group(sym, k) == 0}
    b = {k for k in strength
         if strength[k][0].get("clears_corrected_alpha")
         and strength[k][0]["auc"] > strength[k][1]["auc"]
         and not strength[k][1].get("clears_corrected_alpha")}
    assert a == b, (
        f"the two controls no longer agree: symmetric-only {sorted(a - b)}, "
        f"strength-only {sorted(b - a)}; the claim's 'identical under both' "
        f"is false")
    text = " ".join(DOC.read_text().split())
    assert "identical under both registered controls" in text


def test_the_clinical_null_in_the_claim_block_is_the_runs(scorer):
    import json as _j
    cl = _j.loads(CLINICAL.read_text())
    text = DOC.read_text()
    for m in (cl["structure_max_k10"], cl["degree"]):
        assert f"{m['auc']:.4f}" in text
        assert (f"[{m['auc_ci_lo_corrected']:.4f}, "
                f"{m['auc_ci_hi_corrected']:.4f}]") in text
    pop = cl["population"]
    flat = " ".join(text.split())
    assert f"{pop['positives']} compounds in a trial" in flat
    assert f"{pop['negatives']:,} screened and not" in flat
    assert (f"{pop['positives_invitro_inactive']} of the "
            f"{pop['positives']} positives were measured inactive") in flat
    assert cl["structure_max_k10"]["auc_ci_hi_corrected"] < 0.5, (
        "the clinical cell is no longer wholly below chance and the claim "
        "block says it is")


# --------------------------------------------------------------------------
# THE THIRD CLAIM: a PROSPECTIVE measurement, added to this section on
# 2026-11-20. Everything above it in the document is cross-sectional --
# the labels were split by publication or by assay, never by time -- so
# this is the only block in the claim section that answers "would this
# have worked", and it is the one a reader is most likely to quote.
#
# The document carried the opposite statement until today. Its temporal
# decomposition said "No temporal cell has an AUC interval excluding 0.5
# on this release", which is true of THAT table (temporal_split.py's
# metapath scorers) and false of the project from 2026-11-19 onward.
# --------------------------------------------------------------------------

TEMPORAL = RESULTS_DIR / "temporal_structure_scorer.json"
PROSPECTIVE_MARKER = "<!-- BLOCK: claim_prospective"

#: The verdict word each group prints. `clears BELOW` is NOT a third
#: flavour of failure: it is a cell whose corrected interval excludes
#: chance on the WRONG side, which is a finding. Conflating it with "no
#: clear" is the error the agent's refusal was making until 2026-11-20.
PROSPECTIVE_MARKS = {0: "clears", 1: "no clear", 2: "clears BELOW"}


@pytest.fixture(scope="module")
def prospective() -> dict:
    if not TEMPORAL.exists():
        pytest.skip("temporal_structure_scorer.json absent")
    return json.loads(TEMPORAL.read_text())


def _prospective_group(taxon: str, cell: dict, clears: set[str]) -> int:
    p = cell["primary"]
    if taxon in clears:
        return 0
    if p.get("clears_corrected_alpha") and p.get("direction") == "below":
        return 2
    return 1


def _prospective_expected(doc: dict) -> dict:
    """{virus name: (mark, auc, ci_lo, ci_hi, n_pos, n_neg, degree)}."""
    clears = set(doc["clears"])
    out = {}
    for taxon, cell in doc["cells"].items():
        p = cell.get("primary") or {}
        if not isinstance(p.get("auc"), (int, float)):
            continue
        dg = cell.get("degree") or {}
        degree = (f"{dg['auc']:.4f}"
                  if isinstance(dg.get("auc"), (int, float)) and dg.get("usable")
                  else "--")
        out[cell["name"]] = (
            PROSPECTIVE_MARKS[_prospective_group(taxon, cell, clears)],
            f"{p['auc']:.4f}",
            f"{p['auc_ci_lo_corrected']:.4f}",
            f"{p['auc_ci_hi_corrected']:.4f}",
            str(p["n_pos"]), str(p["n_neg"]), degree)
    return out


@pytest.fixture(scope="module")
def prospective_block() -> dict:
    """The block as the document prints it: {virus: tuple of fields}."""
    text = DOC.read_text()
    assert PROSPECTIVE_MARKER in text, "the prospective claim block is gone"
    body = text.split("to produce one. -->", 1)[1]
    body = body.split("**", 1)[0]
    rows = {}
    for ln in body.splitlines():
        if not ln.startswith("    "):
            continue
        mark = ln[4:18].strip()
        virus = ln[18:32].strip()
        rest = ln[32:].replace("[", "").replace("]", "").replace(",", " ")
        auc, lo, hi, n_pos, _slash, n_neg, degree = rest.split()
        rows[virus] = (mark, auc, lo, hi, n_pos, n_neg, degree)
    assert rows, "the prospective block parsed no rows"
    return rows


def test_the_prospective_block_holds_every_evaluable_cell(
        prospective, prospective_block):
    """NOT THE SIX THAT CLEAR. A block of the clearing cells alone would be
    the selection failure this project is built against, and the two cells
    that clear BELOW chance are the ones a reader most needs."""
    expected = _prospective_expected(prospective)
    assert set(prospective_block) == set(expected), (
        f"block prints {sorted(prospective_block)}, the run evaluated "
        f"{sorted(expected)}")
    assert len(expected) == prospective["n_evaluable"], (
        "the run's own n_evaluable disagrees with the cells carrying an AUC")


@pytest.mark.parametrize("field,label", list(enumerate(
    ("verdict", "AUC", "corrected ci_lo", "corrected ci_hi",
     "n_pos", "n_neg", "degree comparator"))))
def test_every_field_in_the_prospective_block_is_the_runs(
        prospective, prospective_block, field, label):
    expected = _prospective_expected(prospective)
    wrong = [(v, prospective_block[v][field], expected[v][field])
             for v in expected
             if prospective_block[v][field] != expected[v][field]]
    assert not wrong, f"{label} differs: {wrong}"


def test_the_block_prints_the_corrected_interval_not_the_95_percent_one():
    """The registration decides at alpha 0.05/10, z = 2.8070. Printing the
    95% interval would show cells excluding chance that the registered
    family size says do not -- the document claiming more than the study
    it is reporting, out of the same file."""
    doc = json.loads(TEMPORAL.read_text())
    disagree = [(t, c["primary"]) for t, c in doc["cells"].items()
                if isinstance((c.get("primary") or {}).get("auc_ci_lo"), float)
                and c["primary"]["auc_ci_lo"]
                != c["primary"]["auc_ci_lo_corrected"]]
    assert disagree, "no cell distinguishes the two intervals on this release"
    text = DOC.read_text()
    body = text.split("to produce one. -->", 1)[1].split("**", 1)[0]
    for _taxon, p in disagree:
        assert f"{p['auc_ci_lo']:.4f}" not in body, (
            f"the block prints the uncorrected lower bound "
            f"{p['auc_ci_lo']:.4f}")


def test_the_counts_in_the_prospective_sentences_are_the_groups(prospective):
    expected = _prospective_expected(prospective)
    counts = {m: sum(1 for v in expected.values() if v[0] == m)
              for m in PROSPECTIVE_MARKS.values()}
    flat = " ".join(DOC.read_text().split())
    assert (f"**{counts['clears']} of the {len(expected)} evaluable cells "
            f"clear above chance, and {counts['clears BELOW']} clear BELOW "
            f"it**") in flat, (
        f"the sentence does not state the derived groups {counts}")
    usable = len(prospective["comparator_usable"])
    assert (f"**{usable} of the {len(expected)} cells have a usable degree "
            f"comparator.**") in flat


def test_a_cell_that_clears_below_chance_is_not_called_a_failure_to_clear(
        prospective, prospective_block):
    """The distinction the whole block turns on. A cell whose interval
    excludes chance on the wrong side ranks actives reliably LAST, which
    is information; "does not clear" says the opposite of information."""
    expected = _prospective_expected(prospective)
    below = [v for v, f in expected.items() if f[0] == "clears BELOW"]
    assert below, "no cell clears below chance, so this asserts nothing"
    for v in below:
        assert prospective_block[v][0] == "clears BELOW", v
    flat = " ".join(DOC.read_text().split())
    assert "exclude chance on the wrong side" in flat
    assert "inverting them would hold forward" in flat


# --------------------------------------------------------------------------
# AND THE QUALIFICATION ON ALL OF IT: an AUC is not a hit rate.
#
# POST-HOC, and written into the document as a qualification rather than a
# claim. scripts/hit_rate_at_k.py deliberately writes nothing into
# data/results, so this test RE-RUNS it and compares -- one implementation
# of precision@k, not a second copy living in a test.
# --------------------------------------------------------------------------

HIT_RATE_MARKER = "<!-- BLOCK: claim_hit_rate"


def test_the_hit_rate_block_is_marked_as_machine_checked():
    text = DOC.read_text()
    assert HIT_RATE_MARKER in text, "the hit-rate qualification lost its marker"
    head = text.split(HIT_RATE_MARKER, 1)[1].split("-->", 1)[0]
    # It must say what it is, because a post-hoc statistic presented
    # beside registered ones acquires their standing by adjacency.
    assert "POST-HOC" in head
    assert "not in data/results" in head


@pytest.mark.slow
def test_the_hit_rate_figures_in_the_document_are_the_scripts():
    """Re-runs scripts/hit_rate_at_k.py and compares its plate ladder
    against the three percentages the document prints.

    SLOW because it fingerprints the plate population three times, which
    is the price of not keeping a second implementation of precision@k in
    a test file.
    """
    import subprocess
    import tempfile

    script = ROOT / "scripts" / "hit_rate_at_k.py"
    if not script.exists():
        pytest.skip("hit_rate_at_k.py absent")
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "hr.json"
        r = subprocess.run([sys.executable, str(script), "--out", str(out),
                            "--k", "10"],
                           capture_output=True, text=True, cwd=ROOT,
                           timeout=1800)
        assert r.returncode == 0, r.stderr
        doc = json.loads(out.read_text())

    plate = [c for c in doc["plate_cells"]
             if c["population"] == "secondary_ncats"]
    assert len(plate) == 1, f"expected one cytopathic plate cell, got {plate}"
    cell = plate[0]

    text = " ".join(DOC.read_text().split())
    assert f"base rate {cell['base_rate']:.2%}" in text, (
        f"the document does not state the measured base rate "
        f"{cell['base_rate']:.2%}")

    k_cell = cell["k_cell"]
    for label, control in (("no exclusion", "no exclusion"),
                           ("symmetric k = 10", "k10"),
                           (f"strength-matched (k_cell {k_cell})",
                            f"k_cell={k_cell}")):
        got = cell["controls"][control]["per_k"][str(10)] \
            if str(10) in cell["controls"][control]["per_k"] \
            else cell["controls"][control]["per_k"][10]
        pct = f"{got['precision_expected']:.1%}"
        assert f"{label} {pct}" in text, (
            f"the document does not print {label} at {pct}")


@pytest.mark.slow
def test_the_qualified_cell_is_one_the_claim_block_says_clears():
    """The point of the qualification. If the cell whose head holds no
    actives were one that fails its alpha, there would be nothing to
    qualify -- the force of it is that a CLEARING cell behaves this way,
    and the AUC cannot show it."""
    sym = _scorer_cells(SCORER_SYM, "structure_max_k10")
    key = ("plate CPE", "SARS-CoV-2")
    assert key in sym, "the cytopathic plate cell is gone from the family"
    s, _g = sym[key]
    assert s.get("clears_corrected_alpha"), (
        "the qualified cell no longer clears its alpha; the paragraph's "
        "argument has changed and the wording must change with it")
    text = " ".join(DOC.read_text().split())
    assert "one of the cells that clears" in text


def test_the_scope_statement_carries_the_hit_rate_limit_without_copying_it():
    """docs/scope_and_limits.md is where a limit belongs, and it must not
    keep its own copy of the figures.

    A limit restated with its own numbers is a second thing to go stale,
    and this document has had that failure twice. So: the limit must be
    there, it must name the mechanism, and the three percentages must NOT
    be, because `claim_hit_rate` in evaluation.md is their one home.
    """
    scope = ROOT / "docs" / "scope_and_limits.md"
    text = scope.read_text()
    flat = " ".join(text.split())

    assert "An AUC above chance is not an expected hit rate" in text, (
        "the scope statement does not carry the hit-rate limit")
    assert "claim_hit_rate" in flat, "it does not point at the figures' home"
    assert "no actives at all" in flat, "it does not state the finding"

    # The figures the claim block owns, absent from THIS SECTION.
    #
    # SCOPED, AND MATCHED AS WHOLE FIGURES -- the first version of this
    # check was wrong twice over. It searched the whole file, where an
    # unrelated label-census line reads "5,200 active (20.0%)"; and it
    # tested `"0.0%" in flat`, so that 20.0% matched by substring. A test
    # that reports an unrelated figure as a duplication gets relaxed
    # until it reports nothing.
    section = text.split("An AUC above chance is not an expected hit rate",
                         1)[1].split("\n## ", 1)[0]
    doc = DOC.read_text()
    block = doc.split(HIT_RATE_MARKER, 1)[1].split("A precision that", 1)[0]
    owned = set(re.findall(r"\d+\.\d%", block))
    assert owned, "the claim block prints no percentages to protect"
    here = set(re.findall(r"\d+\.\d%", section))
    duplicated = sorted(owned & here)
    assert not duplicated, (
        f"the hit-rate section of scope_and_limits.md repeats "
        f"{duplicated} from the claim block; those have one home and this "
        f"is a second place to drift")
