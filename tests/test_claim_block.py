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
    for name in CHANNELS + ["degree"]:
        m = per[name]
        lo, hi = auc_interval(m["auc"], m["n_pos"], m["n_neg"])
        out[name] = {"auc": m["auc"], "lo": lo, "hi": hi,
                     "n_pos": m["n_pos"],
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


def test_degree_is_still_the_only_scorer_clear_of_chance(chans):
    """The block's own claim about the baseline. If a reasoning channel
    ever clears, this fails -- and that failure is the result, not a test
    to relax."""
    clear = [n for n in CHANNELS
             if chans[n]["lo"] > 0.5 or chans[n]["hi"] < 0.5]
    assert not clear, (
        f"{clear} now has an interval excluding chance; the block says "
        f"degree is the only scorer clear of 0.5")
    d = chans["degree"]
    assert d["lo"] > 0.5, (
        f"degree reads {fmt_auc(d)}, which contains chance -- the block "
        f"calls it the only scorer clear of 0.5")


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
