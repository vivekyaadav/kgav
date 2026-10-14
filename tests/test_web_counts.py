"""The one real claim the prototype UI makes must still be true.

web/client/src/pages/Home.tsx is scrupulous about the drawing: "ILLUSTRATIVE
DATA ONLY", "hand-drawn sample, not from v0.1", and every node carrying a
description that says so. Nothing on the canvas needs checking.

It makes exactly one factual claim:

    The layer counts below are real figures from release {RELEASE}

and those five numbers are typed into the file by hand. They have already
been wrong once -- the comment above them records the previous values as
"invented, wrong by up to elevenfold, and summing to 78,586 against a real
node count of 66,007". A figure copied into a second place and corrected
once is a figure that will drift again; this is what stops it silently.

Skips when no release is on disk, like the other tests that need one.
"""
import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
HOME = ROOT / "web" / "client" / "src" / "pages" / "Home.tsx"
MANIFEST = ROOT / "data" / "releases" / "v0.1" / "MANIFEST.json"

needs_release = pytest.mark.skipif(
    not MANIFEST.exists(), reason="no assembled release on disk")
needs_ui = pytest.mark.skipif(not HOME.exists(), reason="web UI not present")

# The UI shows display names; the release counts Biolink classes. Only the
# two that differ are mapped, so a renamed class fails loudly rather than
# being quietly treated as absent.
DISPLAY_TO_CLASS = {"Virus": "OrganismTaxon", "Compound": "SmallMolecule"}


def _declared() -> dict[str, int]:
    """The counts typed into Home.tsx, by display label."""
    src = HOME.read_text()
    block = re.search(r"const filters[^=]*=\s*\[(.*?)\];", src, re.S)
    assert block, "could not find the `filters` array in Home.tsx"
    out = {}
    for label, count in re.findall(
            r'label:\s*"([A-Za-z]+)"\s*,\s*count:\s*([0-9_]+)', block.group(1)):
        out[label] = int(count.replace("_", ""))
    assert out, "found the filters array but parsed no entries from it"
    return out


@needs_ui
@needs_release
def test_the_ui_layer_counts_match_the_release_it_names():
    """'real figures from release v0.1' has to mean this release, today."""
    by_class = json.loads(MANIFEST.read_text())["nodes_by_class"]
    declared = _declared()

    wrong = []
    for label, shown in declared.items():
        cls = DISPLAY_TO_CLASS.get(label, label)
        actual = by_class.get(cls)
        if actual is None:
            wrong.append(f"{label}: no class {cls!r} in the release "
                         f"(has {sorted(by_class)})")
        elif actual != shown:
            wrong.append(f"{label} ({cls}): UI says {shown:,}, "
                         f"release has {actual:,}")
    assert not wrong, (
        "Home.tsx states these are real figures from the release:\n  "
        + "\n  ".join(wrong)
        + "\nThese five numbers are hand-typed and have drifted once before.")


@needs_ui
@needs_release
def test_the_ui_names_the_release_it_quotes():
    """A correct count against the wrong release name is still wrong."""
    src = HOME.read_text()
    m = re.search(r'const RELEASE\s*=\s*"([^"]+)"', src)
    assert m, "could not find the RELEASE constant in Home.tsx"
    assert m.group(1) == json.loads(MANIFEST.read_text()).get(
        "release", "v0.1"), "Home.tsx quotes a different release"


@needs_ui
def test_every_drawn_node_says_it_is_drawn():
    """No release needed: this is about the file, not the data. The canvas
    may show anything it likes as long as the page never lets it be read as
    retrieved -- which is the whole reason the counts above are the only
    checkable claim here."""
    src = HOME.read_text()
    assert "ILLUSTRATIVE DATA ONLY" in src
    assert "not from {RELEASE}" in src or "not from v0.1" in src
