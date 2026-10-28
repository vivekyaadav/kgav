"""Gate for the ChEMBL join diagnostic.

It answers the one question the register cannot: where ChEMBL actually
filed the data, as opposed to where UniProt curates it. It read the
`surveyed:` block only -- and the nine flaviviruses it was written for
were promoted into `viruses:` and the block emptied, so for weeks it
exited with "no `surveyed:` entries" and answered nothing.
"""
from __future__ import annotations

import sys
from importlib import import_module
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

D = import_module("diagnose_flavivirus_join")

CFG = {
    "defaults": {"family": "Coronaviridae"},
    "viruses": [
        {"name": "DENV-1", "taxon": "11053", "family": "Flaviviridae",
         "full_name": "Dengue virus type 1", "aliases": ["dengue 1"]},
        {"name": "SARS-CoV-2", "taxon": "2697049",
         "full_name": "Severe acute respiratory syndrome coronavirus 2"},
    ],
    "surveyed": [{"name": "Staged", "taxon": "999", "family": "Flaviviridae"}],
}


def _cfg(tmp_path):
    p = tmp_path / "v.yaml"
    p.write_text(yaml.safe_dump(CFG))
    return p


def test_it_reads_the_built_viruses(tmp_path):
    """The block it used to read is empty. Exiting with "no surveyed
    entries" is not an answer to "where did ChEMBL put the data"."""
    names = {e["name"] for e in D._entries(_cfg(tmp_path))}
    assert "DENV-1" in names and "SARS-CoV-2" in names


def test_a_staged_entry_is_still_covered(tmp_path):
    """surveyed: is still read, so a virus staged but not yet built
    does not fall out of the diagnostic that decides whether to build
    it."""
    assert "Staged" in {e["name"] for e in D._entries(_cfg(tmp_path))}


def test_family_uses_the_registers_defaults(tmp_path):
    """Coronaviruses inherit Coronaviridae from defaults. A filter
    reading the per-entry key alone sees None for all of them."""
    corona = D._entries(_cfg(tmp_path), "Coronaviridae")
    assert [e["name"] for e in corona] == ["SARS-CoV-2"]


def test_hints_come_from_the_register(tmp_path):
    """The hardcoded list was eight flavivirus phrases, written when
    this covered `surveyed:`. Three of the nine viruses with no
    direct-acting compound are coronaviruses it could never match."""
    h = D.name_hints(D._entries(_cfg(tmp_path)))
    assert "dengue virus type 1" in h
    assert "severe acute respiratory syndrome coronavirus 2" in h
    assert "dengue 1" in h, "aliases are where ChEMBL's naming often lands"


def test_the_genus_words_survive_as_extras(tmp_path):
    """ChEMBL files some targets under a genus with no serotype, which
    no single entry's name contains."""
    h = D.name_hints(D._entries(_cfg(tmp_path)))
    assert "dengue" in h and "flavivirus" in h and "coronavirus" in h


def test_a_very_short_name_is_not_a_hint(tmp_path):
    """A two-letter fragment LIKE-matches most of target_dictionary."""
    h = D.name_hints([{"name": "AB", "taxon": "1"}])
    assert "ab" not in h


def test_the_real_register_yields_every_virus_and_both_families():
    """16 built viruses, and hints that reach the coronaviruses."""
    cfg = ROOT / "config/viruses.yaml"
    assert len(D._entries(cfg)) == 16
    assert len(D._entries(cfg, "Flaviviridae")) == 9
    assert len(D._entries(cfg, "Coronaviridae")) == 7
    h = D.name_hints(D._entries(cfg))
    assert any("coronavirus" in x for x in h)
    assert any("dengue" in x for x in h)
