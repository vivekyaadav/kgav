"""Validation for config/viruses.yaml.

THE RULE THIS MODULE EXISTS FOR. A virus whose mature-peptide chains live
under a descendant taxon MUST declare that taxon in isolate_taxa. Nothing in
this pipeline may resolve chains by walking down a lineage.

WHY IT IS A RULE RATHER THAN A CONVENTION. UniProt's `taxonomy_id:` search
includes descendants, so querying a species silently returns whatever strain
happens to be curated underneath it. That is convenient and it is exactly the
failure this project has already paid for four times: the discarded identifier
that produces no error, no missing-data signal, just a graph that quietly
lacks the targets everything downstream is supposed to point at. MERS lost
Mpro and RdRp that way (species 1335626 is TrEMBL-only, curation is under
1263720), HKU1 the same (290028 -> 443239).

Measured on 2026-09-17, SEVEN of nine surveyed flaviviruses are curated under
a strain taxon. A lineage-walking lookup would have "worked" for all of them
and left no record of which strain it picked -- so a rebuild months later
could pick a different one and nothing would say so.

The register therefore states the fact, and this validates that the statement
is internally consistent. It does not hit the network: chain_source.verified
records when a human checked, and staleness is a curation question, not
something to paper over with a live lookup at validation time.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

# Lists in viruses.yaml that hold virus entries. `viruses` is what the spine
# ingest builds; `surveyed` is measured but not yet built. Both are validated
# by the same rules, so an entry cannot acquire a problem by being promoted.
ENTRY_LISTS = ("viruses", "surveyed")

REQUIRED_CHAIN_SOURCE = ("accession", "taxon", "chains", "verified")


@dataclass(frozen=True)
class RegisterViolation:
    code: str
    message: str
    ref: str = ""

    def __str__(self) -> str:  # pragma: no cover - display only
        return f"[{self.code}] {self.ref}: {self.message}" if self.ref else \
            f"[{self.code}] {self.message}"


def load(path: str | Path | None = None) -> dict[str, Any]:
    if path is None:
        path = Path(__file__).resolve().parents[2] / "config" / "viruses.yaml"
    with open(path, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def _entries(doc: dict) -> list[tuple[str, dict]]:
    out = []
    for key in ENTRY_LISTS:
        for entry in doc.get(key) or []:
            out.append((key, entry))
    return out


def validate(doc: dict) -> list[RegisterViolation]:
    """Every rule the register must satisfy. Empty list == valid."""
    out: list[RegisterViolation] = []
    seen_taxa: dict[str, str] = {}

    for section, v in _entries(doc):
        name = v.get("name") or "<unnamed>"
        ref = f"{section}:{name}"
        taxon = str(v.get("taxon") or "")

        if not taxon:
            out.append(RegisterViolation("NO_TAXON", "entry has no taxon", ref))
            continue
        if taxon in seen_taxa:
            out.append(RegisterViolation(
                "DUPLICATE_TAXON",
                f"taxon {taxon} is already registered to {seen_taxa[taxon]}", ref))
        seen_taxa[taxon] = ref

        isolates = [str(x) for x in (v.get("isolate_taxa") or [])]
        for iso in isolates:
            if iso == taxon:
                out.append(RegisterViolation(
                    "ISOLATE_IS_SPECIES",
                    f"isolate_taxa lists {iso}, which is the species taxon itself; "
                    f"that records nothing and hides whether a strain is in use", ref))

        cs = v.get("chain_source")
        if not isinstance(cs, dict):
            out.append(RegisterViolation(
                "NO_CHAIN_SOURCE",
                "no chain_source: the taxon whose PRO_ chains this virus uses must "
                "be stated, never discovered by descendant lookup", ref))
            continue
        missing = [k for k in REQUIRED_CHAIN_SOURCE if cs.get(k) in (None, "")]
        if missing:
            out.append(RegisterViolation(
                "INCOMPLETE_CHAIN_SOURCE",
                f"chain_source is missing {', '.join(missing)}", ref))

        chain_taxon = str(cs.get("taxon") or "")
        # THE RULE. A descendant taxon may be used, but only if it is declared.
        if chain_taxon and chain_taxon != taxon and chain_taxon not in isolates:
            out.append(RegisterViolation(
                "UNDECLARED_CHAIN_TAXON",
                f"chains come from taxon {chain_taxon} but the species taxon is "
                f"{taxon} and {chain_taxon} is not in isolate_taxa "
                f"{isolates or '[]'}. Resolving this by descendant lookup is the "
                f"MERS trap: it works, records nothing, and a rebuild may silently "
                f"pick a different strain", ref))

        chains = cs.get("chains")
        if isinstance(chains, int) and chains <= 0:
            out.append(RegisterViolation(
                "NO_CHAINS",
                f"chain_source.chains is {chains}: a polyprotein with no mature "
                f"peptides gives the direct-acting channel nothing to point at", ref))

        if "mondo" not in v:
            out.append(RegisterViolation(
                "MONDO_UNSTATED",
                "mondo must be present, explicitly null if not curated; absent is "
                "indistinguishable from forgotten", ref))
    return out


def declared_taxa(doc: dict, section: str = "viruses") -> dict[str, str]:
    """Every taxid that may appear in a source -> the species taxon to use.

    The same shape the ingests build ad hoc from isolate_taxa, derived here so
    chain_source cannot drift out of the map.
    """
    out: dict[str, str] = {}
    for entry in doc.get(section) or []:
        taxon = str(entry.get("taxon"))
        out[taxon] = taxon
        for iso in entry.get("isolate_taxa") or []:
            out[str(iso)] = taxon
        cs = entry.get("chain_source") or {}
        if cs.get("taxon"):
            out[str(cs["taxon"])] = taxon
    return out
