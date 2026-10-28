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

        # A DECLARED MAPPING MUST POINT SOMEWHERE THIS BUILD CARRIES.
        # ppi_accession_map attributes an interaction measured against one
        # strain's polyprotein to the polyprotein this build holds. If the
        # target is not an accession this entry brings in, the mapping is a
        # claim that silently resolves to nothing and the rows stay lost
        # while the register says they were handled.
        pmap = v.get("ppi_accession_map") or {}
        if pmap:
            carried = set(v.get("extra_accessions") or [])
            if cs.get("accession"):
                carried.add(cs["accession"])
            for src_acc, dst in sorted(pmap.items()):
                if src_acc == dst:
                    out.append(RegisterViolation(
                        "PPI_MAP_SELF",
                        f"ppi_accession_map maps {src_acc} to itself, which "
                        f"records nothing", ref))
                if dst not in carried:
                    out.append(RegisterViolation(
                        "PPI_MAP_TARGET_NOT_CARRIED",
                        f"ppi_accession_map sends {src_acc} to {dst}, which "
                        f"this entry does not carry (chain_source.accession "
                        f"and extra_accessions give {sorted(carried) or '[]'})"
                        f". The mapping would resolve to nothing", ref))
                if src_acc in carried:
                    out.append(RegisterViolation(
                        "PPI_MAP_SOURCE_IS_CARRIED",
                        f"ppi_accession_map rewrites {src_acc}, but this entry "
                        f"carries that accession as a node of its own. "
                        f"Rewriting it would discard the real protein", ref))
            if "ppi_accession_note" not in v:
                out.append(RegisterViolation(
                    "PPI_MAP_UNEXPLAINED",
                    "ppi_accession_map must carry ppi_accession_note saying on "
                    "what basis the accessions were judged the same gene "
                    "product; an undocumented identity claim is the thing this "
                    "register exists to prevent", ref))

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


# Below this length a name is not offered for question matching. The agent's
# resolver skips short strings to avoid spurious substring hits, and a
# register entry that could never be matched is worse than absent: it reads
# as support the agent does not have. "wnv", "yfv" and "jev" are therefore
# NOT declared as aliases -- they would be inert.
MIN_SEARCH_NAME = 4


def family_of(entry: dict, doc: dict) -> str | None:
    """A virus's family, with the register's defaults applied.

    THE RULE HAS FIVE COPIES AND THREE WERE WRONG. Coronaviruses carry
    no per-entry `family`: they inherit it from defaults, which exists
    precisely so the shared Coronaviridae properties are stated once.
    Nine flaviviruses state theirs because they differ from the default.

    So `entry.get("family")` returns None for every coronavirus, and a
    filter written that way covers flaviviruses only. fetch_vhppi and
    diagnose_unresolved_ids were both written that way this week:
    `--family Coronaviridae` matched nothing and exited 2, which reads
    as "the register has no coronaviruses".

    ingest_spine already had it right -- `cfg.get("family") or
    defaults["family"]`, with a comment about the surveyed: block -- and
    so did three diagnostics. Four correct copies did not stop the fifth
    being wrong, which is the argument for there being one.
    """
    return entry.get("family") or (doc.get("defaults") or {}).get("family")


def families(doc: dict, section: str = "viruses") -> dict[str, str]:
    """{taxon: family} for every entry, defaults applied."""
    return {str(v["taxon"]): family_of(v, doc) or ""
            for v in doc.get(section) or []}


def ppi_accession_map(doc: dict, section: str = "viruses") -> dict[str, str]:
    """Measured accession -> the accession this build attributes it to.

    VirHostNet curates some flavivirus interactions against an unreviewed,
    strain-specific polyprotein entry rather than the reviewed polyprotein
    this build carries. The protein is the same gene product of the same
    species -- identical length, same organism -- but it is NOT the same
    UniProt entry, and SAME_AS would say it is.

    So this is an ATTRIBUTION, not an identity: the edge lands on the
    polyprotein this build holds and records the accession actually
    measured in a `measured_accession` qualifier. Nothing in the graph
    claims the two entries are one.

    Keys are bare accessions as the register states them; callers that work
    in node ids prefix them. Collisions across viruses raise, because one
    accession attributed to two polyproteins is a curation error, not a
    precedence question.
    """
    out: dict[str, str] = {}
    owner: dict[str, str] = {}
    for entry in doc.get(section) or []:
        for src, dst in (entry.get("ppi_accession_map") or {}).items():
            if src in out and out[src] != dst:
                raise ValueError(
                    f"{src} is mapped to {out[src]} by {owner[src]} and to "
                    f"{dst} by {entry.get('name')}")
            out[src] = dst
            owner[src] = entry.get("name") or "<unnamed>"
    return out


def search_names(doc: dict, section: str = "viruses") -> dict[str, str]:
    """Every string a person might type for a virus -> its taxon CURIE.

    WHY THE AGENT NEEDS THIS. Entity resolution ran off the release's
    NAMES.json, which carries the register's short `name` and nothing else.
    So "dengue virus 2" matched nothing, no taxon resolved, and the question
    was silently answered about the default virus -- SARS-CoV-2 -- under a
    caveat asserting the question had named no virus. It had. Every virus
    in this register whose common name is not its short form was affected,
    which on this release is every flavivirus.

    `full_name` was already curated here and simply unused.

    Collisions raise rather than resolving to whichever entry happened to be
    last: a string that two taxa both claim is a curation error, and silently
    picking one would reintroduce the same class of wrong-virus answer this
    function exists to end.
    """
    out: dict[str, str] = {}
    for entry in doc.get(section) or []:
        taxon = str(entry.get("taxon") or "")
        if not taxon:
            continue
        curie = f"NCBITaxon:{taxon}"
        candidates = [entry.get("name"), entry.get("full_name")]
        candidates += list(entry.get("aliases") or [])
        for raw in candidates:
            if not raw:
                continue
            key = " ".join(str(raw).lower().split())
            if len(key) < MIN_SEARCH_NAME:
                continue
            if key in out and out[key] != curie:
                raise ValueError(
                    f"register name {key!r} is claimed by both {out[key]} "
                    f"and {curie}; a question using it could not be answered "
                    f"about either without guessing")
            out[key] = curie
    return out
