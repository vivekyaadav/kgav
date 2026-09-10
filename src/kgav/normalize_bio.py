"""Biological identifier normalisation: proteins, taxa, diseases.

Unlike chemicals, these are lookup problems, not algorithmic ones. The failure
mode is not a wrong answer — it is a SILENT one: a source cites a secondary
UniProt accession, the lookup misses, and the edge vanishes with no error.
Every function here therefore reports WHY it failed, never returns None bare.
"""
from __future__ import annotations

import re
import tarfile
from dataclasses import dataclass
from pathlib import Path

UNIPROT_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)


@dataclass(frozen=True)
class Resolved:
    curie: str | None
    status: str          # primary | secondary | demerged | deleted | unmapped | malformed
    original: str
    note: str = ""

    @property
    def ok(self) -> bool:
        return self.curie is not None


# --------------------------------------------------------------------- protein
class UniProtResolver:
    """Resolve UniProt accessions against the secondary and deleted lists.

    sec_ac.txt   secondary -> primary  (accession was merged or renamed)
    delac_sp.txt deleted accessions    (entry withdrawn entirely)

    A demerged accession maps to SEVERAL primaries. Picking one silently is
    wrong; we report it and let the caller decide, because the right choice
    depends on which isoform the source meant.
    """

    def __init__(self, sec_ac: Path | None = None, delac: Path | None = None):
        self.secondary: dict[str, list[str]] = {}
        self.deleted: set[str] = set()
        if sec_ac and Path(sec_ac).exists():
            self._load_secondary(Path(sec_ac))
        if delac and Path(delac).exists():
            self._load_deleted(Path(delac))

    def _load_secondary(self, path: Path) -> None:
        started = False
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not started:
                    if line.startswith("Secondary AC") or set(line.strip()) == {"_"}:
                        started = True
                    continue
                parts = line.split()
                if len(parts) == 2 and UNIPROT_RE.match(parts[0]):
                    self.secondary.setdefault(parts[0], []).append(parts[1])

    def _load_deleted(self, path: Path) -> None:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                acc = line.strip()
                if UNIPROT_RE.match(acc):
                    self.deleted.add(acc)

    def resolve(self, accession: str) -> Resolved:
        acc = (accession or "").strip().upper()
        acc = acc.split(":")[-1]
        base = acc.split("-")[0]          # strip isoform suffix P12345-2

        if not UNIPROT_RE.match(base):
            return Resolved(None, "malformed", accession, "not a UniProt accession pattern")

        if base in self.deleted:
            return Resolved(None, "deleted", accession, "entry withdrawn from UniProt")

        if base in self.secondary:
            primaries = self.secondary[base]
            if len(primaries) > 1:
                return Resolved(None, "demerged", accession,
                                f"maps to {len(primaries)} primaries: {','.join(primaries)}")
            return Resolved(f"UniProtKB:{primaries[0]}", "secondary", accession,
                            f"secondary accession -> {primaries[0]}")

        return Resolved(f"UniProtKB:{base}", "primary", accession)


# ----------------------------------------------------------------------- taxon
class TaxonResolver:
    """Resolve NCBI taxon ids, following merged.dmp so old ids still land.

    Reads directly from taxdump.tar.gz — no extraction step, no stale
    intermediate files on disk.
    """

    def __init__(self, taxdump: Path | None = None):
        self.merged: dict[str, str] = {}
        self.names: dict[str, str] = {}
        self.deleted: set[str] = set()
        self.valid: set[str] = set()
        if taxdump and Path(taxdump).exists():
            self._load(Path(taxdump))

    def _load(self, path: Path) -> None:
        with tarfile.open(path, "r:gz") as tar:
            for member, handler in (("merged.dmp", self._merged),
                                    ("delnodes.dmp", self._delnodes),
                                    ("nodes.dmp", self._nodes),
                                    ("names.dmp", self._names)):
                try:
                    fh = tar.extractfile(member)
                except KeyError:
                    continue
                if fh is not None:
                    handler(fh)

    def _merged(self, fh) -> None:
        for raw in fh:
            parts = raw.decode().split("\t|\t")
            if len(parts) >= 2:
                self.merged[parts[0].strip()] = parts[1].replace("\t|", "").strip()

    def _delnodes(self, fh) -> None:
        for raw in fh:
            self.deleted.add(raw.decode().replace("\t|", "").strip())

    def _nodes(self, fh) -> None:
        for raw in fh:
            self.valid.add(raw.decode().split("\t|\t", 1)[0].strip())

    def _names(self, fh) -> None:
        for raw in fh:
            parts = raw.decode().split("\t|\t")
            if len(parts) >= 4 and "scientific name" in parts[3]:
                self.names[parts[0].strip()] = parts[1].strip()

    def resolve(self, taxid: str | int) -> Resolved:
        t = str(taxid).strip().split(":")[-1]
        if not t.isdigit():
            return Resolved(None, "malformed", str(taxid), "taxon id must be numeric")
        if t in self.deleted:
            return Resolved(None, "deleted", str(taxid), "node deleted from NCBI Taxonomy")
        if t in self.merged:
            new = self.merged[t]
            return Resolved(f"NCBITaxon:{new}", "secondary", str(taxid),
                            f"merged into {new} ({self.names.get(new, '?')})")
        if self.valid and t not in self.valid:
            return Resolved(None, "unmapped", str(taxid), "taxon id not present in NCBI Taxonomy")
        return Resolved(f"NCBITaxon:{t}", "primary", str(taxid), self.names.get(t, ""))


# --------------------------------------------------------------------- disease
class MondoResolver:
    """Map a cross-reference to a MONDO term.

    XREFS ONLY. Never label string matching — 'influenza' matches the disease,
    the virus genus, and a handful of drug trade names, and a label-matched
    graph is one where nobody can tell which.
    """

    def __init__(self, mondo_json: Path | None = None):
        self.xref_to_mondo: dict[str, list[str]] = {}
        self.obsolete: set[str] = set()
        self.labels: dict[str, str] = {}
        if mondo_json and Path(mondo_json).exists():
            self._load(Path(mondo_json))

    def _load(self, path: Path) -> None:
        import json
        with open(path, encoding="utf-8") as fh:
            doc = json.load(fh)
        for graph in doc.get("graphs", []):
            for node in graph.get("nodes", []):
                nid = node.get("id", "")
                if "MONDO_" not in nid:
                    continue
                curie = "MONDO:" + nid.rsplit("MONDO_", 1)[1]
                meta = node.get("meta") or {}
                if meta.get("deprecated"):
                    self.obsolete.add(curie)
                    continue
                if node.get("lbl"):
                    self.labels[curie] = node["lbl"]
                for xref in meta.get("xrefs") or []:
                    val = (xref.get("val") or "").strip()
                    if val:
                        self.xref_to_mondo.setdefault(val.upper(), []).append(curie)

    def resolve(self, xref: str) -> Resolved:
        x = (xref or "").strip().upper()
        if not x:
            return Resolved(None, "malformed", xref, "empty xref")
        if x.startswith("MONDO:"):
            if x in self.obsolete:
                return Resolved(None, "deleted", xref, "MONDO term is deprecated")
            return Resolved(x, "primary", xref, self.labels.get(x, ""))
        hits = self.xref_to_mondo.get(x, [])
        if not hits:
            return Resolved(None, "unmapped", xref, "no MONDO term declares this xref")
        if len(set(hits)) > 1:
            return Resolved(None, "demerged", xref,
                            f"xref claimed by {len(set(hits))} MONDO terms: {','.join(sorted(set(hits))[:4])}")
        return Resolved(hits[0], "primary", xref, self.labels.get(hits[0], ""))
