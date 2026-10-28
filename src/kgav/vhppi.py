"""Day 5: virus-host physical interactions from PSI-MITAB 2.7.

This is the layer that makes the graph antiviral rather than generic. It is the
only bridge between the viral spine and the host interactome, so metapath M3
(drug -> host protein -> viral protein -> virus) exists entirely because of it.

WHAT THE EVIDENCE ACTUALLY SAYS. The SARS-CoV-2 detection-method histogram is
dominated by proximity-dependent biotin identification (MI:1314) and affinity
capture (MI:0400) -- roughly 6,600 of ~7,000 rows. BioID labels anything within
about 10 nm, so those edges mean "in the same complex neighbourhood", not
"binds". Binary methods (two-hybrid, pull-down) are a rounding error by
comparison.

That distinction is preserved rather than averaged away:

    detection_method  PSI-MI code from column 7 -- the assay
    interaction_type  PSI-MI code from column 12 -- MI:0407 direct interaction
                      vs MI:0915 physical association (co-complex)
    source_score      VirHostNet miscore from column 15

All of these edges are evidence_tier 1. Tier means HOW THE FACT WAS ESTABLISHED
(curated experimental / curated inferred / computational / text-mined), not how
good the assay was. BioID is a curated experimental result; demoting it to tier
2 would put it alongside GO's curated-inferred annotations, which it has nothing
in common with, and would make Day 12's tier ablation measure an incoherent
mixture. Filtering to direct interactions is a query against interaction_type,
not a re-ingest.

Consequence to carry into the write-up: paths through this layer mean the drug
target sits in the same complex neighbourhood as a viral protein. That is
weaker than binding, and "physically interacts with" should not be allowed to
imply more than the assay supports.
"""
from __future__ import annotations

import json
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

MI_RE = re.compile(r'psi-mi:"?(MI:\d+)"?')
PRO_RE = re.compile(r"(PRO_\d+)")
ACC_RE = re.compile(r"uniprotkb:([A-Z0-9]+)(?:-PRO_\d+)?", re.IGNORECASE)
TAXID_RE = re.compile(r"taxid:(-?\d+)")
PMID_RE = re.compile(r"pubmed:(\d+)")
SCORE_RE = re.compile(r"miscore:([0-9.]+)")
# MITAB column 12, "Source database(s)": psi-mi:"MI:0471"(MINT). The name
# in parentheses is the curator's own, so the label is DERIVED from the
# data rather than from a MI-code table that would drift.
SOURCE_DB_RE = re.compile(r'psi-mi:"?MI:\d+"?\(([^)]+)\)')

# MI:0407 "direct interaction" asserts a binary contact. MI:0915 "physical
# association" and MI:0914 "association" only assert co-membership of a complex.
DIRECT_INTERACTION = "MI:0407"
# Assays that read out a binary contact rather than co-purification/proximity.
BINARY_METHODS = {"MI:0018", "MI:0096", "MI:0900", "MI:0397", "MI:0055"}


def _first(pattern: re.Pattern, text: str) -> str | None:
    m = pattern.search(text or "")
    return m.group(1) if m else None


def parse_identifier(col: str) -> str | None:
    """MITAB identifier column -> a graph node id.

    A chain reference like uniprotkb:P0DTD1-PRO_0000449624 resolves to the
    CHAIN, not the parent: the source is telling you which mature peptide it
    measured, and collapsing to the polyprotein throws that away.
    """
    if not col or col == "-":
        return None
    chain = _first(PRO_RE, col)
    if chain:
        return f"UniProtKB:{chain}"
    acc = _first(ACC_RE, col)
    return f"UniProtKB:{acc}" if acc else None


def parse_row(line: str) -> dict | None:
    f = line.rstrip("\n").split("\t")
    if len(f) < 15:
        return None
    return {
        "a": parse_identifier(f[0]),
        "b": parse_identifier(f[1]),
        "detection_method": _first(MI_RE, f[6]),
        "pmid": _first(PMID_RE, f[8]),
        "taxid_a": _first(TAXID_RE, f[9]),
        "taxid_b": _first(TAXID_RE, f[10]),
        "interaction_type": _first(MI_RE, f[11]),
        # WHO CURATED THIS ROW, which is not always who served the file.
        # IntAct's release carries IMEx partner records: of 12,496 viral
        # rows, 1,059 are MINT, 110 DIP, and two others -- 9.4% of the
        # layer attributing someone else's curation to IntAct.
        "source_db": _first(SOURCE_DB_RE, f[12]) if len(f) > 12 else None,
        "score": _first(SCORE_RE, f[14]),
    }


def source_infores(name: str | None) -> str | None:
    """A curator's self-reported name -> an infores id.

    Slugged from the name the file gives, not mapped from the MI code.
    A code table is a second place for the truth to live and drift, and
    the six names in this corpus -- virhostnet, IntAct, MINT, DIP,
    HPIDb, matrixdb -- all slug cleanly and unambiguously.
    """
    slug_name = re.sub(r"[^a-z0-9]+", "", (name or "").lower())
    return f"infores:{slug_name}" if slug_name else None


@dataclass(frozen=True)
class RowVerdict:
    """Where one MITAB row dies, or that it survives to be emitted.

    WHY THIS IS A FUNCTION RATHER THAN A COMMENT. The gates run in a fixed
    order and a row is counted at the FIRST one it fails, so a reader who
    wants to know what `unresolved_id` contains cannot simply re-test
    resolution: rows that would also fail it may have died at `unparseable`
    or `no_detection_method` first and been counted there.

    diagnose_flavivirus_join.py already carries a paragraph about this --
    "ingest_mitab checks the taxon map LAST, after unparseable,
    no_detection_method, unresolved_id, self_interaction and viral_viral" --
    written after an earlier version of that script asserted rows were lost
    at a gate they never reached. A second implementation of this order is
    how that claim became wrong; there is one implementation now.

    `gate` is None for a row that survives every gate here. Deduplication is
    NOT a gate: it needs state across rows, so ingest_mitab keeps it.
    """
    gate: str | None
    row: dict | None = None
    a: str | None = None
    b: str | None = None
    host: str | None = None
    viral: str | None = None
    viral_tax: str | None = None
    #: For an unresolved_id row: which side(s) failed to resolve, and the
    #: taxid column of each, verbatim from the file.
    unresolved: tuple[tuple[str, str | None, str | None], ...] = ()
    #: The accession the experiment used, when ppi_accession_map attributed
    #: this row's VIRAL side to a different one. None when nothing was
    #: rewritten, so an edge without it was measured on the entry it names.
    measured: str | None = None


def classify_row(line: str, index, taxon_map: dict[str, str],
                 accession_map: dict[str, str] | None = None) -> RowVerdict:
    """Run one MITAB line through the ingest's gates, in the ingest's order.

    accession_map attributes a measured accession to the one this build
    carries -- see virus_register.ppi_accession_map. It is applied BEFORE
    resolution, because the whole point is that the measured accession does
    not resolve, and the substitution is recorded on the verdict so the
    edge can carry it rather than losing what was measured.
    """
    r = parse_row(line)
    if not r or not r["a"] or not r["b"]:
        return RowVerdict("unparseable", r)
    if not r["detection_method"]:
        return RowVerdict("no_detection_method", r)

    # UNCONDITIONAL, NOT A FALLBACK. Applying it only when resolution fails
    # would make the graph's shape depend on which proteomes happened to be
    # loaded: the same row would attribute to the reference in one build and
    # to a strain entry in another, and nothing would say which.
    measured = {}
    if accession_map:
        for side in ("a", "b"):
            acc = (r[side] or "").removeprefix("UniProtKB:")
            dst = accession_map.get(acc)
            if dst and dst != acc:
                measured[side] = acc
                r[side] = f"UniProtKB:{dst}"

    a, b = index.resolve(r["a"]), index.resolve(r["b"])
    if a is None or b is None:
        # Recorded per side, because "which endpoint failed" is the question
        # the counter alone cannot answer: an unresolved HOST id is lost
        # interactome, an unresolved VIRAL id is a virus this build does not
        # carry, and the single tally reads the same either way.
        bad = tuple((raw, tax, side)
                    for raw, tax, side, got in
                    ((r["a"], r["taxid_a"], "a", a), (r["b"], r["taxid_b"], "b", b))
                    if got is None)
        return RowVerdict("unresolved_id", r, a, b, unresolved=bad)
    if a == b:
        return RowVerdict("self_interaction", r, a, b)

    va, vb = index.is_viral(a), index.is_viral(b)
    if va is None or vb is None:
        return RowVerdict("not_a_protein", r, a, b)
    if va == vb:
        return RowVerdict("viral_viral" if va else "host_host", r, a, b)

    host, viral = (b, a) if va else (a, b)
    viral_tax = r["taxid_a"] if va else r["taxid_b"]
    if viral_tax not in taxon_map:
        return RowVerdict("unmapped_viral_taxon", r, a, b, host, viral, viral_tax)
    return RowVerdict(None, r, a, b, host, viral, viral_tax,
                      measured=measured.get("a" if va else "b"))


def slug(name: str) -> str:
    """A register name -> the stem its raw file carries.

    ONE IMPLEMENTATION. fetch_vhppi writes vhn_<slug>.tab27 and the ingest
    has to read the virus back off the filename to say which virus a file
    belongs to. Two copies of this rule drift, and this repo has already
    paid for that: five hand-built copies of the taxon map, two of the
    trustworthy rule, two of the gate order.
    """
    return re.sub(r"[^a-z0-9]+", "", name.lower())


#: Filename prefixes that predate the sidecar, and the service each one
#: came from. Declared rather than derived: these files were placed by
#: hand and their names are the only record of where they came from, so
#: the mapping is a curation claim and belongs where it can be read.
LEGACY_PREFIX_SOURCE = {
    "vhn": "infores:virhostnet",
    "intact": "infores:intact",
}


def file_source(path) -> str | None:
    """Which service produced this raw file, or None if nothing says.

    THE SIDECAR IS AUTHORITATIVE. ingest_vhppi used to decide this from a
    glob -- vhn_*.tab27 meant VirHostNet -- which is fine while one
    service is ever fetched and silently wrong the moment another is:
    fetch_vhppi takes --base and will happily query MINT, and its file
    would have landed under VirHostNet's name. The tier ablation reads
    these labels, so a mislabelled row is not a cosmetic error.

    None rather than a guess. A file whose provenance nothing records is
    refused by the caller; inventing a source here is the failure this
    function exists to stop.
    """
    side = Path(str(path) + ".source.json")
    if side.exists():
        try:
            src = (json.loads(side.read_text()) or {}).get("source")
        except (ValueError, OSError):
            src = None
        if src:
            return str(src)
    stem = Path(path).stem
    return LEGACY_PREFIX_SOURCE.get(stem.split("_", 1)[0])


def file_taxon(path, by_stem: dict[str, str],
               taxon_map: dict[str, str] | None = None) -> str | None:
    """Which virus a raw MITAB file was fetched for, or None.

    ONE ATTRIBUTION RULE, because two callers need it and they must not
    disagree: the ingest, to say which virus lost a row that died before
    the taxon map, and diagnose_unresolved_ids, to narrow itself to that
    virus. A diagnostic that attributed files differently from the report
    that sent you to it would be worse than no diagnostic.

    The sidecar fetch_vhppi writes names the virus outright and is read
    first, because it survives a rename. The filename is the fallback for
    the hand-placed coronavirus files that predate the fetcher.

    Returns None for a file that covers no single virus -- an
    intact_*.tab27 holds many -- so the caller can report it unattributed
    rather than pin it on one by guesswork.
    """
    taxon = None
    side = Path(str(path) + ".source.json")
    if side.exists():
        try:
            taxon = str(json.loads(side.read_text()).get("taxon") or "")
        except (ValueError, OSError):
            taxon = None
    if not taxon:
        taxon = by_stem.get(Path(path).stem.removeprefix("vhn_"))
    if not taxon:
        return None
    return (taxon_map or {}).get(taxon, taxon)


def ingest_mitab(em, paths: list[Path], index, taxon_map: dict[str, str],
                 source: str, default_date: str,
                 per_file: dict | None = None,
                 accession_map: dict[str, str] | None = None) -> Counter:
    """Emit host<->viral PHYSICALLY_INTERACTS_WITH edges.

    taxon_map maps every viral taxid that may appear (species AND isolate) to
    the species taxid the graph uses. Swiss-Prot and VirHostNet both file
    several coronaviruses under isolate taxa; without the map those rows look
    like an unknown organism and vanish.

    per_file, if given, is filled with {Path: Counter} carrying the same keys
    per file. A row that dies at `unresolved_id` has no taxon -- resolution
    runs before the taxon map -- so the ONLY thing that can say which virus
    it belonged to is the file it came from, which is why the per-file tally
    exists and a per-taxon one cannot replace it. Deduplication is counted
    against the file holding the later copy.
    """
    stats: Counter = Counter()
    seen: set[tuple] = set()

    for path in paths:
        fstats: Counter = (per_file.setdefault(Path(path), Counter())
                           if per_file is not None else Counter())
        for line in Path(path).read_text(errors="replace").splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            stats["rows"] += 1
            fstats["rows"] += 1
            # Every gate, in order, lives in classify_row -- including the
            # ones whose comments used to live here: detection_method is a
            # required qualifier, so an edge without it is dropped rather
            # than defaulted; and viral-viral and host-host pairs are real
            # but out of scope for a BRIDGE layer (host-host belongs to the
            # STRING pass).
            v = classify_row(line, index, taxon_map, accession_map)
            if v.gate:
                stats[v.gate] += 1
                fstats[v.gate] += 1
                continue
            r, host, viral = v.row, v.host, v.viral

            key = (host, viral, r["detection_method"], r["pmid"])
            if key in seen:
                stats["duplicate"] += 1
                fstats["duplicate"] += 1
                continue
            seen.add(key)

            # No host_taxon qualifier: it would duplicate taxon_id on the
            # host protein node, and a fact stored twice is a fact that can
            # disagree with itself.
            quals = {"detection_method": r["detection_method"]}
            if v.measured:
                # WHAT WAS MEASURED, not what it was attributed to. Without
                # this the edge is indistinguishable from one measured on
                # the reference entry, and the attribution becomes an
                # unrecorded claim -- exactly the silent identity this
                # layer refused to make with SAME_AS.
                quals["measured_accession"] = v.measured
                stats["attributed"] += 1
                fstats["attributed"] += 1
            # interaction_type is NOT stored. VirHostNet normalises every row
            # to MI:0915 "physical association" regardless of assay -- all
            # 8,252 SARS-CoV-2 rows, including all 153 two-hybrid ones. A
            # qualifier with one value across every edge is noise that looks
            # like signal. The assay distinction survives in detection_method.

            # THE CURATOR, NOT THE SERVER. `source` is the service the
            # file came from, which is retrieval provenance and lives in
            # the sidecar. primary_knowledge_source is a claim about who
            # established the fact, and the row says so itself. A row
            # with no readable source column falls back to the file's
            # service and is counted, because a silent fallback is how
            # 1,171 MINT and DIP rows came to be filed as IntAct.
            curator = source_infores(r.get("source_db")) or source
            stats[f"src:{curator}"] += 1
            fstats[f"src:{curator}"] += 1
            if not r.get("source_db"):
                stats["source_db_absent"] += 1
            em.edge(host, "PHYSICALLY_INTERACTS_WITH", viral,
                    source=curator, date=default_date, tier=1, quals=quals,
                    pmids=[f"PMID:{r['pmid']}"] if r["pmid"] else None,
                    score=float(r["score"]) if r["score"] else None)
            stats["edges"] += 1
            fstats["edges"] += 1
            stats["direct" if r["interaction_type"] == DIRECT_INTERACTION else "co_complex"] += 1
            if r["detection_method"] in BINARY_METHODS:
                stats["binary_method"] += 1

    return stats
