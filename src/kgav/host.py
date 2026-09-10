"""Day 4: build the host layer.

Human reviewed proteome -> Protein/Gene nodes; UniProt's own Reactome and GO
cross-references -> Pathway nodes and PARTICIPATES_IN edges; STRING -> the host
interactome.

THE STRING THRESHOLD DECISION (recorded here because it changes every downstream
path count):

    combined_score >= 700 is the conventional cut. On the v12.0 human file it
    keeps 474k edges -- but 284k of those (60%) have experimental < 400. They
    clear the bar on text-mining and coexpression. Meanwhile 120k edges with
    strong experimental support fall BELOW 700 and would be discarded.

    Text-mining co-occurrence measures how often two proteins are written about
    together, which tracks study bias, not binding. For a network-proximity
    baseline -- where the whole quantity of interest is distance through the
    interactome -- that inflates closeness around well-studied proteins and
    produces a proximity score that partly measures literature volume.

    So the filter is experimental >= 400 (~309k edges). Every channel score is
    kept as an edge qualifier, so the combined>=700 variant is reconstructible
    by query rather than re-ingest, and both can be reported as a sensitivity
    analysis.

THE IDENTIFIER PROBLEM: STRING keys on Ensembl protein ids, the graph keys on
UniProt. The alias file maps ENSP -> UniProt_AC, but one ENSP frequently maps to
SEVERAL accessions (isoforms, paralogs, obsolete entries). Picking arbitrarily
invents edges; picking none drops real ones. Resolution: keep only accessions
present in the reviewed proteome node set, and report every ENSP that stays
ambiguous or unmapped rather than silently dropping it.
"""
from __future__ import annotations

import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

EXPERIMENTAL_MIN = 400
HUMAN_TAXON = "NCBITaxon:9606"

# STRING's experimental channel aggregates evidence from IMEx, BioGRID and
# others; the per-edge detection method is not preserved in the download. The
# schema wants a PSI-MI CURIE here, so rather than assert one we did not
# observe, the source is named explicitly and the aggregate is flagged.
STRING_DETECTION = "string:experimental_channel_aggregate"


def _sha(s: str) -> str:
    return hashlib.sha256((s or "").encode()).hexdigest()[:32]


def _open(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else open(path, encoding="utf-8")


# ------------------------------------------------------------------ proteome
def load_human_proteome(path: Path) -> list[dict]:
    raw = path.read_bytes()
    if path.suffix == ".gz":
        raw = gzip.decompress(raw)
    return json.loads(raw).get("results", [])


def ingest_proteome(em, entries: list[dict], prov: dict) -> Counter:
    """Protein and Gene nodes, plus Reactome/GO pathways from UniProt xrefs."""
    stats = Counter()

    for e in entries:
        acc = e["primaryAccession"]
        pcurie = f"UniProtKB:{acc}"
        date = (e.get("entryAudit") or {}).get("firstPublicDate", "1970-01-01")[:10]

        symbol = None
        for g in e.get("genes") or []:
            if (g.get("geneName") or {}).get("value"):
                symbol = g["geneName"]["value"]
                break

        name = ((e.get("proteinDescription") or {}).get("recommendedName") or {}) \
            .get("fullName", {}).get("value")

        em.node(pcurie, "Protein",
                gene_symbol=symbol,
                taxon_id=HUMAN_TAXON,
                is_viral=False,
                sequence_hash=_sha((e.get("sequence") or {}).get("value", "")),
                reviewed=True,
                protein_family=name)
        stats["proteins"] += 1

        gene_id = None
        for xr in e.get("uniProtKBCrossReferences") or []:
            db = xr.get("database")

            if db == "GeneID" and gene_id is None:
                gene_id = xr.get("id")

            elif db == "Reactome":
                pw = f"REACT:{xr['id']}"
                label = next((p["value"] for p in xr.get("properties") or []
                              if p.get("key") == "PathwayName"), xr["id"])
                em.node(pw, "Pathway", label=label, source_ontology="reactome")
                em.edge(pcurie, "PARTICIPATES_IN", pw,
                        source=prov["pathway_source"], date=date, tier=1)
                stats["reactome_edges"] += 1

            elif db == "GO":
                props = {p.get("key"): p.get("value") for p in xr.get("properties") or []}
                term = props.get("GoTerm", "")
                # GO evidence codes: IEA is electronically inferred with no
                # curator review and dominates the annotation set. Stored per
                # edge so the noisy tier can be filtered by query later.
                evidence = (props.get("GoEvidenceType") or "").split(":")[0]
                # UniProt prefixes GO terms with aspect: P: process, F: function,
                # C: component. Only biological process is a pathway.
                if not term.startswith("P:"):
                    continue
                pw = xr["id"] if xr["id"].startswith("GO:") else f"GO:{xr['id']}"
                em.node(pw, "Pathway", label=term[2:], source_ontology="go")
                em.edge(pcurie, "PARTICIPATES_IN", pw,
                        source=prov["pathway_source"], date=date,
                        tier=3 if evidence == "IEA" else 2,
                        quals={"role": evidence} if evidence else None)
                stats["go_edges"] += 1

        if gene_id and symbol:
            gcurie = f"NCBIGene:{gene_id}"
            em.node(gcurie, "Gene", symbol=symbol, taxon_id=HUMAN_TAXON, is_viral=False)
            em.edge(pcurie, "ENCODED_BY", gcurie,
                    source=prov["proteome_source"], date=date, tier=1)
            stats["genes"] += 1

    return stats


# ------------------------------------------------------------------- mapping
def build_ensp_map(aliases: Path, valid_accessions: set[str], notes: Counter) -> dict[str, str]:
    """ENSP -> UniProt accession, disambiguated against the proteome node set.

    Multi-mapping is the norm, not the exception. An ENSP that still resolves to
    several reviewed accessions after filtering is reported and dropped: a wrong
    interaction edge is worse than a missing one, because it is indistinguishable
    from real evidence downstream.
    """
    candidates: dict[str, set[str]] = defaultdict(set)
    with _open(aliases) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[2] != "UniProt_AC":
                continue
            ensp, acc = parts[0], parts[1]
            if acc in valid_accessions:
                candidates[ensp].add(acc)

    mapping: dict[str, str] = {}
    for ensp, accs in candidates.items():
        if len(accs) == 1:
            mapping[ensp] = next(iter(accs))
        else:
            notes["ensp_ambiguous_after_filter"] += 1
    return mapping


# -------------------------------------------------------------------- STRING
def ingest_string(em, links: Path, mapping: dict[str, str], prov: dict,
                  experimental_min: int = EXPERIMENTAL_MIN) -> Counter:
    """Stream the links file; emit deduped, UniProt-keyed interaction edges."""
    stats = Counter()
    seen: set[tuple[str, str]] = set()

    with _open(links) as fh:
        header = fh.readline().split()
        col = {name: i for i, name in enumerate(header)}
        i_exp, i_comb = col["experimental"], col["combined_score"]
        channels = ["neighborhood", "fusion", "cooccurence", "coexpression",
                    "experimental", "database", "textmining"]

        for line in fh:
            f = line.split()
            stats["rows"] += 1
            if int(f[i_exp]) < experimental_min:
                continue
            stats["passed_filter"] += 1

            a, b = mapping.get(f[0]), mapping.get(f[1])
            if a is None or b is None:
                stats["unmapped_endpoint"] += 1
                continue
            if a == b:
                stats["self_loop"] += 1
                continue

            # STRING lists every pair twice (A-B and B-A). Canonical ordering
            # halves the edge count; without it every interaction is double
            # counted and degree is inflated twofold.
            key = (a, b) if a < b else (b, a)
            if key in seen:
                stats["reciprocal_duplicate"] += 1
                continue
            seen.add(key)

            quals = {"detection_method": STRING_DETECTION,
                     "string_score": float(f[i_comb]),
                     "channel": ",".join(f"{c}={f[col[c]]}" for c in channels if int(f[col[c]]) > 0)}
            em.edge(f"UniProtKB:{key[0]}", "PHYSICALLY_INTERACTS_WITH", f"UniProtKB:{key[1]}",
                    source=prov["ppi_source"], date=prov["string_release_date"],
                    tier=1, quals=quals)
            stats["edges"] += 1

    return stats
