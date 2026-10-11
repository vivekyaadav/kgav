"""Day 3: build the taxonomy and viral proteome spine.

Emits OrganismTaxon, Gene (viral), Protein (viral) and Disease nodes plus the
edges between them, into a release directory as nodes.jsonl / edges.jsonl.

The one piece of real work here is POLYPROTEIN EXPANSION. Coronavirus replicase
is translated as one polyprotein and cleaved into nsp1-nsp16. The antiviral
targets that matter -- Mpro (nsp5), RdRp (nsp12), helicase (nsp13) -- are chains
INSIDE P0DTD1, not separate accessions. Without expanding them the direct-acting
channel has nothing to point at.

Two traps this handles:
  * pp1a and pp1ab share nsp1-nsp11. Deduped on (taxon, chain label), keeping
    the 1ab parent, or you get sixteen phantom duplicate nodes per virus.
  * Chain features without a PRO_ id cannot be given a stable node key, so they
    are reported and skipped rather than minted with a synthetic id that would
    change between UniProt releases.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

PRO_RE = re.compile(r"^PRO_\d+$")
# nsp11 is a short peptide with no known enzymatic role; kept for completeness.
NSP_RE = re.compile(r"\bnsp(\d+)\b", re.IGNORECASE)

NSP_NUM_RE = re.compile(r"non-structural protein\s+(\d+)\b", re.IGNORECASE)

# Alphacoronaviruses (229E, NL63) name mature peptides by FUNCTION where the
# betacoronaviruses use the nsp number. Without this table, nsp5 (Mpro) and
# nsp13 (helicase) carry no family label for those viruses, and cross-viral
# metapaths silently cannot transfer -- no error, just an evaluation that
# quietly underperforms.
FUNCTION_TO_NSP = {
    "papain-like protease": 3,
    "papain-like proteinase": 3,
    "3c-like proteinase": 5,
    "3c-like protease": 5,
    "main protease": 5,
    "3cl protease": 5,
    "rna-directed rna polymerase": 12,
    "rna-dependent rna polymerase": 12,
    "helicase": 13,
    "exoribonuclease": 14,
    "guanine-n7 methyltransferase": 14,
    "uridylate-specific endoribonuclease": 15,
    "endornase": 15,
    "2'-o-methyl transferase": 16,
    "2'-o-methyltransferase": 16,
    "viral protein genome-linked": 9,
}


def _nsp_number(label):
    """Resolve a chain description to an nsp number.

    Three routes, most explicit first: an nspN token, a "non-structural
    protein N" phrase, then a functional-name lookup. Returns None rather than
    guessing -- an unlabelled chain beats a wrongly labelled one.
    """
    low = (label or "").lower()
    m = NSP_RE.search(low)
    if m:
        return int(m.group(1))
    m = NSP_NUM_RE.search(low)
    if m:
        return int(m.group(1))
    for phrase, num in FUNCTION_TO_NSP.items():
        if phrase in low:
            return num
    return None


# FLAVIVIRUS CHAIN NAMES, from the real annotations rather than a convention.
#
# Read from rest.uniprot.org on 2026-10-11 for Q32ZE1 (Zika, a component of
# UP000054557) and P29990 (DENV-2 strain Thailand/16681/1984): thirteen Chain
# features each, with the SAME thirteen descriptions character for character.
# UniProt annotates the family from one controlled vocabulary, so this table
# covers every flavivirus in the register and there is no per-virus curation.
#
# Exact strings for the structural chains, because they collide on prefixes:
# "Protein prM" and "Peptide pr" are different chains, and a regex on "pr"
# claims both.
FLAVI_STRUCTURAL = {
    "capsid protein c": "C",
    "protein prm": "prM",
    "peptide pr": "pr",
    "small envelope protein m": "M",
    "envelope protein e": "E",
}
# "Serine protease subunit NS2B", "Serine protease NS3",
# "RNA-directed RNA polymerase/Methyltransferase NS5".
FLAVI_TOKEN_RE = re.compile(r"\bns(\d[ab]?)\b", re.IGNORECASE)
# "Non-structural protein 1", "... 2A", "... 4A", "... 4B". The [ab] group is
# the point: a bare (\d+) reads "Non-structural protein 2A" as NS2, which is a
# different chain that also exists.
FLAVI_PHRASE_RE = re.compile(r"non-structural protein (\d[ab]?)", re.IGNORECASE)


def _flavivirus_chain(label: str) -> str | None:
    """Resolve a flavivirus chain description to its own protein name.

    "Genome polyprotein" returns None deliberately: it is the parent product,
    not a mature peptide, and giving it a family key would make the whole
    polyprotein a join target competing with the chains cut from it -- which
    is the undifferentiated-node failure this naming exists to end.
    """
    low = (label or "").lower().strip()
    if low in FLAVI_STRUCTURAL:
        return FLAVI_STRUCTURAL[low]
    m = FLAVI_TOKEN_RE.search(low) or FLAVI_PHRASE_RE.search(low)
    if m:
        return f"NS{m.group(1).upper()}"
    return None


def chain_family(label: str, family: str | None) -> str | None:
    """What protein_family a chain gets, chosen by the virus's family.

    A BRANCH RATHER THAN A MERGED TABLE, on purpose. Every one of
    SARS-CoV-2's sixteen pp1ab chains carries an explicit nspN token
    ("Host translation inhibitor nsp1", "3C-like proteinase nsp5",
    "RNA-directed RNA polymerase nsp12" -- UniProt P0DTD1, read 2026-10-11),
    so the coronavirus routes are untouched here and coronavirus naming cannot
    move: non-Flaviviridae input reaches exactly the code it reached before.

    The two namespaces are kept apart because folds.py:41 refuses to link a
    cysteine protease to a serine one -- "Cys != Ser, no edge, no exceptions".
    nsp5 is a cysteine protease and NS2B-NS3 a serine protease, so a shared
    functional key would assert upstream, and silently, the equivalence that
    gate exists to deny. build_viral_target_index is keyed per taxon, so the
    two never meet.
    """
    if family == "Flaviviridae":
        return _flavivirus_chain(label)
    n = _nsp_number(label)
    return f"nsp{n}" if n else None




@dataclass
class Emit:
    """Accumulates nodes and edges, deduping node ids."""
    nodes: dict[str, dict] = field(default_factory=dict)
    edges: list[dict] = field(default_factory=list)
    notes: Counter = field(default_factory=Counter)

    def node(self, nid: str, cls: str, **props: Any) -> str:
        props = {k: v for k, v in props.items() if v is not None}
        if nid in self.nodes:
            self.nodes[nid]["properties"].update(props)
            self.notes[f"node_merged:{cls}"] += 1
        else:
            self.nodes[nid] = {"id": nid, "class": cls, "properties": props}
        return nid

    def edge(self, subj: str, pred: str, obj: str, *, source: str, date: str,
             tier: int = 1, quals: dict | None = None, pmids: list | None = None) -> None:
        e = {
            "subject": subj, "predicate": pred, "object": obj,
            "qualifiers": quals or {},
            "primary_knowledge_source": source,
            "evidence_tier": tier,
            "first_asserted_date": date,
        }
        if pmids:
            e["publications"] = pmids
        self.edges.append(e)


def _sha(s: str) -> str:
    return hashlib.sha256((s or "").encode()).hexdigest()[:32]


def _entry_date(entry: dict) -> str:
    """First public date for a UniProt entry; the spine's first_asserted_date."""
    d = (entry.get("entryAudit") or {}).get("firstPublicDate")
    return d[:10] if d else "1970-01-01"


def _rec_name(entry: dict) -> str:
    pd = entry.get("proteinDescription") or {}
    for key in ("recommendedName", "submissionNames"):
        v = pd.get(key)
        if isinstance(v, list):
            v = v[0] if v else None
        if v and v.get("fullName", {}).get("value"):
            return v["fullName"]["value"]
    return entry.get("primaryAccession", "?")


def _gene_symbol(entry: dict) -> str | None:
    genes = entry.get("genes") or []
    for g in genes:
        name = (g.get("geneName") or {}).get("value")
        if name:
            return name
        for syn in g.get("orfNames") or []:
            if syn.get("value"):
                return syn["value"]
    return None


def _chain_label(feat: dict) -> str:
    return (feat.get("description") or "").strip()


def _chain_id(feat: dict) -> str | None:
    fid = feat.get("featureId") or ""
    return fid if PRO_RE.match(fid) else None


def _chain_start(feat: dict) -> int | None:
    try:
        return int(feat["location"]["start"]["value"])
    except (KeyError, TypeError, ValueError):
        return None


def load_proteome(path: Path) -> list[dict]:
    raw = path.read_bytes()
    if path.suffix == ".gz":
        raw = gzip.decompress(raw)
    return json.loads(raw).get("results", [])


def ingest_virus(em: Emit, cfg: dict, defaults: dict, prov: dict, proteome_dir: Path) -> dict:
    """Ingest one virus: taxon, genes, proteins, mature peptides, disease."""
    tax, name = cfg["taxon"], cfg["name"]
    taxon_curie = f"NCBITaxon:{tax}"
    stats = Counter()

    # THE ENTRY'S OWN FIELDS WIN OVER THE DEFAULTS.
    #
    # `defaults` in config/viruses.yaml is family: Coronaviridae, and this read
    # `defaults["family"]` unconditionally. The register's surveyed: block
    # already carries nine flaviviruses, each declaring family: Flaviviridae,
    # and its stated promotion path is "an entry graduates by being moved into
    # viruses: above". Move one up and it would have been emitted as a
    # CORONAVIRUS -- the one field that matters for cross-family reasoning,
    # wrong, on the virus the whole v2 expansion is for.
    #
    # It would have validated cleanly: family is required on OrganismTaxon and
    # would have been present, just false. Nothing downstream re-derives it.
    #
    # The other three defaults happen to be right for Flaviviridae too
    # (+ssRNA, enveloped), which is exactly why this would not have announced
    # itself: one wrong field among three correct ones.
    vfamily = cfg.get("family") or defaults["family"]
    em.node(taxon_curie, "OrganismTaxon",
            label=name,
            family=vfamily,
            baltimore_class=cfg.get("baltimore_class") or defaults["baltimore_class"],
            genome_type=cfg.get("genome_type") or defaults.get("genome_type"),
            is_enveloped=defaults["is_enveloped"] if cfg.get("is_enveloped") is None
                         else cfg["is_enveloped"])

    pfile = proteome_dir / f"{cfg['proteome_id']}.json.gz"
    if not pfile.exists():
        em.notes[f"proteome_missing:{name}"] += 1
        return stats

    entries = load_proteome(pfile)

    # Swiss-Prot sometimes curates a reference ISOLATE rather than the species
    # node, leaving the species proteome TrEMBL-only and therefore without any
    # Chain annotation. Where config pins explicit accessions, they are loaded
    # alongside and their taxon_id is rewritten to the species, so cross-viral
    # metapaths resolve against one taxon per virus.
    extra = proteome_dir / f"EXTRA_{tax}.json.gz"
    if cfg.get("extra_accessions") and extra.exists():
        curated = load_proteome(extra)
        keep = set(cfg["extra_accessions"])
        curated = [e for e in curated if e.get("primaryAccession") in keep]
        # Curated entries REPLACE the proteome rather than augmenting it.
        # Accession-level merging leaves both copies, since the TrEMBL and
        # Swiss-Prot entries for the same protein have different accessions.
        entries = curated
        em.notes[f"curated_accessions_substituted:{name}"] += len(curated)
    seen_chains: dict[str, str] = {}          # normalised label -> node id
    parent_len: dict[str, int] = {}

    # Longest polyprotein first, so pp1ab claims the shared nsp1-nsp11 chains
    # before pp1a can. Ties broken on accession for determinism.
    def _order(e: dict) -> tuple:
        return (-len([f for f in e.get("features", []) if f.get("type") == "Chain"]),
                e.get("primaryAccession", ""))

    for entry in sorted(entries, key=_order):
        acc = entry["primaryAccession"]
        pcurie = f"UniProtKB:{acc}"
        seq = (entry.get("sequence") or {}).get("value", "")
        date = _entry_date(entry)
        symbol = _gene_symbol(entry)

        em.node(pcurie, "Protein",
                gene_symbol=symbol,
                taxon_id=taxon_curie,
                is_viral=True,
                sequence_hash=_sha(seq),
                reviewed=entry.get("entryType", "").startswith("UniProtKB reviewed"),
                protein_family=_rec_name(entry))
        stats["proteins"] += 1

        if symbol:
            gcurie = f"KGAV:GENE_{tax}_{symbol}"
            em.node(gcurie, "Gene", symbol=symbol, taxon_id=taxon_curie, is_viral=True)
            em.edge(pcurie, "ENCODED_BY", gcurie,
                    source=prov["proteome_source"], date=date, tier=prov["evidence_tier"])
            em.edge(gcurie, "BELONGS_TO", taxon_curie,
                    source=prov["proteome_source"], date=date, tier=prov["evidence_tier"])
            stats["genes"] += 1

        chains = [f for f in entry.get("features", []) if f.get("type") == "Chain"]
        parent_len[acc] = len(chains)

        # A single-chain entry IS the mature protein, so expanding it would
        # duplicate the same molecule. But sources cite the chain id anyway --
        # VirHostNet references ORF7a and ORF9b by PRO_ id, and those edges
        # silently failed to join while this branch just skipped. Emit an alias
        # to the parent instead of dropping the identifier.
        if len(chains) < 2:
            for feat in chains:
                cid = _chain_id(feat)
                if not cid:
                    continue
                em.node(f"UniProtKB:{cid}", "Protein",
                        taxon_id=taxon_curie,
                        is_viral=True,
                        sequence_hash=_sha(f"{acc}:{cid}"),
                        reviewed=entry.get("entryType", "").startswith("UniProtKB reviewed"))
                em.edge(f"UniProtKB:{cid}", "SAME_AS", pcurie,
                        source=prov["proteome_source"], date=date,
                        tier=prov["evidence_tier"],
                        quals={"merge_rule": "single_chain_equals_parent_protein"})
                em.notes["single_chain_aliased_to_parent"] += 1
            continue

        for feat in chains:
            label = _chain_label(feat)
            cid = _chain_id(feat)
            if not cid:
                em.notes[f"chain_without_pro_id:{name}"] += 1
                continue

            key = label.lower() or cid
            if key in seen_chains:
                # Same physical protein, different chain id: pp1a and pp1ab both
                # annotate nsp1-nsp11. The winner is kept as the node, but the
                # loser's id must remain resolvable -- VirHostNet and IntAct
                # cite pp1a chain ids freely, and without this alias those
                # interactions silently fail to join.
                # Minimal stub so the alias edge resolves. Deliberately NOT
                # given mature_peptide or protein_family: it is an identifier,
                # not a second copy of the protein, and must never be counted
                # as a distinct target or appear in a metapath.
                em.node(f"UniProtKB:{cid}", "Protein",
                        taxon_id=taxon_curie,
                        is_viral=True,
                        sequence_hash=_sha(f"{acc}:{cid}"),
                        reviewed=entry.get("entryType", "").startswith("UniProtKB reviewed"))
                em.edge(f"UniProtKB:{cid}", "SAME_AS", seen_chains[key],
                        source=prov["proteome_source"], date=date,
                        tier=prov["evidence_tier"],
                        quals={"merge_rule": "shared_chain_between_polyproteins"})
                em.notes["chain_shared_between_polyproteins"] += 1
                continue

            ccurie = f"UniProtKB:{cid}"
            seen_chains[key] = ccurie
            fam = chain_family(label, vfamily)

            em.node(ccurie, "Protein",
                    gene_symbol=symbol,
                    taxon_id=taxon_curie,
                    is_viral=True,
                    sequence_hash=_sha(f"{acc}:{cid}"),
                    reviewed=entry.get("entryType", "").startswith("UniProtKB reviewed"),
                    mature_peptide=label or cid,
                    protein_family=fam)

            if symbol:
                em.edge(ccurie, "ENCODED_BY", f"KGAV:GENE_{tax}_{symbol}",
                        source=prov["proteome_source"], date=date, tier=prov["evidence_tier"],
                        quals={"mature_peptide": label or cid,
                               "cleavage_order": _chain_start(feat)})
            stats["mature_peptides"] += 1

    mondo = cfg.get("mondo")
    if mondo:
        em.node(mondo, "Disease", label=f"{name} disease", is_infectious=True,
                causative_taxon=taxon_curie)
        em.edge(taxon_curie, "CAUSES", mondo,
                source=prov["disease_source"], date="1970-01-01",
                tier=prov["evidence_tier"], quals={"is_primary_agent": True})
        stats["diseases"] += 1
    else:
        em.notes[f"no_mondo_curated:{name}"] += 1

    return stats


def build(config_path: Path, proteome_dir: Path, out_dir: Path) -> Emit:
    cfg = yaml.safe_load(config_path.read_text())
    em = Emit()

    for v in cfg["viruses"]:
        s = ingest_virus(em, v, cfg["defaults"], cfg["provenance"], proteome_dir)
        print(f"  {v['name']:<12} proteins={s['proteins']:<4} "
              f"mature={s['mature_peptides']:<4} genes={s['genes']:<4} "
              f"disease={s['diseases']}")

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "nodes.jsonl").write_text(
        "\n".join(json.dumps(n) for n in em.nodes.values()))
    (out_dir / "edges.jsonl").write_text(
        "\n".join(json.dumps(e) for e in em.edges))
    return em


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser()
    root = Path(__file__).resolve().parents[1]
    ap.add_argument("--config", type=Path, default=root / "config" / "viruses.yaml")
    ap.add_argument("--proteomes", type=Path, default=root / "data" / "raw" / "proteomes")
    ap.add_argument("--out", type=Path, default=root / "data" / "releases" / "v0.1-spine")
    args = ap.parse_args()

    print("building spine")
    em = build(args.config, args.proteomes, args.out)

    print(f"\n{len(em.nodes):,} nodes, {len(em.edges):,} edges -> {args.out}")
    by_class = Counter(n["class"] for n in em.nodes.values())
    for c, n in by_class.most_common():
        print(f"  {c:<16} {n:>5}")
    if em.notes:
        print("\nnotes:")
        for k, v in em.notes.most_common():
            print(f"  {k:<44} {v:>5}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
