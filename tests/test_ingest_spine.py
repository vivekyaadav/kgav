"""Day 3 gate: every viral protein resolves to a virus, no orphans, no phantom
duplicate mature peptides from the pp1a/pp1ab overlap.
"""
import gzip
import json
import sys
from pathlib import Path

import pytest
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from ingest_spine import build, ingest_virus

from kgav.emit import Emit
from kgav.schema import load_schema

ROOT = Path(__file__).resolve().parents[1]
PROTEOMES = ROOT / "data" / "raw" / "proteomes"
CONFIG = ROOT / "config" / "viruses.yaml"
needs_data = pytest.mark.skipif(not PROTEOMES.exists(), reason="proteomes not downloaded")


def _chain(pid, desc, start):
    return {"type": "Chain", "featureId": pid, "description": desc,
            "location": {"start": {"value": start}, "end": {"value": start + 100}}}


@pytest.fixture
def fake_proteome(tmp_path):
    """pp1ab (16 chains) + pp1a (12 overlapping) + spike (1 chain)."""
    pp1ab = [_chain(f"PRO_00004496{10+i}", f"Non-structural protein {i+1} (nsp{i+1})", 1 + i * 180)
             for i in range(16)]
    pp1a = [_chain(f"PRO_00004497{10+i}", f"Non-structural protein {i+1} (nsp{i+1})", 1 + i * 180)
            for i in range(12)]
    entries = [
        {"primaryAccession": "P0DTD1", "entryType": "UniProtKB reviewed (Swiss-Prot)",
         "proteinDescription": {"recommendedName": {"fullName": {"value": "Replicase polyprotein 1ab"}}},
         "genes": [{"geneName": {"value": "rep"}}], "sequence": {"value": "MESLVPGFNEK"},
         "entryAudit": {"firstPublicDate": "2020-04-22"}, "features": pp1ab},
        {"primaryAccession": "P0DTC1", "entryType": "UniProtKB reviewed (Swiss-Prot)",
         "proteinDescription": {"recommendedName": {"fullName": {"value": "Replicase polyprotein 1a"}}},
         "genes": [{"geneName": {"value": "rep"}}], "sequence": {"value": "MESLVPGFNE"},
         "entryAudit": {"firstPublicDate": "2020-04-22"}, "features": pp1a},
        {"primaryAccession": "P0DTC2", "entryType": "UniProtKB reviewed (Swiss-Prot)",
         "proteinDescription": {"recommendedName": {"fullName": {"value": "Spike glycoprotein"}}},
         "genes": [{"geneName": {"value": "S"}}], "sequence": {"value": "MFVFLVLLPLV"},
         "entryAudit": {"firstPublicDate": "2020-04-22"},
         "features": [_chain("PRO_0000449647", "Spike glycoprotein", 1)]},
    ]
    d = tmp_path / "proteomes"
    d.mkdir()
    (d / "UP000464024.json.gz").write_bytes(gzip.compress(json.dumps({"results": entries}).encode()))
    return d


@pytest.fixture
def one_virus_config(tmp_path):
    cfg = yaml.safe_load(CONFIG.read_text())
    cfg["viruses"] = [v for v in cfg["viruses"] if v["taxon"] == "2697049"]
    p = tmp_path / "viruses.yaml"
    p.write_text(yaml.safe_dump(cfg))
    return p


@pytest.fixture
def spine(fake_proteome, one_virus_config, tmp_path):
    return build(one_virus_config, fake_proteome, tmp_path / "out")


def test_shared_chains_are_not_duplicated(spine):
    """pp1a and pp1ab share nsp1-nsp11. Sixteen mature peptides, not 28."""
    mature = [n for n in spine.nodes.values()
              if n["class"] == "Protein" and n["properties"].get("mature_peptide")]
    assert len(mature) == 16, f"got {len(mature)} mature peptides"
    assert spine.notes["chain_shared_between_polyproteins"] == 12


def test_single_chain_entries_become_aliases_not_targets(spine):
    """Spike is one chain: the chain IS the protein, so it must not become a
    second target. But sources cite the chain id (VirHostNet references ORF7a
    and ORF9b that way), so the identifier has to resolve. It exists as an
    alias -- no mature_peptide, no protein_family, SAME_AS to the parent."""
    cid = "UniProtKB:PRO_0000449647"
    assert cid in spine.nodes, "chain id must resolve, or source edges silently drop"
    props = spine.nodes[cid]["properties"]
    assert not props.get("mature_peptide"), "alias must not count as a target"
    assert not props.get("protein_family")
    aliases = {(e["subject"], e["object"]) for e in spine.edges if e["predicate"] == "SAME_AS"}
    assert (cid, "UniProtKB:P0DTC2") in aliases


def test_nsp_family_is_extracted(spine):
    fams = {n["properties"].get("protein_family") for n in spine.nodes.values()}
    assert "nsp5" in fams, "Mpro (nsp5) not labelled — direct-acting targets unreachable"
    assert "nsp12" in fams, "RdRp (nsp12) not labelled"


def test_every_viral_protein_reaches_its_virus(spine):
    """No orphans: each protein has a path to the taxon via its gene."""
    genes = {e["object"] for e in spine.edges if e["predicate"] == "ENCODED_BY"}
    to_taxon = {e["subject"] for e in spine.edges if e["predicate"] == "BELONGS_TO"}
    assert genes and genes <= to_taxon, "gene(s) with no BELONGS_TO edge"


def test_all_proteins_flagged_viral(spine):
    for n in spine.nodes.values():
        if n["class"] == "Protein":
            assert n["properties"]["is_viral"] is True


def test_output_validates_against_schema(spine):
    schema = load_schema()
    violations = schema.validate_batch(spine.nodes.values(), spine.edges)
    assert violations == [], [str(v) for v in violations[:10]]


def test_no_dangling_edges(spine):
    ids = set(spine.nodes)
    for e in spine.edges:
        assert e["subject"] in ids and e["object"] in ids, e


def test_every_edge_has_provenance_and_date(spine):
    for e in spine.edges:
        assert e["primary_knowledge_source"]
        assert e["first_asserted_date"]
        assert e["evidence_tier"] in (1, 2, 3, 4)


def test_build_is_deterministic(fake_proteome, one_virus_config, tmp_path):
    a = build(one_virus_config, fake_proteome, tmp_path / "a")
    b = build(one_virus_config, fake_proteome, tmp_path / "b")
    assert list(a.nodes) == list(b.nodes)
    assert a.edges == b.edges


# ------------------------------------------------------------- real data
@needs_data
def test_real_spine_validates(tmp_path):
    spine = build(CONFIG, PROTEOMES, tmp_path / "real")
    schema = load_schema()
    violations = schema.validate_batch(spine.nodes.values(), spine.edges)
    assert violations == [], [str(v) for v in violations[:10]]


@needs_data
def test_real_spine_has_all_seven_viruses(tmp_path):
    spine = build(CONFIG, PROTEOMES, tmp_path / "real")
    taxa = [n for n in spine.nodes.values() if n["class"] == "OrganismTaxon"]
    # DERIVED FROM THE REGISTER, not hardcoded. This read `== 7` and broke
    # when DENV-2 and Zika were promoted -- correctly, but dryrun_promote's
    # breakage scan looked only for `len(doc["viruses"]) == N` and missed both
    # assertions in this file, so the promotion landed reporting one breakage
    # when there were three.
    expected = len(yaml.safe_load(CONFIG.read_text())["viruses"])
    assert len(taxa) == expected, (
        f"register declares {expected} built viruses, spine built "
        f"{len(taxa)}: {sorted(n['properties']['label'] for n in taxa)}")


@needs_data
def test_real_spine_exposes_mpro_and_rdrp(tmp_path):
    """The two targets every coronavirus antiviral programme aims at."""
    spine = build(CONFIG, PROTEOMES, tmp_path / "real")
    fams: dict[str, int] = {}
    for n in spine.nodes.values():
        f = n["properties"].get("protein_family")
        if f:
            fams[f] = fams.get(f, 0) + 1
    assert fams.get("nsp5", 0) >= 1, "no Mpro found in any proteome"
    assert fams.get("nsp12", 0) >= 1, "no RdRp found in any proteome"


@needs_data
def test_every_coronavirus_exposes_the_core_nsps(tmp_path):
    """Every CORONAVIRUS must label the same nsps.

    SCOPED BY FAMILY. This asserted it of every built virus and failed the
    moment DENV-2 was promoted, with "DENV-2 missing ['nsp12', 'nsp13', ...]".
    That failure was correct and the assertion was wrong: flaviviruses have
    no nsps. They carry NS1-NS5 -- the protease is NS3 with an NS2B cofactor,
    the polymerase is NS5 -- and nsp5 is a coronavirus main protease. An
    invariant about coronavirus naming cannot be asserted of a family that
    does not use the naming.

    Alphacoronaviruses name mature peptides by function ("3C-like proteinase")
    where betacoronaviruses use the number. If the alias table regresses, this
    produces no error anywhere -- just a graph where cross-viral metapaths
    cannot transfer, and a leave-one-virus-out evaluation that underperforms
    for reasons nothing reports.
    """
    spine = build(CONFIG, PROTEOMES, tmp_path / "nsp")
    core = {"nsp3", "nsp5", "nsp12", "nsp13", "nsp14", "nsp15", "nsp16"}
    by_taxon: dict[str, set[str]] = {}
    corona: dict[str, str] = {}
    labels: dict[str, str] = {}
    for n in spine.nodes.values():
        if n["class"] == "OrganismTaxon":
            labels[n["id"]] = n["properties"]["label"]
            if n["properties"].get("family") == "Coronaviridae":
                corona[n["id"]] = n["properties"]["label"]
        elif n["class"] == "Protein" and n["properties"].get("protein_family"):
            by_taxon.setdefault(n["properties"]["taxon_id"], set()).add(
                n["properties"]["protein_family"])

    assert len(labels) == len(yaml.safe_load(CONFIG.read_text())["viruses"])
    assert corona, "no Coronaviridae taxon carries family metadata"
    for tid, name in corona.items():
        missing = core - by_taxon.get(tid, set())
        assert not missing, f"{name} missing {sorted(missing)}"


@needs_data
def test_every_built_virus_builds_at_least_one_protein(tmp_path):
    """A virus declared built must produce proteins, or it is an orphan taxon.

    This is what two orphan OrganismTaxon nodes in v0.1 looked like from the
    test suite: the register promised a virus, the spine delivered a label,
    and nothing in between reported it except a note in a build log. It was
    xfail(strict=True) for exactly as long as DENV-2 and Zika sat in
    `viruses:` without cached proteomes -- roughly an hour -- and the xfail is
    gone because the register is, not because the gap was papered over.

    It now guards forward: promote a virus whose proteome is not in
    data/raw/proteomes and this fails BEFORE ingest_chembl emits 1,881
    DANGLING edges and assemble_release invalidates every results file.
    dryrun_promote.py Section 0 checks the same precondition on demand; this
    is the copy that runs whether or not anyone remembers to ask.
    """
    spine = build(CONFIG, PROTEOMES, tmp_path / "proteins")
    counts: dict[str, int] = {}
    labels: dict[str, str] = {}
    for n in spine.nodes.values():
        if n["class"] == "OrganismTaxon":
            labels[n["id"]] = n["properties"]["label"]
            counts.setdefault(n["id"], 0)
        elif n["class"] == "Protein":
            tid = n["properties"].get("taxon_id")
            if tid:
                counts[tid] = counts.get(tid, 0) + 1

    empty = sorted(labels[t] for t, c in counts.items() if c == 0)
    assert not empty, f"declared built but zero proteins: {empty}"


@needs_data
def test_every_built_virus_can_be_reached_from_its_proteins(tmp_path):
    """The M1 TAIL, which the proteins test above does not cover.

    M1 is SmallMolecule -INHIBITS-> Protein -ENCODED_BY-> Gene -BELONGS_TO->
    OrganismTaxon. A virus with proteins but no Gene node has no ENCODED_BY
    and no BELONGS_TO, so nothing can walk from a compound to the taxon and
    the virus is absent from the direct-acting channel entirely.

    THIS IS NOT HYPOTHETICAL. The 2026-10-12 Zika build printed

        Zika   proteins=2   mature=13   genes=0   disease=0

    with BELONGS_TO still at 70 and ENCODED_BY unchanged at 19,370. Zika had
    thirteen mature peptides, 3,653 assay rows resolving to NS3, schema PASS
    -- and was unreachable on the one channel the promotion was for. Q32ZE1
    has no `genes` field, ingest_virus gates the whole Gene block on
    `if symbol`, and nothing reported it as a problem.

    The proteins test passed, because 2 > 0. Proteins are necessary and not
    sufficient; this asserts the rest of the path.
    """
    spine = build(CONFIG, PROTEOMES, tmp_path / "proteins")

    labels = {n["id"]: n["properties"]["label"]
              for n in spine.nodes.values() if n["class"] == "OrganismTaxon"}
    genes: dict[str, int] = {t: 0 for t in labels}
    for n in spine.nodes.values():
        if n["class"] == "Gene":
            tid = n["properties"].get("taxon_id")
            if tid in genes:
                genes[tid] += 1
    reached = {e["object"] for e in spine.edges
               if e["predicate"] == "BELONGS_TO"}

    no_gene = sorted(labels[t] for t, c in genes.items() if c == 0)
    assert not no_gene, (
        f"built with proteins but no Gene node, so unreachable by M1: "
        f"{no_gene}")

    unreached = sorted(labels[t] for t in labels if t not in reached)
    assert not unreached, (
        f"no BELONGS_TO edge, so no metapath ending at the taxon can reach "
        f"them: {unreached}")


# --------------------------------------------------------------------------
# The register's promotion path: "an entry graduates by being moved into
# viruses: above". It must not graduate into the wrong family.
# --------------------------------------------------------------------------
def test_a_graduated_flavivirus_keeps_its_own_family(tmp_path):
    """defaults.family is Coronaviridae, and this read it unconditionally.

    config/viruses.yaml already surveys nine flaviviruses, each declaring
    family: Flaviviridae, and the documented way to build one is to move its
    entry into `viruses:`. Doing that would have emitted a Dengue taxon as
    family=Coronaviridae -- the single field cross-family reasoning depends on,
    wrong, on the viruses the v2 expansion exists for. It validates cleanly:
    family is required and would have been present, just false.
    """
    from ingest_spine import Emit, ingest_virus

    defaults = {"family": "Coronaviridae", "baltimore_class": "IV",
                "genome_type": "ssRNA(+)", "is_enveloped": True}
    prov = {"proteome_source": "infores:uniprot", "disease_source": "infores:mondo",
            "evidence_tier": 1}

    em = Emit()
    ingest_virus(em, {"taxon": "11060", "name": "DENV-2", "family": "Flaviviridae",
                      "proteome_id": "UP000095876"},
                 defaults, prov, tmp_path)          # no proteome file: taxon only
    props = em.nodes["NCBITaxon:11060"]["properties"]
    assert props["family"] == "Flaviviridae"
    # The defaults that are genuinely shared still apply.
    assert props["baltimore_class"] == "IV"
    assert props["is_enveloped"] is True


def test_a_coronavirus_without_its_own_family_still_takes_the_default(tmp_path):
    """The override must not become 'every entry must restate everything'."""
    from ingest_spine import Emit, ingest_virus

    defaults = {"family": "Coronaviridae", "baltimore_class": "IV",
                "genome_type": "ssRNA(+)", "is_enveloped": True}
    prov = {"proteome_source": "infores:uniprot", "disease_source": "infores:mondo",
            "evidence_tier": 1}
    em = Emit()
    ingest_virus(em, {"taxon": "2697049", "name": "SARS-CoV-2",
                      "proteome_id": "UP000464024"}, defaults, prov, tmp_path)
    assert em.nodes["NCBITaxon:2697049"]["properties"]["family"] == "Coronaviridae"


# ------------- _nsp_number is coronavirus-only and partly matches flaviviruses
def test_nsp_number_mislabels_flavivirus_chains_and_misses_the_protease():
    """What the spine would assign to DENV-2's chains, measured against
    P29990's real 13 CHAIN labels (UniProt, 2026-10-11):

        nsp1    Non-structural protein 1                      <- it is NS1
        nsp12   RNA-directed RNA polymerase/Methyltransferase NS5
        None    Serine protease NS3                           <- THE PROBLEM
        None    Serine protease subunit NS2B
        None    Non-structural protein 2A / 4A / 4B
        None    Capsid protein C, prM, pr, M, Envelope protein E

    THE TWO HALVES DO NOT FAIL THE SAME WAY, and that is the finding:

    NS5 WORKS BY ACCIDENT. _nsp_number reads the label as nsp12 and
    chembl.resolve_nsp reads "NS5 RdRp" as nsp12 too, so
    viral_targets["11060"]["nsp12"] is DENV-2's NS5 chain and the edge lands
    on the right protein through a coronavirus family name. The index is
    keyed per taxon, so there is no cross-family contamination -- only a
    wrong name.

    NS3 IS UNREACHABLE. "Serine protease NS3" matches none of _nsp_number's
    three routes, so the chain gets protein_family=None and never enters
    viral_targets. A "NS2B-NS3 protease" description resolves to None and a
    "NS3 helicase" one to nsp13, and neither has an entry -- so both fall
    through to UniProtKB:P29990, the undifferentiated polyprotein. Every
    flavivirus protease measurement would attach to the whole genome product,
    which is the chloroquine-on-spike failure in a new place: one node
    standing in for thirteen proteins.

    PINNED, NOT FIXED. The fix is family-scoped chain naming on both sides at
    once -- NS1..NS5 from the labels AND from the descriptions -- and a
    half-fix is worse than none, because changing one side alone breaks the
    accident that currently makes NS5 work. The ChEMBL description
    distribution that should drive it is measured by
    diagnose_flavivirus_targets.py.
    """
    from ingest_spine import _nsp_number

    # the accidental match that currently works
    assert _nsp_number("RNA-directed RNA polymerase/Methyltransferase NS5") == 12

    # the mislabelled one
    assert _nsp_number("Non-structural protein 1") == 1

    # and the protease, which is the one that matters and is absent
    for label in ("Serine protease NS3", "Serine protease subunit NS2B",
                  "Non-structural protein 2A", "Non-structural protein 4A",
                  "Non-structural protein 4B", "Envelope protein E",
                  "Capsid protein C", "Genome polyprotein"):
        assert _nsp_number(label) is None, label

    # coronavirus naming, which must not move
    assert _nsp_number("3C-like proteinase nsp5") == 5
    assert _nsp_number("Non-structural protein 3") == 3


# Verbatim from rest.uniprot.org, read 2026-10-11 for BOTH polyproteins.
FLAVIVIRUS_CHAINS = (
    "Genome polyprotein", "Capsid protein C", "Protein prM", "Peptide pr",
    "Small envelope protein M", "Envelope protein E",
    "Non-structural protein 1", "Non-structural protein 2A",
    "Serine protease subunit NS2B", "Serine protease NS3",
    "Non-structural protein 4A", "Non-structural protein 4B",
    "RNA-directed RNA polymerase/Methyltransferase NS5",
)


def test_zika_and_denv2_chains_share_one_vocabulary():
    """Read on 2026-10-11: Q32ZE1 (Zika, 64320, a component of UP000054557)
    and P29990 (DENV-2 strain Thailand/16681/1984, 31634, no proteome
    cross-reference at all) each carry thirteen Chain features with the SAME
    thirteen descriptions, character for character.

    The test above measured P29990 alone, which left open whether each virus
    would need its own curated patterns -- nine entries in `surveyed:` would
    have meant nine curation jobs. It does not: UniProt annotates the family
    from one controlled vocabulary, so ONE chain-naming rule covers Zika,
    DENV-2 and the seven flaviviruses still at zero protein rows.

    This is the cheap half of the gate. Pinned here so that if a later commit
    writes per-virus patterns, this says they are not needed.
    """
    from ingest_spine import _nsp_number

    assert len(FLAVIVIRUS_CHAINS) == 13 == len(set(FLAVIVIRUS_CHAINS))

    keyed = {c: _nsp_number(c) for c in FLAVIVIRUS_CHAINS}
    assert sum(1 for v in keyed.values() if v is None) == 11, keyed
    assert [c for c, v in keyed.items() if v is not None] == [
        "Non-structural protein 1",                           # wrong: NS1
        "RNA-directed RNA polymerase/Methyltransferase NS5",  # right by accident
    ]


def test_the_chembl_join_can_reach_the_flavivirus_protease():
    """DECIDED, and the marker deleted as its own reason instructed.

    This was xfail(strict=True) while the key scheme was an open question.
    The question is settled by src/kgav/folds.py:41 -- "NUCLEOPHILE mode:
    Identical residue required. Cys != Ser, no edge, no exceptions." Flavivirus
    NS2B-NS3 is a serine protease and coronavirus nsp5 a cysteine protease, so
    a shared functional key such as "protease" would assert inside
    protein_family the very cross-family equivalence that gate refuses, and
    assert it upstream of the gate. Each family therefore keeps its own
    namespace, and folds/TargetClass stays the one place a cross-family claim
    is made and checked.

    What it was worth fixing for: 315 of the 355 Zika rows that reach the
    chain lookup name a protease, against 9 for the RdRp
    (scripts/diagnose_flavivirus_targets.py, 2026-10-11). Those 315 resolved
    to nothing and landed on the undifferentiated polyprotein.

    _nsp_number itself is deliberately unchanged -- it remains the
    coronavirus-only function, which is why the test above still measures it
    returning None here. chain_family is the family-aware entry point, and
    tests/test_flavivirus_naming.py holds the two-sided detail.
    """
    from ingest_spine import chain_family

    key = chain_family("Serine protease NS3", "Flaviviridae")
    assert key == "NS3", "the protease must be keyable"
    assert not key.lower().startswith("nsp"), \
        f"keyed as {key}, a coronavirus protein's name"


def test_the_declared_gene_symbol_is_looked_up_per_accession():
    """The MECHANISM, tested on a synthetic map rather than on Zika's.

    The first fix for Zika's genes=0 declared ONE symbol for the virus and
    handed it to every entry lacking a gene name, so the register could not
    say different things about different accessions even where it should.
    Per-accession lookup is what makes the claim explicit and auditable.

    Zika's own two accessions now BOTH map to POLY, on UniProt's evidence --
    see test_the_register_declares_one_gene_for_zikas_two_accessions. That is
    a fact about Zika, not about the lookup, so it is tested there and the
    mechanism is tested here with symbols that do differ.
    """
    from ingest_spine import _gene_symbol

    declared = {"ACC_A": "GENE_A", "ACC_B": "GENE_B"}
    assert _gene_symbol({"primaryAccession": "ACC_A"}, declared) == "GENE_A"
    assert _gene_symbol({"primaryAccession": "ACC_B"}, declared) == "GENE_B"
    assert _gene_symbol({"primaryAccession": "ACC_A"}, declared) \
        != _gene_symbol({"primaryAccession": "ACC_B"}, declared), \
        "the register must be able to distinguish two accessions"

    # An accession the register does not name gets nothing. A protein with no
    # gene is unreachable by M1, which is honest; a protein on the wrong gene
    # is a false path.
    assert _gene_symbol({"primaryAccession": "P0DTG0"}, declared) is None
    # And a real UniProt gene name always wins over any declaration.
    assert _gene_symbol(
        {"primaryAccession": "ACC_A",
         "genes": [{"geneName": {"value": "rep"}}]}, declared) == "rep"


def test_the_register_declares_one_gene_for_zikas_two_accessions():
    """ONE GENE, TWO PRODUCTS, on UniProt's own evidence.

    This test used to assert the opposite -- that the two accessions must
    NOT share a symbol -- on the reasoning that "a uORF is by definition a
    different reading frame from the polyprotein". That reasoning was wrong
    twice over: a uORF is defined positionally and the definition says
    nothing about frame, and UniProt states the gene relationship outright.

    Read from data/raw/proteomes/UP000054557.json.gz, every annotation
    evidenced ECO:0000269 against PubMed:39394194:

      MISCELLANEOUS        "Product of an upstream open reading frame of the
                            genome polyprotein GENE"
      ALTERNATIVE PRODUCTS  events ["Alternative initiation"], isoforms
                            P0DXN9-1 "uORF protein" and Q32ZE1-1 "Genome
                            polyprotein"

    Two independent statements that this is one gene. So they share a Gene
    node, and the node_merged:Gene an earlier commit read as a defect is the
    correct outcome.
    """
    import yaml

    cfg = yaml.safe_load(CONFIG.read_text())
    zika = next(v for v in cfg["viruses"] if v["name"] == "Zika")
    declared = zika["gene_symbol"]
    assert isinstance(declared, dict), "per accession, not one string"
    assert declared == {"Q32ZE1": "POLY", "P0DXN9": "POLY"}
    assert len(set(declared.values())) == 1, \
        "UniProt lists these as alternative-initiation isoforms of one gene"


def test_no_script_defines_its_own_emit():
    """ONE ACCUMULATOR. ingest_spine.py carried its own `class Emit`,
    identical to src/kgav/emit.py except that the shared one had since
    gained a `score` parameter -- and the comment explaining why that field
    is not called confidence. The fork missed it in silence.

    Emit.node merges with `properties.update(props)`: last write wins, with
    only a node_merged counter. That semantics decides what a repeated node
    assignment means, and a second copy of it is a second place the rule can
    diverge. It is not hypothetical -- the same merge is how a TargetClass
    holding a cysteine protease and a serine protease came to take whichever
    member was emitted last.

    This is the fifth instance in this repository of one thing living in two
    files: the percentile arithmetic, `trustworthy`, hits_range, the
    metapath constraints, and now the accumulator itself.
    """
    import ast
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    offenders = []
    for p in list((root / "scripts").glob("*.py")) + list(
            (root / "src/kgav").rglob("*.py")):
        if p == root / "src/kgav/emit.py":
            continue
        for n in ast.walk(ast.parse(p.read_text())):
            if isinstance(n, ast.ClassDef) and n.name == "Emit":
                offenders.append(str(p.relative_to(root)))
    assert not offenders, f"Emit redefined outside src/kgav/emit.py: {offenders}"


# ------------------------------- the proteome guard checked the wrong thing
#
# DENV-2's register entry says in as many words that "UP000095876.json.gz must
# still exist, because the pfile.exists() guard returns early before the EXTRA
# file is consulted" -- and then the substitution line throws its contents
# away. On 2026-10-23 UniProt returned an empty result set for UP000095876,
# while an identically formed query for West Nile's UP000008600 returned its
# entry the same minute, so DENV-2 became unbuildable: blocked on a file
# nothing reads.

def _entry(acc, name="Genome polyprotein", chains=2, gene="POLY"):
    return {"primaryAccession": acc,
            "entryType": "UniProtKB reviewed (Swiss-Prot)",
            "proteinDescription": {"recommendedName": {"fullName":
                                                       {"value": name}}},
            "genes": [{"geneName": {"value": gene}}],
            "sequence": {"value": "MNNQRKKARN"},
            "entryAudit": {"firstPublicDate": "2008-01-01"},
            "features": [_chain(f"PRO_900000{i}", f"Protein NS{i+1}", 1 + i * 60)
                         for i in range(chains)]}


def _write(d: Path, stem: str, entries: list) -> None:
    (d / f"{stem}.json.gz").write_bytes(
        gzip.compress(json.dumps({"results": entries}).encode()))


@pytest.fixture
def curated_cfg():
    return {"taxon": "11060", "name": "DENV-2", "family": "Flaviviridae",
            "full_name": "Dengue virus type 2",
            "proteome_id": "UP000095876",
            "extra_accessions": ["P29990"],
            "isolate_taxa": ["31634"]}


# The real blocks, read from config/viruses.yaml rather than invented, so a
# field added there fails here instead of being quietly absent.
_CFG = yaml.safe_load(CONFIG.read_text())
DEFAULTS = _CFG["defaults"]
PROV = _CFG["provenance"]


def _run(cfg, d):
    em = Emit()
    ingest_virus(em, cfg, DEFAULTS, PROV, d)
    return em


def test_curated_accessions_build_the_virus_without_the_proteome(
        curated_cfg, tmp_path):
    """THE FIX. The proteome file is absent and unobtainable; the curated
    accessions are present and are what the build actually consumes."""
    d = tmp_path / "p"
    d.mkdir()
    _write(d, "EXTRA_11060", [_entry("P29990")])
    em = _run(curated_cfg, d)
    assert em.notes["proteome_missing:DENV-2"] == 0
    assert em.notes["curated_accessions_substituted:DENV-2"] == 1
    assert em.notes["proteome_absent_curated_only:DENV-2"] == 1
    assert any(n["id"] == "UniProtKB:P29990" for n in em.nodes.values())


def test_the_absence_is_recorded_not_silent(curated_cfg, tmp_path):
    """Building from a substitute is a different provenance from building
    from the declared proteome, and the run must be able to say so."""
    d = tmp_path / "p"
    d.mkdir()
    _write(d, "EXTRA_11060", [_entry("P29990")])
    em = _run(curated_cfg, d)
    assert em.notes["proteome_absent_curated_only:DENV-2"] == 1


def test_a_present_proteome_does_not_raise_the_absent_note(curated_cfg,
                                                           tmp_path):
    d = tmp_path / "p"
    d.mkdir()
    _write(d, "UP000095876", [_entry("A0A000", chains=0)])
    _write(d, "EXTRA_11060", [_entry("P29990")])
    em = _run(curated_cfg, d)
    assert em.notes["curated_accessions_substituted:DENV-2"] == 1
    assert em.notes["proteome_absent_curated_only:DENV-2"] == 0


def test_curated_entries_still_replace_the_proteome(curated_cfg, tmp_path):
    """Unchanged behaviour: substitution, not augmentation. The TrEMBL and
    Swiss-Prot copies of one protein have different accessions, so merging
    would leave both."""
    d = tmp_path / "p"
    d.mkdir()
    _write(d, "UP000095876", [_entry("A0A000", name="Polyprotein (TrEMBL)")])
    _write(d, "EXTRA_11060", [_entry("P29990")])
    em = _run(curated_cfg, d)
    accs = {n["id"] for n in em.nodes.values() if n["id"].startswith("UniProtKB:")}
    assert "UniProtKB:P29990" in accs
    assert "UniProtKB:A0A000" not in accs


def test_neither_file_is_still_proteome_missing(curated_cfg, tmp_path):
    """The guard must not have been loosened into never firing."""
    d = tmp_path / "p"
    d.mkdir()
    em = _run(curated_cfg, d)
    assert em.notes["proteome_missing:DENV-2"] == 1
    assert not [n for n in em.nodes.values()
                if n["id"].startswith("UniProtKB:")]


def test_an_extra_file_holding_none_of_the_pinned_accessions(curated_cfg,
                                                             tmp_path):
    """A LATENT BUG THIS EXPOSED. With the proteome present and EXTRA holding
    no pinned accession, the old code set entries = [] and built the virus
    with zero proteins out of two files that were both there."""
    d = tmp_path / "p"
    d.mkdir()
    _write(d, "UP000095876", [_entry("A0A000")])
    _write(d, "EXTRA_11060", [_entry("Q99999")])      # not P29990
    em = _run(curated_cfg, d)
    assert em.notes["curated_accessions_absent:DENV-2"] == 1
    accs = {n["id"] for n in em.nodes.values() if n["id"].startswith("UniProtKB:")}
    assert "UniProtKB:A0A000" in accs, "falls back to the proteome, not to nothing"


def test_no_proteome_and_an_unmatched_extra_is_missing_not_empty(curated_cfg,
                                                                 tmp_path):
    d = tmp_path / "p"
    d.mkdir()
    _write(d, "EXTRA_11060", [_entry("Q99999")])
    em = _run(curated_cfg, d)
    assert em.notes["curated_accessions_absent:DENV-2"] == 1
    assert em.notes["proteome_missing:DENV-2"] == 1


def test_a_virus_with_no_extra_accessions_is_unaffected(tmp_path):
    """West Nile declares none and must behave exactly as before."""
    cfg = {"taxon": "11082", "name": "West Nile", "family": "Flaviviridae",
           "full_name": "West Nile virus", "proteome_id": "UP000008600"}
    d = tmp_path / "p"
    d.mkdir()
    em = _run(cfg, d)
    assert em.notes["proteome_missing:West Nile"] == 1
    _write(d, "UP000008600", [_entry("P06935")])
    em2 = _run(cfg, d)
    assert em2.notes["proteome_missing:West Nile"] == 0
    assert em2.notes["curated_accessions_substituted:West Nile"] == 0
