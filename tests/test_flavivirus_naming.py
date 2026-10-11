"""Two-sided flavivirus protein naming, against real strings on both sides.

WHY THE TWO SIDES ARE ONE COMMIT

A chain's protein_family (written by ingest_spine) and an assay description's
resolved family (read by chembl.resolve_nsp) are the two halves of one join
key: build_viral_target_index maps taxon -> {family: chain node id}, and
ingest_activities looks the description's answer up in it. Change one side
alone and the key stops matching -- which would have broken the single
flavivirus join that works today, NS5/RdRp resolving to nsp12 on both sides by
coincidence. So both move here, together.

WHY THE NAMESPACES STAY SEPARATE

src/kgav/folds.py:41 -- "NUCLEOPHILE mode: Identical residue required.
Cys != Ser, no edge, no exceptions." Coronavirus nsp5 is a cysteine protease;
the flavivirus NS2B-NS3 protease is a serine protease. A shared functional key
like "protease" would assert that equivalence inside protein_family, silently,
upstream of the gate written to refuse it. folds/TargetClass stays the one
place a cross-family claim is made and checked.

PROVENANCE OF EVERY STRING HERE

Chain descriptions: rest.uniprot.org, 2026-10-11.
  Q32ZE1  Zika virus (64320), component of UP000054557 -- 13 Chain features
  P29990  DENV-2 strain Thailand/16681/1984 (31634), no proteome xref -- 13
  P0DTD1  SARS-CoV-2 (2697049) pp1ab -- 16 Chain features
Assay descriptions: quoted from scripts/diagnose_flavivirus_targets.py output
on the real ChEMBL database, same date.

Nothing here is a string I composed to make a pattern pass.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

from kgav.chembl import FLAVI_PATTERNS, NSP_PATTERNS, resolve_nsp

ROOT = Path(__file__).resolve().parents[1]


def _spine():
    spec = importlib.util.spec_from_file_location(
        "ingest_spine_naming", ROOT / "scripts" / "ingest_spine.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m          # @dataclass resolves via sys.modules
    spec.loader.exec_module(m)
    return m


spine = _spine()
chain_family = spine.chain_family

# Q32ZE1 and P29990, verbatim and in genome order. Identical on both.
FLAVI_CHAINS = (
    ("Genome polyprotein", None),          # the parent product, not a peptide
    ("Capsid protein C", "C"),
    ("Protein prM", "prM"),
    ("Peptide pr", "pr"),
    ("Small envelope protein M", "M"),
    ("Envelope protein E", "E"),
    ("Non-structural protein 1", "NS1"),
    ("Non-structural protein 2A", "NS2A"),
    ("Serine protease subunit NS2B", "NS2B"),
    ("Serine protease NS3", "NS3"),
    ("Non-structural protein 4A", "NS4A"),
    ("Non-structural protein 4B", "NS4B"),
    ("RNA-directed RNA polymerase/Methyltransferase NS5", "NS5"),
)

# P0DTD1's sixteen, verbatim. Every one carries an explicit nspN token, which
# is why the coronavirus routes need no change: they never relied on the
# function-name table for this polyprotein.
CORONA_CHAINS = (
    ("Replicase polyprotein 1ab", None),
    ("Host translation inhibitor nsp1", "nsp1"),
    ("Non-structural protein 2", "nsp2"),
    ("Papain-like protease nsp3", "nsp3"),
    ("Non-structural protein 4", "nsp4"),
    ("3C-like proteinase nsp5", "nsp5"),
    ("Non-structural protein 6", "nsp6"),
    ("Non-structural protein 7", "nsp7"),
    ("Non-structural protein 8", "nsp8"),
    ("Viral protein genome-linked nsp9", "nsp9"),
    ("Non-structural protein 10", "nsp10"),
    ("RNA-directed RNA polymerase nsp12", "nsp12"),
    ("Helicase nsp13", "nsp13"),
    ("Guanine-N7 methyltransferase nsp14", "nsp14"),
    ("Uridylate-specific endoribonuclease nsp15", "nsp15"),
    ("2'-O-methyltransferase nsp16", "nsp16"),
)


# --------------------------------------------------------------------------
# the chain side
# --------------------------------------------------------------------------
@pytest.mark.parametrize("label,expected", FLAVI_CHAINS,
                         ids=[c[0] for c in FLAVI_CHAINS])
def test_every_flavivirus_chain_gets_its_own_name(label, expected):
    assert chain_family(label, "Flaviviridae") == expected


def test_all_twelve_mature_flavivirus_chains_are_now_keyed():
    """Eleven of thirteen were unkeyed before this commit, and chembl.py skips
    a chain with no protein_family -- so they could not be joined at all."""
    keyed = {c: chain_family(c, "Flaviviridae") for c, _ in FLAVI_CHAINS}
    assert sum(1 for v in keyed.values() if v is None) == 1, keyed
    assert keyed["Genome polyprotein"] is None, "the parent stays unkeyed"
    assert len(set(v for v in keyed.values() if v)) == 12, \
        "twelve distinct chains, no two collapsing onto one key"


def test_prm_and_pr_do_not_collapse_onto_each_other():
    """"Protein prM" and "Peptide pr" are different chains at different
    positions (123-290 and 123-215 on Q32ZE1). A prefix regex claims both,
    which is why these five are matched as exact strings."""
    assert chain_family("Protein prM", "Flaviviridae") == "prM"
    assert chain_family("Peptide pr", "Flaviviridae") == "pr"


def test_ns2a_is_not_read_as_ns2():
    """A bare (\\d+) capture reads "Non-structural protein 2A" as NS2 -- a
    different chain that also exists in the vocabulary. The [ab] group is the
    whole reason this is not a one-line regex."""
    assert chain_family("Non-structural protein 2A", "Flaviviridae") == "NS2A"
    assert chain_family("Non-structural protein 4A", "Flaviviridae") == "NS4A"
    assert chain_family("Non-structural protein 4B", "Flaviviridae") == "NS4B"


@pytest.mark.parametrize("label,expected", CORONA_CHAINS,
                         ids=[c[0] for c in CORONA_CHAINS])
def test_coronavirus_chain_naming_is_unchanged(label, expected):
    """NO REGRESSION, on the real sixteen. The release must come out of this
    commit identical for the seven built viruses."""
    assert chain_family(label, "Coronaviridae") == expected


def test_a_coronavirus_label_is_never_read_by_the_flavivirus_table():
    """The branch is the guarantee, not the patterns. Non-Flaviviridae input
    reaches exactly the code it reached before this commit, so no flavivirus
    pattern can touch a coronavirus row however it is written later."""
    for label, expected in CORONA_CHAINS:
        n = spine._nsp_number(label)
        assert chain_family(label, "Coronaviridae") == (f"nsp{n}" if n else None)
        assert chain_family(label, None) == expected, "default is coronavirus"


def test_the_two_namespaces_are_disjoint():
    """If a key appeared in both, build_viral_target_index being taxon-scoped
    would be the only thing preventing a false identity. Belt and braces: no
    flavivirus key is spelled like a coronavirus one."""
    flavi = {v for _, v in FLAVI_CHAINS if v}
    corona = {v for _, v in CORONA_CHAINS if v}
    assert not flavi & corona
    assert not any(k.lower().startswith("nsp") for k in flavi)


# --------------------------------------------------------------------------
# the activity side
# --------------------------------------------------------------------------
# Quoted from diagnose_flavivirus_targets.py against the real ChEMBL DB.
REAL_ASSAY_DESCRIPTIONS = [
    ("Binding affinity to Zika virus NS2B-NS3 protease", "NS3", "protease"),
    ("Inhibition of Zika virus NS2B/NS3 protease", "NS3", "protease"),
    ("Inhibition of C-terminal His6-tagged linked Zika virus NS2B-NS3 "
     "protease domain expressed", "NS3", "protease"),
    ("Inhibition of N-terminal hexa-histidine tagged Zika virus NS2B/NS3 "
     "protease expressed in E", "NS3", "protease"),
    ("Inhibition of Zika virus unlinked NS2B-NS3 protease preincubated for "
     "5 mins followed by be", "NS3", "protease"),
    ("Inhibition of Dengue virus 2 NS3 helicase incubated for 2 mins "
     "followed by ATP addition", "NS3", "helicase"),
    ("Non-competitive inhibition of Dengue virus 4 NS3 helicase",
     "NS3", "helicase"),
    ("Inhibition of RNA unwinding activity of thioredoxin/His6 tagged "
     "Dengue virus 2 TSV01 NS3 h", "NS3", "helicase"),
    ("Inhibition of Dengue virus type 2 NS5 RNA methyltransferase SAM site",
     "NS5", "methyltransferase"),
    ("Inhibition of Dengue virus type 2 NS5 RNA methyltransferase at 8 nM "
     "enzyme concentration", "NS5", "methyltransferase"),
    ("Inhibition of ZIKV 2'-O-MTase assessed as reduction in [3H]AdoMet "
     "transfer on to capped Gp", "NS5", "methyltransferase"),
    ("PUBCHEM_BIOASSAY: Primary and Confirmatory Screening for Flavivirus "
     "Genomic Capping Enzyme", "NS5", "methyltransferase"),
    ("Inhibition of Zika virus RNA-dependent RNA polymerase using "
     "poly(rC)/rG13 as substrate aft", "NS5", "rdrp"),
    ("Binding affinity to ZIKV RdRP by SPR assay", "NS5", "rdrp"),
    ("Binding affinity to N-terminal His-tagged recombinant ZIKV ZG-01 RdRp "
     "(272 to 903 residues", "NS5", "rdrp"),
    ("Inhibition of Zika virus wild type RNA-dependent RNA polymerase using "
     "ATP as substrate", "NS5", "rdrp"),
]


@pytest.mark.parametrize("desc,family,domain", REAL_ASSAY_DESCRIPTIONS,
                         ids=[d[0][:48] for d in REAL_ASSAY_DESCRIPTIONS])
def test_real_flavivirus_assay_descriptions_resolve(desc, family, domain):
    assert resolve_nsp(desc, "Flaviviridae") == (family, domain)


def test_the_protease_rows_now_reach_ns3_rather_than_the_polyprotein():
    """The number this commit exists to move. 315 of the 355 Zika rows that
    reach the chain lookup name a protease and resolved to nothing before;
    they now resolve to NS3, which is a chain the spine emits."""
    protease_rows = [d for d, f, dom in REAL_ASSAY_DESCRIPTIONS
                     if dom == "protease"]
    assert protease_rows
    for d in protease_rows:
        assert resolve_nsp(d, "Flaviviridae")[0] == "NS3"
        # What it did before: nothing at all.
        assert resolve_nsp(d)[0] is None, \
            "precondition: the coronavirus table could not read this"


def test_the_two_bifunctional_chains_are_split_by_domain_not_by_node():
    """NS3 is protease + helicase; NS5 is methyltransferase + RdRp. UniProt
    annotates one chain each, so the domain qualifier carries the distinction
    rather than inventing nodes the source does not have."""
    assert resolve_nsp("Zika virus NS2B-NS3 protease", "Flaviviridae") \
        == ("NS3", "protease")
    assert resolve_nsp("Dengue virus 2 NS3 helicase", "Flaviviridae") \
        == ("NS3", "helicase")
    assert resolve_nsp("DENV NS5 methyltransferase", "Flaviviridae") \
        == ("NS5", "methyltransferase")
    assert resolve_nsp("ZIKV RdRp", "Flaviviridae") == ("NS5", "rdrp")


def test_the_methyltransferase_rule_precedes_the_rdrp_rule():
    """Order matters and is load-bearing: the bare \\bns5\\b fallback in the
    RdRp rule would otherwise claim "NS5 RNA methyltransferase" and file an
    MTase measurement as a polymerase one. DENV-2's largest single assay --
    1,013 rows, the PubChem capping-enzyme screen -- is MTase."""
    pats = [p.pattern for p, _, _ in FLAVI_PATTERNS]
    mt = next(i for i, p in enumerate(pats) if "mtase" in p)
    rd = next(i for i, p in enumerate(pats) if "rdrp" in p)
    assert mt < rd
    assert resolve_nsp("NS5 RNA methyltransferase SAM site", "Flaviviridae") \
        == ("NS5", "methyltransferase")


def test_coronavirus_descriptions_are_unchanged_by_the_new_argument():
    """resolve_nsp gained an argument. Defaulted or given a coronavirus
    family, every existing caller must get exactly what it got before."""
    for desc in ("Inhibition of SARS-CoV-2 3CL protease",
                 "Inhibition of SARS-CoV-2 Mpro",
                 "Inhibition of SARS-CoV-2 RdRp",
                 "Inhibition of SARS-CoV-2 nsp13 helicase",
                 "Inhibition of SARS-CoV PLpro",
                 "Inhibition of nsp14 exoribonuclease",
                 "Inhibition of SARS-CoV-2 nsp16 2'-O-methyltransferase"):
        assert resolve_nsp(desc) == resolve_nsp(desc, "Coronaviridae")
        assert resolve_nsp(desc) == resolve_nsp(desc, None)
        assert resolve_nsp(desc)[0] is not None, desc


def test_an_unknown_family_falls_back_to_the_coronavirus_table():
    """Not a guess -- the pre-existing behaviour. A taxon with no family in
    the graph gets what it got before this commit, never a flavivirus read."""
    assert resolve_nsp("Inhibition of SARS-CoV-2 Mpro", "Picornaviridae") \
        == resolve_nsp("Inhibition of SARS-CoV-2 Mpro")
    assert len(NSP_PATTERNS) == 9, "the coronavirus table is untouched"


# --------------------------------------------------------------------------
# the two sides agree, which is the point
# --------------------------------------------------------------------------
def test_the_join_key_matches_end_to_end():
    """Simulates build_viral_target_index: key the real chains, then look up
    what the real descriptions resolve to. This is the test that a one-sided
    change breaks."""
    index = {chain_family(label, "Flaviviridae"): f"UniProtKB:PRO_{i}"
             for i, (label, _) in enumerate(FLAVI_CHAINS)
             if chain_family(label, "Flaviviridae")}

    for desc, expected_family, _domain in REAL_ASSAY_DESCRIPTIONS:
        fam, _ = resolve_nsp(desc, "Flaviviridae")
        assert fam == expected_family
        assert fam in index, (
            f"{desc!r} resolves to {fam!r}, which no chain supplies -- the "
            f"row would fall through to the undifferentiated polyprotein")


def test_ns5_no_longer_depends_on_a_coronavirus_coincidence():
    """Before this commit NS5 joined only because "RNA-directed RNA
    polymerase" appears in both families' vocabularies and both sides
    independently said nsp12. That accident is now an explicit rule on each
    side, and neither says nsp12."""
    assert chain_family(
        "RNA-directed RNA polymerase/Methyltransferase NS5", "Flaviviridae") \
        == "NS5"
    assert resolve_nsp("Zika virus RNA-dependent RNA polymerase",
                       "Flaviviridae") == ("NS5", "rdrp")
    # The coincidence, still there in the coronavirus table where it belongs.
    assert spine._nsp_number(
        "RNA-directed RNA polymerase/Methyltransferase NS5") == 12


def test_the_diagnostics_chain_vocabulary_matches_what_the_spine_emits():
    """A SECOND HARDCODED LIST, pinned to the first.

    diagnose_flavivirus_targets.CHAIN_VOCABULARY is the set of keys a chain
    can supply, and the diagnostic uses it to separate "resolved" from
    "would actually join" -- a description resolving to a key no chain carries
    leaves the row on the polyprotein exactly as an unresolved one does, so
    counting them together would overstate the fix.

    That set must equal what chain_family really returns. This repository has
    already had four hand-built taxon maps drift apart before declared_taxa
    replaced them; this is the same shape of mistake one layer down, and the
    cost of it here is a diagnostic that reports a join rate the ingest will
    not reproduce.
    """
    spec = importlib.util.spec_from_file_location(
        "diag_flavi", ROOT / "scripts" / "diagnose_flavivirus_targets.py")
    d = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = d
    spec.loader.exec_module(d)

    emitted = {chain_family(label, "Flaviviridae")
               for label, _ in FLAVI_CHAINS} - {None}
    assert d.CHAIN_VOCABULARY == emitted, (
        f"drifted\n  diagnostic: {sorted(d.CHAIN_VOCABULARY)}\n"
        f"  spine:      {sorted(emitted)}")


def test_every_key_the_activity_side_can_return_is_one_a_chain_supplies():
    """The join's completeness, from the other direction. If a FLAVI_PATTERNS
    entry returns a key no chain carries, every row it matches resolves and
    then silently fails to join -- which looks like success in the resolve
    counters and like nothing in the graph."""
    returnable = {nsp for _p, nsp, _d in FLAVI_PATTERNS}
    emitted = {chain_family(label, "Flaviviridae")
               for label, _ in FLAVI_CHAINS} - {None}
    orphans = returnable - emitted
    assert not orphans, (
        f"these keys resolve but no chain supplies them: {sorted(orphans)}")


# --------------------------------------------------------------------------
# the silent-zero bug: zero rows for a real taxon must not look like zero
# rows for no taxon
# --------------------------------------------------------------------------
def test_query_activities_refuses_an_empty_taxon_list():
    """WHAT ACTUALLY WENT WRONG on 2026-10-11, pinned at the source.

    With no taxa the WHERE clause became `tax_id IN ()`. Standard SQL rejects
    that; SQLite accepts it as always-false. So the query returned zero rows,
    and diagnose_flavivirus_targets printed "DENV-2: 0 activity rows, nothing
    to resolve" for a virus it had itself measured at 8,748 rows an hour
    earlier. No exception, no warning -- a wrong measurement presented as a
    finding, which is the worst failure mode a diagnostic has.

    Raising here fixes it for every caller, not just the one that hit it.
    """
    from kgav.chembl import query_activities

    with pytest.raises(ValueError, match="no taxa"):
        query_activities(Path("/nonexistent/never-opened.db"), [])


def test_declared_taxa_still_defaults_to_the_built_section():
    """The fix that was NOT applied, and why.

    Making declared_taxa include `surveyed:` would have fixed the diagnostic
    in one line and broken the build: ingest_chembl, ingest_vhppi and
    ingest_selectivity all call it, and widening it makes them query DENV-2
    and Zika activity and emit edges to taxa with no protein nodes -- the
    1,881 DANGLING edges of the premature promotion, reached from the other
    direction. The caller passes the section it means.
    """
    import yaml

    from kgav.virus_register import declared_taxa

    cfg = yaml.safe_load((ROOT / "config" / "viruses.yaml").read_text())
    built = declared_taxa(cfg)
    assert len(built) == 11, "the ingests' map must stay the built one"
    for surveyed_taxon in ("11060", "64320", "31634"):
        assert surveyed_taxon not in built, (
            f"{surveyed_taxon} is surveyed; in the ingests' map it would "
            f"produce edges to a taxon with no nodes")

    surveyed = declared_taxa(cfg, "surveyed")
    assert {"11060", "31634", "64320"} <= set(surveyed)


def test_every_flavivirus_in_the_register_resolves_at_least_one_taxon():
    """The invariant the diagnostic now checks at runtime, checked here too.

    A virus that resolves no taxa cannot be measured, and before this it was
    reported as having no data. Covers both sections, because the diagnostic
    reads both.
    """
    import yaml

    from kgav.virus_register import declared_taxa

    cfg = yaml.safe_load((ROOT / "config" / "viruses.yaml").read_text())
    default_family = (cfg.get("defaults") or {}).get("family")
    maps = {sec: declared_taxa(cfg, sec) for sec in ("viruses", "surveyed")}

    for section in ("viruses", "surveyed"):
        for v in cfg.get(section) or []:
            if (v.get("family") or default_family) != "Flaviviridae":
                continue
            taxa = [t for t, sp in maps[section].items()
                    if sp == str(v["taxon"])]
            assert taxa, (
                f"{v['name']} ({section}) resolves no taxa, so it would be "
                f"reported as having no activity data")


# --------------------------------------------------------------------------
# the residue the real database left, and the restraint it required
# --------------------------------------------------------------------------
@pytest.mark.parametrize("desc,expected", [
    ("Inhibitory activity against CF40.NS3pro from Dengue virus type 2",
     ("NS3", "protease")),
    ("Inhibition of Zika virus MTase (4 to 278 residues) expressed in "
     "Escherichia coli T7", ("NS5", "methyltransferase")),
    ("Binding affinity to Zika virus Methyltransferase",
     ("NS5", "methyltransferase")),
    ("Inhibition of West Nile virus protease", ("NS3", "protease")),
])
def test_the_nineteen_rows_the_real_database_left_unresolved(desc, expected):
    """Measured, not guessed. After the two-sided naming landed, the real
    ChEMBL database left exactly 19 unresolved protein rows across the three
    flaviviruses that have any -- DENV-2 12, Zika 6, West Nile 1 -- and all 19
    were these three phrasings (diagnose_flavivirus_targets.py, 2026-10-11).

    Safe only because of a fact about the family: a flavivirus encodes exactly
    ONE viral protease, the NS2B-NS3 heterodimer, and exactly ONE
    methyltransferase, NS5's N-terminal domain. The bare function name is
    unambiguous within Flaviviridae. It is NOT within Coronaviridae, which has
    nsp3 PLpro and nsp5 3CLpro -- one more reason the tables stay apart.
    """
    assert resolve_nsp(desc, "Flaviviridae") == expected


@pytest.mark.parametrize("desc", [
    "Inhibition of furin protease",
    "Inhibition of human cathepsin B protease",
    "Inhibition of TMPRSS2 protease activity",
    "Inhibition of trypsin-like protease",
])
def test_host_proteases_still_fall_through(desc):
    """THE RESTRAINT, pinned. A bare \\bprotease\\b rule would have gained one
    row -- "Inhibition of West Nile virus protease" -- and started claiming
    host protease assays as the viral NS3. Host proteases are real drug
    targets in flavivirus papers (furin cleaves prM), so those rows exist, and
    filing them on NS3 would put an INHIBITS edge on a viral protein for a
    measurement made on a human one.

    Requiring "virus" or "viral" adjacent keeps the one real row and leaves
    these alone. If someone later widens this to a bare word, these fail.
    """
    assert resolve_nsp(desc, "Flaviviridae") == (None, None)


def test_the_specific_rules_still_win_over_the_new_general_ones():
    """Order again. The new bare-ish rules sit after the specific ones, so a
    description naming both a protein and a function keeps the protein."""
    assert resolve_nsp("Inhibition of Dengue virus type 2 NS5 RNA "
                       "methyltransferase SAM site", "Flaviviridae") \
        == ("NS5", "methyltransferase")
    assert resolve_nsp("Inhibition of Dengue virus 2 NS3 helicase",
                       "Flaviviridae") == ("NS3", "helicase")
    assert resolve_nsp("Inhibition of Zika virus NS2B-NS3 protease",
                       "Flaviviridae") == ("NS3", "protease")


def test_the_new_patterns_do_not_touch_coronavirus_rows():
    """nsp3 is a papain-like protease and nsp5 a 3C-like protease, so within
    Coronaviridae a bare "viral protease" is genuinely ambiguous -- which is
    why these rules live in the flavivirus table only."""
    for desc in ("Inhibition of SARS-CoV-2 viral protease",
                 "Inhibition of SARS-CoV-2 virus methyltransferase"):
        before = resolve_nsp(desc)
        assert before == resolve_nsp(desc, "Coronaviridae")
        assert before[0] != "NS3" and before[0] != "NS5"
