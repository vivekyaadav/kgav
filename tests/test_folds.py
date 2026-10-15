"""Gate for phase 2: the two-mode comparison, and the rule that a shared fold
does not license transfer across a different nucleophile.
"""

from dataclasses import dataclass
from dataclasses import field as dc_field
from pathlib import Path

from kgav.emit import Emit
from kgav.folds import NUCLEOPHILES, ChainSpec, ResolvedChain, emit_classes, gate, motif_c_present


def _chain(label, family, role, cat_type, nucleophile, classes):
    spec = ChainSpec(label=label, accession="X", patterns=("x",), family=family,
                     role=role, virus_taxon="NCBITaxon:1")
    # chain_pro_id is part of being resolved: without it .ok is False and a
    # fixture would silently stop being usable by anything that filters on
    # it. A made-up but well-formed PRO_ id keeps these realistic.
    return ResolvedChain(spec=spec, chain_name=label, start=1, end=100,
                         chain_pro_id=f"PRO_{abs(hash(label)) % 10**7:07d}",
                         catalytic_type=cat_type, nucleophile=nucleophile,
                         classes=[{"accession": a, "label": n, "db": "interpro",
                                   "level": "homologous_superfamily",
                                   "start": 1, "end": 100}
                                  for a, n in classes])


# Measured in phase 0: all three families share SSF50494 and G3DSA:2.40.10.10.
SHARED_PROTEASE_FOLD = [("SSF50494", "Trypsin-like serine protease"),
                        ("G3DSA:2.40.10.10", "Trypsin-like serine protease")]
SHARED_RDRP = [("IPR043502", "DNA/RNA polymerase superfamily"),
               ("IPR007094", "RNA-directed RNA polymerase, catalytic domain")]

COV_MPRO = _chain("SARS-CoV-2 Mpro", "Coronaviridae", "protease",
                  "cysteine", "C3408", SHARED_PROTEASE_FOLD)
PICO_3C = _chain("Poliovirus 3C", "Picornaviridae", "protease",
                 "cysteine", "C1712", SHARED_PROTEASE_FOLD)
HCV_NS3 = _chain("HCV NS3", "Flaviviridae", "protease",
                 "serine", "S1165", SHARED_PROTEASE_FOLD)
COV_RDRP = _chain("SARS-CoV-2 nsp12", "Coronaviridae", "polymerase",
                  "none", "", SHARED_RDRP)
PICO_3D = _chain("Poliovirus 3D", "Picornaviridae", "polymerase",
                 "none", "", SHARED_RDRP)


# ------------------------------------------------------- nucleophile mode
def test_same_nucleophile_and_shared_fold_is_allowed():
    """Coronavirus and picornavirus proteases: cysteine both, shared fold.
    This is the transfer the conservation hypothesis predicts."""
    g = gate(COV_MPRO, PICO_3C)
    assert g.allowed and g.mode == "nucleophile" and g.strength == "strong"
    assert "SSF50494" in g.shared


def test_different_nucleophile_is_refused_despite_identical_fold():
    """THE central check. Phase 0 measured that SARS-CoV-2 Mpro and HCV NS3
    share SSF50494 and G3DSA:2.40.10.10 -- by fold they are indistinguishable.
    Without this gate M8 recommends nirmatrelvir for hepatitis C."""
    g = gate(COV_MPRO, HCV_NS3)
    assert not g.allowed
    assert g.catalytic_type_match is False
    assert "cysteine" in g.reason and "serine" in g.reason
    # the shared fold is still reported, so the refusal is explicable
    assert "SSF50494" in g.shared


def test_same_nucleophile_without_shared_fold_is_refused():
    lonely = _chain("other", "Other", "protease", "cysteine", "C10",
                    [("IPR999999", "unrelated")])
    g = gate(COV_MPRO, lonely)
    assert not g.allowed and g.strength == "none"


# ------------------------------------------------------------- motif mode
def test_polymerases_are_compared_in_motif_mode_and_marked_weak():
    """Polymerases have no catalytic nucleophile, and motif C (GDD/SDD) is
    present in all three families -- it does not discriminate the way Cys vs
    Ser does. The result must say so."""
    g = gate(COV_RDRP, PICO_3D)
    assert g.allowed and g.mode == "motif"
    assert g.strength == "weak"
    assert "WEAK" in g.reason


def test_motif_mode_reports_catalytic_domain_entries_when_present():
    g = gate(COV_RDRP, PICO_3D)
    assert "IPR007094" in g.shared


# ------------------------------------------------------------- mixed mode
def test_protease_versus_polymerase_is_refused_outright():
    """A protease and a polymerase are not comparable targets; any fold
    similarity between them is an artifact."""
    g = gate(COV_MPRO, COV_RDRP)
    assert not g.allowed and g.mode == "refused"
    assert "not comparable" in g.reason


def test_refusal_is_symmetric():
    assert gate(COV_MPRO, COV_RDRP).allowed is gate(COV_RDRP, COV_MPRO).allowed
    assert gate(COV_MPRO, HCV_NS3).allowed is gate(HCV_NS3, COV_MPRO).allowed


# ---------------------------------------------------------------- motif C
def test_motif_c_is_a_check_not_a_detector():
    """The pattern also matches CDD in the coronavirus NiRAN domain, which is
    not motif C. It may only be used on a region already known to be a
    polymerase."""
    seq = "AAAAGDDAAAA"
    assert motif_c_present(seq, 1, len(seq)) == ["GDD@5"]
    # SDD, the coronavirus variant, is matched
    assert motif_c_present("AAAASDDAAAA", 1, 11) == ["SDD@5"]
    # CDD is NOT matched by the check pattern
    assert motif_c_present("AAAACDDAAAA", 1, 11) == []


# ------------------------------------------------------------ definitions
def test_nucleophile_table_maps_to_schema_enum_values():
    from kgav.schema import load_schema
    allowed = set(load_schema().enums["catalytic_type"])
    assert set(NUCLEOPHILES.values()) <= allowed


def test_every_phase0_protease_pair_behaves_as_measured():
    """The full phase 0 matrix, encoded so a regression is visible."""
    cys = [COV_MPRO, PICO_3C]
    ser = [HCV_NS3]
    for a in cys:
        for b in cys:
            assert gate(a, b).allowed, f"{a.chain_name} vs {b.chain_name}"
    for a in cys:
        for b in ser:
            assert not gate(a, b).allowed, f"{a.chain_name} vs {b.chain_name}"


# ------------------------------- a class whose members disagree is not a gate
@dataclass
class _Chain:
    node_id: str
    catalytic_type: str
    ok: bool = True
    classes: list = dc_field(default_factory=list)


_SHARED = [{"accession": "IPR009003", "label": "Peptidase S1, PA clan",
            "level": "homologous_superfamily", "start": 1, "end": 200}]


def _emit(order):
    em = Emit()
    stats = emit_classes(em, order, "infores:interpro")
    return em.nodes["INTERPRO:IPR009003"]["properties"], stats


def test_a_mixed_nucleophile_class_is_marked_unknown_not_guessed():
    """Emit.node merges with properties.update -- last write wins, silently,
    recording only a node_merged counter. So a class holding a cysteine
    protease and a serine protease took whichever chain was written last,
    and reversing the order flipped it from serine to cysteine.

    That is the field M8's gate reads. The schema warns that without the
    gate M8 "traverses from SARS-CoV-2 Mpro to HCV NS3 -- same superfamily,
    opposite nucleophile -- and recommends nirmatrelvir for hepatitis C",
    and those two land in the same InterPro superfamily. A gate whose value
    depends on emission order is not a gate.
    """
    mpro = _Chain("UniProtKB:MPRO", "cysteine", classes=list(_SHARED))
    ns3 = _Chain("UniProtKB:NS3", "serine", classes=list(_SHARED))
    fwd, s1 = _emit([mpro, ns3])
    rev, s2 = _emit([ns3, mpro])
    assert fwd["catalytic_type"] == "unknown"
    assert rev["catalytic_type"] == "unknown"
    assert s1["classes_mixed_catalytic_type"] == 1
    assert s2["classes_mixed_catalytic_type"] == 1


def test_emission_order_cannot_change_a_class_catalytic_type():
    """The property the previous code lacked, stated directly."""
    a = _Chain("UniProtKB:A", "cysteine", classes=list(_SHARED))
    b = _Chain("UniProtKB:B", "serine", classes=list(_SHARED))
    c = _Chain("UniProtKB:C", "metallo", classes=list(_SHARED))
    assert _emit([a, b, c])[0] == _emit([c, b, a])[0]


def test_an_agreeing_class_keeps_its_catalytic_type():
    """So "unknown" marks disagreement rather than being a blanket refusal --
    a gate that never opens is as useless as one that never closes."""
    props, stats = _emit([_Chain("UniProtKB:A", "cysteine", classes=list(_SHARED)),
                          _Chain("UniProtKB:B", "cysteine", classes=list(_SHARED))])
    assert props["catalytic_type"] == "cysteine"
    assert stats["classes_mixed_catalytic_type"] == 0


def test_unknown_is_a_declared_value_of_the_enum():
    """Marking a class unknown must be sayable in the schema, not a string
    the validator will reject the moment the folds layer is assembled."""
    from pathlib import Path

    import yaml
    root = Path(__file__).resolve().parents[1]
    d = yaml.safe_load((root / "schema/kg_schema.yaml").read_text())
    assert "unknown" in d["enums"]["catalytic_type"]


# ------------------------------------------------- the id that joins the spine
def test_node_id_is_the_spine_s_pro_id_not_a_coordinate_key():
    """THE REASON THIS LAYER NEVER REACHED THE GRAPH. node_id was
    UniProtKB:{accession}_{start}_{end} while ingest_spine keys the same
    chain on its PRO_ id. Two schemes for one protein, so every
    MEMBER_OF_CLASS edge pointed at a node that does not exist -- and
    assemble_release returns 1 on a dangling edge, which is why `folds` is
    absent from LAYER_PRECEDENCE.
    """
    spec = ChainSpec(label="Mpro", accession="P0DTD1", patterns=("3C-like",),
                     family="Coronaviridae", role="protease",
                     virus_taxon="NCBITaxon:2697049")
    r = ResolvedChain(spec=spec, chain_name="3C-like proteinase",
                      start=3264, end=3569, chain_pro_id="PRO_0000449623")
    assert r.ok
    assert r.node_id == "UniProtKB:PRO_0000449623"
    assert "3264" not in r.node_id, "coordinates must not appear in the key"


def test_a_chain_without_a_pro_id_is_not_resolved():
    """Same rule the spine applies: no PRO_ id, no stable node key. Being
    part of `ok` matters because `ok` is what every caller filters on -- a
    chain that is merely 'not emitted later' would still be compared,
    counted and reported as usable."""
    spec = ChainSpec(label="X", accession="Q1", patterns=("x",),
                     family="Flaviviridae", role="protease",
                     virus_taxon="NCBITaxon:64320")
    r = ResolvedChain(spec=spec, chain_name="x", start=1, end=50)
    assert not r.ok


def test_reading_node_id_without_a_pro_id_raises_rather_than_inventing_one():
    """Returning a coordinate key here is what the previous version did, and
    it failed silently all the way through to a dangling edge."""
    import pytest as _pytest

    spec = ChainSpec(label="X", accession="Q1", patterns=("x",),
                     family="Flaviviridae", role="protease",
                     virus_taxon="NCBITaxon:64320")
    r = ResolvedChain(spec=spec, chain_name="x", start=1, end=50)
    with _pytest.raises(ValueError, match="no PRO_ id"):
        _ = r.node_id


def test_find_chain_keeps_the_feature_id_uniprot_returns():
    """Verified against rest.uniprot.org for P0DTD1 on 2026-10-15 with this
    exact field selection -- a field present in a proteome download is not
    evidence about a single-entry query, so it was checked separately."""
    from kgav import folds

    payload = {"features": [
        {"type": "Chain", "featureId": "PRO_0000449618",
         "description": "Replicase polyprotein 1ab",
         "location": {"start": {"value": 1}, "end": {"value": 7096}}},
        {"type": "Chain", "featureId": "PRO_0000449623",
         "description": "3C-like proteinase nsp5",
         "location": {"start": {"value": 3264}, "end": {"value": 3569}}},
    ]}
    real, folds._get = folds._get, lambda _u: payload
    try:
        name, s, e, note, pro = folds.find_chain("P0DTD1", ("3C-like",))
    finally:
        folds._get = real
    assert (s, e, note, pro) == (3264, 3569, "", "PRO_0000449623")
    assert "3C-like" in name


def test_a_matched_chain_with_no_pro_id_is_reported_not_keyed():
    from kgav import folds

    payload = {"features": [
        {"type": "Chain", "description": "mystery chain",
         "location": {"start": {"value": 1}, "end": {"value": 9}}}]}
    real, folds._get = folds._get, lambda _u: payload
    try:
        _name, s, e, note, pro = folds.find_chain("Q1", ("mystery",))
    finally:
        folds._get = real
    assert (s, e, pro) == (1, 9, "")
    assert "no PRO_ id" in note


def test_role_decides_mixed_mode_not_whether_a_nucleophile_resolved():
    """THE GATE HELD BY LUCK. MIXED was decided by a_nuc != b_nuc -- whether
    a nucleophile was FOUND, a property of the sequence scan. SARS-CoV-2
    PLpro resolves catalytic_type "none" because its spec matches UniProt's
    "Papain-like protease nsp3", the whole 1,945-residue multi-domain chain
    rather than the protease domain. So a PROTEASE reached MOTIF mode beside
    the polymerases, where an allowed weak transfer protease->polymerase was
    possible -- the one comparison the docstring calls an artifact.

    It never fired on v0.1 only because PLpro shares no catalytic-domain
    entry with any cross-family polymerase.
    """
    plpro = _chain("SARS-CoV-2 PLpro", "Coronaviridae", "protease",
                   "none", "", SHARED_RDRP)      # same classes as an RdRp
    g = gate(plpro, PICO_3D)
    assert not g.allowed, "a protease must not transfer to a polymerase"
    assert g.mode == "refused"
    assert "different target roles" in g.reason

    # And the old rule alone would have let it through: both lack a
    # nucleophile, and they share catalytic-domain entries.
    assert plpro.catalytic_type == PICO_3D.catalytic_type == "none"
    assert set(c["accession"] for c in plpro.classes) & \
           set(c["accession"] for c in PICO_3D.classes)


def test_matching_roles_are_still_compared():
    """The role check must not refuse everything: a gate that never opens is
    as useless as one that never closes."""
    assert gate(COV_MPRO, PICO_3C).allowed
    assert gate(COV_RDRP, PICO_3D).mode == "motif"


# ------------------------------------------- which rule read the active sites
def _machinery(sites, seq_len=6000):
    """catalytic_machinery against a stubbed UniProt response.

    `sites` is [(position, residue, description)]; the sequence is built so
    each position holds its residue.
    """
    from kgav import folds

    seq = ["A"] * seq_len
    for pos, aa, _d in sites:
        seq[pos - 1] = aa
    payload = {"sequence": {"value": "".join(seq)},
               "features": [{"location": {"start": {"value": p}},
                             "description": d} for p, _aa, d in sites]}
    real, folds._get = folds._get, lambda _u: payload
    try:
        return folds.catalytic_machinery("X", 1, seq_len)
    finally:
        folds._get = real


def test_an_annotated_nucleophile_decides_and_says_so():
    _res, nuc, typ, basis = _machinery(
        [(100, "H", "Proton acceptor; for 3CL-PRO activity"),
         (200, "C", "Nucleophile; for 3CL-PRO activity")])
    assert (nuc, typ, basis) == ("C200", "cysteine", "annotated")


def test_composition_decides_a_chymotrypsin_like_triad():
    """No annotation, but His-...-Cys with the nucleophile last is still
    readable."""
    _res, nuc, typ, basis = _machinery(
        [(100, "H", "For protease activity"),
         (200, "C", "For protease activity")])
    assert (nuc, typ, basis) == ("C200", "cysteine", "composition")


def test_a_papain_like_triad_is_read_by_composition():
    """THE PLPRO CASE. Its active sites are C1674, H1835, D1849, all
    described "For PL-PRO activity" with no nucleophile annotation. A
    papain-like protease runs Cys-His-Asp with the nucleophile FIRST, and
    the composition rule used to test the LAST residue -- so a cysteine
    protease read "none" and the gate could not evaluate it.
    """
    _res, nuc, typ, basis = _machinery(
        [(1674, "C", "For PL-PRO activity"),
         (1835, "H", "For PL-PRO activity"),
         (1849, "D", "For PL-PRO activity")])
    assert (nuc, typ, basis) == ("C1674", "cysteine", "composition")


def test_a_serine_without_a_histidine_is_still_refused():
    """THE CONDITION THAT KEEPS THE WIDENING HONEST. SARS-CoV-2 RdRp's
    sites are S5151, D5152, D5153: a serine, no histidine. Dropping the
    His test to catch PLpro would read "serine" for a polymerase -- one
    step from the "aspartic D5153" misclassification the explicit branch
    was fixed to prevent.
    """
    _res, nuc, typ, _b = _machinery(
        [(5151, "S", "For polymerase activity"),
         (5152, "D", "For polymerase activity"),
         (5153, "D", "For polymerase activity")], seq_len=6000)
    assert (nuc, typ) == ("", "none")


def test_two_nucleophile_candidates_are_refused_rather_than_picked():
    """Two C/S sites mean the set spans more than one active site and
    nothing here can say which belongs to the protease."""
    _res, nuc, typ, _b = _machinery(
        [(10, "C", "For activity A"), (20, "H", "For activity A"),
         (30, "S", "For activity B")])
    assert (nuc, typ) == ("", "none")


def test_the_widened_rule_moves_only_plpro():
    """THE BLAST RADIUS, measured before the rule changed and pinned here.

    Every chain in SPECS that has active sites, with the residues exactly
    as ingest_folds reported them on 2026-10-15. Fourteen of these were
    already decided by composition; the point of the table is that the
    position-independent rule leaves all fourteen bit-identical and moves
    exactly one chain.
    """
    cases = [
        ("SARS-CoV 3CLpro", "H3281, C3385", "cysteine"),
        ("MERS-CoV 3CLpro", "H3288, C3395", "cysteine"),
        ("HCoV-229E 3CLpro", "H3006, C3109", "cysteine"),
        ("Dengue 2 NS3pro", "H1526, D1550, S1610", "serine"),
        ("Zika NS3pro", "H1549, D1573, S1633", "serine"),
        ("West Nile NS3pro", "H1552, D1576, S1636", "serine"),
        ("HCV NS3pro", "H1083, D1107, S1165", "serine"),
        ("Yellow fever NS3pro", "H1537, D1561, S1622", "serine"),
        ("JEV NS3pro", "H1555, D1579, S1639", "serine"),
        ("EV-A71 3C", "H1588, E1619, C1695", "cysteine"),
        ("Poliovirus 3C", "H1605, E1636, C1712", "cysteine"),
        ("HAV 3C", "H1563, D1603, C1691", "cysteine"),
        ("CVB3 3C", "H1580, E1611, C1687", "cysteine"),
        ("Rhinovirus A 3C", "H1561, E1592, C1668", "cysteine"),
        # the four polymerases that must stay unreadable
        ("SARS-CoV-2 RdRp", "S5151, D5152, D5153", "none"),
        ("SARS-CoV RdRp", "S5128, D5129, D5130", "none"),
        ("Dengue 2 NS5", "K2552, D2637, K2672, E2708", "none"),
        ("Zika NS5", "K2577, D2662, K2698, E2734", "none"),
        # and the one that moves
        ("SARS-CoV-2 PLpro", "N1062, G1070, C1674, H1835, D1849", "cysteine"),
    ]
    for label, residues, expected in cases:
        sites = [(int(r[1:]), r[0], "For activity") for r in residues.split(", ")]
        _res, _nuc, typ, _b = _machinery(sites, seq_len=6000)
        assert typ == expected, f"{label}: {typ} != {expected}"


def test_no_active_sites_at_all_is_distinguishable():
    assert _machinery([])[3] == "none"


# ------------------------------------------- the cache stores input, not verdict
def _ingest():
    import sys
    from importlib import import_module
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    return import_module("ingest_folds")


#: SARS-CoV-2 PLpro, exactly as the 2026-10-15 run printed it. None of the
#: five sites is annotated as a nucleophile, so only the composition rule
#: can read it -- which makes it the chain that moves when that rule moves.
PLPRO_SITES = [(1062, "N", "For ADP-ribosylhydrolase activity"),
               (1070, "G", "For ADP-ribosylhydrolase activity"),
               (1674, "C", "For PL-PRO activity"),
               (1835, "H", "For PL-PRO activity"),
               (1849, "D", "For PL-PRO activity")]


def test_a_cached_verdict_is_recomputed_not_replayed():
    """THE BUG, stated as a test.

    The widened composition rule shipped, the suite went green, and the
    `--offline` run printed "composition 14, none 9, annotated 1" --
    character-identical to the run before it, PLpro still listed under "no
    rule could read them". The rule was never executed: the cache held the
    whole ResolvedChain, derived fields included, and the offline branch
    rebuilt the object from JSON without calling catalytic_machinery.

    So a cache carrying the OLD verdict for PLpro must come back carrying
    the NEW one.
    """
    from kgav.folds import ChainSpec, ResolvedChain

    stale = {"spec": ChainSpec(label="SARS-CoV-2 PLpro", accession="P0DTD1",
                               patterns=("Papain-like proteinase",),
                               family="Coronaviridae", role="protease",
                               virus_taxon="NCBITaxon:2697049").__dict__,
             "chain_name": "Papain-like proteinase", "start": 819, "end": 2763,
             "chain_pro_id": "PRO_0000449623",
             "active_sites": [list(s) for s in PLPRO_SITES],
             # what the previous rule decided, written into the file
             "catalytic_type": "none", "nucleophile": "", "catalytic_basis": "none",
             "catalytic_residues": []}

    loaded, = _ingest().load_cache([stale])
    assert isinstance(loaded, ResolvedChain)
    assert loaded.catalytic_type == "cysteine"
    assert loaded.nucleophile == "C1674"
    assert loaded.catalytic_basis == "composition"


def test_the_cache_file_does_not_contain_a_verdict_to_go_stale():
    """Belt to the braces above: nothing derived is written in the first
    place, so there is no stale value even for code that forgets to
    recompute."""
    from kgav.folds import DERIVED_FIELDS, ChainSpec, ResolvedChain

    c = ResolvedChain(spec=ChainSpec("L", "P1", ("p",), "F", "protease", "T:1"),
                      chain_name="x", start=1, end=9, chain_pro_id="PRO_1",
                      active_sites=[list(s) for s in PLPRO_SITES])
    c.redecide()
    assert c.catalytic_type == "cysteine"          # it IS decided in memory
    written, = _ingest().cache_payload([c])
    assert DERIVED_FIELDS.isdisjoint(written)      # and not on disk
    assert written["active_sites"]                 # the input is


def test_the_completeness_guard_does_not_demand_the_fields_it_dropped():
    """The two halves have to agree: a cache the writer just produced must
    satisfy the reader's guard. They were separate expressions, and this is
    the assertion that keeps them one."""
    from kgav.folds import ChainSpec, ResolvedChain

    ing = _ingest()
    c = ResolvedChain(spec=ChainSpec("L", "P1", ("p",), "F", "protease", "T:1"),
                      chain_name="x", start=1, end=9, chain_pro_id="PRO_1",
                      active_sites=[])
    assert ing.missing_inputs(ing.cache_payload([c])) == []


def test_a_cache_predating_active_sites_is_refused_not_defaulted():
    """An old cache has no sites to re-decide from, so it cannot be
    corrected -- it must be refused rather than silently read as a chain
    with no catalytic machinery."""
    ing = _ingest()
    old = {"spec": {"label": "L", "accession": "P1", "patterns": ["p"],
                    "family": "F", "role": "protease", "virus_taxon": "T:1"},
           "chain_name": "x", "start": 1, "end": 9, "chain_pro_id": "PRO_1",
           "classes": [], "note": ""}
    # Every input field added since, not just the first one. The guard is
    # derived from the dataclass precisely so this list grows by itself.
    assert ing.missing_inputs([old]) == ["active_sites", "domains",
                                         "window_end", "window_note",
                                         "window_start"]


def test_a_failed_fetch_and_an_empty_fetch_stay_different():
    """None is "the network did not answer", [] is "it answered, there are
    no active sites". Collapsing them would report an outage as biology."""
    from kgav.folds import decide_machinery

    assert decide_machinery(None)[2] == "unknown"
    assert decide_machinery([])[2] == "none"


def test_the_rule_is_reachable_without_the_network():
    """decide_machinery must not touch _get -- if it does, the cache path
    is making network calls and the split has quietly undone itself."""
    from kgav import folds

    def _boom(_u):
        raise AssertionError("decide_machinery reached the network")

    real, folds._get = folds._get, _boom
    try:
        assert folds.decide_machinery(PLPRO_SITES)[1] == "C1674"
    finally:
        folds._get = real


# --------------------------------------------- the window the analysis reads
#: P0DTD1's Domain annotations inside nsp3 (819-2763), copied from the
#: UniProt probe run on 2026-10-15. Nine domains in one chain.
NSP3_DOMAINS = [(821, 929, "Ubiquitin-like 1"), (1025, 1194, "Macro 1"),
                (1231, 1359, "Macro 2"), (1367, 1494, "Macro 3"),
                (1496, 1561, "DPUP"), (1565, 1620, "Ubiquitin-like 2"),
                (1634, 1898, "Peptidase C16"),
                (1911, 2021, "Nucleic acid-binding"), (2046, 2155, "G2M"),
                (2247, 2317, "3Ecto"), (2395, 2763, "CoV Nsp3 Y")]


def test_the_window_defaults_to_the_whole_chain():
    """Narrowing is opt-in. A chain that is one enzyme must be unaffected."""
    c = ResolvedChain(spec=ChainSpec("L", "P1", ("p",), "F", "protease", "T:1"),
                      chain_name="x", start=3264, end=3569,
                      chain_pro_id="PRO_1")
    assert c.window == (3264, 3569)
    assert c.window_is_chain


def test_a_named_domain_narrows_the_window_to_itself():
    from kgav.folds import pick_domain

    assert pick_domain(NSP3_DOMAINS, "Peptidase C16") == (1634, 1898, "")


def test_narrowing_cuts_the_macrodomain_sites_out_of_the_catalytic_input():
    """WHY THE WINDOW MATTERS, in the one case measured end to end.

    PLpro's sites over the whole nsp3 chain are N1062, G1070, C1674, H1835,
    D1849 -- and the first two are the MACRODOMAIN's ADP-ribosylhydrolase
    sites, 600 residues from the protease. The rule returned cysteine over
    that set, but only because neither N nor G is a nucleophile candidate.
    Over the real protease domain the same answer no longer depends on
    that.
    """
    from kgav.folds import decide_machinery, pick_domain

    lo, hi, why = pick_domain(NSP3_DOMAINS, "Peptidase C16")
    assert not why
    inside = [s for s in PLPRO_SITES if lo <= s[0] <= hi]
    assert [s[0] for s in inside] == [1674, 1835, 1849]      # the triad alone
    residues, nuc, typ, basis = decide_machinery(inside)
    assert (nuc, typ, basis) == ("C1674", "cysteine", "composition")
    assert residues == ["C1674", "H1835", "D1849"]
    # and the un-narrowed set carried two residues from another enzyme
    assert decide_machinery(PLPRO_SITES)[0][:2] == ["N1062", "G1070"]


def test_a_requested_narrowing_that_fails_does_not_fall_back_silently():
    """The whole chain is the wrong window; reverting to it without saying
    so is how this survived unnoticed in the first place."""
    from kgav.folds import pick_domain

    s, e, why = pick_domain(NSP3_DOMAINS, "Peptidase C30")
    assert (s, e) == (0, 0)
    assert "no Domain matched" in why
    assert "Peptidase C16" in why          # says what IS there

    s, e, why = pick_domain(None, "Peptidase C16")
    assert (s, e) == (0, 0) and "fetch failed" in why


def test_two_matching_domains_are_refused_rather_than_picked():
    doubled = NSP3_DOMAINS + [(2500, 2600, "Peptidase C16")]
    from kgav.folds import pick_domain

    s, e, why = pick_domain(doubled, "Peptidase C16")
    assert (s, e) == (0, 0) and "ambiguous" in why


def test_the_window_is_cached_because_it_is_an_input_not_a_verdict():
    """Everything downstream is a function of the window, so a cache that
    omitted it would re-derive classes and catalytic type against the wrong
    region -- the previous bug with a different field."""
    from kgav.folds import DERIVED_FIELDS

    assert DERIVED_FIELDS.isdisjoint({"window_start", "window_end", "domains"})
    ing = _ingest()
    c = ResolvedChain(spec=ChainSpec("L", "P1", ("p",), "F", "protease", "T:1"),
                      chain_name="x", start=819, end=2763, chain_pro_id="PRO_1",
                      window_start=1634, window_end=1898, active_sites=[])
    written, = ing.cache_payload([c])
    assert written["window_start"] == 1634 and written["window_end"] == 1898
    assert ing.missing_inputs([written]) == []
    assert ing.load_cache([written])[0].window == (1634, 1898)


def test_the_docstring_no_longer_claims_a_split_that_does_not_exist():
    """It said the NS3 protease and helicase were "declared separately
    below, with the boundary taken from the domain annotation". There is
    one NS3 spec per flavivirus and, until this commit, no code read a
    Domain annotation anywhere."""
    src = (Path(__file__).resolve().parents[1]
           / "scripts/ingest_folds.py").read_text()
    head = src.split('"""')[1]
    assert "Both are declared\nseparately below" not in head
    assert "Neither\nclause was ever true" in head


# --------------------------------------- the declared windows, pinned to data
#: Chain -> the UniProt Domain its window narrows to, and the catalytic
#: residues the 2026-10-15 run read. Every entry is from that run's output,
#: not from recall. The residues are here so a wrong domain name cannot be
#: adopted quietly: a window that excludes a chain's own catalytic residues
#: is the wrong window whatever it is called.
#: label -> (domain, residues that MUST be inside the window,
#:           residues deliberately left outside and why)
#: Every number is from a run's printed output, not from recall.
DECLARED_WINDOWS = {
    "SARS-CoV-2 PLpro":    ("Peptidase C16", (1634, 1898),
                            [1674, 1835, 1849], []),
    "Dengue 2 NS3pro":     ("Peptidase S7", (1476, 1653),
                            [1526, 1550, 1610], []),
    "Zika NS3pro":         ("Peptidase S7", (1499, 1676),
                            [1549, 1573, 1633], []),
    "West Nile NS3pro":    ("Peptidase S7", (1502, 1679),
                            [1552, 1576, 1636], []),
    "HCV NS3pro":          ("Peptidase S29", (1027, 1208),
                            [1083, 1107, 1165], []),
    "Yellow fever NS3pro": ("Peptidase S7", (1485, 1665),
                            [1537, 1561, 1622], []),
    "JEV NS3pro":          ("Peptidase S7", (1505, 1682),
                            [1555, 1579, 1639], []),
    # the seven polymerases: "RdRp catalytic" is the only domain name
    # UniProt annotates on all of them
    "SARS-CoV-2 RdRp":     ("RdRp catalytic", (5004, 5166),
                            [5151, 5152, 5153], []),
    "SARS-CoV RdRp":       ("RdRp catalytic", (4981, 5143),
                            [5128, 5129, 5130], []),
    "HCV NS5B":            ("RdRp catalytic", (2634, 2752), [], []),
    # NS5's four annotated sites are the METHYLTRANSFERASE's -- the run
    # placed every one inside "mRNA cap 0-1 NS5-type MT". Excluding them is
    # the point, so they are recorded as required-outside rather than
    # quietly dropped from the table.
    "Dengue 2 NS5":        ("RdRp catalytic", (3020, 3169), [],
                            [2552, 2637, 2672, 2708]),
    "Zika NS5":            ("RdRp catalytic", (3045, 3195), [],
                            [2577, 2662, 2698, 2734]),
    "Poliovirus 3D":       ("RdRp catalytic", (1975, 2090), [], []),
    "EV-A71 3D":           ("RdRp catalytic", (1958, 2073), [], []),
}


def test_every_narrowed_spec_narrows_to_the_domain_it_was_measured_for():
    """A BULK EDIT PUT "Peptidase S7" ON A PICORNAVIRUS.

    Six specs needed the same two-token addition, so I did it with a regex.
    It matched the wrong lines: two NS3pro specs were skipped, JEV got HCV's
    S29, and EV-A71 3C -- a picornavirus 3C protease with no S7 domain
    anywhere in its entry -- was given domain="Peptidase S7". pick_domain
    would have refused it at runtime and printed WINDOW NOT NARROWED, so it
    would not have corrupted a result; but it would have reached the user.

    The table is the record of what each window was MEASURED to be. Nothing
    may carry a domain that is not in it, and nothing in it may lose one.
    """
    specs = {s.label: s for s in _ingest().SPECS}
    for label, spec in specs.items():
        want = DECLARED_WINDOWS.get(label, ("",))[0]
        assert spec.domain == want, (
            f"{label}: domain={spec.domain!r}, measured {want!r}")


def test_no_narrowed_window_excludes_its_own_catalytic_residues():
    """The name alone is not enough: a plausible name can resolve to a
    region that does not hold the enzyme's triad. That window is wrong even
    though the narrowing "worked"."""
    for label, (dom, (lo, hi), inside, _out) in DECLARED_WINDOWS.items():
        missing = [p for p in inside if not lo <= p <= hi]
        assert not missing, f"{label}: {missing} outside {dom} {lo}-{hi}"


def test_residues_recorded_as_excluded_really_are_excluded():
    """NS5's four annotated sites belong to its methyltransferase, and
    narrowing to the polymerase is supposed to drop them. If a future
    boundary quietly pulled them back in, the polymerase's catalytic type
    would once again be read from a capping enzyme."""
    checked = 0
    for label, (dom, (lo, hi), _in, out) in DECLARED_WINDOWS.items():
        for p in out:
            assert not lo <= p <= hi, (
                f"{label}: {p} was recorded as deliberately outside "
                f"{dom} {lo}-{hi}, but falls inside it")
            checked += 1
    assert checked == 8, "both NS5 chains' four MTase sites"


def test_the_polymerase_windows_all_use_the_one_shared_domain_name():
    """Narrowing coronaviruses to the 568 aa "Nsp12 RNA-dependent RNA
    polymerase" and picornaviruses to the 116 aa "RdRp catalytic" would
    compare unlike regions as if alike. "RdRp catalytic" was chosen because
    intersecting the seven entries' domain-name sets leaves exactly it."""
    poly = [s for s in _ingest().SPECS if s.role == "polymerase"]
    assert len(poly) == 7
    assert {s.domain for s in poly} == {"RdRp catalytic"}
    spans = [DECLARED_WINDOWS[s.label][1] for s in poly]
    widths = [hi - lo + 1 for lo, hi in spans]
    assert max(widths) / min(widths) < 1.5, dict(zip(
        [s.label for s in poly], widths))


def test_a_spec_cannot_request_a_domain_no_chain_of_its_family_has():
    """EV-A71 3C with domain="Peptidase S7" is the shape of the bug: a
    plausible name from the wrong family. Nothing outside the measured table
    requests a narrowing at all."""
    for s in _ingest().SPECS:
        if s.domain:
            assert s.label in DECLARED_WINDOWS, (
                f"{s.label} requests domain={s.domain!r} but no run has "
                f"measured that domain on its entry")


# ------------------------------------- the report accounts for every chain
def _rc(label, role, start, end, *, domains=None, sites=None,
        window=None, classes=0):
    c = ResolvedChain(
        spec=ChainSpec(label, "P1", ("p",), "F", role, "T:1"),
        chain_name=label, start=start, end=end, chain_pro_id="PRO_1",
        domains=domains, active_sites=sites,
        classes=[{"accession": f"X{i}"} for i in range(classes)])
    if window:
        c.window_start, c.window_end = window
    return c


def test_the_chain_line_reports_the_window_it_actually_read():
    """It printed the CHAIN's coordinates beside the WINDOW's class count:
    "[1476-2093] ... 11 classes" when the 11 came from 1476-1653. Two
    numbers from different regions, side by side."""
    lab = _ingest().window_label
    assert lab(_rc("3CLpro", "protease", 3264, 3569)) == "[3264-3569]"
    narrowed = _rc("Dengue 2 NS3pro", "protease", 1476, 2093,
                   window=(1476, 1653))
    assert lab(narrowed) == "[1476-1653 of 1476-2093]"


def test_a_chain_with_no_annotated_domain_is_not_dropped_from_the_table():
    """HAV 3C: three catalytic residues, zero annotated Domains. The filter
    required `c.domains`, so the one report whose job is to say where every
    active site sits omitted it entirely -- and said nothing about doing
    so."""
    hav = _rc("HAV 3C", "protease", 1520, 1738, domains=[],
              sites=[(1563, "H", ""), (1603, "D", ""), (1691, "C", "")])
    rows = _ingest().active_site_rows([hav])
    assert [r[2] for r in rows] == ["H1563", "D1603", "C1691"]
    assert all(r[3] == "(entry annotates no Domain)" for r in rows)


def test_a_site_outside_every_domain_reads_differently_from_no_domains():
    """"The entry annotates nothing" and "the entry annotates domains and
    this residue is in none of them" are different findings about the
    curation, and the second is the one that suggests a wrong window."""
    c = _rc("X", "protease", 100, 400,
            domains=[(100, 200, "Peptidase Q1")],
            sites=[(150, "C", ""), (350, "H", "")])
    rows = _ingest().active_site_rows([c])
    assert rows[0][3] == "Peptidase Q1"
    assert rows[1][3] == "(no annotated domain)"


def test_nested_domains_are_all_named_not_just_the_first():
    """SARS-CoV-2 RdRp's S5151 is in "Nsp12 RNA-dependent RNA polymerase"
    (4757-5324) AND "RdRp catalytic" (5004-5166). Printing one would hide
    which narrowing is available."""
    c = _rc("SARS-CoV-2 RdRp", "polymerase", 4393, 5324,
            domains=[(4399, 4653, "NiRAN"),
                     (4757, 5324, "Nsp12 RNA-dependent RNA polymerase"),
                     (5004, 5166, "RdRp catalytic")],
            sites=[(5151, "S", "")])
    assert _ingest().active_site_rows([c])[0][3] == (
        "Nsp12 RNA-dependent RNA polymerase + RdRp catalytic")


def test_sites_outside_the_role_filter_are_still_excluded():
    """The helicase is ingested, not scored. Widening the filter to include
    domain-less chains must not also have widened it to ineligible roles."""
    hel = _rc("SARS-CoV-2 helicase", "helicase", 5325, 5925, domains=[],
              sites=[(5400, "C", "")])
    assert _ingest().active_site_rows([hel]) == []


# ------------------------------- what the shared classes actually distinguish
def _pair(mode, shared, allowed=True):
    spec = ChainSpec("x", "P", ("p",), "F", "protease", "T:1")
    c = ResolvedChain(spec=spec)
    from kgav.folds import GateResult
    return (c, c, GateResult(allowed, mode, True, list(shared), "strong", ""))


def test_a_class_in_every_allowed_pair_is_counted_as_such():
    """IPR007094 is "RNA-directed RNA polymerase, catalytic domain". Two
    polymerases sharing it says they are both polymerases -- a constant
    across the mode, which cannot rank one comparator above another."""
    res = ([_pair("motif", ["IPR007094", "IPR043502"])] * 10
           + [_pair("motif", ["IPR001205", "IPR007094"])] * 4)
    n, counts = _ingest().shared_class_census(res)["motif"]
    assert n == 14
    assert counts["IPR007094"] == 14          # universal
    assert counts["IPR043502"] == 10          # discriminating
    assert counts["IPR001205"] == 4


def test_refused_pairs_are_not_counted():
    """The census describes the evidence behind transfers that WERE made.
    A refusal's shared fold is the gate's story, not the transfer's."""
    res = [_pair("nucleophile", ["IPR009003"]),
           _pair("nucleophile", ["IPR009003"], allowed=False)]
    n, counts = _ingest().shared_class_census(res)["nucleophile"]
    assert (n, counts["IPR009003"]) == (1, 1)


def test_modes_are_counted_separately():
    """A class universal in one mode need not be in another, and merging
    them would hide both facts."""
    res = ([_pair("nucleophile", ["IPR009003"])] * 3
           + [_pair("motif", ["IPR007094"])] * 5)
    c = _ingest().shared_class_census(res)
    assert c["nucleophile"][0] == 3 and c["motif"][0] == 5
    assert set(c) == {"nucleophile", "motif"}


def test_a_mode_whose_every_class_is_universal_is_detectable():
    """The condition the report calls out: identical fold evidence on every
    pair, so the hop is an eligibility filter and not a similarity score --
    which is what M8's DWPC would be weighting."""
    res = [_pair("nucleophile", ["G3DSA:2.40.10.10", "IPR009003"])] * 20
    n, counts = _ingest().shared_class_census(res)["nucleophile"]
    assert all(k == n for k in counts.values())
    assert len(counts) == 2
