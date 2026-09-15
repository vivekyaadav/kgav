"""Gate for phase 2: the two-mode comparison, and the rule that a shared fold
does not license transfer across a different nucleophile.
"""

from kgav.folds import NUCLEOPHILES, ChainSpec, ResolvedChain, gate, motif_c_present


def _chain(label, family, role, cat_type, nucleophile, classes):
    spec = ChainSpec(label=label, accession="X", patterns=("x",), family=family,
                     role=role, virus_taxon="NCBITaxon:1")
    return ResolvedChain(spec=spec, chain_name=label, start=1, end=100,
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
