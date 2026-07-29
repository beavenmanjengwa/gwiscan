"""Tests for the ProtParam property computation."""

from gwiscan.features.protparam import COLUMNS, properties


def test_columns():
    assert COLUMNS == ["protein_id", "length_aa", "molecular_weight",
                       "isoelectric_point", "negatively_charged_residues",
                       "positively_charged_residues", "instability_index",
                       "gravy", "aliphatic_index",
                       "ec_cystines", "ec_reduced"]


def test_charged_residue_counts():
    p = properties("DDEERRKK")
    assert p["negatively_charged_residues"] == 4   # 2 Asp + 2 Glu
    assert p["positively_charged_residues"] == 4   # 2 Arg + 2 Lys


def test_extinction_coefficient():
    # WWYYCC: W=2, Y=2, C=2. Expasy/Biopython 280 nm coefficients:
    #   reduced  = 2*5500 + 2*1490            = 13980
    #   cystines = reduced + (2 // 2) * 125   = 14105
    p = properties("WWYYCC")
    assert p["ec_reduced"] == 13980
    assert p["ec_cystines"] == 14105


def test_x_sequence_is_all_dashes():
    p = properties("ACDEFGHXKL")
    assert set(p.values()) == {"-"}
    assert set(p) == set(COLUMNS[1:])


def test_aliphatic_index_scale():
    # amino_acids_percent is mole-percent (0-100): poly-Ala -> AI == 100.
    assert properties("A" * 40)["aliphatic_index"] == 100.0
    # One each of A,V,I,L (+ others) -> A%=V%=I%=L% for the aliphatic residues.
    # ACDEFGHIKL: A=10, V=0, I=10, L=10 -> 10 + 0 + 3.9*(10+10) = 88.0
    assert properties("ACDEFGHIKL")["aliphatic_index"] == 88.0


def test_basic_properties_are_numeric():
    p = properties("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ")
    assert p["length_aa"] == 33
    assert isinstance(p["molecular_weight"], float) and p["molecular_weight"] > 0
    assert isinstance(p["isoelectric_point"], float)
    assert isinstance(p["gravy"], float)
