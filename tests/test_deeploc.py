#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# test_deeploc.py - DeepLoc 2.1 results-parsing tests.                                             #
#                                                                                                  #
# Checks that only protein_id and the summary columns are kept (the per-compartment probability    #
# columns stay in the raw file), the ID column is renamed so it joins the pipeline, and DeepLoc's  #
# Lysosome/Vacuole class collapses to Vacuole in either slash order.                               #
#                                                                                                  #
####################################################################################################
"""

from gwiscan.features.deeploc import parse_results

# DeepLoc-style results CSV: summary columns then per-compartment probabilities
# (including a 'Lysosome/Vacuole' probability column, which must be dropped).
RESULTS_CSV = (
    "Protein_ID,Localizations,Signals,Membrane types,Cytoplasm,Lysosome/Vacuole,Plastid\n"
    "XP_1,Cytoplasm,Peroxisomal targeting signal,Soluble,0.90,0.01,0.02\n"
    "XP_2,Lysosome/Vacuole,,Peripheral,0.10,0.85,0.00\n"
    "XP_3,Vacuole/Lysosome,,Soluble,0.10,0.80,0.00\n"
    "XP_4,Plastid,,Soluble,0.00,0.00,0.95\n"
)


def _write(tmp_path):
    p = tmp_path / "results.csv"
    p.write_text(RESULTS_CSV)
    return p


def test_keeps_only_summary_columns(tmp_path):
    df = parse_results(_write(tmp_path))
    assert list(df.columns) == ["protein_id", "Localizations", "Signals", "Membrane types"]
    # probability columns (including the 'Lysosome/Vacuole' probability) dropped
    assert "Cytoplasm" not in df.columns
    assert "Lysosome/Vacuole" not in df.columns


def test_lysosome_vacuole_collapsed(tmp_path):
    df = parse_results(_write(tmp_path))
    loc = dict(zip(df["protein_id"], df["Localizations"]))
    assert loc["XP_2"] == "Vacuole"     # Lysosome/Vacuole -> Vacuole
    assert loc["XP_3"] == "Vacuole"     # Vacuole/Lysosome -> Vacuole (either order)
    assert loc["XP_1"] == "Cytoplasm"   # unaffected
    assert loc["XP_4"] == "Plastid"


def test_id_column_renamed(tmp_path):
    df = parse_results(_write(tmp_path))
    assert "protein_id" in df.columns
    assert set(df["protein_id"]) == {"XP_1", "XP_2", "XP_3", "XP_4"}
