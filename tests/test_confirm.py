"""Tests for the InterProScan confirmation gate."""

import pandas as pd

from gwiscan.confirm import confirmed_rows

FAMILY_TO_PFAM = {"Legume": "PF00139", "GNA": "PF01453", "CRA": "PF00704"}
# EUL intentionally absent (no Pfam).

# candidates_merged: protein_id + assigned family
MERGED = pd.DataFrame([
    {"protein_id": "P_leg", "family": "Legume"},   # has PF00139 -> confirmed
    {"protein_id": "P_gna", "family": "GNA"},      # lacks PF01453 -> dropped
    {"protein_id": "P_cra", "family": "CRA"},      # has PF00704 -> confirmed
    {"protein_id": "P_eul", "family": "EUL"},      # no Pfam family -> kept
])

# InterProScan Pfam matches (protein_id, sig_acc)
INTERPRO = pd.DataFrame([
    {"protein_id": "P_leg", "sig_acc": "PF00139"},
    {"protein_id": "P_leg", "sig_acc": "PF00069"},   # kinase, irrelevant
    {"protein_id": "P_gna", "sig_acc": "PF00069"},   # NOT its family's PF01453
    {"protein_id": "P_cra", "sig_acc": "PF00704"},
])


def test_gate_keeps_confirmed_and_eul_drops_unconfirmed():
    final = confirmed_rows(MERGED, INTERPRO, FAMILY_TO_PFAM)
    kept = dict(zip(final["protein_id"], final["family"]))

    assert kept == {"P_leg": "Legume", "P_cra": "CRA", "P_eul": "EUL"}
    assert "P_gna" not in kept          # GNA hit with no PF01453 -> false positive dropped


def test_eul_kept_even_with_no_interpro_hit():
    # EUL protein has no InterProScan rows at all -> still kept (no Pfam to gate).
    final = confirmed_rows(MERGED, INTERPRO, FAMILY_TO_PFAM)
    assert "P_eul" in set(final["protein_id"])
