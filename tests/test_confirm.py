#!/usr/bin/env python3
"""
#####################################################################################################
#                                                                                                   #
# test_confirm.py - InterProScan confirmation-gate tests.                                           #
#                                                                                                   #
# The gate keeps a (protein, family) row when the protein carries ANY of the family's confirmation  #
# accessions, across whichever member database (Pfam / CDD / PANTHER / ...) each prefix names. A    #
# family with no accessions (a custom-HMM family) is not gated and is kept.                         #
#                                                                                                   #
####################################################################################################
"""

import pandas as pd

from gwiscan.confirm import confirmed_rows

# family -> set of confirmation accessions (what io.family_confirm_accessions returns).
FAMILY_TO_ACCS = {
    "Legume": {"PF00139"},                 # Pfam
    "GNA": {"PF01453"},                    # Pfam
    "CRA": {"cd02879"},                    # no Pfam -> confirms on a CDD accession
    "EUL": {"PTHR31257"},                  # no Pfam -> confirms on a PANTHER accession
    "MULTI": {"PF01476", "PTHR27007"},     # confirms on either
    "CUSTOM": set(),                       # custom-HMM family: not gated -> always kept
}

MERGED = pd.DataFrame([
    {"protein_id": "P_leg", "family": "Legume"},
    {"protein_id": "P_gna", "family": "GNA"},
    {"protein_id": "P_cra", "family": "CRA"},
    {"protein_id": "P_eul", "family": "EUL"},
    {"protein_id": "P_multi", "family": "MULTI"},
    {"protein_id": "P_cust", "family": "CUSTOM"},
])

# InterProScan matches (protein_id, sig_acc), any member database.
INTERPRO = pd.DataFrame([
    {"protein_id": "P_leg", "sig_acc": "PF00139"},
    {"protein_id": "P_leg", "sig_acc": "PF00069"},      # kinase, irrelevant
    {"protein_id": "P_gna", "sig_acc": "PF00069"},      # NOT its family's PF01453
    {"protein_id": "P_cra", "sig_acc": "cd02879"},      # CDD confirms CRA
    {"protein_id": "P_eul", "sig_acc": "PTHR31257"},    # PANTHER confirms EUL
    {"protein_id": "P_multi", "sig_acc": "PTHR27007"},  # one of MULTI's two
])


def test_gate_confirms_on_any_application():
    final = confirmed_rows(MERGED, INTERPRO, FAMILY_TO_ACCS)
    kept = dict(zip(final["protein_id"], final["family"]))
    assert kept == {
        "P_leg": "Legume", "P_cra": "CRA", "P_eul": "EUL",
        "P_multi": "MULTI", "P_cust": "CUSTOM",
    }
    assert "P_gna" not in kept          # GNA hit without PF01453 -> false positive dropped


def test_multi_accession_confirms_on_either():
    # MULTI needs PF01476 OR PTHR27007; the protein carries only PTHR27007.
    final = confirmed_rows(MERGED, INTERPRO, FAMILY_TO_ACCS)
    assert "P_multi" in set(final["protein_id"])


def test_non_pfam_family_gated_not_auto_passed():
    # CRA/EUL have no Pfam but DO get a real gate now: a CRA protein missing cd02879
    # is dropped rather than auto-confirmed.
    interpro_no_cra = INTERPRO[INTERPRO["protein_id"] != "P_cra"]
    final = confirmed_rows(MERGED, interpro_no_cra, FAMILY_TO_ACCS)
    assert "P_cra" not in set(final["protein_id"])


def test_custom_family_kept_without_any_hit():
    # CUSTOM has no confirmation accessions (and no InterProScan rows) -> kept.
    final = confirmed_rows(MERGED, INTERPRO, FAMILY_TO_ACCS)
    assert "P_cust" in set(final["protein_id"])
