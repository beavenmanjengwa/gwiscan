#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# test_domain_bed.py - Lectin domain BED builder tests.                                            #
#                                                                                                  #
# Covers the family-driven logic: one row per domain match, tandem domains as separate rows, CRA   #
# PF00704 pulled only for DIAMOND-identified CRA proteins (not every PF00704 hit), incidental non- #
# lectin domains excluded, and EUL skipped.                                                        #
#                                                                                                  #
####################################################################################################
"""

import pandas as pd

from gwiscan.domain_bed import BED_HEADER, build_bed

FAMILY_TO_PFAM = {"Legume": "PF00139", "GNA": "PF01453", "CRA": "PF00704"}
# EUL intentionally absent (no Pfam).

# candidates_merged: protein_id + assigned family (family).
CANDIDATES = pd.DataFrame([
    {"protein_id": "P_leg", "family": "Legume"},   # one legume domain
    {"protein_id": "P_tandem", "family": "GNA"},   # two GNA domains
    {"protein_id": "P_cra", "family": "CRA"},      # CRA -> PF00704 wanted
    {"protein_id": "P_chit", "family": "Legume"},  # legume protein that ALSO has PF00704 (a chitinase domain)
    {"protein_id": "P_eul", "family": "EUL"},      # no Pfam -> skipped
])

# InterProScan Pfam matches (protein_id, sig_acc, start, end).
INTERPRO = pd.DataFrame([
    {"protein_id": "P_leg", "sig_acc": "PF00139", "start": 34, "end": 280},
    {"protein_id": "P_leg", "sig_acc": "PF00069", "start": 300, "end": 400},   # kinase, not lectin
    {"protein_id": "P_tandem", "sig_acc": "PF01453", "start": 10, "end": 110},
    {"protein_id": "P_tandem", "sig_acc": "PF01453", "start": 130, "end": 230},
    {"protein_id": "P_cra", "sig_acc": "PF00704", "start": 25, "end": 300},
    {"protein_id": "P_chit", "sig_acc": "PF00139", "start": 40, "end": 285},
    {"protein_id": "P_chit", "sig_acc": "PF00704", "start": 320, "end": 590},   # chitinase, NOT a CRA lectin
    {"protein_id": "P_eul", "sig_acc": "PF00139", "start": 5, "end": 200},      # stray hit, but EUL has no pfam
])


def _bed():
    return build_bed(CANDIDATES, FAMILY_TO_PFAM, INTERPRO)


def test_header():
    assert BED_HEADER == ["protein_id", "pfam_id", "start", "end"]


def test_one_row_per_domain_and_tandem():
    rows = _bed()
    tandem = [r for r in rows if r[0] == "P_tandem"]
    assert tandem == [["P_tandem", "PF01453", 10, 110], ["P_tandem", "PF01453", 130, 230]]


def test_incidental_nonlectin_excluded():
    rows = _bed()
    leg = [r for r in rows if r[0] == "P_leg"]
    assert leg == [["P_leg", "PF00139", 34, 280]]        # PF00069 kinase dropped


def test_cra_pf00704_only_for_cra_protein():
    rows = _bed()
    # CRA protein: its PF00704 domain IS included
    assert ["P_cra", "PF00704", 25, 300] in rows
    # Legume protein that also has a PF00704 chitinase domain: PF00704 NOT included,
    # only its assigned-family (Legume/PF00139) domain
    chit = [r for r in rows if r[0] == "P_chit"]
    assert chit == [["P_chit", "PF00139", 40, 285]]
    assert not any(r[1] == "PF00704" and r[0] == "P_chit" for r in rows)


def test_eul_skipped():
    rows = _bed()
    assert not any(r[0] == "P_eul" for r in rows)


# --- custom-HMM family: coordinates from its own hmmscan hit --------------------

CUSTOM_CANDIDATES = pd.DataFrame([
    {"protein_id": "P_cra1", "family": "CRA", "method": "hmm",
     "start": 25, "end": 300},
    {"protein_id": "P_cra1", "family": "CRA", "method": "hmm",   # tandem
     "start": 340, "end": 610},
    {"protein_id": "P_cra2", "family": "CRA", "method": "blast",  # blast row: no coords
     "start": "-", "end": "-"},
    {"protein_id": "P_gna", "family": "GNA", "method": "hmm",
     "start": 10, "end": 110},
])


def test_custom_hmm_uses_hmmscan_coords_not_interpro():
    # CRA is a custom-HMM family; GNA stays a Pfam family.
    rows = build_bed(CUSTOM_CANDIDATES, {"GNA": "PF01453"},
                     INTERPRO, custom_families={"CRA"})
    cra = [r for r in rows if r[0].startswith("P_cra")]
    # hmmscan coords used; pfam_id is the family name; tandem kept; blast row skipped.
    assert cra == [
        ["P_cra1", "CRA", 25, 300],
        ["P_cra1", "CRA", 340, 610],
    ]
    # The blast-only CRA row contributes no coordinates.
    assert not any(r[0] == "P_cra2" for r in rows)
