#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# test_family_map.py - family_map schema-loader tests (io.family_records).                         #
#                                                                                                  #
# Locks the one correctness-critical rule: a family may carry a Pfam accession for                 #
# confirmation/coordinates while not being pressed as an identifying HMM. CRA's PF00704 also       #
# matches active chitinases, so pressing it would mis-identify every GH18 protein as CRA;          #
# HmmPress=no prevents that.                                                                       #
#                                                                                                  #
####################################################################################################
"""

from pathlib import Path

from gwiscan import io
from gwiscan.domain_bed import family_pfam_map

SYNThetic = (
    "Family\tPfamModel\tBlastModel\tHmmPress\n"
    "GNA\tPF01453\tgna.fasta\n"          # normal: pressed
    "CRA\tPF00704\tcra.fasta\tno\n"      # confirm/coords only, NOT pressed
    "EUL\t-\teul.fasta\n"                # no Pfam: blast-only
    "Legume\tPF00139.27\tleg.fasta\n"    # versioned accession -> bare
)


def _write(tmp_path, text):
    p = tmp_path / "family.tsv"
    p.write_text(text)
    return p


def test_records_fields(tmp_path):
    recs = {r["family"]: r for r in io.family_records(_write(tmp_path, SYNThetic))}

    assert recs["GNA"]["pfam_model"] == "PF01453"
    assert recs["GNA"]["hmm_press"] is True

    # CRA: has the Pfam, but is explicitly not pressed.
    assert recs["CRA"]["pfam_model"] == "PF00704"
    assert recs["CRA"]["hmm_press"] is False

    # EUL: no Pfam at all -> nothing to press.
    assert recs["EUL"]["pfam_model"] == ""
    assert recs["EUL"]["hmm_press"] is False

    # Version suffix is stripped to the bare accession.
    assert recs["Legume"]["pfam_model"] == "PF00139"


def test_pfam_to_family_includes_unpressed_cra(tmp_path):
    # CRA is present for confirmation/coordinate mapping even though unpressed.
    m = io.pfam_to_family(_write(tmp_path, SYNThetic))
    assert m == {"PF01453": "GNA", "PF00704": "CRA", "PF00139": "Legume"}


def test_family_pfam_map_inverse(tmp_path):
    m = family_pfam_map(_write(tmp_path, SYNThetic))
    assert m["CRA"] == "PF00704"
    assert m["GNA"] == "PF01453"
    assert "EUL" not in m


SUPER_MAP = (
    "Superfamily\tFamily\tPfamModel\tBlastModel\n"
    "Lectin\tGNA\tPF01453\tgna.fasta\n"
    "Lectin\tLegume\tPF00139\tleg.fasta\n"
    "Lectin\tCRA\tPF00704\tcra.fasta\n"
)


def test_family_to_superfamily(tmp_path):
    m = io.family_to_superfamily(_write(tmp_path, SUPER_MAP))
    assert m == {"GNA": "Lectin", "Legume": "Lectin", "CRA": "Lectin"}


def test_family_mode_has_no_superfamily(tmp_path):
    # Flat family_map (no Superfamily column) -> empty mapping (family mode).
    assert io.family_to_superfamily(_write(tmp_path, SYNThetic)) == {}


CUSTOM_MAP = (
    "Family\tPfamModel\tBlastModel\tHmmPress\n"
    "GNA\tPF01453\tgna.fasta\n"          # Pfam accession -> download + confirm key
    "CRA\tCRA.hmm\tcra.fasta\n"          # custom HMM -> self-contained, no Pfam
    "EUL\t-\teul.fasta\n"                # DIAMOND-only
)


def test_custom_hmm_family_has_no_pfam_but_is_pressed(tmp_path):
    recs = {r["family"]: r for r in io.family_records(_write(tmp_path, CUSTOM_MAP))}

    # Pfam route: accession is the confirm key; pressed file is <acc>.hmm.
    assert recs["GNA"]["pfam_model"] == "PF01453"
    assert recs["GNA"]["hmm_file"] == "PF01453.hmm"
    assert recs["GNA"]["hmm_is_custom"] is False
    assert recs["GNA"]["hmm_press"] is True

    # Custom HMM: identifies itself, no Pfam / InterProScan, pressed as given.
    assert recs["CRA"]["pfam_model"] == ""          # no confirm Pfam
    assert recs["CRA"]["hmm_file"] == "CRA.hmm"
    assert recs["CRA"]["hmm_is_custom"] is True
    assert recs["CRA"]["hmm_press"] is True

    # DIAMOND-only family: no HMM at all.
    assert recs["EUL"]["hmm_file"] == ""
    assert recs["EUL"]["hmm_press"] is False


def test_custom_hmm_absent_from_pfam_maps(tmp_path):
    p = _write(tmp_path, CUSTOM_MAP)
    # A custom-HMM family has no Pfam, so it never appears in the Pfam maps
    # (it is not InterProScan-confirmed; it confirms via its own HMM).
    assert "CRA" not in io.pfam_to_family(p).values()
    assert "PF00704" not in io.pfam_to_family(p)
    assert io.custom_hmm_families(p) == {"CRA"}


def test_shipped_family_map_cra_is_custom_hmm():
    shipped = Path(__file__).resolve().parents[1] / "config" / "family.tsv"
    recs = {r["family"]: r for r in io.family_records(shipped)}
    # CRA uses a custom HMM (hmmscan + BLAST, no InterProScan): no Pfam accession.
    assert recs["CRA"]["hmm_is_custom"] is True
    assert recs["CRA"]["pfam_model"] == ""
    assert recs["CRA"]["hmm_file"] == "CRA.hmm"
    # Every family has a BlastModel (DIAMOND runs for all).
    assert all(r["blast_model"] for r in recs.values())
