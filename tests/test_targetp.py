#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# test_targetp.py - TargetP 2.0 short-summary parsing tests.                                       #
#                                                                                                  #
# The full summary is kept for reproducibility (type + per-class probabilities + cleavage site).   #
# Parsing reads the summary's own header, so it stays correct for the plant layout (-org pl, extra #
# cTP/luTP columns) and keeps the space-containing CS field intact.                                #
#                                                                                                  #
####################################################################################################
"""

from gwiscan.features.targetp import parse_summary, resolve_targetp

# Plant short summary: header '#'-lines, then ID, Prediction, noTP, SP, mTP,
# cTP, luTP, CS Position. Tab-separated; the CS field contains spaces.
PLANT_SUMMARY = (
    "# TargetP-2.0\tOrganism: Plant\tTimestamp: 20260710\n"
    "# ID\tPrediction\tnoTP\tSP\tmTP\tcTP\tluTP\tCS Position\n"
    "Medtr_glob.1\tnoTP\t0.9987\t0.0001\t0.0007\t0.0003\t0.0002\t\n"
    "Medtr_secr.1\tSP\t0.0012\t0.9980\t0.0005\t0.0002\t0.0001\tCS pos: 24-25. Pr: 0.89\n"
    "Medtr_chlo.1\tcTP\t0.0100\t0.0100\t0.0100\t0.9600\t0.0100\tCS pos: 55-56. Pr: 0.70\n"
)

NONPLANT_SUMMARY = (
    "# TargetP-2.0\tOrganism: Non-plant\n"
    "# ID\tPrediction\tnoTP\tSP\tmTP\tCS Position\n"
    "Prot_a\tSP\t0.02\t0.97\t0.01\tCS pos: 20-21. Pr: 0.8\n"
)


def test_header_renamed_and_full_columns(tmp_path):
    f = tmp_path / "targetp_summary.targetp2"
    f.write_text(PLANT_SUMMARY)
    header, rows = parse_summary(f)

    # ID -> protein_id, Prediction -> targetp_type, CS Position -> cs_position;
    # the plant class columns are kept verbatim.
    assert header == ["protein_id", "targetp_type", "noTP", "SP", "mTP",
                      "cTP", "luTP", "cs_position"]
    assert len(rows) == 3
    assert all(len(r) == len(header) for r in rows)


def test_values_preserved(tmp_path):
    f = tmp_path / "targetp_summary.targetp2"
    f.write_text(PLANT_SUMMARY)
    header, rows = parse_summary(f)
    by_id = {r[0]: dict(zip(header, r)) for r in rows}

    secr = by_id["Medtr_secr.1"]
    assert secr["targetp_type"] == "SP"
    assert secr["SP"] == "0.9980"                       # probability retained
    assert secr["cs_position"] == "CS pos: 24-25. Pr: 0.89"  # spaces intact

    glob = by_id["Medtr_glob.1"]
    assert glob["targetp_type"] == "noTP"
    assert glob["cs_position"] == ""                    # empty CS for noTP

    chlo = by_id["Medtr_chlo.1"]
    assert chlo["targetp_type"] == "cTP"                # plant-only class
    assert chlo["cTP"] == "0.9600"


def test_resolve_targetp_install_dir(tmp_path):
    # Given the unpacked install directory, resolve to <dir>/bin/targetp.
    (tmp_path / "bin").mkdir()
    exe = tmp_path / "bin" / "targetp"
    exe.write_text("#!/bin/sh\n")
    assert resolve_targetp(str(tmp_path)) == str(exe)
    # bin/targetp under any given path is preferred.
    assert resolve_targetp(str(tmp_path / "bin" / "targetp")) == str(exe)


def test_resolve_targetp_plain_name():
    # A bare PATH name with no bin/targetp beside it is returned unchanged.
    assert resolve_targetp("targetp") == "targetp"


def test_nonplant_layout(tmp_path):
    f = tmp_path / "targetp_summary.targetp2"
    f.write_text(NONPLANT_SUMMARY)
    header, rows = parse_summary(f)
    # No cTP/luTP columns for non-plant; parsing follows the header either way.
    assert header == ["protein_id", "targetp_type", "noTP", "SP", "mTP", "cs_position"]
    assert dict(zip(header, rows[0]))["SP"] == "0.97"
