#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# test_deeptmhmm.py - DeepTMHMM parsing tests.                                                     #
#                                                                                                  #
# Coordinates come from TMRs.gff3 (the structured coordinate output); the per-protein class label  #
# comes from the predicted_topologies.3line header.                                                #
#                                                                                                  #
####################################################################################################
"""

import pytest

from gwiscan import io
from gwiscan.config import Config
from gwiscan.features import deeptmhmm
from gwiscan.features.deeptmhmm import (
    OUT_HEADER,
    build_rows,
    parse_gff,
    parse_type_labels,
)

# Real TMRs.gff3 shape (tab-separated region rows, '#'/ '//' separators). Includes
# trailing tabs like the actual output, plus a beta-barrel record.
GFF = (
    "##gff-version 3\n"
    "# Medtr0002s0020.1 Length: 288\n"
    "# Medtr0002s0020.1 Number of predicted TMRs: 0\n"
    "Medtr0002s0020.1\tinside\t1\t288\t\t\t\n"
    "//\n"
    "# Medtr0015s0030.1 Length: 662\n"
    "# Medtr0015s0030.1 Number of predicted TMRs: 1\n"
    "Medtr0015s0030.1\tsignal\t1\t20\t\t\t\n"
    "Medtr0015s0030.1\toutside\t21\t264\t\t\t\n"
    "Medtr0015s0030.1\tTMhelix\t265\t285\t\t\t\n"
    "Medtr0015s0030.1\tinside\t286\t662\t\t\t\n"
    "//\n"
    "# Barrel.1 Length: 50\n"
    "# Barrel.1 Number of predicted TMRs: 2\n"
    "Barrel.1\tBeta sheet\t10\t20\t\t\t\n"
    "Barrel.1\tperiplasm\t21\t30\t\t\t\n"
    "Barrel.1\tBeta sheet\t31\t40\t\t\t\n"
    "//\n"
)

THREE_LINE_HEADERS = (
    ">Medtr0002s0020.1 | GLOB\nAAA\nIII\n"
    ">Medtr0015s0030.1 | SP+TM\nAAA\nIII\n"
    ">Barrel.1 | BETA\nAAA\nIII\n"
)


def test_parse_gff_coordinates(tmp_path):
    f = tmp_path / "TMRs.gff3"
    f.write_text(GFF)
    coords = parse_gff(f)

    glob = coords["Medtr0002s0020.1"]
    assert glob["length"] == 288
    assert glob["signal"] == [] and glob["tm"] == [] and glob["beta"] == []

    sptm = coords["Medtr0015s0030.1"]
    assert sptm["length"] == 662
    assert sptm["signal"] == [(1, 20)]
    assert sptm["tm"] == [(265, 285)]

    barrel = coords["Barrel.1"]
    assert barrel["beta"] == [(10, 20), (31, 40)]
    assert barrel["tm"] == []


def test_type_labels(tmp_path):
    f = tmp_path / "predicted_topologies.3line"
    f.write_text(THREE_LINE_HEADERS)
    labels = parse_type_labels(f)
    assert labels == {
        "Medtr0002s0020.1": "GLOB",
        "Medtr0015s0030.1": "SP+TM",
        "Barrel.1": "BETA",
    }


def test_build_rows_combines_both(tmp_path):
    gff = tmp_path / "TMRs.gff3"
    three = tmp_path / "predicted_topologies.3line"
    gff.write_text(GFF)
    three.write_text(THREE_LINE_HEADERS)

    rows = {r[0]: dict(zip(OUT_HEADER, r)) for r in build_rows(gff, three)}

    sptm = rows["Medtr0015s0030.1"]
    assert sptm["topology"] == "SP+TM"            # from .3line
    assert sptm["signal_peptide"] == "1-20"       # from gff, matches TMRs.gff3
    assert sptm["tm_regions"] == "265-285"
    assert sptm["n_tm_regions"] == 1
    assert sptm["beta_regions"] == "-"

    glob = rows["Medtr0002s0020.1"]
    assert glob["topology"] == "GLOB"
    assert glob["signal_peptide"] == "-"
    assert glob["tm_regions"] == "-"
    assert glob["n_tm_regions"] == 0

    # A beta-barrel gets its topology label AND its β-strand ranges (general
    # families can be β-barrels -- transporters, porins).
    barrel = rows["Barrel.1"]
    assert barrel["topology"] == "BETA"
    assert barrel["tm_regions"] == "-"
    assert barrel["n_tm_regions"] == 0
    assert barrel["beta_regions"] == "10-20;31-40"


def test_build_rows_without_3line(tmp_path):
    gff = tmp_path / "TMRs.gff3"
    gff.write_text(GFF)
    rows = {r[0]: dict(zip(OUT_HEADER, r)) for r in build_rows(gff, None)}
    # Coordinates still present; topology label just blank without the .3line.
    assert rows["Medtr0015s0030.1"]["tm_regions"] == "265-285"
    assert rows["Medtr0015s0030.1"]["topology"] == ""


def test_local_run_ignores_stale_biolib_output(tmp_path, monkeypatch):
    cfg = Config(root=tmp_path, DEEPTMHMM_MODE="local", DEEPTMHMM_DIR=str(tmp_path / "tool"))
    (tmp_path / "tool").mkdir()
    (tmp_path / "tool" / "predict.py").write_text("# stub\n")
    cfg.result("final_candidates.fasta").parent.mkdir(parents=True, exist_ok=True)
    cfg.result("final_candidates.fasta").write_text(">new\nAAAA\n")
    stale = cfg.result("deeptmhmm_out") / "biolib_results"
    stale.mkdir(parents=True)
    (stale / "TMRs.gff3").write_text("# stale Length: 1\nstale\tTMhelix\t1\t1\n")

    def fake_run(cmd, **kwargs):
        out = cfg.result("deeptmhmm_out") / "predict"
        out.mkdir(parents=True)
        (out / "TMRs.gff3").write_text("# fresh Length: 1\nfresh\tTMhelix\t1\t1\n")

    monkeypatch.setattr(deeptmhmm.external, "run", fake_run)
    deeptmhmm.run(cfg)
    rows = io.read_tsv(cfg.result("deeptmhmm.tsv"))
    assert rows["protein_id"].tolist() == ["fresh"]


def test_unknown_deeptmhmm_mode_fails_explicitly(tmp_path):
    cfg = Config(root=tmp_path, DEEPTMHMM_MODE="loacl")
    cfg.result("final_candidates.fasta").parent.mkdir(parents=True, exist_ok=True)
    cfg.result("final_candidates.fasta").write_text(">p\nAAAA\n")
    with pytest.raises(ValueError, match="must be 'local' or 'biolib'"):
        deeptmhmm.run(cfg)
