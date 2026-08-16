#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# test_clipkit.py - ClipKIT trimming-stage tests, and IQ-TREE's use of its output.                 #
#                                                                                                  #
####################################################################################################
"""

import pytest

from gwiscan import clipkit, iqtree
from gwiscan.config import Config


# --- trimmed_path ------------------------------------------------------------

def test_trimmed_path_mirrors_aligned_name(tmp_path):
    aln = tmp_path / "msa" / "GNA_aligned.fasta"
    assert clipkit.trimmed_path(aln) == tmp_path / "msa" / "GNA_trimmed.fasta"


def test_trimmed_path_handles_mature_track(tmp_path):
    aln = tmp_path / "msa" / "GNA_mature_aligned.fasta"
    assert clipkit.trimmed_path(aln).name == "GNA_mature_trimmed.fasta"


# --- run(): builds the right ClipKIT command per alignment -------------------

def test_run_trims_each_aligned_fasta(tmp_path, monkeypatch):
    msa = tmp_path / "intermediate" / "msa"
    msa.mkdir(parents=True)
    (msa / "GNA_aligned.fasta").write_text(">a\nAAAA\n>b\nAAAA\n")
    (msa / "EUL_aligned.fasta").write_text(">a\nCCCC\n>b\nCCCC\n")
    # a stray trimmed file must NOT be re-trimmed (glob is *_aligned.fasta only)
    (msa / "GNA_trimmed.fasta").write_text(">a\nAA\n")

    cmds = []
    monkeypatch.setattr(clipkit.external, "require", lambda b: None)
    monkeypatch.setattr(clipkit.external, "run", lambda cmd, **kw: cmds.append([str(c) for c in cmd]))

    clipkit.run(Config(root=tmp_path))

    assert len(cmds) == 2                              # GNA + EUL, not the stray trimmed
    for cmd in cmds:
        assert cmd[0] == "clipkit"
        assert "-m" in cmd and cmd[cmd.index("-m") + 1] == "smart-gap"
        assert "-o" in cmd
        out = cmd[cmd.index("-o") + 1]
        assert out.endswith("_trimmed.fasta")


def test_run_honours_clipkit_mode(tmp_path, monkeypatch):
    msa = tmp_path / "intermediate" / "msa"
    msa.mkdir(parents=True)
    (msa / "GNA_aligned.fasta").write_text(">a\nAAAA\n>b\nAAAA\n")

    cmds = []
    monkeypatch.setattr(clipkit.external, "require", lambda b: None)
    monkeypatch.setattr(clipkit.external, "run", lambda cmd, **kw: cmds.append([str(c) for c in cmd]))

    clipkit.run(Config(root=tmp_path, CLIPKIT_MODE="kpic-smart-gap"))

    assert cmds[0][cmds[0].index("-m") + 1] == "kpic-smart-gap"


def test_run_warns_when_no_alignments(tmp_path, monkeypatch, capsys):
    (tmp_path / "intermediate" / "msa").mkdir(parents=True)
    monkeypatch.setattr(clipkit.external, "require", lambda b: None)
    monkeypatch.setattr(clipkit.external, "run", lambda *a, **k: pytest.fail("should not run"))
    clipkit.run(Config(root=tmp_path))
    assert "No" in capsys.readouterr().out


# --- iqtree consumes trimmed if present, else the raw alignment --------------

def test_iqtree_prefers_trimmed_alignment(tmp_path):
    msa = tmp_path / "msa"
    msa.mkdir()
    aln = msa / "GNA_aligned.fasta"
    aln.write_text(">a\nAAAA\n")
    trimmed = msa / "GNA_trimmed.fasta"
    trimmed.write_text(">a\nAA\n")
    assert iqtree._tree_input(aln) == trimmed


def test_iqtree_falls_back_to_untrimmed_when_no_trimmed(tmp_path):
    msa = tmp_path / "msa"
    msa.mkdir()
    aln = msa / "GNA_aligned.fasta"
    aln.write_text(">a\nAAAA\n")
    # no trimmed sibling -> use the raw alignment (ClipKIT optional / skipped)
    assert iqtree._tree_input(aln) == aln
