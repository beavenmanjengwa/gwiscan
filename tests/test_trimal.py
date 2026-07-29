"""Tests for the trimAl trimming stage and IQ-TREE's use of its output."""

import pytest

from gwiscan import iqtree, trimal
from gwiscan.config import Config


# --- trimmed_path / method flag ----------------------------------------------

def test_trimmed_path_mirrors_aligned_name(tmp_path):
    aln = tmp_path / "msa" / "GNA_aligned.fasta"
    assert trimal.trimmed_path(aln) == tmp_path / "msa" / "GNA_trimmed.fasta"


def test_trimmed_path_handles_mature_track(tmp_path):
    aln = tmp_path / "msa" / "GNA_mature_aligned.fasta"
    assert trimal.trimmed_path(aln).name == "GNA_mature_trimmed.fasta"


def test_method_flag_default_is_automated1(tmp_path):
    assert trimal._method_flag(Config(root=tmp_path)) == "-automated1"


def test_method_flag_prefixes_dash(tmp_path):
    assert trimal._method_flag(Config(root=tmp_path, TRIMAL_METHOD="gappyout")) == "-gappyout"


def test_method_flag_passes_through_leading_dash(tmp_path):
    assert trimal._method_flag(Config(root=tmp_path, TRIMAL_METHOD="-strict")) == "-strict"


# --- run(): builds the right trimAl command per alignment --------------------

def test_run_trims_each_aligned_fasta(tmp_path, monkeypatch):
    msa = tmp_path / "intermediate" / "msa"
    msa.mkdir(parents=True)
    (msa / "GNA_aligned.fasta").write_text(">a\nAAAA\n>b\nAAAA\n")
    (msa / "EUL_aligned.fasta").write_text(">a\nCCCC\n>b\nCCCC\n")
    # a stray trimmed file must NOT be re-trimmed (glob is *_aligned.fasta only)
    (msa / "GNA_trimmed.fasta").write_text(">a\nAA\n")

    cmds = []
    monkeypatch.setattr(trimal.external, "require", lambda b: None)
    monkeypatch.setattr(trimal.external, "run", lambda cmd, **kw: cmds.append([str(c) for c in cmd]))

    trimal.run(Config(root=tmp_path))

    assert len(cmds) == 2                              # GNA + EUL, not the stray trimmed
    for cmd in cmds:
        assert cmd[0] == "trimal"
        assert "-automated1" in cmd
        assert "-in" in cmd and "-out" in cmd
        out = cmd[cmd.index("-out") + 1]
        assert out.endswith("_trimmed.fasta")


def test_run_warns_when_no_alignments(tmp_path, monkeypatch, capsys):
    (tmp_path / "intermediate" / "msa").mkdir(parents=True)
    monkeypatch.setattr(trimal.external, "require", lambda b: None)
    monkeypatch.setattr(trimal.external, "run", lambda *a, **k: pytest.fail("should not run"))
    trimal.run(Config(root=tmp_path))
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
    # no trimmed sibling -> use the raw alignment (trimAl optional / skipped)
    assert iqtree._tree_input(aln) == aln
