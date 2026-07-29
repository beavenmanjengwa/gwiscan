"""Tests for the IQ-TREE per-family command builder."""

from gwiscan.config import Config
from gwiscan.iqtree import _iqtree_cmd


def test_cmd_has_core_options(tmp_path):
    cfg = Config(root=tmp_path, IQTREE_BIN="iqtree", IQTREE_MODEL="MFP",
                 IQTREE_BOOTSTRAP=1000, THREADS=8)
    cmd = _iqtree_cmd(cfg, tmp_path / "GNA_aligned.fasta", tmp_path / "trees" / "GNA", n_seqs=20)
    assert cmd[0] == "iqtree"
    assert cmd[cmd.index("-s") + 1].endswith("GNA_aligned.fasta")
    assert cmd[cmd.index("-T") + 1] == "8"
    assert cmd[cmd.index("-m") + 1] == "MFP"
    assert cmd[cmd.index("--prefix") + 1].endswith("GNA")
    assert "-redo" in cmd
    assert cmd[cmd.index("-seed") + 1] == "12345"  # fixed seed -> reproducible
    assert cmd[cmd.index("-B") + 1] == "1000"      # bootstrap on with enough taxa


def test_bootstrap_skipped_when_too_few_taxa(tmp_path):
    cfg = Config(root=tmp_path, IQTREE_BOOTSTRAP=1000)
    cmd = _iqtree_cmd(cfg, tmp_path / "X_aligned.fasta", tmp_path / "X", n_seqs=3)
    assert "-B" not in cmd                          # UFBoot needs >=4 sequences


def test_bootstrap_off_when_zero(tmp_path):
    cfg = Config(root=tmp_path, IQTREE_BOOTSTRAP=0)
    cmd = _iqtree_cmd(cfg, tmp_path / "X_aligned.fasta", tmp_path / "X", n_seqs=50)
    assert "-B" not in cmd
