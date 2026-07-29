#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# iqtree.py - Per-family maximum-likelihood phylogenetic trees with IQ-TREE (the `iqtree` stage).  #
#                                                                                                  #
# Builds a tree for each family from its MAFFT alignment and writes the IQ-TREE output             #
# (.treefile plus the model/report files) to intermediate/trees/. Model selection uses ModelFinder      #
# (IQTREE_MODEL, default MFP) with ultrafast bootstrap (IQTREE_BOOTSTRAP, default 1000; 0 off).    #
#                                                                                                  #
# Install IQ-TREE via conda so iqtree resolves by name: conda install bioconda::iqtree             #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

from Bio import SeqIO

from . import external, trimal
from .config import Config


def _tree_input(aligned):
    """The alignment IQ-TREE should build from: the trimAl-trimmed file if the
    `trim` stage produced one, otherwise the raw MAFFT alignment (trimAl is
    optional / may have been skipped, so fall back gracefully)."""
    trimmed = trimal.trimmed_path(aligned)
    return trimmed if trimmed.exists() else aligned


def _bootstrap_replicates(cfg: Config) -> int:
    """Ultrafast bootstrap replicate count from config (0 if unset/invalid)."""
    try:
        return int(cfg.IQTREE_BOOTSTRAP or 0)
    except (TypeError, ValueError):
        return 0


def _iqtree_cmd(cfg: Config, aln, prefix, n_seqs: int) -> list:
    """IQ-TREE command for one family alignment. Ultrafast bootstrap is added only
    when enabled and there are enough sequences for it (>=4)."""
    cmd = [cfg.IQTREE_BIN, "-s", str(aln), "-T", str(cfg.THREADS),
           "--prefix", str(prefix), "-seed", str(cfg.IQTREE_SEED), "-redo"]
    if cfg.IQTREE_MODEL:
        cmd += ["-m", cfg.IQTREE_MODEL]
    if _bootstrap_replicates(cfg) > 0 and n_seqs >= 4:
        cmd += ["-B", str(cfg.IQTREE_BOOTSTRAP)]
    return cmd


def run(cfg: Config) -> None:
    """Build a per-family ML tree from each MAFFT alignment."""
    cfg.ensure_dirs()
    external.require(cfg.IQTREE_BIN)
    msa_dir = cfg.result("msa")
    trees_dir = cfg.result("trees")
    trees_dir.mkdir(parents=True, exist_ok=True)

    alignments = sorted(msa_dir.glob("*_aligned.fasta"))
    if not alignments:
        external.log(f"[WARN] No *_aligned.fasta in {msa_dir}; run msa first.")
        return

    min_seqs = 4 if _bootstrap_replicates(cfg) > 0 else 3
    external.log("[iqtree] Building per-family trees...")
    for aln in alignments:
        family = aln.name[: -len("_aligned.fasta")]
        tree_in = _tree_input(aln)          # trimmed alignment if trimAl ran, else raw
        n_seqs = sum(1 for _ in SeqIO.parse(str(tree_in), "fasta"))
        if n_seqs < min_seqs:
            external.log(f"[SKIP] {family}: {n_seqs} sequence(s), need >={min_seqs} for a tree")
            continue
        src = "trimmed" if tree_in is not aln else "untrimmed"
        external.log(f"[iqtree] {family}: {n_seqs} sequences ({src} alignment)...")
        external.run(_iqtree_cmd(cfg, tree_in, trees_dir / family, n_seqs))
        external.log(f"[OK] {family}: tree -> {family}.treefile")

    external.log("[iqtree] IQ-TREE step done.")
