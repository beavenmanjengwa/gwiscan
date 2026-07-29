#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# trimal.py - Trim each MAFFT alignment with trimAl before tree-building (the `trim` stage).       #
#                                                                                                  #
# Runs after msa and before iqtree: for every msa/{family}_aligned.fasta it writes a cleaned       #
# msa/{family}_trimmed.fasta using trimAl's heuristic column filter (TRIMAL_METHOD, default        #
# -automated1). Only IQ-TREE consumes the trimmed alignment (poorly-aligned columns hurt tree      #
# inference); WebLogo and MEME deliberately keep using the full alignment so logos/motifs show     #
# the complete domain. trimAl is optional: if it isn't installed the `run` stage auto-skips and    #
# iqtree falls back to the untrimmed alignment.                                                    #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

from . import external
from .config import Config

_ALIGNED_SUFFIX = "_aligned.fasta"
_TRIMMED_SUFFIX = "_trimmed.fasta"


def trimmed_path(aligned):
    """msa/{family}_aligned.fasta -> msa/{family}_trimmed.fasta (same directory)."""
    return aligned.with_name(aligned.name[: -len(_ALIGNED_SUFFIX)] + _TRIMMED_SUFFIX)


def _method_flag(cfg: Config) -> str:
    """trimAl heuristic as a CLI flag: 'automated1' -> '-automated1'. A method
    already written with a leading '-' is passed through unchanged."""
    method = str(cfg.TRIMAL_METHOD or "automated1").strip()
    return method if method.startswith("-") else f"-{method}"


def run(cfg: Config) -> None:
    """Trim every family's MAFFT alignment with trimAl."""
    cfg.ensure_dirs()
    external.require(cfg.TRIMAL_BIN)
    msa_dir = cfg.result("msa")

    alignments = sorted(msa_dir.glob(f"*{_ALIGNED_SUFFIX}"))
    if not alignments:
        external.log(f"[WARN] No *{_ALIGNED_SUFFIX} in {msa_dir}; run msa first.")
        return

    flag = _method_flag(cfg)
    external.log(f"[trim] Trimming {len(alignments)} alignment(s) with trimAl ({flag})...")
    for aln in alignments:
        out = trimmed_path(aln)
        external.run([cfg.TRIMAL_BIN, "-in", aln, "-out", out, flag])
        external.log(f"[OK] {aln.name} -> {out.name}")

    external.log("[trim] trimAl step done.")
