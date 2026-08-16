#!/usr/bin/env python3
"""
#####################################################################################################
#                                                                                                   #
# clipkit.py - Trim each MAFFT alignment with ClipKIT before tree-building (the `trim` stage).      #
#                                                                                                   #
# Runs after msa and before iqtree: for every msa/{family}_aligned.fasta it writes a cleaned        #
# msa/{family}_trimmed.fasta using ClipKIT's gap-aware column filter (CLIPKIT_MODE, default         #
# smart-gap). Only IQ-TREE consumes the trimmed alignment (poorly-aligned columns hurt tree         #
# inference); WebLogo and MEME deliberately keep using the full alignment so logos/motifs show      #
# the complete domain. ClipKIT is optional: if it isn't installed the `run` stage auto-skips and    #
# iqtree falls back to the untrimmed alignment.                                                     #
#                                                                                                   #
#####################################################################################################
"""

from __future__ import annotations

from . import external
from .config import Config

_ALIGNED_SUFFIX = "_aligned.fasta"
_TRIMMED_SUFFIX = "_trimmed.fasta"


def trimmed_path(aligned):
    """msa/{family}_aligned.fasta -> msa/{family}_trimmed.fasta (same directory)."""
    return aligned.with_name(aligned.name[: -len(_ALIGNED_SUFFIX)] + _TRIMMED_SUFFIX)


def _clipkit_cmd(cfg: Config, aln, out) -> list:
    """ClipKIT command: the input alignment is a positional argument, the cleaned
    alignment is written with -o in the same (FASTA) format, and -m sets the
    trimming mode (CLIPKIT_MODE, default smart-gap)."""
    mode = str(cfg.CLIPKIT_MODE or "smart-gap").strip()
    return [cfg.CLIPKIT_BIN, str(aln), "-m", mode, "-o", str(out)]


def run(cfg: Config) -> None:
    """Trim every family's MAFFT alignment with ClipKIT."""
    cfg.ensure_dirs()
    external.require(cfg.CLIPKIT_BIN)
    msa_dir = cfg.result("msa")

    alignments = sorted(msa_dir.glob(f"*{_ALIGNED_SUFFIX}"))
    if not alignments:
        external.log(f"[WARN] No *{_ALIGNED_SUFFIX} in {msa_dir}; run msa first.")
        return

    mode = str(cfg.CLIPKIT_MODE or "smart-gap").strip()
    external.log(f"[trim] Trimming {len(alignments)} alignment(s) with ClipKIT (-m {mode})...")
    for aln in alignments:
        out = trimmed_path(aln)
        external.run(_clipkit_cmd(cfg, aln, out))
        external.log(f"[OK] {aln.name} -> {out.name}")

    external.log("[trim] ClipKIT step done.")
