#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# msa.py - Per-family multiple sequence alignment with MAFFT (the `msa` stage).                    #
#                                                                                                  #
# Aligns two independent tracks with `mafft --auto`:                                               #
#   * domain_sequences/{family}_domains.fasta -> msa/{family}_aligned.fasta                        #
#   * mature_sequences/{family}_mature.fasta  -> msa/{family}_mature_aligned.fasta                 #
# The domain track aligns the extracted family domains; the mature track aligns the full TargetP   #
# mature sequences (presequence cleaved). Families with fewer than two sequences are skipped.      #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

from Bio import SeqIO

from . import external
from .config import Config


def _aligned_name(fasta) -> str:
    """{family}_domains.fasta -> {family}_aligned.fasta;
    {family}_mature.fasta  -> {family}_mature_aligned.fasta."""
    return fasta.stem.replace("_domains", "") + "_aligned.fasta"


def run(cfg: Config) -> None:
    """Align each family's domain and mature sequences with MAFFT."""
    cfg.ensure_dirs()
    external.require("mafft")
    msa_dir = cfg.result("msa")
    msa_dir.mkdir(parents=True, exist_ok=True)

    targets = (sorted(cfg.result("domain_sequences").glob("*_domains.fasta"))
               + sorted(cfg.result("mature_sequences").glob("*_mature.fasta")))
    if not targets:
        external.log("[WARN] No *_domains.fasta or *_mature.fasta; "
                     "run extract-domains / extract-mature first.")
        return

    external.log("[msa] Running MAFFT per family...")
    for fasta in targets:
        n_seqs = sum(1 for _ in SeqIO.parse(str(fasta), "fasta"))
        if n_seqs < 2:
            external.log(f"[SKIP] {fasta.name}: only {n_seqs} sequence(s), need >=2 for MSA")
            continue

        out_msa = msa_dir / _aligned_name(fasta)
        external.log(f"[msa] {fasta.name}: aligning {n_seqs} sequences...")
        # MAFFT writes the alignment to stdout; keep it out of the log stream.
        external.run(
            ["mafft", "--auto", "--thread", cfg.THREADS, "--reorder", fasta],
            stdout_path=out_msa,
        )
        external.log(f"[OK] {fasta.name}: MSA written -> {out_msa.name}")

    external.log("[msa] MAFFT step done.")
