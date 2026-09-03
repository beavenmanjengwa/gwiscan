#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# msa.py - Per-family multiple sequence alignment with MAFFT (the `msa` stage).                    #
#                                                                                                  #
# Aligns two independent tracks with MAFFT (strategy from MAFFT_ALGORITHM, default --auto):         #
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


# MAFFT strategy presets (MAFFT_ALGORITHM) -> the flags that select them. `auto`
# lets MAFFT choose by input size; linsi/ginsi/einsi are the accurate iterative
# strategies (L/G/E-INS-i) for divergent sequences; fftns2 is the fast progressive
# strategy for large sets. Aliases accept both the short and canonical spellings.
_MAFFT_STRATEGIES = {
    "auto": ["--auto"],
    "linsi": ["--localpair", "--maxiterate", "1000"],
    "l-ins-i": ["--localpair", "--maxiterate", "1000"],
    "ginsi": ["--globalpair", "--maxiterate", "1000"],
    "g-ins-i": ["--globalpair", "--maxiterate", "1000"],
    "einsi": ["--genafpair", "--maxiterate", "1000"],
    "e-ins-i": ["--genafpair", "--maxiterate", "1000"],
    "fftns2": ["--retree", "2"],
    "fft-ns-2": ["--retree", "2"],
    "fftnsi": ["--retree", "2", "--maxiterate", "2"],
}


def _mafft_strategy_args(cfg: Config) -> list:
    """The MAFFT flags for the configured strategy (MAFFT_ALGORITHM). Unknown names
    fall back to --auto with a warning rather than passing an invalid flag."""
    key = str(cfg.MAFFT_ALGORITHM or "auto").strip().lower()
    if key not in _MAFFT_STRATEGIES:
        external.log(f"[WARN] unknown MAFFT_ALGORITHM {cfg.MAFFT_ALGORITHM!r}; using --auto "
                     f"(valid: {', '.join(sorted(_MAFFT_STRATEGIES))}).")
        key = "auto"
    return _MAFFT_STRATEGIES[key]


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

    strategy_args = _mafft_strategy_args(cfg)
    external.log(f"[msa] Running MAFFT per family (strategy: {cfg.MAFFT_ALGORITHM})...")
    for fasta in targets:
        n_seqs = sum(1 for _ in SeqIO.parse(str(fasta), "fasta"))
        if n_seqs < 2:
            external.log(f"[SKIP] {fasta.name}: only {n_seqs} sequence(s), need >=2 for MSA")
            continue

        out_msa = msa_dir / _aligned_name(fasta)
        external.log(f"[msa] {fasta.name}: aligning {n_seqs} sequences...")
        # MAFFT writes the alignment to stdout; keep it out of the log stream.
        external.run(
            ["mafft", *strategy_args, "--thread", cfg.THREADS, "--reorder", fasta],
            stdout_path=out_msa,
        )
        external.log(f"[OK] {fasta.name}: MSA written -> {out_msa.name}")

    external.log("[msa] MAFFT step done.")
