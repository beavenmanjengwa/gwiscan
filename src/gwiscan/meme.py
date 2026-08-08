#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# meme.py - Per-family motif discovery with the MEME Suite (the `meme` stage).                     #
#                                                                                                  #
# Runs MEME on each family's clean domain FASTA and writes the standard MEME output                #
# (meme.html / meme.xml / meme.txt + motif logos) to intermediate/meme/{family}/. Only the number       #
# of motifs is opinionated (MEME_NMOTIFS, default 15); everything else uses MEME's own defaults.   #
#                                                                                                  #
# MEME's install puts binaries in two folders, so add BOTH to PATH so `meme` resolves by name:     #
#   export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-<version>:$PATH                             #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

from Bio import SeqIO

from . import external
from .config import Config


def run(cfg: Config) -> None:
    """Run MEME motif discovery on each family's domain FASTA."""
    cfg.ensure_dirs()
    external.require(cfg.MEME_BIN)
    domain_dir = cfg.result("domain_sequences")
    meme_dir = cfg.result("meme")
    meme_dir.mkdir(parents=True, exist_ok=True)

    domain_fastas = sorted(domain_dir.glob("*_domains.fasta"))
    if not domain_fastas:
        external.log(f"[WARN] No *_domains.fasta in {domain_dir}; run extract-domains first.")
        return

    external.log(f"[meme] Running MEME per family (nmotifs={cfg.MEME_NMOTIFS})...")
    for domain_fasta in domain_fastas:
        family = domain_fasta.name[: -len("_domains.fasta")]
        n_seqs = sum(1 for _ in SeqIO.parse(str(domain_fasta), "fasta"))
        if n_seqs < 2:
            external.log(f"[SKIP] {family}: only {n_seqs} sequence(s), need >=2 for MEME")
            continue

        out_family = meme_dir / family
        external.log(f"[meme] {family}: {n_seqs} sequences -> {out_family.name}/")
        # MEME is a leaf stage (nothing consumes its output), and its post-run image
        # conversion shells out to Ghostscript, which fails on some environments
        # (e.g. a project path containing spaces: "undefinedfilename"). A MEME failure
        # must not sink the whole species -- warn and move to the next family.
        try:
            external.run([
                cfg.MEME_BIN, domain_fasta,
                "-protein",
                "-nmotifs", cfg.MEME_NMOTIFS,
                "-oc", out_family,
            ])
        except RuntimeError as e:
            external.log(f"[WARN] {family}: MEME failed, skipping (motifs are optional). {e}")
            continue
        external.log(f"[OK] {family}: MEME written -> {out_family.name}/")

    external.log("[meme] MEME step done.")
