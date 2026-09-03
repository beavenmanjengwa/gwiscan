#!/usr/bin/env python3
"""
###########################################################################################################
#                                                                                                         #
# meme.py - Per-family motif discovery with the MEME Suite (the `meme` stage).                            #
#                                                                                                         #
# Runs MEME on each family's extracted domains AND its TargetP mature sequences (the same                 #
# two tracks the msa / iqtree stages use) and writes the standard MEME output (meme.html /                #
# meme.xml / meme.txt + motif logos) to intermediate/meme/{family}/ and {family}_mature/. Only the number #
# of motifs is opinionated (MEME_NMOTIFS, default 15); everything else uses MEME's own defaults.          #
#                                                                                                         #
# MEME's install puts binaries in two folders, so add BOTH to PATH so `meme` resolves by name:            #
#   export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-<version>:$PATH                                    #
#                                                                                                         #
###########################################################################################################
"""

from __future__ import annotations

from Bio import SeqIO

from . import external
from .config import Config


def run(cfg: Config) -> None:
    """Run MEME motif discovery on each family's domain and mature FASTA."""
    cfg.ensure_dirs()
    external.require(cfg.MEME_BIN)
    meme_dir = cfg.result("meme")
    meme_dir.mkdir(parents=True, exist_ok=True)

    # Both tracks, like the msa / iqtree stages: the extracted family domains and the
    # TargetP mature sequences. Output folders are {family}/ and {family}_mature/.
    targets = (sorted(cfg.result("domain_sequences").glob("*_domains.fasta"))
               + sorted(cfg.result("mature_sequences").glob("*_mature.fasta")))
    if not targets:
        external.log("[WARN] No *_domains.fasta or *_mature.fasta; "
                     "run extract-domains / extract-mature first.")
        return

    width = ""
    if str(cfg.MEME_MINW).strip() or str(cfg.MEME_MAXW).strip():
        width = f", width {cfg.MEME_MINW or '6'}-{cfg.MEME_MAXW or '50'}"
    external.log(f"[meme] Running MEME per family (nmotifs={cfg.MEME_NMOTIFS}{width})...")
    for fasta in targets:
        label = fasta.stem.replace("_domains", "")   # {family} or {family}_mature
        n_seqs = sum(1 for _ in SeqIO.parse(str(fasta), "fasta"))
        if n_seqs < 2:
            external.log(f"[SKIP] {label}: only {n_seqs} sequence(s), need >=2 for MEME")
            continue

        out_family = meme_dir / label
        external.log(f"[meme] {label}: {n_seqs} sequences -> {out_family.name}/")
        # MEME is a leaf stage (nothing consumes its output), and its post-run image
        # conversion shells out to Ghostscript, which fails on some environments
        # (e.g. a project path containing spaces: "undefinedfilename"). A MEME failure
        # must not sink the whole species -- warn and move to the next family.
        cmd = [cfg.MEME_BIN, fasta, "-protein", "-nmotifs", cfg.MEME_NMOTIFS, "-oc", out_family]
        # Optional motif-width bounds; empty => MEME's own defaults (minw 6, maxw 50).
        if str(cfg.MEME_MINW).strip():
            cmd += ["-minw", str(cfg.MEME_MINW).strip()]
        if str(cfg.MEME_MAXW).strip():
            cmd += ["-maxw", str(cfg.MEME_MAXW).strip()]
        try:
            external.run(cmd)
        except RuntimeError as e:
            external.log(f"[WARN] {label}: MEME failed, skipping (motifs are optional). {e}")
            continue
        external.log(f"[OK] {label}: MEME written -> {out_family.name}/")

    external.log("[meme] MEME step done.")
