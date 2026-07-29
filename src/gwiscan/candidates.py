#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# candidates.py - Merge HMM and DIAMOND hits into the candidate set (the `merge` stage).           #
#                                                                                                  #
# Concatenates hmm_hits.tsv and blast_hits.tsv, sorts best-hit-first per protein/family, and       #
# drops exact duplicate domain rows (keeping every family a protein hits). Writes the merged       #
# table and the candidate protein sequences pulled from the proteome.                              #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

import pandas as pd
from Bio import SeqIO

from . import external, io
from .config import Config


def run(cfg: Config) -> None:
    """Merge HMM + DIAMOND hits and extract the candidate sequences."""
    cfg.ensure_dirs()
    hmm_hits = cfg.result("hmm_hits.tsv")
    blast_hits = cfg.result("blast_hits.tsv")
    merged_tsv = cfg.result("candidates_merged.tsv")
    candidates_fasta = cfg.result("candidates.fasta")

    frames = []
    if hmm_hits.exists():
        hmm_df = io.read_tsv(hmm_hits)
        frames.append(hmm_df)
        external.log(f"[OK] HMM hits loaded: {len(hmm_df)} domain hits")
    else:
        external.log("[WARN] hmm_hits.tsv not found, skipping.")

    if blast_hits.exists():
        blast_df = io.read_tsv(blast_hits)
        if len(blast_df) > 0:
            frames.append(blast_df)
            external.log(f"[OK] BLAST hits loaded: {len(blast_df)} hits")
        else:
            external.log("[WARN] blast_hits.tsv is empty, skipping.")
    else:
        external.log("[WARN] blast_hits.tsv not found, skipping.")

    if not frames:
        raise RuntimeError("No hits found from either search. Run a search stage first.")

    merged = pd.concat(frames, ignore_index=True)
    merged["evalue"] = pd.to_numeric(merged["evalue"], errors="coerce")
    merged["bitscore"] = pd.to_numeric(merged["bitscore"], errors="coerce")

    # Best hit first per protein/family; keep all families, drop exact duplicate rows.
    merged.sort_values(["protein_id", "family", "evalue"], inplace=True)
    merged.drop_duplicates(
        subset=["protein_id", "family", "start", "end"], inplace=True
    )

    io.write_df(merged, merged_tsv, "tsv")
    external.log(
        f"[OK] Merged candidates: {len(merged)} domain hits across "
        f"{merged['protein_id'].nunique()} proteins -> candidates_merged.tsv"
    )

    # Extract the unique candidate sequences from the proteome. Order-preserving
    # dedup (merged is already deterministically sorted) so candidates.fasta has a
    # stable, reproducible row order run-to-run.
    candidate_ids = list(dict.fromkeys(merged["protein_id"].astype(str)))
    proteome = SeqIO.to_dict(SeqIO.parse(str(cfg.proteome), "fasta"))
    found = [pid for pid in candidate_ids if pid in proteome]
    missing = [pid for pid in candidate_ids if pid not in proteome]

    SeqIO.write([proteome[pid] for pid in found], str(candidates_fasta), "fasta")
    external.log(f"[OK] Extracted {len(found)} candidate sequences -> candidates.fasta")

    if missing:
        external.log(f"[WARN] {len(missing)} protein IDs not found in proteome.fasta:")
        for pid in missing[:10]:
            external.log(f"       {pid}")
        if len(missing) > 10:
            external.log(f"       ... and {len(missing) - 10} more")
