#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# confirm.py - The InterProScan-confirmed FINAL candidate set (the `confirm` gate).                #
#                                                                                                  #
# Search gives initial candidates; InterProScan is the gate that turns them into final ones. A      #
# merged (protein, family) row is kept when the family's Pfam is reported on that protein, or when   #
# the family has no Pfam (EUL / custom-HMM), which is kept as already confirmed. Outputs the        #
# confirmed table + sequences that the per-candidate annotators read.                              #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

from collections import defaultdict

from Bio import SeqIO

from . import external, io
from .config import Config
from .domain_bed import family_pfam_map


def _bare(acc) -> str:
    return str(acc).split(".")[0]


def confirmed_rows(merged_df, interpro_df, family_to_pfam):
    """Keep (protein, family) rows confirmed by InterProScan; families with no Pfam kept."""
    pfams_by_protein = defaultdict(set)
    for _, r in interpro_df.iterrows():
        pfams_by_protein[str(r["protein_id"])].add(_bare(r["sig_acc"]))

    keep = []
    for idx, r in merged_df.iterrows():
        pfam = family_to_pfam.get(r["family"])   # None for EUL / custom-HMM families
        if pfam is None or pfam in pfams_by_protein.get(str(r["protein_id"]), ()):
            keep.append(idx)
    return merged_df.loc[keep]


def run(cfg: Config) -> None:
    """Gate the merged candidates through InterProScan into the final set."""
    cfg.ensure_dirs()
    merged = cfg.result("candidates_merged.tsv")
    interpro = cfg.result("interproscan.tsv")
    cand_fasta = cfg.result("candidates.fasta")
    out_tsv = cfg.result("final_candidates.tsv")
    out_fasta = cfg.result("final_candidates.fasta")

    if not merged.exists():
        raise FileNotFoundError(f"candidates_merged.tsv not found: {merged} (run `merge` first)")
    if not interpro.exists():
        raise FileNotFoundError(f"interproscan.tsv not found: {interpro} (run `interpro` first)")

    merged_df = io.read_tsv(merged)
    interpro_df = io.read_tsv(interpro)
    family_to_pfam = family_pfam_map(cfg.family_map)

    final_df = confirmed_rows(merged_df, interpro_df, family_to_pfam)
    io.write_df(final_df, out_tsv, "tsv")

    proteome = SeqIO.to_dict(SeqIO.parse(str(cand_fasta), "fasta"))
    ids = list(dict.fromkeys(final_df["protein_id"].astype(str)))   # unique, order-preserving
    records = [proteome[i] for i in ids if i in proteome]
    SeqIO.write(records, str(out_fasta), "fasta")

    dropped = merged_df["protein_id"].nunique() - len(ids)
    external.log(
        f"[OK] Final candidates: {len(ids)} proteins confirmed "
        f"({len(final_df)} rows; {dropped} proteins dropped) -> final_candidates.tsv/.fasta"
    )
