#!/usr/bin/env python3
"""
######################################################################################################
#                                                                                                    #
# confirm.py - The InterProScan-confirmed FINAL candidate set (the `confirm` gate).                  #
#                                                                                                    #
# Search gives initial candidates; InterProScan is the gate that turns them into final ones. A       #
# merged (protein, family) row is kept when the protein carries any of the family's confirmation      #
# accessions -- the InterProModel column (Pfam, CDD or PANTHER), or the PfamModel accession when      #
# that column is absent. A family with no accessions (a custom-HMM family)                            #
# is kept as already confirmed. Outputs the confirmed table + sequences the annotators read.          #
#                                                                                                    #
######################################################################################################
"""

from __future__ import annotations

from collections import defaultdict

from Bio import SeqIO

from . import external, io
from .config import Config


def _bare(acc) -> str:
    return str(acc).split(".")[0]


def confirmed_rows(merged_df, interpro_df, family_to_accs):
    """Keep (protein, family) rows confirmed by InterProScan.

    A family's confirmation accessions come from its InterProModel column, or the
    PfamModel accession if that column is absent (see io.family_confirm_accessions).
    A row is kept when the protein carries ANY of the family's accessions -- across
    whichever member database (Pfam, CDD or PANTHER) those accessions name -- so a
    family with no usable Pfam still confirms on a CDD or PANTHER id. A family with
    no accessions at all (a custom-HMM family) is not gated here and is kept.

    An accession absent from ``interpro_df`` is treated as a genuine non-match, not a
    missing annotation: the interpro stage annotates every candidate completely (the
    API path fails the run if any chunk does not finish, and the local path runs the
    whole set in one call), so an absent accession means "InterProScan saw this
    protein and did not call it", never "the annotation is missing".
    """
    sigs_by_protein = defaultdict(set)
    for _, r in interpro_df.iterrows():
        sigs_by_protein[str(r["protein_id"])].add(_bare(r["sig_acc"]))

    keep = []
    for idx, r in merged_df.iterrows():
        accs = family_to_accs.get(r["family"], set())   # empty for custom-HMM families
        if not accs or (accs & sigs_by_protein.get(str(r["protein_id"]), set())):
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
    interpro_df = io.read_interpro_tsv(interpro)
    family_to_accs = io.family_confirm_accessions(cfg.family_map)

    final_df = confirmed_rows(merged_df, interpro_df, family_to_accs)
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
