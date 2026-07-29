#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# domain_bed.py - BED-like table of family domain coordinates (the `domain-bed` stage).            #
#                                                                                                  #
# For each candidate, emits rows [protein_id, pfam_id, start, end]. Pfam families take their       #
# coordinates from the InterProScan match of their own family Pfam; custom-HMM families take       #
# them from their own hmmscan hit (start/end already in candidates_merged), with no InterProScan.  #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

from . import external, io
from .config import Config

BED_HEADER = ["protein_id", "pfam_id", "start", "end"]


def _versionless(acc) -> str:
    return str(acc).split(".")[0]


def family_pfam_map(family_map_path) -> dict:
    """family -> bare Pfam accession, for families that have a PfamModel."""
    return {r["family"]: r["pfam_model"]
            for r in io.family_records(family_map_path) if r["pfam_model"]}


def build_bed(candidates_df, family_to_pfam, interpro_df, custom_families=()) -> list:
    """Rows [protein_id, pfam_id, start, end].

    Two coordinate sources, by how the family is identified:
      * Pfam families — one row per InterProScan match of the candidate's own
        family Pfam (PF00704 for CRA on the Pfam route, etc.).
      * custom-HMM families — one row per hmmscan hit of that family, using the
        model's own alignment coordinates (start/end already in
        candidates_merged). No InterProScan; pfam_id is the family name.
    """
    custom_families = set(custom_families)

    # Index InterProScan Pfam matches by (protein_id, bare pfam) -> [(start, end)].
    matches = {}
    for _, r in interpro_df.iterrows():
        key = (str(r["protein_id"]), _versionless(r["sig_acc"]))
        matches.setdefault(key, []).append((int(r["start"]), int(r["end"])))

    rows = []
    seen = set()
    for _, c in candidates_df.iterrows():
        pid = str(c["protein_id"])
        family = c["family"]

        # Custom-HMM family: coordinates from its own hmmscan hit.
        if family in custom_families:
            if str(c.get("method")) != "hmm":
                continue
            try:
                start, end = int(c["start"]), int(c["end"])
            except (ValueError, TypeError, KeyError):
                continue
            if (pid, family, start, end) in seen:
                continue
            seen.add((pid, family, start, end))
            rows.append([pid, family, start, end])
            continue

        # Pfam family: coordinates from InterProScan matches of the family Pfam.
        pfam = family_to_pfam.get(family)
        if not pfam or (pid, pfam) in seen:
            continue
        seen.add((pid, pfam))
        for start, end in matches.get((pid, pfam), []):
            rows.append([pid, pfam, start, end])
    return rows


def run(cfg: Config) -> None:
    """Emit each candidate's family domain coordinates to intermediate/domains.bed."""
    cfg.ensure_dirs()
    merged = cfg.result("candidates_merged.tsv")
    interpro = cfg.result("interproscan.tsv")
    out_bed = cfg.result("domains.bed")

    if not merged.exists():
        raise FileNotFoundError(f"candidates_merged.tsv not found: {merged} (run `merge` first)")
    if not interpro.exists():
        raise FileNotFoundError(f"interproscan.tsv not found: {interpro} (run `interpro` first)")

    candidates_df = io.read_tsv(merged)
    interpro_df = io.read_tsv(interpro)
    family_to_pfam = family_pfam_map(cfg.family_map)
    custom_families = io.custom_hmm_families(cfg.family_map)

    rows = build_bed(candidates_df, family_to_pfam, interpro_df, custom_families)
    io.write_tsv(out_bed, BED_HEADER, rows)

    n_prot = len({r[0] for r in rows})
    external.log(f"[OK] Family domains: {len(rows)} domain rows "
                 f"({n_prot} proteins) -> domains.bed")
