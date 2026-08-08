#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# domains.py - Extract and name per-family domain sequences (the `extract-domains` stage).         #
#                                                                                                  #
# Reads intermediate/domains.bed, slices each domain out of candidates.fasta, and writes per family     #
# into intermediate/domain_sequences/: clean systematic ids (for MSA/trees) and provenance-headed       #
# copies, plus domain_map.tsv (the domain_id/protein/family/pfam/coord crosswalk).                 #
#                                                                                                  #
# Clean domain id = {species_prefix}{Family}{serial}.{domain}, e.g. AthGNA001.1. Coordinates are   #
# 1-based (InterProScan for Pfam families, the family's own hmmscan hit for custom-HMM families).  #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import external, io
from .config import Config

MAP_HEADER = ["domain_id", "protein_id", "family", "pfam_id", "start", "end"]


def assign_domain_ids(bed_df, pfam_to_family, species_prefix):
    """Assign systematic domain ids. Returns a list of dicts with MAP_HEADER keys.

    Grouped by family; within a family, proteins are sorted and numbered (serial);
    within a protein, domains are ordered by start and numbered.
    """
    by_family = defaultdict(lambda: defaultdict(list))
    for _, r in bed_df.iterrows():
        pfam = str(r["pfam_id"])
        family = pfam_to_family.get(pfam, pfam)
        by_family[family][str(r["protein_id"])].append((int(r["start"]), int(r["end"]), pfam))

    records = []
    for family in sorted(by_family):
        for serial, pid in enumerate(sorted(by_family[family]), start=1):
            for d_idx, (start, end, pfam) in enumerate(sorted(by_family[family][pid]), start=1):
                records.append({
                    "domain_id": f"{species_prefix}{family}{serial:03d}.{d_idx}",
                    "protein_id": pid,
                    "family": family,
                    "pfam_id": pfam,
                    "start": start,
                    "end": end,
                })
    return records


def run(cfg: Config) -> None:
    """Slice and name each family's domain sequences from candidates.fasta."""
    cfg.ensure_dirs()
    bed = cfg.result("domains.bed")
    cand_fasta = cfg.result("candidates.fasta")
    out_dir = cfg.result("domain_sequences")
    out_dir.mkdir(parents=True, exist_ok=True)

    if not bed.exists():
        raise FileNotFoundError(f"domains.bed not found: {bed} (run `gwiscan domain-bed` first)")
    if not cand_fasta.exists():
        raise FileNotFoundError(f"candidates.fasta not found: {cand_fasta}")

    bed_df = io.read_tsv(bed)
    proteome = SeqIO.to_dict(SeqIO.parse(str(cand_fasta), "fasta"))
    # Architecture mode: the BED's pfam_id column already holds the architecture name
    # (coordinates come from the members' own hmmscan hits), so no map is needed.
    pfam_to_family = {} if cfg.is_architecture else io.pfam_to_family(cfg.family_map)

    records = assign_domain_ids(bed_df, pfam_to_family, cfg.SPECIES_PREFIX)

    clean = defaultdict(list)
    annotated = defaultdict(list)
    skipped = 0
    for rec in records:
        pid = rec["protein_id"]
        if pid not in proteome:
            skipped += 1
            continue
        full = str(proteome[pid].seq)
        start = max(0, rec["start"] - 1)     # 1-based inclusive -> 0-based slice
        end = min(len(full), rec["end"])
        if end <= start:
            skipped += 1
            continue
        seq = Seq(full[start:end])
        family = rec["family"]
        clean[family].append(SeqRecord(seq, id=rec["domain_id"], description=""))
        annotated[family].append(
            SeqRecord(seq, id=f'{pid}:{rec["start"]}-{rec["end"]}', description=f'| {rec["pfam_id"]}')
        )

    for family in sorted(clean):
        SeqIO.write(clean[family], str(out_dir / f"{family}_domains.fasta"), "fasta")
        SeqIO.write(annotated[family], str(out_dir / f"{family}_domains.annotated.fasta"), "fasta")

    io.write_tsv(out_dir / "domain_map.tsv", MAP_HEADER,
                 [[rec[c] for c in MAP_HEADER] for rec in records])

    external.log(
        f"[OK] Extracted {len(records) - skipped} domains across {len(clean)} families "
        + (f"({skipped} skipped) " if skipped else "")
        + "-> domain_sequences/ (clean + annotated FASTA, domain_map.tsv)"
    )
    for family in sorted(clean):
        external.log(f"     {family}: {len(clean[family])} domains")
