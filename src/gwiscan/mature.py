#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# mature.py - Per-family mature-sequence FASTAs from TargetP (the `extract-mature` stage).         #
#                                                                                                  #
# TargetP runs with -mature and writes intermediate/targetp_raw/*_mature.fasta, the full sequences with #
# the signal/transit presequence cleaved off. This groups those mature sequences by family (from   #
# final_candidates) into intermediate/mature_sequences/{family}_mature.fasta, so msa + iqtree build a   #
# mature full-sequence tree per family.                                                            #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

from collections import defaultdict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from . import external, io
from .config import Config


def find_mature_fasta(cfg: Config):
    """TargetP's mature-sequence FASTA under intermediate/targetp_raw/, or None."""
    raw_dir = cfg.result("targetp_raw")
    if not raw_dir.exists():
        return None
    matches = sorted(raw_dir.glob("*mature*.fasta"))
    return matches[0] if matches else None


def group_by_family(final_df) -> dict:
    """{family: [protein_id, ...]} from final_candidates, unique and ordered."""
    by_family = defaultdict(list)
    for pid, family in zip(final_df["protein_id"].astype(str),
                           final_df["family"].astype(str)):
        if pid not in by_family[family]:
            by_family[family].append(pid)
    return by_family


def run(cfg: Config) -> None:
    """Write per-family FASTAs of the TargetP mature sequences."""
    cfg.ensure_dirs()
    src = find_mature_fasta(cfg)
    if src is None:
        external.log("[WARN] No TargetP mature FASTA in intermediate/targetp_raw/; run targetp first.")
        return

    final = cfg.result("final_candidates.tsv")
    if not final.exists():
        raise FileNotFoundError(
            f"final_candidates.tsv not found: {final} (run `gwiscan confirm` first)")

    final_df = io.read_tsv(final)
    mature = SeqIO.to_dict(SeqIO.parse(str(src), "fasta"))
    out_dir = cfg.result("mature_sequences")
    out_dir.mkdir(parents=True, exist_ok=True)

    total, missing = 0, 0
    for family, protein_ids in sorted(group_by_family(final_df).items()):
        records = []
        for pid in protein_ids:
            rec = mature.get(pid)
            if rec is None:
                missing += 1
                continue
            records.append(SeqRecord(rec.seq, id=pid, description=""))
        if records:
            SeqIO.write(records, str(out_dir / f"{family}_mature.fasta"), "fasta")
            total += len(records)
            external.log(f"[OK] {family}: {len(records)} mature sequences")

    if missing:
        external.log(f"[WARN] {missing} candidate(s) absent from the TargetP mature FASTA")
    external.log(f"[OK] Mature sequences -> mature_sequences/ ({total} sequences)")
