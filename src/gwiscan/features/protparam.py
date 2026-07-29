#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# protparam.py - Physicochemical properties via Biopython Bio.SeqUtils.ProtParam (`protparam`).    #
#                                                                                                  #
# Per candidate: length, molecular weight (kDa), pI, negatively (Asp+Glu) and positively (Arg+Lys) #
# charged-residue counts, instability index, GRAVY, aliphatic index (Ikai 1980, from mole-percent  #
# amino_acids_percent), and molar extinction coefficient at 280 nm (both Expasy values). Sequences #
# with X or unweighable residues (B/Z/U/O) get '-' in every field. Writes PROTPARAM_FORMATS        #
# (default tsv + xlsx; tsv feeds the join).                                                        #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from .. import external, io
from ..config import Config

COLUMNS = ["protein_id", "length_aa", "molecular_weight", "isoelectric_point",
           "negatively_charged_residues", "positively_charged_residues",
           "instability_index", "gravy", "aliphatic_index",
           "ec_cystines", "ec_reduced"]

_DASHES = {c: "-" for c in COLUMNS[1:]}

# Valid output formats (config `protparam_formats`). Writing (with camelCase
# headers) is delegated to io.write_df.
_FORMATS = ("tsv", "csv", "xlsx")


def properties(seq: str) -> dict:
    """Physicochemical properties for one sequence, or dashes if not analysable."""
    if "X" in seq:
        return dict(_DASHES)
    try:
        pa = ProteinAnalysis(seq)
        aa = pa.amino_acids_percent          # mole percent (0-100) in this Biopython
        aliphatic_index = (aa.get("A", 0)
                           + 2.9 * aa.get("V", 0)
                           + 3.9 * (aa.get("I", 0) + aa.get("L", 0)))
        # Molar extinction coefficient at 280 nm (M-1 cm-1). Biopython returns
        # (reduced, cystines): reduced = all Cys reduced; cystines = all Cys
        # paired into cystines (Expasy's two values).
        ext_reduced, ext_cystines = pa.molar_extinction_coefficient()
        return {
            "length_aa": len(seq),
            "molecular_weight": round(pa.molecular_weight() / 1000, 2),
            "isoelectric_point": round(pa.isoelectric_point(), 2),
            "negatively_charged_residues": seq.count("D") + seq.count("E"),  # Asp + Glu
            "positively_charged_residues": seq.count("R") + seq.count("K"),  # Arg + Lys
            "instability_index": round(pa.instability_index(), 2),
            "gravy": round(pa.gravy(), 3),
            "aliphatic_index": round(aliphatic_index, 2),
            "ec_cystines": ext_cystines,   # EC, all Cys as cystines
            "ec_reduced": ext_reduced,     # EC, all Cys reduced
        }
    except Exception as e:  # noqa: BLE001 - non-standard residue (B/Z/U/O), etc.
        external.log(f"[WARN] ProtParam failed for a sequence: {e}")
        return dict(_DASHES)


def run(cfg: Config) -> None:
    cfg.ensure_dirs()
    cand_fasta = cfg.result("final_candidates.fasta")
    if not cand_fasta.exists():
        raise FileNotFoundError(
            f"final_candidates.fasta not found: {cand_fasta} (run `gwiscan confirm` first)"
        )

    rows = []
    for rec in SeqIO.parse(str(cand_fasta), "fasta"):
        seq = str(rec.seq).upper().replace("*", "").strip()
        if not seq:
            external.log(f"[WARN] skipping empty sequence: {rec.id}")
            continue
        rows.append({"protein_id": rec.id, **properties(seq)})

    df = pd.DataFrame(rows, columns=COLUMNS)

    written = []
    for fmt in cfg.PROTPARAM_FORMATS:
        if fmt not in _FORMATS:
            external.log(f"[WARN] unknown protparam format '{fmt}' (use tsv/xlsx/csv)")
            continue
        io.write_df(df, cfg.result(f"protparam.{fmt}"), fmt)
        written.append(f"protparam.{fmt}")
    external.log(f"[OK] ProtParam: {len(df)} sequences -> {', '.join(written)}")
