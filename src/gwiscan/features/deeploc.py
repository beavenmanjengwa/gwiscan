#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# deeploc.py - Subcellular localization with DeepLoc 2.1 (the `deeploc` stage).                     #
#                                                                                                  #
# DeepLoc 2.1 is a user-installed academic download; point at it via DEEPLOC_BIN (PATH name or      #
# path). Runs -m DEEPLOC_MODEL, default Accurate (ProtT5, ~32GB RAM, downloaded on first use);      #
# set Fast (ESM1b) on a memory-limited machine. The full results CSV is kept under                 #
# intermediate/deeploc_out/; deeploc.tsv carries protein_id, Localizations, Signals, Membrane types.    #
# DeepLoc's Lysosome/Vacuole class is collapsed to Vacuole (plants have vacuoles).                 #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path

import pandas as pd

from .. import external, io
from ..config import Config

# DeepLoc results columns. The ID column plus these three are kept in
# deeploc.tsv; the rest (per-compartment probabilities) stay in the raw CSV.
KEEP_COLUMNS = ["Localizations", "Signals", "Membrane types"]

# DeepLoc's "Lysosome/Vacuole" class -> "Vacuole" (either order, any case/spacing).
_LYSO_RE = re.compile(r"lysosome\s*/\s*vacuole|vacuole\s*/\s*lysosome", re.IGNORECASE)


def _clean_localization(value):
    if not isinstance(value, str):
        return value
    return _LYSO_RE.sub("Vacuole", value)


def parse_results(path) -> pd.DataFrame:
    """Keep protein_id + the DeepLoc summary columns; collapse Lysosome/Vacuole."""
    df = pd.read_csv(path).rename(columns={"Protein_ID": "protein_id"})
    df = df[["protein_id", *KEEP_COLUMNS]]
    df["Localizations"] = df["Localizations"].map(_clean_localization)
    return df


def run(cfg: Config) -> None:
    cfg.ensure_dirs()
    cand_fasta = cfg.result("final_candidates.fasta")
    out_dir = cfg.result("deeploc_out")
    out_tsv = cfg.result("deeploc.tsv")
    out_dir.mkdir(parents=True, exist_ok=True)

    if not cand_fasta.exists():
        raise FileNotFoundError(
            f"final_candidates.fasta not found: {cand_fasta} (run `gwiscan confirm` first)"
        )

    # deeploc_bin may be a PATH name or an absolute path to the executable.
    if shutil.which(cfg.DEEPLOC_BIN) is None and not Path(cfg.DEEPLOC_BIN).exists():
        raise FileNotFoundError(
            f"DeepLoc executable not found: {cfg.DEEPLOC_BIN!r}. DeepLoc 2.1 requires a "
            f"free academic license — download and install it, then set deeploc_bin / "
            f"DEEPLOC_BIN to the 'deeploc2' executable. "
            f"https://services.healthtech.dtu.dk/services/DeepLoc-2.1/"
        )

    # Clear our output folder first so the only CSV left afterwards is this run's.
    for old in out_dir.glob("*.csv"):
        old.unlink()

    external.log(f"[deeploc] Running DeepLoc 2.1 (-m {cfg.DEEPLOC_MODEL})...")
    external.run([
        cfg.DEEPLOC_BIN,
        "-f", cand_fasta.resolve(),
        "-o", out_dir,
        "-m", cfg.DEEPLOC_MODEL,
    ])

    csvs = list(out_dir.glob("*.csv"))
    if not csvs:
        raise FileNotFoundError(
            f"DeepLoc produced no results CSV in {out_dir} (check the output above)"
        )

    df = parse_results(csvs[0])
    io.write_df(df, out_tsv, "tsv")
    external.log(f"[OK] DeepLoc: {len(df)} entries -> deeploc.tsv")
