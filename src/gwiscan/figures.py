#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# figures.py - Per-family ProtParam distribution figures + summary stats (the `figures` stage).    #
#                                                                                                  #
# Runs after compile. Calls the bundled R/ggplot2 script on final_results/gwiscan_results.tsv to    #
# write, into final_results/: a faceted boxplot of the ProtParam properties (600 dpi PNG/TIFF +     #
# PDF), a per-family summary table, and the list of Tukey outliers. Optional: needs Rscript on      #
# PATH; if it isn't installed the `run` stage auto-skips.                                           #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

from pathlib import Path

from . import external
from .config import Config

# The R script ships inside the package (scripts/), so it is found whether gwiscan
# is installed editable or from a built package.
_R_SCRIPT = Path(__file__).parent / "scripts" / "protparam_boxplots.R"


def run(cfg: Config) -> None:
    """Render the ProtParam boxplots + stats from the compiled results table."""
    cfg.ensure_dirs()
    external.require(cfg.RSCRIPT_BIN)

    results = cfg.final_dir / "gwiscan_results.tsv"
    if not results.exists():
        external.log(f"[WARN] {results} not found; run compile first. Skipping figures.")
        return

    # The R script reads gwiscan_results.tsv and writes its outputs in the working
    # directory, so run it inside final_results/ (per species in multi-species mode).
    external.log(f"[figures] Rendering ProtParam distributions from {results.name}...")
    external.run([cfg.RSCRIPT_BIN, "--vanilla", str(_R_SCRIPT)], cwd=cfg.final_dir)
    external.log("[OK] ProtParam figures + stats written to final_results/ "
                 "(protparam_boxplots.*, protparam_stats.csv, protparam_outliers.csv)")
