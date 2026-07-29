#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# logos.py - Per-family sequence-conservation logos with WebLogo (the `weblogo` stage).            #
#                                                                                                  #
# Generates PDF + PNG + logodata (TXT) for each *_aligned.fasta produced by the msa stage.         #
# Requires the optional weblogo dependency (pip install weblogo).                                  #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

from . import external
from .config import Config

_FORMATS = [("pdf", "pdf"), ("png", "png"), ("txt", "logodata")]


def _weblogo(cfg, msa_fasta, out_file, fmt, title) -> bool:
    """Run WebLogo for one alignment/format; return True on success."""
    result = external.run([
        cfg.WEBLOGO_BIN,
        "--fin", str(msa_fasta),
        "--fout", str(out_file),
        "--format", fmt,
        "--sequence-type", "protein",
        "--units", "bits",
        "--color-scheme", "chemistry",
        "--resolution", "600",
        "--size", "large",
        "--title", title,
        "--fineprint", "",
        "--errorbars", "YES",
    ], check=False)
    return result.returncode == 0


def run(cfg: Config) -> None:
    """Draw a conservation logo per family from its MSA."""
    cfg.ensure_dirs()
    external.require(cfg.WEBLOGO_BIN)
    msa_dir = cfg.result("msa")
    logo_dir = cfg.result("weblogo")
    logo_dir.mkdir(parents=True, exist_ok=True)

    # Domain alignments only -- a conservation logo is a per-domain figure, so the
    # mature full-sequence alignments are not logo'd.
    msa_files = [f for f in sorted(msa_dir.glob("*_aligned.fasta"))
                 if not f.name.endswith("_mature_aligned.fasta")]
    if not msa_files:
        external.log(f"[WARN] No aligned FASTA in {msa_dir}; run the msa stage first.")
        return

    external.log(f"[OK] Found {len(msa_files)} MSA file(s)")
    for msa_fasta in msa_files:
        family = msa_fasta.name[: -len("_aligned.fasta")]
        title = f"{family} domain conservation"
        external.log(f"Generating logos for {family}...")

        made = [
            label.upper()
            for label, fmt in _FORMATS
            if _weblogo(cfg, msa_fasta, logo_dir / f"{family}_logo.{label}", fmt, title)
        ]
        if made:
            external.log(f"[OK] {family}: {' + '.join(made)} generated")
        else:
            external.log(f"[WARN] {family}: all formats failed")

    external.log("\n[OK] WebLogo step complete. Logos in: intermediate/weblogo/")
