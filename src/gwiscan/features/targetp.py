#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# targetp.py - N-terminal presequence prediction with TargetP 2.0 (the `targetp` stage).           #
#                                                                                                  #
# TargetP 2.0 is a user-installed academic download (not redistributable); point at it via         #
# TARGETP_BIN (PATH name, executable path, or install dir -> <dir>/bin/targetp; see                #
# resolve_targetp). Runs once with -org pl -format short -gff3 -mature. Keeps the full short-       #
# summary (type + all class probabilities + cleavage site), columns named from the summary's own   #
# header (plant + non-plant layouts); ID -> protein_id. Raw GFF3 + mature FASTA go to              #
# intermediate/targetp_raw/.                                                                             #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

import shutil
from pathlib import Path

from .. import external, io
from ..config import Config

SUMMARY_SUFFIX = "_summary.targetp2"


def resolve_targetp(targetp_bin) -> str:
    """Resolve the TargetP executable.

    Accepts a PATH name (``targetp``), a direct path to the executable, or the
    unpacked install directory. TargetP always keeps the executable at
    ``<dir>/bin/targetp`` (with ``<dir>/lib`` beside it so it finds its
    libraries), so if a directory (or any path) has ``bin/targetp`` under it,
    that is used.
    """
    nested = Path(targetp_bin) / "bin" / "targetp"
    if nested.exists():
        return str(nested)
    return str(targetp_bin)

# Summary header names -> pipeline column names (others kept verbatim, e.g. the
# standard TargetP class columns noTP / SP / mTP / cTP / luTP).
_RENAME = {"ID": "protein_id", "Prediction": "targetp_type", "CS Position": "cs_position"}


def _split(line: str) -> list:
    """Tab-split when tabs are present (keeps the space-containing CS field
    intact), else fall back to whitespace."""
    return line.split("\t") if "\t" in line else line.split()


def parse_summary(path):
    """Parse a TargetP 2.0 short summary into ``(header, rows)``.

    Reads the ``# ID Prediction ...`` comment as the column header and keeps all
    data columns. Returns the normalised header and one row per protein.
    """
    header = None
    data = []
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith("#"):
                fields = _split(line.lstrip("#").strip())
                if header is None and fields[:2] == ["ID", "Prediction"]:
                    header = [_RENAME.get(f, f) for f in fields]
                continue
            data.append(_split(line))

    if header is None:
        # No recognisable header — fall back to generic names so nothing is lost.
        width = max((len(r) for r in data), default=2)
        header = ["protein_id", "targetp_type"] + [f"col_{i}" for i in range(3, width + 1)]

    n = len(header)
    rows = [(r + [""] * n)[:n] for r in data]   # normalise each row to header width
    return header, rows


def run(cfg: Config) -> None:
    cfg.ensure_dirs()
    cand_fasta = cfg.result("final_candidates.fasta")
    out_dir = cfg.result("targetp_raw")
    out_tsv = cfg.result("targetp.tsv")
    out_dir.mkdir(parents=True, exist_ok=True)

    if not cand_fasta.exists():
        raise FileNotFoundError(
            f"final_candidates.fasta not found: {cand_fasta} (run `gwiscan confirm` first)"
        )

    # targetp_bin may be a PATH name, a path to the executable, or the install
    # directory (which contains bin/targetp).
    targetp = resolve_targetp(cfg.TARGETP_BIN)
    if shutil.which(targetp) is None and not Path(targetp).exists():
        raise FileNotFoundError(
            f"TargetP executable not found: {cfg.TARGETP_BIN!r}. TargetP 2.0 requires a "
            f"free academic license — download and install it, then set targetp_bin / "
            f"TARGETP_BIN to the 'targetp' executable or its install directory (which "
            f"contains bin/targetp). https://services.healthtech.dtu.dk/services/TargetP-2.0/"
        )

    prefix = out_dir / "targetp"
    external.log(
        f"[targetp] Running TargetP 2.0 (-org {cfg.TARGETP_ORGANISM} -format short, "
        f"-gff3 -mature)..."
    )
    external.run([
        targetp,
        "-fasta", cand_fasta.resolve(),
        "-org", cfg.TARGETP_ORGANISM,
        "-format", "short",
        "-prefix", prefix,
        "-gff3",      # feature file of processed sequences
        "-mature",    # FASTA of mature sequences (presequence cleaved off)
    ])

    summary = out_dir / f"targetp{SUMMARY_SUFFIX}"
    if not summary.exists():
        matches = sorted(out_dir.glob(f"*{SUMMARY_SUFFIX}"))
        if not matches:
            raise FileNotFoundError(
                f"TargetP produced no *{SUMMARY_SUFFIX} file in {out_dir} "
                f"(check the TargetP output above)"
            )
        summary = matches[0]

    header, rows = parse_summary(summary)
    io.write_tsv(out_tsv, header, rows)
    external.log(f"[OK] TargetP: {len(rows)} entries -> targetp.tsv")
