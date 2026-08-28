#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# provenance.py - Write run provenance (the `provenance` stage).                                   #
#                                                                                                  #
# Records the GWIscan version, the full effective parameter set, external tool versions, and the   #
# SHA-256 checksums of the input proteome and family table, to final_results/provenance.txt.       #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

import hashlib
import shutil
import subprocess
from datetime import datetime
from pathlib import Path

from . import __version__, external
from .config import DEFAULTS, Config


def _available(tool: str) -> bool:
    return shutil.which(tool) is not None or Path(tool).exists()


def _tool_version(cmd) -> str:
    """First non-empty line of a binary's version/help output, or 'unknown'."""
    try:
        out = subprocess.run([str(c) for c in cmd], stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, text=True, timeout=30)
        for line in out.stdout.splitlines():
            if line.strip():
                return line.strip()
    except Exception:  # noqa: BLE001
        pass
    return "unknown"


def _sha256(path) -> str:
    """SHA-256 of a file, streamed in chunks."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run(cfg: Config) -> None:
    """Write final_results/provenance.txt for the run."""
    cfg.final_dir.mkdir(parents=True, exist_ok=True)
    out = cfg.final_dir / "provenance.txt"

    lines = [
        "GWIscan run provenance",
        f"Generated: {datetime.now()}",
        f"GWIscan version: {__version__}",
        "",
        "--- Parameters (effective config) ---",
    ]
    for key in DEFAULTS:
        lines.append(f"{key}={getattr(cfg, key)}")

    lines += ["", "--- Tool versions ---"]
    probes = [
        ("hmmsearch", ["hmmsearch", "-h"]),
        ("diamond", ["diamond", "version"]),
        ("seqkit", ["seqkit", "version"]),
        ("mafft", ["mafft", "--version"]),
        ("weblogo", ["weblogo", "--version"]),
        ("meme", ["meme", "-version"]),
        ("iqtree", [cfg.IQTREE_BIN, "--version"]),
        ("targetp", [cfg.TARGETP_BIN, "-h"]),
        ("deeploc", [cfg.DEEPLOC_BIN, "--help"]),
        ("biolib", ["biolib", "--version"]),
    ]
    if str(cfg.INTERPRO_MODE).lower() == "local":
        probes.append(("interproscan", [cfg.INTERPROSCAN_BIN, "--version"]))
    for label, cmd in probes:
        if _available(str(cmd[0])):
            lines.append(f"{label}: {_tool_version(cmd)}")
    if str(cfg.INTERPRO_MODE).lower() == "api":
        lines.append(f"interproscan: EBI API (appl={cfg.INTERPRO_APPL})")

    lines += ["", "--- Input checksums ---"]
    for path in (cfg.proteome, cfg.family_map):
        if path.exists():
            lines.append(f"sha256  {path.name}  {_sha256(path)}")

    out.write_text("\n".join(lines) + "\n")
    external.log(f"[OK] Provenance written: {out}")
