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


def tool_versions(cfg: Config) -> list:
    """(tool, version) for every tool the pipeline can use, in run order.

    The single source of truth for tool versioning, shared by provenance.txt and
    the results workbook's Tool_versions sheet (compile). A tool absent from this
    environment is reported as '(not found)' rather than omitted, so the list is
    always complete — every tool carries a version. Tools that live in their own
    conda env are probed via cfg's configured executable (an absolute path in the
    config still resolves); the ones without a --version flag (DeepTMHMM's pinned
    app tag, the EBI InterProScan release) are read from config / the run manifest."""
    rows = [("GWIscan", __version__)]

    # External binaries, probed with the configured executable.
    probes = [
        ("hmmsearch", ["hmmsearch", "-h"]),
        ("diamond", ["diamond", "version"]),
        ("seqkit", ["seqkit", "version"]),
        ("mafft", ["mafft", "--version"]),
        ("targetp", [cfg.TARGETP_BIN, "-h"]),
        ("deeploc", [cfg.DEEPLOC_BIN, "--help"]),
        ("clipkit", [cfg.CLIPKIT_BIN, "--version"]),
        ("weblogo", [cfg.WEBLOGO_BIN, "--version"]),
        ("meme", [cfg.MEME_BIN, "-version"]),
        ("iqtree", [cfg.IQTREE_BIN, "--version"]),
        ("Rscript", [cfg.RSCRIPT_BIN, "--version"]),
    ]
    for label, cmd in probes:
        rows.append((label, _tool_version(cmd) if _available(str(cmd[0])) else "(not found)"))

    # Python packages that do annotation work directly (no external binary):
    # Biopython (ProtParam physicochemistry + FASTA parsing) and pandas (tables).
    for label, module in (("biopython (ProtParam)", "Bio"), ("pandas", "pandas")):
        try:
            rows.append((label, getattr(__import__(module), "__version__", "unknown")))
        except Exception:  # noqa: BLE001
            rows.append((label, "(not installed)"))

    # DeepTMHMM: the predictor version is the pinned app tag (DEEPTMHMM_VERSION),
    # not a binary --version. Record it plus how it runs.
    if str(cfg.DEEPTMHMM_MODE).lower() == "local":
        rows.append(("deeptmhmm", f"DeepTMHMM (local install, DEEPTMHMM_DIR="
                                  f"{cfg.DEEPTMHMM_DIR or 'unset'})"))
    else:
        biolib_v = _tool_version(["biolib", "--version"]) if _available("biolib") else "(not found)"
        rows.append(("deeptmhmm", f"DeepTMHMM {cfg.DEEPTMHMM_VERSION} (biolib CLI: {biolib_v})"))

    # InterProScan: a local install's --version, or the EBI REST service with the
    # authoritative InterProScan/InterPro release captured in the run manifest.
    if str(cfg.INTERPRO_MODE).lower() == "local":
        v = _tool_version([cfg.INTERPROSCAN_BIN, "--version"]) if _available(cfg.INTERPROSCAN_BIN) else "(not found)"
        rows.append(("interproscan", f"local install {v}"))
    else:
        from .features import interpro
        try:
            version = interpro._api_version(cfg)
            apps = interpro._appl_for_version(interpro._run_applications(cfg), version)
            endpoint = interpro.IPRSCAN_ENDPOINTS[version]
            rows.append(("interproscan", f"EBI InterProScan {version} API "
                                         f"(endpoint={endpoint}, appl={','.join(apps)})"))
        except Exception as e:  # noqa: BLE001
            rows.append(("interproscan", f"EBI API (version resolution failed: {e})"))
        manifest = cfg.result("interproscan.manifest.txt")
        if manifest.exists():
            for ln in manifest.read_text().splitlines():
                if ln.startswith("interproscan-version:"):
                    rows.append(("interproscan release", ln.split(":", 1)[1].strip()))
    return rows


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
    for tool, version in tool_versions(cfg):
        lines.append(f"{tool}: {version}")

    lines += ["", "--- Input checksums ---"]
    for path in (cfg.proteome, cfg.family_map):
        if path.exists():
            lines.append(f"sha256  {path.name}  {_sha256(path)}")

    out.write_text("\n".join(lines) + "\n")
    external.log(f"[OK] Provenance written: {out}")
