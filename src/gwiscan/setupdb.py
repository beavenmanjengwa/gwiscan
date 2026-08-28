#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# setupdb.py - Build the HMM database and the DIAMOND proteome database (`setup-db`).              #
#                                                                                                  #
# setup_shared: fetch each family's Pfam HMM from the InterPro API (or reuse a provided/custom     #
# db/hmm/<name>.hmm), concatenate into one .hmm (hmmsearch reads it directly, no hmmpress), and    #
# verify every family's BlastModel FASTA.                                                          #
# setup_proteome: build the DIAMOND database of this proteome (both DIAMOND rounds search it).     #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

import gzip
from datetime import datetime

from . import external, hmm, io, net
from .config import Config

# One HMM per Pfam accession (gzipped) from the InterPro API, not the full
# Pfam-A.hmm. A locally provided db/hmm/<acc>.hmm is used instead if present.
_HMM_URL = "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{acc}?annotation=hmm"


def _ensure_hmm(cfg: Config, acc: str):
    """Return db/hmm/<acc>.hmm, downloading + gunzipping it if not already present."""
    dest = cfg.hmm_dir / f"{acc}.hmm"
    if dest.exists():
        external.log(f"[OK] HMM present (provided): {dest.name}")
        return dest
    url = _HMM_URL.format(acc=acc)
    external.log(f"[..] Fetching HMM {acc} from InterPro: {url}")
    try:
        # net.fetch adds a User-Agent and retries transient failures with backoff
        # (same hardening as the go-basic.obo download).
        data = net.fetch(url, timeout=60)
    except Exception as e:  # noqa: BLE001 - surface any network/HTTP error clearly
        raise RuntimeError(f"failed to download HMM for {acc} from {url}: {e}")
    if not data:
        raise RuntimeError(f"empty HMM download for {acc} from {url}")
    try:
        data = gzip.decompress(data)          # API serves gzip; tolerate plain too
    except (OSError, gzip.BadGzipFile):
        pass
    dest.write_bytes(data)
    external.log(f"[OK] Downloaded HMM: {dest.name}")
    return dest


def _build_hmm_db(cfg: Config, hmm_models, db_path) -> None:
    """Concatenate the collected HMM files into db_path. hmmsearch reads this .hmm
    directly, so no hmmpress step is needed."""
    if not hmm_models:
        external.log(f"[WARN] no models for {db_path.name}; HMM search will be empty")
    external.log(f"[{datetime.now()}] Building HMM database {db_path.name} ({len(hmm_models)} models)...")
    with open(db_path, "wb") as out:
        for hmm_path in hmm_models:
            out.write(hmm_path.read_bytes())
    external.log(f"[OK] HMM database built: {db_path.name}")


def setup_architecture(cfg: Config) -> None:
    """Architecture-mode setup: fetch and build two Pfam HMM databases from the rules.
    A PRIMARY-only db seeds the genome-wide search (pass 1); a primary+required db
    searches the candidates (pass 2). No DIAMOND / BLAST -- this mode is hmmsearch-only
    and InterProScan annotates the final candidates afterwards."""
    from . import architecture

    cfg.hmm_dir.mkdir(parents=True, exist_ok=True)
    rules = architecture.read_rules(cfg.architecture_map)
    primary = architecture.primary_accessions(rules)
    every = architecture.all_accessions(rules)
    external.log(f"[{datetime.now()}] Architecture mode: {len(rules)} rule(s); "
                 f"primary {primary}; primary+required {every}")

    _build_hmm_db(cfg, [_ensure_hmm(cfg, acc) for acc in primary], cfg.primary_hmm_db)
    _build_hmm_db(cfg, [_ensure_hmm(cfg, acc) for acc in every], cfg.hmm_db)


def setup_shared(cfg: Config) -> None:
    """Species-independent setup: build the identifying HMM db and verify every
    family's BlastModel query FASTA. Safe to run once and reuse across species
    (the HMM db and model FASTAs are shared, read-only references).
    """
    if cfg.is_architecture:
        setup_architecture(cfg)
        return

    cfg.hmm_dir.mkdir(parents=True, exist_ok=True)
    cfg.blast_dir.mkdir(parents=True, exist_ok=True)
    external.require("diamond")

    records = io.family_records(cfg.family_map)

    # Collect the identifying HMM models (hmm_press families only).
    # Pfam accession -> downloaded (cached). Custom HMM -> must be provided in db/hmm/.
    hmm_models = []
    for r in records:
        if not r["hmm_press"]:
            continue
        if r["hmm_is_custom"]:
            dest = cfg.hmm_dir / r["hmm_file"]
            if not dest.exists():
                raise FileNotFoundError(
                    f"custom HMM for family '{r['family']}' not found: {dest} "
                    f"(build it and place it in db/hmm/)"
                )
            # A custom HMM without GA thresholds makes the later `hmmsearch --cut_ga`
            # abort the whole search, so refuse it now with the same guidance
            # preflight gives.
            ga_error = hmm.custom_hmm_ga_error(dest, r["family"])
            if ga_error:
                raise ValueError(ga_error)
            external.log(f"[OK] custom HMM present: {dest.name}")
        else:
            dest = _ensure_hmm(cfg, r["pfam_model"])
        hmm_models.append(dest)

    if not hmm_models:
        external.log("[WARN] no families have an identifying HMM; HMM search will be empty")

    external.log(f"[{datetime.now()}] Building HMM database ({len(hmm_models)} models)...")
    with open(cfg.hmm_db, "wb") as out:
        for hmm_path in hmm_models:
            out.write(hmm_path.read_bytes())
    external.log(f"[OK] HMM database built: {cfg.hmm_db.name}")

    # Every family is DIAMOND-searched, so every BlastModel query must exist.
    for r in records:
        if not r["blast_model"]:
            raise FileNotFoundError(f"family '{r['family']}' has no BlastModel in the family table")
        model_fasta = cfg.blast_dir / r["blast_model"]
        if not model_fasta.exists():
            raise FileNotFoundError(
                f"Missing BLAST model FASTA (round-1 query) for {r['family']}: {model_fasta}"
            )
        external.log(f"[OK] blast model query present: {r['blast_model']}")


def setup_proteome(cfg: Config) -> None:
    """Per-species setup: build the DIAMOND database of this species' proteome
    (both DIAMOND rounds search it). Writes to cfg.proteome_db, which is
    species-namespaced in multi-species mode.
    """
    cfg.blast_dir.mkdir(parents=True, exist_ok=True)
    external.require("diamond")
    if not cfg.proteome.exists():
        raise FileNotFoundError(f"proteome not found: {cfg.proteome}")
    external.log(f"[{datetime.now()}] Building DIAMOND database of the proteome...")
    external.run([
        "diamond", "makedb",
        "--in", cfg.proteome,
        "--db", cfg.proteome_db,
        "--threads", cfg.THREADS,
    ])
    external.log(f"[OK] DIAMOND proteome db built: {cfg.proteome_db.name}.dmnd")


def run(cfg: Config) -> None:
    """Build the shared HMM database and this proteome's DIAMOND database.

    Architecture mode is hmmsearch-only, so it builds just the component HMM db and
    skips the DIAMOND proteome database entirely."""
    external.log(f"[{datetime.now()}] Starting database setup")
    setup_shared(cfg)
    if not cfg.is_architecture:
        setup_proteome(cfg)
    external.log(f"[{datetime.now()}] Database setup complete.")
