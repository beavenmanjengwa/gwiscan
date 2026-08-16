#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# preflight.py - Verify required tools, Python packages, and inputs (the `preflight` stage).       #
#                                                                                                  #
# Fails fast with a clear message rather than letting a long run die halfway on a missing tool.    #
# Checks the InterProScan requirement for the chosen mode (EBI_EMAIL for api, the binary local).   #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

import importlib
import shutil
from pathlib import Path

from . import external, hmm, io
from .config import Config

# (binary, version command) for tools that must be present.
_REQUIRED_BINS = [
    ("hmmscan", ["hmmscan", "-h"]),
    ("hmmpress", ["hmmpress", "-h"]),
    ("diamond", ["diamond", "version"]),
    ("seqkit", ["seqkit", "version"]),
    ("mafft", ["mafft", "--version"]),
]

_REQUIRED_PKGS = ["Bio", "pandas", "openpyxl", "requests", "yaml"]


def _check_family_reference_files(cfg: Config) -> bool:
    """Verify every reference file the family table names is present. Returns True
    if anything is missing (so the caller can fail preflight).

    For each family: its BlastModel FASTA (db/blast/<name>.fasta -- required for
    all, since DIAMOND runs for every family) and, for a family identified by a
    user-provided custom HMM (not a downloadable Pfam accession), that HMM in
    db/hmm/<name>.hmm. Pfam-accession HMMs are NOT checked here -- setup-db fetches
    those from InterPro on demand, so their absence is expected and fine."""
    try:
        records = io.family_records(cfg.family_map)
    except Exception as e:  # noqa: BLE001 - a malformed table should fail preflight, not crash it
        external.log(f"[MISSING] family table could not be parsed: {e}")
        return True

    missing = False
    for r in records:
        family = r["family"]
        blast_model = r["blast_model"]
        if not blast_model:
            external.log(f"[MISSING] family '{family}' has no BlastModel in the family table "
                         f"(DIAMOND runs for every family, so it is required)")
            missing = True
        elif not (cfg.blast_dir / blast_model).exists():
            external.log(f"[MISSING] BlastModel FASTA for '{family}' not found: "
                         f"{cfg.blast_dir / blast_model}")
            missing = True

        if r["hmm_is_custom"]:
            custom_hmm = cfg.hmm_dir / r["hmm_file"]
            if not custom_hmm.exists():
                external.log(f"[MISSING] custom HMM for '{family}' not found: {custom_hmm} "
                             f"(build it and place it in db/hmm/)")
                missing = True
            elif r["hmm_press"]:
                # An identifying custom HMM must carry GA thresholds, or hmmscan
                # --cut_ga aborts the whole search once it is pressed in. Catch it
                # here rather than mid-run with HMMER's terse "GA unavailable" error.
                ga_error = hmm.custom_hmm_ga_error(custom_hmm, family)
                if ga_error:
                    external.log(f"[MISSING] {ga_error}")
                    missing = True

    if not missing:
        external.log(f"[OK] family reference files present ({len(records)} families)")
    return missing


def _duplicate_fasta_ids(path) -> list:
    """Sequence ids (the first header token) that appear more than once in a FASTA.

    merge/confirm/architecture index the proteome with SeqIO.to_dict, which raises
    on a duplicate id partway through the run with a terse Biopython message;
    preflight uses this to surface the offending ids up front instead."""
    from collections import Counter
    ids = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                parts = line[1:].split()
                ids.append(parts[0] if parts else "")
    return sorted(i for i, count in Counter(ids).items() if count > 1)


def _check_one_proteome(path: Path, label: str) -> bool:
    """Validate a single proteome FASTA: exists, has records, unique ids. Returns
    True on any problem, logging a labelled message for each."""
    if not path.exists():
        external.log(f"[MISSING] {label} proteome not found: {path}")
        return True
    with open(path, encoding="utf-8") as handle:
        n = sum(1 for line in handle if line.startswith(">"))
    if n == 0:
        external.log(f"[ERROR] {label} proteome contains no FASTA records: {path}")
        return True
    external.log(f"[OK] {label} proteome: {n} sequences")
    dup = _duplicate_fasta_ids(path)
    if dup:
        external.log(f"[MISSING] {label} proteome has {len(dup)} duplicate sequence id(s) "
                     f"(e.g. {', '.join(dup[:5])}); ids must be unique because the pipeline "
                     f"indexes the proteome by id. Deduplicate the FASTA.")
        return True
    return False


def _check_proteomes(cfg: Config) -> bool:
    """Validate the run's input proteome(s). Returns True on any problem.

    Multi-species runs (config/species.tsv present, and not already inside a single
    species' sub-run) have no single input/proteome.fasta -- each proteome comes from
    the manifest -- so every manifest species' proteome is validated instead, plus a
    warning for any missing annotation. Otherwise the single input proteome is checked."""
    if cfg.species_manifest.exists() and not cfg.SPECIES:
        from . import species as species_mod
        try:
            manifest = species_mod.read_species_manifest(cfg.species_manifest)
        except Exception as e:  # noqa: BLE001 - a bad manifest should fail preflight cleanly
            external.log(f"[MISSING] species manifest could not be parsed: {e}")
            return True
        bad = False
        for entry in manifest:
            proteome = Path(entry["proteome"])
            proteome = proteome if proteome.is_absolute() else (cfg.root / proteome)
            bad = _check_one_proteome(proteome, entry["prefix"]) or bad
            annotation = entry.get("annotation")
            if annotation:
                ann = Path(annotation)
                ann = ann if ann.is_absolute() else (cfg.root / ann)
                if not ann.exists():
                    external.log(f"[WARN] {entry['prefix']}: annotation not found: {ann} "
                                 f"(this species runs without chromosomal coordinates)")
        if not bad:
            external.log(f"[OK] {len(manifest)} species proteome(s) present "
                         f"({', '.join(e['prefix'] for e in manifest)})")
        return bad
    return _check_one_proteome(cfg.proteome, "input")


def _check_config_values(cfg: Config) -> bool:
    """Sanity-check numeric knobs so a typo (identity 300, evalue "1e-5x") fails at
    preflight with a clear message, not hours later inside a tool invocation with
    whatever cryptic error that tool emits. Returns True on any problem.

    Each entry: (label, predicate over the value, rule text). A value that fails
    its caster (e.g. non-numeric threads) is caught and reported the same way."""
    checks = [
        ("THREADS", cfg.THREADS, lambda v: int(v) >= 1, "must be an integer >= 1"),
        ("DIAMOND_IDENTITY", cfg.DIAMOND_IDENTITY, lambda v: 0 <= float(v) <= 100,
         "must be a percentage between 0 and 100"),
        ("DIAMOND_COVERAGE_R2", cfg.DIAMOND_COVERAGE_R2, lambda v: 0 <= float(v) <= 100,
         "must be a percentage between 0 and 100"),
        ("DIAMOND_EVALUE", cfg.DIAMOND_EVALUE, lambda v: float(v) >= 0,
         'must be a non-negative number (e.g. "1e-5")'),
        ("DIAMOND_BSR", cfg.DIAMOND_BSR, lambda v: 0 <= float(v) <= 1,
         "must be a Blast Score Ratio between 0 and 1"),
        ("CONCORDANCE_MIN", cfg.CONCORDANCE_MIN, lambda v: 0 <= float(v) <= 1,
         "must be a Jaccard fraction between 0 and 1"),
        ("IQTREE_BOOTSTRAP", cfg.IQTREE_BOOTSTRAP, lambda v: int(v) >= 0,
         "must be an integer >= 0 (0 = no bootstrap)"),
        ("MEME_NMOTIFS", cfg.MEME_NMOTIFS, lambda v: int(v) >= 1,
         "must be an integer >= 1"),
        ("SPECIES_PARALLEL", cfg.SPECIES_PARALLEL, lambda v: int(v) >= 0,
         "must be an integer >= 0 (0 = auto)"),
        ("DEEPTMHMM_MODE", cfg.DEEPTMHMM_MODE,
         lambda v: str(v).strip().lower() in {"biolib", "local"},
         "must be 'biolib' or 'local'"),
    ]
    bad = False
    for label, value, ok, rule in checks:
        try:
            valid = ok(value)
        except (TypeError, ValueError):
            valid = False
        if not valid:
            external.log(f"[MISSING] {label}={value!r} is invalid — {rule}")
            bad = True
    if not bad:
        external.log("[OK] config parameters within valid ranges")
    return bad


def run(cfg: Config) -> None:
    """Executes environment diagnostic validations to assert system and workflow integrity."""
    cfg.ensure_dirs()
    # Validate before either mode-specific pre-flight path: architecture has its
    # own shorter pre-flight and would otherwise not reach _check_config_values.
    if str(cfg.DEEPTMHMM_MODE).strip().lower() not in {"biolib", "local"}:
        raise RuntimeError(
            f"invalid DEEPTMHMM_MODE {cfg.DEEPTMHMM_MODE!r}; use 'biolib' or 'local'"
        )
    # Architecture mode is hmmscan-only with its own inputs (the rules table, no
    # DIAMOND / InterProScan); it has its own lighter pre-flight.
    if cfg.is_architecture:
        from . import architecture
        architecture.preflight(cfg)
        return
    fail = False
    external.log("[preflight] GWIscan pre-flight check")
    external.log("--------------------------------------")

    for binary, _ in _REQUIRED_BINS:
        if shutil.which(binary) is None:
            external.log(f"[MISSING] {binary} not found on PATH")
            fail = True
        else:
            external.log(f"[OK] {binary}")

    # Optional tools: warn only (skippable via flags / separate steps). Each may be
    # a PATH name or an absolute path (weblogo/meme/iqtree/targetp/deeploc bins), so
    # external.available accepts either.
    for opt in (cfg.TARGETP_BIN, cfg.DEEPLOC_BIN, "biolib",
                cfg.WEBLOGO_BIN, cfg.MEME_BIN, cfg.CLIPKIT_BIN, cfg.IQTREE_BIN,
                cfg.RSCRIPT_BIN):
        if external.available(opt):
            external.log(f"[OK] {opt} found")
        else:
            external.log(f"[WARN] optional tool not found: {opt}")

    for pkg in _REQUIRED_PKGS:
        try:
            importlib.import_module(pkg)
            external.log(f"[OK] python package: {pkg}")
        except ImportError:
            external.log(f"[MISSING] python package not importable: {pkg}")
            fail = True

    # In multi-species mode the proteomes come from config/species.tsv, not the
    # single-species input/proteome.fasta; _check_proteomes handles both.
    fail = _check_proteomes(cfg) or fail

    if not cfg.family_map.exists():
        external.log(f"[MISSING] family table not found: {cfg.family_map} "
                     f"(expected for MODE={cfg.MODE})")
        fail = True
    else:
        external.log(f"[OK] family table present: {cfg.family_map.name}")
        # The reference files the family table POINTS AT must exist too, or the run
        # dies later inside setup-db (after hmmpress has already run) on a missing
        # BlastModel FASTA or custom HMM. Check them here so a typo fails fast.
        fail = _check_family_reference_files(cfg) or fail

    # InterProScan is mandatory; check the requirement for the chosen mode.
    ipr_mode = str(cfg.INTERPRO_MODE).lower()
    if ipr_mode == "local":
        if external.available(cfg.INTERPROSCAN_BIN):
            external.log(f"[OK] InterProScan (local): {cfg.INTERPROSCAN_BIN}")
        else:
            external.log(f"[MISSING] local InterProScan not found: {cfg.INTERPROSCAN_BIN}")
            fail = True
    elif ipr_mode == "api":
        if cfg.EBI_EMAIL:
            external.log("[OK] InterProScan (api): EBI_EMAIL set")
        else:
            external.log("[MISSING] INTERPRO_MODE=api needs EBI_EMAIL "
                         "(set it, or use INTERPRO_MODE=local)")
            fail = True
    else:
        external.log(f"[MISSING] unknown INTERPRO_MODE {cfg.INTERPRO_MODE!r} (use 'api' or 'local')")
        fail = True

    # Numeric parameter ranges (identity/coverage/e-value/BSR/... ), so a typo
    # fails here rather than deep inside a tool call.
    fail = _check_config_values(cfg) or fail

    # Isoform advisory: the pipeline never collapses isoforms — the input must be
    # one protein per gene, asserted via primary_transcript.
    if cfg.PRIMARY_TRANSCRIPT:
        external.log("[OK] primary_transcript=true — input treated as one protein per gene")
    else:
        external.log(
            "[WARN] primary_transcript=false — if the proteome contains splice isoforms, "
            "family counts and domain trees will be inflated. Supply a primary-transcript "
            "proteome (one protein per gene) and set primary_transcript=true."
        )

    if fail:
        raise RuntimeError("pre-flight check failed — see messages above")
    external.log("[PASS] Pre-flight check passed.")
