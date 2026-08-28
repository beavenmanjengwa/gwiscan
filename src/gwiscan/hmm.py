#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# hmm.py - hmmsearch against the HMM database and domtbl parsing (the `search-hmm` stage).          #
#                                                                                                  #
# Uses hmmsearch (the profiles are the query, the proteome the target DB): far faster than hmmscan  #
# for a genome-wide search and it saturates the cores. In hmmsearch --domtblout the target is the   #
# protein and the query is the HMM, so column 1 is the protein id and columns 4-5 the HMM name/acc. #
# The column indices below are named accordingly and covered by tests/test_hmm_parse.py.            #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

import re
from datetime import datetime

from . import external, io
from .config import Config
from .schema import HIT_HEADER

# A GA (gathering) threshold line in an HMMER profile header, e.g. "GA  27.0 27.0;".
_GA_RE = re.compile(r"^GA\s")

# Guidance shown when a custom identifying HMM lacks GA thresholds. Kept in one
# place so preflight and setup-db give the user the same, correct instructions.
GA_THRESHOLD_HELP = (
    "The HMM search runs `hmmsearch --cut_ga`, which applies each model's own GA "
    "(gathering) cutoff, so every identifying model must declare one. A profile from "
    "hmmbuild has no GA line unless its source alignment carried one. Add a GA cutoff "
    "to the profile (build it from a Stockholm alignment with a `#=GF GA <seq> <dom>` "
    "line, or insert a `GA  <seq> <dom>;` line into the .hmm header), or identify this "
    "family by a Pfam accession instead."
)


def has_ga_thresholds(hmm_path) -> bool:
    """True when every model in an HMMER profile file declares a GA threshold.

    An .hmm file may hold several concatenated models; ``hmmsearch --cut_ga`` aborts
    the whole search if any one of them lacks a GA cutoff, so this requires one GA
    line per model (and at least one model)."""
    n_models = n_ga = 0
    with open(hmm_path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("HMMER3/"):     # each model begins with the format line
                n_models += 1
            elif _GA_RE.match(line):
                n_ga += 1
    return n_models >= 1 and n_ga == n_models


def custom_hmm_ga_error(hmm_path, family: str) -> str | None:
    """An explanatory error string if a custom identifying HMM lacks GA thresholds,
    else None. Shared by preflight (logs it) and setup-db (raises it)."""
    if has_ga_thresholds(hmm_path):
        return None
    return f"custom HMM for '{family}' has no GA (gathering) thresholds: {hmm_path}. {GA_THRESHOLD_HELP}"


# hmmsearch --domtblout fixed columns (0-based); see HMMER Userguide. hmmsearch's
# target is the protein and its query is the HMM (the reverse of hmmscan).
_COL = {
    "protein": 0,       # target NAME: protein id from proteome.fasta
    "hmm_name": 3,      # query NAME:  HMM model name (e.g. Lectin_legB)
    "hmm_acc": 4,       # query ACC:   HMM accession (e.g. PF00139.27)
    "i_evalue": 12,     # domain independent E-value
    "dom_score": 13,    # domain bit score
    "ali_from": 17,     # alignment start on the protein (target sequence)
    "ali_to": 18,       # alignment end on the protein
}
_MIN_COLS = 22


def run(cfg: Config) -> int:
    """Run hmmsearch on the proteome and parse the hits."""
    cfg.ensure_dirs()
    if not cfg.proteome.exists():
        raise FileNotFoundError(f"proteome not found: {cfg.proteome}")
    if not cfg.hmm_db.exists():
        raise FileNotFoundError(
            f"HMM database not found: {cfg.hmm_db} (run `gwiscan setup-db` first)"
        )
    external.require("hmmsearch")

    domtbl = cfg.result("hmm_domtbl.txt")
    external.log(f"[{datetime.now()}] Running hmmsearch with --cut_ga...")
    external.run([
        "hmmsearch",
        "--cut_ga",
        "--noali",                              # the domtbl is parsed; skip the big alignment dump
        "--cpu", cfg.THREADS,
        "--domtblout", domtbl,
        "--tblout", cfg.result("hmm_tbl.txt"),
        "-o", cfg.result("hmm_full.out"),
        cfg.hmm_db,
        cfg.proteome,
    ])
    external.log(f"[{datetime.now()}] hmmsearch complete.")

    n = parse_domtbl(cfg, domtbl)
    external.log(f"[{datetime.now()}] hmmsearch step done.")
    return n


def parse_domtbl(cfg: Config, domtbl=None) -> int:
    """Parse an hmmsearch --domtblout file into intermediate/hmm_hits.tsv."""
    domtbl = domtbl or cfg.result("hmm_domtbl.txt")
    if not domtbl.exists():
        raise FileNotFoundError(f"domtblout not found: {domtbl} (run hmmsearch first)")

    # Family/multi-family modes map each Pfam to its curated family. Architecture
    # mode has no per-domain family — the label IS the HMM's own domain name (e.g.
    # B_lectin, Pkinase), so the map is empty and every hit uses that name (the
    # fallback below). architecture.run then combines the domains per protein.
    pfam_to_family = {} if cfg.is_architecture else io.pfam_to_family(cfg.family_map)
    rows, unmapped, skipped = [], set(), 0

    with open(domtbl) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split()
            if len(cols) < _MIN_COLS:                    # short/truncated row (e.g. a killed job)
                skipped += 1
                continue

            hmm_acc = cols[_COL["hmm_acc"]]              # e.g. PF00139.27, or '-' for a custom HMM
            pfam_base = hmm_acc.split(".")[0]            # e.g. PF00139
            family = pfam_to_family.get(pfam_base)
            if family is None:
                family = cols[_COL["hmm_name"]]          # fall back to raw HMM name
                # A custom HMM legitimately has no accession ('-'); its name IS the
                # family, so that is not an "unmapped Pfam" to warn about.
                if pfam_base != "-":
                    unmapped.add(pfam_base)

            rows.append([
                cols[_COL["protein"]],
                family,
                hmm_acc,
                cols[_COL["i_evalue"]],
                cols[_COL["dom_score"]],
                cols[_COL["ali_from"]],
                cols[_COL["ali_to"]],
                "hmm",
            ])

    io.write_tsv(cfg.result("hmm_hits.tsv"), HIT_HEADER, rows)
    n_prot = len({r[0] for r in rows})
    n_fam = len({r[1] for r in rows})
    external.log(
        f"[OK] Parsed {len(rows)} domain hits "
        f"({n_prot} proteins, {n_fam} families) -> hmm_hits.tsv"
    )
    if skipped:
        external.log(
            f"[WARN] {skipped} malformed/short domtbl line(s) skipped -- the domtbl "
            f"may be truncated (a killed hmmsearch?); results could be incomplete."
        )
    if unmapped and not cfg.is_architecture:
        external.log(
            f"[WARN] {len(unmapped)} Pfam accession(s) not in the family table, "
            f"used raw HMM name instead: {', '.join(sorted(unmapped))}"
        )
    return len(rows)
