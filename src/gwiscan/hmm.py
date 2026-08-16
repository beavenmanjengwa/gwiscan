#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# hmm.py - hmmscan against the HMM database and domtbl parsing (the `search-hmm` stage).           #
#                                                                                                  #
# In hmmscan --domtblout the target is the HMM model and the query is the protein: column 1 is     #
# the HMM name, column 4 the protein id. The column indices below are named accordingly and        #
# covered by tests/test_hmm_parse.py.                                                              #
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
    "The HMM search runs `hmmscan --cut_ga`, which applies each model's own GA "
    "(gathering) cutoff, so every identifying model must declare one. A profile from "
    "hmmbuild has no GA line unless its source alignment carried one. Add a GA cutoff "
    "to the profile (build it from a Stockholm alignment with a `#=GF GA <seq> <dom>` "
    "line, or insert a `GA  <seq> <dom>;` line into the .hmm header), or identify this "
    "family by a Pfam accession instead."
)


def has_ga_thresholds(hmm_path) -> bool:
    """True when every model in an HMMER profile file declares a GA threshold.

    An .hmm file may hold several concatenated models; ``hmmscan --cut_ga`` aborts
    the whole search if any one of them lacks a GA cutoff, so this requires one GA
    line per model (and at least one model)."""
    n_models = n_ga = 0
    with open(hmm_path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("HMMER3"):      # each model begins with the format line
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


# hmmscan --domtblout fixed columns (0-based); see HMMER Userguide.
_COL = {
    "target_name": 0,   # HMM model NAME  (e.g. Lectin_legB)
    "target_acc": 1,    # HMM ACC         (e.g. PF00139.27)
    "query_name": 3,    # protein id from proteome.fasta
    "i_evalue": 12,     # domain independent E-value
    "dom_score": 13,    # domain bit score
    "ali_from": 17,     # alignment start on the query (protein)
    "ali_to": 18,       # alignment end on the query (protein)
}
_MIN_COLS = 22


def run(cfg: Config) -> int:
    """Run hmmscan on the proteome and parse the hits."""
    cfg.ensure_dirs()
    if not cfg.proteome.exists():
        raise FileNotFoundError(f"proteome not found: {cfg.proteome}")
    if not cfg.hmm_db.exists():
        raise FileNotFoundError(
            f"pressed HMM database not found: {cfg.hmm_db} (run `gwiscan setup-db` first)"
        )
    external.require("hmmscan")

    domtbl = cfg.result("hmm_domtbl.txt")
    external.log(f"[{datetime.now()}] Running hmmscan with --cut_ga...")
    external.run([
        "hmmscan",
        "--cut_ga",
        "--cpu", cfg.THREADS,
        "--domtblout", domtbl,
        "--tblout", cfg.result("hmm_tbl.txt"),
        "-o", cfg.result("hmm_full.out"),
        cfg.hmm_db,
        cfg.proteome,
    ])
    external.log(f"[{datetime.now()}] hmmscan complete.")

    n = parse_domtbl(cfg, domtbl)
    external.log(f"[{datetime.now()}] hmmscan step done.")
    return n


def parse_domtbl(cfg: Config, domtbl=None) -> int:
    """Parse a hmmscan --domtblout file into intermediate/hmm_hits.tsv."""
    domtbl = domtbl or cfg.result("hmm_domtbl.txt")
    if not domtbl.exists():
        raise FileNotFoundError(f"domtblout not found: {domtbl} (run hmmscan first)")

    # Family/superfamily modes map each Pfam to its curated family. Architecture
    # mode has no per-domain family — the label IS the HMM's own domain name (e.g.
    # B_lectin, Pkinase), so the map is empty and every hit uses that name (the
    # fallback below). architecture.run then combines the domains per protein.
    pfam_to_family = {} if cfg.is_architecture else io.pfam_to_family(cfg.family_map)
    rows, unmapped = [], set()

    with open(domtbl) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split()
            if len(cols) < _MIN_COLS:
                continue

            hmm_acc = cols[_COL["target_acc"]]           # e.g. PF00139.27
            pfam_base = hmm_acc.split(".")[0]            # e.g. PF00139
            family = pfam_to_family.get(pfam_base)
            if family is None:
                family = cols[_COL["target_name"]]       # fall back to raw HMM name
                unmapped.add(pfam_base)

            rows.append([
                cols[_COL["query_name"]],
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
    if unmapped and not cfg.is_architecture:
        external.log(
            f"[WARN] {len(unmapped)} Pfam accession(s) not in the family table, "
            f"used raw HMM name instead: {', '.join(sorted(unmapped))}"
        )
    return len(rows)
