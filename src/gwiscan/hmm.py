#!/usr/bin/env python3
"""
#####################################################################################################
#                                                                                                   #
# hmm.py - hmmsearch against the HMM database and domtbl parsing (the `search-hmm` stage).          #
#                                                                                                   #
# Uses hmmsearch (the profiles are the query, the proteome the target DB): far faster than hmmscan  #
# for a genome-wide search and it saturates the cores. In hmmsearch --domtblout the target is the   #
# protein and the query is the HMM, so column 1 is the protein id and columns 4-5 the HMM name/acc. #
# The column indices below are named accordingly and covered by tests/test_hmm_parse.py.            #
#                                                                                                   #
#####################################################################################################
"""

from __future__ import annotations

import re
from datetime import datetime

from . import external, io
from .config import Config
from .schema import HIT_HEADER

# Pfam's three curated per-model bit-score cutoffs, any of which hmmsearch can use
# with --cut_<kind>: GA (gathering, the balanced default), TC (trusted, strictest:
# the weakest included member), NC (noise, loosest: the best excluded non-member).
_CUTOFF_NAMES = {"ga": "gathering", "tc": "trusted", "nc": "noise"}


def cutoff_kind(cfg) -> str:
    """The model cutoff the run uses (HMM_CUTOFF): 'ga' (default), 'tc' or 'nc'.
    Ignored when HMM_EVALUE is set (that switches to an E-value cutoff instead)."""
    kind = str(getattr(cfg, "HMM_CUTOFF", "ga") or "ga").strip().lower()
    if kind not in _CUTOFF_NAMES:
        raise RuntimeError(
            f"invalid HMM_CUTOFF {getattr(cfg, 'HMM_CUTOFF', '')!r}; use ga, tc or nc "
            f"(gathering / trusted / noise)."
        )
    return kind


def cutoff_help(kind: str) -> str:
    """Guidance shown when a custom identifying HMM lacks the chosen cutoff line.
    Kept in one place so preflight and setup-db give identical instructions."""
    tag, name = kind.upper(), _CUTOFF_NAMES[kind]
    return (
        f"The HMM search runs `hmmsearch --cut_{kind}`, which applies each model's own "
        f"{tag} ({name}) cutoff, so every identifying model must declare one. A profile "
        f"from hmmbuild has no {tag} line unless its source alignment carried one. Add a "
        f"{tag} cutoff to the profile (build it from a Stockholm alignment with a "
        f"`#=GF {tag} <seq> <dom>` line, or insert a `{tag}  <seq> <dom>;` line into the "
        f".hmm header), or identify this family by a Pfam accession instead."
    )


# Back-compat: the GA-specific help constant other code/tests reference.
GA_THRESHOLD_HELP = cutoff_help("ga")


def has_cutoff_thresholds(hmm_path, kind: str = "ga") -> bool:
    """True when every model in an HMMER profile file declares the given cutoff line
    (GA/TC/NC).

    An .hmm file may hold several concatenated models; ``hmmsearch --cut_<kind>``
    aborts the whole search if any one of them lacks that cutoff, so this requires
    one such line per model (and at least one model)."""
    line_re = re.compile(rf"^{kind.upper()}\s")
    n_models = n_cut = 0
    with open(hmm_path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("HMMER3/"):     # each model begins with the format line
                n_models += 1
            elif line_re.match(line):
                n_cut += 1
    return n_models >= 1 and n_cut == n_models


def has_ga_thresholds(hmm_path) -> bool:
    """True when every model declares a GA threshold (has_cutoff_thresholds for GA)."""
    return has_cutoff_thresholds(hmm_path, "ga")


def threshold_args(cfg) -> list:
    """The hmmsearch reporting/inclusion cutoff for this run: a model cutoff line by
    default (``--cut_ga`` / ``--cut_tc`` / ``--cut_nc``, per HMM_CUTOFF), or a
    sequence E-value (``-E``) when HMM_EVALUE is set. A model cutoff requires every
    model to declare that line; an E-value works for custom HMMs that carry none and
    gives a more sensitive search."""
    evalue = str(getattr(cfg, "HMM_EVALUE", "") or "").strip()
    return ["-E", evalue] if evalue else [f"--cut_{cutoff_kind(cfg)}"]


def uses_ga(cfg) -> bool:
    """True when the run scores HMMs by a model cutoff line (no HMM_EVALUE set), i.e.
    when every identifying model must declare that GA/TC/NC cutoff. (Named for the GA
    default; use cutoff_kind(cfg) for which one.)"""
    return not str(getattr(cfg, "HMM_EVALUE", "") or "").strip()


def custom_hmm_cutoff_error(hmm_path, family: str, cfg) -> str | None:
    """An explanatory error string if a custom identifying HMM lacks the cutoff line
    the run needs (the HMM_CUTOFF kind), else None. None also when HMM_EVALUE is set
    (no cutoff line required). Shared by preflight (logs it) and setup-db (raises)."""
    if not uses_ga(cfg):
        return None
    kind = cutoff_kind(cfg)
    if has_cutoff_thresholds(hmm_path, kind):
        return None
    return (f"custom HMM for '{family}' has no {kind.upper()} ({_CUTOFF_NAMES[kind]}) "
            f"thresholds: {hmm_path}. {cutoff_help(kind)}")


def custom_hmm_ga_error(hmm_path, family: str) -> str | None:
    """GA-specific form of custom_hmm_cutoff_error (kept for back-compat)."""
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
    cutoff = threshold_args(cfg)
    external.log(f"[{datetime.now()}] Running hmmsearch with {' '.join(map(str, cutoff))}...")
    external.run([
        "hmmsearch",
        *cutoff,                                # --cut_ga (default) or -E <HMM_EVALUE>
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
