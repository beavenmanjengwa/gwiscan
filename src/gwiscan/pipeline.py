#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# pipeline.py - Full-pipeline orchestration (the linear `run` for one proteome).                   #
#                                                                                                  #
# Runs the stages in order, teeing each stage's stdout to logs/<stage>.log. `run` supports         #
# resuming: --from-stage skips everything before a given stage (its outputs must already be on     #
# disk from an earlier run), --until stops after one, --skip drops named stages entirely (e.g. an  #
# optional tool that isn't installed) without aborting the rest. `gwiscan run --list-stages`        #
# prints the valid keys. For a DAG with automatic resume/parallelism instead, use the Snakemake     #
# workflow, which calls the same per-stage subcommands.                                            #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

import contextlib
import sys

from . import (
    candidates,
    confirm,
    coords,
    diamond,
    domain_bed,
    domains,
    external,
    hmm,
    iqtree,
    logos,
    mature,
    meme,
    msa,
    preflight,
    figures,
    provenance,
    score,
    setupdb,
    trimal,
)
from . import compile as compile_stage
from .config import Config
from .features import deeploc, deeptmhmm, interpro, protparam, targetp

# The full stage order as (key, label, func, logname). key is what --from-stage,
# --until and --skip name a stage by; it doubles as the `gwiscan <key>` subcommand
# name wherever one exists. "setup-shared" is dropped from the list run() actually
# uses when include_shared_setup=False (multi-species: done once, not per species).
#
# func is a lambda that looks up e.g. preflight.run at CALL time rather than the
# module-level function object directly -- STAGES is built once at import, and
# tests (and anyone scripting against the module) monkeypatch stage functions by
# setting them on the stage's module (monkeypatch.setattr(pipeline.preflight,
# "run", ...)); a bound reference captured here at import time would never see
# that patch.
STAGES = [
    ("preflight", "00 Pre-flight checks", lambda cfg: preflight.run(cfg), "00_preflight"),
    ("setup-shared", "00 Setup shared databases", lambda cfg: setupdb.setup_shared(cfg), "00_setup_shared"),
    ("setup-db", "00 Setup proteome db", lambda cfg: setupdb.setup_proteome(cfg), "00_setup_proteome"),
    ("search-hmm", "01 HMMscan", lambda cfg: hmm.run(cfg), "01_hmmscan"),
    ("search-diamond", "02 DIAMOND BLASTp", lambda cfg: diamond.run(cfg), "02_diamond"),
    ("merge", "03 Merge candidates", lambda cfg: candidates.run(cfg), "03_merge_candidates"),
    ("score", "04 Family detectability", lambda cfg: score.run(cfg), "04_detectability"),
    ("interpro", "05 InterProScan", lambda cfg: interpro.run(cfg), "05_interproscan"),
    ("confirm", "06 Confirm final candidates", lambda cfg: confirm.run(cfg), "06_confirm"),
    ("protparam", "07 ProtParam", lambda cfg: protparam.run(cfg), "07_protparam"),
    ("targetp", "08 TargetP 2.0", lambda cfg: targetp.run(cfg), "08_targetp"),
    ("deeptmhmm", "09 DeepTMHMM", lambda cfg: deeptmhmm.run(cfg), "09_deeptmhmm"),
    ("deeploc", "10 DeepLoc 2.1", lambda cfg: deeploc.run(cfg), "10_deeploc"),
    ("coords", "11 Genomic coordinates", lambda cfg: coords.run(cfg), "11_coords"),
    ("compile", "12 Compile results", lambda cfg: compile_stage.run(cfg), "12_compile_results"),
    ("domain-bed", "13 Domain BED", lambda cfg: domain_bed.run(cfg), "13_domain_bed"),
    ("extract-domains", "14 Extract domains", lambda cfg: domains.run(cfg), "14_extract_domains"),
    ("extract-mature", "15 Extract mature sequences", lambda cfg: mature.run(cfg), "15_extract_mature"),
    ("msa", "16 MAFFT MSA", lambda cfg: msa.run(cfg), "16_mafft"),
    ("trim", "17 Trim alignments (trimAl)", lambda cfg: trimal.run(cfg), "17_trimal"),
    ("weblogo", "18 WebLogo", lambda cfg: logos.run(cfg), "18_weblogo"),
    ("meme", "19 MEME motifs", lambda cfg: meme.run(cfg), "19_meme"),
    ("iqtree", "20 IQ-TREE", lambda cfg: iqtree.run(cfg), "20_iqtree"),
    ("figures", "21 ProtParam figures", lambda cfg: figures.run(cfg), "21_figures"),
    ("provenance", "22 Provenance", lambda cfg: provenance.run(cfg), "22_provenance"),
]
STAGE_KEYS = [key for key, *_ in STAGES]

# Stages backed by an OPTIONAL external tool, mapped to the Config attribute
# naming that tool's binary. If the tool isn't installed `run` auto-skips the
# stage with a notice instead of aborting the whole pipeline -- preflight already
# reports these as optional [WARN]s, so a hard abort here contradicted that.
# weblogo/meme/iqtree are leaf stages (nothing consumes their output). `trim` is
# safe to include even though iqtree reads its output, because iqtree explicitly
# falls back to the untrimmed alignment when the trimmed file is absent. Feature
# stages whose output IS consumed with no fallback (targetp/deeploc/deeptmhmm ->
# compile) are deliberately NOT here: skipping them would break the final table.
OPTIONAL_TOOL_STAGES = {
    "trim": "TRIMAL_BIN",     # iqtree falls back to the untrimmed alignment if skipped
    "weblogo": "WEBLOGO_BIN",
    "meme": "MEME_BIN",
    "iqtree": "IQTREE_BIN",
    "figures": "RSCRIPT_BIN",  # ProtParam boxplots via R; skipped if Rscript absent
}


class _Tee:
    """Write to several streams at once (console + per-stage log file)."""

    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for s in self.streams:
            s.write(data)

    def flush(self):
        for s in self.streams:
            s.flush()


def _stage(cfg, label, func, logname):
    log_path = cfg.logs / f"{logname}.log"
    external.log("")
    external.log(f"[STAGE] {label}")
    external.log("-" * 55)
    with open(log_path, "w") as fh, contextlib.redirect_stdout(_Tee(sys.stdout, fh)):
        func(cfg)
    external.log(f"[DONE] {label} (log: {log_path})")


def plan_stages(stages, from_stage="", until_stage="", skip=()) -> list:
    """The ordered slice of ``stages`` (key, label, func, logname) tuples to run,
    given --from-stage / --until / --skip. Pure and side-effect-free, so the
    resume/stop/skip logic is testable without invoking any external tool.

    --from-stage / --until name a stage to start at / stop after (inclusive); both
    must be one of ``stages``' keys. --skip names stages the caller should omit
    from execution (validated here against ``stages`` so a typo fails fast, but
    still returned in the slice -- the caller decides what "skip" means, e.g. still
    logging a [SKIP] line). Raises ValueError, naming the bad key and the valid
    ones, on any unrecognised key, or if --until falls at or before --from-stage
    (nothing would run)."""
    keys = [key for key, *_ in stages]

    def _index(value, flag):
        if not value:
            return None
        if value not in keys:
            raise ValueError(
                f"{flag} {value!r} is not a stage in this run. Valid stages: {', '.join(keys)}"
            )
        return keys.index(value)

    start_idx = _index(from_stage, "--from-stage")
    until_idx = _index(until_stage, "--until")
    start = start_idx if start_idx is not None else 0
    end = until_idx + 1 if until_idx is not None else len(keys)

    if start >= end:
        raise ValueError(
            f"--until {until_stage!r} occurs at or before --from-stage {from_stage!r}; nothing to run"
        )

    unknown_skips = sorted(set(skip) - set(keys))
    if unknown_skips:
        raise ValueError(
            f"--skip has unknown stage(s) {', '.join(unknown_skips)}. Valid stages: {', '.join(keys)}"
        )

    return stages[start:end]


def run(cfg: Config, include_shared_setup: bool = True) -> None:
    """Run the full pipeline for one proteome.

    include_shared_setup presses the shared HMM database and validates the BLAST
    model FASTAs. The multi-species driver does that once up front, then calls this
    per species with include_shared_setup=False -- each species still builds its own
    proteome db and writes into its own intermediate/<SPECIES>/ tree, staying independent.

    cfg.FROM_STAGE / cfg.UNTIL_STAGE / cfg.SKIP_STAGES (--from-stage / --until /
    --skip) resume, stop early, or drop named stages -- see plan_stages().
    """
    cfg.ensure_dirs()
    tag = f" [{cfg.SPECIES}]" if cfg.SPECIES else ""
    external.log("=" * 56)
    external.log(f" GWIscan - genome-wide identification and annotation pipeline for gene families or superfamilies{tag}")
    external.log("=" * 56)
    external.log(f" Threads         : {cfg.THREADS}")
    external.log(f" DIAMOND e-value : {cfg.DIAMOND_EVALUE}")
    external.log(f" Project dir     : {cfg.root}")
    if cfg.SPECIES:
        external.log(f" Species         : {cfg.SPECIES}  ({cfg.proteome})")
    external.log("=" * 56)

    stages = [s for s in STAGES if s[0] != "setup-shared" or include_shared_setup]
    skip = set(cfg.SKIP_STAGES or ())
    try:
        selected = plan_stages(stages, cfg.FROM_STAGE, cfg.UNTIL_STAGE, skip)
    except ValueError as e:
        raise RuntimeError(str(e))

    # Auto-skip optional-tool leaf stages whose tool isn't installed (explicit
    # --skip already covers the rest). Explicit --skip wins the message if both.
    auto_skip = {
        key: getattr(cfg, attr)
        for key, attr in OPTIONAL_TOOL_STAGES.items()
        if key not in skip and not external.available(getattr(cfg, attr))
    }

    if cfg.FROM_STAGE:
        external.log(f"[resume] Starting from stage '{cfg.FROM_STAGE}' "
                     f"-- earlier stages' outputs must already be on disk.")
    if cfg.UNTIL_STAGE:
        external.log(f"[resume] Stopping after stage '{cfg.UNTIL_STAGE}'.")

    for key, label, func, logname in selected:
        if key in skip:
            external.log(f"\n[SKIP] {label} (--skip {key})")
            continue
        if key in auto_skip:
            external.log(
                f"\n[AUTO-SKIP] {label}: optional tool '{auto_skip[key]}' not found. "
                f"Install it (or set its *_BIN config/flag) to enable this stage; "
                f"pass --skip {key} to silence this notice."
            )
            continue
        _stage(cfg, label, func, logname)

    external.log("\n" + "=" * 56)
    external.log(" GWIscan complete.")
    external.log(f" Results:    {cfg.final_dir}")
    external.log("=" * 56)
