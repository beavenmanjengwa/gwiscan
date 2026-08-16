#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# cli.py - GWIscan command-line interface.                                                         #
#                                                                                                  #
# One entry point (`gwiscan`) with one subcommand per pipeline stage, plus `run` for the whole     #
# pipeline. Every subcommand takes the same global options (project dir, config file, parameter    #
# overrides) and dispatches to a stage's run(cfg) function.                                        #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

import argparse

from . import (
    __version__,
    architecture,
    candidates,
    confirm,
    coords,
    diamond,
    domains,
    external,
    figures,
    hmm,
    iqtree,
    clipkit,
    domain_bed,
    logos,
    mature,
    meme,
    msa,
    pipeline,
    preflight,
    provenance,
    score,
    setupdb,
    species,
)
from . import compile as compile_stage
from .config import Config
from .features import deeploc, deeptmhmm, interpro, protparam, targetp

# command name -> (help text, stage callable taking a Config)
COMMANDS = {
    "preflight": ("Check tools, packages and inputs", preflight.run),
    "setup-db": ("Build the HMM + DIAMOND databases", setupdb.run),
    "search-hmm": ("Run hmmscan and parse hits", hmm.run),
    "search-diamond": ("Run two-round DIAMOND BLASTp", diamond.run),
    "merge": ("Merge HMM + DIAMOND candidates", candidates.run),
    "architecture": ("Identify proteins with a given domain architecture", architecture.run),
    "score": ("Per-family detectability profile", score.run),
    "protparam": ("Physicochemical properties", protparam.run),
    "targetp": ("Signal/transit peptides", targetp.run),
    "deeptmhmm": ("Transmembrane topology", deeptmhmm.run),
    "deeploc": ("Subcellular localization", deeploc.run),
    "interpro": ("Domain + GO annotation", interpro.run),
    "confirm": ("Confirm final candidates via InterProScan", confirm.run),
    "domain-bed": ("BED of family domain coordinates", domain_bed.run),
    "coords": ("Genomic coordinates of members", coords.run),
    "compile": ("Join tables into final TSV + XLSX", compile_stage.run),
    "extract-domains": ("Extract per-family domain sequences", domains.run),
    "extract-mature": ("Per-family TargetP mature sequences", mature.run),
    "msa": ("Per-family MAFFT alignment", msa.run),
    "trim": ("Trim alignments before trees", clipkit.run),
    "weblogo": ("Per-family conservation logos", logos.run),
    "meme": ("Per-family MEME motif discovery", meme.run),
    "iqtree": ("Per-family ML trees", iqtree.run),
    "figures": ("ProtParam distribution figures + stats", figures.run),
    "provenance": ("Write run provenance", provenance.run),
    "run": ("Run the full pipeline", None),
}


def _run(cfg: Config) -> None:
    """`run` dispatch: multi-species if a species manifest is present and no single
    species was pinned with --species; otherwise a single-proteome run."""
    if cfg.species_manifest.exists() and not cfg.SPECIES:
        species.run(cfg)
    else:
        pipeline.run(cfg)

_OVERRIDE_KEYS = (
    "OUTPUT",
    "THREADS",
    "VERBOSE",
    "MODE",
    "DIAMOND_EVALUE",
    "DIAMOND_IDENTITY",
    "DIAMOND_COVERAGE_R2",
    "DIAMOND_SENSITIVE_R2",
    "DIAMOND_BSR",
    "PRIMARY_TRANSCRIPT",
    "ANNOTATION",
    "EBI_EMAIL",
    "INTERPRO_APPL",
    "INTERPRO_MODE",
    "INTERPROSCAN_BIN",
    "INTERPRO_LOOKUP",
    "TARGETP_BIN",
    "TARGETP_ORGANISM",
    "DEEPLOC_BIN",
    "DEEPLOC_MODEL",
    "DEEPTMHMM_VERSION",
    "PROTPARAM_FORMATS",
    "SPECIES_PREFIX",
    "MEME_NMOTIFS",
    "WEBLOGO_BIN",
    "MEME_BIN",
    "CLIPKIT_BIN",
    "CLIPKIT_MODE",
    "RSCRIPT_BIN",
    "IQTREE_BIN",
    "IQTREE_MODEL",
    "IQTREE_BOOTSTRAP",
    "IQTREE_SEED",
    "SPECIES_PARALLEL",
    "SPECIES",
    "PROTEOME",
    "FROM_STAGE",
    "UNTIL_STAGE",
    "SKIP_STAGES",
    "ADD_STAGES",
    "SPECIES_ONLY",
    "RETRY_FAILED",
)


def _add_global_opts(p: argparse.ArgumentParser) -> None:
    p.add_argument("-C", "--workdir", default=None,
                   help="Project directory: inputs/config/db (default: current directory)")
    p.add_argument("-o", "--output", dest="OUTPUT", default=None,
                   help="Where final_results/, intermediate/ and logs/ go (default: the project directory)")
    p.add_argument("--config", default=None, help="Path to config.yaml")
    p.add_argument("-p", "--threads", dest="THREADS", type=int, default=None,
                   help="Worker threads / CPUs (default 4)")
    p.add_argument("--verbose", dest="VERBOSE", action=argparse.BooleanOptionalAction, default=None,
                   help="Stream each stage's full output to the terminal. Default: off — "
                        "only stage start/done is shown; full detail always goes to logs/.")
    p.add_argument("--mode", dest="MODE", default=None,
                   choices=["family", "superfamily", "architecture"],
                   help="Mode: family (flat), superfamily (grouped rollup), or "
                        "architecture (domain-combination, hmmscan-only)")
    p.add_argument("--evalue", dest="DIAMOND_EVALUE", default=None,
                   help="DIAMOND E-value cutoff")
    p.add_argument("--identity", dest="DIAMOND_IDENTITY", type=int, default=None,
                   help="DIAMOND round-2 min %% identity (round 1 is E-value only)")
    p.add_argument("--coverage-r2", dest="DIAMOND_COVERAGE_R2", type=int, default=None,
                   help="DIAMOND round-2 min %% query coverage")
    p.add_argument("--diamond-sensitive-r2", dest="DIAMOND_SENSITIVE_R2",
                   action=argparse.BooleanOptionalAction, default=None,
                   help="Round-2 --sensitive (default off)")
    p.add_argument("--diamond-bsr", dest="DIAMOND_BSR", type=float, default=None,
                   help="Blast Score Ratio seed cutoff for no-HMM families (default 0.4)")
    p.add_argument("--primary-transcript", dest="PRIMARY_TRANSCRIPT",
                   action=argparse.BooleanOptionalAction, default=None,
                   help="Input is one protein per gene (default true; --no-primary-transcript if it has isoforms)")
    p.add_argument("--annotation", dest="ANNOTATION", default=None,
                   help="Genome annotation (GTF/GFF3) for the coords stage; overrides input/annotation.*")
    p.add_argument("--ebi-email", dest="EBI_EMAIL", default=None,
                   help="Email for the EBI InterProScan API")
    p.add_argument("--interpro-appl", dest="INTERPRO_APPL", default=None,
                   help="InterProScan 6 member db(s), comma-separated (default Pfam)")
    p.add_argument("--interpro-mode", dest="INTERPRO_MODE", default=None,
                   choices=["api", "local"],
                   help="InterProScan: api (EBI REST) or local (installed interproscan.sh)")
    p.add_argument("--interproscan-bin", dest="INTERPROSCAN_BIN", default=None,
                   help="Path/name of the local InterProScan executable (local mode)")
    p.add_argument("--interpro-lookup", dest="INTERPRO_LOOKUP",
                   action=argparse.BooleanOptionalAction, default=None,
                   help="Local mode: use the online precalc lookup (default off = offline)")
    p.add_argument("--targetp-bin", dest="TARGETP_BIN", default=None,
                   help="Path or PATH name of the TargetP 2.0 executable")
    p.add_argument("--targetp-organism", dest="TARGETP_ORGANISM", default=None,
                   help="TargetP organism group (pl = plant)")
    p.add_argument("--deeploc-bin", dest="DEEPLOC_BIN", default=None,
                   help="Path or PATH name of the DeepLoc 2.1 executable")
    p.add_argument("--deeploc-model", dest="DEEPLOC_MODEL", default=None,
                   choices=["Accurate", "Fast"],
                   help="DeepLoc model: Accurate (ProtT5, ~32GB RAM) or Fast (ESM1b)")
    p.add_argument("--deeptmhmm-version", dest="DEEPTMHMM_VERSION", default=None,
                   help="Pinned DeepTMHMM biolib version (default 1.0.24, runs locally)")
    p.add_argument("--protparam-formats", dest="PROTPARAM_FORMATS", default=None,
                   type=lambda s: [x.strip() for x in s.split(",") if x.strip()],
                   help="ProtParam output formats, comma-separated: tsv,xlsx,csv (default tsv,xlsx)")
    p.add_argument("--species-prefix", dest="SPECIES_PREFIX", default=None,
                   help="Prefix for systematic domain ids, e.g. Ath -> AthGNA001.1")
    p.add_argument("--meme-nmotifs", dest="MEME_NMOTIFS", type=int, default=None,
                   help="Number of MEME motifs per family (default 15)")
    p.add_argument("--weblogo-bin", dest="WEBLOGO_BIN", default=None,
                   help="WebLogo executable: PATH name or absolute path (default weblogo)")
    p.add_argument("--meme-bin", dest="MEME_BIN", default=None,
                   help="MEME executable: PATH name or absolute path (default meme)")
    p.add_argument("--clipkit-bin", dest="CLIPKIT_BIN", default=None,
                   help="ClipKIT executable: PATH name or absolute path (default clipkit)")
    p.add_argument("--clipkit-mode", dest="CLIPKIT_MODE", default=None,
                   help="ClipKIT trimming mode: smart-gap (default), gappy, kpic, kpic-smart-gap, ...")
    p.add_argument("--rscript-bin", dest="RSCRIPT_BIN", default=None,
                   help="Rscript executable for the figures stage (default Rscript)")
    p.add_argument("--iqtree-bin", dest="IQTREE_BIN", default=None,
                   help="IQ-TREE executable (default iqtree; install: conda install bioconda::iqtree)")
    p.add_argument("--iqtree-model", dest="IQTREE_MODEL", default=None,
                   help="IQ-TREE substitution model (default MFP = ModelFinder Plus)")
    p.add_argument("--iqtree-bootstrap", dest="IQTREE_BOOTSTRAP", type=int, default=None,
                   help="IQ-TREE ultrafast bootstrap replicates (default 1000; 0 = off)")
    p.add_argument("--iqtree-seed", dest="IQTREE_SEED", type=int, default=None,
                   help="IQ-TREE RNG seed for reproducible trees (default 12345)")
    p.add_argument("--species-parallel", dest="SPECIES_PARALLEL", type=int, default=None,
                   help="Multi-species: concurrent species (0 = auto = cores // threads)")
    p.add_argument("--species", dest="SPECIES", default=None,
                   help="Run a single species by prefix (namespaces outputs to <prefix>/)")
    p.add_argument("--proteome", dest="PROTEOME", default=None,
                   help="Proteome path override (used with --species)")


def _add_run_opts(p: argparse.ArgumentParser) -> None:
    """Resume/stop/skip flags for `run` only -- meaningless for a single-stage
    subcommand, so kept off the shared global options."""
    p.add_argument("--from-stage", dest="FROM_STAGE", default=None,
                   help="Resume from this stage onward, skipping everything before it "
                        "(that stage's inputs -- the earlier stages' outputs -- must "
                        "already be on disk). See --list-stages.")
    p.add_argument("--until", dest="UNTIL_STAGE", default=None,
                   help="Stop after this stage (inclusive). See --list-stages.")
    p.add_argument("--skip", dest="SKIP_STAGES", default=None,
                   type=lambda s: [x.strip() for x in s.split(",") if x.strip()],
                   help="Comma-separated stage(s) to skip entirely (e.g. an optional "
                        "tool you don't have installed), continuing past them instead "
                        "of aborting the run. See --list-stages.")
    p.add_argument("--add", dest="ADD_STAGES", default=None,
                   type=lambda s: [x.strip() for x in s.split(",") if x.strip()],
                   help="Comma-separated stage(s) that are off by default to turn on "
                        "for this run: the phylogeny workflow trim + iqtree. "
                        "'--add iqtree' also runs its trim prerequisite. See --list-stages.")
    p.add_argument("--list-stages", action="store_true",
                   help="Print the ordered stage keys (for --from-stage/--until/--skip) and exit.")
    p.add_argument("--only-species", dest="SPECIES_ONLY", default=None,
                   type=lambda s: [x.strip() for x in s.split(",") if x.strip()],
                   help="Multi-species: run only these manifest prefixes (comma-separated).")
    p.add_argument("--retry-failed", dest="RETRY_FAILED", action="store_true", default=None,
                   help="Multi-species: re-run only the species marked FAILED in the previous "
                        "run's logs/species_status.tsv.")


def _config_from_args(args) -> Config:
    overrides = {k: getattr(args, k, None) for k in _OVERRIDE_KEYS}
    return Config.load(root=args.workdir, configfile=args.config, overrides=overrides)


# Shown by `gwiscan --help`. Kept separate from the module docstring above: that
# docstring is the source-file header comment (its ASCII box mangles once argparse
# rewraps it into one paragraph), not user-facing help text.
_CLI_DESCRIPTION = (
    "Genome-wide identification and annotation pipeline for gene families or superfamilies. "
    "One subcommand per pipeline stage, plus `run` for the whole pipeline. "
    "Run `gwiscan <stage> --help` for a stage's options."
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="gwiscan", description=_CLI_DESCRIPTION)
    parser.add_argument("--version", action="version",
                        version=f"gwiscan {__version__}")
    sub = parser.add_subparsers(dest="command", metavar="<stage>")
    for name, (help_text, _) in COMMANDS.items():
        sp = sub.add_parser(name, help=help_text)
        _add_global_opts(sp)
        if name == "run":
            _add_run_opts(sp)
    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not args.command:
        parser.print_help()
        return 2

    if args.command == "run" and getattr(args, "list_stages", False):
        for key in pipeline.STAGE_KEYS:
            print(key)
        return 0

    cfg = _config_from_args(args)
    func = COMMANDS[args.command][1] or _run   # "run" dispatches via _run
    try:
        func(cfg)
    except (FileNotFoundError, RuntimeError) as e:
        external.log(f"[ERROR] {e}")
        return 1
    return 0
