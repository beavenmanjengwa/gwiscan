#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# config.py - Runtime configuration and project directory layout.                                  #
#                                                                                                  #
# A single Config object carries the layout and every tunable parameter, resolved once in order:   #
# built-in defaults < config.yaml < environment variables < CLI flags. Paths resolve relative to   #
# a project directory (--workdir, default: current dir), so the tool is installed once and run     #
# against any dataset.                                                                             #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

import yaml


def _csv_list(value):
    """Parse a comma-separated string into a clean list (for env/CLI overrides)."""
    if isinstance(value, list):
        return value
    return [x.strip() for x in str(value).split(",") if x.strip()]


def _to_bool(value):
    """Parse a truthy env string ('true', '1', 'yes', 'on') into a bool."""
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in ("1", "true", "yes", "on")


# Intermediate working files are grouped into per-tool subfolders so the working
# directory is navigable instead of one flat dump. Config.result() routes a
# filename to its subfolder via _result_subdir(); everything read/written through
# result() moves together, and single-file stage outputs (targetp.tsv, deeploc.tsv,
# domains.bed, chromosome_map.tsv, ...) stay at the intermediate root.
_CANDIDATE_FILES = {
    "candidates.fasta", "candidates_merged.tsv", "candidates_primary.fasta",
    "final_candidates.tsv", "final_candidates.fasta",
}


def _result_subdir(name: str) -> str:
    """The intermediate subfolder for a result filename ('' = the root)."""
    if name.startswith("interproscan") or name == "interpro_input.fasta":
        return "interproscan"
    if name in _CANDIDATE_FILES:
        return "candidates"
    if name.startswith("hmm_") or name.startswith("arch_"):
        return "hmm"
    if name.startswith("diamond_") or name in ("blast_hits.tsv", "family_detectability.tsv"):
        return "diamond"
    # Per-protein prediction annotators (TargetP signal/transit peptide, DeepLoc
    # subcellular localization, DeepTMHMM membrane topology).
    if name in ("targetp.tsv", "deeploc.tsv", "deeptmhmm.tsv"):
        return "prediction"
    return ""

# Parameter defaults. Keys are UPPER_CASE — one name per setting across the
# config.yaml key, the environment variable, and the Config field.
DEFAULTS = {
    "THREADS": 4,
    # Where final_results/, intermediate/ and logs/ are written. Empty = the
    # project directory (--workdir), so outputs sit next to the inputs by default.
    # Set it (CLI -o/--output or OUTPUT here) to send outputs to a separate
    # directory; inputs, config, and db are still read from the project directory.
    "OUTPUT": "",
    # Reporting mode:
    #   "family"       — flat families (one identifying model per family).
    #   "superfamily"  — families grouped under a Superfamily column (superfamily.tsv);
    #                    adds a superfamily rollup. e.g. Lectin over GNA, Legume, CRA...
    #   "architecture" — domain-COMBINATION mode. A protein is classified by the SET
    #                    of Pfam domains it carries (e.g. LecRLK = a lectin domain AND
    #                    a kinase domain on the same chain). hmmscan only — the
    #                    co-occurrence of the component domains is the criterion, so no
    #                    DIAMOND and no InterProScan. Rules live in architecture.tsv.
    "MODE": "family",
    "DIAMOND_EVALUE": "1e-5",
    "DIAMOND_IDENTITY": 30,        # round 2 min % identity (native seeds; round 1 is E-value only)
    "DIAMOND_COVERAGE_R2": 70,     # round 2 query coverage % (native seeds)
    "DIAMOND_SENSITIVE_R2": False, # round 2 --sensitive; default off (native-to-native)
    # Round-2 seeds: HMM-validated (round-1 hits that also pass hmmscan) when a
    # family has an HMM; otherwise Blast Score Ratio -- keep round-1 subjects whose
    # bitscore / model self-bitscore >= this (length-normalized, Rasko et al. 2005).
    "DIAMOND_BSR": 0.4,
    # Per-family detectability (score stage): a family counts as "concordant" when
    # hmmscan and the independent DIAMOND round 1 agree at least this much
    # (Jaccard). Reporting cutoff only -- the raw metrics are always written.
    "CONCORDANCE_MIN": 0.7,
    # The input proteome is one protein per gene (primary transcript). Default
    # true — that is the expected input. Set false only if your proteome still
    # contains splice isoforms, to get the inflation advisory (see preflight).
    # The pipeline never collapses isoforms itself.
    # Genome annotation (GTF or GFF3) for chromosomal coordinates. Supply the
    # REPRESENTATIVE annotation matching the input proteome (one transcript per
    # gene), so every protein resolves to exactly one locus. Empty = look for
    # input/annotation.gtf, then .gff3/.gff; absent entirely = the coords stage
    # is skipped and the run continues without chromosomal columns.
    "ANNOTATION": "",
    "PRIMARY_TRANSCRIPT": True,
    "EBI_EMAIL": "",
    # InterProScan 6 member database(s), comma-separated. Default Pfam + CDD;
    # InterProScan runs all requested DBs on the small candidate set for the full
    # domain architecture + GO. Add more (SMART, PROSITEProfiles, SUPERFAMILY, ...).
    "INTERPRO_APPL": "Pfam,CDD,PROSITE-profiles",
    # How InterProScan runs:
    #   "api"   — EBI InterProScan 6 REST API (needs EBI_EMAIL, internet).
    #   "local" — a local install (INTERPROSCAN_BIN), offline-capable, no email.
    "INTERPRO_MODE": "api",
    # Local mode: install InterProScan system-wide (on PATH) so "interproscan.sh"
    # resolves by name and no directory has to be given. Set a full path here only
    # for a non-PATH install.
    "INTERPROSCAN_BIN": "interproscan.sh",
    # Local mode only: use InterProScan's precalculated-match lookup (an online
    # service). Default off so local runs work fully offline (adds -dp).
    "INTERPRO_LOOKUP": False,
    # Gene Ontology: InterProScan reports GO ids only, so go-basic.obo is fetched
    # once into db/go/ to add readable GO names (a copy already there is reused).
    "GO_OBO_URL": "https://purl.obolibrary.org/obo/go/go-basic.obo",
    "TARGETP_BIN": "targetp",
    "TARGETP_ORGANISM": "pl",
    "DEEPLOC_BIN": "deeploc2",
    # DeepLoc model: 'Accurate' (ProtT5, high quality, ~32GB RAM, downloaded on
    # first use) or 'Fast' (ESM1b, lighter/faster). See features/deeploc.py.
    "DEEPLOC_MODEL": "Accurate",
    # DeepTMHMM biolib version. 1.0.24 runs LOCALLY (seconds); newer tags (e.g.
    # 1.0.57) route to biolib's cloud and queue for minutes. See deeptmhmm.py.
    "DEEPTMHMM_VERSION": "1.0.24",
    # ProtParam output formats — any of tsv, xlsx, csv. tsv feeds the pipeline join.
    "PROTPARAM_FORMATS": ["tsv", "xlsx"],
    # Prefix for systematic domain ids (e.g. "Ath" -> AthGNA001.1). One species
    # per run, so this is a single constant that keeps ids unique across runs.
    "SPECIES_PREFIX": "",
    # MEME motif discovery: number of motifs per family. Everything else uses
    # MEME's own defaults.
    "MEME_NMOTIFS": 15,
    # Optional-tool executables, like TARGETP_BIN/DEEPLOC_BIN/IQTREE_BIN: a PATH
    # name (resolved by name) or an absolute path to the binary. So WebLogo/MEME can
    # live in a separate conda env without being on the PATH `gwiscan` runs from.
    "WEBLOGO_BIN": "weblogo",
    "MEME_BIN": "meme",
    # trimAl: trims each MAFFT alignment before tree-building. TRIMAL_METHOD is a
    # trimAl heuristic (automated1 [default], gappyout, strict, strictplus);
    # written as -automated1 on the command line. Only IQ-TREE uses the trimmed
    # alignment (WebLogo/MEME keep the full one). Optional -- if trimAl isn't
    # installed the `trim` stage is skipped and iqtree uses the untrimmed alignment.
    "TRIMAL_BIN": "trimal",
    "TRIMAL_METHOD": "automated1",
    # IQ-TREE per-family phylogenetic trees. Install via conda so iqtree is on
    # PATH: conda install bioconda::iqtree
    # ProtParam distribution figures + stats (figures stage). Runs the bundled
    # R/ggplot2 script; optional, auto-skipped if Rscript isn't installed.
    "RSCRIPT_BIN": "Rscript",
    "IQTREE_BIN": "iqtree",
    "IQTREE_MODEL": "MFP",      # ModelFinder Plus (auto model selection)
    "IQTREE_BOOTSTRAP": 1000,   # ultrafast bootstrap replicates (-B); 0 = no bootstrap
    "IQTREE_SEED": 12345,       # RNG seed (-seed) so trees/support values are reproducible
    # Multi-species mode: how many species pipelines run concurrently.
    #   0 (default) = auto = cores // THREADS, so species × per-species threads
    #                 totals ~one thread per core (full CPU use, no oversubscription,
    #                 bounded RAM / InterProScan-API load).
    #   N > 0       = run exactly N species at once (raise it if the machine has
    #                 the RAM/cores; lower it to be gentler on the EBI API).
    # Each species is still independent — the smallest proteomes finish first.
    "SPECIES_PARALLEL": 0,
    # `run` resume/stop/skip (see pipeline.plan_stages; `gwiscan run --list-stages`
    # prints the valid keys). FROM_STAGE/UNTIL_STAGE resume from or stop after a
    # named stage (inclusive); earlier/later stages' outputs must already be on
    # disk. SKIP_STAGES (comma-separated) drops named stages entirely and keeps
    # going -- e.g. an optional tool (meme, weblogo, iqtree) you don't have
    # installed, so it no longer aborts the whole run.
    "FROM_STAGE": "",
    "UNTIL_STAGE": "",
    "SKIP_STAGES": [],
    # Multi-species subset selection (see species.py). SPECIES_ONLY
    # (--only-species) runs only the named manifest prefixes; RETRY_FAILED
    # (--retry-failed) re-runs just the species marked FAILED in the previous
    # run's persisted logs/species_status.tsv. Both are ignored in single-species
    # mode. RETRY_FAILED wins if both are given.
    "SPECIES_ONLY": [],
    "RETRY_FAILED": False,
}

# Type casters for settings coming from strings (env vars, CLI). The env var
# name IS the config key (both UPPER_CASE); keys without an entry stay as-is.
_CASTERS = {
    "THREADS": int,
    "DIAMOND_IDENTITY": int,
    "DIAMOND_COVERAGE_R2": int,
    "DIAMOND_SENSITIVE_R2": _to_bool,
    "DIAMOND_BSR": float,
    "CONCORDANCE_MIN": float,
    "PRIMARY_TRANSCRIPT": _to_bool,
    "PROTPARAM_FORMATS": _csv_list,
    "MEME_NMOTIFS": int,
    "SPECIES_PARALLEL": int,
    "INTERPRO_LOOKUP": _to_bool,
    "IQTREE_BOOTSTRAP": int,
    "IQTREE_SEED": int,
    "SKIP_STAGES": _csv_list,
    "SPECIES_ONLY": _csv_list,
    "RETRY_FAILED": _to_bool,
}


@dataclass
class Config:
    root: Path
    THREADS: int = 4
    OUTPUT: str = ""
    MODE: str = "family"
    DIAMOND_EVALUE: str = "1e-5"
    DIAMOND_IDENTITY: int = 30
    DIAMOND_COVERAGE_R2: int = 70
    DIAMOND_SENSITIVE_R2: bool = False
    DIAMOND_BSR: float = 0.4
    CONCORDANCE_MIN: float = 0.7
    ANNOTATION: str = ""
    PRIMARY_TRANSCRIPT: bool = True
    EBI_EMAIL: str = ""
    INTERPRO_APPL: str = "Pfam,CDD,PROSITE-profiles"
    INTERPRO_MODE: str = "api"
    INTERPROSCAN_BIN: str = "interproscan.sh"
    INTERPRO_LOOKUP: bool = False
    GO_OBO_URL: str = "https://purl.obolibrary.org/obo/go/go-basic.obo"
    TARGETP_BIN: str = "targetp"
    TARGETP_ORGANISM: str = "pl"
    DEEPLOC_BIN: str = "deeploc2"
    DEEPLOC_MODEL: str = "Accurate"
    DEEPTMHMM_VERSION: str = "1.0.24"
    PROTPARAM_FORMATS: list = field(default_factory=lambda: ["tsv", "xlsx"])
    SPECIES_PREFIX: str = ""
    MEME_NMOTIFS: int = 15
    WEBLOGO_BIN: str = "weblogo"
    MEME_BIN: str = "meme"
    TRIMAL_BIN: str = "trimal"
    TRIMAL_METHOD: str = "automated1"
    RSCRIPT_BIN: str = "Rscript"
    IQTREE_BIN: str = "iqtree"
    IQTREE_MODEL: str = "MFP"
    IQTREE_BOOTSTRAP: int = 1000
    IQTREE_SEED: int = 12345
    SPECIES_PARALLEL: int = 0
    # `run` resume/stop/skip -- see DEFAULTS above and pipeline.plan_stages.
    FROM_STAGE: str = ""
    UNTIL_STAGE: str = ""
    SKIP_STAGES: list = field(default_factory=list)
    # Multi-species subset selection -- see DEFAULTS above and species.py.
    SPECIES_ONLY: list = field(default_factory=list)
    RETRY_FAILED: bool = False
    # Multi-species runtime (set per species by the multi-species driver; empty in
    # single-species mode). SPECIES namespaces the outputs into
    # final_results/<SPECIES>/, intermediate/<SPECIES>/ and logs/<SPECIES>/;
    # PROTEOME overrides the input proteome path for that species. Both empty =>
    # ordinary single-species run.
    SPECIES: str = ""
    PROTEOME: str = ""

    # --- Directory / file layout (single source of truth) ---
    @property
    def input_dir(self) -> Path:
        return self.root / "input"

    @property
    def proteome(self) -> Path:
        """The input proteome. In multi-species mode PROTEOME points at that
        species' proteome (resolved relative to the project root if not absolute);
        otherwise the default input/proteome.fasta."""
        if self.PROTEOME:
            p = Path(self.PROTEOME)
            return p if p.is_absolute() else (self.root / p)
        return self.input_dir / "proteome.fasta"

    @property
    def annotation(self):
        """The genome annotation accompanying the proteome, or None if absent.
        ANNOTATION names it explicitly (per species in multi-species mode);
        otherwise input/annotation.{gtf,gff3,gff} is used, whichever is present."""
        if self.ANNOTATION:
            p = Path(self.ANNOTATION)
            return p if p.is_absolute() else (self.root / p)
        for suffix in ("gtf", "gff3", "gff"):
            candidate = self.input_dir / f"annotation.{suffix}"
            if candidate.exists():
                return candidate
        return None

    @property
    def config_dir(self) -> Path:
        return self.root / "config"

    @property
    def is_architecture(self) -> bool:
        """True when running the domain-combination mode (MODE: architecture)."""
        return str(self.MODE).lower() == "architecture"

    @property
    def family_map(self) -> Path:
        """The family table: config/superfamily.tsv in superfamily mode, else
        config/family.tsv. The filename matches MODE."""
        name = "superfamily.tsv" if str(self.MODE).lower() == "superfamily" else "family.tsv"
        return self.config_dir / name

    @property
    def architecture_map(self) -> Path:
        """The architecture rules table (config/architecture.tsv), used only in
        MODE: architecture. One row per domain-combination class."""
        return self.config_dir / "architecture.tsv"

    @property
    def species_manifest(self) -> Path:
        """config/species.tsv (Prefix, Proteome). Its presence switches `run`
        into multi-species mode."""
        return self.config_dir / "species.tsv"

    @property
    def db_dir(self) -> Path:
        return self.root / "db"

    @property
    def hmm_dir(self) -> Path:
        return self.db_dir / "hmm"

    @property
    def hmm_db(self) -> Path:
        return self.hmm_dir / "all_models.hmm"

    @property
    def primary_hmm_db(self) -> Path:
        """Architecture mode: the pressed db of just the PRIMARY domain HMM(s), used
        for the genome-wide seed pass. hmm_db holds primary+required for pass 2."""
        return self.hmm_dir / "primary_models.hmm"

    @property
    def blast_dir(self) -> Path:
        return self.db_dir / "blast"

    @property
    def go_obo(self) -> Path:
        """Cached Gene Ontology file (go-basic.obo) for GO id -> name mapping."""
        return self.db_dir / "go" / "go-basic.obo"

    @property
    def proteome_db(self) -> Path:
        """DIAMOND db of the proteome (both rounds search it). Per-species in
        multi-species mode so runs never share/overwrite each other's db."""
        name = f"{self.SPECIES}_proteome_db" if self.SPECIES else "proteome_db"
        return self.blast_dir / name

    @property
    def output_root(self) -> Path:
        """Where final_results/, intermediate/ and logs/ are written. OUTPUT if set (absolute, or
        relative to the current directory), else the project directory. Inputs,
        config, and db always stay under the project directory (root)."""
        if self.OUTPUT:
            out = Path(self.OUTPUT)
            return out if out.is_absolute() else (Path.cwd() / out)
        return self.root

    @property
    def results(self) -> Path:
        """Intermediate / working files: intermediate/ (namespaced to
        intermediate/<SPECIES>/ in multi-species mode). Every stage's working
        outputs land here; the headline deliverables go to final_dir instead.
        Method name kept `results` (and result()) so call sites are unchanged."""
        base = self.output_root / "intermediate"
        return base / self.SPECIES if self.SPECIES else base

    @property
    def final_dir(self) -> Path:
        """The headline deliverables — final TSV/XLSX, members GFF3, provenance —
        in a TOP-LEVEL final_results/ (a sibling of intermediate/, not nested
        inside it) so the main output isn't buried. Namespaced to
        final_results/<SPECIES>/ in multi-species mode."""
        base = self.output_root / "final_results"
        return base / self.SPECIES if self.SPECIES else base

    @property
    def logs(self) -> Path:
        base = self.output_root / "logs"
        return base / self.SPECIES if self.SPECIES else base

    def result(self, name: str) -> Path:
        """Path to a file inside the intermediate/ working directory.

        Grouped into a per-tool subfolder (interproscan/, candidates/, hmm/,
        diamond/) when the filename maps to one, else the intermediate root. The
        chosen directory is created so writers (including external tools that will
        not make it themselves) can write straight to the returned path."""
        sub = _result_subdir(name)
        directory = self.results / sub if sub else self.results
        directory.mkdir(parents=True, exist_ok=True)
        return directory / name

    @property
    def protparam_dir(self) -> Path:
        """intermediate/<species>/protparam/ — every ProtParam working file (the
        table, and the figures stage's stats + outlier CSVs) grouped in one folder."""
        return self.results / "protparam"

    def ensure_dirs(self) -> None:
        for d in (self.results, self.logs):
            d.mkdir(parents=True, exist_ok=True)

    @classmethod
    def load(cls, root=None, configfile=None, overrides=None) -> "Config":
        """Build a Config from defaults < config.yaml < env < CLI overrides."""
        root = Path(root or Path.cwd()).resolve()
        values = dict(DEFAULTS)

        cfg_path = Path(configfile) if configfile else (root / "config.yaml")
        if cfg_path.exists():
            with open(cfg_path) as fh:
                loaded = yaml.safe_load(fh) or {}
            values.update({k: v for k, v in loaded.items() if k in DEFAULTS})

        # Env var name == config key (both UPPER_CASE); cast string values.
        for key in DEFAULTS:
            raw = os.environ.get(key)
            if raw is not None:
                values[key] = _CASTERS.get(key, str)(raw)

        if overrides:
            values.update({k: v for k, v in overrides.items() if v is not None})

        return cls(root=root, **values)
