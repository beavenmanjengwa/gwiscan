# GWIscan (Genome-Wide Identification Scan)

GWIscan is a pipeline for the genome-wide identification and annotation of gene
families and superfamilies which combines BLAST, profile hidden Markov models, and
InterProScan. It runs in family, superfamily, or architecture modes, on a single
proteome or across multiple species. It provides detailed annotation features
including gene coordinates, domain architecture, physicochemical properties,
signal-peptide and transmembrane topology, subcellular localization, and GO terms.
It can also perform multiple sequence alignment, MEME motif discovery, and
phylogenetic tree generation.

## Installation

Install from Bioconda (pulls in the search and alignment tools automatically):

```bash
conda create -n gwiscan -c bioconda -c conda-forge gwiscan
conda activate gwiscan
```

Or build the environment from source:

```bash
mamba env create -f environment.yml      # or 'conda env create' (mamba resolves dependencies faster)
conda activate gwiscan
pip install -e .                         # the `gwiscan` command
```

Either way you get `HMMER`, `diamond`, `seqkit`, `mafft`, `clipkit`, `iqtree`,
`meme`, `weblogo`, and R. DeepTMHMM runs through
`pybiolib` (`pip install pybiolib`; already included in the source environment). InterProScan uses the EBI web service
by default (set `EBI_EMAIL`); to run it offline, install InterProScan and set
`INTERPRO_MODE: local` and `INTERPROSCAN_BIN` if it's not executable at system level.

[TargetP 2.0](https://services.healthtech.dtu.dk/services/TargetP-2.0/) and
[DeepLoc 2.1](https://services.healthtech.dtu.dk/services/DeepLoc-2.1/) should be
installed manually, as they need an academic license from DTU. Install each one
(its own conda environment is easiest) and add its path to `TARGETP_BIN` and
`DEEPLOC_BIN` (in the config file, an env var, or `--targetp-bin` /
`--deeploc-bin`).

Before a long run, check the tools, packages, and inputs with:

```bash
gwiscan preflight -C project/
```

## Project layout

Run `gwiscan` from a project folder, or point at one with `-C /path`:

```
project/
├── config.yaml                 # settings (THREADS, MODE, EBI_EMAIL, ...)
├── config/
│   ├── family.tsv              # family table (or superfamily.tsv for MODE: superfamily)
│   └── species.tsv             # species list (multi-species only)
├── input/proteome.fasta        # the proteome to scan (single species)
├── db/
│   ├── hmm/                    # custom HMMs; Pfam HMMs download here
│   └── blast/                  # BlastModel FASTAs (one per family)
├── final_results/              # the main results: final TSV/XLSX, GFF3, provenance
├── intermediate/               # working files from each stage
└── logs/
```

More detail is in `docs/project-layout.txt`.

## The family table

Copy `config/family.tsv.example` to `config/family.tsv` and edit it (or
`superfamily.tsv` for `MODE: superfamily`, which adds a `Superfamily` column).
The columns are:

| Column | Meaning |
|--------|---------|
| `Family` | Family name. Used in the outputs and domain ids (for example `GNA` gives `AthGNA001.1`). |
| `PfamModel` | A Pfam accession (`PF01453`), a custom HMM file in `db/hmm/` (`name.hmm`), or `-` for none. A custom identifying HMM must declare GA (gathering) thresholds, because the HMM search runs `hmmscan --cut_ga`; preflight checks this and tells you how to fix it. |
| `BlastModel` | The DIAMOND query FASTA in `db/blast/`. Every family needs one. |

```
Family    PfamModel   BlastModel
GNA       PF01453     AAA33346.1.fasta
CRA       CRA.hmm     ABL98074.1.fasta
EUL       -           ABW73993.1.fasta
```

## Architecture mode (domain combinations)

This mode identifies genes that encode proteins with a given domain architecture: a
**primary** domain that defines and
seeds the family, plus one or more **required** domains that must also be present.
Both are Pfam HMMs.

The search runs in two hmmscan passes:

1. hmmscan the whole proteome against the **primary** HMM(s) to get candidates.
2. hmmscan only those candidates against the primary+required HMMs, keeping the
   ones that also carry every required domain. Those are the final candidates.

The final candidates then go through InterProScan, ProtParam, TargetP,
DeepTMHMM, DeepLoc, coordinate mapping, and result compilation. `family` holds
the architecture name in every output.

Set `MODE: architecture` in `config.yaml`, then copy
`config/architecture.tsv.example` to `config/architecture.tsv`. One row per
architecture:

| Column | Meaning |
|--------|---------|
| `Architecture` | Family name, used in the outputs. |
| `Primary` | The defining Pfam domain that seeds the genome-wide search. One slot; alternatives allowed with `\|`. Pick the domain that best defines the family. |
| `Required` | Pfam domain(s) that must also be present, `+`-separated (AND). Any slot may list alternatives with `\|` (OR), so `PF00069\|PF07714` accepts either Pkinase or Pkinase_Tyr. |
| `Class` | Optional rollup label grouping architectures. |

```
Architecture   Primary   Required          Class
G-LecRLK       PF01453   PF00069|PF07714   LecRLK
L-LecRLK       PF00139   PF00069|PF07714   LecRLK
C-LecRLK       PF00059   PF00069|PF07714   LecRLK
```

Run it with `gwiscan run -C project/ --mode architecture` (set `EBI_EMAIL`, or
`INTERPRO_MODE: local`, since InterProScan annotates the final candidates). The
results are the usual `final_results/gwiscan_results.tsv`/`.xlsx`, with `family`
holding the architecture and the domain architecture filled in from InterProScan.

## Usage

```bash
gwiscan run -C project/ --threads 8 --ebi-email you@example.com
```

The main results go to `project/final_results/` and the working files to
`project/intermediate/`. Send both (and `logs/`) somewhere else with `-o/--output`.
A command cheat-sheet is in [`docs/COMMANDS.md`](docs/COMMANDS.md).

You can resume a run, stop after a stage, skip a stage, or add an off-by-default one:

```bash
gwiscan run --list-stages                      # show the stage names
gwiscan run -C project/ --from-stage merge     # start from a stage
gwiscan run -C project/ --until compile        # stop after a stage
gwiscan run -C project/ --skip weblogo         # skip a stage
gwiscan run -C project/ --add iqtree           # add the off-by-default tree
```

A default run stops at the annotation table. Phylogenetic tree building is off by
default; enable it with `--add iqtree` (or `ADD_STAGES: [iqtree]` in `config.yaml`).

`weblogo` and `meme` are optional and skipped if their tool is not installed. Any
flag above also works as a `config.yaml` setting or env var.

To scan several species, add `config/species.tsv` (`Prefix`, `Proteome`):

```
Prefix   Proteome
Ath      genomes/Athaliana.fasta
Gma      genomes/Gmax.fasta
```

Each species has its own `final_results/<Prefix>/` and `intermediate/<Prefix>/`
folders. Species are scheduled stage by stage: every species finishes a stage before
any starts the next, and within a stage the species run in parallel.

Some stages run one species at a time instead: `deeptmhmm`, `targetp`, and `deeploc`
(listed in `SERIAL_STAGES`), because each needs a resource that several species
cannot share at once, such as a single GPU, a large in-memory model, or a
license-limited tool. The per-species thread count is rebalanced for each stage, so a stage running
one species at a time still uses every core.

Three settings control how many species run at once, from a global default to a
per-stage override:

| Setting | Type | What it controls |
|---------|------|------------------|
| `SPECIES_PARALLEL` | one number | the default number of species run at once, for every stage (0 = auto from cores/threads) |
| `SERIAL_STAGES` | list of stage names | stages forced to one species at a time |
| `STAGE_PARALLEL` | map `{stage: N}` | a per-stage number that overrides the default and `SERIAL_STAGES` |

Per stage the count is `STAGE_PARALLEL[stage]` if set, else `1` if the stage is in
`SERIAL_STAGES`, else `SPECIES_PARALLEL`. For example:

```yaml
SPECIES_PARALLEL: 4              # default: every other stage runs 4 species at once
SERIAL_STAGES: [deeptmhmm]       # deeptmhmm runs 1 species at a time
STAGE_PARALLEL: {interpro: 2}    # interpro runs 2 species at a time
```

One species failing does not stop the others; use `--only-species Ath,Gma` or
`--retry-failed` to rerun a subset.

After the species finish, a cross-species summary is written to the top-level
`final_results/`:

- `all_species_summary.tsv` / `.xlsx` — a families (rows) by species (columns)
  matrix of the member protein count, with per-family and per-species totals. Every
  configured family has a row, so families absent in a species show as 0.
- `all_species_members.tsv` — every species' members stacked into one table with a
  `species` column.

This works in every mode (in architecture mode the rows are the architectures).

Each stage can also run on its own (`gwiscan <stage>`); `gwiscan <stage> --help`
lists the options. Settings are read in this order: built-in defaults, then
`config.yaml`, then env vars, then command-line flags. Any setting can be given
as an environment variable named `GWISCAN_<KEY>` (e.g. `GWISCAN_THREADS`); the
bare `<KEY>` name still works, but the `GWISCAN_` prefix is preferred for the few
generic names (`THREADS`, `OUTPUT`, `MODE`, `ANNOTATION`) so they cannot pick up an
unrelated value already in your environment.

## Outputs

Deliverables go to `final_results/`; working files to `intermediate/`. In
multi-species runs each species has its own `final_results/<Prefix>/` and
`intermediate/<Prefix>/`.

| Path | What it is |
|------|------------|
| `final_results/gwiscan_results.tsv` / `.xlsx` | The annotation table, one row per member (the `.xlsx` groups columns by source). Genomic coordinates are added when a GFF or GTF annotation is provided. |
| `final_results/gwiscan_members.gff3` | The members as genome features. |
| `final_results/provenance.txt` | Tool versions, settings, and input checksums. |
| `intermediate/candidates/candidates_merged.tsv` / `.fasta` | Merged candidates from the HMM and BLAST searches. |
| `intermediate/interproscan/interproscan.tsv` | Domain and GO annotations. |
| `intermediate/protparam/protparam.tsv` | Physicochemical properties. |
| `intermediate/prediction/{targetp,deeptmhmm,deeploc}.tsv` | Signal peptide, transmembrane topology, and subcellular localization. |
| `intermediate/msa/{Family}_aligned.fasta` | Per-family alignment. |
| `intermediate/weblogo/`, `meme/`, `trees/` | Per-family logos, motifs, and trees (`trees/` only with `--add iqtree`). |

Each stage also writes a log under `logs/` (one folder per species in multi-species runs).

## Development

```bash
pip install -e ".[dev]"
pytest
```
