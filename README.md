# GWIscan (Genome-Wide Identification Scan)

Genome-wide identification and annotation pipeline for gene families or superfamilies. You
list the families in a config table, so it works for any families. The bundled
example is a set of plant lectins.

For each family it runs DIAMOND BLASTp, and for families that have a profile HMM
it also runs HMMER `hmmscan`. It then checks the candidates with InterProScan and
annotates them with physicochemical properties, signal and transit peptides,
transmembrane topology, subcellular localization, and domain and GO terms.
Finally it aligns each family and builds conservation logos, MEME motifs, and
IQ-TREE trees. It is one Python package with a single command, `gwiscan`, whose
subcommands are the pipeline stages. You can run one proteome, or many species at
once.

## Installation

```bash
conda env create -f environment.yml      # conda tools + Python deps
conda activate gwiscan
pip install -e .                         # the `gwiscan` command
```

This puts all the conda tools in the `gwiscan` environment, ready to use: HMMER,
`diamond`, `seqkit`, `mafft`, `trimal`, `iqtree`, `meme`, `weblogo`, and
`pybiolib` (which runs DeepTMHMM). InterProScan uses the EBI web service by
default (nothing to install; you just set `EBI_EMAIL`). To run it offline, install
InterProScan and set `INTERPRO_MODE: local` and `INTERPROSCAN_BIN`.

[TargetP 2.0](https://services.healthtech.dtu.dk/services/TargetP-2.0/) and
[DeepLoc 2.1](https://services.healthtech.dtu.dk/services/DeepLoc-2.1/) should be
installed manually, as they need an academic license from DTU. Install each one
(its own conda environment is easiest) and add its path to `TARGETP_BIN` and
`DEEPLOC_BIN` (in the config file, an env var, or `--targetp-bin` /
`--deeploc-bin`). Any path works.

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
| `PfamModel` | A Pfam accession (`PF01453`), a custom HMM file in `db/hmm/` (`name.hmm`), or `-` for none. |
| `BlastModel` | The DIAMOND query FASTA in `db/blast/`. Every family needs one. |

```
Family    PfamModel   BlastModel
GNA       PF01453     AAA33346.1.fasta
CRA       CRA.hmm     ABL98074.1.fasta
EUL       -           ABW73993.1.fasta
```

## Usage

```bash
gwiscan run -C project/ --threads 8 --ebi-email you@example.com
```

The main results go to `project/final_results/` and the working files to
`project/intermediate/`. Send both (and `logs/`) somewhere else with `-o/--output`.

You can resume, stop, or skip stages, so a long run that fails part way does not
have to start over:

```bash
gwiscan run --list-stages                      # show the stage names
gwiscan run -C project/ --from-stage merge     # start from a stage
gwiscan run -C project/ --until compile        # stop after a stage
gwiscan run -C project/ --skip meme,weblogo    # skip stages and keep going
```

If an optional tool (`trim`, `weblogo`, `meme`, `iqtree`) is not installed, that
stage is skipped instead of stopping the run. You can point to any of these tools
with its `*_BIN` setting. All the flags above also work as `config.yaml` settings
or env vars.

To scan several species at once, add `config/species.tsv` (`Prefix`, `Proteome`).
Each species runs on its own in parallel, into `final_results/<Prefix>/` and
`intermediate/<Prefix>/`:

```
Prefix   Proteome
Ath      genomes/Athaliana.fasta
Gma      genomes/Gmax.fasta
```

`SPECIES_PARALLEL` sets how many run at once (0 means auto). One species failing
does not stop the others; re-run a subset with `--only-species Ath,Gma` or just
the failed ones with `--retry-failed`.

There is also a Snakemake workflow that runs the same stages as a resumable,
parallel pipeline:

```bash
snakemake -s workflow/Snakefile --cores 8
```

Each stage can also run on its own (`gwiscan <stage>`); `gwiscan <stage> --help`
lists the options. Settings are read in this order: built-in defaults, then
`config.yaml`, then env vars, then command-line flags.

## Outputs

| File | What it is |
|------|-------------|
| `final_results/gwiscan_results.tsv` / `.xlsx` | The final results table (the `.xlsx` groups the columns by source tool) |
| `final_results/gwiscan_members.gff3` | The members as genome features |
| `final_results/provenance.txt` | Tool versions, settings, input checksums |
| `intermediate/candidates_merged.tsv` / `.fasta` | The merged HMM and DIAMOND candidates |
| `intermediate/interproscan.tsv` | Domain and GO annotations |
| `intermediate/{protparam,targetp,deeptmhmm,deeploc}.tsv` | The per-tool annotation tables |
| `intermediate/msa/{Family}_aligned.fasta` | Per-family MAFFT alignment |
| `intermediate/weblogo/`, `meme/`, `trees/` | Per-family logos, motifs, and trees |

Each stage writes a log under `logs/` (one folder per species in multi-species runs).

## Development

```bash
pip install -e ".[dev]"
pytest
```
