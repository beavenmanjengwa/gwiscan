# GWIscan (Genome-Wide Identification scan)

GWIscan is a pipeline for the genome-wide identification and in silico characterization of
gene families. It combines homology-based search, profile hidden Markov model search, and 
InterPro family signatures for identification. It runs in family, multi-family, or architecture
modes, on a single proteome or across multiple species. It provides annotation including gene  
coordinates, domain architecture, physicochemical properties, signal peptides and transmembrane 
topology, subcellular localization, and GO terms. It can also perform multiple sequence alignment,
MEME motif discovery, and phylogenetic tree construction.

## Installation

Install from Bioconda (pulls in the search and alignment tools automatically):

```bash
mamba create -n gwiscan -c bioconda -c conda-forge gwiscan
# Or use conda
conda create -n gwiscan -c bioconda -c conda-forge gwiscan
conda activate gwiscan
```

Or build the environment from source:

```bash
mamba env create -f environment.yml      
# or 'conda env create', mamba can resolves dependencies faster
conda activate gwiscan
pip install -e .                         # the `gwiscan` command
```

Either way you get `HMMER`, `diamond`, `seqkit`, `mafft`, `clipkit`, `iqtree`,
`meme`, `weblogo`, and R. DeepTMHMM runs through
`pybiolib` (`pip install pybiolib`; already included in the source environment). InterProScan uses the EBI web service
by default (set `EBI_EMAIL`); to run it offline, install InterProScan and set
`INTERPRO_MODE: local` and `INTERPROSCAN_BIN` if it's not executable at system level.

The EBI web service has two API versions, chosen with `IPRSCAN_VERSION` (config /
env) or `--iprscan-version {5,6}` (CLI). **Default is `5`**, the stable
InterProScan 5 service. `6` targets InterProScan 6, whose match-lookup step has a
known server-side fault that fails every job; if you hit it, GWIscan stops and
tells you to rerun with `--iprscan-version 5`. The two versions name the same
databases differently (v5 `PfamA`/`Panther`, v6 `Pfam`/`PANTHER`); GWIscan derives
the applications from the family table and maps them to the chosen version
automatically, so you never set application names by hand. The InterProScan and
InterPro release that produced a run are recorded in
`intermediate/<species>/interproscan/interproscan.manifest.txt`, in the
`Tool_versions` sheet of `gwiscan_results.xlsx`, and in `provenance.txt`.

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
│   ├── family.tsv              # family table (or multi-family.tsv for MODE: multi-family)
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
`multi-family.tsv` for `MODE: multi-family`, which adds a `Multifamily` column).
The columns are:

| Column | Meaning |
|--------|---------|
| `Family` | Family name. Used in the outputs and domain ids (for example `GNA` gives `AthGNA001.1`). |
| `PfamModel` | How the family is IDENTIFIED: a Pfam accession (`PF01453`), a custom HMM file in `db/hmm/` (`name.hmm`), or `-` for none (DIAMOND only). A custom identifying HMM must declare the model cutoff the run uses, because the HMM search runs `hmmsearch --cut_<HMM_CUTOFF>` by default (`ga` gathering, `tc` trusted, or `nc` noise); preflight checks this (or set `HMM_EVALUE` to use an E-value cutoff instead). |
| `BlastModel` | The DIAMOND query FASTA in `db/blast/`. Every family needs one. |
| `InterProModel` | Optional. Accession(s) that CONFIRM a candidate via InterProScan: a Pfam (`PF...`), CDD (`cd...`), or PANTHER (`PTHR...`) id, pipe-separated for alternatives (`PF01476\|PTHR27007`). GWIscan enables the matching applications automatically. Blank falls back to the `PfamModel` accession. A family with no usable Pfam model (`PfamModel` = `-`) can still be confirmed by giving a CDD or PANTHER accession here. |

```
Family    PfamModel   BlastModel          InterProModel
GNA       PF01453     AAA33346.1.fasta
EUL       -           ABW73993.1.fasta    PTHR31257
```

## Architecture mode

This mode identifies genes that encode proteins with a given domain architecture: a
**primary** domain that defines and
seeds the family, plus one or more **required** domains that must also be present.
Both are Pfam HMMs.

The search runs in two HMMER hmmsearch passes:

1. hmmsearch the whole proteome against the **primary** HMM(s) to get candidates.
2. hmmsearch only those candidates against the primary+required HMMs.

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

To scan several species, add `config/species.tsv` (`Prefix`, `Proteome`,
`Annotation`):

```
Prefix   Proteome                  Annotation
Ath      genomes/Athaliana.fasta   genomes/Athaliana.gtf
Gma      genomes/Gmax.fasta        genomes/Gmax.gtf
```

`Annotation` is optional. Provide the GTF/GFF3 that matches the proteome, and each member's 
chromosome, start, end, and strand will be included. Leave the cell empty and the species is
still analysed, just without chromosomal coordinates.

Each species has its own `final_results/<Prefix>/` and `intermediate/<Prefix>/`
folders. Species are scheduled stage by stage: every species finishes a stage before
any starts the next, and within a stage the species run in parallel.

Some stages run one species at a time instead: `deeptmhmm`, `targetp`, and `deeploc`
(listed in `SERIAL_STAGES`), because each needs a resource that several species
cannot share at once, such as a single GPU, a large in-memory model, or a
license-limited tool. The per-species thread count is rebalanced for each stage, so a stage running
one species at a time still uses every core.

The annotation stage is the most demanding stage, since `DeepTMHMM`, `TargetP`, and `DeepLoc` run neural 
network models whose runtime scales with the number of confirmed members and the available hardware.
In multi-species runs these stages process one species at a time by default, a behaviour controlled 
by SERIAL_STAGES, while the per-species thread count is rebalanced within each stage to use the 
available cores, and the number of species running in parallel is limited by SPECIES_PARALLEL. 

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

- `all_species_summary.tsv` / `.xlsx`: a families (rows) by species (columns)
  matrix of the member protein count, with per-family and per-species totals. Every
  configured family has a row, so families absent in a species show as 0.
- `all_species_members.tsv`: every species' members stacked into one table with a
  `species` column.

This works in every mode (in architecture mode the rows are the architectures).

Each stage can also run on its own (`gwiscan <stage>`); `gwiscan <stage> --help`
lists the options. Settings are read in this order: built-in defaults, then
`config.yaml`, then environment variables, then command-line flags. Any setting can be
set as an environment variable, either as `<KEY>` or as `GWISCAN_<KEY>` (for example
`GWISCAN_THREADS`). For most settings either form works. For a few generic names
(`THREADS`, `OUTPUT`, `MODE`, `ANNOTATION`) use the `GWISCAN_` prefix, since your
shell or another program may already use the bare name for something else, and GWIscan 
would otherwise pick up that unrelated value.

## Outputs

Compiled results go to `final_results/`; working files to `intermediate/`. In
multi-species runs each species has its own `final_results/<Prefix>/` and
`intermediate/<Prefix>/`.

| Path | What it is |
|------|------------|
| `final_results/gwiscan_results.tsv` / `.xlsx` | The annotation table, one row per member (the `.xlsx` groups columns by source). Genomic coordinates are added when a GFF or GTF annotation is provided. |
| `final_results/gwiscan_members.gff3` | The members as genome features. |
| `final_results/provenance.txt` | Tool versions, settings, and input checksums. |
| `intermediate/candidates/candidates_merged.tsv` / `.fasta` | Merged candidates from the HMM and DIAMOND searches. |
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