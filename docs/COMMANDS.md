# gwiscan commands

## Install

```bash
# from bioconda (pulls in the tools automatically)
conda create -n gwiscan -c bioconda -c conda-forge gwiscan
conda activate gwiscan
```

Or from source:

```bash
# go to the gwiscan repo (the folder with environment.yml and pyproject.toml)
cd /path/to/gwiscan

# create the env with all tools + python deps (mamba preferred: much faster solve)
mamba env create -f environment.yml

# activate it
conda activate gwiscan

# install gwiscan into the env; source edits take effect immediately
pip install -e .
```

```bash
# if the env already exists, update it from environment.yml instead of creating
mamba env update -n gwiscan -f environment.yml

# remove the env entirely
conda env remove -n gwiscan
```

## Activate (every session)

```bash
conda activate gwiscan
```

## Help

```bash
# list all stages/subcommands
gwiscan --help

# all options for the full pipeline
gwiscan run --help

# options for one stage, e.g. gwiscan interpro --help
gwiscan <stage> --help

# ordered stage names (for --from-stage/--until/--skip)
gwiscan run --list-stages

# version
gwiscan --version
```

## Full runs

```bash
# multi-family level
gwiscan run -C /path/to/ProjectDir --mode multi-family -p 8

# architecture mode
gwiscan run -C /path/to/ProjectDir --mode architecture -p 8

# family level
gwiscan run -C /path/to/ProjectDir --mode family -p 8
```

## Part of the pipeline

```bash
# stop after the merge stage
gwiscan run -C ProjectDir --until merge

# resume from interpro
gwiscan run -C ProjectDir --from-stage interpro

# skip a reporting stage
gwiscan run -C ProjectDir --skip meme,weblogo

# add the off-by-default phylogeny workflow (trim + iqtree)
gwiscan run -C ProjectDir --add iqtree
```

`trim` and `iqtree` are off by default (slow, and absent from the final table).
`--add iqtree` runs both: adding the tree pulls in its `trim` prerequisite, so you
never trim without building the tree. `gwiscan trim` / `gwiscan iqtree` still run
on their own.

## Multi-species

```bash
# rerun only FAILED species
gwiscan run -C ProjectDir --retry-failed

# run only these species
gwiscan run -C ProjectDir --only-species Ath,Osa

# one species, resumed from a stage
gwiscan run -C /path/to/ProjectDir --only-species Stu --from-stage deeptmhmm

# default width: 2 species at a time; each stage receives a matching thread budget
gwiscan run -C ProjectDir -p 8 --species-parallel 2
```

## Single stage

```bash
# check tools/inputs
gwiscan preflight -C ProjectDir

# build HMM + DIAMOND db
gwiscan setup-db -C ProjectDir

# join tables -> final TSV/XLSX
gwiscan compile -C ProjectDir
```

## Stage order

| Phase | Stages (in order) |
|-------|-------------------|
| Setup | `preflight`, `setup-shared`, `setup-db` |
| Identify | `search-hmm` + `search-diamond`, `merge`, `score` |
| Annotate | `interpro`, `confirm`, `protparam`, `targetp`, `deeptmhmm`, `deeploc`, `coords` |
| Compile | `compile` |
| Domains | `domain-bed`, `extract-domains`, `extract-mature` |
| Align | `msa` |
| Logos/motifs | `weblogo`, `meme` |
| Tree | `trim`, `iqtree` |
| Finish | `figures`, `provenance` |

`trim` and `iqtree` are off by default; add them with `--add iqtree`. Print the
live list with `gwiscan run --list-stages`.

## Common run flags

```
-C, --workdir DIR       project dir (inputs/config/db)
-o, --output DIR        where results/logs go
--config FILE           path to config.yaml
-p, --threads N         threads (default 4)
--mode MODE             family | multi-family | architecture
--from-stage NAME       resume from this stage
--until NAME            stop after this stage
--skip NAMES            skip stage(s), comma-separated
--add NAMES             turn on off-by-default stage(s) (trim/iqtree), comma-separated
--list-stages           print stage order and exit
--retry-failed          rerun only FAILED species
--only-species PREFIXES run only these species, comma-separated
--species-parallel N    concurrent species
```
