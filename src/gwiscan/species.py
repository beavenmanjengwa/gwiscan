#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# species.py - Multi-species driver: run each proteome as an independent parallel job.             #
#                                                                                                  #
# Reads config/species.tsv (Prefix, Proteome), presses the shared databases once, then schedules   #
# each pipeline stage across species. Each species writes into intermediate/<Prefix>/ and is fully #
# isolated; one species failing does not abort the others.                                         #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

import csv
import dataclasses
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from . import external, io, pipeline, setupdb
from .config import Config

# Per-species OK/FAILED outcomes are persisted here after every multi-species run
# so `--retry-failed` can re-run just the ones that failed, and so the record
# survives the process exiting. Written under the run's log dir (not namespaced to
# any single species).
_STATUS_FILENAME = "species_status.tsv"


def _status_path(cfg: Config) -> Path:
    return cfg.output_root / "logs" / _STATUS_FILENAME


def _load_status(path: Path) -> dict:
    """{prefix -> 'OK'|'FAILED'} from a prior run, or {} if none/unreadable."""
    if not path.exists():
        return {}
    status = {}
    for row in path.read_text().splitlines():
        parts = row.split("\t")
        if len(parts) == 2 and parts[0] and parts[0] != "Prefix":
            status[parts[0]] = parts[1]
    return status


def _write_status(path: Path, status: dict) -> None:
    """Persist {prefix -> status}, merging is the caller's job (this overwrites)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["Prefix\tStatus"]
    lines += [f"{prefix}\t{state}" for prefix, state in sorted(status.items())]
    path.write_text("\n".join(lines) + "\n")


def _select_species(cfg: Config, species: list) -> list:
    """Apply --retry-failed / --only-species to the manifest rows. Returns the
    subset to run (all of them if neither is set). Raises if a named prefix isn't
    in the manifest; returns [] (with a message) if --retry-failed has nothing to
    retry."""
    all_prefixes = {s["prefix"] for s in species}

    if cfg.RETRY_FAILED:
        prior = _load_status(_status_path(cfg))
        wanted = {p for p, state in prior.items() if state == "FAILED"}
        if not wanted:
            external.log("[multi] --retry-failed: no previously FAILED species recorded; nothing to do.")
            return []
        flag = "--retry-failed"
    elif cfg.SPECIES_ONLY:
        wanted = set(cfg.SPECIES_ONLY)
        flag = "--only-species"
    else:
        return species

    unknown = wanted - all_prefixes
    if unknown:
        raise RuntimeError(
            f"{flag}: prefix(es) not in {cfg.species_manifest.name}: {', '.join(sorted(unknown))}"
        )
    subset = [s for s in species if s["prefix"] in wanted]
    external.log(f"[multi] {flag}: running {len(subset)} of {len(species)} species: "
                 f"{', '.join(s['prefix'] for s in subset)}")
    return subset


def _total_mem_gb():
    """Total physical RAM in GB, or None if it can't be read portably (e.g. on a
    platform without these sysconf names). Used only to keep the auto concurrency
    from oversubscribing memory."""
    try:
        return os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE") / 1e9
    except (ValueError, AttributeError, OSError):
        return None


def _mem_per_species_gb(cfg: Config, key: str) -> float | None:
    """Peak RAM estimate for a particular heavy stage, or None if it is not
    resource-limited by a known in-memory model."""
    if key == "deeploc" and str(cfg.DEEPLOC_MODEL).lower() == "accurate":
        return 32.0
    return None


def _concurrency(cfg: Config, n_species: int) -> int:
    """How many species to run at once.

    SPECIES_PARALLEL > 0 forces that many. 0 (default) auto-sizes to cores // THREADS
    so the outer (species) x inner (per-species THREADS) parallelism totals ~one
    thread per core -- full CPU use without oversubscription, and without launching
    every species at once (which can exhaust RAM or hammer the InterProScan API).

    Never more than the number of species, never less than 1. Stage-specific RAM
    caps are applied by _stage_concurrency, where the actual tool is known.
    """
    if cfg.SPECIES_PARALLEL > 0:
        return min(cfg.SPECIES_PARALLEL, n_species)
    cores = os.cpu_count() or 1
    auto = max(1, cores // max(1, cfg.THREADS))

    return min(auto, n_species)


def read_species_manifest(path) -> list:
    """Parse config/species.tsv into [{'prefix', 'proteome', 'annotation'}, ...].

    Accepts Prefix, Proteome and an optional Annotation (case-insensitive headers);
    blank lines and lines starting with '#' are ignored. Annotation is that species'
    representative GTF/GFF3, giving its members chromosomal coordinates. Raises on
    missing columns, empty values, or duplicate prefixes (which would collide on
    intermediate/<Prefix>/).
    """
    with open(path) as fh:
        rows = [ln for ln in fh if ln.strip() and not ln.lstrip().startswith("#")]
    reader = csv.DictReader(rows, delimiter="\t")
    lower = {(c or "").strip().lower(): c for c in (reader.fieldnames or [])}
    if "prefix" not in lower or "proteome" not in lower:
        raise RuntimeError(
            f"{path}: species manifest needs 'Prefix' and 'Proteome' columns "
            f"(found: {reader.fieldnames})"
        )

    out, seen = [], set()
    for row in reader:
        prefix = (row[lower["prefix"]] or "").strip()
        proteome = (row[lower["proteome"]] or "").strip()
        if not prefix or not proteome:
            raise RuntimeError(f"{path}: every row needs a non-empty Prefix and Proteome (got {row})")
        if prefix in seen:
            raise RuntimeError(f"{path}: duplicate prefix '{prefix}' (each must be unique)")
        seen.add(prefix)
        annotation = (row.get(lower.get("annotation", ""), "") or "").strip()
        out.append({"prefix": prefix, "proteome": proteome, "annotation": annotation})
    if not out:
        raise RuntimeError(f"{path}: no species rows found")
    return out


def _species_config(cfg: Config, prefix: str, proteome: str, annotation: str = "") -> Config:
    """A per-species Config: namespaces outputs to intermediate/<prefix>/ and points at
    that species' proteome and annotation; the prefix also drives the domain ids."""
    return dataclasses.replace(cfg, SPECIES=prefix, SPECIES_PREFIX=prefix,
                               PROTEOME=proteome, ANNOTATION=annotation)


def _configured_families(cfg: Config) -> list:
    """Every family the study defines, so the cross-species matrix has a row for each
    even when a species has zero of it. Architecture names in MODE: architecture,
    else the family/multi-family table's families."""
    try:
        if cfg.is_architecture:
            from . import architecture
            return architecture.arch_names(architecture.read_rules(cfg.architecture_map))
        return [r["family"] for r in io.family_records(cfg.family_map)]
    except (OSError, ValueError, KeyError):
        return []


def write_combined_summary(cfg: Config, prefixes: list) -> None:
    """One cross-species summary over the per-species results.

    Reads each species' final_results/<prefix>/gwiscan_results.tsv and writes, to the
    top-level final_results/:
      * all_species_summary.tsv  -- families (rows) x species (columns) matrix of the
        member protein count, with per-family and per-species Totals; every configured
        family has a row even at zero.
      * all_species_summary.xlsx -- that matrix, a Multifamily matrix (multi-family mode),
        and every species' members stacked into one Members sheet (a Species column added).
      * all_species_members.tsv  -- the stacked members as TSV.
    Species without a results file (never ran / failed early) are skipped.
    """
    import pandas as pd

    final_root = cfg.output_root / "final_results"
    frames, found = [], []
    for prefix in prefixes:
        res = final_root / prefix / "gwiscan_results.tsv"
        if not res.exists():
            continue
        df = io.read_tsv(res, low_memory=False)
        if "family" not in df.columns or "protein_id" not in df.columns:
            continue
        df.insert(0, "species", prefix)
        frames.append(df)
        found.append(prefix)

    if not frames:
        external.log("[multi] No per-species results found; skipping combined summary.")
        return

    combined = pd.concat(frames, ignore_index=True)

    def _matrix(group_col: str) -> "pd.DataFrame":
        m = (combined.groupby([group_col, "species"])["protein_id"].nunique()
             .unstack("species").reindex(columns=found).fillna(0).astype(int))
        if group_col == "family":                       # show every configured family
            m = m.reindex(sorted(set(m.index) | set(_configured_families(cfg))), fill_value=0)
        else:
            m = m.sort_index()
        m["Total"] = m.sum(axis=1)
        m.loc["Total"] = m.sum(axis=0)
        return m.reset_index().rename(columns={group_col: group_col.capitalize()})

    fam_matrix = _matrix("family")
    out_tsv = final_root / "all_species_summary.tsv"
    out_xlsx = final_root / "all_species_summary.xlsx"
    out_members = final_root / "all_species_members.tsv"
    final_root.mkdir(parents=True, exist_ok=True)

    # The matrix columns are species prefixes (Ath, Csa, ...) and Family/Total, so
    # they are written verbatim -- NOT camelCased. The stacked members table keeps
    # the same camelCase headers as the per-species gwiscan_results.tsv.
    fam_matrix.to_csv(out_tsv, sep="\t", index=False)
    io.write_df(combined, out_members, "tsv")

    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        fam_matrix.to_excel(writer, sheet_name="Family_matrix", index=False)
        if "multifamily" in combined.columns:
            _matrix("multifamily").to_excel(writer, sheet_name="Multifamily_matrix", index=False)
        combined.rename(columns=io.to_camel).to_excel(writer, sheet_name="Members", index=False)

    external.log(
        f"[multi] Combined summary over {len(found)} species ({', '.join(found)}): "
        f"{combined['protein_id'].nunique()} members -> {out_tsv.name}, "
        f"{out_xlsx.name}, {out_members.name}"
    )


def _stage_concurrency(cfg: Config, key: str, n_alive: int) -> int:
    """How many species run stage ``key`` at once.

    Use STAGE_PARALLEL[key] if it is set; otherwise run one species at a time when
    the stage is listed in SERIAL_STAGES; otherwise use the default width from
    _concurrency (SPECIES_PARALLEL, or auto from cores/threads). Capped at the
    species still alive, and at least 1."""
    override = (cfg.STAGE_PARALLEL or {}).get(key)
    if override is not None and int(override) > 0:
        want = int(override)
    elif key in set(cfg.SERIAL_STAGES or ()):
        want = 1
    else:
        want = _concurrency(cfg, n_alive)
        per_species_gb = _mem_per_species_gb(cfg, key)
        mem_gb = _total_mem_gb()
        if per_species_gb and mem_gb is not None:
            want = min(want, max(1, int(mem_gb // per_species_gb)))
    return max(1, min(want, n_alive))


def _stage_threads(cfg: Config, cores: int, conc: int) -> int:
    """Per-species thread budget for a stage running ``conc`` species at once.

    It is exactly the machine budget divided by the number of concurrent species,
    so concurrent tools cannot collectively exceed the available cores. A
    serialised stage receives every core."""
    return max(1, cores // max(1, conc))


def _validate_stage_policy(cfg: Config, available_stages: list) -> None:
    """Fail before work when policy names or widths are invalid for this MODE.

    ``available_stages`` must be the complete mode-specific catalog, not the
    current --from-stage/--until slice: a policy for an earlier or later stage is
    still valid during a resumed run.
    """
    valid = {key for key, *_ in available_stages}
    unknown_serial = sorted(set(cfg.SERIAL_STAGES or ()) - valid)
    unknown_parallel = sorted(set((cfg.STAGE_PARALLEL or {})) - valid)
    if unknown_serial or unknown_parallel:
        details = []
        if unknown_serial:
            details.append(f"SERIAL_STAGES: {', '.join(unknown_serial)}")
        if unknown_parallel:
            details.append(f"STAGE_PARALLEL: {', '.join(unknown_parallel)}")
        raise RuntimeError(
            f"unknown stage policy key(s) ({'; '.join(details)}). "
            f"Valid stages: {', '.join(key for key, *_ in available_stages)}"
        )
    for key, value in (cfg.STAGE_PARALLEL or {}).items():
        try:
            width = int(value)
        except (TypeError, ValueError):
            raise RuntimeError(f"STAGE_PARALLEL[{key!r}] must be a positive integer, got {value!r}")
        if width < 1:
            raise RuntimeError(f"STAGE_PARALLEL[{key!r}] must be a positive integer, got {value!r}")


def _lookup_stage(cfg: Config, key: str) -> tuple:
    """The (label, func, logname) for stage ``key`` in this run's MODE. Resolved
    from pipeline.stages_for INSIDE the worker so the stage's lambda is rebuilt
    locally and never has to be pickled across the process boundary."""
    for k, label, func, logname in pipeline.stages_for(cfg):
        if k == key:
            return label, func, logname
    raise KeyError(f"stage {key!r} not found for this run's mode")


def _run_stage(cfg: Config, key: str) -> tuple:
    """Worker: run ONE stage (named by ``key``) for one species (shared setup
    already done). Returns (prefix, None) on success or (prefix, error_message) on
    failure, so a species failing a stage never aborts the others -- it is dropped
    from the remaining stages instead.

    Only the picklable Config and the stage key cross the process boundary; the
    stage function (a lambda) is looked up locally via _lookup_stage."""
    # This worker runs in its own process but shares the terminal with the other
    # species; tag every line so interleaved output stays readable. (Each species'
    # authoritative per-stage log still lives in its own logs/<prefix>/ tree.)
    external.set_line_prefix(f"[{cfg.SPECIES}] ")
    try:
        label, func, logname = _lookup_stage(cfg, key)
        pipeline._stage(cfg, label, func, logname)
        return (cfg.SPECIES, None)
    except Exception as e:  # noqa: BLE001 - report per-species, keep the rest running
        return (cfg.SPECIES, f"{type(e).__name__}: {e}")
    finally:
        external.set_line_prefix("")


def run(cfg: Config) -> None:
    """Run config/species.tsv with stage-barrier scheduling: every species clears a
    stage before any moves on, so per-stage concurrency can differ. Most stages run
    the species in parallel; SERIAL_STAGES / STAGE_PARALLEL pin resource-bound stages
    (DeepTMHMM, TargetP, DeepLoc) to fewer species at once, and the per-species thread
    budget is rebalanced so a serialised stage still uses the whole machine.

    --only-species / --retry-failed narrow it to a subset; every run's per-species
    OK/FAILED outcomes are persisted to logs/species_status.tsv (feeding
    --retry-failed). A species that fails one stage is dropped from later stages but
    keeps whatever it already produced (so the combined summary still includes it)."""
    manifest = read_species_manifest(cfg.species_manifest)
    external.log(f"[multi] {len(manifest)} species: {', '.join(s['prefix'] for s in manifest)}")

    species = _select_species(cfg, manifest)
    if not species:
        return   # --retry-failed with nothing to retry; already logged

    # Shared, once: press the HMM db + verify BLAST model FASTAs (read-only refs).
    external.log("[multi] Building shared databases (HMM + model FASTAs)...")
    setupdb.setup_shared(cfg)

    configs = {s["prefix"]: _species_config(cfg, s["prefix"], s["proteome"], s.get("annotation", ""))
               for s in species}
    # Create each species' intermediate/<prefix>/ and logs/<prefix>/ up front:
    # pipeline._stage opens the per-stage log before the stage runs, so the tree
    # must exist even though the stage functions also ensure their own dirs.
    for c in configs.values():
        c.ensure_dirs()

    # The plan (stage list + skip/auto-skip) is identical across species -- same
    # MODE and same --from-stage/--until/--skip -- so resolve it once from any
    # config. Shared setup was just done, so it is excluded here.
    sample = next(iter(configs.values()))
    selected, skip, auto_skip, off = pipeline.resolve_plan(sample, include_shared_setup=False)
    available_stages = [s for s in pipeline.stages_for(sample) if s[0] != "setup-shared"]
    _validate_stage_policy(cfg, available_stages)

    cores = os.cpu_count() or 1
    serial = ", ".join(cfg.SERIAL_STAGES) or "none"
    external.log(
        f"[multi] {len(configs)} species x {len(selected)} stage(s), stage-barrier "
        f"scheduling (cores={cores}, default width {_concurrency(cfg, len(configs))}). "
        f"Serialised stages: {serial}."
    )
    if cfg.FROM_STAGE:
        external.log(f"[multi] [resume] starting from stage '{cfg.FROM_STAGE}'.")
    if cfg.UNTIL_STAGE:
        external.log(f"[multi] [resume] stopping after stage '{cfg.UNTIL_STAGE}'.")

    alive = dict(configs)     # prefix -> cfg for species still running
    failed = {}               # prefix -> first error message

    for key, label, _func, _logname in selected:
        notice = pipeline.skip_notice(key, label, skip, auto_skip, off)
        if notice:
            external.log("\n[multi] " + notice)
            continue
        if not alive:
            break

        conc = _stage_concurrency(cfg, key, len(alive))
        threads = _stage_threads(cfg, cores, conc)
        external.log(
            f"\n[multi] {label} -- {len(alive)} species, {conc} at once, {threads} threads each"
        )
        # Rebalance the per-species thread budget for this stage so it fills the
        # machine (a serialised stage gets every core). Each species still writes
        # into its own intermediate/<prefix>/ tree, staying independent.
        stage_cfgs = {p: dataclasses.replace(c, THREADS=threads) for p, c in alive.items()}
        stage_failures = {}
        with ProcessPoolExecutor(max_workers=conc) as pool:
            futures = {pool.submit(_run_stage, c, key): p
                       for p, c in stage_cfgs.items()}
            for fut in as_completed(futures):
                prefix, err = fut.result()
                if err:
                    stage_failures[prefix] = err
        # Drop species that failed this stage; they keep their earlier outputs.
        for prefix, err in stage_failures.items():
            failed[prefix] = err
            alive.pop(prefix, None)
            external.log(f"[multi] [FAIL] {prefix} at stage '{key}': {err}")

    for prefix in alive:
        external.log(f"[multi] [DONE] {prefix} -> {cfg.output_root / 'final_results' / prefix}")

    # Persist outcomes, MERGING with any prior status so species not in this run
    # (e.g. when only a subset ran) keep their previous OK/FAILED state -- that is
    # what makes a later --retry-failed see the right set.
    status = _load_status(_status_path(cfg))
    for prefix in configs:
        status[prefix] = "FAILED" if prefix in failed else "OK"
    _write_status(_status_path(cfg), status)

    # Cross-species summary over every manifest species whose results are on disk
    # (a subset run, or a species that failed only in a late optional stage, still
    # contributes its finished results). Runs before the failure-raise so it is
    # always produced.
    try:
        write_combined_summary(cfg, [s["prefix"] for s in manifest])
    except Exception as e:  # noqa: BLE001 - the summary is a convenience, never fatal
        external.log(f"[multi] [WARN] could not write combined summary: {type(e).__name__}: {e}")

    external.log("=" * 56)
    external.log(f"[multi] Complete: {len(alive)}/{len(configs)} species succeeded")
    external.log(f"[multi] Status written: {_status_path(cfg)}")
    if failed:
        ordered = [p for p in configs if p in failed]
        external.log(f"[multi] Failed: {', '.join(ordered)}")
        external.log("[multi] Re-run just these with:  gwiscan run --retry-failed")
        raise RuntimeError(f"{len(failed)} species failed: {', '.join(ordered)}")
