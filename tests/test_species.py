#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# test_species.py - Multi-species mode tests: Config path namespacing + manifest parsing.          #
#                                                                                                  #
# Each species is an independent run written into its own results tree (no mixing); the prefix is  #
# the project key and also drives the domain ids.                                                  #
#                                                                                                  #
####################################################################################################
"""

import pytest

from gwiscan.config import Config
from gwiscan.species import (
    _concurrency,
    _stage_concurrency,
    _stage_threads,
    _validate_stage_policy,
    read_species_manifest,
    write_combined_summary,
)
from gwiscan import io


# --- Config path namespacing -------------------------------------------------

def test_single_species_paths_unchanged(tmp_path):
    cfg = Config(root=tmp_path)
    assert cfg.results == tmp_path / "intermediate"
    assert cfg.final_dir == tmp_path / "final_results"
    assert cfg.logs == tmp_path / "logs"
    assert cfg.proteome == tmp_path / "input" / "proteome.fasta"
    assert cfg.proteome_db == tmp_path / "db" / "blast" / "proteome_db"
    # intermediate files are grouped into per-tool subfolders (hmm/, candidates/, ...)
    assert cfg.result("hmm_hits.tsv") == tmp_path / "intermediate" / "hmm" / "hmm_hits.tsv"


def test_multi_species_namespaced_paths(tmp_path):
    cfg = Config(root=tmp_path, SPECIES="Ath", PROTEOME="genomes/ath.fasta")
    assert cfg.results == tmp_path / "intermediate" / "Ath"
    assert cfg.logs == tmp_path / "logs" / "Ath"
    # final deliverables: top-level final_results/<Prefix>/, not nested in intermediate/
    assert cfg.final_dir == tmp_path / "final_results" / "Ath"
    assert cfg.result("hmm_hits.tsv") == tmp_path / "intermediate" / "Ath" / "hmm" / "hmm_hits.tsv"
    # per-species proteome db so species never share/overwrite each other
    assert cfg.proteome_db == tmp_path / "db" / "blast" / "Ath_proteome_db"
    # relative proteome resolved against the project root
    assert cfg.proteome == tmp_path / "genomes" / "ath.fasta"
    # shared references stay un-namespaced
    assert cfg.hmm_db == tmp_path / "db" / "hmm" / "all_models.hmm"
    assert cfg.family_map == tmp_path / "config" / "family.tsv"


def test_absolute_proteome_override(tmp_path):
    abs_path = tmp_path / "elsewhere" / "abs.fasta"
    cfg = Config(root=tmp_path, SPECIES="Gma", PROTEOME=str(abs_path))
    assert cfg.proteome == abs_path


# --- Concurrency sizing ------------------------------------------------------

def test_concurrency_explicit_cap(tmp_path):
    cfg = Config(root=tmp_path, SPECIES_PARALLEL=3)
    assert _concurrency(cfg, n_species=10) == 3      # honour the explicit cap
    assert _concurrency(cfg, n_species=2) == 2       # never exceed #species


def test_concurrency_auto_bounds_by_cores_over_threads(tmp_path, monkeypatch):
    # 0 = auto = cores // THREADS (bounded nested parallelism, no oversubscription).
    # Hold RAM large so this isolates the cores/threads math from the memory cap.
    monkeypatch.setattr("gwiscan.species.os.cpu_count", lambda: 16)
    monkeypatch.setattr("gwiscan.species._total_mem_gb", lambda: 1024.0)
    cfg = Config(root=tmp_path, SPECIES_PARALLEL=0, THREADS=4)
    assert _concurrency(cfg, n_species=10) == 4      # 16 // 4
    # capped by #species
    assert _concurrency(cfg, n_species=2) == 2


def test_concurrency_auto_never_below_one(tmp_path, monkeypatch):
    # THREADS larger than the whole machine still runs one species at a time.
    monkeypatch.setattr("gwiscan.species.os.cpu_count", lambda: 4)
    monkeypatch.setattr("gwiscan.species._total_mem_gb", lambda: 1024.0)
    cfg = Config(root=tmp_path, SPECIES_PARALLEL=0, THREADS=8)
    assert _concurrency(cfg, n_species=5) == 1


def test_concurrency_does_not_apply_deeploc_memory_to_other_stages(tmp_path, monkeypatch):
    # The scheduler applies model-RAM limits to the DeepLoc stage, not globally.
    monkeypatch.setattr("gwiscan.species.os.cpu_count", lambda: 32)
    monkeypatch.setattr("gwiscan.species._total_mem_gb", lambda: 64.0)
    cfg = Config(root=tmp_path, SPECIES_PARALLEL=0, THREADS=1, DEEPLOC_MODEL="Accurate")
    assert _concurrency(cfg, n_species=10) == 10


def test_concurrency_light_model_is_bounded_by_cores(tmp_path, monkeypatch):
    monkeypatch.setattr("gwiscan.species.os.cpu_count", lambda: 8)
    monkeypatch.setattr("gwiscan.species._total_mem_gb", lambda: 64.0)
    cfg = Config(root=tmp_path, SPECIES_PARALLEL=0, THREADS=1, DEEPLOC_MODEL="Fast")
    assert _concurrency(cfg, n_species=10) == 8


def test_concurrency_explicit_cap_is_not_memory_limited(tmp_path, monkeypatch):
    monkeypatch.setattr("gwiscan.species._total_mem_gb", lambda: 8.0)
    cfg = Config(root=tmp_path, SPECIES_PARALLEL=6, DEEPLOC_MODEL="Accurate")
    assert _concurrency(cfg, n_species=10) == 6


def test_concurrency_unknown_memory_falls_back_to_core_bound(tmp_path, monkeypatch):
    # If RAM can't be read, behave exactly as before (cores/threads only).
    monkeypatch.setattr("gwiscan.species.os.cpu_count", lambda: 16)
    monkeypatch.setattr("gwiscan.species._total_mem_gb", lambda: None)
    cfg = Config(root=tmp_path, SPECIES_PARALLEL=0, THREADS=4, DEEPLOC_MODEL="Accurate")
    assert _concurrency(cfg, n_species=10) == 4


# --- Per-stage concurrency (stage-barrier scheduling) ------------------------

def test_stage_concurrency_serialises_named_stages(tmp_path, monkeypatch):
    # Default width would allow several, but SERIAL_STAGES pins DeepTMHMM to 1.
    monkeypatch.setattr("gwiscan.species.os.cpu_count", lambda: 32)
    monkeypatch.setattr("gwiscan.species._total_mem_gb", lambda: 1024.0)
    cfg = Config(root=tmp_path, SPECIES_PARALLEL=4)
    assert _stage_concurrency(cfg, "deeptmhmm", n_alive=4) == 1   # serialised by default
    assert _stage_concurrency(cfg, "targetp", n_alive=4) == 1
    assert _stage_concurrency(cfg, "search-hmm", n_alive=4) == 4  # ordinary stage = full width


def test_stage_concurrency_explicit_override_wins(tmp_path):
    # STAGE_PARALLEL beats both SERIAL_STAGES and the default width.
    cfg = Config(root=tmp_path, SPECIES_PARALLEL=4, STAGE_PARALLEL={"deeptmhmm": 2})
    assert _stage_concurrency(cfg, "deeptmhmm", n_alive=4) == 2
    # never exceed the species still alive
    assert _stage_concurrency(cfg, "deeptmhmm", n_alive=1) == 1


def test_stage_concurrency_custom_serial_list(tmp_path):
    cfg = Config(root=tmp_path, SPECIES_PARALLEL=4, SERIAL_STAGES=["interpro"])
    assert _stage_concurrency(cfg, "interpro", n_alive=4) == 1
    # deeptmhmm no longer serialised because the list was overridden
    assert _stage_concurrency(cfg, "deeptmhmm", n_alive=4) == 4


def test_stage_concurrency_applies_accurate_deeploc_memory_cap(tmp_path, monkeypatch):
    monkeypatch.setattr("gwiscan.species.os.cpu_count", lambda: 32)
    monkeypatch.setattr("gwiscan.species._total_mem_gb", lambda: 64.0)
    cfg = Config(root=tmp_path, THREADS=1, SERIAL_STAGES=[])
    assert _stage_concurrency(cfg, "deeploc", n_alive=10) == 2
    assert _stage_concurrency(cfg, "interpro", n_alive=10) == 10


def test_stage_threads_fills_machine_when_serial(tmp_path):
    cfg = Config(root=tmp_path, THREADS=8)
    # one species at a time -> the stage gets every core
    assert _stage_threads(cfg, cores=32, conc=1) == 32
    # shared four ways -> cores/conc
    assert _stage_threads(cfg, cores=32, conc=4) == 8
    # never oversubscribes the core budget
    assert _stage_threads(cfg, cores=6, conc=1) == 6
    assert _stage_threads(cfg, cores=32, conc=8) == 4


def test_stage_policy_rejects_unknown_or_invalid_settings(tmp_path):
    selected = [("search-hmm", "HMM", None, "01_hmmscan")]
    with pytest.raises(RuntimeError, match="SERIAL_STAGES: typo"):
        _validate_stage_policy(Config(root=tmp_path, SERIAL_STAGES=["typo"]), selected)
    with pytest.raises(RuntimeError, match="positive integer"):
        _validate_stage_policy(
            Config(root=tmp_path, SERIAL_STAGES=[], STAGE_PARALLEL={"search-hmm": 0}), selected
        )


def test_stage_policy_accepts_valid_stages_outside_a_resume_slice(tmp_path):
    catalog = [
        ("targetp", "TargetP", None, "08_targetp"),
        ("deeptmhmm", "DeepTMHMM", None, "09_deeptmhmm"),
        ("deeploc", "DeepLoc", None, "10_deeploc"),
    ]
    # A resumed deeptmhmm run still carries the default targetp policy; that is
    # valid because targetp exists in the mode's full stage catalog.
    _validate_stage_policy(Config(root=tmp_path), catalog)


# --- Manifest parsing --------------------------------------------------------

def _write(tmp_path, text):
    p = tmp_path / "species.tsv"
    p.write_text(text)
    return p


def test_manifest_basic(tmp_path):
    m = _write(tmp_path, "Prefix\tProteome\nAth\tgenomes/ath.fasta\nGma\tgenomes/gma.fasta\n")
    assert read_species_manifest(m) == [
        {"prefix": "Ath", "proteome": "genomes/ath.fasta", "annotation": ""},
        {"prefix": "Gma", "proteome": "genomes/gma.fasta", "annotation": ""},
    ]


def test_manifest_optional_annotation_column(tmp_path):
    # Annotation is optional per row: a species without one is still analysed,
    # only without chromosomal coordinates.
    m = _write(tmp_path, "Prefix\tProteome\tAnnotation\n"
                         "Ath\tgenomes/ath.fasta\tgenomes/ath.gtf\n"
                         "Gma\tgenomes/gma.fasta\t\n")
    assert read_species_manifest(m) == [
        {"prefix": "Ath", "proteome": "genomes/ath.fasta", "annotation": "genomes/ath.gtf"},
        {"prefix": "Gma", "proteome": "genomes/gma.fasta", "annotation": ""},
    ]


def test_manifest_ignores_comments_and_blanks(tmp_path):
    m = _write(tmp_path, "# a note\nPrefix\tProteome\n\nAth\tath.fasta\n# skip\nGma\tgma.fasta\n")
    recs = read_species_manifest(m)
    assert [r["prefix"] for r in recs] == ["Ath", "Gma"]


def test_manifest_duplicate_prefix_rejected(tmp_path):
    m = _write(tmp_path, "Prefix\tProteome\nAth\ta.fasta\nAth\tb.fasta\n")
    with pytest.raises(RuntimeError, match="duplicate prefix"):
        read_species_manifest(m)


def test_manifest_missing_column_rejected(tmp_path):
    m = _write(tmp_path, "Name\tProteome\nAth\ta.fasta\n")
    with pytest.raises(RuntimeError, match="Prefix.*Proteome"):
        read_species_manifest(m)


def test_manifest_empty_value_rejected(tmp_path):
    m = _write(tmp_path, "Prefix\tProteome\nAth\t\n")
    with pytest.raises(RuntimeError, match="non-empty"):
        read_species_manifest(m)


# --- subset selection: --only-species / --retry-failed -----------------------

from gwiscan import species as species_mod


def test_stage_barrier_orders_stages_and_drops_failed_species(tmp_path, monkeypatch):
    """The next stage starts only after the prior one, and excludes failures."""
    manifest = tmp_path / "config" / "species.tsv"
    manifest.parent.mkdir()
    manifest.write_text("Prefix\tProteome\nAth\ta.fasta\nGma\tg.fasta\n")
    cfg = Config(root=tmp_path, SPECIES_PARALLEL=2, SERIAL_STAGES=["second"])
    selected = [
        ("first", "First", None, "01_first"),
        ("second", "Second", None, "02_second"),
    ]
    calls, widths = [], []

    class Future:
        def __init__(self, value):
            self.value = value

        def result(self):
            return self.value

    class ImmediatePool:
        def __init__(self, max_workers):
            widths.append(max_workers)

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def submit(self, func, *args):
            return Future(func(*args))

    def fake_run_stage(stage_cfg, key):
        calls.append((key, stage_cfg.SPECIES, stage_cfg.THREADS))
        if key == "first" and stage_cfg.SPECIES == "Gma":
            return "Gma", "intentional failure"
        return stage_cfg.SPECIES, None

    monkeypatch.setattr(species_mod, "ProcessPoolExecutor", ImmediatePool)
    monkeypatch.setattr(species_mod, "as_completed", lambda futures: list(futures))
    monkeypatch.setattr(species_mod.os, "cpu_count", lambda: 8)
    monkeypatch.setattr(species_mod.setupdb, "setup_shared", lambda _cfg: None)
    monkeypatch.setattr(species_mod.pipeline, "resolve_plan", lambda *_args, **_kwargs: (selected, set(), {}, set()))
    monkeypatch.setattr(species_mod.pipeline, "stages_for", lambda _cfg: selected)
    monkeypatch.setattr(species_mod, "_run_stage", fake_run_stage)
    monkeypatch.setattr(species_mod, "write_combined_summary", lambda *_args: None)

    with pytest.raises(RuntimeError, match="1 species failed: Gma"):
        species_mod.run(cfg)

    assert widths == [2, 1]
    assert calls == [("first", "Ath", 4), ("first", "Gma", 4), ("second", "Ath", 8)]
    assert species_mod._load_status(species_mod._status_path(cfg)) == {"Ath": "OK", "Gma": "FAILED"}

_MANIFEST = [
    {"prefix": "Ath", "proteome": "a.fasta", "annotation": ""},
    {"prefix": "Gma", "proteome": "g.fasta", "annotation": ""},
    {"prefix": "Vvi", "proteome": "v.fasta", "annotation": ""},
]


def test_select_species_no_filter_returns_all(tmp_path):
    cfg = Config(root=tmp_path)
    assert species_mod._select_species(cfg, _MANIFEST) == _MANIFEST


def test_select_species_only_species_subset(tmp_path):
    cfg = Config(root=tmp_path, SPECIES_ONLY=["Gma", "Vvi"])
    got = [s["prefix"] for s in species_mod._select_species(cfg, _MANIFEST)]
    assert got == ["Gma", "Vvi"]


def test_select_species_only_species_unknown_prefix_raises(tmp_path):
    cfg = Config(root=tmp_path, SPECIES_ONLY=["Nope"])
    with pytest.raises(RuntimeError, match=r"--only-species: prefix\(es\) not in .*Nope"):
        species_mod._select_species(cfg, _MANIFEST)


def test_select_species_retry_failed_reads_status(tmp_path):
    cfg = Config(root=tmp_path)
    # simulate a prior run: Ath OK, Gma FAILED, Vvi FAILED
    species_mod._write_status(species_mod._status_path(cfg),
                              {"Ath": "OK", "Gma": "FAILED", "Vvi": "FAILED"})
    cfg2 = Config(root=tmp_path, RETRY_FAILED=True)
    got = sorted(s["prefix"] for s in species_mod._select_species(cfg2, _MANIFEST))
    assert got == ["Gma", "Vvi"]


def test_select_species_retry_failed_nothing_to_retry(tmp_path, capsys):
    cfg = Config(root=tmp_path)
    species_mod._write_status(species_mod._status_path(cfg), {"Ath": "OK", "Gma": "OK"})
    cfg2 = Config(root=tmp_path, RETRY_FAILED=True)
    assert species_mod._select_species(cfg2, _MANIFEST) == []
    assert "nothing to do" in capsys.readouterr().out


def test_select_species_retry_failed_no_status_file(tmp_path, capsys):
    cfg = Config(root=tmp_path, RETRY_FAILED=True)
    assert species_mod._select_species(cfg, _MANIFEST) == []


def test_retry_failed_takes_precedence_over_only_species(tmp_path):
    cfg0 = Config(root=tmp_path)
    species_mod._write_status(species_mod._status_path(cfg0), {"Gma": "FAILED"})
    cfg = Config(root=tmp_path, RETRY_FAILED=True, SPECIES_ONLY=["Ath"])
    got = [s["prefix"] for s in species_mod._select_species(cfg, _MANIFEST)]
    assert got == ["Gma"]   # retry-failed wins


def test_status_roundtrip_and_merge(tmp_path):
    path = species_mod._status_path(Config(root=tmp_path))
    species_mod._write_status(path, {"Ath": "OK", "Gma": "FAILED"})
    loaded = species_mod._load_status(path)
    assert loaded == {"Ath": "OK", "Gma": "FAILED"}
    # merge semantics: update one, keep the rest
    loaded["Gma"] = "OK"
    loaded["Vvi"] = "FAILED"
    species_mod._write_status(path, loaded)
    assert species_mod._load_status(path) == {"Ath": "OK", "Gma": "OK", "Vvi": "FAILED"}


def test_status_load_missing_file_is_empty(tmp_path):
    assert species_mod._load_status(tmp_path / "nope.tsv") == {}


# --- Cross-species combined summary ------------------------------------------

def _write_results(cfg, prefix, rows):
    """Write a minimal final_results/<prefix>/gwiscan_results.tsv (camelCase headers)."""
    d = cfg.output_root / "final_results" / prefix
    d.mkdir(parents=True, exist_ok=True)
    io.write_tsv(d / "gwiscan_results.tsv", ["protein_id", "family"], rows)


def test_combined_summary_matrix(tmp_path):
    cfg = Config(root=tmp_path)
    _write_results(cfg, "Ath", [["a1", "GNA"], ["a2", "GNA"], ["a3", "Legume"]])
    _write_results(cfg, "Gma", [["g1", "GNA"], ["g2", "Legume"], ["g3", "Legume"]])

    write_combined_summary(cfg, ["Ath", "Gma"])

    tsv = tmp_path / "final_results" / "all_species_summary.tsv"
    lines = [ln.split("\t") for ln in tsv.read_text().splitlines()]
    header, body = lines[0], {r[0]: r for r in lines[1:]}
    # species columns kept verbatim (not camelCased), plus Family and Total
    assert header == ["Family", "Ath", "Gma", "Total"]
    assert body["GNA"] == ["GNA", "2", "1", "3"]
    assert body["Legume"] == ["Legume", "1", "2", "3"]
    assert body["Total"] == ["Total", "3", "3", "6"]
    # stacked members carry a species column
    members = io.read_tsv(tmp_path / "final_results" / "all_species_members.tsv")
    assert set(members["species"]) == {"Ath", "Gma"}
    assert len(members) == 6


def test_combined_summary_skips_missing(tmp_path):
    cfg = Config(root=tmp_path)
    _write_results(cfg, "Ath", [["a1", "GNA"]])
    # Gma has no results on disk -> skipped, not an error
    write_combined_summary(cfg, ["Ath", "Gma"])
    tsv = tmp_path / "final_results" / "all_species_summary.tsv"
    assert tsv.exists()
    assert "Gma" not in tsv.read_text().splitlines()[0]
