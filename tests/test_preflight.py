"""Tests for preflight's family-reference-file validation.

The gap this locks: preflight used to check only that config/family.tsv EXISTS,
not that the db/blast/*.fasta and db/hmm/*.hmm files it names are present -- so a
typo or missing model FASTA sailed through preflight and only blew up hours later
inside setup-db. These tests exercise _check_family_reference_files directly.
"""

import pytest

from gwiscan import preflight
from gwiscan.config import Config


def _project(tmp_path, family_tsv):
    """A minimal project dir with a family.tsv; returns a Config rooted there."""
    (tmp_path / "config").mkdir()
    (tmp_path / "config" / "family.tsv").write_text(family_tsv)
    (tmp_path / "db" / "blast").mkdir(parents=True)
    (tmp_path / "db" / "hmm").mkdir(parents=True)
    return Config(root=tmp_path)


def _touch(cfg, relpath):
    p = cfg.root / relpath
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(">x\nACDEFGHIKL\n")
    return p


def test_all_reference_files_present_passes(tmp_path):
    cfg = _project(
        tmp_path,
        "Family\tPfamModel\tBlastModel\n"
        "GNA\tPF01453\tgna.fasta\n"
        "EUL\t-\teul.fasta\n",           # blast-only family still needs its FASTA
    )
    _touch(cfg, "db/blast/gna.fasta")
    _touch(cfg, "db/blast/eul.fasta")
    assert preflight._check_family_reference_files(cfg) is False  # False = nothing missing


def test_missing_blast_model_is_detected(tmp_path, capsys):
    cfg = _project(
        tmp_path,
        "Family\tPfamModel\tBlastModel\n"
        "GNA\tPF01453\tgna.fasta\n",
    )
    # gna.fasta deliberately NOT created
    assert preflight._check_family_reference_files(cfg) is True
    out = capsys.readouterr().out
    assert "BlastModel FASTA for 'GNA' not found" in out
    assert "gna.fasta" in out


def test_missing_custom_hmm_is_detected(tmp_path, capsys):
    cfg = _project(
        tmp_path,
        "Family\tPfamModel\tBlastModel\n"
        "CRA\tCRA.hmm\tcra.fasta\n",       # custom HMM (not a Pfam accession)
    )
    _touch(cfg, "db/blast/cra.fasta")      # blast model present...
    # ...but db/hmm/CRA.hmm is missing
    assert preflight._check_family_reference_files(cfg) is True
    out = capsys.readouterr().out
    assert "custom HMM for 'CRA' not found" in out


def test_pfam_accession_hmm_not_required_on_disk(tmp_path):
    # A Pfam-accession family's HMM is downloaded by setup-db, so its absence from
    # db/hmm/ must NOT fail preflight -- only the BlastModel FASTA is required.
    cfg = _project(
        tmp_path,
        "Family\tPfamModel\tBlastModel\n"
        "GNA\tPF01453\tgna.fasta\n",
    )
    _touch(cfg, "db/blast/gna.fasta")
    # no db/hmm/PF01453.hmm on disk
    assert preflight._check_family_reference_files(cfg) is False


def test_full_preflight_fails_when_blast_model_missing(tmp_path, monkeypatch):
    # End-to-end: run() must raise (not pass) when a referenced file is missing,
    # independent of tool availability.
    cfg = _project(
        tmp_path,
        "Family\tPfamModel\tBlastModel\n"
        "GNA\tPF01453\tgna.fasta\n",
    )
    # a proteome so that check passes and we isolate the family-file failure
    (tmp_path / "input").mkdir()
    (tmp_path / "input" / "proteome.fasta").write_text(">p1\nACDEFGHIKL\n")
    cfg = Config(root=tmp_path, EBI_EMAIL="x@y.z")
    monkeypatch.setattr(preflight.shutil, "which", lambda *_a, **_k: "/usr/bin/tool")
    monkeypatch.setattr(preflight.importlib, "import_module", lambda *_a, **_k: object())

    with pytest.raises(RuntimeError, match="pre-flight check failed"):
        preflight.run(cfg)


# --- numeric config validation (_check_config_values) ------------------------

def test_config_values_defaults_are_valid(tmp_path):
    assert preflight._check_config_values(Config(root=tmp_path)) is False


def test_config_values_identity_out_of_range(tmp_path, capsys):
    cfg = Config(root=tmp_path, DIAMOND_IDENTITY=300)
    assert preflight._check_config_values(cfg) is True
    assert "DIAMOND_IDENTITY" in capsys.readouterr().out


def test_config_values_coverage_negative(tmp_path):
    assert preflight._check_config_values(Config(root=tmp_path, DIAMOND_COVERAGE_R2=-5)) is True


def test_config_values_evalue_not_a_number(tmp_path, capsys):
    cfg = Config(root=tmp_path, DIAMOND_EVALUE="1e-5x")
    assert preflight._check_config_values(cfg) is True
    assert "DIAMOND_EVALUE" in capsys.readouterr().out


def test_config_values_evalue_valid_exponent_ok(tmp_path):
    assert preflight._check_config_values(Config(root=tmp_path, DIAMOND_EVALUE="1e-10")) is False


def test_config_values_bsr_above_one(tmp_path):
    assert preflight._check_config_values(Config(root=tmp_path, DIAMOND_BSR=1.5)) is True


def test_config_values_concordance_above_one(tmp_path):
    assert preflight._check_config_values(Config(root=tmp_path, CONCORDANCE_MIN=2.0)) is True


def test_config_values_threads_zero(tmp_path):
    assert preflight._check_config_values(Config(root=tmp_path, THREADS=0)) is True


def test_config_values_bootstrap_zero_is_allowed(tmp_path):
    # 0 = no bootstrap, a valid choice
    assert preflight._check_config_values(Config(root=tmp_path, IQTREE_BOOTSTRAP=0)) is False


def test_config_values_nmotifs_zero_invalid(tmp_path):
    assert preflight._check_config_values(Config(root=tmp_path, MEME_NMOTIFS=0)) is True


def test_config_values_reports_multiple_problems(tmp_path, capsys):
    cfg = Config(root=tmp_path, DIAMOND_IDENTITY=300, DIAMOND_BSR=9, THREADS=0)
    assert preflight._check_config_values(cfg) is True
    out = capsys.readouterr().out
    assert "DIAMOND_IDENTITY" in out and "DIAMOND_BSR" in out and "THREADS" in out
