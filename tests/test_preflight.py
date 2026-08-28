#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# test_preflight.py - preflight family-reference-file validation tests.                            #
#                                                                                                  #
# preflight checks not only that config/family.tsv exists but that the db/blast/*.fasta and        #
# db/hmm/*.hmm files it names are present, so a typo or a missing model FASTA is caught up front   #
# rather than inside setup-db. These tests exercise _check_family_reference_files directly.        #
#                                                                                                  #
####################################################################################################
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


def test_check_proteomes_multispecies_all_present(tmp_path):
    # Multi-species: proteomes come from species.tsv, NOT input/proteome.fasta.
    (tmp_path / "config").mkdir()
    (tmp_path / "config" / "species.tsv").write_text(
        "Prefix\tProteome\nAth\tinput/a.fa\nGma\tinput/b.fa\n")
    (tmp_path / "input").mkdir()
    (tmp_path / "input" / "a.fa").write_text(">p1\nACDE\n")
    (tmp_path / "input" / "b.fa").write_text(">p2\nFGHI\n")
    cfg = Config(root=tmp_path)
    assert preflight._check_proteomes(cfg) is False   # no single input/proteome.fasta needed


def test_check_proteomes_multispecies_reports_missing(tmp_path, capsys):
    (tmp_path / "config").mkdir()
    (tmp_path / "config" / "species.tsv").write_text(
        "Prefix\tProteome\nAth\tinput/a.fa\nGma\tinput/missing.fa\n")
    (tmp_path / "input").mkdir()
    (tmp_path / "input" / "a.fa").write_text(">p1\nACDE\n")
    cfg = Config(root=tmp_path)
    assert preflight._check_proteomes(cfg) is True
    out = capsys.readouterr().out
    assert "Gma" in out and "missing.fa" in out


def test_check_proteomes_single_species(tmp_path):
    (tmp_path / "input").mkdir()
    (tmp_path / "input" / "proteome.fasta").write_text(">p1\nACDE\n>p2\nFGHI\n")
    cfg = Config(root=tmp_path)   # no config/species.tsv -> single-species path
    assert preflight._check_proteomes(cfg) is False


def test_duplicate_fasta_ids_found(tmp_path):
    # Duplicate ids would crash the proteome index (SeqIO.to_dict) mid-run.
    p = tmp_path / "proteome.fasta"
    p.write_text(">P1 kinase\nACDE\n>P2\nFGHI\n>P1 other\nKLMN\n>P3\nPQRS\n>P2\nTVWY\n")
    assert preflight._duplicate_fasta_ids(p) == ["P1", "P2"]


def test_duplicate_fasta_ids_none_when_unique(tmp_path):
    p = tmp_path / "proteome.fasta"
    p.write_text(">P1\nACDE\n>P2 desc\nFGHI\n>P3\nKLMN\n")
    assert preflight._duplicate_fasta_ids(p) == []


_HMM_WITH_GA = (
    "HMMER3/f [3.4]\nNAME  CRA\nLENG  100\nGA    25.00 25.00;\nHMM  A  C\n//\n"
)
_HMM_NO_GA = "HMMER3/f [3.4]\nNAME  CRA\nLENG  100\nHMM  A  C\n//\n"


def test_custom_hmm_without_ga_thresholds_is_detected(tmp_path, capsys):
    # An identifying custom HMM lacking GA cutoffs would make hmmsearch --cut_ga
    # abort mid-run; preflight must catch it up front.
    cfg = _project(
        tmp_path,
        "Family\tPfamModel\tBlastModel\n"
        "CRA\tCRA.hmm\tcra.fasta\n",
    )
    _touch(cfg, "db/blast/cra.fasta")
    (cfg.root / "db" / "hmm" / "CRA.hmm").write_text(_HMM_NO_GA)
    assert preflight._check_family_reference_files(cfg) is True
    out = capsys.readouterr().out
    assert "no GA (gathering) thresholds" in out
    assert "CRA" in out


def test_custom_hmm_with_ga_thresholds_passes(tmp_path):
    cfg = _project(
        tmp_path,
        "Family\tPfamModel\tBlastModel\n"
        "CRA\tCRA.hmm\tcra.fasta\n",
    )
    _touch(cfg, "db/blast/cra.fasta")
    (cfg.root / "db" / "hmm" / "CRA.hmm").write_text(_HMM_WITH_GA)
    assert preflight._check_family_reference_files(cfg) is False


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


def test_config_values_deeptmhmm_mode_invalid(tmp_path):
    assert preflight._check_config_values(Config(root=tmp_path, DEEPTMHMM_MODE="loacl")) is True


def test_config_values_reports_multiple_problems(tmp_path, capsys):
    cfg = Config(root=tmp_path, DIAMOND_IDENTITY=300, DIAMOND_BSR=9, THREADS=0)
    assert preflight._check_config_values(cfg) is True
    out = capsys.readouterr().out
    assert "DIAMOND_IDENTITY" in out and "DIAMOND_BSR" in out and "THREADS" in out
