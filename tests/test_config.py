#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# test_config.py - Configuration precedence and path-layout tests.                                 #
#                                                                                                  #
# Settings use one UPPER_CASE name across the config.yaml key, the environment variable, and the   #
# Config field; derived paths stay snake_case properties.                                          #
#                                                                                                  #
####################################################################################################
"""

import textwrap

from gwiscan.config import Config


def test_defaults(tmp_path):
    cfg = Config.load(root=tmp_path)
    assert cfg.THREADS == 4
    assert cfg.DIAMOND_EVALUE == "1e-5"
    assert cfg.PRIMARY_TRANSCRIPT is True
    assert cfg.proteome == tmp_path / "input" / "proteome.fasta"
    assert cfg.hmm_db == tmp_path / "db" / "hmm" / "all_models.hmm"


def test_config_file_then_override(tmp_path):
    (tmp_path / "config.yaml").write_text(textwrap.dedent("""
        THREADS: 16
        DIAMOND_EVALUE: "1e-10"
        INTERPRO_APPL: "Pfam,SMART"
    """))
    cfg = Config.load(root=tmp_path)
    assert cfg.THREADS == 16
    assert cfg.DIAMOND_EVALUE == "1e-10"
    assert cfg.INTERPRO_APPL == "Pfam,SMART"

    # CLI overrides win over the config file; None values are ignored.
    cfg2 = Config.load(root=tmp_path, overrides={"THREADS": 32, "DIAMOND_EVALUE": None})
    assert cfg2.THREADS == 32
    assert cfg2.DIAMOND_EVALUE == "1e-10"


def test_output_defaults_to_project_dir(tmp_path):
    # No OUTPUT -> intermediate/logs sit under the project directory; the final
    # deliverables land in a top-level final_results/ (not nested inside).
    cfg = Config.load(root=tmp_path)
    assert cfg.results == tmp_path / "intermediate"
    assert cfg.logs == tmp_path / "logs"
    assert cfg.final_dir == tmp_path / "final_results"


def test_output_dir_separates_results_from_inputs(tmp_path):
    out = tmp_path / "elsewhere"
    cfg = Config.load(root=tmp_path, overrides={"OUTPUT": str(out)})
    # intermediate/final_results/logs go to OUTPUT...
    assert cfg.results == out / "intermediate"
    assert cfg.logs == out / "logs"
    assert cfg.final_dir == out / "final_results"    # top-level, sibling of intermediate/
    # ...but inputs, config, and db stay under the project directory.
    assert cfg.proteome == tmp_path / "input" / "proteome.fasta"
    assert cfg.config_dir == tmp_path / "config"
    assert cfg.db_dir == tmp_path / "db"


def test_output_dir_multispecies_namespaced(tmp_path):
    out = tmp_path / "out"
    cfg = Config.load(root=tmp_path, overrides={"OUTPUT": str(out), "SPECIES": "Ath"})
    assert cfg.results == out / "intermediate" / "Ath"
    assert cfg.final_dir == out / "final_results" / "Ath"
    assert cfg.logs == out / "logs" / "Ath"


def test_env_override(tmp_path, monkeypatch):
    monkeypatch.setenv("THREADS", "8")
    monkeypatch.setenv("EBI_EMAIL", "me@example.org")
    monkeypatch.setenv("PRIMARY_TRANSCRIPT", "false")
    cfg = Config.load(root=tmp_path)
    assert cfg.THREADS == 8
    assert cfg.EBI_EMAIL == "me@example.org"
    assert cfg.PRIMARY_TRANSCRIPT is False


def test_optional_tool_bins_default_to_path_names(tmp_path):
    # WebLogo/MEME are overridable like TARGETP_BIN/DEEPLOC_BIN/IQTREE_BIN; the
    # defaults are the bare PATH names.
    cfg = Config.load(root=tmp_path)
    assert cfg.WEBLOGO_BIN == "weblogo"
    assert cfg.MEME_BIN == "meme"


def test_optional_tool_bins_overridable(tmp_path, monkeypatch):
    # config.yaml, env, and CLI override all resolve, same as any other setting.
    (tmp_path / "config.yaml").write_text("WEBLOGO_BIN: /envs/x/bin/weblogo\n")
    monkeypatch.setenv("MEME_BIN", "/opt/meme/bin/meme")
    cfg = Config.load(root=tmp_path)
    assert cfg.WEBLOGO_BIN == "/envs/x/bin/weblogo"
    assert cfg.MEME_BIN == "/opt/meme/bin/meme"

    cfg2 = Config.load(root=tmp_path, overrides={"WEBLOGO_BIN": "/cli/weblogo"})
    assert cfg2.WEBLOGO_BIN == "/cli/weblogo"   # CLI wins over config.yaml


def test_unknown_config_key_warns_and_is_ignored(tmp_path, capsys):
    # A misspelled key must not be applied silently with the default; the user is
    # told which settings were unrecognised.
    (tmp_path / "config.yaml").write_text("THREDS: 16\nDIAMOND_EVALUE: \"1e-8\"\n")
    cfg = Config.load(root=tmp_path)
    out = capsys.readouterr().out
    assert "unknown setting" in out
    assert "THREDS" in out
    assert cfg.THREADS == 4              # typo ignored, default kept
    assert cfg.DIAMOND_EVALUE == "1e-8"  # valid key still applied


def test_all_known_keys_produce_no_warning(tmp_path, capsys):
    (tmp_path / "config.yaml").write_text("THREADS: 8\nMODE: multi-family\n")
    Config.load(root=tmp_path)
    assert "unknown setting" not in capsys.readouterr().out


def test_prefixed_env_var_is_read(tmp_path, monkeypatch):
    monkeypatch.setenv("GWISCAN_THREADS", "12")
    cfg = Config.load(root=tmp_path)
    assert cfg.THREADS == 12


def test_prefixed_env_var_wins_over_bare(tmp_path, monkeypatch):
    monkeypatch.setenv("GWISCAN_THREADS", "12")
    monkeypatch.setenv("THREADS", "3")
    cfg = Config.load(root=tmp_path)
    assert cfg.THREADS == 12


def test_bare_generic_env_var_still_works_but_warns(tmp_path, monkeypatch, capsys):
    # Back-compat: the bare name is honoured, with a heads-up because it is generic.
    monkeypatch.delenv("GWISCAN_THREADS", raising=False)
    monkeypatch.setenv("THREADS", "7")
    cfg = Config.load(root=tmp_path)
    assert cfg.THREADS == 7
    out = capsys.readouterr().out
    assert "GWISCAN_THREADS" in out and "generic" in out


def test_bare_specific_env_var_is_silent(tmp_path, monkeypatch, capsys):
    # A gwiscan-specific name is unlikely to clash, so no warning for the bare form.
    monkeypatch.setenv("DIAMOND_EVALUE", "1e-12")
    cfg = Config.load(root=tmp_path)
    assert cfg.DIAMOND_EVALUE == "1e-12"
    assert "generic" not in capsys.readouterr().out
