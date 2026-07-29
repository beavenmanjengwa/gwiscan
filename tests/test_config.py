"""Tests for configuration precedence and path layout.

Settings use one UPPER_CASE name across the config.yaml key, the environment
variable, and the Config field; derived paths stay snake_case properties.
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
