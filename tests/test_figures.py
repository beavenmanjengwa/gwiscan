"""Tests for the figures stage (ProtParam boxplots + stats via the bundled R script)."""

from pathlib import Path

from gwiscan import figures
from gwiscan.config import Config


def test_bundled_r_script_exists():
    # The R script must ship inside the package so it is found from any install.
    assert figures._R_SCRIPT.exists()
    assert figures._R_SCRIPT.name == "protparam_boxplots.R"


def test_run_invokes_rscript_in_final_dir(tmp_path, monkeypatch):
    cfg = Config(root=tmp_path)
    cfg.final_dir.mkdir(parents=True)
    (cfg.final_dir / "gwiscan_results.tsv").write_text("proteinId\tfamily\n")

    calls = {}
    monkeypatch.setattr(figures.external, "require", lambda b: None)
    def _run(cmd, **kw):
        calls["cmd"] = [str(c) for c in cmd]
        calls["cwd"] = kw.get("cwd")
    monkeypatch.setattr(figures.external, "run", _run)

    figures.run(cfg)

    assert calls["cmd"][0] == "Rscript"
    assert any(c.endswith("protparam_boxplots.R") for c in calls["cmd"])
    # intermediate/protparam/ is passed as the last arg (stats + outliers go there)
    assert Path(calls["cmd"][-1]) == cfg.protparam_dir
    # runs inside final_results/ so the script finds gwiscan_results.tsv there
    assert Path(calls["cwd"]) == cfg.final_dir


def test_run_skips_when_no_results(tmp_path, monkeypatch, capsys):
    cfg = Config(root=tmp_path)
    monkeypatch.setattr(figures.external, "require", lambda b: None)
    monkeypatch.setattr(figures.external, "run",
                        lambda *a, **k: (_ for _ in ()).throw(AssertionError("must not run")))
    figures.run(cfg)                       # no gwiscan_results.tsv -> warn and return
    assert "not found" in capsys.readouterr().out


def test_rscript_bin_configurable(tmp_path, monkeypatch):
    cfg = Config(root=tmp_path, RSCRIPT_BIN="/opt/R/bin/Rscript")
    cfg.final_dir.mkdir(parents=True)
    (cfg.final_dir / "gwiscan_results.tsv").write_text("x\n")
    seen = {}
    monkeypatch.setattr(figures.external, "require", lambda b: None)
    monkeypatch.setattr(figures.external, "run", lambda cmd, **kw: seen.update(bin=str(cmd[0])))
    figures.run(cfg)
    assert seen["bin"] == "/opt/R/bin/Rscript"
