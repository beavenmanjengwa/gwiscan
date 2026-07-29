"""Tests for external.py -- the subprocess wrapper every stage funnels through.

A regression here (wrong check field, swallowed stderr, bad stringification)
would silently affect every external-tool stage at once, so the wrapper's own
behaviour is worth locking independently of any one stage.
"""

import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

from gwiscan import external


def _fake_run(monkeypatch, returncode=0, stdout="out", stderr="err"):
    """Patch subprocess.run to record its call and return a canned result."""
    seen = {}

    def _run(cmd, **kwargs):
        seen["cmd"] = cmd
        seen["kwargs"] = kwargs
        return SimpleNamespace(returncode=returncode, stdout=stdout, stderr=stderr)

    monkeypatch.setattr(external.subprocess, "run", _run)
    return seen


def test_run_stringifies_command_elements(monkeypatch, capsys):
    seen = _fake_run(monkeypatch)
    external.run(["tool", Path("/a/b"), 42])
    # every element passed to subprocess.run is a str (Paths/ints coerced)
    assert seen["cmd"] == ["tool", "/a/b", "42"]
    assert all(isinstance(c, str) for c in seen["cmd"])


def test_run_merges_stdout_stderr_when_no_stdout_path(monkeypatch, capsys):
    _fake_run(monkeypatch, stdout="hello-merged")
    external.run(["tool"])
    # merged mode: stdout+stderr combined onto our stdout, echoed
    assert "hello-merged" in capsys.readouterr().out


def test_run_routes_stdout_to_file_and_echoes_only_stderr(monkeypatch, capsys, tmp_path):
    seen = _fake_run(monkeypatch, stdout="TO_FILE", stderr="ONLY_STDERR")
    out_file = tmp_path / "out.txt"
    external.run(["tool"], stdout_path=out_file)
    # with stdout_path, subprocess.run is told to write stdout to the file handle
    assert "stdout" in seen["kwargs"]
    # and only stderr is echoed to our stdout
    printed = capsys.readouterr().out
    assert "ONLY_STDERR" in printed
    assert "TO_FILE" not in printed


def test_run_raises_on_nonzero_with_command_in_message(monkeypatch, capsys):
    _fake_run(monkeypatch, returncode=2)
    with pytest.raises(RuntimeError, match=r"exit 2.*failing-tool"):
        external.run(["failing-tool", "--flag"])


def test_run_check_false_does_not_raise(monkeypatch, capsys):
    _fake_run(monkeypatch, returncode=2)
    result = external.run(["failing-tool"], check=False)  # must not raise
    assert result.returncode == 2


def test_available_true_for_path_name(monkeypatch):
    monkeypatch.setattr(external.shutil, "which", lambda b: "/usr/bin/" + b)
    assert external.available("mafft") is True


def test_available_true_for_existing_absolute_path(monkeypatch, tmp_path):
    monkeypatch.setattr(external.shutil, "which", lambda b: None)  # not on PATH
    real = tmp_path / "weblogo"
    real.write_text("#!/bin/sh\n")
    assert external.available(str(real)) is True          # but exists as a file


def test_available_false_when_neither(monkeypatch):
    monkeypatch.setattr(external.shutil, "which", lambda b: None)
    assert external.available("/nope/not/here") is False


def test_require_raises_when_unavailable(monkeypatch):
    monkeypatch.setattr(external.shutil, "which", lambda b: None)
    with pytest.raises(FileNotFoundError, match="required tool not found on PATH: ghost"):
        external.require("ghost")


def test_require_ok_when_available(monkeypatch):
    monkeypatch.setattr(external.shutil, "which", lambda b: "/usr/bin/x")
    external.require("x")   # must not raise


def test_run_real_subprocess_roundtrip():
    # One un-mocked call to catch signature drift against the real subprocess API.
    result = external.run(["true"])
    assert result.returncode == 0
    with pytest.raises(RuntimeError):
        external.run(["false"])
    # sanity: subprocess is really the stdlib module
    assert external.subprocess is subprocess


# --- per-process log line prefix (multi-species interleaving) -----------------

def test_log_prefix_applied_and_cleared(capsys):
    external.set_line_prefix("[Ath] ")
    try:
        external.log("HMMscan done")
        assert capsys.readouterr().out == "[Ath] HMMscan done\n"
    finally:
        external.set_line_prefix("")
    external.log("no prefix")
    assert capsys.readouterr().out == "no prefix\n"


def test_log_prefix_tags_each_line_but_not_blanks(capsys):
    external.set_line_prefix("[Gma] ")
    try:
        external.log("line1\n\nline2")
    finally:
        external.set_line_prefix("")
    out = capsys.readouterr().out
    # each non-empty physical line tagged; the blank separator stays blank
    assert out == "[Gma] line1\n\n[Gma] line2\n"


def test_log_no_prefix_by_default(capsys):
    external.log("plain")
    assert capsys.readouterr().out == "plain\n"
