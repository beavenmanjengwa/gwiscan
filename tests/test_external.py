#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# test_external.py - Tests for external.py, the subprocess wrapper every stage funnels through.    #
#                                                                                                  #
# A regression here (wrong check field, swallowed stderr, bad stringification) would affect every  #
# external-tool stage at once, so the wrapper's own behaviour is verified independently of any one #
# stage.                                                                                           #
#                                                                                                  #
####################################################################################################
"""

import subprocess
import sys
from pathlib import Path

import pytest

from gwiscan import external

# Real subprocesses (via this interpreter) rather than a mocked subprocess: the
# wrapper now streams a live pipe, so exercising the real API is both simpler and
# a stronger guard than asserting Popen was called a certain way.
_PY = sys.executable


def _emit(out="", err="", code=0):
    """A portable command that writes given text to stdout/stderr and exits code."""
    script = (
        f"import sys\n"
        f"sys.stdout.write({out!r})\n"
        f"sys.stderr.write({err!r})\n"
        f"sys.exit({code})\n"
    )
    return [_PY, "-c", script]


def test_run_stringifies_command_elements(capsys):
    # A Path and an int in the command must be coerced to str (Popen rejects ints);
    # the child echoes them back so we can see they arrived intact.
    script = "import sys; print(sys.argv[1]); print(sys.argv[2])"
    result = external.run([_PY, "-c", script, Path("/a/b"), 42])
    assert "/a/b" in result.stdout
    assert "42" in result.stdout


def test_run_merges_stdout_stderr_when_no_stdout_path(capsys):
    result = external.run(_emit(out="hello-out\n", err="hello-err\n"))
    printed = capsys.readouterr().out
    # merged mode: both streams echoed live AND captured on result.stdout
    assert "hello-out" in printed and "hello-err" in printed
    assert "hello-out" in result.stdout and "hello-err" in result.stdout


def test_run_routes_stdout_to_file_and_echoes_only_stderr(capsys, tmp_path, monkeypatch):
    # Markers come from the environment so they appear only in the child's runtime
    # output, not in the echoed [CMD] line (which prints the script source).
    monkeypatch.setenv("GW_OUT", "TO_FILE")
    monkeypatch.setenv("GW_ERR", "ONLY_STDERR")
    script = ("import os, sys\n"
              "sys.stdout.write(os.environ['GW_OUT'] + '\\n')\n"
              "sys.stderr.write(os.environ['GW_ERR'] + '\\n')\n")
    out_file = tmp_path / "out.txt"
    result = external.run([_PY, "-c", script], stdout_path=out_file)
    # stdout went to the file, not the console; stderr streamed to the console
    assert "TO_FILE" in out_file.read_text()
    printed = capsys.readouterr().out
    assert "ONLY_STDERR" in printed
    assert "TO_FILE" not in printed
    assert result.stderr == "ONLY_STDERR\n"


def test_run_raises_on_nonzero_with_command_in_message(capsys):
    with pytest.raises(RuntimeError, match=r"exit 3"):
        external.run(_emit(code=3))


def test_run_check_false_does_not_raise(capsys):
    result = external.run(_emit(code=2), check=False)   # must not raise
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
    # Un-mocked calls to catch signature drift against the real subprocess API.
    assert external.run(_emit(code=0)).returncode == 0
    with pytest.raises(RuntimeError):
        external.run(_emit(code=1))
    # sanity: subprocess is really the stdlib module
    assert external.subprocess is subprocess


def test_run_streams_all_lines_in_order(capsys):
    # Every line of a multi-line output is echoed, in order (live streaming, not a
    # single dump), and fully captured on the result.
    script = "import sys\nfor i in range(5): print('line', i)\n"
    result = external.run([_PY, "-c", script])
    printed = capsys.readouterr().out
    for i in range(5):
        assert f"line {i}" in printed
    assert printed.index("line 0") < printed.index("line 4")
    assert result.stdout.count("line ") == 5


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
