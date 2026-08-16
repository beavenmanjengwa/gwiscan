#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# external.py - Helpers for invoking external bioinformatics binaries.                             #
#                                                                                                  #
# Every stage that shells out to a tool (hmmscan, hmmpress, diamond, mafft, targetp, biolib,       #
# deeploc2, weblogo, iqtree) goes through run(), so command logging and failure handling are       #
# uniform across stages. Stages print to stdout/stderr; the CLI/pipeline layer captures each       #
# stage's output into logs/<stage>.log.                                                            #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path


# Optional per-process line prefix. In multi-species mode each species runs in
# its own process (ProcessPoolExecutor) but they share the terminal, so their
# lines interleave; the species worker sets this to "[<prefix>] " so every line
# is self-identifying. Empty (the default) => no prefix, single-species output
# is unchanged. Module-global is safe: it's per-process, never shared/raced.
_LINE_PREFIX = ""


def set_line_prefix(prefix: str) -> None:
    """Prefix every subsequent log() line (this process only). '' clears it."""
    global _LINE_PREFIX
    _LINE_PREFIX = prefix


def log(msg: str) -> None:
    if _LINE_PREFIX:
        # Prefix every physical line so multi-line messages stay tagged and
        # blank separator lines aren't turned into stray "[Ath] " noise.
        msg = "\n".join(_LINE_PREFIX + ln if ln else ln for ln in msg.split("\n"))
    print(msg, flush=True)


def available(binary: str) -> bool:
    """True if a tool is usable: resolvable on PATH by name, or given as an
    existing file path (an absolute/relative path to a user-installed binary)."""
    return shutil.which(binary) is not None or Path(binary).exists()


def require(binary: str) -> None:
    """Raise if a required external tool is neither on PATH nor an existing path."""
    if not available(binary):
        raise FileNotFoundError(f"required tool not found on PATH: {binary}")


def _stream(pipe) -> str:
    """Echo a text pipe to our stdout line by line as it arrives, and return the
    full captured text. Streaming (rather than reading once at the end) gives live
    progress for the long-running tools -- hmmscan, IQ-TREE, local InterProScan --
    instead of a silent wait until the tool exits."""
    captured = []
    for line in pipe:
        sys.stdout.write(line)
        sys.stdout.flush()
        captured.append(line)
    return "".join(captured)


def run(cmd, check: bool = True, stdout_path=None, cwd=None) -> subprocess.CompletedProcess:
    """Run an external command, echoing it and its output live.

    cmd elements are str()-cast, so Paths and ints are fine. If stdout_path is
    given, the child's stdout is written there and only stderr is echoed;
    otherwise stdout and stderr are merged and echoed. Output is streamed as it is
    produced. The returned CompletedProcess carries the captured echoed stream
    (its stdout in merged mode, its stderr in stdout_path mode). Raises
    RuntimeError on non-zero exit when check.
    """
    cmd = [str(c) for c in cmd]
    cwd = str(cwd) if cwd is not None else None
    log(f"[CMD] {' '.join(cmd)}")

    if stdout_path is not None:
        # Child stdout goes straight to the file; only stderr is streamed to us.
        with open(stdout_path, "w") as out:
            proc = subprocess.Popen(cmd, stdout=out, stderr=subprocess.PIPE,
                                    text=True, cwd=cwd)
            stderr_text = _stream(proc.stderr)
            proc.wait()
        result = subprocess.CompletedProcess(cmd, proc.returncode,
                                             stdout=None, stderr=stderr_text)
    else:
        # Merge stderr into stdout and stream the single combined pipe.
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True, cwd=cwd)
        stdout_text = _stream(proc.stdout)
        proc.wait()
        result = subprocess.CompletedProcess(cmd, proc.returncode,
                                             stdout=stdout_text, stderr=None)

    if check and result.returncode != 0:
        raise RuntimeError(f"command failed (exit {result.returncode}): {' '.join(cmd)}")
    return result
