#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# test_net.py - Tests for net.fetch, the resilient single-file HTTP GET.                           #
#                                                                                                  #
# Locks the two behaviours the download call sites depend on: a User-Agent is always sent          #
# (purl.obolibrary.org 403s requests without one), and transient failures (connection errors,      #
# 429/5xx) are retried with backoff while a genuine non-transient error (403/404) raises           #
# immediately.                                                                                     #
#                                                                                                  #
####################################################################################################
"""

import urllib.error

import pytest

from gwiscan import net


class _Resp:
    def __init__(self, data):
        self._data = data
    def read(self):
        return self._data
    def __enter__(self):
        return self
    def __exit__(self, *a):
        return False


def _patch_urlopen(monkeypatch, side_effects):
    """side_effects: list of either bytes (success) or an Exception to raise.
    Records the Request objects seen. Returns the record dict."""
    calls = {"requests": [], "n": 0}
    seq = list(side_effects)

    def _fake_urlopen(req, timeout=None):
        calls["requests"].append(req)
        calls["n"] += 1
        effect = seq.pop(0)
        if isinstance(effect, Exception):
            raise effect
        return _Resp(effect)

    monkeypatch.setattr(net.urllib.request, "urlopen", _fake_urlopen)
    monkeypatch.setattr(net.time, "sleep", lambda _s: None)  # no real waiting
    return calls


def test_fetch_sends_user_agent(monkeypatch):
    calls = _patch_urlopen(monkeypatch, [b"DATA"])
    assert net.fetch("https://example.org/x") == b"DATA"
    req = calls["requests"][0]
    assert req.get_header("User-agent") == net.USER_AGENT


def test_fetch_merges_extra_headers(monkeypatch):
    calls = _patch_urlopen(monkeypatch, [b"D"])
    net.fetch("https://example.org/x", headers={"Accept": "text/plain"})
    req = calls["requests"][0]
    assert req.get_header("User-agent") == net.USER_AGENT
    assert req.get_header("Accept") == "text/plain"


def test_fetch_retries_transient_5xx_then_succeeds(monkeypatch):
    err = urllib.error.HTTPError("u", 503, "busy", {}, None)
    calls = _patch_urlopen(monkeypatch, [err, err, b"OK"])
    assert net.fetch("https://example.org/x", retries=5) == b"OK"
    assert calls["n"] == 3   # two failures + one success


def test_fetch_retries_connection_error(monkeypatch):
    err = urllib.error.URLError("connection reset")
    calls = _patch_urlopen(monkeypatch, [err, b"OK"])
    assert net.fetch("https://example.org/x", retries=3) == b"OK"
    assert calls["n"] == 2


def test_fetch_does_not_retry_403(monkeypatch):
    err = urllib.error.HTTPError("u", 403, "forbidden", {}, None)
    calls = _patch_urlopen(monkeypatch, [err, b"OK"])
    with pytest.raises(urllib.error.HTTPError):
        net.fetch("https://example.org/x", retries=5)
    assert calls["n"] == 1   # raised immediately, no retry


def test_fetch_does_not_retry_404(monkeypatch):
    err = urllib.error.HTTPError("u", 404, "missing", {}, None)
    _patch_urlopen(monkeypatch, [err])
    with pytest.raises(urllib.error.HTTPError):
        net.fetch("https://example.org/x")


def test_fetch_exhausts_retries_and_raises_runtimeerror(monkeypatch):
    err = urllib.error.HTTPError("u", 500, "boom", {}, None)
    calls = _patch_urlopen(monkeypatch, [err, err, err])
    with pytest.raises(RuntimeError, match="after 3 attempts"):
        net.fetch("https://example.org/x", retries=2)   # 1 initial + 2 retries = 3
    assert calls["n"] == 3


# --- the download call sites route through net.fetch --------------------------

def test_setupdb_ensure_hmm_uses_net_fetch(monkeypatch, tmp_path):
    import gzip
    from gwiscan import setupdb
    from gwiscan.config import Config

    payload = gzip.compress(b"HMMER3/f model")
    seen = {}
    def _fake_fetch(url, **kw):
        seen["url"] = url
        return payload
    monkeypatch.setattr(setupdb.net, "fetch", _fake_fetch)

    cfg = Config(root=tmp_path)
    cfg.hmm_dir.mkdir(parents=True)
    dest = setupdb._ensure_hmm(cfg, "PF01453")
    assert dest.read_bytes() == b"HMMER3/f model"   # gunzipped
    assert "PF01453" in seen["url"]


def test_go_ensure_obo_uses_net_fetch(monkeypatch, tmp_path):
    from gwiscan import go
    from gwiscan.config import Config

    monkeypatch.setattr(go.net, "fetch", lambda url, **kw: b"format-version: 1.2\n")
    cfg = Config(root=tmp_path)
    dest = go.ensure_obo(cfg)
    assert dest.exists()
    assert dest.read_bytes().startswith(b"format-version")
