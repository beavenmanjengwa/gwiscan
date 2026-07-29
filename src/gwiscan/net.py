#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# net.py - Small resilient HTTP GET for the urllib-based downloads (Pfam HMMs, go-basic.obo).      #
#                                                                                                  #
# Both downloads fetch a single file from an EBI/OBO endpoint with the stdlib urllib. Two things   #
# they both need, learned the hard way: a User-Agent (purl.obolibrary.org 403s requests without    #
# one), and a few retries with backoff (a public endpoint has transient blips). fetch() centralises #
# both so neither call site repeats the pattern. (The InterProScan REST client uses requests +      #
# urllib3 Retry directly; this is only for the plain single-file GETs.)                            #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

import time
import urllib.error
import urllib.request

from . import __version__, external

USER_AGENT = f"gwiscan/{__version__} (https://github.com/; bioinformatics pipeline)"
# HTTP statuses worth retrying (transient): rate-limit + the 5xx family.
_RETRY_STATUSES = frozenset({429, 500, 502, 503, 504})
_DEFAULT_RETRIES = 5
_DEFAULT_BACKOFF = 2.0   # seconds: sleep = backoff * 2**attempt -> 2, 4, 8, 16


def fetch(url: str, *, timeout: float = 60, retries: int = _DEFAULT_RETRIES,
          backoff: float = _DEFAULT_BACKOFF, headers: dict | None = None) -> bytes:
    """GET ``url`` and return its bytes, sending a User-Agent and retrying transient
    failures (connection errors, timeouts, HTTP 429/5xx) with exponential backoff.

    A non-retryable HTTP error (e.g. 403/404) raises immediately -- retrying a
    consistent client/policy rejection only wastes time. After ``retries`` transient
    failures the last error is raised. ``retries`` is the number of RETRIES, so the
    URL is requested at most ``retries + 1`` times."""
    hdrs = {"User-Agent": USER_AGENT, **(headers or {})}
    last_exc: Exception | None = None

    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url, headers=hdrs)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except urllib.error.HTTPError as e:
            if e.code not in _RETRY_STATUSES:
                raise                       # 403/404/... -> not transient, don't retry
            last_exc = e
        except (urllib.error.URLError, TimeoutError) as e:
            last_exc = e                    # connection reset, DNS, timeout, ...

        if attempt < retries:
            wait = backoff * (2 ** attempt)
            external.log(f"[..] fetch failed ({last_exc}); retry {attempt + 1}/{retries} in {wait:.0f}s")
            time.sleep(wait)

    raise RuntimeError(f"failed to fetch {url} after {retries + 1} attempts: {last_exc}")
