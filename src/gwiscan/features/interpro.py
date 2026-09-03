#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# interpro.py - Domain / GO annotation via InterProScan (the `interpro` stage).                    #
#                                                                                                  #
# Two modes (INTERPRO_MODE), both producing the same outputs so downstream is identical:           #
#   * api   - EBI InterProScan REST API. IPRSCAN_VERSION picks the service: 5 (default, stable)    #
#             or 6. Chunks, polls, fetches TSV + GFF. Needs EBI_EMAIL and internet.                #
#   * local - a local install (INTERPROSCAN_BIN, e.g. interproscan.sh, system-wide on PATH).       #
#             Offline by default (adds -dp) unless INTERPRO_LOOKUP is set.                         #
#                                                                                                  #
# API v5 vs v6: the two services name the SAME member databases differently on the /run `appl`     #
# field (proven against the live API: v5 requires PfamA/Panther and rejects Pfam/PANTHER with      #
# HTTP 400; v6 requires Pfam/PANTHER). CDD is spelled the same in both. GWIscan derives its        #
# canonical names (Pfam/CDD/PANTHER) from the family table and maps them to the chosen version     #
# here, so the two must never be hand-mixed. v6 additionally carries a known server-side defect    #
# in its Matches-API lookup step ("Null key for a Map", modules/lookup/main.nf) that fails every   #
# job regardless of applications or batch size; a failed v6 chunk is diagnosed and the run stops   #
# with a clear instruction to rerun with --iprscan-version 5.                                      #
#                                                                                                  #
# Runs only on Pfam-accession-family candidates. Writes the combined interproscan.tsv (the file    #
# the pipeline reads; a short '#'-commented provenance block sits above its header), per-member-    #
# database splits (interproscan.pfam.tsv, interproscan.cdd.tsv), the merged interproscan.gff3,     #
# and interproscan.manifest.txt (service, endpoint, client, InterProScan version, applications,    #
# per-chunk job ids) so a reviewer can reconstruct which InterPro release produced the inventory.  #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

import re
import tempfile
import time
from datetime import datetime
from functools import lru_cache
from pathlib import Path

import requests
from Bio import SeqIO
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .. import __version__, external, io
from ..config import Config

# EBI InterProScan REST endpoints, by API version (IPRSCAN_VERSION).
IPRSCAN_ENDPOINTS = {
    5: "https://www.ebi.ac.uk/Tools/services/rest/iprscan5",
    6: "https://www.ebi.ac.uk/Tools/services/rest/iprscan6",
}

# GWIscan's canonical member-database names (from io.application_for: Pfam / CDD /
# PANTHER) -> the exact string each API version accepts on the /run `appl` field.
# Verified against the live service on 2026-09-03: v5 rejects "Pfam" and "PANTHER"
# with HTTP 400 ('Value for "appl" is not valid') and requires "PfamA" / "Panther";
# v6 requires "Pfam" / "PANTHER". CDD is identical in both. A canonical name not in
# a version's map is passed through unchanged, so a future database still submits.
_APPL_BY_VERSION = {
    5: {"Pfam": "PfamA", "PANTHER": "Panther", "CDD": "CDD"},
    6: {"Pfam": "Pfam", "PANTHER": "PANTHER", "CDD": "CDD"},
}

CHUNK_SIZE = 1000     # EBI REST API max sequences per submitted job
MAX_CONCURRENT_JOBS = 5   # EBI fair-use: max jobs kept in flight at once (cap is 30)
POLL_INTERVAL = 30    # seconds between status checks
MAX_WAIT = 3600       # max seconds to wait per job
# Network resilience: a multi-hour poll against a public REST service will hit the
# occasional transient failure (a dropped connection, a momentary 5xx, rate
# limiting). Retry those a few times with exponential backoff so one blip doesn't
# discard all the upstream HMM/DIAMOND work for that species. allowed_methods=None
# retries POST (the submit) too, not just idempotent GETs; raise_on_status=False
# lets the caller's own raise_for_status() produce the clean HTTP error message on
# a genuine, non-transient 4xx after retries are exhausted.
_RETRY_TOTAL = 5
_RETRY_BACKOFF = 2.0   # 0s, 2s, 4s, 8s, 16s between attempts
_RETRY_STATUSES = (429, 500, 502, 503, 504)

# Signatures of the known InterProScan 6 server-side lookup defect. When a v6 job
# FAILUREs, its error/log carries these; matching any one classifies the failure as
# the defect (not a data or pipeline problem) so the stop message can point at v5.
_V6_LOOKUP_SIGNATURES = (
    "Null key for a Map",
    "modules/lookup/main.nf",
    "INTERPROSCAN:LOOKUP:GET_MATCHES",
    "querying the Matches API",
)


@lru_cache(maxsize=1)
def _session() -> requests.Session:
    """A requests Session whose adapter retries transient network/HTTP failures
    with exponential backoff. Built once (cached) and reused for every call."""
    retry = Retry(
        total=_RETRY_TOTAL,
        connect=_RETRY_TOTAL,
        read=_RETRY_TOTAL,
        status=_RETRY_TOTAL,
        backoff_factor=_RETRY_BACKOFF,
        status_forcelist=_RETRY_STATUSES,
        allowed_methods=None,       # retry POST (submit) as well as GET
        raise_on_status=False,      # let callers raise_for_status() for a clean message
    )
    session = requests.Session()
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session
# GO-term lookup on; pathways off. Both are sent explicitly: the API defaults
# 'pathways' to true, so keeping it off requires setting it to "false".
GOTERMS = "true"
PATHWAYS = "false"

# InterProScan TSV columns, in order (the API result has NO header row). The
# trailing go_terms and pathways columns are OPTIONAL — present only when the
# --goterms / --pathways options are on (per the InterProScan output-format
# docs). So we name columns by how many the data actually has, rather than
# assuming a fixed count. compile reads protein_id / analysis / ipr_acc /
# ipr_desc / go_terms.
BASE_COLUMNS = ["protein_id", "seq_md5", "length", "analysis", "sig_acc",
                "sig_desc", "start", "end", "score", "status", "date",
                "ipr_acc", "ipr_desc"]
OPTIONAL_COLUMNS = ["go_terms", "pathways"]  # appended in this order when present


def _header_for(ncols: int) -> list:
    """Column names for a result row of ``ncols`` fields."""
    names = BASE_COLUMNS + OPTIONAL_COLUMNS
    if ncols <= len(names):
        return names[:ncols]
    # More columns than expected: pad generically so nothing is dropped.
    return names + [f"col_{i}" for i in range(len(names) + 1, ncols + 1)]


def _appl_list(appl) -> list:
    """Member-database ids as a clean list. The API takes one ``appl`` field per db,
    so we send a list (requests -> repeated appl= params)."""
    if isinstance(appl, list):
        return [a.strip() for a in appl if a.strip()]
    return [a.strip() for a in str(appl).split(",") if a.strip()]


def _run_applications(cfg: Config) -> list:
    """The InterProScan applications to run, in GWIscan's CANONICAL names (Pfam /
    CDD / PANTHER), enabled automatically from the family table: whichever member
    databases the confirmation accessions name, unioned across families.
    Architecture mode identifies on Pfam domains, so it uses Pfam. Falls back to
    Pfam if nothing is derived. Map to a version's submit names with
    _appl_for_version()."""
    if cfg.is_architecture:
        return ["Pfam"]
    return list(io.interpro_applications(cfg.family_map)) or ["Pfam"]


def _api_version(cfg: Config) -> int:
    """The EBI REST API version to submit to (IPRSCAN_VERSION): 5 or 6."""
    try:
        version = int(cfg.IPRSCAN_VERSION)
    except (TypeError, ValueError):
        version = -1
    if version not in IPRSCAN_ENDPOINTS:
        raise RuntimeError(
            f"unsupported IPRSCAN_VERSION {cfg.IPRSCAN_VERSION!r}; use 5 (default, "
            f"stable) or 6."
        )
    return version


def _appl_for_version(apps, version) -> list:
    """Map canonical application names to the strings the given API version accepts
    on /run (v5: PfamA/Panther, v6: Pfam/PANTHER; CDD both). Unknown names pass
    through unchanged."""
    mapping = _APPL_BY_VERSION[version]
    return [mapping.get(a, a) for a in apps]


def _submit(base, email, sequences_fasta, appl):
    """POST one chunk to <base>/run. `appl` is the already-version-mapped list."""
    resp = _session().post(f"{base}/run", data={
        "email": email,
        "title": "GWIscan_batch",
        "sequence": sequences_fasta,
        "appl": _appl_list(appl),   # one appl= per member db, version-mapped
        "goterms": GOTERMS,
        "pathways": PATHWAYS,
        "stype": "p",          # protein (explicit rather than relying on default)
    }, timeout=60)
    resp.raise_for_status()
    return resp.text.strip()


def _poll(base, job_id):
    elapsed = 0
    while elapsed < MAX_WAIT:
        resp = _session().get(f"{base}/status/{job_id}", timeout=30)
        resp.raise_for_status()
        status = resp.text.strip()
        # EBI's terminal failure status is FAILURE (not FAILED); include both so a
        # failed job returns immediately instead of polling until MAX_WAIT.
        if status in ("FINISHED", "FAILURE", "FAILED", "ERROR", "NOT_FOUND"):
            return status
        time.sleep(POLL_INTERVAL)
        elapsed += POLL_INTERVAL
        external.log(f"  [{job_id}] status: {status} ({elapsed}s elapsed)")
    return "TIMEOUT"


def _result_types(base, job_id):
    """Available result-type identifiers for a finished job (self-documenting,
    like /parameters) — so we don't guess the GFF type's id."""
    resp = _session().get(f"{base}/resulttypes/{job_id}", timeout=30)
    resp.raise_for_status()
    return re.findall(r"<identifier>([^<]+)</identifier>", resp.text)


def _fetch_result(base, job_id, rtype):
    resp = _session().get(f"{base}/result/{job_id}/{rtype}", timeout=120)
    resp.raise_for_status()
    return resp.text


_VERSION_RE = re.compile(r'interproscan-version["\s:=]+([0-9][0-9.\-]+)')


def _version_from_text(text) -> str:
    """The 'interproscan-version' string (e.g. '5.78-109.0' = software-release) if
    present in an InterProScan XML/JSON/GFF result, else ''."""
    m = _VERSION_RE.search(text or "")
    return m.group(1) if m else ""


def _job_version(base, job_id, gff_text="") -> str:
    """Best-effort authoritative version for a FINISHED job: the
    'interproscan-version' attribute (software + InterPro release, e.g.
    '5.78-109.0'). Tries the GFF already in hand first (no extra request), then a
    small streamed read of the XML head (the attribute sits in the root element),
    then JSON. '' if none is retrievable."""
    found = _version_from_text(gff_text)
    if found:
        return found
    # Stream the XML head only — the attribute is in the opening root tag.
    try:
        with _session().get(f"{base}/result/{job_id}/xml", timeout=60, stream=True) as r:
            r.raise_for_status()
            head = ""
            for chunk in r.iter_content(chunk_size=1024):
                head += chunk.decode("utf-8", "ignore") if isinstance(chunk, bytes) else chunk
                found = _version_from_text(head)
                if found or len(head) > 8192:
                    break
            if found:
                return found
    except Exception:  # noqa: BLE001
        pass
    for rtype in ("json", "out"):
        try:
            found = _version_from_text(_fetch_result(base, job_id, rtype))
        except Exception:  # noqa: BLE001
            continue
        if found:
            return found
    return ""


def _failure_diagnostics(base, job_id) -> str:
    """Best-effort error/log text for a failed job (for diagnosis). '' if none."""
    for rtype in ("error", "log", "out"):
        try:
            txt = _fetch_result(base, job_id, rtype)
        except Exception:  # noqa: BLE001
            continue
        if txt and txt.strip():
            return txt
    return ""


def _is_v6_lookup_defect(diag) -> bool:
    return any(sig in diag for sig in _V6_LOOKUP_SIGNATURES)


def _fail(cfg, version, base, i, n_chunks, job_id, status, n_seqs, diag) -> RuntimeError:
    """The error a non-FINISHED chunk raises. On the known InterProScan 6 lookup
    defect, name it plainly and tell the user to rerun with --iprscan-version 5 —
    the tool deciding to stop rather than silently dropping the chunk's members."""
    if version == 6 and _is_v6_lookup_defect(diag):
        return RuntimeError(
            f"InterProScan 6 chunk {i + 1}/{n_chunks} (job {job_id}) failed: {status}. "
            f"Known InterProScan 6 server-side lookup failure. "
            f"Rerun with --iprscan-version 5."
        )
    switch = " Try --iprscan-version 5." if version == 6 else ""
    return RuntimeError(
        f"InterProScan {version} chunk {i + 1}/{n_chunks} (job {job_id}) did not finish: "
        f"{status}.{switch} Re-run once EBI is reachable "
        f"(gwiscan run --from-stage interpro ...), or use INTERPRO_MODE=local."
    )


def _combine_gff(gff_texts) -> str:
    """Merge per-chunk GFF3 into one valid feature GFF3: a single version pragma,
    all feature/pragma lines, and the trailing ##FASTA sequence blocks dropped
    (the sequences are already in candidates.fasta)."""
    out = ["##gff-version 3"]
    for text in gff_texts:
        for line in text.splitlines():
            if line.startswith("##gff-version"):
                continue
            if line.strip() == "##FASTA":
                break
            out.append(line)
    return "\n".join(out) + "\n"


def _run_api(cfg: Config, cand_fasta) -> tuple:
    """EBI InterProScan REST API (IPRSCAN_VERSION 5 or 6): chunk, submit, poll,
    fetch TSV + GFF. Returns (tsv_lines, gff_texts, meta)."""
    if not cfg.EBI_EMAIL:
        raise RuntimeError(
            "InterProScan API mode needs an email: set EBI_EMAIL, or use "
            "INTERPRO_MODE=local with a local install."
        )

    version = _api_version(cfg)
    base = IPRSCAN_ENDPOINTS[version]
    apps = _run_applications(cfg)
    submit_apps = _appl_for_version(apps, version)

    records = list(SeqIO.parse(str(cand_fasta), "fasta"))
    external.log(f"[OK] Loaded {len(records)} candidate sequences")
    n_chunks = (len(records) + CHUNK_SIZE - 1) // CHUNK_SIZE
    external.log(f"[OK] Submitting {n_chunks} job(s) to EBI InterProScan {version} API "
                 f"(applications: {', '.join(submit_apps)})...")

    meta = {
        "service": f"InterProScan {version} (EBI REST API)",
        "api_version": version,
        "endpoint": base,
        "client": f"gwiscan {__version__} (requests REST client)",
        "applications_canonical": list(apps),
        "applications_submitted": list(submit_apps),
        "date": datetime.now().isoformat(timespec="seconds"),
        "chunks": [],
        "interproscan_versions": [],
    }

    # Submit in waves of at most MAX_CONCURRENT_JOBS jobs: submit the wave, wait for
    # and fetch every job in it, then start the next wave. EMBL-EBI's Job Dispatcher
    # asks that you keep concurrent jobs modest (its cap is 30) and not submit more
    # until the running ones finish; capping a species at 5 in flight keeps several
    # species running in parallel (SPECIES_PARALLEL) within that fair-use limit while
    # still overlapping jobs for speed. Do NOT raise the cap or submit every chunk up
    # front -- that stacks a species' whole chunk count as concurrent jobs and,
    # multiplied across parallel species, breaches the limit.
    all_lines, gff_texts = [], []
    for wave_start in range(0, n_chunks, MAX_CONCURRENT_JOBS):
        wave = list(range(wave_start, min(wave_start + MAX_CONCURRENT_JOBS, n_chunks)))

        # Submit the whole wave up front, so its jobs run on EBI concurrently.
        jobs = []
        for i in wave:
            chunk = records[i * CHUNK_SIZE:(i + 1) * CHUNK_SIZE]
            chunk_fasta = "\n".join(f">{rec.id}\n{rec.seq}" for rec in chunk)
            external.log(f"\n[Chunk {i + 1}/{n_chunks}] Submitting {len(chunk)} sequences...")
            job_id = _submit(base, cfg.EBI_EMAIL, chunk_fasta, submit_apps)
            external.log(f"[Chunk {i + 1}/{n_chunks}] Job ID: {job_id}")
            jobs.append((i, job_id, chunk))

        # Wait for and fetch each job in the wave before starting the next wave.
        for i, job_id, chunk in jobs:
            status = _poll(base, job_id)
            external.log(f"[Chunk {i + 1}/{n_chunks}] Final status: {status}")
            # A chunk that does not FINISH leaves its sequences unannotated. The confirm
            # stage keeps a Pfam family's members only when their Pfam is reported here,
            # so a missing annotation is indistinguishable from a genuine non-match: the
            # members would be dropped and the run would still exit 0 with a plausible
            # but incomplete table. Stop with a clear, resumable error instead. (The
            # session already retries transient blips, so reaching here means a real
            # failure, not a momentary one.)
            if status != "FINISHED":
                diag = _failure_diagnostics(base, job_id)
                raise _fail(cfg, version, base, i, n_chunks, job_id, status, len(chunk), diag)

            # Discover result-type ids for this job, then fetch TSV (pipeline table)
            # and GFF (raw feature output) using the ids the API actually offers.
            types = _result_types(base, job_id)
            tsv_id = "tsv" if "tsv" in types else None
            gff_id = next((t for t in types if t.lower() in ("gff", "gff3")), None)

            # No TSV on a FINISHED job is as damaging as an unfinished chunk (same
            # silent-drop), so it is likewise fatal. A missing GFF only costs the
            # auxiliary feature file, which downstream stages do not read, so it stays
            # a warning.
            if not tsv_id:
                raise RuntimeError(
                    f"InterProScan chunk {i + 1}/{n_chunks} (job {job_id}) finished but offers "
                    f"no 'tsv' result type ({types}); its {len(chunk)} sequence(s) cannot be "
                    f"annotated. Re-run the interpro stage, or use INTERPRO_MODE=local."
                )
            lines = [ln for ln in _fetch_result(base, job_id, tsv_id).strip().split("\n") if ln]
            all_lines.extend(lines)
            external.log(f"[Chunk {i + 1}/{n_chunks}] Retrieved {len(lines)} TSV rows.")

            gff_text = ""
            if gff_id:
                gff_text = _fetch_result(base, job_id, gff_id)
                gff_texts.append(gff_text)
            else:
                external.log(f"[WARN] Job {job_id}: no GFF result type ({types}).")

            # Authoritative per-chunk version (software + InterPro release) — recorded
            # for the manifest and used to enforce one version across the whole run.
            ipr_version = _job_version(base, job_id, gff_text)
            meta["chunks"].append({
                "chunk": f"{i + 1}/{n_chunks}", "job_id": job_id,
                "n_sequences": len(chunk), "interproscan_version": ipr_version or "unknown",
                "tsv_rows": len(lines),
            })
            if ipr_version and ipr_version not in meta["interproscan_versions"]:
                meta["interproscan_versions"].append(ipr_version)

    _enforce_single_version(meta)
    return all_lines, gff_texts, meta


def _enforce_single_version(meta) -> None:
    """Hard check (not a warning): refuse to combine annotations produced under more
    than one InterProScan/InterPro version in a single run. Mixed releases make the
    inventory unreconstructable, so the run stops before any combined output is
    written."""
    seen = list(meta.get("interproscan_versions", []))
    if len(seen) > 1:
        rows = "; ".join(f"{c['chunk']} (job {c['job_id']}): {c['interproscan_version']}"
                         for c in meta["chunks"])
        raise RuntimeError(
            "InterProScan results in this run came from MORE THAN ONE version: "
            f"{sorted(seen)}. Combining annotations from different InterPro releases "
            "would make the inventory impossible to reconstruct, so the combined "
            "output is NOT written. This usually means EBI rolled a release mid-run. "
            "Re-run the interpro stage so every chunk is annotated by one release.\n"
            f"  Per chunk: {rows}"
        )


def _local_cmd(cfg: Config, cand_fasta, outdir) -> list:
    """interproscan.sh command line for local mode. Writes TSV + GFF3 into
    outdir; adds -dp (disable the online precalc lookup) unless INTERPRO_LOOKUP."""
    cmd = [
        cfg.INTERPROSCAN_BIN,
        "-i", str(cand_fasta),
        "-f", "TSV,GFF3",
        "-d", str(outdir),
        "-appl", ",".join(_run_applications(cfg)),
        "-goterms",
        "-cpu", str(cfg.THREADS),
    ]
    if not cfg.INTERPRO_LOOKUP:
        cmd.append("-dp")     # offline: no precalculated-match lookup service
    return cmd


def _run_local(cfg: Config, cand_fasta) -> tuple:
    """Local InterProScan install: run interproscan.sh, read its TSV + GFF3.
    Returns (tsv_lines, gff_texts, meta). The local TSV has the same columns as the
    API result, so downstream parsing is identical."""
    external.require(cfg.INTERPROSCAN_BIN)
    apps = _run_applications(cfg)
    with tempfile.TemporaryDirectory(prefix="gwiscan_ipr_") as tmp:
        external.log(f"[OK] Running local InterProScan: {cfg.INTERPROSCAN_BIN}")
        external.run(_local_cmd(cfg, cand_fasta, tmp))
        tmp_p = Path(tmp)
        all_lines = []
        for f in sorted(tmp_p.glob("*.tsv")):
            all_lines.extend(ln for ln in f.read_text().splitlines() if ln.strip())
        gff_texts = [f.read_text() for f in sorted(tmp_p.glob("*.gff3"))]
    external.log(f"[OK] Local InterProScan: {len(all_lines)} TSV rows")
    ipr_version = ""
    for text in gff_texts:
        ipr_version = _version_from_text(text)
        if ipr_version:
            break
    meta = {
        "service": "InterProScan (local install)",
        "api_version": "local",
        "endpoint": str(cfg.INTERPROSCAN_BIN),
        "client": f"gwiscan {__version__} (local subprocess)",
        "applications_canonical": list(apps),
        "applications_submitted": list(apps),
        "date": datetime.now().isoformat(timespec="seconds"),
        "chunks": [{"chunk": "1/1", "job_id": "local", "n_sequences": None,
                    "interproscan_version": ipr_version or "unknown",
                    "tsv_rows": len(all_lines)}],
        "interproscan_versions": [ipr_version] if ipr_version else [],
    }
    return all_lines, gff_texts, meta


_ANALYSIS_COL = BASE_COLUMNS.index("analysis")   # member-database column (Pfam, CDD, ...)


def _provenance_preamble(meta) -> list:
    """The '#'-commented provenance lines written above interproscan.tsv's header,
    so a reviewer opening the TSV sees which service and InterProScan/InterPro
    release produced it. io.read_interpro_tsv() strips these before parsing."""
    ipr_versions = meta.get("interproscan_versions") or ["unknown"]
    return [
        f"# GWIscan interpro stage — {meta['service']}",
        f"# endpoint: {meta['endpoint']}",
        f"# client: {meta['client']}",
        f"# applications: {', '.join(meta['applications_submitted'])}",
        f"# interproscan-version (software-InterPro release): {', '.join(ipr_versions)}",
        f"# generated: {meta['date']}",
    ]


def _write_tsv(path, header, lines, preamble=None) -> None:
    with open(path, "w") as out:
        for pre in (preamble or []):
            out.write(pre + "\n")
        out.write(header + "\n")
        for line in lines:
            out.write(line + "\n")


def _write_outputs(all_lines, gff_texts, out_tsv, out_gff, preamble=None) -> None:
    """Write the combined interproscan.tsv (the file the pipeline reads) and the
    merged interproscan.gff3, plus one per-member-database split file
    (interproscan_Pfam.tsv, interproscan_CDD.tsv, ...) for convenience. Each TSV
    carries the same '#'-commented provenance preamble above its header. Same for
    API and local mode."""
    ncols = max((len(line.split("\t")) for line in all_lines), default=len(BASE_COLUMNS))
    header = "\t".join(io.to_camel(c) for c in _header_for(ncols))

    # Combined table — the canonical file every downstream stage reads.
    _write_tsv(out_tsv, header, all_lines, preamble)

    # Per-member-database split (by the 'analysis' column) alongside the combined,
    # named interproscan.<db>.tsv (e.g. interproscan.pfam.tsv, interproscan.cdd.tsv).
    by_analysis = {}
    for line in all_lines:
        cols = line.split("\t")
        analysis = cols[_ANALYSIS_COL] if len(cols) > _ANALYSIS_COL else "unknown"
        by_analysis.setdefault(analysis, []).append(line)
    for analysis, lines in by_analysis.items():
        db = re.sub(r"[^A-Za-z0-9]+", "", analysis).lower() or "unknown"
        _write_tsv(out_tsv.with_name(f"{out_tsv.stem}.{db}{out_tsv.suffix}"), header, lines, preamble)

    # Always write the GFF (empty -> just the version pragma) so it is a stable,
    # trackable output alongside the TSV.
    out_gff.write_text(_combine_gff(gff_texts))


def _write_manifest(path, meta) -> None:
    """The run manifest: service, endpoint, client, applications, InterProScan
    version(s), and every chunk's job id — the reviewer-facing record of which
    InterPro release produced the inventory."""
    lines = [
        "GWIscan InterProScan run manifest",
        f"service: {meta['service']}",
        f"api_version: {meta['api_version']}",
        f"endpoint: {meta['endpoint']}",
        f"client: {meta['client']}",
        f"applications (canonical): {', '.join(meta['applications_canonical'])}",
        f"applications (submitted): {', '.join(meta['applications_submitted'])}",
        f"interproscan-version: {', '.join(meta.get('interproscan_versions') or ['unknown'])}",
        f"generated: {meta['date']}",
        "",
        "chunks:",
    ]
    for c in meta["chunks"]:
        n = c["n_sequences"]
        n_str = f"{n} seqs, " if n is not None else ""
        lines.append(f"  chunk {c['chunk']}: job {c['job_id']} "
                     f"({n_str}{c['tsv_rows']} rows, version {c['interproscan_version']})")
    Path(path).write_text("\n".join(lines) + "\n")


def pfam_candidate_ids(cfg: Config) -> set:
    """Protein ids that are candidates of a family identified by a real Pfam
    accession. Only these are sent to InterProScan; custom-HMM families (hmmscan +
    BLAST) and BLAST-only families are confirmed by their own search and never
    touch InterProScan."""
    merged = cfg.result("candidates_merged.tsv")
    if not merged.exists():
        return set()
    df = io.read_tsv(merged)
    # Architecture mode: every final candidate is annotated (no per-family Pfam gate).
    if cfg.is_architecture:
        return set(df["protein_id"].astype(str))
    pfam_families = {r["family"] for r in io.family_records(cfg.family_map) if r["pfam_model"]}
    return set(df.loc[df["family"].isin(pfam_families), "protein_id"].astype(str))


def run(cfg: Config) -> None:
    cfg.ensure_dirs()
    cand_fasta = cfg.result("candidates.fasta")
    out_tsv = cfg.result("interproscan.tsv")
    out_gff = cfg.result("interproscan.gff3")
    out_manifest = cfg.result("interproscan.manifest.txt")

    if not cand_fasta.exists():
        raise FileNotFoundError(f"candidates.fasta not found: {cand_fasta}")

    # Submit only Pfam-accession families; custom-HMM / BLAST-only families skip
    # InterProScan entirely (they are confirmed by their own hmmscan/BLAST).
    keep_ids = pfam_candidate_ids(cfg)
    all_records = list(SeqIO.parse(str(cand_fasta), "fasta"))
    records = [r for r in all_records if r.id in keep_ids]
    if not records:
        external.log("[OK] No Pfam-accession candidates; skipping InterProScan.")
        _write_outputs([], [], out_tsv, out_gff)
        return

    ipr_fasta = cfg.result("interpro_input.fasta")
    SeqIO.write(records, str(ipr_fasta), "fasta")
    # Only mention the skipped ones when there actually are any.
    n_skipped = len(all_records) - len(records)
    skip_note = f" ({n_skipped} custom-HMM / BLAST-only candidates skip it)" if n_skipped else ""
    external.log(f"[OK] {len(records)} Pfam-family candidates sent to InterProScan{skip_note}")

    mode = str(cfg.INTERPRO_MODE).lower()
    if mode == "local":
        all_lines, gff_texts, meta = _run_local(cfg, ipr_fasta)
    elif mode == "api":
        all_lines, gff_texts, meta = _run_api(cfg, ipr_fasta)
    else:
        raise RuntimeError(
            f"unknown INTERPRO_MODE {cfg.INTERPRO_MODE!r} (use 'api' or 'local')"
        )

    _write_outputs(all_lines, gff_texts, out_tsv, out_gff, _provenance_preamble(meta))
    _write_manifest(out_manifest, meta)
    external.log(
        f"\n[OK] InterProScan complete ({mode}, {meta['service']}): {len(all_lines)} "
        f"annotations -> interproscan.tsv + interproscan.gff3 + interproscan.manifest.txt"
    )
