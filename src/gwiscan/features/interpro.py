#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# interpro.py - Domain / GO annotation via InterProScan 6 (the `interpro` stage).                  #
#                                                                                                  #
# Two modes (INTERPRO_MODE), both producing the same outputs so downstream is identical:           #
#   * api   - EBI InterProScan 6 REST API (chunks of 30, polls, fetches TSV + GFF3). Needs         #
#             EBI_EMAIL and internet.                                                              #
#   * local - a local install (INTERPROSCAN_BIN, e.g. interproscan.sh, system-wide on PATH).       #
#             Offline by default (adds -dp) unless INTERPRO_LOOKUP is set.                         #
#                                                                                                  #
# Runs only on Pfam-accession-family candidates. Writes the combined interproscan.tsv (the file    #
# the pipeline reads), per-member-database splits (interproscan.pfam.tsv, interproscan.cdd.tsv),   #
# and the merged interproscan.gff3.                                                                #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

import re
import tempfile
import time
from functools import lru_cache
from pathlib import Path

import requests
from Bio import SeqIO
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .. import external, io
from ..config import Config

BASE_URL = "https://www.ebi.ac.uk/Tools/services/rest/iprscan6"
CHUNK_SIZE = 30       # EBI limit per submission
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
    """Member-database ids as a clean list. INTERPRO_APPL is comma-separated
    (default "Pfam"; opt into CDD with "Pfam,CDD"); the API takes one ``appl``
    field per db, so we send a list (requests -> repeated appl= params)."""
    if isinstance(appl, list):
        return [a.strip() for a in appl if a.strip()]
    return [a.strip() for a in str(appl).split(",") if a.strip()]


def _submit(email, sequences_fasta, appl):
    resp = _session().post(f"{BASE_URL}/run", data={
        "email": email,
        "title": "GWIscan_batch",
        "sequence": sequences_fasta,
        "appl": _appl_list(appl),   # one appl= per member db; default ["Pfam"]
        "goterms": GOTERMS,
        "pathways": PATHWAYS,
        "stype": "p",          # protein (explicit rather than relying on default)
    }, timeout=60)
    resp.raise_for_status()
    return resp.text.strip()


def _poll(job_id):
    elapsed = 0
    while elapsed < MAX_WAIT:
        resp = _session().get(f"{BASE_URL}/status/{job_id}", timeout=30)
        resp.raise_for_status()
        status = resp.text.strip()
        if status in ("FINISHED", "FAILED", "ERROR", "NOT_FOUND"):
            return status
        time.sleep(POLL_INTERVAL)
        elapsed += POLL_INTERVAL
        external.log(f"  [{job_id}] status: {status} ({elapsed}s elapsed)")
    return "TIMEOUT"


def _result_types(job_id):
    """Available result-type identifiers for a finished job (self-documenting,
    like /parameters) — so we don't guess the GFF type's id."""
    resp = _session().get(f"{BASE_URL}/resulttypes/{job_id}", timeout=30)
    resp.raise_for_status()
    return re.findall(r"<identifier>([^<]+)</identifier>", resp.text)


def _fetch_result(job_id, rtype):
    resp = _session().get(f"{BASE_URL}/result/{job_id}/{rtype}", timeout=120)
    resp.raise_for_status()
    return resp.text


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
    """EBI InterProScan 6 REST API: chunk, submit, poll, fetch TSV + GFF3.
    Returns (tsv_lines, gff_texts)."""
    if not cfg.EBI_EMAIL:
        raise RuntimeError(
            "InterProScan API mode needs an email: set EBI_EMAIL, or use "
            "INTERPRO_MODE=local with a local install."
        )

    records = list(SeqIO.parse(str(cand_fasta), "fasta"))
    external.log(f"[OK] Loaded {len(records)} candidate sequences")
    n_chunks = (len(records) + CHUNK_SIZE - 1) // CHUNK_SIZE
    external.log(f"[OK] Submitting {n_chunks} job(s) to EBI InterProScan 6 API...")

    # One job in flight at a time: submit a chunk, wait for it, fetch it, then move
    # to the next. EMBL-EBI's Job Dispatcher asks that no more than a handful of jobs
    # run concurrently and that you not submit more until the running ones finish;
    # keeping a single job open per species means several species running in parallel
    # (SPECIES_PARALLEL) still stay well within that fair-use limit. Do NOT change this
    # to submit every chunk up front -- that stacks a species' whole chunk count as
    # concurrent jobs and, multiplied across parallel species, breaches the limit.
    all_lines, gff_texts = [], []
    for i in range(n_chunks):
        chunk = records[i * CHUNK_SIZE:(i + 1) * CHUNK_SIZE]
        chunk_fasta = "\n".join(f">{rec.id}\n{rec.seq}" for rec in chunk)

        external.log(f"\n[Chunk {i + 1}/{n_chunks}] Submitting {len(chunk)} sequences...")
        job_id = _submit(cfg.EBI_EMAIL, chunk_fasta, cfg.INTERPRO_APPL)
        external.log(f"[Chunk {i + 1}/{n_chunks}] Job ID: {job_id}")

        status = _poll(job_id)
        external.log(f"[Chunk {i + 1}/{n_chunks}] Final status: {status}")
        # A chunk that does not FINISH leaves its sequences unannotated. The confirm
        # stage keeps a Pfam family's members only when their Pfam is reported here,
        # so a missing annotation is indistinguishable from a genuine non-match: the
        # members would be dropped and the run would still exit 0 with a plausible
        # but incomplete table. Stop with a clear, resumable error instead. (The
        # session already retries transient blips, so reaching here means a real
        # failure, not a momentary one.)
        if status != "FINISHED":
            raise RuntimeError(
                f"InterProScan chunk {i + 1}/{n_chunks} (job {job_id}) did not finish: "
                f"status {status}. Its {len(chunk)} sequence(s) would be left unannotated. "
                f"Rather than silently drop those members, the run stops here. Re-run this "
                f"stage once EBI is reachable (gwiscan run --from-stage interpro ...), or "
                f"switch to a local install (INTERPRO_MODE=local)."
            )

        # Discover result-type ids for this job, then fetch TSV (pipeline table)
        # and GFF3 (raw feature output) using the ids the API actually offers.
        types = _result_types(job_id)
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
        lines = [ln for ln in _fetch_result(job_id, tsv_id).strip().split("\n") if ln]
        all_lines.extend(lines)
        external.log(f"[Chunk {i + 1}/{n_chunks}] Retrieved {len(lines)} TSV rows.")

        if gff_id:
            gff_texts.append(_fetch_result(job_id, gff_id))
        else:
            external.log(f"[WARN] Job {job_id}: no GFF result type ({types}).")

        if i < n_chunks - 1:
            time.sleep(5)  # be polite to EBI between jobs

    return all_lines, gff_texts


def _local_cmd(cfg: Config, cand_fasta, outdir) -> list:
    """interproscan.sh command line for local mode. Writes TSV + GFF3 into
    outdir; adds -dp (disable the online precalc lookup) unless INTERPRO_LOOKUP."""
    cmd = [
        cfg.INTERPROSCAN_BIN,
        "-i", str(cand_fasta),
        "-f", "TSV,GFF3",
        "-d", str(outdir),
        "-appl", ",".join(_appl_list(cfg.INTERPRO_APPL)),
        "-goterms",
        "-cpu", str(cfg.THREADS),
    ]
    if not cfg.INTERPRO_LOOKUP:
        cmd.append("-dp")     # offline: no precalculated-match lookup service
    return cmd


def _run_local(cfg: Config, cand_fasta) -> tuple:
    """Local InterProScan install: run interproscan.sh, read its TSV + GFF3.
    Returns (tsv_lines, gff_texts). The local TSV has the same columns as the
    API result, so downstream parsing is identical."""
    external.require(cfg.INTERPROSCAN_BIN)
    with tempfile.TemporaryDirectory(prefix="gwiscan_ipr_") as tmp:
        external.log(f"[OK] Running local InterProScan: {cfg.INTERPROSCAN_BIN}")
        external.run(_local_cmd(cfg, cand_fasta, tmp))
        tmp_p = Path(tmp)
        all_lines = []
        for f in sorted(tmp_p.glob("*.tsv")):
            all_lines.extend(ln for ln in f.read_text().splitlines() if ln.strip())
        gff_texts = [f.read_text() for f in sorted(tmp_p.glob("*.gff3"))]
    external.log(f"[OK] Local InterProScan: {len(all_lines)} TSV rows")
    return all_lines, gff_texts


_ANALYSIS_COL = BASE_COLUMNS.index("analysis")   # member-database column (Pfam, CDD, ...)


def _write_tsv(path, header, lines) -> None:
    with open(path, "w") as out:
        out.write(header + "\n")
        for line in lines:
            out.write(line + "\n")


def _write_outputs(all_lines, gff_texts, out_tsv, out_gff) -> None:
    """Write the combined interproscan.tsv (the file the pipeline reads) and the
    merged interproscan.gff3, plus one per-member-database split file
    (interproscan_Pfam.tsv, interproscan_CDD.tsv, ...) for convenience. Same for
    API and local mode."""
    ncols = max((len(line.split("\t")) for line in all_lines), default=len(BASE_COLUMNS))
    header = "\t".join(io.to_camel(c) for c in _header_for(ncols))

    # Combined table — the canonical file every downstream stage reads.
    _write_tsv(out_tsv, header, all_lines)

    # Per-member-database split (by the 'analysis' column) alongside the combined,
    # named interproscan.<db>.tsv (e.g. interproscan.pfam.tsv, interproscan.cdd.tsv).
    by_analysis = {}
    for line in all_lines:
        cols = line.split("\t")
        analysis = cols[_ANALYSIS_COL] if len(cols) > _ANALYSIS_COL else "unknown"
        by_analysis.setdefault(analysis, []).append(line)
    for analysis, lines in by_analysis.items():
        db = re.sub(r"[^A-Za-z0-9]+", "", analysis).lower() or "unknown"
        _write_tsv(out_tsv.with_name(f"{out_tsv.stem}.{db}{out_tsv.suffix}"), header, lines)

    # Always write the GFF (empty -> just the version pragma) so it is a stable,
    # trackable output alongside the TSV.
    out_gff.write_text(_combine_gff(gff_texts))


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
        all_lines, gff_texts = _run_local(cfg, ipr_fasta)
    elif mode == "api":
        all_lines, gff_texts = _run_api(cfg, ipr_fasta)
    else:
        raise RuntimeError(
            f"unknown INTERPRO_MODE {cfg.INTERPRO_MODE!r} (use 'api' or 'local')"
        )

    _write_outputs(all_lines, gff_texts, out_tsv, out_gff)
    external.log(
        f"\n[OK] InterProScan complete ({mode}): {len(all_lines)} annotations "
        f"-> interproscan.tsv + interproscan.gff3"
    )
