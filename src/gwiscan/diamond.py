#!/usr/bin/env python3
"""
#######################################################################################################
#                                                                                                     #
# diamond.py - Two-round DIAMOND BLASTp, run for every family (the `search-diamond` stage).           #
#                                                                                                     #
# Round 1: the family's foreign model vs the proteome, DIAMOND_SENSITIVITY, E-value and coverage.      #
#   * Maximaize sensitivity to pick distant candidates candidates. Round-2 SEEDS are chosen as:       #
#   * HMM-validated  — round-1 subjects that ALSO pass the family hmmscan (two independent methods    #
#                      agreeing = confident members, no false-positive seeds to amplify); used        #
#                      whenever the family has an HMM.                                                #
#   * Blast Score Ratio (BSR) — fallback for families with no HMM: round-1 subjects whose             #
#                      bitscore / model self-bitscore >= DIAMOND_BSR (length-normalized membership,   #
#                      Rasko et al. 2005).                                                            #
# Round 2: those native seeds vs the proteome (identity + stricter coverage) recover the family.      #
#                                                                                                     #
#######################################################################################################
"""

from __future__ import annotations

import csv
import tempfile
from datetime import datetime
from pathlib import Path

from . import external, io
from .config import Config
from .schema import HIT_HEADER

# DIAMOND --outfmt 6 column layout, named so field lookups are self-explanatory.
DIAMOND_FIELDS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]
_COL = {name: i for i, name in enumerate(DIAMOND_FIELDS)}


def _sensitivity_flag(cfg) -> str | None:
    """The DIAMOND sensitivity mode flag (e.g. ``--ultra-sensitive``) from
    DIAMOND_SENSITIVITY, applied to both rounds. Returns None for ``fast``/default
    (DIAMOND's fast mode takes no flag). ``ultra-sensitive`` (the default) matches
    NCBI BLASTP sensitivity."""
    mode = str(cfg.DIAMOND_SENSITIVITY or "").strip().lstrip("-")
    if not mode or mode in ("fast", "default", "none"):
        return None
    return f"--{mode}"


def _round(cfg, query, db, out_path, sensitivity=None, identity=None, coverage=None) -> int:
    """One DIAMOND blastp round. Always filtered by E-value; identity and
    query-coverage are added only when given. ``sensitivity`` is a DIAMOND mode
    flag (e.g. ``--ultra-sensitive``) or None."""
    external.require("diamond")
    cmd = [
        "diamond", "blastp",
        "--db", db,
        "--query", query,
        "--out", out_path,
        "--outfmt", "6", *DIAMOND_FIELDS,
        "--evalue", cfg.DIAMOND_EVALUE,
        "--threads", cfg.THREADS,
    ]
    if sensitivity:
        cmd.append(sensitivity)
    if identity is not None:
        cmd += ["--id", identity]
    if coverage is not None:
        cmd += ["--query-cover", coverage]
    external.run(cmd)
    with open(out_path) as fh:
        return sum(1 for _ in fh)


def _round1_subjects(r1_out) -> set:
    """Unique subject protein ids from a round-1 hit table."""
    subjects = set()
    with open(r1_out) as fh:
        for line in fh:
            if line.strip():
                subjects.add(line.split("\t")[_COL["sseqid"]])
    return subjects


def _hmm_validated_seeds(r1_out, hmm_ids) -> list:
    """Seeds = round-1 subjects that ALSO pass the family HMM (hmmscan). Two
    independent methods agreeing gives high-confidence members and no false-positive
    seeds to amplify in round 2."""
    return sorted(_round1_subjects(r1_out) & set(hmm_ids))


def _self_bitscores(cfg, model_fasta, tmp_dir) -> dict:
    """Per-model-sequence self-alignment bitscore, for the Blast Score Ratio."""
    external.require("diamond")
    self_db = tmp_dir / "self_db"
    external.run(["diamond", "makedb", "--in", model_fasta, "--db", self_db,
                  "--threads", cfg.THREADS])
    self_out = tmp_dir / "self.tsv"
    external.run(["diamond", "blastp", "--db", self_db, "--query", model_fasta,
                  "--out", self_out, "--outfmt", "6", *DIAMOND_FIELDS,
                  "--threads", cfg.THREADS])
    scores = {}
    with open(self_out) as fh:
        for line in fh:
            if not line.strip():
                continue
            c = line.split("\t")
            if c[_COL["qseqid"]] == c[_COL["sseqid"]]:      # the self hit
                q = c[_COL["qseqid"]]
                scores[q] = max(scores.get(q, 0.0), float(c[_COL["bitscore"]]))
    return scores


def _bsr_seeds(r1_out, self_scores, threshold) -> list:
    """Seeds for a family with no HMM: round-1 subjects whose Blast Score Ratio
    (bitscore / the model's self-bitscore) is >= threshold — a length-normalized,
    absolute membership criterion."""
    best = {}
    with open(r1_out) as fh:
        for line in fh:
            if not line.strip():
                continue
            c = line.split("\t")
            self_sc = self_scores.get(c[_COL["qseqid"]])
            if not self_sc:
                continue
            bsr = float(c[_COL["bitscore"]]) / self_sc
            member = c[_COL["sseqid"]]
            best[member] = max(best.get(member, 0.0), bsr)
    return sorted(m for m, bsr in best.items() if bsr >= threshold)


def _seqkit_extract(pattern_file, fasta_in, fasta_out) -> None:
    external.require("seqkit")
    external.run([
        "seqkit", "grep",
        "--pattern-file", pattern_file,
        fasta_in,
        "--out-file", fasta_out,
    ])


def _best_per_member(r2_out) -> dict:
    """Best HSP (highest bitscore) per proteome family member (sseqid)."""
    best = {}
    with open(r2_out) as fh:
        for line in fh:
            if not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            member = c[_COL["sseqid"]]
            bitscore = float(c[_COL["bitscore"]])
            if member not in best or bitscore > best[member][1]:
                best[member] = (c[_COL["evalue"]], bitscore,
                                c[_COL["sstart"]], c[_COL["send"]])
    return best


def _process_family(cfg, family, blast_model, proteome_db, tmp_dir, writer, hmm_ids) -> None:
    model_fasta = cfg.blast_dir / blast_model
    if not model_fasta.exists():
        raise FileNotFoundError(f"model FASTA not found for {family}: {model_fasta}")

    # Round 1: foreign model vs proteome, DIAMOND_SENSITIVITY, E-value only.
    external.log(f"[{datetime.now()}] [{family}] Round 1: model vs proteome "
                 f"({cfg.DIAMOND_SENSITIVITY}, E<={cfg.DIAMOND_EVALUE})...")
    r1_out = cfg.result(f"diamond_{family}_r1.tsv")
    r1_count = _round(cfg, model_fasta, proteome_db, r1_out, sensitivity=_sensitivity_flag(cfg))
    external.log(f"[OK] {family} Round 1: {r1_count} hits")
    if r1_count == 0:
        external.log(f"[WARN] {family}: no Round 1 hits, skipping Round 2.")
        return

    # Seeds: HMM-validated (round-1 subjects that also pass hmmscan) when the family
    # has an HMM; otherwise Blast Score Ratio.
    if hmm_ids:
        seeds = _hmm_validated_seeds(r1_out, hmm_ids)
        seed_desc = "HMM-validated"
    else:
        self_scores = _self_bitscores(cfg, model_fasta, tmp_dir)
        seeds = _bsr_seeds(r1_out, self_scores, cfg.DIAMOND_BSR)
        seed_desc = f"BSR>={cfg.DIAMOND_BSR}"
    external.log(f"[OK] {family} Round 1: {len(seeds)} seed sequences ({seed_desc})")
    if not seeds:
        external.log(f"[WARN] {family}: no confident seeds ({seed_desc}); skipping Round 2.")
        return

    seed_ids = tmp_dir / f"{family}_seed_ids.txt"
    seed_ids.write_text("\n".join(seeds) + "\n")
    seed_fasta = tmp_dir / f"{family}_seeds.fasta"
    _seqkit_extract(seed_ids, cfg.proteome, seed_fasta)

    # Round 2: native seeds vs proteome (identity + stricter coverage).
    external.log(f"[{datetime.now()}] [{family}] Round 2: seeds vs proteome "
                 f"({cfg.DIAMOND_SENSITIVITY}, id {cfg.DIAMOND_IDENTITY}%, "
                 f"cover {cfg.DIAMOND_COVERAGE_R2}%)...")
    r2_out = cfg.result(f"diamond_{family}_r2.tsv")
    _round(cfg, seed_fasta, proteome_db, r2_out,
           sensitivity=_sensitivity_flag(cfg),
           identity=cfg.DIAMOND_IDENTITY, coverage=cfg.DIAMOND_COVERAGE_R2)
    # Only the member IDs (+ coords) are needed for blast_hits.tsv; merge extracts
    # the sequences from the proteome downstream, so no FASTA is written here.
    best = _best_per_member(r2_out)
    external.log(f"[OK] {family}: {len(best)} family members")

    for member, (evalue, bitscore, sstart, send) in sorted(best.items()):
        writer.writerow([member, family, "-", evalue, bitscore, sstart, send, "blast"])


def _hmm_hits_by_family(cfg) -> dict:
    """{family: set(protein_ids)} from intermediate/hmm_hits.tsv (empty if search-hmm
    has not run). Used to HMM-validate the DIAMOND round-2 seeds."""
    path = cfg.result("hmm_hits.tsv")
    if not path.exists():
        return {}
    df = io.read_tsv(path)
    return {str(fam): set(sub["protein_id"].astype(str))
            for fam, sub in df.groupby("family")}


def run(cfg: Config) -> None:
    """Run two-round DIAMOND for every family into intermediate/blast_hits.tsv."""
    cfg.ensure_dirs()
    proteome_db = cfg.proteome_db
    if not Path(f"{proteome_db}.dmnd").exists():
        raise FileNotFoundError(
            f"DIAMOND proteome db not found: {proteome_db}.dmnd (run `gwiscan setup-db` first)"
        )

    hmm_by_family = _hmm_hits_by_family(cfg)
    blast_out = cfg.result("blast_hits.tsv")
    with open(blast_out, "w", newline="") as hits_fh, \
            tempfile.TemporaryDirectory(prefix="gwiscan_diamond_") as tmp:
        writer = csv.writer(hits_fh, delimiter="\t")
        writer.writerow([io.to_camel(c) for c in HIT_HEADER])

        external.log(f"[{datetime.now()}] Two-round DIAMOND BLASTp (model -> proteome, seeds -> proteome)...")
        external.log(
            f"[{datetime.now()}] E-value: {cfg.DIAMOND_EVALUE} | "
            f"Sensitivity: {cfg.DIAMOND_SENSITIVITY} (both rounds) | "
            f"Round 1: E-value only | "
            f"Seeds: HMM-validated, else BSR>={cfg.DIAMOND_BSR} | "
            f"Round 2: id {cfg.DIAMOND_IDENTITY}%, cover {cfg.DIAMOND_COVERAGE_R2}%"
        )
        for r in io.family_records(cfg.family_map):
            _process_family(cfg, r["family"], r["blast_model"], proteome_db,
                            Path(tmp), writer, hmm_by_family.get(r["family"], set()))

    external.log(f"[{datetime.now()}] DIAMOND step done.")
