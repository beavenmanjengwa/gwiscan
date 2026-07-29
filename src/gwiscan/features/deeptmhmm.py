#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# deeptmhmm.py - Transmembrane topology via DeepTMHMM (biolib CLI) (the `deeptmhmm` stage).         #
#                                                                                                  #
# DEEPTMHMM_VERSION pins where it runs: 1.0.24 runs locally in seconds; newer tags queue on         #
# biolib's cloud (same predictor). The app takes only --fasta and writes into biolib_results/ in    #
# the cwd, so it runs with cwd set to the stage output dir. Parsing: predicted_topologies.3line ->  #
# per-protein topology class (GLOB/SP/TM/SP+TM/BETA); TMRs.gff3 -> signal-peptide, TM-helix         #
# (+ n_tm_regions), and beta-strand coordinate ranges. Raw files kept in intermediate/deeptmhmm_out/.    #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

import re

from .. import external, io
from ..config import Config

OUT_HEADER = ["protein_id", "topology", "signal_peptide", "tm_regions",
              "n_tm_regions", "beta_regions"]

TOPOLOGY_FILE = "predicted_topologies.3line"
GFF_FILE = "TMRs.gff3"

_LENGTH_RE = re.compile(r"#\s+(\S+)\s+Length:\s+(\d+)")


def _fmt(segments: list) -> str:
    """Render coordinate ranges as '265-285' (';'-joined), or '-' if none."""
    return ";".join(f"{s}-{e}" for s, e in segments) or "-"


def parse_type_labels(path) -> dict:
    """Map protein_id -> DeepTMHMM class label from the .3line headers.

    Header form is ``>{protein_id} | {TYPE}`` (e.g. ``>Medtr0015s0030.1 | SP+TM``).
    """
    labels = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            pid = line[1:].split("|")[0].split()[0]
            labels[pid] = line.split("|", 1)[1].strip() if "|" in line else ""
    return labels


def parse_gff(path) -> dict:
    """Parse TMRs.gff3 into per-protein coordinate segments.

    Returns ``{protein_id: {"length": int, "signal": [(s, e)], "tm": [...],
    "beta": [...]}}``. Region rows are tab-separated ``id, type, start, end``;
    the type may contain a space (e.g. "Beta sheet"), so we split on tabs.
    """
    data = {}

    def record(pid):
        return data.setdefault(pid, {"length": 0, "signal": [], "tm": [], "beta": []})

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line == "//":
                continue
            if line.startswith("#"):
                m = _LENGTH_RE.match(line)
                if m:
                    record(m.group(1))["length"] = int(m.group(2))
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                continue
            pid, rtype = cols[0], cols[1].strip().lower()
            seg = (int(cols[2]), int(cols[3]))
            rec = record(pid)
            if rtype == "signal":
                rec["signal"].append(seg)
            elif rtype == "tmhelix":
                rec["tm"].append(seg)
            elif "beta" in rtype:
                rec["beta"].append(seg)
    return data


def build_rows(gff_path, three_line_path=None) -> list:
    """Combine gff coordinates with 3line type labels into output rows."""
    coords = parse_gff(gff_path)
    labels = parse_type_labels(three_line_path) if three_line_path else {}
    rows = []
    for pid, rec in coords.items():
        rows.append([
            pid,
            labels.get(pid, ""),
            _fmt(rec["signal"]),
            _fmt(rec["tm"]),
            len(rec["tm"]),
            _fmt(rec["beta"]),
        ])
    return rows


def run(cfg: Config) -> None:
    cfg.ensure_dirs()
    cand_fasta = cfg.result("final_candidates.fasta")
    out_dir = cfg.result("deeptmhmm_out")
    out_tsv = cfg.result("deeptmhmm.tsv")
    out_dir.mkdir(parents=True, exist_ok=True)

    if not cand_fasta.exists():
        raise FileNotFoundError(f"final_candidates.fasta not found: {cand_fasta} (run `confirm` first)")
    external.require("biolib")

    app = "DTU/DeepTMHMM"
    if cfg.DEEPTMHMM_VERSION:
        app += f":{cfg.DEEPTMHMM_VERSION}"

    external.log(f"[deeptmhmm] Running {app}...")
    # No --output flag; biolib writes into biolib_results/ under cwd, so run in
    # out_dir. Pass an absolute FASTA path since cwd is changed.
    external.run(
        ["biolib", "run", app, "--fasta", cand_fasta.resolve()],
        cwd=out_dir,
    )

    gff = sorted(out_dir.rglob(GFF_FILE))
    if not gff:
        raise FileNotFoundError(
            f"DeepTMHMM produced no {GFF_FILE} under {out_dir} "
            f"(check the biolib output above)"
        )
    three = sorted(out_dir.rglob(TOPOLOGY_FILE))

    rows = build_rows(gff[0], three[0] if three else None)
    io.write_tsv(out_tsv, OUT_HEADER, rows)
    external.log(f"[OK] DeepTMHMM parsed: {len(rows)} entries -> deeptmhmm.tsv")
