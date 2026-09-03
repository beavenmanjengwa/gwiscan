#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# score.py - Per-family detectability profile (the `score` stage).                                 #
#                                                                                                  #
# Families differ in how findable they are: a conserved family (Hevein) is recovered by pairwise   #
# BLAST and by a profile HMM alike, while a divergent one (GNA/bulb-type) is recovered by the HMM  #
# and missed by BLAST entirely. Running both searches genome-wide makes that difference            #
# measurable, so it is reported per family instead of assumed.                                     #
#                                                                                                  #
# The two hit sets compared are two independent methods, each as its full result:                   #
#   H - hmmsearch hits (profile HMM vs proteome).                                                    #
#   B - DIAMOND BLAST members (blast_hits.tsv), the full two-round result.                           #
#                                                                                                  #
# Writes intermediate/family_detectability.tsv (raw counts and fractions, always) and logs a per-family #
# verdict naming the concrete action where identification is weak.                                 #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

from . import external, io
from .config import Config

HEADER = [
    "family", "has_hmm", "n_hmm", "n_blast", "n_both", "n_union",
    "jaccard", "hmm_only_frac", "blast_only_frac", "blast_recall_of_hmm",
    "detectability",
]


def blast_hits_by_family(cfg: Config) -> dict:
    """{family: {protein_id, ...}} from blast_hits.tsv — every DIAMOND BLAST member
    (the full two-round result)."""
    path = cfg.result("blast_hits.tsv")
    if not path.exists():
        return {}
    hits = io.read_tsv(path)
    return {str(family): set(sub["protein_id"].astype(str))
            for family, sub in hits.groupby("family")}


def hmm_hits_by_family(cfg: Config) -> dict:
    """{family: {protein_id, ...}} from hmm_hits.tsv."""
    path = cfg.result("hmm_hits.tsv")
    if not path.exists():
        return {}
    hits = io.read_tsv(path)
    return {str(family): set(sub["protein_id"].astype(str))
            for family, sub in hits.groupby("family")}


def family_metrics(hmm_ids, blast_ids, has_hmm: bool, concordance_min: float) -> dict:
    """Detectability metrics and class for one family from its two hit sets.

    ``concordance_min`` is a reporting cutoff on the Jaccard index; every raw
    count and fraction is returned as well, so the classification can be
    reproduced or re-thresholded from the written table.
    """
    hmm_set, blast_set = set(hmm_ids), set(blast_ids)
    both, union = hmm_set & blast_set, hmm_set | blast_set
    n_union = len(union)

    jaccard = len(both) / n_union if n_union else 0.0
    hmm_only = len(hmm_set - blast_set) / n_union if n_union else 0.0
    blast_only = len(blast_set - hmm_set) / n_union if n_union else 0.0
    blast_recall = len(both) / len(hmm_set) if hmm_set else 0.0

    if n_union == 0:
        detectability = "not_detected"
    elif not has_hmm:
        detectability = "blast_only_family"
    elif jaccard >= concordance_min:
        detectability = "concordant"
    elif hmm_only > blast_only:
        detectability = "hmm_dominant"
    else:
        detectability = "blast_dominant"

    return {
        "n_hmm": len(hmm_set),
        "n_blast": len(blast_set),
        "n_both": len(both),
        "n_union": n_union,
        "jaccard": round(jaccard, 4),
        "hmm_only_frac": round(hmm_only, 4),
        "blast_only_frac": round(blast_only, 4),
        "blast_recall_of_hmm": round(blast_recall, 4),
        "detectability": detectability,
    }


def verdict(family: str, metrics: dict) -> str:
    """A log line stating what the family's detectability class means for the run."""
    detectability = metrics["detectability"]
    if detectability == "hmm_dominant":
        return (f"[VERDICT] {family}: HMM-dominant — DIAMOND recovered "
                f"{metrics['blast_recall_of_hmm']:.0%} of the {metrics['n_hmm']} hmmsearch hits. "
                f"A BLAST-only survey of this family would under-report it.")
    if detectability == "blast_dominant":
        return (f"[VERDICT] {family}: BLAST-dominant — {metrics['n_blast'] - metrics['n_both']} "
                f"member(s) found by BLAST alone; the profile may be too strict for this family.")
    if detectability == "blast_only_family":
        return (f"[VERDICT] {family}: no HMM — {metrics['n_blast']} member(s) from BLAST only. "
                f"A curated family HMM (PfamModel = <family>.hmm) would add profile evidence.")
    if detectability == "not_detected":
        return f"[VERDICT] {family}: not detected by either search."
    return (f"[OK] {family}: concordant (Jaccard {metrics['jaccard']:.2f}) — "
            f"both searches recover the same family.")


def run(cfg: Config) -> None:
    """Profile every family's detectability and write family_detectability.tsv."""
    cfg.ensure_dirs()
    external.log("[score] Per-family detectability: hmmsearch vs DIAMOND BLAST")

    by_family = hmm_hits_by_family(cfg)
    blast_by_family = blast_hits_by_family(cfg)
    rows, counts = [], {}

    for record in io.family_records(cfg.family_map):
        family = record["family"]
        metrics = family_metrics(
            by_family.get(family, set()),
            blast_by_family.get(family, set()),
            bool(record["hmm_press"]),
            cfg.CONCORDANCE_MIN,
        )
        rows.append([family, record["hmm_press"]] + [metrics[c] for c in HEADER[2:]])
        counts[metrics["detectability"]] = counts.get(metrics["detectability"], 0) + 1
        external.log(verdict(family, metrics))

    io.write_tsv(cfg.result("family_detectability.tsv"), HEADER, rows)
    summary = ", ".join(f"{n} {klass}" for klass, n in sorted(counts.items()))
    external.log(f"[OK] Detectability profile: {len(rows)} families ({summary}) "
                 f"-> family_detectability.tsv")
