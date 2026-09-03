#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# architecture.py - Domain-COMBINATION identification (MODE: architecture).                         #
#                                                                                                  #
# The family/multi-family modes call a protein per domain hit. Architecture mode classifies a       #
# protein by a COMBINATION of Pfam domains on one chain: a PRIMARY domain that defines and seeds    #
# the family, plus one or more REQUIRED domains that must also be present. Both are Pfam HMMs.       #
#                                                                                                  #
# The search is two hmmsearch passes, so the required domains are never scanned genome-wide:          #
#   1. hmmsearch the whole proteome against the PRIMARY HMM(s)  -> candidate proteins.                #
#   2. hmmsearch only those candidates against the primary+required HMMs -> keep the ones that also   #
#      carry every required domain. Those are the final candidates.                                #
#                                                                                                  #
# The final candidates then flow into the ordinary annotation pipeline unchanged (InterProScan for  #
# the complete annotation, ProtParam, TargetP, DeepLoc, coordinates, compile, ...). family here is  #
# the Architecture name, so every downstream stage groups by it.                                   #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

import re
from collections import defaultdict
from datetime import datetime

from Bio import SeqIO

from . import external, hmm, io
from .config import Config
from .schema import HIT_HEADER

_PFAM_RE = re.compile(r"^PF\d{4,6}$", re.IGNORECASE)

# hmmsearch --domtblout fixed columns (0-based); mirror hmm.py. hmmsearch's target
# is the protein and its query is the HMM (the reverse of hmmscan).
_COL = {"protein": 0, "hmm_name": 3, "hmm_acc": 4,
        "i_evalue": 12, "dom_score": 13, "ali_from": 17, "ali_to": 18}
_MIN_COLS = 22


def _bare(acc) -> str:
    """PF00069.27 -> PF00069 (versionless, upper-cased, whitespace-stripped)."""
    return str(acc).strip().split(".")[0].upper()


def _slot(token: str) -> list[str]:
    """One requirement slot ('PF00069' or 'PF00069|PF07714') -> its bare Pfam
    alternatives. A slot is satisfied when ANY of its accessions is present."""
    return [_bare(a) for a in token.split("|") if a.strip()]


def _slots(cell: str) -> list[list[str]]:
    """A '+'-separated cell -> a list of OR-slots (AND across the slots)."""
    return [s for s in (_slot(tok) for tok in str(cell).split("+") if tok.strip()) if s]


def read_rules(path) -> list[dict]:
    """Parse config/architecture.tsv into rule dicts.

    Columns (TAB-separated; '#' comment lines and blanks ignored):
      * Architecture — the family name, verbatim in outputs. Required, unique.
      * Primary      — the defining Pfam domain hmmsearch seeds the search with. One
                       slot; alternatives allowed with '|'. Required.
      * Required     — the Pfam domain(s) that must ALSO be present, '+'-separated
                       (AND); each slot may list '|' alternatives (OR). Required.
      * Class        — optional rollup label grouping architectures. Defaults to the
                       Architecture name.

    Both Primary and Required are Pfam accessions (PFxxxxx) — the HMMs hmmsearch
    searches with.
    """
    rules, seen = [], set()
    for row in io.read_family_map(path):
        name = (str(row.get("Architecture") or row.get("Name") or "")).strip()
        if not name:
            raise ValueError(f"architecture table {path}: a row has no Architecture name")
        if name in seen:
            raise ValueError(f"architecture table {path}: duplicate Architecture {name!r}")
        seen.add(name)

        primary = _slots(row.get("Primary") or "")
        if len(primary) != 1:
            raise ValueError(
                f"architecture table {path}: {name!r} needs exactly one Primary domain "
                f"slot (got {row.get('Primary')!r})"
            )
        required = _slots(row.get("Required") or "")
        if not required:
            raise ValueError(
                f"architecture table {path}: {name!r} has no Required domains "
                f"(a combination needs the primary plus at least one required domain)"
            )

        cls = (str(row.get("Class") or "").strip() or "-")
        cls = name if cls == "-" else cls

        for acc in primary[0] + [a for slot in required for a in slot]:
            if not _PFAM_RE.match(acc):
                raise ValueError(
                    f"architecture table {path}: {name!r} references {acc!r}, which is not a "
                    f"Pfam accession (PFxxxxx). Primary and Required must name Pfam accessions."
                )

        rules.append({"architecture": name, "class": cls,
                      "primary": primary[0], "required": required})
    if not rules:
        raise ValueError(f"architecture table {path}: no architecture rows found")
    return rules


def _accessions(rules, which) -> list[str]:
    """Ordered, de-duplicated bare Pfam accessions for a slot kind."""
    seen, out = set(), []
    for rule in rules:
        slots = [rule["primary"]] if which == "primary" else rule["required"]
        if which == "all":
            slots = [rule["primary"]] + rule["required"]
        for acc in [a for slot in slots for a in slot]:
            if acc not in seen:
                seen.add(acc)
                out.append(acc)
    return out


def primary_accessions(rules) -> list[str]:
    """Pfam accessions used to SEED the genome-wide search (pass 1)."""
    return _accessions(rules, "primary")


def all_accessions(rules) -> list[str]:
    """Every Pfam accession (primary + required) — searched on candidates (pass 2)."""
    return _accessions(rules, "all")


def arch_names(rules) -> list[str]:
    return [r["architecture"] for r in rules]


def fam_to_pfam(rules) -> dict:
    """Architecture name -> its primary Pfam accession (for compile's architecture flag)."""
    return {r["architecture"]: r["primary"][0] for r in rules}


def _parse_domtbl(path):
    """hmmsearch --domtblout -> list of hit dicts (protein/name/pfam/evalue/bitscore/start/end)."""
    hits = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c = line.split()
            if len(c) < _MIN_COLS:
                continue
            hits.append({
                "protein": c[_COL["protein"]],
                "name": c[_COL["hmm_name"]],
                "pfam": _bare(c[_COL["hmm_acc"]]),
                "acc": c[_COL["hmm_acc"]],
                "evalue": c[_COL["i_evalue"]],
                "bitscore": c[_COL["dom_score"]],
                "start": int(c[_COL["ali_from"]]),
                "end": int(c[_COL["ali_to"]]),
            })
    return hits


def _hmmsearch(cfg: Config, hmm_db, query_fasta, tag):
    """Run hmmsearch (the hmm_db profiles vs the query_fasta sequences) with the run's
    cutoff (Pfam GA thresholds by default, or -E when HMM_EVALUE is set); return the
    domtbl path. hmmsearch saturates the cores, unlike hmmscan, on a genome-wide search."""
    domtbl = cfg.result(f"arch_{tag}_domtbl.txt")
    external.log(f"[{datetime.now()}] hmmsearch ({tag}) ...")
    external.run([
        "hmmsearch", *hmm.threshold_args(cfg), "--noali", "--cpu", cfg.THREADS,
        "--domtblout", domtbl,
        "-o", cfg.result(f"arch_{tag}.out"),
        hmm_db, query_fasta,
    ])
    return domtbl


def classify(candidate_hits, rules):
    """From the pass-2 hits on candidates, return HIT_HEADER rows for each matched
    (protein, architecture): the member's primary+required domain hits, relabeled
    with family = the architecture. A protein matches a rule when its Pfam set holds
    the primary and every required slot."""
    by_protein = defaultdict(list)
    for h in candidate_hits:
        by_protein[h["protein"]].append(h)

    rows = []
    for pid, hits in by_protein.items():
        present = {h["pfam"] for h in hits}
        for rule in rules:
            slots = [rule["primary"]] + rule["required"]
            if not all(any(a in present for a in slot) for slot in slots):
                continue
            wanted = {a for slot in slots for a in slot}
            for h in sorted(hits, key=lambda x: x["start"]):
                if h["pfam"] in wanted:
                    rows.append([pid, rule["architecture"], h["acc"], h["evalue"],
                                 h["bitscore"], h["start"], h["end"], "hmm"])
    return rows


def preflight(cfg: Config) -> None:
    """Architecture-mode pre-flight: proteome, hmmsearch, a parseable rules table
    naming only Pfam accessions, and the InterProScan requirement (it annotates the
    final candidates)."""
    external.log("[preflight] GWIscan architecture-mode pre-flight check")
    external.log("--------------------------------------")
    fail = False

    for binary in ("hmmsearch",):
        if external.available(binary):
            external.log(f"[OK] {binary}")
        else:
            external.log(f"[MISSING] {binary} not found on PATH")
            fail = True

    # Proteome(s): multi-species runs validate every manifest species; single-species
    # validates input/proteome.fasta. Shared with the family/multi-family pre-flight.
    from . import preflight
    fail = preflight._check_proteomes(cfg) or fail

    if not cfg.architecture_map.exists():
        external.log(f"[MISSING] architecture table not found: {cfg.architecture_map}")
        fail = True
    else:
        try:
            rules = read_rules(cfg.architecture_map)
            external.log(f"[OK] architecture table: {len(rules)} rule(s); "
                         f"primary {primary_accessions(rules)}, "
                         f"required {_accessions(rules, 'required')}")
        except (ValueError, OSError) as e:
            external.log(f"[MISSING] architecture table invalid: {e}")
            fail = True

    # InterProScan annotates the final candidates, so its requirement is checked here.
    ipr = str(cfg.INTERPRO_MODE).lower()
    if ipr == "local" and not external.available(cfg.INTERPROSCAN_BIN):
        external.log(f"[MISSING] local InterProScan not found: {cfg.INTERPROSCAN_BIN}")
        fail = True
    elif ipr == "api" and not cfg.EBI_EMAIL:
        external.log("[MISSING] INTERPRO_MODE=api needs EBI_EMAIL (or use INTERPRO_MODE=local)")
        fail = True
    elif ipr not in ("local", "api"):
        external.log(f"[MISSING] unknown INTERPRO_MODE {cfg.INTERPRO_MODE!r} (use 'api' or 'local')")
        fail = True
    else:
        external.log("[OK] InterProScan requirement satisfied (annotates final candidates)")

    if fail:
        raise RuntimeError("pre-flight check failed — see messages above")
    external.log("[PASS] Pre-flight check passed.")


def run(cfg: Config) -> None:
    """Two-stage hmmsearch identification of the domain-combination members.

    Pass 1 seeds on the primary domain across the whole proteome; pass 2 searches
    only those candidates against the primary+required HMMs and keeps the ones that
    carry every required domain. Writes the standard intermediate files (family =
    the architecture) so the ordinary annotation pipeline runs on the members:
      candidates.fasta / final_candidates.fasta  — member sequences
      candidates_merged.tsv / final_candidates.tsv — member domain hits (HIT_HEADER)
    """
    cfg.ensure_dirs()
    external.require("hmmsearch")
    if not cfg.proteome.exists():
        raise FileNotFoundError(f"proteome not found: {cfg.proteome}")
    if not cfg.primary_hmm_db.exists() or not cfg.hmm_db.exists():
        raise FileNotFoundError("HMM databases not built (run `setup-db` first)")

    rules = read_rules(cfg.architecture_map)
    # Index the proteome (file offsets only) rather than loading every record into
    # memory -- a plant genome can carry 100k+ proteins.
    proteome = SeqIO.index(str(cfg.proteome), "fasta")
    try:
        # Pass 1: seed on the primary domain across the whole proteome.
        primary_hits = _parse_domtbl(_hmmsearch(cfg, cfg.primary_hmm_db, cfg.proteome, "primary"))
        candidate_ids = list(dict.fromkeys(h["protein"] for h in primary_hits))
        external.log(f"[OK] Pass 1 (primary seed): {len(candidate_ids)} candidate proteins")

        cand_seed_fasta = cfg.result("candidates_primary.fasta")
        SeqIO.write((proteome[i] for i in candidate_ids if i in proteome),
                    str(cand_seed_fasta), "fasta")

        # Pass 2: search candidates against primary+required, then keep full combinations.
        if candidate_ids:
            candidate_hits = _parse_domtbl(_hmmsearch(cfg, cfg.hmm_db, cand_seed_fasta, "required"))
        else:
            candidate_hits = []
        rows = classify(candidate_hits, rules)

        io.write_tsv(cfg.result("candidates_merged.tsv"), HIT_HEADER, rows)
        io.write_tsv(cfg.result("final_candidates.tsv"), HIT_HEADER, rows)

        member_ids = list(dict.fromkeys(r[0] for r in rows))
        members = [proteome[i] for i in member_ids if i in proteome]
        # candidates.fasta == final members: InterProScan and domain extraction both
        # operate on the final candidates only.
        SeqIO.write(members, str(cfg.result("candidates.fasta")), "fasta")
        SeqIO.write(members, str(cfg.result("final_candidates.fasta")), "fasta")
    finally:
        proteome.close()

    per_arch = defaultdict(set)
    for r in rows:
        per_arch[r[1]].add(r[0])
    external.log(f"[OK] Final candidates: {len(member_ids)} proteins across "
                 f"{len(per_arch)} architecture(s) -> final_candidates.tsv/.fasta")
    for arch in sorted(per_arch):
        external.log(f"     {arch}: {len(per_arch[arch])} proteins")
