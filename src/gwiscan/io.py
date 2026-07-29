#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# io.py - Shared file I/O: family table, TSV/CSV/XLSX read+write, FASTA loading.                   #
#                                                                                                  #
# Output files use camelCase headers (proteinId, ...); in memory the pipeline keeps snake_case,    #
# so casing is converted only at the read/write boundary (like camelCase JSON on the wire,         #
# snake_case in the code).                                                                         #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Columns whose casing must be preserved verbatim (TargetP class abbreviations
# that fit neither snake_case nor camelCase).
_PRESERVE: set[str] = {"SP", "mTP", "cTP", "luTP", "noTP"}


def to_camel(name: str) -> str:
    """snake_case / 'Title case' -> camelCase, for output headers."""
    name_str = str(name)
    if name_str in _PRESERVE:
        return name_str
    parts = re.split(r"[_ ]+", name_str)
    first, *rest = parts
    return first[:1].lower() + first[1:] + "".join(p[:1].upper() + p[1:] for p in rest)


def to_snake(name: str) -> str:
    """camelCase -> snake_case, to normalise headers read back from files."""
    name_str = str(name)
    if name_str in _PRESERVE:
        return name_str
    return re.sub(r"([A-Z])", lambda m: "_" + m.group(1).lower(), name_str).lstrip("_")


def read_family_map(path: Path | str) -> list[dict[str, str]]:
    """Return the family table (config/family.tsv, or config/superfamily.tsv in
    superfamily mode) as a list of row dicts (not casing-converted).

    Blank lines and lines starting with '#' are ignored, so the family table can
    carry header comments documenting its columns.
    """
    with open(Path(path), "r", encoding="utf-8") as fh:
        rows = [ln for ln in fh if ln.strip() and not ln.lstrip().startswith("#")]
    return list(csv.DictReader(rows, delimiter="\t"))


def _clean(value: Any) -> str:
    """Strip a cell; treat '-' and '' as empty (the map's 'no value' marker)."""
    v = (str(value) if value is not None else "").strip()
    return "" if v == "-" else v


_PFAM_RE = re.compile(r"^PF\d{4,6}$", re.IGNORECASE)


def family_records(path: Path | str) -> list[dict[str, Any]]:
    """Normalised family-table rows — the single source of truth for the schema.

    The ``PfamModel`` cell is one of three things, told apart by its form:
      * a **Pfam accession** (``PF00139``) — the HMM is downloaded from InterPro,
        pressed for identification, and the *same accession* is the InterProScan
        confirmation / coordinate key.
      * a **custom HMM file** (``*.hmm``, e.g. ``CRA.hmm``, provided in db/hmm/) —
        a user-built, family-specific profile. It identifies on its own; its
        hmmscan hit gives the domain coordinates and *is* the confirmation, so it
        needs no Pfam and never touches InterProScan. Use this for families whose
        Pfam is not specific enough to identify by (e.g. CRA's PF00704, the shared
        GH18 enzyme domain, which also matches active chitinases).
      * **empty** (``-``) — no HMM at all; the family is DIAMOND-only (EUL).

    Each record has:
      * ``family`` / ``superfamily`` — names (superfamily present only in that mode).
      * ``pfam_model``   — bare Pfam accession, or '' (custom HMM and DIAMOND-only
                           families have no confirm/coordinate Pfam).
      * ``hmm_file``     — filename pressed into the HMM db (``<acc>.hmm`` for an
                           accession, the given file for a custom HMM), or ''.
      * ``hmm_is_custom``— True when ``hmm_file`` is a user-provided profile (not
                           downloaded, no InterProScan involvement).
      * ``blast_model``  — DIAMOND round-1 query FASTA (DIAMOND runs for every
                           family, so this is required for all).
      * ``hmm_press``    — whether ``hmm_file`` is pressed as an identifying model:
                           has an HMM and ``HmmPress`` is not 'no'. (A Pfam kept
                           only for confirmation, like PF00704, sets HmmPress=no.)
    """
    out = []
    for row in read_family_map(path):
        model = _clean(row.get("PfamModel"))
        press_flag = _clean(row.get("HmmPress")).lower() not in ("no", "false", "0")
        family = (str(row.get("Family") or row.get("FamilyName") or "")).strip()

        if model.lower().endswith(".hmm"):
            pfam, hmm_file, is_custom = "", model, True
        elif model:
            pfam = model.split(".")[0]              # versionless accession
            hmm_file, is_custom = f"{pfam}.hmm", False
        else:
            pfam, hmm_file, is_custom = "", "", False

        out.append({
            "family": family,
            "superfamily": _clean(row.get("Superfamily")) or None,
            "pfam_model": pfam,
            "hmm_file": hmm_file,
            "hmm_is_custom": is_custom,
            "blast_model": _clean(row.get("BlastModel")),
            "hmm_press": bool(hmm_file) and press_flag,
        })
    return out


def pfam_to_family(path: Path | str) -> dict[str, str]:
    """Bare Pfam accession (PF00139) -> family (Legume), for families whose
    PfamModel is a Pfam accession. Includes CRA's PF00704 when CRA is kept on the
    Pfam route (confirmation/coordinates); custom-HMM families have no accession
    and are absent here (they confirm + get coordinates via their own hmmscan hit).
    """
    return {
        str(r["pfam_model"]): str(r["family"])
        for r in family_records(path) 
        if r["pfam_model"]
    }


def custom_hmm_families(path: Path | str) -> set[str]:
    """Families identified by a user-provided custom HMM (no Pfam / InterProScan;
    coordinates come from their own hmmscan hit).
    """
    return {str(r["family"]) for r in family_records(path) if r["hmm_is_custom"]}


def family_to_superfamily(path: Path | str) -> dict[str, str | None]:
    """family -> superfamily, for families that declare one (superfamily mode)."""
    return {
        str(r["family"]): r["superfamily"]
        for r in family_records(path) 
        if r["superfamily"]
    }


def write_tsv(path: Path | str, header: list[str], rows: list[list[Any]]) -> None:
    """Write a TSV; the header (snake_case) is emitted as camelCase."""
    file_path = Path(path)
    file_path.parent.mkdir(parents=True, exist_ok=True)
    with open(file_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([to_camel(c) for c in header])
        writer.writerows(rows)


def write_df(df: pd.DataFrame, path: Path | str, fmt: str) -> None:
    """Write a DataFrame with camelCase headers as tsv / csv / xlsx."""
    out = df.rename(columns=to_camel)
    if fmt == "tsv":
        out.to_csv(Path(path), sep="\t", index=False)
    elif fmt == "csv":
        out.to_csv(Path(path), index=False)
    elif fmt == "xlsx":
        out.to_excel(Path(path), index=False, engine="openpyxl")
    else:
        raise ValueError(f"unknown output format: {fmt!r}")


def read_tsv(path: Path | str, **kwargs: Any) -> pd.DataFrame:
    """Read a pipeline TSV, normalising camelCase headers back to snake_case."""
    df = pd.read_csv(Path(path), sep="\t", **kwargs)
    return df.rename(columns=to_snake)


def load_proteome(path: Path | str) -> dict[str, SeqRecord]:
    """Load a FASTA into an id -> SeqRecord dict."""
    return SeqIO.to_dict(SeqIO.parse(str(path), "fasta"))