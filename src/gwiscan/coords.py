#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# coords.py - Genomic coordinates for family members (the `coords` stage).                         #
#                                                                                                  #
# Reads the species annotation (GTF or GFF3) that accompanies the input proteome and gives every   #
# confirmed member its chromosome, start, end, strand and gene structure (intron count).           #
# Supply the REPRESENTATIVE annotation -- the same one-transcript-per-gene set the proteome comes   #
# from -- so each protein resolves to one locus and no isoform ambiguity enters the coordinates.    #
#                                                                                                  #
# Writes intermediate/chromosome_map.tsv (feeds chromosomal distribution figures such as MapChart) and  #
# final_results/gwiscan_members.gff3 (the members as genomic features, for browsers and synteny    #
# tools). Each locus carries the annotation feature id it came from alongside the protein id.      #
#                                                                                                  #
####################################################################################################
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from . import external, io
from .config import Config

# Genomic coordinates, named apart from the protein-space start/end (hit schema)
# so the two coordinate systems are never read as one. intron_count is the gene
# structure read from the annotation (intron count = exon count - 1; the exon count
# itself is just intron_count + 1, so only one of the two is stored). Intron number is
# reported because an intronless gene reads as a clear 0.
MAP_HEADER = ["protein_id", "gene_id", "family", "chrom", "gene_start", "gene_end",
              "strand", "intron_count"]

# Attribute keys naming a feature, most specific first: a protein id in the FASTA
# matches a transcript/protein id before it matches the gene it belongs to. `Name`
# is included because Phytozome GFF3 puts the clean transcript id there
# (ID=Cucsa.000200.1.v1.122;Name=Cucsa.000200.1), and that is what proteomes carry.
ID_KEYS = ("transcript_id", "protein_id", "ID", "Name", "gene_id", "Parent")

# Ensembl-style prefixes on GFF3 identifiers (transcript:AT1G01010.1).
ID_PREFIXES = ("transcript", "gene", "cds", "protein", "mrna")


def parse_attributes(field: str) -> dict:
    """Column 9 of a GTF (``key "value";``) or GFF3 (``key=value;``) line."""
    attrs = {}
    for part in str(field).strip().strip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            key, _, value = part.partition("=")
        else:
            key, _, value = part.partition(" ")
        key, value = key.strip(), value.strip().strip('"')
        if key and value:
            attrs.setdefault(key, value)
    return attrs


def strip_prefix(identifier: str) -> str:
    """``transcript:AT1G01010.1`` -> ``AT1G01010.1``; other ids pass through."""
    prefix, sep, rest = identifier.partition(":")
    return rest if sep and prefix.lower() in ID_PREFIXES and rest else identifier


def gene_of(feature_type: str, attrs: dict) -> str:
    """The gene a line belongs to: its ``gene_id`` (GTF), its own ``ID`` when the
    line IS the gene, or its ``Parent`` (GFF3 transcript). '' when none is stated."""
    if attrs.get("gene_id"):
        return strip_prefix(attrs["gene_id"])
    if str(feature_type).lower() == "gene" and attrs.get("ID"):
        return strip_prefix(attrs["ID"])
    if attrs.get("Parent"):
        return strip_prefix(attrs["Parent"].split(",")[0])
    return ""


def parse_annotation(path: Path | str) -> tuple:
    """({feature id: (chrom, start, end, strand)}, {feature id: gene id}).

    Every identifier a line carries is indexed, and a feature's span grows to cover
    each line naming it, so a transcript resolves whether the file states it as one
    ``transcript``/``mRNA`` record or only as its exon/CDS children. The second map
    records which gene each identifier belongs to, so a transcript-level match still
    reports the gene."""
    spans, genes = {}, {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            # Simple 4-column coordinate file: chrom, feature id, start, end (no
            # strand or attributes). Some prepared annotations ship in this form.
            if len(fields) == 4:
                chrom = fields[0]
                ident = strip_prefix(fields[1].strip())
                try:
                    start, end = int(fields[2]), int(fields[3])
                except ValueError:
                    continue
                if ident:
                    previous = spans.get(ident)
                    if previous is None:
                        spans[ident] = (chrom, start, end, ".")
                    elif previous[0] == chrom:
                        spans[ident] = (chrom, min(previous[1], start),
                                        max(previous[2], end), previous[3])
                    genes.setdefault(ident, ident.rsplit(".", 1)[0] if "." in ident else ident)
                continue
            if len(fields) < 9:
                continue
            chrom, strand = fields[0], fields[6]
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                continue
            attrs = parse_attributes(fields[8])
            gene = gene_of(fields[2], attrs)
            for key in ID_KEYS:
                for value in str(attrs.get(key, "")).split(","):
                    identifier = strip_prefix(value.strip())
                    if not identifier:
                        continue
                    previous = spans.get(identifier)
                    if previous is None:
                        spans[identifier] = (chrom, start, end, strand)
                    elif previous[0] == chrom:
                        spans[identifier] = (chrom, min(previous[1], start),
                                             max(previous[2], end), strand)
                    if gene:
                        genes.setdefault(identifier, gene)
    return spans, genes


# Feature types that ARE a transcript (their id is what exon/CDS children name via
# Parent), and the attribute keys under which a transcript may also be known.
_TX_FEATURES = ("mrna", "transcript")
_TX_ALIAS_KEYS = ("ID", "Name", "transcript_id")


def count_exons(path: Path | str) -> dict:
    """{transcript id: exon count} from a GTF/GFF3, for gene-structure reporting.

    Counts ``exon`` features per transcript (a transcript with no exon lines but with
    ``CDS`` lines falls back to its CDS count, so annotations that ship only CDS still
    yield a structure count). A transcript is named by its children through
    ``transcript_id`` (GTF) or ``Parent`` (GFF3); the count is then also propagated to
    the transcript's other identifiers (``Name`` in particular), so a proteome keyed on
    the Phytozome ``Name`` still resolves. Intron number is ``exon count - 1``."""
    exon, cds = defaultdict(int), defaultdict(int)
    aliases: dict = defaultdict(set)
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            ftype = str(fields[2]).lower()
            attrs = parse_attributes(fields[8])
            if ftype in _TX_FEATURES:
                own = strip_prefix(str(attrs.get("ID", "")).strip())
                if own:
                    for key in _TX_ALIAS_KEYS:
                        value = attrs.get(key)
                        if value:
                            aliases[own].add(strip_prefix(str(value).split(",")[0].strip()))
                continue
            if ftype not in ("exon", "cds"):
                continue
            transcript = ""
            for key in ("transcript_id", "Parent"):
                if attrs.get(key):
                    transcript = strip_prefix(str(attrs[key]).split(",")[0].strip())
                    break
            if transcript:
                (exon if ftype == "exon" else cds)[transcript] += 1

    counts = {}
    for transcript in set(exon) | set(cds):
        counts[transcript] = exon[transcript] if exon.get(transcript) else cds.get(transcript, 0)
    # Give a transcript's count to its alternate ids too (e.g. Phytozome Name=), so a
    # proteome keyed on one of those still resolves the gene structure.
    for transcript, alias_set in aliases.items():
        count = counts.get(transcript)
        if count:
            for alias in alias_set:
                counts.setdefault(alias, count)
    return counts


def exons_for(pid: str, gene_id: str, counts: dict):
    """Exon count for a protein, matching the same id resolution ``locate`` uses
    (exact id, then isoform-suffix stripped, then the gene id). None if unknown."""
    for key in (pid, pid.rsplit(".", 1)[0], gene_id):
        if key in counts:
            return counts[key]
    return None


def locate(protein_ids, spans, genes=None) -> tuple:
    """({protein_id: (gene_id, chrom, start, end, strand)}, [unresolved ids]).

    A protein id is looked up as given, then with a trailing isoform suffix removed
    (``AT1G01010.1`` -> ``AT1G01010``) for annotations keyed on the gene. A protein
    id matches a TRANSCRIPT, so the gene reported is the one the annotation states
    for that feature, not the identifier the lookup happened to match on."""
    genes = genes or {}
    located, unresolved = {}, []
    for pid in protein_ids:
        for key in (pid, pid.rsplit(".", 1)[0]):
            hit = spans.get(key)
            if hit:
                located[pid] = (genes.get(key, key), *hit)
                break
        else:
            unresolved.append(pid)
    return located, unresolved


def write_gff3(path: Path, rows) -> None:
    """The members as genomic features, one gene line per protein. Rows carry the
    MAP_HEADER columns; the intron count (when known) is added as an attribute."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("##gff-version 3\n")
        for row in rows:
            pid, gene_id, family, chrom, start, end, strand = row[:7]
            attributes = f"ID={gene_id};protein={pid};family={family}"
            intron_count = row[7] if len(row) > 7 else ""
            if intron_count not in ("", None):
                attributes += f";introns={intron_count}"
            handle.write(f"{chrom}\tGWIscan\tgene\t{start}\t{end}\t.\t{strand}\t.\t"
                         f"{attributes}\n")


def run(cfg: Config) -> None:
    """Map every confirmed member onto the genome from the species annotation."""
    cfg.ensure_dirs()
    annotation = cfg.annotation
    if annotation is None or not annotation.exists():
        external.log("[SKIP] No genome annotation (GTF/GFF3) provided; "
                     "chromosomal coordinates are not available for this run.")
        return

    final = cfg.result("final_candidates.tsv")
    if not final.exists():
        raise FileNotFoundError(
            f"final_candidates.tsv not found: {final} (run `gwiscan confirm` first)")

    external.log(f"[..] Reading annotation: {annotation.name}")
    spans, genes = parse_annotation(annotation)
    exon_counts = count_exons(annotation)
    external.log(f"[OK] Annotation features indexed: {len(spans)}")

    final_df = io.read_tsv(final)
    pairs, seen = [], set()
    for pid, family in zip(final_df["protein_id"].astype(str),
                           final_df["family"].astype(str)):
        if (pid, family) not in seen:
            seen.add((pid, family))
            pairs.append((pid, family))

    located, unresolved = locate([p for p, _ in pairs], spans, genes)
    rows = []
    for pid, family in pairs:
        if pid not in located:
            continue
        gene_id, chrom, start, end, strand = located[pid]
        n_exons = exons_for(pid, gene_id, exon_counts)
        # Report intron number (exon count - 1); blank when the annotation gives no
        # exon/CDS features for the transcript. Exon count itself is intron_count + 1.
        intron_count = "" if n_exons is None else max(n_exons - 1, 0)
        rows.append([pid, gene_id, family, chrom, start, end, strand, intron_count])
    rows.sort(key=lambda r: (r[3], int(r[4])))

    io.write_tsv(cfg.result("chromosome_map.tsv"), MAP_HEADER, rows)
    write_gff3(cfg.final_dir / "gwiscan_members.gff3", rows)

    if unresolved:
        external.log(f"[ERROR] {len(set(unresolved))} protein(s) are absent from the "
                     f"annotation (e.g. {sorted(set(unresolved))[0]}): the annotation and "
                     f"the proteome are not the same gene set.")
    external.log(f"[OK] Chromosomal coordinates: {len(rows)} member loci "
                 f"-> chromosome_map.tsv, gwiscan_members.gff3")
