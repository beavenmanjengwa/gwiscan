#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# test_io.py - Casing-boundary tests: camelCase on disk, snake_case in memory.                     #
#                                                                                                  #
####################################################################################################
"""

from gwiscan import io


def test_to_camel_snake_cases():
    assert io.to_camel("protein_id") == "proteinId"
    assert io.to_camel("molecular_weight_kda") == "molecularWeightKda"
    assert io.to_camel("go_terms") == "goTerms"
    assert io.to_camel("gene_start") == "geneStart"
    assert io.to_camel("Membrane types") == "membraneTypes"
    assert io.to_camel("Localizations") == "localizations"
    assert io.to_camel("gravy") == "gravy"


def test_preserved_targetp_classes():
    # TargetP class abbreviations must not be mangled either way.
    for name in ("SP", "mTP", "cTP", "luTP", "noTP"):
        assert io.to_camel(name) == name
        assert io.to_snake(name) == name


def test_roundtrip_for_referenced_columns():
    # Columns the pipeline joins/reads on must survive camel->snake exactly.
    for snake in ["protein_id", "family", "accession", "evalue", "bitscore",
                  "start", "end", "pfam_id", "ipr_acc", "go_terms",
                  "method", "cs_position", "targetp_type"]:
        assert io.to_snake(io.to_camel(snake)) == snake


def _all_pipeline_columns():
    """Every snake_case column name the pipeline actually emits/reads, gathered from
    the schema/header constants so the round-trip is checked against the real set,
    not a hand-copied sample that can drift."""
    from gwiscan.schema import HIT_HEADER
    from gwiscan.coords import MAP_HEADER
    from gwiscan.score import HEADER as SCORE_HEADER
    from gwiscan.features.interpro import BASE_COLUMNS, OPTIONAL_COLUMNS
    from gwiscan.compile import COLUMN_BANDS

    cols = set(HIT_HEADER) | set(MAP_HEADER) | set(SCORE_HEADER)
    cols |= set(BASE_COLUMNS) | set(OPTIONAL_COLUMNS)
    cols |= {c for _, band in COLUMN_BANDS for c in band}
    return cols


def test_every_pipeline_column_roundtrips():
    # The casing boundary must be lossless for every real column: camel on disk,
    # snake in memory, back to the same snake. Preserved names (SP/mTP/...) are
    # covered separately and pass through unchanged here too.
    for snake in _all_pipeline_columns():
        assert io.to_snake(io.to_camel(snake)) == snake, snake


def test_generated_snake_names_roundtrip():
    # A generative sweep over well-formed snake_case names (1-4 lowercase/numeric
    # words) to guard the boundary beyond the columns in use today.
    words = ["protein", "id", "go", "terms", "n", "tm", "regions", "ec", "value2",
             "aa", "kda", "x"]
    from itertools import product
    for n in (1, 2, 3):
        for combo in product(words, repeat=n):
            snake = "_".join(combo)
            assert io.to_snake(io.to_camel(snake)) == snake, snake


def test_write_is_camel_read_is_snake(tmp_path):
    p = tmp_path / "hits.tsv"
    io.write_tsv(p, ["protein_id", "family", "start"], [["X1", "Legume", "34"]])

    # On disk: camelCase header
    first_line = p.read_text().splitlines()[0]
    assert first_line == "proteinId\tfamily\tstart"

    # Read back: snake_case columns for in-memory use
    df = io.read_tsv(p)
    assert list(df.columns) == ["protein_id", "family", "start"]
    assert df.iloc[0]["protein_id"] == "X1"


# --- InterProScan confirmation: InterProModel column & application derivation ----

def test_application_for_infers_member_db_from_prefix():
    assert io.application_for("PF01453") == "Pfam"
    assert io.application_for("cd02879") == "CDD"
    assert io.application_for("PTHR31257") == "PANTHER"
    # Limited to the three gatekeeper databases: anything else is unrecognised.
    assert io.application_for("SM00108") is None
    assert io.application_for("PS50000") is None
    assert io.application_for("XYZ999") is None


def test_accession_list_parses_pipe_and_strips_version():
    assert io._accession_list("PF01476|PTHR27007") == ["PF01476", "PTHR27007"]
    assert io._accession_list("PF00139.27") == ["PF00139"]     # versionless
    assert io._accession_list("-") == []
    assert io._accession_list("") == []


_FAMILY_TSV = (
    "Family\tPfamModel\tBlastModel\tInterProModel\n"
    "GNA\tPF01453\tgna.fasta\t\n"                       # falls back to PfamModel
    "CRA\t-\tcra.fasta\tcd02879\n"                      # no Pfam -> CDD
    "EUL\t-\teul.fasta\tPTHR31257\n"                    # no Pfam -> PANTHER
    "MULTI\tPF01476\tmulti.fasta\tPF01476|PTHR27007\n"  # confirms on either
)


def _family_map(tmp_path):
    p = tmp_path / "family.tsv"
    p.write_text(_FAMILY_TSV)
    return p


def test_family_confirm_accessions_uses_column_or_pfam_fallback(tmp_path):
    accs = io.family_confirm_accessions(_family_map(tmp_path))
    assert accs["GNA"] == {"PF01453"}                  # fallback to PfamModel
    assert accs["CRA"] == {"cd02879"}
    assert accs["EUL"] == {"PTHR31257"}
    assert accs["MULTI"] == {"PF01476", "PTHR27007"}


def test_interpro_applications_is_the_union(tmp_path):
    apps = io.interpro_applications(_family_map(tmp_path))
    # Pfam (GNA/MULTI), CDD (CRA), PANTHER (EUL/MULTI) -- one run over the union.
    assert set(apps) == {"Pfam", "CDD", "PANTHER"}
