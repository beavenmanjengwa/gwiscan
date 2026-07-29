"""Tests for the casing boundary: camelCase on disk, snake_case in memory."""

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
