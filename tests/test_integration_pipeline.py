#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# test_integration_pipeline.py - End-to-end wiring test for the core join chain (merge -> confirm  #
# -> compile).                                                                                     #
#                                                                                                  #
# The unit tests cover each stage in isolation; this one seeds the files the search stages would   #
# have produced and runs the real merge, confirm and compile functions in sequence, so a cross-    #
# stage break (a renamed column, a camel/snake slip at the file boundary, a changed join key or    #
# dedup rule) is caught here rather than only in a live run. No external binary or network is      #
# touched.                                                                                         #
#                                                                                                  #
####################################################################################################
"""

from gwiscan import candidates, compile as compile_stage, confirm, io
from gwiscan.config import Config
from gwiscan.features.interpro import BASE_COLUMNS
from gwiscan.schema import HIT_HEADER


def _project(tmp_path):
    """A minimal single-species project: family table + proteome, ready for the
    intermediate files to be seeded under it."""
    (tmp_path / "config").mkdir()
    (tmp_path / "config" / "family.tsv").write_text(
        "Family\tPfamModel\tBlastModel\n"
        "GNA\tPF01453\tgna.fasta\n"      # Pfam family -> InterProScan confirms it
        "EUL\t-\teul.fasta\n"            # BLAST-only family -> kept without a Pfam
    )
    (tmp_path / "input").mkdir()
    (tmp_path / "input" / "proteome.fasta").write_text(
        ">P_gna kinase\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ\n"   # GNA: hmm + blast
        ">P_gna2\nMSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSG\n"        # GNA: hmm only
        ">P_eul\nMAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ\n"          # EUL: blast only
    )
    return Config(root=tmp_path)


def _seed_search_hits(cfg):
    """hmm_hits.tsv + blast_hits.tsv as the two search stages would write them."""
    io.write_tsv(cfg.result("hmm_hits.tsv"), HIT_HEADER, [
        ["P_gna", "GNA", "PF01453.1", "1e-40", "130.0", "10", "110", "hmm"],
        ["P_gna2", "GNA", "PF01453.1", "1e-25", "95.0", "12", "108", "hmm"],
    ])
    io.write_tsv(cfg.result("blast_hits.tsv"), HIT_HEADER, [
        ["P_gna", "GNA", "-", "1e-30", "120.0", "8", "112", "blast"],
        ["P_eul", "EUL", "-", "1e-18", "80.0", "5", "100", "blast"],
    ])


def _seed_interproscan(cfg):
    """interproscan.tsv with a Pfam PF01453 hit on both GNA proteins (so confirm
    keeps them) and none on the EUL protein. Includes the optional go_terms column,
    as a real --goterms run does."""
    header = BASE_COLUMNS + ["go_terms"]
    def row(pid):
        return [pid, "md5", "120", "Pfam", "PF01453", "Bulb-type lectin domain",
                "10", "110", "1e-30", "T", "2026-01-01", "IPR001480",
                "Bulb-type lectin", "-"]
    io.write_tsv(cfg.result("interproscan.tsv"), header, [row("P_gna"), row("P_gna2")])


def _seed_coords(cfg):
    """chromosome_map.tsv as the coords stage would write it, including the gene
    structure column (n_introns)."""
    from gwiscan.coords import MAP_HEADER
    io.write_tsv(cfg.result("chromosome_map.tsv"), MAP_HEADER, [
        ["P_gna", "g1", "GNA", "Chr1", "100", "500", "+", "2"],
        ["P_gna2", "g2", "GNA", "Chr1", "800", "1200", "+", "0"],
        ["P_eul", "g3", "EUL", "Chr2", "50", "400", "-", "1"],
    ])


def _seed_annotations(cfg):
    """Minimal protein-level annotation tables (one row per confirmed protein)."""
    pids = ["P_gna", "P_gna2", "P_eul"]
    io.write_tsv(cfg.protparam_dir / "protparam.tsv",
                 ["protein_id", "length_aa", "molecular_weight"],
                 [[p, "33", "3700.0"] for p in pids])
    io.write_tsv(cfg.result("targetp.tsv"),
                 ["protein_id", "targetp_type", "cs_position"],
                 [[p, "noTP", "-"] for p in pids])
    io.write_tsv(cfg.result("deeptmhmm.tsv"),
                 ["protein_id", "topology", "n_tm_regions"],
                 [[p, "SP", "0"] for p in pids])
    io.write_tsv(cfg.result("deeploc.tsv"),
                 ["protein_id", "localizations"],
                 [[p, "Extracellular"] for p in pids])


def test_merge_confirm_compile_chain(tmp_path):
    cfg = _project(tmp_path)
    _seed_search_hits(cfg)

    # merge: concatenate the two hit sets + extract candidate sequences.
    candidates.run(cfg)
    merged = io.read_tsv(cfg.result("candidates_merged.tsv"))
    assert set(merged["protein_id"]) == {"P_gna", "P_gna2", "P_eul"}
    assert (cfg.result("candidates.fasta")).exists()

    # confirm: InterProScan gate. GNA proteins need PF01453; EUL (no Pfam) is kept.
    _seed_interproscan(cfg)
    confirm.run(cfg)
    final = io.read_tsv(cfg.result("final_candidates.tsv"))
    assert set(zip(final["protein_id"], final["family"])) == {
        ("P_gna", "GNA"), ("P_gna2", "GNA"), ("P_eul", "EUL")}

    # compile: join every annotation into the final table.
    _seed_coords(cfg)
    _seed_annotations(cfg)
    compile_stage.run(cfg)

    out = cfg.final_dir / "gwiscan_results.tsv"
    assert out.exists()
    assert (cfg.final_dir / "gwiscan_results.xlsx").exists()

    df = io.read_tsv(out)
    # One row per (protein, family), no domain-hit duplication.
    assert len(df) == 3
    for col in ("protein_id", "family", "evidence_level", "evidence_support",
                "domain_architecture", "localizations"):
        assert col in df.columns, col

    support = dict(zip(df["protein_id"], df["evidence_support"]))
    assert support["P_gna"] == "both"          # hmm + blast
    assert support["P_gna2"] == "hmm_only"
    assert support["P_eul"] == "blast_only"

    # Gene structure from the annotation flows into the final table (introns only).
    assert "intron_count" in df.columns and "n_exons" not in df.columns
    introns = dict(zip(df["protein_id"], df["intron_count"]))
    assert int(introns["P_gna"]) == 2 and int(introns["P_gna2"]) == 0

    # The transient hit-level columns are dropped from the final table.
    for dropped in ("evalue", "bitscore", "accession", "method", "start", "end"):
        assert dropped not in df.columns

    # GNA proteins carry the Pfam domain architecture; EUL (no InterProScan) is '-'.
    arch = dict(zip(df["protein_id"], df["domain_architecture"]))
    assert "Bulb-type lectin" in arch["P_gna"]
    assert arch["P_eul"] == "-"
