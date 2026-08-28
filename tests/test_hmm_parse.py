#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# test_hmm_parse.py - hmmsearch --domtblout parsing tests.                                          #
#                                                                                                  #
# Pins the column orientation: in hmmsearch output the target is the protein and the query is the   #
# HMM model. Keeping these mapped correctly keeps protein_id and family straight and the candidate  #
# set populated downstream.                                                                        #
#                                                                                                  #
####################################################################################################
"""

from gwiscan import hmm, io
from gwiscan.config import Config
from gwiscan.schema import HIT_HEADER

FAMILY_MAP = (
    "Family\tPfamModel\tBlastModel\tHmmPress\n"
    "Legume\tPF00139\tAAA33983.1.fasta\n"
    "GNA\tPF01453\tAAA33346.1.fasta\n"
)

# Two realistic hmmsearch domtbl data lines (23 fields) plus comment lines that must
# be ignored. Field order: tname tacc tlen qname qacc qlen  Eval score bias # of
# cEval iEval domSc bias hmmF hmmT aliF aliT envF envT acc desc...
# hmmsearch: tname is the PROTEIN, qname/qacc are the HMM name/accession.
DOMTBL = (
    "# target name accession tlen query name ...\n"
    "#------------------- ---------- ----- \n"
    "sp|P12345|LEC_SOYBN - 285 Lectin_legB PF00139.27 250 1.2e-80 270.5 0.1 1 1 "
    "3.4e-84 6.8e-80 268.9 0.1 1 249 34 280 33 282 0.98 Legume lectin domain\n"
    "GNA_prot_007 - 160 B_lectin PF01453.31 105 2.1e-40 130.2 0.0 1 1 "
    "5.0e-44 9.9e-40 128.0 0.0 1 104 12 118 10 120 0.97 D-mannose binding lectin\n"
)


def _make_project(tmp_path):
    (tmp_path / "config").mkdir()
    (tmp_path / "config" / "family.tsv").write_text(FAMILY_MAP)
    cfg = Config(root=tmp_path)
    # cfg.result() routes the hmm working files into intermediate/hmm/ and creates it.
    cfg.result("hmm_domtbl.txt").write_text(DOMTBL)
    return cfg


def _read_hits(cfg):
    # Files are camelCase on disk; io.read_tsv normalises back to snake_case.
    return io.read_tsv(cfg.result("hmm_hits.tsv")).astype(str).to_dict("records")


def test_parse_produces_shared_schema(tmp_path):
    cfg = _make_project(tmp_path)
    n = hmm.parse_domtbl(cfg)
    assert n == 2
    rows = _read_hits(cfg)
    assert list(rows[0].keys()) == HIT_HEADER


def test_protein_and_family_not_swapped(tmp_path):
    cfg = _make_project(tmp_path)
    hmm.parse_domtbl(cfg)
    by_protein = {r["protein_id"]: r for r in _read_hits(cfg)}

    # The target column is the protein id, not the HMM name.
    assert "sp|P12345|LEC_SOYBN" in by_protein
    assert "GNA_prot_007" in by_protein
    assert "Lectin_legB" not in by_protein  # that's the HMM name, must not be an id

    row = by_protein["sp|P12345|LEC_SOYBN"]
    assert row["family"] == "Legume"        # mapped from PF00139 via family_map
    assert row["accession"] == "PF00139.27"
    assert row["evalue"] == "6.8e-80"  # i-Evalue, not the full-seq E-value
    assert row["bitscore"] == "268.9"
    assert row["start"] == "34"           # ali coords on the protein (target)
    assert row["end"] == "280"
    assert row["method"] == "hmm"


def test_unmapped_accession_falls_back_to_hmm_name(tmp_path):
    cfg = _make_project(tmp_path)
    # A hit whose accession is absent from family_map should keep the raw name
    # rather than being dropped.
    extra = DOMTBL + (
        "orphan_prot - 100 Mystery_dom PF99999.1 90 1e-10 40.0 0.0 1 1 "
        "1e-12 1e-10 38.0 0.0 1 89 5 88 3 90 0.9 unknown\n"
    )
    cfg.result("hmm_domtbl.txt").write_text(extra)
    hmm.parse_domtbl(cfg)
    by_protein = {r["protein_id"]: r for r in _read_hits(cfg)}
    assert by_protein["orphan_prot"]["family"] == "Mystery_dom"


# --- GA (gathering) threshold detection for custom identifying HMMs -----------
#
# hmmsearch --cut_ga applies each model's own GA cutoff and aborts the whole search
# if any pressed model lacks one. A profile from hmmbuild has no GA line unless its
# source alignment carried one, so a user's custom identifying HMM must be checked
# before it is pressed. has_ga_thresholds() is what preflight and setup-db use.

_HMM_WITH_GA = (
    "HMMER3/f [3.4 | Aug 2023]\n"
    "NAME  CRA\n"
    "LENG  100\n"
    "GA    25.00 25.00;\n"
    "HMM          A        C\n"
    "//\n"
)
_HMM_NO_GA = (
    "HMMER3/f [3.4 | Aug 2023]\n"
    "NAME  CRA\n"
    "LENG  100\n"
    "HMM          A        C\n"
    "//\n"
)


def test_has_ga_thresholds_true_when_present(tmp_path):
    p = tmp_path / "cra.hmm"
    p.write_text(_HMM_WITH_GA)
    assert hmm.has_ga_thresholds(p) is True


def test_has_ga_thresholds_false_when_absent(tmp_path):
    p = tmp_path / "cra.hmm"
    p.write_text(_HMM_NO_GA)
    assert hmm.has_ga_thresholds(p) is False


def test_has_ga_thresholds_requires_every_model(tmp_path):
    # Two concatenated models, only the first declaring GA: hmmsearch --cut_ga would
    # still abort, so this must read as missing.
    p = tmp_path / "two.hmm"
    p.write_text(_HMM_WITH_GA + _HMM_NO_GA)
    assert hmm.has_ga_thresholds(p) is False


def test_has_ga_thresholds_ignores_other_ga_prefixed_lines(tmp_path):
    # A line like "GATHER..." is not a GA threshold line ("GA" + whitespace).
    p = tmp_path / "x.hmm"
    p.write_text(_HMM_NO_GA.replace("NAME  CRA\n", "NAME  CRA\nGATHERING not a cutoff\n"))
    assert hmm.has_ga_thresholds(p) is False


def test_custom_hmm_ga_error_names_family_and_gives_guidance(tmp_path):
    p = tmp_path / "cra.hmm"
    p.write_text(_HMM_NO_GA)
    msg = hmm.custom_hmm_ga_error(p, "CRA")
    assert msg is not None
    assert "CRA" in msg and "--cut_ga" in msg


def test_custom_hmm_ga_error_none_when_present(tmp_path):
    p = tmp_path / "cra.hmm"
    p.write_text(_HMM_WITH_GA)
    assert hmm.custom_hmm_ga_error(p, "CRA") is None
