"""Tests for hmmscan --domtblout parsing.

Pins the column orientation: in hmmscan output the *target* is the HMM model and
the *query* is the protein. Keeping these mapped correctly is what keeps
protein_id and family straight and the candidate set populated downstream.
"""

from gwiscan import hmm, io
from gwiscan.config import Config
from gwiscan.schema import HIT_HEADER

FAMILY_MAP = (
    "Family\tPfamModel\tBlastModel\tHmmPress\n"
    "Legume\tPF00139\tAAA33983.1.fasta\n"
    "GNA\tPF01453\tAAA33346.1.fasta\n"
)

# Two realistic domtbl data lines (23 fields) plus comment lines that must be
# ignored. Field order: tname tacc tlen qname qacc qlen  Eval score bias # of
# cEval iEval domSc bias hmmF hmmT aliF aliT envF envT acc desc...
DOMTBL = (
    "# target name accession tlen query name ...\n"
    "#------------------- ---------- ----- \n"
    "Lectin_legB PF00139.27 250 sp|P12345|LEC_SOYBN - 285 1.2e-80 270.5 0.1 1 1 "
    "3.4e-84 6.8e-80 268.9 0.1 1 249 34 280 33 282 0.98 Legume lectin domain\n"
    "B_lectin PF01453.31 105 GNA_prot_007 - 160 2.1e-40 130.2 0.0 1 1 "
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

    # The query column is the protein id, not the HMM name.
    assert "sp|P12345|LEC_SOYBN" in by_protein
    assert "GNA_prot_007" in by_protein
    assert "Lectin_legB" not in by_protein  # that's the HMM name, must not be an id

    row = by_protein["sp|P12345|LEC_SOYBN"]
    assert row["family"] == "Legume"        # mapped from PF00139 via family_map
    assert row["accession"] == "PF00139.27"
    assert row["evalue"] == "6.8e-80"  # i-Evalue, not the full-seq E-value
    assert row["bitscore"] == "268.9"
    assert row["start"] == "34"           # ali coords on the query
    assert row["end"] == "280"
    assert row["method"] == "hmm"


def test_unmapped_accession_falls_back_to_hmm_name(tmp_path):
    cfg = _make_project(tmp_path)
    # A hit whose accession is absent from family_map should keep the raw name
    # rather than being dropped.
    extra = DOMTBL + (
        "Mystery_dom PF99999.1 90 orphan_prot - 100 1e-10 40.0 0.0 1 1 "
        "1e-12 1e-10 38.0 0.0 1 89 5 88 3 90 0.9 unknown\n"
    )
    cfg.result("hmm_domtbl.txt").write_text(extra)
    hmm.parse_domtbl(cfg)
    by_protein = {r["protein_id"]: r for r in _read_hits(cfg)}
    assert by_protein["orphan_prot"]["family"] == "Mystery_dom"
