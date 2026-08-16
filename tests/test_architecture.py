#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# test_architecture.py - MODE: architecture domain-combination identification tests.               #
#                                                                                                  #
# Both the primary and the required domains are Pfam HMMs. The search is a two-pass hmmscan        #
# (primary genome-wide, then required on the candidates); these tests pin the rules grammar, the   #
# pressed-database accession sets, and the classify step that keeps a protein only when it carries #
# the primary and every required domain.                                                           #
#                                                                                                  #
####################################################################################################
"""

import pytest

from gwiscan import architecture as A
from gwiscan.config import Config
from gwiscan.schema import HIT_HEADER

ARCH_TABLE = (
    "Architecture\tPrimary\tRequired\tClass\n"
    "# a comment line that must be ignored\n"
    "G-LecRLK\tPF01453\tPF00069|PF07714\tLecRLK\n"
    "L-LecRLK\tPF00139\tPF00069|PF07714\tLecRLK\n"
)


def _rules(tmp_path):
    (tmp_path / "config").mkdir(exist_ok=True)
    path = tmp_path / "config" / "architecture.tsv"
    path.write_text(ARCH_TABLE)
    return A.read_rules(path)


# candidate-pass hits: dicts as produced by A._parse_domtbl
def _hit(protein, name, pfam, start, end, ev=1e-40, bs=100.0):
    return {"protein": protein, "name": name, "pfam": pfam, "acc": pfam + ".1",
            "evalue": str(ev), "bitscore": str(bs), "start": start, "end": end}


def test_read_rules_grammar(tmp_path):
    rules = _rules(tmp_path)
    g = next(r for r in rules if r["architecture"] == "G-LecRLK")
    assert g["primary"] == ["PF01453"]
    assert g["required"] == [["PF00069", "PF07714"]]   # one required slot, OR inside
    assert g["class"] == "LecRLK"


def test_primary_and_all_accessions(tmp_path):
    rules = _rules(tmp_path)
    assert A.primary_accessions(rules) == ["PF01453", "PF00139"]
    # pass-2 db is primary + required, de-duplicated in first-seen order
    assert A.all_accessions(rules) == ["PF01453", "PF00069", "PF07714", "PF00139"]


def test_primary_required_must_be_pfam(tmp_path):
    (tmp_path / "config").mkdir(exist_ok=True)
    p = tmp_path / "config" / "architecture.tsv"
    p.write_text("Architecture\tPrimary\tRequired\tClass\nBad\tB_lectin\tPF00069\tX\n")
    with pytest.raises(ValueError):
        A.read_rules(p)


def test_missing_required_rejected(tmp_path):
    (tmp_path / "config").mkdir(exist_ok=True)
    p = tmp_path / "config" / "architecture.tsv"
    p.write_text("Architecture\tPrimary\tRequired\tClass\nBad\tPF01453\t-\tX\n")
    with pytest.raises(ValueError):
        A.read_rules(p)


def test_classify_keeps_only_full_combinations(tmp_path):
    rules = _rules(tmp_path)
    hits = [
        # p_g: B_lectin (primary) + Pkinase (required) -> G-LecRLK
        _hit("p_g", "B_lectin", "PF01453", 10, 120),
        _hit("p_g", "Pkinase", "PF00069", 200, 460),
        # p_l: Lectin_legB + Pkinase_Tyr -> L-LecRLK (OR alternative kinase)
        _hit("p_l", "Lectin_legB", "PF00139", 20, 250),
        _hit("p_l", "Pkinase_Tyr", "PF07714", 300, 560),
        # p_no_kinase: lectin only -> dropped (candidate from pass 1, fails required)
        _hit("p_no_kinase", "B_lectin", "PF01453", 10, 120),
    ]
    rows = A.classify(hits, rules)
    by = {(r[0], r[1]) for r in rows}
    assert ("p_g", "G-LecRLK") in by
    assert ("p_l", "L-LecRLK") in by
    assert not any(pid == "p_no_kinase" for pid, _ in by)      # no required kinase
    assert ("p_g", "L-LecRLK") not in by                       # wrong lectin
    # rows carry the standard hit schema, family = the architecture
    assert all(len(r) == len(HIT_HEADER) for r in rows)
    assert all(r[7] == "hmm" for r in rows)


def test_classify_emits_primary_and_required_domains(tmp_path):
    rules = _rules(tmp_path)
    hits = [
        _hit("p_g", "B_lectin", "PF01453", 10, 120),
        _hit("p_g", "Pkinase", "PF00069", 200, 460),
    ]
    rows = [r for r in A.classify(hits, rules) if r[1] == "G-LecRLK"]
    pfams = sorted(r[2].split(".")[0] for r in rows)
    assert pfams == ["PF00069", "PF01453"]     # both the primary and the required domain


def test_config_mode_and_paths(tmp_path):
    cfg = Config(root=tmp_path, MODE="architecture")
    assert cfg.is_architecture
    assert cfg.architecture_map == tmp_path / "config" / "architecture.tsv"
    assert cfg.primary_hmm_db == tmp_path / "db" / "hmm" / "primary_models.hmm"
    assert cfg.hmm_db == tmp_path / "db" / "hmm" / "all_models.hmm"
