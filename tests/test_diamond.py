#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# test_diamond.py - DIAMOND search tests: round options, seed selection, best-per-member.          #
#                                                                                                  #
# Both rounds search the proteome (hit proteins are the subject, sseqid). Round 1 is --very-       #
# sensitive, E-value only. Round-2 seeds are chosen by Blast Score Ratio for every family; DIAMOND #
# is self-contained and never reads hmmsearch results.                                            #
#                                                                                                  #
####################################################################################################
"""

from gwiscan import diamond
from gwiscan.config import Config
from gwiscan.diamond import (
    _best_per_member,
    _bsr_seeds,
)

# outfmt6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
R2 = (
    "Prot_A\tProt_A\t100\t200\t0\t0\t1\t200\t1\t200\t0.0\t400\n"
    "Prot_A\tProt_C\t45\t100\t50\t2\t1\t100\t5\t104\t1e-20\t150\n"
    "Prot_B\tProt_C\t60\t120\t30\t1\t1\t120\t3\t122\t1e-30\t200\n"
    "Prot_B\tProt_B\t100\t180\t0\t0\t1\t180\t1\t180\t0.0\t360\n"
)


def _r1(subjects_and_scores):
    """One round-1 line per (subject, bitscore): 'model' hits each subject."""
    return "".join(
        f"model\t{sub}\t60\t200\t0\t0\t1\t200\t1\t200\t1e-30\t{bs}\n"
        for sub, bs in subjects_and_scores
    )


def test_bsr_seeds_by_ratio(tmp_path):
    r1 = tmp_path / "r1.tsv"
    r1.write_text(_r1([("Prot_A", 400), ("Prot_B", 250), ("Prot_C", 150)]))
    # model self-bitscore 500 -> BSR = 0.80, 0.50, 0.30.
    seeds = _bsr_seeds(r1, {"model": 500.0}, threshold=0.4)
    assert seeds == ["Prot_A", "Prot_B"]      # 0.30 is below 0.4


def test_both_rounds_evalue_only(tmp_path, monkeypatch):
    calls = []

    def fake_run(cmd, **kwargs):
        cmd = [str(c) for c in cmd]
        calls.append(cmd)
        open(cmd[cmd.index("--out") + 1], "w").close()

    monkeypatch.setattr(diamond.external, "run", fake_run)
    monkeypatch.setattr(diamond.external, "require", lambda b: None)
    cfg = Config(root=tmp_path)

    # Both rounds are called E-value + sensitivity only -- no identity/coverage.
    diamond._round(cfg, "q", "db", tmp_path / "r1.tsv", sensitivity="--very-sensitive")
    diamond._round(cfg, "q", "db", tmp_path / "r2.tsv",
                   sensitivity=diamond._sensitivity_flag(cfg))

    r1, r2 = calls
    assert "--very-sensitive" in r1
    assert "--id" not in r1 and "--query-cover" not in r1   # round 1: E-value only
    assert "--id" not in r2 and "--query-cover" not in r2   # round 2: E-value only


def test_sensitivity_flag_from_config(tmp_path):
    # default is ultra-sensitive (matches NCBI BLASTP), applied to both rounds
    assert diamond._sensitivity_flag(Config(root=tmp_path)) == "--ultra-sensitive"
    assert diamond._sensitivity_flag(
        Config(root=tmp_path, DIAMOND_SENSITIVITY="very-sensitive")) == "--very-sensitive"
    # a value written with leading dashes is tolerated
    assert diamond._sensitivity_flag(
        Config(root=tmp_path, DIAMOND_SENSITIVITY="--sensitive")) == "--sensitive"
    # fast / empty means DIAMOND's default mode: no flag
    assert diamond._sensitivity_flag(Config(root=tmp_path, DIAMOND_SENSITIVITY="fast")) is None
    assert diamond._sensitivity_flag(Config(root=tmp_path, DIAMOND_SENSITIVITY="")) is None


def test_best_per_member_keeps_highest_bitscore(tmp_path):
    r2 = tmp_path / "r2.tsv"
    r2.write_text(R2)
    best = _best_per_member(r2)
    assert set(best) == {"Prot_A", "Prot_B", "Prot_C"}
    assert best["Prot_A"][1] == 400.0
    assert best["Prot_C"][1] == 200.0          # Prot_C hit by two seeds -> keep 200
    assert best["Prot_C"][2] == "3" and best["Prot_C"][3] == "122"
