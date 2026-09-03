#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# test_score.py - Per-family detectability-profile tests (score stage).                            #
#                                                                                                  #
# The cases mirror the two extremes: a conserved family both searches recover, and a divergent     #
# family only the profile HMM recovers.                                                            #
#                                                                                                  #
####################################################################################################
"""

from gwiscan.config import Config
from gwiscan.score import HEADER, blast_hits_by_family, family_metrics, verdict

MIN = 0.7


def test_concordant_family():
    # both searches recover the same members (conserved family, e.g. Hevein)
    m = family_metrics({"P1", "P2", "P3"}, {"P1", "P2", "P3"}, has_hmm=True, concordance_min=MIN)
    assert m["jaccard"] == 1.0
    assert m["hmm_only_frac"] == 0.0 and m["blast_only_frac"] == 0.0
    assert m["blast_recall_of_hmm"] == 1.0
    assert m["detectability"] == "concordant"


def test_hmm_dominant_family_blast_finds_nothing():
    # divergent family (e.g. GNA/bulb-type): hmmsearch hits, DIAMOND empty
    m = family_metrics({"P1", "P2", "P3", "P4"}, set(), has_hmm=True, concordance_min=MIN)
    assert m["n_blast"] == 0 and m["jaccard"] == 0.0
    assert m["hmm_only_frac"] == 1.0
    assert m["blast_recall_of_hmm"] == 0.0
    assert m["detectability"] == "hmm_dominant"
    assert "under-report" in verdict("GNA", m)


def test_partial_overlap_below_threshold_is_hmm_dominant():
    # 8 HMM hits, 2 shared, 1 BLAST-only -> jaccard 2/9, HMM side larger
    m = family_metrics({f"P{i}" for i in range(8)}, {"P0", "P1", "X"},
                       has_hmm=True, concordance_min=MIN)
    assert m["n_union"] == 9 and m["n_both"] == 2
    assert m["hmm_only_frac"] > m["blast_only_frac"]
    assert m["detectability"] == "hmm_dominant"


def test_blast_dominant_when_profile_misses_members():
    m = family_metrics({"P1"}, {"P1", "P2", "P3", "P4"}, has_hmm=True, concordance_min=MIN)
    assert m["blast_only_frac"] > m["hmm_only_frac"]
    assert m["detectability"] == "blast_dominant"
    assert "3 member(s) found by BLAST alone" in verdict("Fam", m)


def test_family_without_hmm_is_stated_not_judged():
    # no HMM -> no comparison is possible; the verdict points at the custom-HMM route
    m = family_metrics(set(), {"P1", "P2"}, has_hmm=False, concordance_min=MIN)
    assert m["detectability"] == "blast_only_family"
    assert "curated family HMM" in verdict("EUL", m)


def test_family_found_by_neither_search():
    m = family_metrics(set(), set(), has_hmm=True, concordance_min=MIN)
    assert m["n_union"] == 0 and m["jaccard"] == 0.0
    assert m["detectability"] == "not_detected"


def test_metrics_cover_every_written_column():
    m = family_metrics({"P1"}, {"P1"}, has_hmm=True, concordance_min=MIN)
    # family/has_hmm come from the family table; the rest must come from the metrics
    assert set(HEADER[2:]) <= set(m)


def test_blast_hits_by_family_groups_members(tmp_path):
    # concordance compares hmmsearch vs ALL DIAMOND BLAST members (blast_hits.tsv),
    # grouped per family.
    cfg = Config(root=tmp_path)
    cfg.result("blast_hits.tsv").write_text(
        "proteinId\tfamily\taccession\tevalue\tbitscore\tstart\tend\tmethod\n"
        "P1\tGNA\t-\t1e-40\t150\t1\t100\tblast\n"
        "P2\tGNA\t-\t1e-20\t100\t1\t100\tblast\n"
        "P3\tCRA\t-\t1e-10\t80\t1\t100\tblast\n"
    )
    assert blast_hits_by_family(cfg) == {"GNA": {"P1", "P2"}, "CRA": {"P3"}}
    assert blast_hits_by_family(Config(root=tmp_path / "empty")) == {}
