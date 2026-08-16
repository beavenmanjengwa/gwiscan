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

from gwiscan.score import HEADER, family_metrics, round1_subjects, verdict

MIN = 0.7


def test_concordant_family():
    # both searches recover the same members (conserved family, e.g. Hevein)
    m = family_metrics({"P1", "P2", "P3"}, {"P1", "P2", "P3"}, has_hmm=True, concordance_min=MIN)
    assert m["jaccard"] == 1.0
    assert m["hmm_only_frac"] == 0.0 and m["blast_only_frac"] == 0.0
    assert m["blast_recall_of_hmm"] == 1.0
    assert m["detectability"] == "concordant"


def test_hmm_dominant_family_blast_finds_nothing():
    # divergent family (e.g. GNA/bulb-type): hmmscan hits, DIAMOND round 1 empty
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


def test_round1_subjects_reads_the_sseqid_column(tmp_path):
    r1 = tmp_path / "diamond_GNA_r1.tsv"
    r1.write_text("q1\tP1\t80.0\t100\t0\t0\t1\t100\t1\t100\t1e-40\t150\n"
                  "q1\tP2\t60.0\t100\t0\t0\t1\t100\t1\t100\t1e-20\t100\n"
                  "q2\tP1\t55.0\t100\t0\t0\t1\t100\t1\t100\t1e-15\t90\n")
    assert round1_subjects(r1) == {"P1", "P2"}
    assert round1_subjects(tmp_path / "absent.tsv") == set()
