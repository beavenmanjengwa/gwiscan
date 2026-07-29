"""Tests for the per-candidate evidence mark (compile stage).

The level carries one '+' per method supporting the call, so the tests pin the
mark against the criteria naming it: the two must always agree.
"""

import pandas as pd

from gwiscan.compile import confidence, evidence_support, support_databases


def test_evidence_support_per_protein_and_family():
    merged = pd.DataFrame([
        {"protein_id": "P1", "family": "GNA", "method": "hmm"},
        {"protein_id": "P1", "family": "GNA", "method": "blast"},
        {"protein_id": "P2", "family": "GNA", "method": "hmm"},
        {"protein_id": "P3", "family": "EUL", "method": "blast"},
        # the same protein can be recovered differently in a second family
        {"protein_id": "P1", "family": "EUL", "method": "blast"},
    ])
    support = evidence_support(merged)
    assert support[("P1", "GNA")] == "both"
    assert support[("P2", "GNA")] == "hmm_only"
    assert support[("P3", "EUL")] == "blast_only"
    assert support[("P1", "EUL")] == "blast_only"


def test_support_databases_counts_distinct_member_dbs():
    ipr = pd.DataFrame([
        {"protein_id": "P1", "analysis": "Pfam"},
        {"protein_id": "P1", "analysis": "CDD"},
        {"protein_id": "P1", "analysis": "Pfam"},      # same db twice counts once
        {"protein_id": "P2", "analysis": "Pfam"},
    ])
    counts = support_databases(ipr)
    assert counts["P1"] == 2
    assert counts["P2"] == 1


def test_all_three_methods_give_three_plus():
    crit, level = confidence("both", 2)
    assert crit == "hmmProfile+blastHit+interproDomain"
    assert level == "+++"


def test_hmm_only_with_interpro_domain_is_two_plus():
    # the GNA case: BLAST is blind, but the profile hit is independently annotated
    crit, level = confidence("hmm_only", 1)
    assert crit == "hmmProfile+interproDomain"
    assert level == "++"


def test_blast_only_without_interpro_is_one_plus():
    crit, level = confidence("blast_only", 0)
    assert crit == "blastHit"
    assert level == "+"


def test_no_supporting_method_gives_no_mark():
    crit, level = confidence("-", 0)
    assert crit == "-"
    assert level == "-"


def test_the_mark_always_counts_the_criteria_named_beside_it():
    for support in ("both", "hmm_only", "blast_only", "-"):
        for n_databases in (0, 1, 3):
            crit, level = confidence(support, n_databases)
            named = 0 if crit == "-" else len(crit.split("+"))
            assert len(level.replace("-", "")) == named


def test_empty_inputs_are_handled():
    assert evidence_support(pd.DataFrame()) == {}
    assert support_databases(None) == {}
