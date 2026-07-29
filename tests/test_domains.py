"""Tests for systematic domain-id assignment."""

import pandas as pd

from gwiscan.domains import MAP_HEADER, assign_domain_ids

PFAM_TO_FAMILY = {"PF01453": "GNA", "PF00139": "Legume", "PF00704": "CRA"}

# protein_id, pfam_id, start, end  (as in domains.bed)
BED = pd.DataFrame([
    {"protein_id": "Xp2", "pfam_id": "PF01453", "start": 130, "end": 230},  # GNA, 2nd by coord
    {"protein_id": "Xp2", "pfam_id": "PF01453", "start": 10, "end": 110},   # GNA, 1st by coord (same protein)
    {"protein_id": "Xp1", "pfam_id": "PF01453", "start": 5, "end": 100},    # GNA, different protein
    {"protein_id": "Xa", "pfam_id": "PF00139", "start": 34, "end": 280},    # Legume
    {"protein_id": "Xc", "pfam_id": "PF00704", "start": 25, "end": 300},    # CRA
])


def test_id_format_and_serials():
    recs = assign_domain_ids(BED, PFAM_TO_FAMILY, "Ath")
    ids = {r["domain_id"]: r for r in recs}

    # per-family serial over sorted proteins: Xp1 -> 001, Xp2 -> 002
    assert "AthGNA001.1" in ids                       # Xp1, single domain
    assert ids["AthGNA001.1"]["protein_id"] == "Xp1"

    # tandem domains on Xp2 share serial 002, ordered by start -> .1 then .2
    assert ids["AthGNA002.1"]["start"] == 10
    assert ids["AthGNA002.2"]["start"] == 130
    assert ids["AthGNA002.1"]["protein_id"] == "Xp2"

    # other families numbered independently
    assert "AthLegume001.1" in ids
    assert "AthCRA001.1" in ids
    assert ids["AthCRA001.1"]["pfam_id"] == "PF00704"


def test_no_separators_and_map_fields():
    recs = assign_domain_ids(BED, PFAM_TO_FAMILY, "Ath")
    assert list(recs[0].keys()) == MAP_HEADER
    for r in recs:
        assert "_" not in r["domain_id"]                    # no underscores
        assert r["domain_id"].startswith("Ath")


def test_empty_prefix():
    recs = assign_domain_ids(BED, PFAM_TO_FAMILY, "")
    assert any(r["domain_id"] == "GNA001.1" for r in recs)
