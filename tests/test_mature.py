"""Tests for the mature-sequence track (TargetP -mature -> per-family FASTAs).

This is an independent track from the domain track: full mature sequences
(presequence cleaved) grouped per family, for a mature-sequence tree.
"""

import pandas as pd

from gwiscan.config import Config
from gwiscan.mature import find_mature_fasta, group_by_family
from gwiscan.msa import _aligned_name


def test_group_by_family_unique_and_ordered():
    final_df = pd.DataFrame([
        {"protein_id": "P1", "family": "GNA"},
        {"protein_id": "P1", "family": "GNA"},    # duplicate row (tandem domains)
        {"protein_id": "P2", "family": "GNA"},
        {"protein_id": "P3", "family": "Legume"},
        {"protein_id": "P1", "family": "Legume"},  # chimera: also in Legume
    ])
    grouped = group_by_family(final_df)
    assert grouped["GNA"] == ["P1", "P2"]        # deduped, order preserved
    assert grouped["Legume"] == ["P3", "P1"]     # a protein can be in two families


def test_find_mature_fasta(tmp_path):
    cfg = Config(root=tmp_path)
    (tmp_path / "intermediate").mkdir()
    assert find_mature_fasta(cfg) is None          # no targetp_raw yet

    raw = tmp_path / "intermediate" / "targetp_raw"
    raw.mkdir()
    assert find_mature_fasta(cfg) is None          # dir present, no mature file

    (raw / "targetp_mature.fasta").write_text(">P1\nMKT\n")
    assert find_mature_fasta(cfg).name == "targetp_mature.fasta"


def test_aligned_names_keep_tracks_separate(tmp_path):
    # domain track drops the _domains tag; mature track keeps _mature in the name,
    # so the two alignments (and their trees) never collide.
    assert _aligned_name(tmp_path / "GNA_domains.fasta") == "GNA_aligned.fasta"
    assert _aligned_name(tmp_path / "GNA_mature.fasta") == "GNA_mature_aligned.fasta"
