"""Tests for GO id -> name / namespace mapping from go-basic.obo.

InterProScan reports GO ids only; these map them to names and split by aspect.
"""

from gwiscan.go import clean_ids, names_by_aspect, names_for, parse_obo

# Includes an obsolete term and a [Typedef] stanza, as in the real file.
OBO = """format-version: 1.2
data-version: releases/2026-01-01

[Term]
id: GO:0004672
name: protein kinase activity
namespace: molecular_function
def: "Catalysis of the transfer of a phosphate group." [GOC:x]

[Term]
id: GO:0005537
name: mannose binding
namespace: molecular_function

[Term]
id: GO:0016020
name: membrane
namespace: cellular_component

[Term]
id: GO:0006468
name: protein phosphorylation
namespace: biological_process

[Term]
id: GO:0000005
name: obsolete ribosomal chaperone activity
namespace: molecular_function
is_obsolete: true
consider: GO:0042254

[Typedef]
id: part_of
name: part of
"""


def _obo(tmp_path):
    p = tmp_path / "go-basic.obo"
    p.write_text(OBO, encoding="utf-8")
    return p


def test_parse_captures_id_name_namespace(tmp_path):
    terms = parse_obo(_obo(tmp_path))
    assert terms["GO:0004672"] == {"name": "protein kinase activity",
                                   "namespace": "molecular_function"}
    assert terms["GO:0016020"]["namespace"] == "cellular_component"
    assert terms["GO:0006468"]["namespace"] == "biological_process"
    # obsolete terms still map (they keep their GO name)
    assert terms["GO:0000005"]["name"].startswith("obsolete")
    # [Typedef] ids are not GO ids and must not leak in
    assert "part_of" not in terms
    assert len(terms) == 5


def test_names_for_preserves_order_and_unknowns(tmp_path):
    terms = parse_obo(_obo(tmp_path))
    assert names_for("GO:0004672|GO:0005537", terms) == "protein kinase activity|mannose binding"
    assert names_for("GO:0005537|GO:9999999", terms) == "mannose binding|GO:9999999"
    assert names_for("-", terms) == "-"
    assert names_for("GO:0004672", {}) == "-"


def test_names_by_aspect_splits_three_ways(tmp_path):
    terms = parse_obo(_obo(tmp_path))
    got = names_by_aspect("GO:0004672|GO:0006468|GO:0016020|GO:0005537", terms)
    assert got["molecular_function"] == "protein kinase activity|mannose binding"
    assert got["biological_process"] == "protein phosphorylation"
    assert got["cellular_component"] == "membrane"


def test_names_by_aspect_empty_aspects_dashed(tmp_path):
    terms = parse_obo(_obo(tmp_path))
    got = names_by_aspect("GO:0005537", terms)
    assert got["molecular_function"] == "mannose binding"
    assert got["biological_process"] == "-"
    assert got["cellular_component"] == "-"


def test_clean_ids_strips_interproscan_source_annotation():
    # InterProScan's go_terms field appends the source db to each id, e.g.
    # "GO:0005515(InterPro)"; the bare id is what must be looked up / reported.
    assert clean_ids("GO:0005515(InterPro)|GO:0004553(PANTHER)") == [
        "GO:0005515", "GO:0004553",
    ]
    assert clean_ids("GO:0005537") == ["GO:0005537"]  # already bare, unaffected


def test_names_for_resolves_ids_with_source_annotation(tmp_path):
    terms = parse_obo(_obo(tmp_path))
    assert (names_for("GO:0004672(InterPro)|GO:0005537(InterPro)", terms)
            == "protein kinase activity|mannose binding")
