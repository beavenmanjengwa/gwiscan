"""Tests for genomic coordinates of family members (coords stage).

Both annotation flavours are covered, because the same pipeline must accept the
GTF and the GFF3 a genome portal hands out for the same representative gene set.
"""

from gwiscan.coords import (
    locate,
    parse_annotation,
    parse_attributes,
    strip_prefix,
    write_gff3,
)

GTF = """\
#!genome-build test
Chr1\tphytozome\tgene\t1000\t5000\t.\t+\t.\tgene_id "Medtr1g001"; transcript_id "Medtr1g001.1";
Chr1\tphytozome\texon\t1000\t2000\t.\t+\t.\tgene_id "Medtr1g001"; transcript_id "Medtr1g001.1";
Chr1\tphytozome\texon\t4000\t5000\t.\t+\t.\tgene_id "Medtr1g001"; transcript_id "Medtr1g001.1";
Chr3\tphytozome\tgene\t8000\t9500\t.\t-\t.\tgene_id "Medtr3g077"; transcript_id "Medtr3g077.1";
"""

GFF3 = """\
##gff-version 3
Chr1\taraport\tgene\t1000\t5000\t.\t+\t.\tID=gene:AT1G01010;Name=AT1G01010
Chr1\taraport\tmRNA\t1000\t5000\t.\t+\t.\tID=transcript:AT1G01010.1;Parent=gene:AT1G01010
Chr1\taraport\tCDS\t1200\t2000\t.\t+\t0\tID=CDS:AT1G01010.1;Parent=transcript:AT1G01010.1
"""


def test_parse_attributes_handles_both_flavours():
    gtf = parse_attributes('gene_id "Medtr1g001"; transcript_id "Medtr1g001.1";')
    assert gtf["gene_id"] == "Medtr1g001"
    assert gtf["transcript_id"] == "Medtr1g001.1"

    gff = parse_attributes("ID=transcript:AT1G01010.1;Parent=gene:AT1G01010")
    assert gff["ID"] == "transcript:AT1G01010.1"
    assert gff["Parent"] == "gene:AT1G01010"


def test_strip_prefix_only_strips_known_feature_prefixes():
    assert strip_prefix("transcript:AT1G01010.1") == "AT1G01010.1"
    assert strip_prefix("gene:AT1G01010") == "AT1G01010"
    # a chromosome-style id must survive intact
    assert strip_prefix("Chr1:12345") == "Chr1:12345"


def test_gtf_span_covers_every_line_naming_the_feature(tmp_path):
    path = tmp_path / "rep.gtf"
    path.write_text(GTF, encoding="utf-8")
    spans, genes = parse_annotation(path)
    # exons at 1000-2000 and 4000-5000 -> the transcript spans 1000-5000
    assert spans["Medtr1g001.1"] == ("Chr1", 1000, 5000, "+")
    assert genes["Medtr1g001.1"] == "Medtr1g001"      # transcript -> its gene
    assert spans["Medtr1g001"] == ("Chr1", 1000, 5000, "+")
    assert spans["Medtr3g077.1"] == ("Chr3", 8000, 9500, "-")


def test_gff3_ids_are_indexed_without_their_prefixes(tmp_path):
    path = tmp_path / "rep.gff3"
    path.write_text(GFF3, encoding="utf-8")
    spans, genes = parse_annotation(path)
    assert spans["AT1G01010.1"] == ("Chr1", 1000, 5000, "+")
    assert genes["AT1G01010.1"] == "AT1G01010"        # mRNA Parent -> its gene
    assert spans["AT1G01010"] == ("Chr1", 1000, 5000, "+")


def test_locate_falls_back_to_the_gene_id(tmp_path):
    path = tmp_path / "rep.gtf"
    path.write_text(GTF, encoding="utf-8")
    spans, genes = parse_annotation(path)

    # exact transcript id, gene-keyed id via the isoform-suffix fallback, and a miss
    located, unresolved = locate(["Medtr1g001.1", "Medtr3g077.2", "NotInGenome.1"],
                                 spans, genes)
    # a transcript-level match reports the GENE, never the transcript id
    assert located["Medtr1g001.1"] == ("Medtr1g001", "Chr1", 1000, 5000, "+")
    assert located["Medtr3g077.2"] == ("Medtr3g077", "Chr3", 8000, 9500, "-")
    assert unresolved == ["NotInGenome.1"]


def test_unresolved_proteins_are_reported_not_dropped(tmp_path):
    path = tmp_path / "rep.gtf"
    path.write_text(GTF, encoding="utf-8")
    spans, genes = parse_annotation(path)
    located, unresolved = locate(["Medtr1g001.1", "Ghost1", "Ghost2"], spans, genes)
    assert len(located) == 1
    assert sorted(unresolved) == ["Ghost1", "Ghost2"]


def test_write_gff3_emits_valid_feature_lines(tmp_path):
    out = tmp_path / "members.gff3"
    write_gff3(out, [["Medtr1g001.1", "Medtr1g001", "GNA", "Chr1", 1000, 5000, "+"]])
    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "##gff-version 3"
    fields = lines[1].split("\t")
    assert len(fields) == 9
    assert fields[0] == "Chr1" and fields[1] == "GWIscan" and fields[2] == "gene"
    assert fields[3] == "1000" and fields[4] == "5000" and fields[6] == "+"
    assert fields[8] == "ID=Medtr1g001;protein=Medtr1g001.1;family=GNA"


# --- new formats: 4-column files + Phytozome Name= attribute --------------------

from gwiscan import coords as _coords


def test_four_column_annotation(tmp_path):
    # chrom, feature id, start, end  (no strand/attributes) -- e.g. Arabidopsis.gff
    gff = tmp_path / "simple.gff"
    gff.write_text("atChr1\tAT1G01010.1\t3631\t5899\n"
                   "Gm01\tGlyma.01G000100.2\t78503\t103594\n")
    spans, genes = _coords.parse_annotation(gff)
    located, unresolved = _coords.locate(["AT1G01010.1", "Glyma.01G000100.2"], spans, genes)
    assert located["AT1G01010.1"] == ("AT1G01010", "atChr1", 3631, 5899, ".")
    assert located["Glyma.01G000100.2"][1:4] == ("Gm01", 78503, 103594)
    assert unresolved == []


def test_phytozome_name_attribute_is_indexed(tmp_path):
    # Proteome id matches Name=, not ID= (Phytozome GFF3).
    gff = tmp_path / "phyto.gff3"
    gff.write_text(
        "sc9\tphytozome\tmRNA\t1869\t4028\t.\t-\t.\t"
        "ID=Cucsa.000200.1.v1.122;Name=Cucsa.000200.1;Parent=Cucsa.000200.v1\n")
    spans, _ = _coords.parse_annotation(gff)
    located, unresolved = _coords.locate(["Cucsa.000200.1"], spans)
    assert located["Cucsa.000200.1"][1:] == ("sc9", 1869, 4028, "-")


def test_dot_p_suffix_resolves_against_name(tmp_path):
    # Pvulgaris proteome ids carry a trailing ".p" that the GFF Name does not.
    gff = tmp_path / "pv.gff3"
    gff.write_text(
        "Chr01\tphytozome\tmRNA\t1705\t6715\t.\t-\t.\t"
        "ID=Phvul.001G000400.3.v2.1;Name=Phvul.001G000400.3;Parent=Phvul.001G000400.v2.1\n")
    spans, _ = _coords.parse_annotation(gff)
    located, unresolved = _coords.locate(["Phvul.001G000400.3.p"], spans)
    assert "Phvul.001G000400.3.p" in located
    assert located["Phvul.001G000400.3.p"][1:4] == ("Chr01", 1705, 6715)
