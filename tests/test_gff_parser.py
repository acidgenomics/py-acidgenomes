"""Tests for the core GFF/GTF streaming parser."""

from pathlib import Path

import pandas as pd

from acidgenomes.gff._parser import parse_gff, read_gff_directives

FIXTURES = Path(__file__).parent / "fixtures"


def test_parse_gtf_returns_dataframe():
    df = parse_gff(FIXTURES / "ensembl_genes.gtf")
    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0


def test_parse_gtf_standard_columns():
    df = parse_gff(FIXTURES / "ensembl_genes.gtf")
    for col in ("seqnames", "source", "type", "start", "end", "strand"):
        assert col in df.columns, f"Missing column: {col}"


def test_parse_gtf_gene_attributes_present():
    df = parse_gff(FIXTURES / "ensembl_genes.gtf")
    assert "gene_id" in df.columns
    assert "gene_name" in df.columns


def test_parse_gtf_feature_filter():
    df = parse_gff(FIXTURES / "ensembl_genes.gtf", feature_filter={"gene"})
    assert set(df["type"].unique()) == {"gene"}


def test_parse_gtf_coordinates_are_int():
    df = parse_gff(FIXTURES / "ensembl_genes.gtf")
    assert df["start"].dtype.kind == "i"
    assert df["end"].dtype.kind == "i"


def test_read_directives_ensembl():
    directives = read_gff_directives(FIXTURES / "ensembl_genes.gtf")
    assert isinstance(directives, dict)
    assert "genome-build" in directives
    assert directives["genome-build"] == "GRCh38"


def test_read_directives_wormbase():
    directives = read_gff_directives(FIXTURES / "wormbase_genes.gtf")
    assert "genebuild-version" in directives
    assert directives["genebuild-version"] == "WS289"


def test_parse_flybase_gtf():
    df = parse_gff(FIXTURES / "flybase_genes.gtf")
    assert "gene_id" in df.columns
    fbgn = df[df["type"] == "gene"]["gene_id"]
    assert all(g.startswith("FBgn") for g in fbgn.dropna())


def test_parse_wormbase_gtf():
    df = parse_gff(FIXTURES / "wormbase_genes.gtf")
    wb = df[df["type"] == "gene"]["gene_id"]
    assert all(g.startswith("WBGene") for g in wb.dropna())
