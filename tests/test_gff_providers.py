"""Tests for provider-specific GFF feature extraction."""

from pathlib import Path

from acidgenomes.gff._parser import parse_gff

FIXTURES = Path(__file__).parent / "fixtures"


def test_ensembl_extract_genes():
    from acidgenomes.gff._ensembl import extract_genes

    df = parse_gff(FIXTURES / "ensembl_genes.gtf")
    genes = extract_genes(df, "GTF")
    assert len(genes) > 0
    assert "gene_id" in genes.columns
    assert "gene_id_version" in genes.columns
    # All rows should be genes.
    assert set(genes["type"].unique()) == {"gene"}


def test_ensembl_extract_transcripts():
    from acidgenomes.gff._ensembl import extract_transcripts

    df = parse_gff(FIXTURES / "ensembl_genes.gtf")
    txs = extract_transcripts(df, "GTF")
    assert len(txs) > 0
    assert "transcript_id" in txs.columns
    assert "gene_id" in txs.columns


def test_ensembl_extract_exons():
    from acidgenomes.gff._ensembl import extract_exons

    df = parse_gff(FIXTURES / "ensembl_genes.gtf")
    exons = extract_exons(df, "GTF")
    assert len(exons) > 0
    assert "exon_id" in exons.columns


def test_gencode_extract_genes():
    from acidgenomes.gff._gencode import extract_genes

    df = parse_gff(FIXTURES / "gencode_genes.gtf")
    genes = extract_genes(df, "GTF")
    assert len(genes) > 0
    assert "gene_id" in genes.columns
    assert "gene_id_version" in genes.columns
    # Unversioned IDs should not have a dot.
    assert all("." not in g for g in genes["gene_id"].dropna())
    # Versioned IDs should have a dot.
    assert all("." in g for g in genes["gene_id_version"].dropna())


def test_flybase_extract_genes():
    from acidgenomes.gff._flybase import extract_genes

    df = parse_gff(FIXTURES / "flybase_genes.gtf")
    genes = extract_genes(df, "GTF")
    assert len(genes) > 0
    assert all(g.startswith("FBgn") for g in genes["gene_id"].dropna())


def test_flybase_extract_transcripts():
    from acidgenomes.gff._flybase import extract_transcripts

    df = parse_gff(FIXTURES / "flybase_genes.gtf")
    txs = extract_transcripts(df, "GTF")
    assert len(txs) > 0
    assert "transcript_id" in txs.columns


def test_wormbase_extract_genes():
    from acidgenomes.gff._wormbase import extract_genes

    df = parse_gff(FIXTURES / "wormbase_genes.gtf")
    genes = extract_genes(df, "GTF")
    assert len(genes) > 0
    assert all(g.startswith("WBGene") for g in genes["gene_id"].dropna())


def test_wormbase_gene_id_validation():
    import pandas as pd

    from acidgenomes.gff._wormbase import extract_genes

    bad_df = pd.DataFrame(
        {
            "type": ["gene", "gene"],
            "gene_id": ["WBGene00000001", "BADID12345"],
            "seqnames": ["I", "I"],
            "start": [1, 100],
            "end": [50, 200],
            "strand": ["+", "+"],
        }
    )
    import warnings

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        genes = extract_genes(bad_df, "GTF")
    assert len(genes) == 1
    assert genes.iloc[0]["gene_id"] == "WBGene00000001"
    assert any("malformed" in str(warning.message).lower() for warning in w)
