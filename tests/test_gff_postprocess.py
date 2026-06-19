"""Tests for GFF post-processing: broadClass, version handling, column normalization."""

import pandas as pd

from acidgenomes.gff._postprocess import (
    add_broad_class,
    drop_all_na_columns,
    handle_versions,
    remove_unwanted_biotypes,
    sort_ranges,
    standardize_column_names,
)


def test_broad_class_coding():
    df = pd.DataFrame(
        {
            "gene_biotype": ["protein_coding"],
            "seqnames": ["1"],
            "gene_name": ["TP53"],
        }
    )
    result = add_broad_class(df)
    assert result["broad_class"].iloc[0] == "coding"


def test_broad_class_mito_by_chromosome():
    df = pd.DataFrame(
        {
            "gene_biotype": ["protein_coding"],
            "seqnames": ["MT"],
            "gene_name": ["MT-ND1"],
        }
    )
    result = add_broad_class(df)
    assert result["broad_class"].iloc[0] == "mito"


def test_broad_class_mito_by_gene_name():
    df = pd.DataFrame(
        {
            "gene_biotype": ["protein_coding"],
            "seqnames": ["1"],
            "gene_name": ["mt-Nd1"],
        }
    )
    result = add_broad_class(df)
    assert result["broad_class"].iloc[0] == "mito"


def test_broad_class_pseudo():
    df = pd.DataFrame(
        {
            "gene_biotype": ["pseudogene"],
            "seqnames": ["1"],
            "gene_name": ["GENE1P"],
        }
    )
    result = add_broad_class(df)
    assert result["broad_class"].iloc[0] == "pseudo"


def test_broad_class_small_rna():
    for biotype in ("miRNA", "snoRNA", "snRNA", "rRNA"):
        df = pd.DataFrame(
            {
                "gene_biotype": [biotype],
                "seqnames": ["1"],
                "gene_name": ["GENE1"],
            }
        )
        result = add_broad_class(df)
        assert result["broad_class"].iloc[0] == "small", f"Expected 'small' for {biotype}"


def test_broad_class_ig():
    df = pd.DataFrame(
        {
            "gene_biotype": ["ig_c_gene"],
            "seqnames": ["14"],
            "gene_name": ["IGHG1"],
        }
    )
    result = add_broad_class(df)
    assert result["broad_class"].iloc[0] == "ig"


def test_remove_unwanted_biotypes_ensembl():
    df = pd.DataFrame(
        {
            "gene_id": ["ENSG1", "ENSG2", "ENSG3"],
            "gene_biotype": ["protein_coding", "TEC", "LRG_gene"],
        }
    )
    result = remove_unwanted_biotypes(df, provider="Ensembl")
    assert len(result) == 1
    assert result.iloc[0]["gene_biotype"] == "protein_coding"


def test_remove_unwanted_biotypes_gencode_keeps_tec():
    df = pd.DataFrame(
        {
            "gene_id": ["ENSG1", "ENSG2"],
            "gene_biotype": ["protein_coding", "TEC"],
        }
    )
    result = remove_unwanted_biotypes(df, provider="GENCODE")
    assert len(result) == 2


def test_handle_versions_ignore_false():
    df = pd.DataFrame(
        {
            "gene_id": ["ENSG00000001"],
            "gene_id_version": ["ENSG00000001.5"],
        }
    )
    result = handle_versions(
        df, id_col="gene_id", version_col="gene_id_version", ignore_version=False
    )
    assert result["gene_id"].iloc[0] == "ENSG00000001.5"
    assert result["gene_id_no_version"].iloc[0] == "ENSG00000001"


def test_handle_versions_ignore_true():
    df = pd.DataFrame(
        {
            "gene_id": ["ENSG00000001"],
            "gene_id_version": ["ENSG00000001.5"],
        }
    )
    result = handle_versions(
        df, id_col="gene_id", version_col="gene_id_version", ignore_version=True
    )
    assert result["gene_id"].iloc[0] == "ENSG00000001"
    assert "gene_id_no_version" not in result.columns


def test_standardize_column_names_gene_type():
    df = pd.DataFrame({"gene_type": ["protein_coding"]})
    result = standardize_column_names(df)
    assert "gene_biotype" in result.columns
    assert "gene_type" not in result.columns


def test_standardize_column_names_symbol():
    df = pd.DataFrame({"symbol": ["TP53"]})
    result = standardize_column_names(df)
    assert "gene_name" in result.columns
    assert "symbol" not in result.columns


def test_sort_ranges():
    df = pd.DataFrame(
        {
            "seqnames": ["2", "1", "1"],
            "start": [100, 300, 50],
            "gene_id": ["g3", "g2", "g1"],
        }
    )
    result = sort_ranges(df)
    assert result.iloc[0]["gene_id"] == "g1"
    assert result.iloc[1]["gene_id"] == "g2"
    assert result.iloc[2]["gene_id"] == "g3"


def test_drop_all_na_columns():
    df = pd.DataFrame(
        {
            "gene_id": ["ENSG1"],
            "empty_col": [None],
        }
    )
    result = drop_all_na_columns(df)
    assert "empty_col" not in result.columns
    assert "gene_id" in result.columns
