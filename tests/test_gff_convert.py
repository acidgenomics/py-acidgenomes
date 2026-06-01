"""Tests for DataFrame-to-GenomicRanges conversion."""

import pandas as pd
import pytest
from genomicranges import GenomicRanges

from acidgenomes.gff._convert import dataframe_to_granges


def _make_gene_df() -> pd.DataFrame:
    return pd.DataFrame({
        "seqnames": ["1", "1", "MT"],
        "start": [65419, 450703, 3307],
        "end": [71585, 451697, 4262],
        "strand": ["+", "-", "+"],
        "gene_id": ["ENSG00000186092", "ENSG00000284733", "ENSG00000198888"],
        "gene_name": ["OR4F5", "OR4F29", "MT-ND1"],
    })


def test_returns_genomic_ranges():
    df = _make_gene_df()
    gr = dataframe_to_granges(df, names_col="gene_id")
    assert isinstance(gr, GenomicRanges)


def test_correct_length():
    df = _make_gene_df()
    gr = dataframe_to_granges(df, names_col="gene_id")
    assert len(gr) == 3


def test_names_set():
    df = _make_gene_df()
    gr = dataframe_to_granges(df, names_col="gene_id")
    assert gr.names is not None
    assert "ENSG00000186092" in list(gr.names)


def test_metadata_attached():
    df = _make_gene_df()
    meta = {"provider": "Ensembl", "organism": "Homo sapiens"}
    gr = dataframe_to_granges(df, names_col="gene_id", metadata=meta)
    assert gr.metadata["provider"] == "Ensembl"


def test_raises_on_empty_df():
    df = pd.DataFrame()
    with pytest.raises(ValueError, match="empty"):
        dataframe_to_granges(df, names_col="gene_id")


def test_raises_on_missing_names_col():
    df = _make_gene_df()
    with pytest.raises(ValueError, match="gene_bad_col"):
        dataframe_to_granges(df, names_col="gene_bad_col")


def test_raises_on_missing_required_cols():
    df = pd.DataFrame({"gene_id": ["ENSG1"]})
    with pytest.raises(ValueError, match="seqnames"):
        dataframe_to_granges(df, names_col="gene_id")


def test_mcols_present():
    df = _make_gene_df()
    gr = dataframe_to_granges(df, names_col="gene_id")
    mcols = gr.mcols
    assert mcols is not None
    assert "gene_name" in mcols.column_names
