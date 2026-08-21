"""Tests for _genomicranges_to_pandas_safe in _classes.py.

Covers the workaround for a confirmed bug in ``genomicranges``
(version 0.8.4): ``GenomicRanges.to_pandas()`` overwrites the range-columns
DataFrame's index with ``gr.names`` (e.g. Ensembl gene IDs) *before*
concatenating it with ``gr.mcols.to_pandas()``, whose index is left at the
default positional ``RangeIndex``. Because the two indices are then
disjoint, ``pd.concat(..., axis=1)`` outer-joins on mismatched labels
instead of stacking columns positionally, silently doubling the row count
with every row missing either its coordinate columns or its metadata
columns.

Confirmed 2026-08-20 against a live Ensembl GTF parse via
``make_ensembl_genes_from_gtf``: both the human and mouse gene annotation
tables built from this package's output were silently corrupted at exactly
2x the true row count.
"""

import pandas as pd

from acidgenomes._classes import _genomicranges_to_pandas_safe
from acidgenomes.gff._convert import dataframe_to_granges


def _make_gene_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "seqnames": ["1", "1", "MT"],
            "start": [65419, 450703, 3307],
            "end": [71585, 451697, 4262],
            "strand": ["+", "-", "+"],
            "gene_id": ["ENSG00000186092", "ENSG00000284733", "ENSG00000198888"],
            "gene_name": ["OR4F5", "OR4F29", "MT-ND1"],
        }
    )


def test_upstream_to_pandas_reproduces_the_bug() -> None:
    """Confirm the bug is real: GenomicRanges.to_pandas() doubles the rows.

    Uses dataframe_to_granges -- the real construction path
    make_granges_from_gff uses for every provider -- not a hand-built
    GenomicRanges, so this reflects the actual production shape.
    """
    df = _make_gene_df()
    gr = dataframe_to_granges(df, names_col="gene_id")
    buggy = gr.to_pandas()
    assert len(buggy) == 2 * len(gr)
    assert buggy["gene_id"].isna().sum() == len(gr)
    assert buggy["seqnames"].isna().sum() == len(gr)


def test_safe_conversion_returns_one_row_per_range() -> None:
    """The workaround returns exactly len(gr) rows, not double."""
    df = _make_gene_df()
    gr = dataframe_to_granges(df, names_col="gene_id")
    out = _genomicranges_to_pandas_safe(gr)
    assert len(out) == len(gr)


def test_safe_conversion_has_no_nan_split_artifacts() -> None:
    """No row is missing either its identity columns or its coordinate columns."""
    df = _make_gene_df()
    gr = dataframe_to_granges(df, names_col="gene_id")
    out = _genomicranges_to_pandas_safe(gr)
    assert bool(out["gene_id"].notna().all())
    assert bool(out["seqnames"].notna().all())
    assert list(out["gene_id"]) == df["gene_id"].tolist()
    assert list(out["gene_name"]) == df["gene_name"].tolist()


def test_safe_conversion_index_matches_names() -> None:
    """The returned index is gr.names, applied only after the positional concat."""
    df = _make_gene_df()
    gr = dataframe_to_granges(df, names_col="gene_id")
    out = _genomicranges_to_pandas_safe(gr)
    assert list(out.index) == df["gene_id"].tolist()
