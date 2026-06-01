"""DataFrame-to-GenomicRanges conversion bridge."""

from typing import Any

import pandas as pd
from biocframe import BiocFrame
from genomicranges import GenomicRanges
from iranges import IRanges


def dataframe_to_granges(
    df: pd.DataFrame,
    *,
    names_col: str,
    metadata: dict[str, Any] | None = None,
) -> GenomicRanges:
    """Convert a GFF DataFrame to a BiocPy GenomicRanges object.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with at minimum ``seqnames``, ``start``, ``end``, ``strand``
        columns (1-based GFF coordinates).
    names_col : str
        Column name to use as range names (e.g. ``"gene_id"``).
    metadata : dict or None
        Arbitrary metadata dict attached to the returned GenomicRanges object.

    Returns
    -------
    GenomicRanges
        BiocPy GenomicRanges object with names from ``names_col`` and all
        remaining DataFrame columns as metadata columns.
    """
    if df.empty:
        msg = "Cannot convert an empty DataFrame to GenomicRanges."
        raise ValueError(msg)

    required = {"seqnames", "start", "end"}
    missing = required - set(df.columns)
    if missing:
        msg = f"DataFrame missing required columns: {', '.join(sorted(missing))}"
        raise ValueError(msg)

    if names_col not in df.columns:
        msg = f"Names column {names_col!r} not found in DataFrame."
        raise ValueError(msg)

    # GFF is 1-based inclusive; IRanges uses 0-based start + width.
    starts = (df["start"].astype(int) - 1).tolist()
    widths = (df["end"].astype(int) - df["start"].astype(int) + 1).tolist()
    iranges = IRanges(start=starts, width=widths)

    seqnames = df["seqnames"].astype(str).tolist()
    strand_col = "strand" if "strand" in df.columns else None
    strand = df[strand_col].astype(str).tolist() if strand_col else None

    names = df[names_col].astype(str).tolist()

    # All columns except the genomic coordinate columns become mcols.
    coord_cols = {
        "seqnames", "start", "end", "strand", "score", "frame", "source", "type", "format"
    }
    mcol_names = [c for c in df.columns if c not in coord_cols]
    mcols: dict[str, list] = {}
    for col in mcol_names:
        mcols[col] = df[col].tolist()

    bioc_mcols = BiocFrame(mcols, number_of_rows=len(df)) if mcols else None

    gr = GenomicRanges(
        seqnames=seqnames,
        ranges=iranges,
        strand=strand,
        mcols=bioc_mcols,
        names=names,
        metadata=metadata or {},
    )
    return gr
