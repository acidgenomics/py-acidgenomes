"""WormBase-specific GTF feature extraction.

WormBase only supports GTF format (GFF3 is on the denylist).
Validates that gene IDs match the ``WBGene[0-9]{8}`` pattern and drops
malformed entries with a warning.

Ports R's .rtracklayerWormbase*() functions from internal-rtracklayer.R.
"""

import re
import warnings

import pandas as pd

_WBGENE_PATTERN = re.compile(r"^WBGene\d{8}$")


def _validate_gene_ids(df: pd.DataFrame, level: str) -> pd.DataFrame:
    """Drop rows where gene_id does not match WBGene pattern."""
    if "gene_id" not in df.columns:
        return df
    valid = df["gene_id"].str.match(_WBGENE_PATTERN, na=False)
    n_dropped = (~valid).sum()
    if n_dropped > 0:
        warnings.warn(
            f"Dropping {n_dropped} {level} entries with malformed WormBase gene_id.",
            stacklevel=3,
        )
    return pd.DataFrame(df[valid]).reset_index(drop=True)


def extract_genes(df: pd.DataFrame, gff_format: str) -> pd.DataFrame:
    """Extract gene-level features from a WormBase GTF DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Full parsed GFF DataFrame.
    gff_format : str
        Must be ``"GTF"`` (GFF3 not supported).

    Returns
    -------
    pd.DataFrame
        Gene-level rows with validated WBGene IDs.
    """
    genes = pd.DataFrame(df[df["type"] == "gene"]).copy()
    if genes.empty:
        msg = "No gene features found in WormBase GTF."
        raise ValueError(msg)
    genes = _validate_gene_ids(genes, "gene")
    return genes


def extract_transcripts(df: pd.DataFrame, gff_format: str) -> pd.DataFrame:
    """Extract transcript-level features from a WormBase GTF DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Full parsed GFF DataFrame.
    gff_format : str
        Must be ``"GTF"``.

    Returns
    -------
    pd.DataFrame
        Transcript-level rows with gene metadata joined.
    """
    genes = extract_genes(df, gff_format)
    gene_cols_prefix = [c for c in genes.columns if c.startswith("gene_")]

    txs = pd.DataFrame(df[df["type"] == "transcript"]).copy()
    if txs.empty:
        msg = "No transcript features found in WormBase GTF."
        raise ValueError(msg)
    txs = _validate_gene_ids(txs, "transcript")

    # Join gene metadata.
    gene_extra = [c for c in gene_cols_prefix if c not in txs.columns and c != "gene_id"]
    join_cols = ["gene_id", *gene_extra]
    txs = txs.merge(genes[join_cols], on="gene_id", how="left")
    return txs.reset_index(drop=True)


def extract_exons(df: pd.DataFrame, gff_format: str) -> pd.DataFrame:
    """Extract exon-level features from a WormBase GTF DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Full parsed GFF DataFrame.
    gff_format : str
        Must be ``"GTF"``.

    Returns
    -------
    pd.DataFrame
        Exon-level rows with validated WBGene IDs.
    """
    exons = pd.DataFrame(df[df["type"] == "exon"]).copy()
    if exons.empty:
        msg = "No exon features found in WormBase GTF."
        raise ValueError(msg)
    exons = _validate_gene_ids(exons, "exon")
    return exons
