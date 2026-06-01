"""FlyBase-specific GTF feature extraction.

FlyBase only supports GTF format (GFF3 is on the denylist).
Gene IDs follow the FBgn pattern; transcript type filtering uses
``pseudogene`` or ``*RNA`` types.

Ports R's .rtracklayerFlybase*() functions from internal-rtracklayer.R.
"""

import re

import pandas as pd

_TRANSCRIPT_TYPE_PATTERN = re.compile(r"^pseudogene$|RNA$", re.IGNORECASE)


def _is_transcript_type(type_val: str) -> bool:
    return bool(_TRANSCRIPT_TYPE_PATTERN.search(type_val))


def extract_genes(df: pd.DataFrame, gff_format: str) -> pd.DataFrame:
    """Extract gene-level features from a FlyBase GTF DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Full parsed GFF DataFrame.
    gff_format : str
        Must be ``"GTF"`` (GFF3 is not supported for FlyBase).

    Returns
    -------
    pd.DataFrame
        Gene-level rows.
    """
    genes = df[df["type"] == "gene"].copy()
    if genes.empty:
        msg = "No gene features found in FlyBase GTF."
        raise ValueError(msg)
    return genes.reset_index(drop=True)


def extract_transcripts(df: pd.DataFrame, gff_format: str) -> pd.DataFrame:
    """Extract transcript-level features from a FlyBase GTF DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Full parsed GFF DataFrame.
    gff_format : str
        Must be ``"GTF"``.

    Returns
    -------
    pd.DataFrame
        Transcript rows (pseudogene and *RNA types).
    """
    if "type" not in df.columns:
        msg = "FlyBase GTF missing 'type' column."
        raise ValueError(msg)
    keep = df["type"].apply(_is_transcript_type)
    txs = df[keep].copy()
    if txs.empty:
        msg = "No transcript features found in FlyBase GTF."
        raise ValueError(msg)
    return txs.reset_index(drop=True)


def extract_exons(df: pd.DataFrame, gff_format: str) -> pd.DataFrame:
    """Extract exon-level features from a FlyBase GTF DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Full parsed GFF DataFrame.
    gff_format : str
        Must be ``"GTF"``.

    Returns
    -------
    pd.DataFrame
        Exon-level rows.
    """
    exons = df[df["type"] == "exon"].copy()
    if exons.empty:
        msg = "No exon features found in FlyBase GTF."
        raise ValueError(msg)
    return exons.reset_index(drop=True)
