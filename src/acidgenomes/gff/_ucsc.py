"""UCSC-specific GTF feature extraction.

UCSC GTF files use a source column value of ``knownGene``, ``ncbiRefSeq``,
``refGene``, or ``ensGene``. The bcbio-nextgen ``ref-transcripts.gtf`` is
also treated as UCSC.

Only GTF is supported for UCSC in this package.
"""

import pandas as pd


def extract_genes(df: pd.DataFrame, gff_format: str) -> pd.DataFrame:
    """Extract gene-level features from a UCSC GTF DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Full parsed GFF DataFrame.
    gff_format : str
        Must be ``"GTF"``.

    Returns
    -------
    pd.DataFrame
        Gene-level rows (type == "transcript" since UCSC GTF lacks ``gene``
        feature type; deduplicated by gene_id).
    """
    # UCSC GTF files typically lack explicit "gene" feature rows.
    # Derive gene-level information from transcript rows deduplicated by gene_id.
    if "type" in df.columns and "gene" in df["type"].values:
        genes = pd.DataFrame(df[df["type"] == "gene"]).copy()
    else:
        # Collapse transcripts to gene level.
        txs = (
            pd.DataFrame(df[df["type"] == "transcript"]).copy()
            if "type" in df.columns
            else df.copy()
        )
        if "gene_id" not in txs.columns:
            msg = "UCSC GTF missing 'gene_id' column."
            raise ValueError(msg)
        genes = txs.drop_duplicates(subset=["gene_id"]).copy()
    if genes.empty:
        msg = "No gene features found in UCSC GTF."
        raise ValueError(msg)
    return genes.reset_index(drop=True)


def extract_transcripts(df: pd.DataFrame, gff_format: str) -> pd.DataFrame:
    """Extract transcript-level features from a UCSC GTF DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Full parsed GFF DataFrame.
    gff_format : str
        Must be ``"GTF"``.

    Returns
    -------
    pd.DataFrame
        Transcript-level rows.
    """
    if "type" not in df.columns:
        msg = "UCSC GTF missing 'type' column."
        raise ValueError(msg)
    txs = pd.DataFrame(df[df["type"] == "transcript"]).copy()
    if txs.empty:
        msg = "No transcript features found in UCSC GTF."
        raise ValueError(msg)
    return txs.reset_index(drop=True)


def extract_exons(df: pd.DataFrame, gff_format: str) -> pd.DataFrame:
    """Extract exon-level features from a UCSC GTF DataFrame.

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
    exons = pd.DataFrame(df[df["type"] == "exon"]).copy()
    if exons.empty:
        msg = "No exon features found in UCSC GTF."
        raise ValueError(msg)
    return exons.reset_index(drop=True)
