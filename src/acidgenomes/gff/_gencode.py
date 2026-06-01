"""GENCODE-specific GFF3/GTF feature extraction.

GENCODE uses ``gene_type`` instead of ``gene_biotype``, and primary identifier
fields already contain the version (e.g. ``ENSG00000000003.14``). The
unversioned form is derived by stripping the version suffix.

Ports the logic from R's .rtracklayerGencode*() functions in
internal-rtracklayer.R.
"""

import pandas as pd

from acidgenomes._strip_versions import (
    strip_exon_versions,
    strip_gene_versions,
    strip_transcript_versions,
)


def extract_genes(df: pd.DataFrame, gff_format: str) -> pd.DataFrame:
    """Extract gene-level features from a GENCODE GFF/GTF DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Full parsed GFF DataFrame.
    gff_format : str
        ``"GTF"`` or ``"GFF3"``.

    Returns
    -------
    pd.DataFrame
        Gene-level rows.
    """
    if gff_format == "GTF":
        return _genes_gtf(df)
    return _genes_gff3(df)


def extract_transcripts(df: pd.DataFrame, gff_format: str) -> pd.DataFrame:
    """Extract transcript-level features from a GENCODE GFF/GTF DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Full parsed GFF DataFrame.
    gff_format : str
        ``"GTF"`` or ``"GFF3"``.

    Returns
    -------
    pd.DataFrame
        Transcript-level rows.
    """
    if gff_format == "GTF":
        return _transcripts_gtf(df)
    return _transcripts_gff3(df)


def extract_exons(df: pd.DataFrame, gff_format: str) -> pd.DataFrame:
    """Extract exon-level features from a GENCODE GFF/GTF DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Full parsed GFF DataFrame.
    gff_format : str
        ``"GTF"`` or ``"GFF3"``.

    Returns
    -------
    pd.DataFrame
        Exon-level rows.
    """
    if gff_format == "GTF":
        return _exons_gtf(df)
    return _exons_gff3(df)


def _genes_gtf(df: pd.DataFrame) -> pd.DataFrame:
    genes = df[df["type"] == "gene"].copy()
    if genes.empty:
        msg = "No gene features found in GENCODE GTF."
        raise ValueError(msg)
    # gene_id already contains versioned ID in GENCODE.
    genes["gene_id_version"] = genes["gene_id"]
    genes["gene_id"] = strip_gene_versions(genes["gene_id"].tolist())
    for col in ("ont",):
        if col in genes.columns:
            genes = genes.drop(columns=[col])
    return genes.reset_index(drop=True)


def _genes_gff3(df: pd.DataFrame) -> pd.DataFrame:
    genes = df[df["type"] == "gene"].copy()
    if genes.empty:
        msg = "No gene features found in GENCODE GFF3."
        raise ValueError(msg)
    # In GENCODE GFF3, the ID attribute holds the versioned gene_id.
    if "ID" in genes.columns:
        genes["gene_id_version"] = genes["ID"]
        genes["gene_id"] = strip_gene_versions(genes["gene_id_version"].tolist())
    else:
        genes["gene_id_version"] = genes["gene_id"]
        genes["gene_id"] = strip_gene_versions(genes["gene_id"].tolist())
    for col in ("ID", "Parent", "ont"):
        if col in genes.columns:
            genes = genes.drop(columns=[col])
    return genes.reset_index(drop=True)


def _transcripts_gtf(df: pd.DataFrame) -> pd.DataFrame:
    txs = df[df["type"] == "transcript"].copy()
    if txs.empty:
        msg = "No transcript features found in GENCODE GTF."
        raise ValueError(msg)
    # transcript_id and gene_id already versioned.
    txs["transcript_id_version"] = txs["transcript_id"]
    txs["transcript_id"] = strip_transcript_versions(txs["transcript_id"].tolist())
    txs["gene_id_version"] = txs["gene_id"]
    txs["gene_id"] = strip_gene_versions(txs["gene_id"].tolist())
    for col in ("ont",):
        if col in txs.columns:
            txs = txs.drop(columns=[col])
    return txs.reset_index(drop=True)


def _transcripts_gff3(df: pd.DataFrame) -> pd.DataFrame:
    txs = df[df["type"] == "transcript"].copy()
    if txs.empty:
        msg = "No transcript features found in GENCODE GFF3."
        raise ValueError(msg)
    if "ID" in txs.columns:
        txs["transcript_id_version"] = txs["ID"]
        txs["transcript_id"] = strip_transcript_versions(txs["transcript_id_version"].tolist())
    else:
        txs["transcript_id_version"] = txs["transcript_id"]
        txs["transcript_id"] = strip_transcript_versions(txs["transcript_id"].tolist())
    if "Parent" in txs.columns:
        txs["gene_id_version"] = txs["Parent"].astype(str)
        txs["gene_id"] = strip_gene_versions(txs["gene_id_version"].tolist())
    else:
        txs["gene_id_version"] = txs["gene_id"]
        txs["gene_id"] = strip_gene_versions(txs["gene_id"].tolist())
    for col in ("ID", "Parent", "ont"):
        if col in txs.columns:
            txs = txs.drop(columns=[col])
    return txs.reset_index(drop=True)


def _exons_gtf(df: pd.DataFrame) -> pd.DataFrame:
    exons = df[df["type"] == "exon"].copy()
    if exons.empty:
        msg = "No exon features found in GENCODE GTF."
        raise ValueError(msg)
    exons["exon_id_version"] = exons["exon_id"]
    exons["exon_id"] = strip_exon_versions(exons["exon_id"].tolist())
    exons["transcript_id_version"] = exons["transcript_id"]
    exons["transcript_id"] = strip_transcript_versions(exons["transcript_id"].tolist())
    exons["gene_id_version"] = exons["gene_id"]
    exons["gene_id"] = strip_gene_versions(exons["gene_id"].tolist())
    for col in ("ont",):
        if col in exons.columns:
            exons = exons.drop(columns=[col])
    return exons.reset_index(drop=True)


def _exons_gff3(df: pd.DataFrame) -> pd.DataFrame:
    exons = df[df["type"] == "exon"].copy()
    if exons.empty:
        msg = "No exon features found in GENCODE GFF3."
        raise ValueError(msg)
    exons["exon_id_version"] = exons["exon_id"]
    exons["exon_id"] = strip_exon_versions(exons["exon_id"].tolist())
    exons["transcript_id_version"] = exons["transcript_id"]
    exons["transcript_id"] = strip_transcript_versions(exons["transcript_id"].tolist())
    exons["gene_id_version"] = exons["gene_id"]
    exons["gene_id"] = strip_gene_versions(exons["gene_id"].tolist())
    for col in ("ID", "Parent", "ont"):
        if col in exons.columns:
            exons = exons.drop(columns=[col])
    return exons.reset_index(drop=True)
