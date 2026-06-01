"""Ensembl-specific GFF3/GTF feature extraction.

Ports the logic from R's .rtracklayerEnsembl*() functions in
internal-rtracklayer.R.
"""

import pandas as pd


def extract_genes(df: pd.DataFrame, gff_format: str) -> pd.DataFrame:
    """Extract gene-level features from an Ensembl GFF/GTF DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Full parsed GFF DataFrame.
    gff_format : str
        ``"GTF"`` or ``"GFF3"``.

    Returns
    -------
    pd.DataFrame
        Gene-level rows with standardized columns.
    """
    if gff_format == "GTF":
        return _genes_gtf(df)
    return _genes_gff3(df)


def extract_transcripts(df: pd.DataFrame, gff_format: str) -> pd.DataFrame:
    """Extract transcript-level features from an Ensembl GFF/GTF DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Full parsed GFF DataFrame.
    gff_format : str
        ``"GTF"`` or ``"GFF3"``.

    Returns
    -------
    pd.DataFrame
        Transcript-level rows with gene metadata joined.
    """
    if gff_format == "GTF":
        return _transcripts_gtf(df)
    return _transcripts_gff3(df)


def extract_exons(df: pd.DataFrame, gff_format: str) -> pd.DataFrame:
    """Extract exon-level features from an Ensembl GFF/GTF DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Full parsed GFF DataFrame.
    gff_format : str
        ``"GTF"`` or ``"GFF3"``.

    Returns
    -------
    pd.DataFrame
        Exon-level rows with transcript metadata joined.
    """
    if gff_format == "GTF":
        return _exons_gtf(df)
    return _exons_gff3(df)


def _genes_gtf(df: pd.DataFrame) -> pd.DataFrame:
    genes = df[df["type"] == "gene"].copy()
    if genes.empty:
        msg = "No gene features found in Ensembl GTF."
        raise ValueError(msg)
    # Build versioned ID from gene_id + gene_version columns.
    if "gene_version" in genes.columns:
        genes["gene_id_version"] = genes["gene_id"] + "." + genes["gene_version"].astype(str)
        genes = genes.drop(columns=["gene_version"])
    else:
        genes["gene_id_version"] = genes["gene_id"]
    return genes.reset_index(drop=True)


def _genes_gff3(df: pd.DataFrame) -> pd.DataFrame:
    # Ensembl GFF3: gene rows have a non-NA gene_id attribute.
    if "gene_id" in df.columns:
        genes = df[df["gene_id"].notna()].copy()
    else:
        genes = df[df["type"] == "gene"].copy()
    if genes.empty:
        msg = "No gene features found in Ensembl GFF3."
        raise ValueError(msg)
    # Rename columns to match GTF conventions.
    if "Name" in genes.columns and "gene_name" not in genes.columns:
        genes = genes.rename(columns={"Name": "gene_name"})
    if "biotype" in genes.columns and "gene_biotype" not in genes.columns:
        genes = genes.rename(columns={"biotype": "gene_biotype"})
    # Build versioned ID.
    if "version" in genes.columns:
        genes["gene_id_version"] = genes["gene_id"] + "." + genes["version"].astype(str)
        genes = genes.drop(columns=["version"])
    else:
        genes["gene_id_version"] = genes["gene_id"]
    for col in ("Alias", "ID", "Parent"):
        if col in genes.columns:
            genes = genes.drop(columns=[col])
    return genes.reset_index(drop=True)


def _transcripts_gtf(df: pd.DataFrame) -> pd.DataFrame:
    genes = _genes_gtf(df)
    txs = df[df["type"] == "transcript"].copy()
    if txs.empty:
        msg = "No transcript features found in Ensembl GTF."
        raise ValueError(msg)
    for key in ("gene", "transcript"):
        id_col = f"{key}_id"
        ver_col = f"{key}_version"
        ver_id_col = f"{key}_id_version"
        if ver_col in txs.columns:
            txs[ver_id_col] = txs[id_col] + "." + txs[ver_col].astype(str)
            txs = txs.drop(columns=[ver_col])
        else:
            txs[ver_id_col] = txs[id_col]
    # Join gene metadata.
    gene_cols = ["gene_id"] + [c for c in genes.columns if c not in txs.columns and c != "gene_id"]
    txs = txs.merge(genes[gene_cols], on="gene_id", how="left")
    return txs.reset_index(drop=True)


def _transcripts_gff3(df: pd.DataFrame) -> pd.DataFrame:
    genes = _genes_gff3(df)
    if "transcript_id" in df.columns:
        txs = df[df["transcript_id"].notna()].copy()
    else:
        txs = pd.DataFrame()
    if txs.empty:
        msg = "No transcript features found in Ensembl GFF3."
        raise ValueError(msg)
    if "Name" in txs.columns and "transcript_name" not in txs.columns:
        txs = txs.rename(columns={"Name": "transcript_name"})
    if "biotype" in txs.columns and "tx_biotype" not in txs.columns:
        txs = txs.rename(columns={"biotype": "tx_biotype"})
    # Extract gene_id from "Parent" field which is prefixed "gene:".
    if "Parent" in txs.columns:
        txs["gene_id"] = txs["Parent"].str.replace("^gene:", "", regex=True)
    if "version" in txs.columns:
        txs["transcript_id_version"] = txs["transcript_id"] + "." + txs["version"].astype(str)
        txs = txs.drop(columns=["version"])
    else:
        txs["transcript_id_version"] = txs["transcript_id"]
    for col in ("Alias", "ID", "Parent"):
        if col in txs.columns:
            txs = txs.drop(columns=[col])
    gene_cols = ["gene_id"] + [c for c in genes.columns if c not in txs.columns and c != "gene_id"]
    txs = txs.merge(genes[gene_cols], on="gene_id", how="left")
    return txs.reset_index(drop=True)


def _exons_gtf(df: pd.DataFrame) -> pd.DataFrame:
    txs = _transcripts_gtf(df)
    exons = df[df["type"] == "exon"].copy()
    if exons.empty:
        msg = "No exon features found in Ensembl GTF."
        raise ValueError(msg)
    for key in ("exon", "gene", "transcript"):
        id_col = f"{key}_id" if key != "transcript" else "transcript_id"
        ver_col = f"{key}_version"
        ver_id_col = f"{id_col}_version"
        if ver_col in exons.columns:
            exons[ver_id_col] = exons[id_col] + "." + exons[ver_col].astype(str)
            exons = exons.drop(columns=[ver_col])
        else:
            exons[ver_id_col] = exons[id_col]
    tx_extra = [c for c in txs.columns if c not in exons.columns and c != "transcript_id"]
    tx_cols = ["transcript_id", *tx_extra]
    exons = exons.merge(txs[tx_cols], on="transcript_id", how="left")
    return exons.reset_index(drop=True)


def _exons_gff3(df: pd.DataFrame) -> pd.DataFrame:
    txs = _transcripts_gff3(df)
    exons = df[df["exon_id"].notna()].copy() if "exon_id" in df.columns else pd.DataFrame()
    if exons.empty:
        msg = "No exon features found in Ensembl GFF3."
        raise ValueError(msg)
    if "Name" in exons.columns and "exon_name" not in exons.columns:
        exons = exons.rename(columns={"Name": "exon_name"})
    # Extract transcript_id from Parent field prefixed "transcript:".
    if "Parent" in exons.columns:
        exons["transcript_id"] = exons["Parent"].str.replace("^transcript:", "", regex=True)
    if "version" in exons.columns:
        exons["exon_id_version"] = exons["exon_id"] + "." + exons["version"].astype(str)
        exons = exons.drop(columns=["version"])
    else:
        exons["exon_id_version"] = exons["exon_id"]
    for col in ("Alias", "ID", "Parent"):
        if col in exons.columns:
            exons = exons.drop(columns=[col])
    tx_extra = [c for c in txs.columns if c not in exons.columns and c != "transcript_id"]
    tx_cols = ["transcript_id", *tx_extra]
    exons = exons.merge(txs[tx_cols], on="transcript_id", how="left")
    return exons.reset_index(drop=True)
