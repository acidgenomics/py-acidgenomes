r"""RefSeq-specific GFF3/GTF feature extraction.

RefSeq stores NCBI gene IDs in the ``Dbxref`` attribute (GFF3) or the
``db_xref`` attribute (GTF). Transcript IDs follow the pattern
``[A-Z]{2}_[0-9]+\.[0-9]+`` (e.g. ``NM_000014.6``).

Ports the logic from R's .rtracklayerRefseq*() functions in
internal-rtracklayer.R.
"""

import re

import pandas as pd

_TRANSCRIPT_ID_PATTERN = re.compile(r"^[A-Z]{2}_[0-9]+\.[0-9]+$")
_GENE_ID_RE = re.compile(r"GeneID:(\d+)")


def _extract_ncbi_gene_id_from_dbxref(dbxref_value: object) -> int | None:
    """Extract NCBI gene ID from a Dbxref attribute value (GFF3)."""
    if dbxref_value is None:
        return None
    entries = dbxref_value if isinstance(dbxref_value, list) else [str(dbxref_value)]
    for entry in entries:
        m = _GENE_ID_RE.search(str(entry))
        if m:
            return int(m.group(1))
    return None


def _extract_ncbi_gene_ids_from_df_gff3(df: pd.DataFrame) -> pd.DataFrame:
    """Build a mapping of ID → ncbi_gene_id from RefSeq GFF3 gene rows."""
    gene_rows = df[df["gbkey"] == "Gene"].copy() if "gbkey" in df.columns else pd.DataFrame()
    if gene_rows.empty:
        return pd.DataFrame(columns=["ID", "ncbi_gene_id"])
    gene_rows = gene_rows[["ID", "Dbxref"]].dropna(subset=["ID"])
    gene_rows["ncbi_gene_id"] = gene_rows["Dbxref"].apply(_extract_ncbi_gene_id_from_dbxref)
    return gene_rows[["ID", "ncbi_gene_id"]].dropna().drop_duplicates(subset=["ID"])


def _extract_ncbi_gene_ids_from_df_gtf(df: pd.DataFrame) -> pd.DataFrame:
    """Build a mapping of gene_id → ncbi_gene_id from RefSeq GTF exon rows."""
    if "db_xref" not in df.columns:
        return pd.DataFrame(columns=["gene_id", "ncbi_gene_id"])
    exon_rows = df[df["type"] == "exon"].copy() if "type" in df.columns else df.copy()
    keep = exon_rows["db_xref"].str.contains("GeneID:", na=False)
    exon_rows = exon_rows[keep][["gene_id", "db_xref"]].dropna()
    ncbi_ids_extracted = exon_rows["db_xref"].str.extract(r"GeneID:(\d+)")
    exon_rows["ncbi_gene_id"] = ncbi_ids_extracted.astype(float).astype("Int64")
    result = exon_rows[["gene_id", "ncbi_gene_id"]].dropna().drop_duplicates(subset=["gene_id"])
    return result


def extract_genes(df: pd.DataFrame, gff_format: str) -> pd.DataFrame:
    """Extract gene-level features from a RefSeq GFF/GTF DataFrame.

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
    """Extract transcript-level features from a RefSeq GFF/GTF DataFrame.

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
    """Extract exon-level features from a RefSeq GFF/GTF DataFrame.

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
    # NCBI gene IDs must be extracted from exon rows before gene filtering.
    ncbi_ids = _extract_ncbi_gene_ids_from_df_gtf(df)
    genes = df[df["type"] == "gene"].copy()
    if genes.empty:
        msg = "No gene features found in RefSeq GTF."
        raise ValueError(msg)
    if not ncbi_ids.empty:
        genes = genes.merge(ncbi_ids, on="gene_id", how="left")
    genes = genes.rename(columns={"gene_id": "parent_gene_id", "gene": "gene_id"})
    for col in ("gbkey", "note"):
        if col in genes.columns:
            genes = genes.drop(columns=[col])
    return genes.reset_index(drop=True)


def _genes_gff3(df: pd.DataFrame) -> pd.DataFrame:
    if "gbkey" not in df.columns:
        msg = "RefSeq GFF3 missing 'gbkey' column."
        raise ValueError(msg)
    genes = df[df["gbkey"] == "Gene"].copy()
    if genes.empty:
        msg = "No gene features found in RefSeq GFF3."
        raise ValueError(msg)
    ncbi_ids = _extract_ncbi_gene_ids_from_df_gff3(df)
    if not ncbi_ids.empty:
        genes = genes.merge(ncbi_ids, on="ID", how="left")
    # ID column → parent_gene_id (strip "gene-" prefix), gene column → gene_id.
    genes = genes.rename(columns={"ID": "parent_gene_id", "gene": "gene_id"})
    genes["parent_gene_id"] = genes["parent_gene_id"].str.replace("^gene-", "", regex=True)
    for col in (
        "Name", "Note", "Parent", "end_range", "experiment", "gbkey", "start_range", "transl_except"
    ):
        if col in genes.columns:
            genes = genes.drop(columns=[col])
    return genes.reset_index(drop=True)


def _transcripts_gtf(df: pd.DataFrame) -> pd.DataFrame:
    genes = _genes_gtf(df)
    gene_meta_cols = ["parent_gene_id"] + [
        c for c in ("description", "gene_biotype", "ncbi_gene_id") if c in genes.columns
    ]
    gene_meta = (
        genes[gene_meta_cols]
        .dropna(subset=["parent_gene_id"])
        .drop_duplicates(subset=["parent_gene_id"])
    )

    txs = df[df["type"] == "transcript"].copy()
    if txs.empty:
        msg = "No transcript features found in RefSeq GTF."
        raise ValueError(msg)
    keep = txs["transcript_id"].str.match(_TRANSCRIPT_ID_PATTERN, na=False)
    txs = txs[keep]
    txs = txs.rename(columns={"gene_id": "parent_gene_id", "gene": "gene_id"})
    txs = txs.merge(gene_meta, on="parent_gene_id", how="left")
    for col in ("experiment", "gbkey", "note"):
        if col in txs.columns:
            txs = txs.drop(columns=[col])
    return txs.reset_index(drop=True)


def _transcripts_gff3(df: pd.DataFrame) -> pd.DataFrame:
    genes = _genes_gff3(df)
    gene_meta_cols = ["parent_gene_id"] + [
        c for c in ("description", "gene_biotype", "ncbi_gene_id") if c in genes.columns
    ]
    gene_meta = (
        genes[gene_meta_cols]
        .dropna(subset=["parent_gene_id"])
        .drop_duplicates(subset=["parent_gene_id"])
    )

    if "transcript_id" in df.columns:
        txs = df[df["transcript_id"].notna()].copy()
    else:
        txs = pd.DataFrame()
    if txs.empty:
        msg = "No transcript features found in RefSeq GFF3."
        raise ValueError(msg)
    keep = txs["transcript_id"].str.match(_TRANSCRIPT_ID_PATTERN, na=False)
    txs = txs[keep]
    # Extract parent_gene_id from Parent attribute: strip "gene-" prefix.
    if "Parent" in txs.columns:
        txs["parent_gene_id"] = txs["Parent"].astype(str).str.replace("^gene-", "", regex=True)
    txs = txs.rename(columns={"gene": "gene_id"})
    txs = txs.merge(gene_meta, on="parent_gene_id", how="left")
    for col in (
        "ID", "Name", "Note", "Parent",
        "end_range", "experiment", "gbkey", "start_range", "transl_except",
    ):
        if col in txs.columns:
            txs = txs.drop(columns=[col])
    return txs.reset_index(drop=True)


def _exons_gff3(df: pd.DataFrame) -> pd.DataFrame:
    exons = df[df["type"] == "exon"].copy()
    if exons.empty:
        msg = "No exon features found in RefSeq GFF3."
        raise ValueError(msg)
    if "transcript_id" in exons.columns:
        keep = exons["transcript_id"].str.match(_TRANSCRIPT_ID_PATTERN, na=False)
    else:
        keep = pd.Series([False] * len(exons))
    exons = exons[keep]
    exons["ncbi_gene_id"] = exons["Dbxref"].apply(_extract_ncbi_gene_id_from_dbxref)
    exons = exons.rename(columns={"gene": "gene_id"})
    for col in ("Dbxref", "ID", "Name", "Note", "Parent", "gbkey"):
        if col in exons.columns:
            exons = exons.drop(columns=[col])
    return exons.reset_index(drop=True)


def _exons_gtf(df: pd.DataFrame) -> pd.DataFrame:
    ncbi_ids = _extract_ncbi_gene_ids_from_df_gtf(df)
    exons = df[df["type"] == "exon"].copy()
    if exons.empty:
        msg = "No exon features found in RefSeq GTF."
        raise ValueError(msg)
    if "transcript_id" in exons.columns:
        keep = exons["transcript_id"].str.match(_TRANSCRIPT_ID_PATTERN, na=False)
    else:
        keep = pd.Series([False] * len(exons))
    exons = exons[keep]
    if not ncbi_ids.empty:
        exons = exons.merge(ncbi_ids, on="gene_id", how="left")
    exons = exons.rename(columns={"gene_id": "parent_gene_id", "gene": "gene_id"})
    for col in ("gbkey", "note"):
        if col in exons.columns:
            exons = exons.drop(columns=[col])
    return exons.reset_index(drop=True)
