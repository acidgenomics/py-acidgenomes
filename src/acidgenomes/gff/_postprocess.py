"""Post-processing steps applied to all parsed GFF DataFrames.

Ports the logic from R's internal-GenomicRanges.R (.makeGRanges pipeline).
"""

import re

import pandas as pd

# Biotypes that are always removed (present in GFF but not in FASTA).
_REMOVE_BIOTYPES_DEFAULT = frozenset({"LRG_gene"})
# Additional biotypes removed for non-GENCODE providers.
_REMOVE_BIOTYPES_NON_GENCODE = frozenset({"artifact", "TEC"})

# Biotypes classified as noncoding long RNA.
_NONCODING_BIOTYPES = frozenset({"known_ncrna", "lincRNA", "lncRNA", "non_coding"})

# Biotypes classified as small RNA.
_SMALL_RNA_BIOTYPES = frozenset(
    {
        "miRNA",
        "misc_RNA",
        "ribozyme",
        "rRNA",
        "scaRNA",
        "scRNA",
        "snoRNA",
        "snRNA",
        "sRNA",
    }
)

# Biotypes classified as NMD/NSD decaying transcripts.
_DECAYING_BIOTYPES = frozenset({"non_stop_decay", "nonsense_mediated_decay"})


def _broad_class_row(biotype: str | None, seqnames: str | None, gene_name: str | None) -> str:  # noqa: PLR0911
    """Classify one feature into a broadClass string."""
    seqnames = str(seqnames) if seqnames is not None else ""
    gene_name = str(gene_name) if gene_name is not None else ""
    biotype = str(biotype) if biotype is not None else ""

    # Mitochondrial: chromosome "MT*" (case-insensitive) or gene name "mt:" / "mt-"
    if re.match(r"^MT", seqnames, re.IGNORECASE):
        return "mito"
    if re.match(r"^mt[:\-]", gene_name, re.IGNORECASE):
        return "mito"

    if biotype == "protein_coding":
        return "coding"
    if biotype in _NONCODING_BIOTYPES:
        return "noncoding"
    if re.search(r"pseudo", biotype, re.IGNORECASE):
        return "pseudo"
    if biotype in _SMALL_RNA_BIOTYPES:
        return "small"
    if biotype in _DECAYING_BIOTYPES:
        return "decaying"
    if re.match(r"^ig_", biotype, re.IGNORECASE):
        return "ig"
    if re.match(r"^tr_", biotype, re.IGNORECASE):
        return "tcr"
    return "other"


def add_broad_class(df: pd.DataFrame) -> pd.DataFrame:
    """Add a broadClass column using biotype, chromosome, and gene name.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with optional columns: ``gene_biotype``/``tx_biotype``,
        ``seqnames``, ``gene_name``.

    Returns
    -------
    pd.DataFrame
        Input DataFrame with a new ``broad_class`` column appended.
    """
    if "gene_biotype" in df.columns:
        biotype_col: str | None = "gene_biotype"
    elif "tx_biotype" in df.columns:
        biotype_col = "tx_biotype"
    else:
        biotype_col = None
    biotypes = df[biotype_col] if biotype_col else pd.Series([None] * len(df))
    seqnames = df["seqnames"] if "seqnames" in df.columns else pd.Series([None] * len(df))
    gene_names = df["gene_name"] if "gene_name" in df.columns else pd.Series([None] * len(df))

    broad = [
        _broad_class_row(b, s, g) for b, s, g in zip(biotypes, seqnames, gene_names, strict=False)
    ]
    result = df.copy()
    result["broad_class"] = pd.Categorical(broad)
    return result


def remove_unwanted_biotypes(df: pd.DataFrame, *, provider: str) -> pd.DataFrame:
    """Remove experimental biotypes not present in FASTA files.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with optional ``gene_biotype`` column.
    provider : str
        Provider name string (e.g. "Ensembl", "GENCODE").

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame.
    """
    if "gene_biotype" not in df.columns:
        return df
    to_remove = set(_REMOVE_BIOTYPES_DEFAULT)
    if provider != "GENCODE":
        to_remove |= _REMOVE_BIOTYPES_NON_GENCODE
    mask = df["gene_biotype"].isin(list(to_remove))
    return pd.DataFrame(df[~mask]).reset_index(drop=True)


def standardize_column_names(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize column names to camelCase conventions from the R package.

    Handles the key renames documented in R's .standardizeMcolsNames():
    - ``gene_type`` → ``gene_biotype`` (GENCODE uses gene_type)
    - ``biotype`` → ``gene_biotype`` (Ensembl GFF3 uses biotype)
    - ``gene_symbol`` / ``symbol`` → ``gene_name``
    - ``tx_type`` → ``tx_biotype``

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame.

    Returns
    -------
    pd.DataFrame
        DataFrame with renamed columns.
    """
    df = df.copy()
    renames: dict[str, str] = {}

    if "gene_type" in df.columns and "gene_biotype" not in df.columns:
        renames["gene_type"] = "gene_biotype"
    if "biotype" in df.columns and "gene_biotype" not in df.columns:
        renames["biotype"] = "gene_biotype"
    if "tx_type" in df.columns and "tx_biotype" not in df.columns:
        renames["tx_type"] = "tx_biotype"
    if "symbol" in df.columns and "gene_name" not in df.columns:
        renames["symbol"] = "gene_name"
    if "gene_symbol" in df.columns and "gene_name" not in df.columns:
        renames["gene_symbol"] = "gene_name"
    # Drop duplicate symbol if gene_name already exists
    if "symbol" in df.columns and "gene_name" in df.columns:
        df = df.drop(columns=["symbol"])
    if "entrez_id" in df.columns and "ncbi_gene_id" not in df.columns:
        renames["entrez_id"] = "ncbi_gene_id"
    if "artif_dupl" in df.columns:
        renames["artif_dupl"] = "artifactual_duplication"

    df = df.rename(columns=renames)

    # Add gene_name fallback from gene_id
    if "gene_id" in df.columns and "gene_name" not in df.columns:
        df["gene_name"] = df["gene_id"]

    # Add tx_name fallback from transcript_id
    if "transcript_id" in df.columns and "tx_name" not in df.columns:
        df["tx_name"] = df["transcript_id"]

    return df


def handle_versions(
    df: pd.DataFrame,
    *,
    id_col: str,
    version_col: str,
    ignore_version: bool,
) -> pd.DataFrame:
    """Manage versioned identifier columns.

    When ``ignore_version=False``, moves the versioned ID into the primary column
    and stores the unversioned ID in ``{id_col}_no_version``.
    When ``ignore_version=True``, keeps ``{id_col}`` unversioned and stores the
    versioned form in ``{id_col}_version``.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame.
    id_col : str
        Primary identifier column, e.g. ``"gene_id"``.
    version_col : str
        Column holding the versioned identifier, e.g. ``"gene_id_version"``.
    ignore_version : bool
        Whether to strip versions from the primary identifier.

    Returns
    -------
    pd.DataFrame
        Modified DataFrame.
    """
    if id_col not in df.columns or version_col not in df.columns:
        return df
    no_version_col = f"{id_col}_no_version"
    df = df.copy()
    if not ignore_version:
        # Primary ID becomes versioned; unversioned stored separately.
        df[no_version_col] = df[id_col]
        df[id_col] = df[version_col]
    # Always keep version column for downstream use.
    return df


def sort_ranges(df: pd.DataFrame) -> pd.DataFrame:
    """Sort a GFF DataFrame by genomic location (seqnames, start).

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with ``seqnames`` and ``start`` columns.

    Returns
    -------
    pd.DataFrame
        Sorted DataFrame with reset index.
    """
    if "seqnames" not in df.columns or "start" not in df.columns:
        return df
    return df.sort_values(["seqnames", "start"]).reset_index(drop=True)


def drop_all_na_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Remove columns that are entirely NA.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame.

    Returns
    -------
    pd.DataFrame
        DataFrame without all-NA columns.
    """
    return df.dropna(axis=1, how="all")
