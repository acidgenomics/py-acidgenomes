"""GFF/GTF parsing subpackage.

Public API
----------
make_granges_from_gff
    Parse any GFF3 or GTF file into a BiocPy GenomicRanges object, with
    automatic provider detection and provider-specific feature extraction.
"""

from pathlib import Path
from typing import Literal

from genomicranges import GenomicRanges

from acidgenomes.gff._convert import dataframe_to_granges
from acidgenomes.gff._metadata import get_gff_metadata
from acidgenomes.gff._parser import parse_gff
from acidgenomes.gff._postprocess import (
    add_broad_class,
    drop_all_na_columns,
    handle_versions,
    remove_unwanted_biotypes,
    sort_ranges,
    standardize_column_names,
)
from acidgenomes.gff._providers import Provider, is_supported_gff

# Map (provider, format) → (module, extract_fn_name)
_PROVIDER_MODULES = {
    Provider.ENSEMBL: "acidgenomes.gff._ensembl",
    Provider.FLYBASE: "acidgenomes.gff._flybase",
    Provider.GENCODE: "acidgenomes.gff._gencode",
    Provider.REFSEQ: "acidgenomes.gff._refseq",
    Provider.UCSC: "acidgenomes.gff._ucsc",
    Provider.WORMBASE: "acidgenomes.gff._wormbase",
}

_LEVEL_TO_NAMES_COL = {
    "genes": "gene_id",
    "transcripts": "transcript_id",
    "exons": "exon_id",
}

# Feature types to pre-filter during parsing per level (optimization).
_LEVEL_FEATURE_FILTER: dict[str, set[str]] = {
    "genes": {"gene", "region"},
    "transcripts": {
        "gene",
        "transcript",
        "mRNA",
        "pseudogene",
        "ncRNA",
        "rRNA",
        "tRNA",
        "snRNA",
        "snoRNA",
        "miRNA",
    },
    "exons": {"gene", "transcript", "mRNA", "exon", "pseudogene", "ncRNA", "rRNA", "tRNA"},
}


def make_granges_from_gff(
    file: str | Path,
    *,
    level: Literal["genes", "transcripts", "exons"] = "genes",
    ignore_version: bool = False,
    extra_mcols: bool = False,
) -> GenomicRanges:
    """Parse a GFF3 or GTF file into a BiocPy GenomicRanges object.

    Automatically detects the annotation provider (Ensembl, GENCODE, RefSeq,
    UCSC, FlyBase, WormBase) from directives and filename, then applies
    provider-specific extraction and post-processing.

    Parameters
    ----------
    file : str or Path
        Path to a GFF3 or GTF file (plain or gzip-compressed).
    level : {"genes", "transcripts", "exons"}
        Feature level to extract.
    ignore_version : bool
        If ``True``, strip version suffixes from identifiers (e.g.
        ``ENSG00000000003.14`` → ``ENSG00000000003``). If ``False`` (default),
        the primary identifier column contains the versioned form and the
        unversioned form is stored in a ``*_no_version`` column.
    extra_mcols : bool
        If ``True``, add ``broad_class`` annotations and (for Ensembl/GENCODE)
        fetch additional metadata from Ensembl FTP.

    Returns
    -------
    GenomicRanges
        BiocPy GenomicRanges object with feature metadata in mcols and file
        metadata attached to the object's ``metadata`` attribute.

    Raises
    ------
    ValueError
        If the file is not a supported GFF format or provider cannot be detected.
    """
    path = Path(file)
    if not is_supported_gff(path):
        msg = f"Unsupported GFF format: {path.name!r}. Use GTF for FlyBase and WormBase."
        raise ValueError(msg)

    meta = get_gff_metadata(path)
    provider = meta.provider
    gff_format = meta.format.value  # "GTF" or "GFF3"

    import importlib

    mod = importlib.import_module(_PROVIDER_MODULES[provider])
    extract_fn = getattr(mod, f"extract_{level}")

    # Parse with feature-type pre-filter for efficiency.
    feature_filter = _LEVEL_FEATURE_FILTER.get(level)
    df = parse_gff(path, feature_filter=feature_filter)

    # Provider-specific feature extraction.
    df = extract_fn(df, gff_format)

    # Shared post-processing.
    df = standardize_column_names(df)

    # Handle identifier versions.
    names_col = _LEVEL_TO_NAMES_COL[level]
    version_col = f"{names_col}_version"
    df = handle_versions(
        df, id_col=names_col, version_col=version_col, ignore_version=ignore_version
    )

    df = remove_unwanted_biotypes(df, provider=str(provider))
    df = sort_ranges(df)
    df = drop_all_na_columns(df)

    if extra_mcols:
        df = add_broad_class(df)

    # Build metadata dict from GffMetadata.
    gr_metadata: dict = {
        "file": str(meta.file),
        "format": str(meta.format),
        "provider": str(meta.provider),
        "level": level,
        "ignore_version": ignore_version,
    }
    for key in ("organism", "genome_build", "release", "gff_version", "md5", "sha256"):
        val = getattr(meta, key, None)
        if val is not None:
            gr_metadata[key] = val

    # Use names_col if present; fall back to transcript_id or gene_id.
    if names_col not in df.columns:
        for fallback in ("transcript_id", "gene_id"):
            if fallback in df.columns:
                names_col = fallback
                break

    return dataframe_to_granges(df, names_col=names_col, metadata=gr_metadata)


__all__ = ["make_granges_from_gff"]
