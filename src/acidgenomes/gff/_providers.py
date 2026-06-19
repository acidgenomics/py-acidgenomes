"""GFF provider definitions and filename pattern matching."""

import re
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path


class Provider(StrEnum):
    """Supported GFF/GTF annotation providers."""

    ENSEMBL = "Ensembl"
    FLYBASE = "FlyBase"
    GENCODE = "GENCODE"
    REFSEQ = "RefSeq"
    UCSC = "UCSC"
    WORMBASE = "WormBase"


class GffFormat(StrEnum):
    """GFF file format."""

    GFF3 = "GFF3"
    GTF = "GTF"


# Filename patterns ported from R's .gffPatterns in AllGlobals.R.
# Each pattern captures provider-specific groups used for metadata extraction.
GFF_FILENAME_PATTERNS: dict[Provider, re.Pattern[str]] = {
    Provider.ENSEMBL: re.compile(
        # e.g. "Homo_sapiens.GRCh38.102.gtf.gz"
        r"^(?:[a-z0-9]+_)?"  # optional BiocFileCache prefix
        r"([A-Z][a-z]+_[a-z]+)"  # organism: "Homo_sapiens"
        r"\.([A-Za-z0-9]+)"  # genome build: "GRCh38"
        r"\.([0-9]+)"  # release: "102"
        r"(?:\.chr_patch_hapl_scaff)?"
        r"\.(gff3|gtf)"
        r"(?:\.gz)?$"
    ),
    Provider.FLYBASE: re.compile(
        # e.g. "dmel-all-r6.37.gtf.gz"
        r"^(?:[a-z0-9]+_)?"
        r"([^-]+)"  # organism abbreviation: "dmel"
        r"-([^-]+)"  # subset: "all"
        r"-(r[0-9]+\.[0-9]+)"  # release: "r6.37"
        r"\.(gff|gtf)"
        r"(?:\.gz)?$"
    ),
    Provider.GENCODE: re.compile(
        # e.g. "gencode.v36.annotation.gtf.gz"
        r"^(?:[a-z0-9]+_)?"
        r"gencode"
        r"\.v([M0-9]+)"  # release: "36" (human) or "M25" (mouse)
        r"(lift37)?"  # optional GRCh37 liftover marker
        r"\.annotation"
        r"\.(gff3|gtf)"
        r"(?:\.gz)?$"
    ),
    Provider.REFSEQ: re.compile(
        # e.g. "GCF_000001405.38_GRCh38.p12_genomic.gff.gz"
        r"^(?:[0-9a-z]+_)?"
        r"(GC[AF]_[0-9]+\.[0-9]+)"  # accession: "GCF_000001405.38"
        r"_([^_]+)"  # genome build: "GRCh38.p12"
        r"_(.+)"  # suffix: "genomic" or "full_analysis_set"
        r"\.(gff|gtf)"
        r"(?:\.gz)?$"
    ),
    Provider.UCSC: re.compile(
        # e.g. "hg38.knownGene.gtf.gz"
        r"^(?:[0-9a-z]+_)?"
        r"([a-z]+[A-Za-z]+[0-9]+)"  # genome build: "hg38"
        r"\.(ensGene|knownGene|ncbiRefSeq|refGene)"
        r"\.gtf"
        r"(?:\.gz)?$"
    ),
    Provider.WORMBASE: re.compile(
        # e.g. "c_elegans.PRJNA13758.WS279.canonical_geneset.gtf.gz"
        r"^(?:[a-z0-9]+_)?"
        r"([a-z]_[a-z]+)"  # organism: "c_elegans"
        r"\.([A-Z0-9]+)"  # project: "PRJNA13758"
        r"\.(WS[0-9]+)"  # release: "WS279"
        r"\.([a-z_]+)"  # gene set: "canonical_geneset"
        r"\.(gff3|gtf)"
        r"(?:\.gz)?$"
    ),
}

# Denylist: FlyBase GFF3 and WormBase GFF3 are not supported; use GTF.
_DENYLIST_PATTERNS: list[re.Pattern[str]] = [
    re.compile(r"^(?:[a-z0-9]+_)?([^-]+)-([^-]+)-(r[0-9]+\.[0-9]+)\.gff(?:\.gz)?$"),
    re.compile(
        r"^(?:[a-z0-9]+_)?([a-z]_[a-z]+)\.([A-Z0-9]+)\.(WS[0-9]+)\.([a-z_]+)\.gff3(?:\.gz)?$"
    ),
]

# UCSC source column values that identify UCSC-origin files.
UCSC_SOURCES: re.Pattern[str] = re.compile(
    r"^(ensGene|knownGene|ncbiRefSeq|ncbiRefSeq\.\d{4}-\d{2}-\d{2}|refGene)$"
)


@dataclass
class GffMetadata:
    """Metadata extracted from a GFF/GTF file."""

    file: Path
    format: GffFormat
    provider: Provider
    organism: str | None = None
    genome_build: str | None = None
    release: str | int | None = None
    gff_version: str | None = None
    directives: dict[str, str] = field(default_factory=dict)
    md5: str | None = None
    sha256: str | None = None


def detect_provider_from_filename(path: Path) -> Provider | None:
    """Detect GFF provider from the filename using regex patterns.

    Parameters
    ----------
    path : Path
        Path to the GFF/GTF file.

    Returns
    -------
    Provider or None
        Detected provider, or None if no pattern matches.
    """
    name = path.name
    for provider, pattern in GFF_FILENAME_PATTERNS.items():
        if pattern.match(name):
            return provider
    return None


def is_supported_gff(path: Path) -> bool:
    """Check if the GFF file is supported (not on the denylist).

    Parameters
    ----------
    path : Path
        Path to the GFF/GTF file.

    Returns
    -------
    bool
        False if the file is a FlyBase GFF3 or WormBase GFF3 (unsupported).
    """
    name = path.name
    return all(not pattern.match(name) for pattern in _DENYLIST_PATTERNS)
