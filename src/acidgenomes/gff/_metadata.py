"""GFF/GTF metadata extraction and provider auto-detection."""

import contextlib
import hashlib
import re
from pathlib import Path

from acidgenomes._detect import detect_organism
from acidgenomes.gff._parser import read_first_data_line, read_gff_directives
from acidgenomes.gff._providers import (
    GFF_FILENAME_PATTERNS,
    UCSC_SOURCES,
    GffFormat,
    GffMetadata,
    Provider,
    detect_provider_from_filename,
)


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _detect_format_from_extension(path: Path) -> GffFormat:
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    if name.endswith(".gtf"):
        return GffFormat.GTF
    return GffFormat.GFF3


def _detect_provider_from_directives(directives: dict[str, str]) -> Provider | None:
    """Detect provider from GFF directive keys/values.

    Logic ported from R's .gffProvider() in internal-gff.R.
    """
    provider_val = directives.get("provider", "")
    anno_source = directives.get("annotation-source", "")
    genebuild_version = directives.get("genebuild-version", "")

    if provider_val == "GENCODE":
        return Provider.GENCODE
    if anno_source.startswith("NCBI"):
        return Provider.REFSEQ
    # Ensembl has a cluster of specific directives.
    ensembl_keys = {
        "genebuild-last-updated",
        "genome-build",
        "genome-build-accession",
        "genome-version",
    }
    if ensembl_keys.issubset(directives.keys()):
        return Provider.ENSEMBL
    # Also accept a subset: genome-build + genome-build-accession + genome-date + genome-version
    ensembl_keys2 = {"genome-build", "genome-build-accession", "genome-date", "genome-version"}
    if ensembl_keys2.issubset(directives.keys()):
        return Provider.ENSEMBL
    if re.match(r"^WS[0-9]+$", genebuild_version):
        return Provider.WORMBASE
    return None


def _detect_provider_from_source_column(path: Path) -> Provider | None:
    """Detect provider by inspecting the source column of the first data line."""
    first = read_first_data_line(path)
    if first is None:
        return None
    source = first.get("source", "")
    if UCSC_SOURCES.match(source):
        return Provider.UCSC
    if source == "FlyBase":
        return Provider.FLYBASE
    if source == "WormBase":
        return Provider.WORMBASE
    return None


def _detect_genome_build_from_directives(
    directives: dict[str, str], provider: Provider
) -> str | None:
    """Extract genome build from directive key-value pairs."""
    if provider == Provider.GENCODE:
        description = directives.get("description", "")
        # GRCh37 liftover: "mapped to GRCh37"
        m = re.search(r"mapped to ([^ ]+)", description)
        if m:
            return m.group(1)
        # Standard: "GENCODE human comprehensive gene annotation on the ... genome (GRCh38)"
        m = re.search(r"genome \(([^)]+)\)", description)
        if m:
            return m.group(1)
    genome_build = directives.get("genome-build", "")
    if genome_build:
        # Ensembl: "genome-build GRCh38 109" → take last token
        parts = genome_build.split()
        return parts[-1] if parts else None
    return None


def _extract_metadata_from_filename(path: Path, provider: Provider) -> dict:
    """Extract organism, genome_build, release from filename.

    Returns a partial dict with whatever can be determined.
    """
    name = path.name
    pattern = GFF_FILENAME_PATTERNS.get(provider)
    if pattern is None:
        return {}
    m = pattern.match(name)
    if m is None:
        return {}
    groups = m.groups()
    result: dict = {}
    if provider == Provider.ENSEMBL:
        # groups: (organism, genome_build, release, suffix?, format, gz?)
        result["organism"] = groups[0].replace("_", " ")
        result["genome_build"] = groups[1]
        result["release"] = int(groups[2])
    elif provider == Provider.FLYBASE:
        # groups: (organism_abbrev, subset, release, format, gz?)
        if groups[0] == "dmel":
            result["organism"] = "Drosophila melanogaster"
        result["release"] = groups[2]  # e.g. "r6.37"
        result["genome_build"] = result["release"]
    elif provider == Provider.GENCODE:
        # groups: (release, lift37?, format, gz?)
        release_str = groups[0]
        result["release"] = int(release_str) if release_str.isdigit() else release_str
    elif provider == Provider.UCSC:
        # groups: (genome_build, source_type, gz?)
        result["genome_build"] = groups[0]
    elif provider == Provider.WORMBASE:
        # groups: (organism, project, release, gene_set, format, gz?)
        if groups[0] == "c_elegans":
            result["organism"] = "Caenorhabditis elegans"
        result["release"] = groups[2]  # e.g. "WS279"
        result["genome_build"] = result["release"]
    return result


def get_gff_metadata(path: Path) -> GffMetadata:
    """Extract metadata from a GFF/GTF file.

    Reads directives, detects provider, genome build, organism, and release.
    Detection order: directives → filename regex → source column.

    Parameters
    ----------
    path : Path
        Path to a GFF3 or GTF file.

    Returns
    -------
    GffMetadata
        Fully-populated metadata object.

    Raises
    ------
    ValueError
        If provider cannot be detected.
    """
    path = Path(path)
    directives = read_gff_directives(path)

    gff_format = _detect_format_from_extension(path)
    # GENCODE GRCh37 GTF incorrectly advertises format as gff3 in directives.
    format_from_directive = directives.get("format", "").upper()
    if format_from_directive == "GTF":
        gff_format = GffFormat.GTF
    elif format_from_directive == "GFF3" and (
        directives.get("provider", "") != "GENCODE" or gff_format == GffFormat.GFF3
    ):
        # Only trust if it's not the known GENCODE GRCh37 bug.
        gff_format = GffFormat.GFF3

    gff_version = directives.get("gff-version")

    # Provider detection: directives → filename → source column
    provider = _detect_provider_from_directives(directives)
    if provider is None:
        provider = detect_provider_from_filename(path)
    if provider is None:
        provider = _detect_provider_from_source_column(path)
    if provider is None and path.name == "ref-transcripts.gtf":
        provider = Provider.UCSC
    if provider is None:
        msg = f"Failed to detect provider from {path.name!r}."
        raise ValueError(msg)

    genome_build = _detect_genome_build_from_directives(directives, provider)

    # Extract additional metadata from filename
    filename_meta = _extract_metadata_from_filename(path, provider)
    organism = filename_meta.get("organism")
    release = filename_meta.get("release")
    if genome_build is None:
        genome_build = filename_meta.get("genome_build")

    # RefSeq: parse release from annotation-source directive
    if provider == Provider.REFSEQ and release is None:
        anno_source = directives.get("annotation-source", "")
        m = re.match(r"^NCBI.+Annotation\sRelease\s([.0-9]+)$", anno_source)
        if m:
            release = m.group(1)

    # WormBase: parse release from genebuild-version directive
    if provider == Provider.WORMBASE and release is None:
        release = directives.get("genebuild-version")
        if genome_build is None and release:
            genome_build = release

    # Organism fallback: try to detect from genome build string
    if organism is None and genome_build:
        with contextlib.suppress(Exception):
            organism = detect_organism([genome_build])

    return GffMetadata(
        file=path,
        format=gff_format,
        provider=provider,
        organism=organism,
        genome_build=genome_build,
        release=release,
        gff_version=gff_version,
        directives=directives,
        md5=_md5(path),
        sha256=_sha256(path),
    )
