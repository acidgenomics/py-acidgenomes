"""Download UCSC genome files."""

from pathlib import Path

from acidgenomes.download._common import _download_file, _make_output_dir, _slugify


def download_ucsc_genome(
    organism: str,
    *,
    genome_build: str | None = None,
    output_dir: Path = Path.cwd(),
    cache: bool = False,
) -> dict[str, Path]:
    """Download UCSC genome FASTA and annotation files.

    Parameters
    ----------
    organism : str
        Latin organism name, e.g. ``"Homo sapiens"``.
    genome_build : str or None
        UCSC genome build (e.g. ``"hg38"``). Auto-detected if ``None``.
    output_dir : Path
        Parent directory for the download.
    cache : bool
        If ``True``, skip downloading files that already exist.

    Returns
    -------
    dict[str, Path]
        Mapping of file type to local path.
    """
    from acidgenomes._genome_build import current_ucsc_genome_build

    if genome_build is None:
        genome_build = current_ucsc_genome_build(organism)

    label = f"{_slugify(organism)}-{genome_build}-ucsc"
    out_dir = _make_output_dir(Path(output_dir), label)

    base = f"https://hgdownload.soe.ucsc.edu/goldenPath/{genome_build}"
    files = {
        "genome": (
            f"{base}/bigZips/{genome_build}.fa.gz",
            out_dir / "genome" / f"{genome_build}.fa.gz",
        ),
        "gtf_known_gene": (
            f"https://hgdownload.soe.ucsc.edu/goldenPath/{genome_build}/bigZips/genes/{genome_build}.knownGene.gtf.gz",
            out_dir / "annotation" / f"{genome_build}.knownGene.gtf.gz",
        ),
        "gtf_ncbi_refseq": (
            f"https://hgdownload.soe.ucsc.edu/goldenPath/{genome_build}/bigZips/genes/{genome_build}.ncbiRefSeq.gtf.gz",
            out_dir / "annotation" / f"{genome_build}.ncbiRefSeq.gtf.gz",
        ),
    }

    import requests

    result: dict[str, Path] = {}
    for key, (url, dest) in files.items():
        if cache and dest.exists():
            result[key] = dest
            continue
        # Suppress only HTTP 404 (file absent for this assembly); propagate real errors.
        try:
            result[key] = _download_file(url, dest)
        except requests.HTTPError as exc:
            if exc.response is not None and exc.response.status_code == 404:
                continue
            raise

    return result
