"""Download GENCODE genome files."""

from pathlib import Path

from acidgenomes.download._common import _download_file, _make_output_dir, _slugify


def download_gencode_genome(
    organism: str,
    *,
    genome_build: str | None = None,
    release: str | None = None,
    output_dir: Path = Path.cwd(),
    cache: bool = False,
) -> dict[str, Path]:
    """Download GENCODE genome FASTA, transcriptome FASTA, GTF, and GFF3.

    Parameters
    ----------
    organism : str
        ``"Homo sapiens"`` or ``"Mus musculus"``.
    genome_build : str or None
        Genome build (e.g. ``"GRCh38"``). Auto-detected if ``None``.
    release : str or None
        GENCODE release (e.g. ``"44"`` for human, ``"M33"`` for mouse).
        Auto-detected if ``None``.
    output_dir : Path
        Parent directory for the download.
    cache : bool
        If ``True``, skip downloading files that already exist.

    Returns
    -------
    dict[str, Path]
        Mapping of file type to local path.
    """
    from acidgenomes._genome_build import current_gencode_genome_build
    from acidgenomes._genome_version import current_gencode_version

    if genome_build is None:
        genome_build = current_gencode_genome_build(organism)
    if release is None:
        release = current_gencode_version(organism)

    is_mouse = organism.lower().startswith("mus")
    if is_mouse and not str(release).startswith("M"):
        release_prefix = f"M{release}"
    else:
        release_prefix = str(release)
    base_url = "https://ftp.ebi.ac.uk/pub/databases/gencode"
    org_path = "Gencode_mouse" if is_mouse else "Gencode_human"
    release_path = f"release_{release_prefix}"

    gtf_name = f"gencode.v{release_prefix}.annotation.gtf.gz"
    gff3_name = f"gencode.v{release_prefix}.annotation.gff3.gz"
    # Use the detected/provided genome_build (not a hardcoded version).
    # e.g. GRCm39 for modern mouse, GRCh38 for human.
    genome_name = f"{genome_build}.primary_assembly.genome.fa.gz"
    tx_name = f"gencode.v{release_prefix}.transcripts.fa.gz"

    label = f"{_slugify(organism)}-{genome_build}-gencode-{release_prefix}"
    out_dir = _make_output_dir(Path(output_dir), label)

    files = {
        "genome": (
            f"{base_url}/{org_path}/{release_path}/{genome_name}",
            out_dir / "genome" / genome_name,
        ),
        "transcriptome": (
            f"{base_url}/{org_path}/{release_path}/{tx_name}",
            out_dir / "transcriptome" / tx_name,
        ),
        "gtf": (
            f"{base_url}/{org_path}/{release_path}/{gtf_name}",
            out_dir / "annotation" / gtf_name,
        ),
        "gff3": (
            f"{base_url}/{org_path}/{release_path}/{gff3_name}",
            out_dir / "annotation" / gff3_name,
        ),
    }

    result: dict[str, Path] = {}
    for key, (url, dest) in files.items():
        if cache and dest.exists():
            result[key] = dest
            continue
        result[key] = _download_file(url, dest)

    return result
