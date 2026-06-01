"""Download Ensembl genome files."""

import re
from pathlib import Path

from acidgenomes.download._common import _download_file, _make_output_dir, _slugify


def download_ensembl_genome(
    organism: str,
    *,
    genome_build: str | None = None,
    release: int | None = None,
    output_dir: Path = Path.cwd(),
    cache: bool = False,
) -> dict[str, Path]:
    """Download Ensembl genome FASTA, transcriptome FASTA, and GTF/GFF3.

    Parameters
    ----------
    organism : str
        Latin organism name, e.g. ``"Homo sapiens"``.
    genome_build : str or None
        Ensembl genome build. Auto-detected if ``None``.
    release : int or None
        Ensembl release number. Auto-detected if ``None``.
    output_dir : Path
        Parent directory for the download.
    cache : bool
        If ``True``, skip downloading files that already exist.

    Returns
    -------
    dict[str, Path]
        Mapping of file type (``"genome"``, ``"transcriptome"``, ``"gtf"``,
        ``"gff3"``) to the downloaded local path.
    """
    from acidgenomes._genome_build import current_ensembl_genome_build
    from acidgenomes._genome_version import current_ensembl_version

    if genome_build is None:
        genome_build = current_ensembl_genome_build(organism)
    if release is None:
        release = current_ensembl_version()

    build_clean = re.sub(r"\.p\d+$", "", genome_build)
    slug = organism.replace(" ", "_")
    slug_lower = slug.lower()
    label = f"{_slugify(organism)}-{build_clean}-ensembl-{release}"
    out_dir = _make_output_dir(Path(output_dir), label)

    base_ftp = f"https://ftp.ensembl.org/pub/release-{release}"
    files = {
        "genome": (
            f"{base_ftp}/fasta/{slug_lower}/dna/"
            f"{slug}.{build_clean}.dna.primary_assembly.fa.gz",
            out_dir / "genome" / f"{slug}.{build_clean}.dna.primary_assembly.fa.gz",
        ),
        "transcriptome": (
            f"{base_ftp}/fasta/{slug_lower}/cdna/"
            f"{slug}.{build_clean}.cdna.all.fa.gz",
            out_dir / "transcriptome" / f"{slug}.{build_clean}.cdna.all.fa.gz",
        ),
        "gtf": (
            f"{base_ftp}/gtf/{slug_lower}/{slug}.{build_clean}.{release}.gtf.gz",
            out_dir / "annotation" / f"{slug}.{build_clean}.{release}.gtf.gz",
        ),
        "gff3": (
            f"{base_ftp}/gff3/{slug_lower}/{slug}.{build_clean}.{release}.gff3.gz",
            out_dir / "annotation" / f"{slug}.{build_clean}.{release}.gff3.gz",
        ),
    }

    result: dict[str, Path] = {}
    for key, (url, dest) in files.items():
        if cache and dest.exists():
            result[key] = dest
            continue
        result[key] = _download_file(url, dest)

    return result
