"""Download RefSeq genome files."""

from pathlib import Path

from acidgenomes.download._common import _download_file, _make_output_dir, _slugify


def download_refseq_genome(
    organism: str,
    *,
    taxonomic_group: str | None = None,
    genome_build: str | None = None,
    output_dir: Path | None = None,
    cache: bool = False,
) -> dict[str, Path]:
    """Download RefSeq genome GFF3 and GTF annotation files.

    Parameters
    ----------
    organism : str
        Latin organism name, e.g. ``"Homo sapiens"``.
    taxonomic_group : str or None
        NCBI taxonomic group (e.g. ``"vertebrate_mammalian"``).
        Auto-detected from internal mapping if ``None``.
    genome_build : str or None
        RefSeq genome build. Auto-detected if ``None``.
    output_dir : Path or None
        Parent directory for the download. Defaults to the current
        working directory at call time if ``None``.
    cache : bool
        If ``True``, skip downloading files that already exist.

    Returns
    -------
    dict[str, Path]
        Mapping of file type to local path.
    """
    from acidgenomes._data import NCBI_TAXONOMIC_GROUPS
    from acidgenomes._genome_build import current_refseq_genome_build

    if taxonomic_group is None:
        org_groups = NCBI_TAXONOMIC_GROUPS.get(organism, {})
        taxonomic_group = org_groups.get("refseq", "vertebrate_mammalian")
    # genome_build used for label only; accession is resolved from assembly_summary.
    if genome_build is None:
        genome_build = current_refseq_genome_build(organism)
    if output_dir is None:
        output_dir = Path.cwd()

    # RefSeq uses NCBI FTP paths. We look up the accession from NCBI assembly_summary.
    import requests

    summary_url = (
        f"https://ftp.ncbi.nlm.nih.gov/genomes/refseq/{taxonomic_group}/"
        f"{organism.replace(' ', '_')}/assembly_summary.txt"
    )
    resp = requests.get(summary_url, timeout=30)
    resp.raise_for_status()
    # Parse assembly_summary: find latest reference genome accession.
    accession = None
    ftp_path = None
    for line in resp.text.splitlines():
        if line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 20:
            continue
        if cols[10] == "reference genome" and cols[11] == "latest":
            accession = cols[0]
            ftp_path = cols[19]
            break
    if accession is None or ftp_path is None:
        msg = f"Could not find reference genome accession for {organism!r} in RefSeq."
        raise ValueError(msg)

    # ftp_path is like "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/.../name"
    base = ftp_path.replace("ftp://", "https://")
    dirname = base.rstrip("/").split("/")[-1]
    label = f"{_slugify(organism)}-{genome_build}-refseq"
    out_dir = _make_output_dir(Path(output_dir), label)

    files = {
        "gff3": (
            f"{base}/{dirname}_genomic.gff.gz",
            out_dir / "annotation" / f"{dirname}_genomic.gff.gz",
        ),
        "gtf": (
            f"{base}/{dirname}_genomic.gtf.gz",
            out_dir / "annotation" / f"{dirname}_genomic.gtf.gz",
        ),
        "genome": (
            f"{base}/{dirname}_genomic.fna.gz",
            out_dir / "genome" / f"{dirname}_genomic.fna.gz",
        ),
        "transcriptome": (
            f"{base}/{dirname}_rna.fna.gz",
            out_dir / "transcriptome" / f"{dirname}_rna.fna.gz",
        ),
    }

    import requests

    result: dict[str, Path] = {}
    for key, (url, dest) in files.items():
        if cache and dest.exists():
            result[key] = dest
            continue
        # Suppress only HTTP 404 (file absent for this organism); propagate real errors.
        try:
            result[key] = _download_file(url, dest)
        except requests.HTTPError as exc:
            if exc.response is not None and exc.response.status_code == 404:
                continue
            raise

    return result
