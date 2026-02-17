"""Fetch current genome builds from provider APIs.

Ported from ``currentGenomeBuild.R``.
"""

from __future__ import annotations

import logging

from acidgenomes._cache import fetch_json, fetch_text

logger = logging.getLogger(__name__)


def current_ensembl_genome_build(organism: str = "Homo sapiens") -> str:
    """Return the current Ensembl genome build for *organism*.

    Parameters
    ----------
    organism : str
        Latin organism name.

    Returns
    -------
    str
        Genome build name (e.g. ``'GRCh38'``).

    Examples
    --------
    >>> current_ensembl_genome_build("Homo sapiens")  # doctest: +SKIP
    'GRCh38'
    """
    species = organism.lower().replace(" ", "_")
    url = f"https://rest.ensembl.org/info/assembly/{species}?content-type=application/json"
    data = fetch_json(url)
    build = data.get("assembly_name") or data.get("default_coord_system_version")
    if build is None:
        raise RuntimeError(f"Failed to detect genome build for '{organism}'.")
    return str(build)


def current_gencode_genome_build(organism: str = "Homo sapiens") -> str:
    """Return the current GENCODE genome build.

    Delegates to Ensembl since GENCODE is built on Ensembl annotations.

    Parameters
    ----------
    organism : str

    Returns
    -------
    str
    """
    return current_ensembl_genome_build(organism)


def current_refseq_genome_build(organism: str = "Homo sapiens") -> str:
    """Return the current RefSeq genome build for *organism*.

    Parses the NCBI assembly summary report.

    Parameters
    ----------
    organism : str

    Returns
    -------
    str
    """
    url = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/assembly_summary.txt"
    text = fetch_text(url)
    for line in text.splitlines():
        if line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 16:
            continue
        org_name = parts[7]  # organism_name
        if org_name.lower() == organism.lower() and parts[4] == "reference genome":
            return parts[15]  # asm_name
    raise RuntimeError(f"Failed to detect RefSeq genome build for '{organism}'.")


def current_ucsc_genome_build(organism: str = "Homo sapiens") -> str:
    """Return the current UCSC genome build for *organism*.

    Parameters
    ----------
    organism : str

    Returns
    -------
    str
    """
    organism.lower().replace(" ", "_")
    url = "https://api.genome.ucsc.edu/list/ucscGenomes"
    data = fetch_json(url)
    genomes = data.get("ucscGenomes", {})
    # Find the most recent genome for this organism.
    matches = []
    for db_name, info in genomes.items():
        sci_name = info.get("scientificName", "")
        if sci_name.lower() == organism.lower():
            order_key = info.get("orderKey", 0)
            matches.append((order_key, db_name))
    if not matches:
        raise RuntimeError(f"No UCSC genome found for '{organism}'.")
    matches.sort(reverse=True)
    return matches[0][1]
