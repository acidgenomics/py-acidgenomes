"""Fetch current genome release version numbers.

Ported from ``currentGenomeVersion.R``.
"""

from __future__ import annotations

import logging
import re

from acidgenomes._cache import fetch_text

logger = logging.getLogger(__name__)


def current_ensembl_version() -> int:
    """Return the current Ensembl release version.

    Returns
    -------
    int

    Examples
    --------
    >>> current_ensembl_version()  # doctest: +SKIP
    112
    """
    url = "https://ftp.ensembl.org/pub/current_README"
    text = fetch_text(url)
    m = re.search(r"Ensembl Release\s+(\d+)", text, re.IGNORECASE)
    if m is None:
        raise RuntimeError("Failed to parse Ensembl release from README.")
    return int(m.group(1))


def current_gencode_version(organism: str = "Homo sapiens") -> str:
    """Return the current GENCODE release version.

    Parameters
    ----------
    organism : str

    Returns
    -------
    str
        Version string (e.g. ``'46'`` for human, ``'M35'`` for mouse).
    """
    short = "mouse" if organism == "Mus musculus" else "human"
    url = f"https://www.gencodegenes.org/{short}/"
    text = fetch_text(url)
    m = re.search(r"Release\s+(M?\d+)", text)
    if m is None:
        raise RuntimeError("Failed to parse GENCODE version.")
    return m.group(1)


def current_refseq_version() -> str:
    """Return the current RefSeq release version.

    Returns
    -------
    str
    """
    url = "https://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER"
    return fetch_text(url).strip()


def current_flybase_version(dmel: bool = True) -> str:
    """Return the current FlyBase release version.

    Parameters
    ----------
    dmel : bool
        If ``True``, returns the *Drosophila melanogaster* annotation version.
        Otherwise returns the FlyBase data release version.

    Returns
    -------
    str
    """
    url = "https://ftp.flybase.org/releases/current/RELEASE_NOTES"
    text = fetch_text(url)
    if dmel:
        m = re.search(r"Dmel\s+Release\s+([\w.]+)", text)
    else:
        m = re.search(r"Release\s+(FB\d{4}_\d+)", text)
    if m is None:
        raise RuntimeError("Failed to parse FlyBase version.")
    return m.group(1)


def current_wormbase_version() -> str:
    """Return the current WormBase release version.

    Returns
    -------
    str
    """
    url = "https://ftp.wormbase.org/pub/wormbase/releases/current-production-release"
    text = fetch_text(url)
    m = re.search(r"(WS\d+)", text)
    if m is None:
        raise RuntimeError("Failed to parse WormBase version.")
    return m.group(1)
