"""GENCODE release history and Ensembl-to-GENCODE release mapping."""

import re

import pandas as pd

from acidgenomes._cache import fetch_text

_GENCODE_HISTORY_URL_HUMAN = "https://www.gencodegenes.org/human/releases.html"
_GENCODE_HISTORY_URL_MOUSE = "https://www.gencodegenes.org/mouse/releases.html"


def gencode_release_history(organism: str = "Homo sapiens") -> pd.DataFrame:
    """Fetch the full GENCODE release history table.

    Scrapes the GENCODE website for a table of all releases and their
    corresponding Ensembl release and genome build.

    Parameters
    ----------
    organism : str
        ``"Homo sapiens"`` (default) or ``"Mus musculus"``.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ``gencode_release``, ``ensembl_release``,
        ``genome_build``.
    """
    if "musculus" in organism.lower():
        url = _GENCODE_HISTORY_URL_MOUSE
    else:
        url = _GENCODE_HISTORY_URL_HUMAN
    html = fetch_text(url)

    # Parse table rows: look for release number, Ensembl version, genome build.
    # Pattern: rows like "47", "GRCh38", "Ensembl 113"
    row_pattern = re.compile(
        r"<tr[^>]*>.*?<td[^>]*>([\dM.]+)</td>.*?<td[^>]*>([^<]+)</td>.*?<td[^>]*>[^<]*(?:Ensembl|ensembl)[^<]*?(\d+)[^<]*</td>",
        re.DOTALL,
    )
    records: list[dict] = []
    for m in row_pattern.finditer(html):
        records.append({
            "gencode_release": m.group(1).strip(),
            "genome_build": m.group(2).strip(),
            "ensembl_release": int(m.group(3).strip()),
        })

    return pd.DataFrame(records)


def map_ensembl_to_gencode(release: int, organism: str = "Homo sapiens") -> str | None:
    """Map an Ensembl release number to the corresponding GENCODE release.

    This is the inverse of :func:`~acidgenomes._mapping.map_gencode_to_ensembl`.

    Parameters
    ----------
    release : int
        Ensembl release number.
    organism : str
        ``"Homo sapiens"`` or ``"Mus musculus"``.

    Returns
    -------
    str or None
        GENCODE release string (e.g. ``"44"`` or ``"M33"``), or ``None`` if
        not found.
    """
    history = gencode_release_history(organism)
    if history.empty:
        return None
    match = history[history["ensembl_release"] == release]
    if match.empty:
        return None
    return str(match.iloc[0]["gencode_release"])
