"""Organism detection from gene identifiers.

Ported from ``detectOrganism.R``.
"""

from __future__ import annotations

import logging

from acidgenomes._data import DETECT_ORGANISM_DATA

logger = logging.getLogger(__name__)

_MAX_SCAN = 50  # Check at most this many identifiers.


def _detect_single(string: str) -> str | None:
    """Try to match a single string against known organism patterns.

    Returns the organism name or ``None``.
    """
    for entry in DETECT_ORGANISM_DATA:
        if entry["pattern"].search(string):
            return entry["organism"]
    return None


def detect_organism(identifiers: list[str]) -> str:
    """Detect the organism from a list of gene/transcript identifiers.

    Parameters
    ----------
    identifiers : list[str]
        Gene or transcript identifiers (e.g. Ensembl IDs).

    Returns
    -------
    str
        Latin organism name.

    Raises
    ------
    ValueError
        If the organism cannot be detected or identifiers map to
        multiple organisms.

    Examples
    --------
    >>> detect_organism(["ENSG00000000003", "ENSG00000000005"])
    'Homo sapiens'
    """
    if not identifiers:
        raise ValueError("identifiers must not be empty.")
    detected: set[str] = set()
    for ident in identifiers[:_MAX_SCAN]:
        org = _detect_single(str(ident))
        if org is not None:
            detected.add(org)
    if len(detected) == 0:
        raise ValueError(
            f"Failed to detect organism from identifiers. First identifier: '{identifiers[0]}'."
        )
    if len(detected) > 1:
        raise ValueError(
            f"Multiple organisms detected: {', '.join(sorted(detected))}. "
            "Input must contain identifiers from a single organism."
        )
    return detected.pop()
