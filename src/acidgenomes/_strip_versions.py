"""Strip version suffixes from Ensembl identifiers.

Ported from ``stripGeneVersions.R``, ``stripTranscriptVersions.R``,
``stripExonVersions.R``.
"""

from __future__ import annotations

import re

# Patterns that match Ensembl-style versioned identifiers.
_ENSEMBL_VERSION_RE = re.compile(r"^(ENS[A-Z]*\d+)\.\d+$")
_FLYBASE_VERSION_RE = re.compile(r"^(FBgn\d+)\.\d+$")
_WORMBASE_VERSION_RE = re.compile(r"^(WBGene\d+)\.\d+$")

_VERSION_PATTERNS = [
    _ENSEMBL_VERSION_RE,
    _FLYBASE_VERSION_RE,
    _WORMBASE_VERSION_RE,
]


def _strip(identifiers: list[str]) -> list[str]:
    """Strip version suffixes from identifiers."""
    result: list[str] = []
    for ident in identifiers:
        stripped = ident
        for pattern in _VERSION_PATTERNS:
            m = pattern.match(ident)
            if m:
                stripped = m.group(1)
                break
        result.append(stripped)
    return result


def strip_gene_versions(identifiers: list[str]) -> list[str]:
    """Strip version suffixes from gene identifiers.

    Parameters
    ----------
    identifiers : list[str]

    Returns
    -------
    list[str]

    Examples
    --------
    >>> strip_gene_versions(["ENSG00000000003.14"])
    ['ENSG00000000003']
    """
    return _strip(identifiers)


def strip_transcript_versions(identifiers: list[str]) -> list[str]:
    """Strip version suffixes from transcript identifiers.

    Parameters
    ----------
    identifiers : list[str]

    Returns
    -------
    list[str]

    Examples
    --------
    >>> strip_transcript_versions(["ENST00000000233.10"])
    ['ENST00000000233']
    """
    return _strip(identifiers)


def strip_exon_versions(identifiers: list[str]) -> list[str]:
    """Strip version suffixes from exon identifiers.

    Parameters
    ----------
    identifiers : list[str]

    Returns
    -------
    list[str]

    Examples
    --------
    >>> strip_exon_versions(["ENSE00000000001.2"])
    ['ENSE00000000001']
    """
    return _strip(identifiers)
