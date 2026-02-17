"""URL download caching and HTTP fetch utilities.

Inspired by BiocFileCache from R/Bioconductor. Provides local file
caching so large downloads are not repeated unnecessarily.
"""

from __future__ import annotations

import hashlib
import logging
import os
from pathlib import Path
from typing import Any

import requests

logger = logging.getLogger(__name__)

_DEFAULT_CACHE_DIR = Path.home() / ".cache" / "acidgenomes"


def get_cache_dir() -> Path:
    """Return the cache directory, creating it if needed.

    Uses ``ACIDGENOMES_CACHE_DIR`` env var if set, otherwise
    ``~/.cache/acidgenomes``.
    """
    cache = Path(os.environ.get("ACIDGENOMES_CACHE_DIR", _DEFAULT_CACHE_DIR))
    cache.mkdir(parents=True, exist_ok=True)
    return cache


def cache_url(url: str, *, force: bool = False) -> Path:
    """Download a URL to the local cache and return the file path.

    Parameters
    ----------
    url : str
        Remote URL to download.
    force : bool
        Re-download even if the cached file already exists.

    Returns
    -------
    Path
        Path to the cached file on disk.
    """
    digest = hashlib.sha256(url.encode()).hexdigest()[:16]
    # Preserve file extension(s) for pandas readers.
    basename = url.rsplit("/", 1)[-1]
    dest = get_cache_dir() / f"{digest}_{basename}"
    if dest.exists() and not force:
        logger.debug("Using cached file: %s", dest)
        return dest
    logger.info("Downloading %s", url)
    resp = requests.get(url, timeout=120, stream=True)
    resp.raise_for_status()
    with open(dest, "wb") as fh:
        for chunk in resp.iter_content(chunk_size=1 << 16):
            fh.write(chunk)
    logger.debug("Cached to %s", dest)
    return dest


def fetch_json(url: str, **kwargs: Any) -> dict:
    """GET *url* and return the parsed JSON body."""
    resp = requests.get(url, timeout=30, **kwargs)
    resp.raise_for_status()
    return resp.json()


def fetch_text(url: str, **kwargs: Any) -> str:
    """GET *url* and return the response body as text."""
    resp = requests.get(url, timeout=30, **kwargs)
    resp.raise_for_status()
    return resp.text
