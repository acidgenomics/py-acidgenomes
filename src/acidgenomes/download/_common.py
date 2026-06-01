"""Shared utilities for genome download functions."""

import hashlib
import re
from pathlib import Path

import requests


def _make_output_dir(base: Path, label: str) -> Path:
    """Create and return a structured output directory."""
    out = base / label
    out.mkdir(parents=True, exist_ok=True)
    (out / "genome").mkdir(exist_ok=True)
    (out / "transcriptome").mkdir(exist_ok=True)
    (out / "annotation").mkdir(exist_ok=True)
    return out


def _slugify(s: str) -> str:
    return re.sub(r"\s+", "_", s).lower()


def _download_file(url: str, dest: Path, *, chunk_size: int = 65536) -> Path:
    """Stream a URL to a local file."""
    with requests.get(url, stream=True, timeout=120) as resp:
        resp.raise_for_status()
        with dest.open("wb") as fh:
            for chunk in resp.iter_content(chunk_size=chunk_size):
                fh.write(chunk)
    return dest


def _md5_file(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()
