"""Core streaming GFF3/GTF file parser."""

import gzip
import re
from pathlib import Path
from typing import IO

import pandas as pd

# GFF3 attribute format: key=value;key=value
# GTF attribute format:  key "value"; key "value";
_GTF_ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')
_GFF3_ATTR_RE = re.compile(r"([^=;]+)=([^;]+)")

_GFF_COLS = ["seqnames", "source", "type", "start", "end", "score", "strand", "frame", "attributes"]


def _is_gzip(path: Path) -> bool:
    with path.open("rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def _open_gff(path: Path) -> IO[str]:
    if _is_gzip(path):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def _parse_gtf_attributes(attr_str: str) -> dict[str, str | list[str]]:
    """Parse GTF attribute string into a dict.

    Parameters
    ----------
    attr_str : str
        Raw GTF attributes column, e.g. ``gene_id "ENSG0001"; gene_name "TP53";``.

    Returns
    -------
    dict[str, str | list[str]]
        Parsed key-value pairs.
    """
    result: dict[str, str | list[str]] = {}
    for key, value in _GTF_ATTR_RE.findall(attr_str):
        if key in result:
            existing = result[key]
            if isinstance(existing, list):
                existing.append(value)
            else:
                result[key] = [existing, value]
        else:
            result[key] = value
    return result


def _parse_gff3_attributes(attr_str: str) -> dict[str, str | list[str]]:
    """Parse GFF3 attribute string into a dict.

    Parameters
    ----------
    attr_str : str
        Raw GFF3 attributes column, e.g. ``ID=gene1;Name=TP53;Dbxref=GeneID:7157``.

    Returns
    -------
    dict[str, str | list[str]]
        Parsed key-value pairs. Multi-valued attributes (comma-separated) are
        returned as lists.
    """
    result: dict[str, str | list[str]] = {}
    for raw_key, raw_value in _GFF3_ATTR_RE.findall(attr_str):
        key = raw_key.strip()
        values = [v.strip() for v in raw_value.split(",")]
        if len(values) == 1:
            result[key] = values[0]
        else:
            result[key] = values
    return result


def _detect_format(attr_str: str) -> str:
    """Detect GFF format from an attribute string.

    Parameters
    ----------
    attr_str : str
        A raw attributes column value.

    Returns
    -------
    str
        ``"GTF"`` or ``"GFF3"``.
    """
    if "=" in attr_str and '"' not in attr_str:
        return "GFF3"
    return "GTF"


def read_gff_directives(path: Path, *, max_lines: int = 200) -> dict[str, str]:
    """Read directive (comment) lines from the head of a GFF/GTF file.

    Matches lines of the form ``#!key value`` or ``##key: value``.

    Parameters
    ----------
    path : Path
        Path to the GFF/GTF file.
    max_lines : int
        Maximum number of lines to read.

    Returns
    -------
    dict[str, str]
        Mapping of directive key to value.
    """
    pattern = re.compile(r"^(#!|#+)([a-z-]+)(?::)?\s+(.+)$")
    directives: dict[str, str] = {}
    with _open_gff(path) as fh:
        for i, raw_line in enumerate(fh):
            if i >= max_lines:
                break
            line = raw_line.rstrip("\n")
            if not line.startswith("#"):
                break
            m = pattern.match(line)
            if m:
                key = m.group(2).strip()
                value = m.group(3).strip()
                directives[key] = value
    return directives


def parse_gff(
    path: Path,
    *,
    feature_filter: set[str] | None = None,
) -> pd.DataFrame:
    """Parse a GFF3 or GTF file into a DataFrame.

    Streams the file line by line, so memory usage is proportional to the
    number of rows that pass the feature filter, not the full file size.

    Parameters
    ----------
    path : Path
        Path to a GFF3 or GTF file (plain or gzip-compressed).
    feature_filter : set[str] or None
        If provided, only rows whose ``type`` column is in this set are kept.
        Pass ``{"gene"}`` to extract only gene features.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns from ``_GFF_COLS`` plus all flattened attribute
        key-value pairs. GFF coordinates are 1-based; this function does NOT
        convert them.
    """
    records: list[dict] = []
    attr_format: str | None = None

    with _open_gff(path) as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t", 8)
            if len(parts) < 9:
                continue
            feat_type = parts[2]
            if feature_filter is not None and feat_type not in feature_filter:
                continue
            attr_str = parts[8]
            if attr_format is None:
                attr_format = _detect_format(attr_str)
            if attr_format == "GTF":
                attrs = _parse_gtf_attributes(attr_str)
            else:
                attrs = _parse_gff3_attributes(attr_str)
            row: dict = {
                "seqnames": parts[0],
                "source": parts[1],
                "type": feat_type,
                "start": int(parts[3]),
                "end": int(parts[4]),
                "score": parts[5],
                "strand": parts[6],
                "frame": parts[7],
            }
            row.update(attrs)
            records.append(row)

    if not records:
        return pd.DataFrame(columns=_GFF_COLS)

    df = pd.DataFrame.from_records(records)
    # Ensure standard columns are present and in front.
    for col in _GFF_COLS[:-1]:  # exclude "attributes" — we expand instead
        if col not in df.columns:
            df[col] = pd.NA
    df["format"] = attr_format or "GTF"
    return df


def read_first_data_line(path: Path) -> dict[str, str] | None:
    """Read the first non-comment data line and return its split fields.

    Parameters
    ----------
    path : Path
        Path to the GFF/GTF file.

    Returns
    -------
    dict[str, str] or None
        Mapping of column name to value, or None if no data lines found.
    """
    with _open_gff(path) as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t", 8)
            if len(parts) >= 2:
                return {"seqnames": parts[0], "source": parts[1]}
    return None
