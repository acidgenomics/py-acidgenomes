"""Parse FASTA headers to build transcript-to-gene mappings.

Supports Ensembl, GENCODE, FlyBase, and WormBase FASTA header formats.
"""

import gzip
import re
from pathlib import Path

import pandas as pd

from acidgenomes._classes import TxToGene
from acidgenomes._strip_versions import strip_gene_versions, strip_transcript_versions

# Ensembl: ">ENST00000000233.10 cdna chromosome:GRCh38:... gene:ENSG00000000003.15 ..."
_ENSEMBL_HEADER = re.compile(r"^>(\S+)\s.*\bgene:(\S+)")

# GENCODE: ">ENST00000000233.10|ENSG00000000003.15|...|gene_symbol|..."
_GENCODE_HEADER = re.compile(r"^>(\S+)\|(\S+)\|")

# FlyBase: ">FBtr0300690 type=mRNA; ... parent=FBgn0031208,..."
_FLYBASE_HEADER = re.compile(r"^>(\S+)\s.*\bparent=([^,;\s]+)")

# WormBase: ">Y74C9A.2a.1 gene=WBGene00022276 ..."
_WORMBASE_HEADER = re.compile(r"^>(\S+)\s.*\bgene=(\S+)")


def make_tx_to_gene_from_fasta(
    file: str | Path,
    *,
    ignore_version: bool = False,
) -> TxToGene:
    """Build a TxToGene mapping by parsing FASTA transcript headers.

    Detects the FASTA format (Ensembl, GENCODE, FlyBase, WormBase) from
    the header structure of the first FASTA record.

    Parameters
    ----------
    file : str or Path
        Path to a FASTA file (plain or gzip-compressed). Only header lines
        are read — sequences are skipped.
    ignore_version : bool
        If ``True``, strip version suffixes from transcript and gene IDs.

    Returns
    -------
    TxToGene
        DataFrame with columns ``tx_id`` and ``gene_id``.
    """
    path = Path(file)
    opener = gzip.open if path.suffix == ".gz" or _is_gzip(path) else open
    records: list[dict[str, str]] = []
    pattern: re.Pattern | None = None

    with opener(path, "rt", encoding="utf-8") as fh:  # type: ignore[call-overload]
        for raw_line in fh:
            if not raw_line.startswith(">"):
                continue
            line = raw_line.rstrip("\n")
            if pattern is None:
                pattern = _detect_pattern(line)
                if pattern is None:
                    continue
            m = pattern.match(line)
            if m:
                tx_id = m.group(1)
                gene_id = m.group(2)
                records.append({"tx_id": tx_id, "gene_id": gene_id})

    if not records:
        msg = f"No transcript-gene pairs found in {path.name!r}."
        raise ValueError(msg)

    df = pd.DataFrame(records).drop_duplicates()

    if ignore_version:
        df["tx_id"] = strip_transcript_versions(df["tx_id"].tolist())
        df["gene_id"] = strip_gene_versions(df["gene_id"].tolist())
        df = df.drop_duplicates()

    return TxToGene(data=df)


def _is_gzip(path: Path) -> bool:
    with path.open("rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def _detect_pattern(header: str) -> re.Pattern | None:
    """Detect the FASTA header format from a single header line."""
    if "|" in header and _GENCODE_HEADER.match(header):
        return _GENCODE_HEADER
    if "gene:" in header and _ENSEMBL_HEADER.match(header):
        return _ENSEMBL_HEADER
    if "parent=" in header and _FLYBASE_HEADER.match(header):
        return _FLYBASE_HEADER
    if "gene=" in header and _WORMBASE_HEADER.match(header):
        return _WORMBASE_HEADER
    return None
