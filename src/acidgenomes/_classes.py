"""S4 class equivalents as Python dataclasses.

The R package defines 28 S4 classes. Here they are represented as
``@dataclass`` wrappers around a ``pandas.DataFrame`` with attached
metadata dictionary.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import pandas as pd

# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------


def _check_required_columns(df: pd.DataFrame, cols: list[str], cls: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{cls} requires column(s): {', '.join(missing)}")


def _check_no_duplicates(df: pd.DataFrame, col: str, cls: str) -> None:
    if df[col].duplicated().any():
        n = df[col].duplicated().sum()
        raise ValueError(f"{cls}: {n} duplicate(s) in column '{col}'.")


def _check_has_rows(df: pd.DataFrame, cls: str) -> None:
    if len(df) == 0:
        raise ValueError(f"{cls}: DataFrame has zero rows.")


def _check_no_na(df: pd.DataFrame, col: str, cls: str) -> None:
    n = df[col].isna().sum()
    if n:
        raise ValueError(f"{cls}: {n} NA value(s) in column '{col}'.")


# ---------------------------------------------------------------------------
# Base annotated DataFrame
# ---------------------------------------------------------------------------


@dataclass
class _AnnotatedDataFrame:
    """A DataFrame with attached metadata dictionary."""

    data: pd.DataFrame
    metadata: dict[str, Any] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.data)

    def __repr__(self) -> str:
        name = type(self).__name__
        nrow, ncol = self.data.shape
        return f"{name}: {nrow} x {ncol}"

    def summary(self) -> str:
        parts = [repr(self)]
        for k, v in self.metadata.items():
            parts.append(f"  {k}: {v}")
        return "\n".join(parts)


# ---------------------------------------------------------------------------
# GRanges-like base (has organism/provider/genome_build metadata)
# ---------------------------------------------------------------------------


@dataclass
class _GRangesBase(_AnnotatedDataFrame):
    """Base for classes that represent genomic ranges."""

    @property
    def organism(self) -> str | None:
        return self.metadata.get("organism")

    @property
    def provider(self) -> str | None:
        return self.metadata.get("provider")

    @property
    def genome_build(self) -> str | None:
        return self.metadata.get("genome_build")


# ---------------------------------------------------------------------------
# Provider-specific bases (factory)
# ---------------------------------------------------------------------------


def _make_provider_base(provider: str) -> type:
    """Create a provider-specific base class."""

    @dataclass
    class _ProviderBase(_GRangesBase):
        def __post_init__(self) -> None:
            self.metadata.setdefault("provider", provider)

    _ProviderBase.__name__ = f"_{provider}Base"
    _ProviderBase.__qualname__ = f"_{provider}Base"
    return _ProviderBase


_EnsemblBase = _make_provider_base("Ensembl")
_FlybaseBase = _make_provider_base("FlyBase")
_GencodeBase = _make_provider_base("GENCODE")
_RefseqBase = _make_provider_base("RefSeq")
_UcscBase = _make_provider_base("UCSC")
_WormbaseBase = _make_provider_base("WormBase")


# ---------------------------------------------------------------------------
# Ensembl
# ---------------------------------------------------------------------------


@dataclass
class EnsemblGenes(_EnsemblBase):
    """Ensembl gene annotations."""

    pass


@dataclass
class EnsemblTranscripts(_EnsemblBase):
    """Ensembl transcript annotations."""

    pass


@dataclass
class EnsemblExons(_EnsemblBase):
    """Ensembl exon annotations."""

    pass


# ---------------------------------------------------------------------------
# FlyBase
# ---------------------------------------------------------------------------


@dataclass
class FlybaseGenes(_FlybaseBase):
    """FlyBase gene annotations."""

    pass


@dataclass
class FlybaseTranscripts(_FlybaseBase):
    """FlyBase transcript annotations."""

    pass


@dataclass
class FlybaseExons(_FlybaseBase):
    """FlyBase exon annotations."""

    pass


# ---------------------------------------------------------------------------
# GENCODE
# ---------------------------------------------------------------------------


@dataclass
class GencodeGenes(_GencodeBase):
    """GENCODE gene annotations."""

    pass


@dataclass
class GencodeTranscripts(_GencodeBase):
    """GENCODE transcript annotations."""

    pass


@dataclass
class GencodeExons(_GencodeBase):
    """GENCODE exon annotations."""

    pass


# ---------------------------------------------------------------------------
# RefSeq
# ---------------------------------------------------------------------------


@dataclass
class RefseqGenes(_RefseqBase):
    """RefSeq gene annotations."""

    pass


@dataclass
class RefseqTranscripts(_RefseqBase):
    """RefSeq transcript annotations."""

    pass


@dataclass
class RefseqExons(_RefseqBase):
    """RefSeq exon annotations."""

    pass


# ---------------------------------------------------------------------------
# UCSC
# ---------------------------------------------------------------------------


@dataclass
class UcscGenes(_UcscBase):
    """UCSC gene annotations."""

    pass


@dataclass
class UcscTranscripts(_UcscBase):
    """UCSC transcript annotations."""

    pass


@dataclass
class UcscExons(_UcscBase):
    """UCSC exon annotations."""

    pass


# ---------------------------------------------------------------------------
# WormBase
# ---------------------------------------------------------------------------


@dataclass
class WormbaseGenes(_WormbaseBase):
    """WormBase gene annotations."""

    pass


@dataclass
class WormbaseTranscripts(_WormbaseBase):
    """WormBase transcript annotations."""

    pass


@dataclass
class WormbaseExons(_WormbaseBase):
    """WormBase exon annotations."""

    pass


# ---------------------------------------------------------------------------
# Gene nomenclature / mapping classes
# ---------------------------------------------------------------------------


@dataclass
class Hgnc(_AnnotatedDataFrame):
    """HGNC (Human Gene Nomenclature Committee) data."""

    pass


@dataclass
class Mgi(_AnnotatedDataFrame):
    """MGI (Mouse Genome Informatics) data."""

    pass


@dataclass
class NcbiGeneInfo(_AnnotatedDataFrame):
    """NCBI gene information."""

    pass


@dataclass
class NcbiGeneHistory(_AnnotatedDataFrame):
    """NCBI gene history."""

    pass


@dataclass
class EnsemblToNcbi(_AnnotatedDataFrame):
    """Ensembl-to-NCBI gene identifier mapping."""

    pass


@dataclass
class NcbiToEnsembl(_AnnotatedDataFrame):
    """NCBI-to-Ensembl gene identifier mapping."""

    pass


@dataclass
class GeneToSymbol(_AnnotatedDataFrame):
    """Gene identifier-to-symbol mapping."""

    pass


@dataclass
class TxToGene(_AnnotatedDataFrame):
    """Transcript-to-gene identifier mapping."""

    pass


@dataclass
class JaxHumanToMouse(_AnnotatedDataFrame):
    """JAX human-to-mouse ortholog mapping."""

    pass


@dataclass
class ProteinToGene(_AnnotatedDataFrame):
    """Protein-to-gene identifier mapping."""

    pass
