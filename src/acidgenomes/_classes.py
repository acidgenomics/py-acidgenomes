"""S4 class equivalents as Python dataclasses.

The R package defines 28 S4 classes. Here they are represented as
``@dataclass`` wrappers around a ``pandas.DataFrame`` (for tabular data) or
a BiocPy ``GenomicRanges`` object (for genomic interval data).

Genomic range classes (all ``_GRangesBase`` subclasses) now store data
internally as a ``GenomicRanges`` object, with a ``.data`` property that
returns a ``pandas.DataFrame`` for backward compatibility.
"""

from dataclasses import dataclass, field
from typing import Any

import pandas as pd
from genomicranges import GenomicRanges

from acidgenomes.gff._convert import dataframe_to_granges

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
    if bool(n):
        raise ValueError(f"{cls}: {n} NA value(s) in column '{col}'.")


# ---------------------------------------------------------------------------
# Base annotated DataFrame (tabular data — HGNC, NCBI, mapping tables)
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
        """Return a human-readable summary string."""
        parts = [repr(self)]
        for k, v in self.metadata.items():
            parts.append(f"  {k}: {v}")
        return "\n".join(parts)


# ---------------------------------------------------------------------------
# GRanges-backed base (genomic interval data)
# ---------------------------------------------------------------------------


@dataclass
class _GRangesWrapper:
    """Base for classes that store data as a BiocPy GenomicRanges object.

    Accepts either a GenomicRanges or a pandas DataFrame on construction.
    When given a DataFrame, it is stored directly in ``_df`` as a fallback;
    when given a GenomicRanges, it is stored in ``_gr``.

    The ``.data`` property always returns a DataFrame (from ``.gr.to_pandas()``
    or the stored DataFrame). The ``.granges`` property returns the underlying
    GenomicRanges, building it lazily from ``_df`` if needed.
    """

    _gr: GenomicRanges | None = field(default=None, init=False, repr=False)
    _df: pd.DataFrame | None = field(default=None, init=False, repr=False)
    metadata: dict[str, Any] = field(default_factory=dict)

    def __init__(
        self,
        data: GenomicRanges | pd.DataFrame,
        metadata: dict[str, Any] | None = None,
    ) -> None:
        self.metadata = metadata or {}
        if isinstance(data, GenomicRanges):
            self._gr = data
            self._df = None
        else:
            self._df = data
            self._gr = None

    def __len__(self) -> int:
        if self._gr is not None:
            return len(self._gr)
        return len(self._df) if self._df is not None else 0

    def __repr__(self) -> str:
        name = type(self).__name__
        return f"{name}: {len(self)} ranges"

    @property
    def granges(self) -> GenomicRanges:
        """Return the underlying GenomicRanges object, building it if needed."""
        if self._gr is not None:
            return self._gr
        if self._df is not None:
            # Try to infer the primary identifier column.
            for col in ("gene_id", "transcript_id", "exon_id"):
                if col in self._df.columns:
                    self._gr = dataframe_to_granges(self._df, names_col=col, metadata=self.metadata)
                    return self._gr
            msg = "Cannot determine names_col for GenomicRanges conversion."
            raise ValueError(msg)
        msg = "No data stored in this object."
        raise ValueError(msg)

    @property
    def data(self) -> pd.DataFrame:
        """Return data as a pandas DataFrame (backward-compatible accessor)."""
        if self._df is not None and self._gr is None:
            return self._df
        gr = self.granges
        return gr.to_pandas()

    @property
    def organism(self) -> str | None:
        """Organism name from metadata."""
        return self.metadata.get("organism")

    @property
    def provider(self) -> str | None:
        """Provider name from metadata."""
        return self.metadata.get("provider")

    @property
    def genome_build(self) -> str | None:
        """Genome build from metadata."""
        return self.metadata.get("genome_build")

    def summary(self) -> str:
        """Return a human-readable summary string."""
        parts = [repr(self)]
        for k, v in self.metadata.items():
            parts.append(f"  {k}: {v}")
        return "\n".join(parts)


# Keep _GRangesBase as an alias for backward compatibility.
_GRangesBase = _GRangesWrapper


# ---------------------------------------------------------------------------
# Provider-specific bases
# ---------------------------------------------------------------------------


class _EnsemblBase(_GRangesWrapper):
    def __init__(
        self, data: GenomicRanges | pd.DataFrame, metadata: dict[str, Any] | None = None
    ) -> None:
        m = dict(metadata) if metadata else {}
        m.setdefault("provider", "Ensembl")
        super().__init__(data, m)


class _FlybaseBase(_GRangesWrapper):
    def __init__(
        self, data: GenomicRanges | pd.DataFrame, metadata: dict[str, Any] | None = None
    ) -> None:
        m = dict(metadata) if metadata else {}
        m.setdefault("provider", "FlyBase")
        super().__init__(data, m)


class _GencodeBase(_GRangesWrapper):
    def __init__(
        self, data: GenomicRanges | pd.DataFrame, metadata: dict[str, Any] | None = None
    ) -> None:
        m = dict(metadata) if metadata else {}
        m.setdefault("provider", "GENCODE")
        super().__init__(data, m)


class _RefseqBase(_GRangesWrapper):
    def __init__(
        self, data: GenomicRanges | pd.DataFrame, metadata: dict[str, Any] | None = None
    ) -> None:
        m = dict(metadata) if metadata else {}
        m.setdefault("provider", "RefSeq")
        super().__init__(data, m)


class _UcscBase(_GRangesWrapper):
    def __init__(
        self, data: GenomicRanges | pd.DataFrame, metadata: dict[str, Any] | None = None
    ) -> None:
        m = dict(metadata) if metadata else {}
        m.setdefault("provider", "UCSC")
        super().__init__(data, m)


class _WormbaseBase(_GRangesWrapper):
    def __init__(
        self, data: GenomicRanges | pd.DataFrame, metadata: dict[str, Any] | None = None
    ) -> None:
        m = dict(metadata) if metadata else {}
        m.setdefault("provider", "WormBase")
        super().__init__(data, m)


# ---------------------------------------------------------------------------
# Ensembl
# ---------------------------------------------------------------------------


class EnsemblGenes(_EnsemblBase):
    """Ensembl gene annotations."""


class EnsemblTranscripts(_EnsemblBase):
    """Ensembl transcript annotations."""


class EnsemblExons(_EnsemblBase):
    """Ensembl exon annotations."""


# ---------------------------------------------------------------------------
# FlyBase
# ---------------------------------------------------------------------------


class FlybaseGenes(_FlybaseBase):
    """FlyBase gene annotations."""


class FlybaseTranscripts(_FlybaseBase):
    """FlyBase transcript annotations."""


class FlybaseExons(_FlybaseBase):
    """FlyBase exon annotations."""


# ---------------------------------------------------------------------------
# GENCODE
# ---------------------------------------------------------------------------


class GencodeGenes(_GencodeBase):
    """GENCODE gene annotations."""


class GencodeTranscripts(_GencodeBase):
    """GENCODE transcript annotations."""


class GencodeExons(_GencodeBase):
    """GENCODE exon annotations."""


# ---------------------------------------------------------------------------
# RefSeq
# ---------------------------------------------------------------------------


class RefseqGenes(_RefseqBase):
    """RefSeq gene annotations."""


class RefseqTranscripts(_RefseqBase):
    """RefSeq transcript annotations."""


class RefseqExons(_RefseqBase):
    """RefSeq exon annotations."""


# ---------------------------------------------------------------------------
# UCSC
# ---------------------------------------------------------------------------


class UcscGenes(_UcscBase):
    """UCSC gene annotations."""


class UcscTranscripts(_UcscBase):
    """UCSC transcript annotations."""


class UcscExons(_UcscBase):
    """UCSC exon annotations."""


# ---------------------------------------------------------------------------
# WormBase
# ---------------------------------------------------------------------------


class WormbaseGenes(_WormbaseBase):
    """WormBase gene annotations."""


class WormbaseTranscripts(_WormbaseBase):
    """WormBase transcript annotations."""


class WormbaseExons(_WormbaseBase):
    """WormBase exon annotations."""


# ---------------------------------------------------------------------------
# Gene nomenclature / mapping classes (tabular data — keep DataFrame-backed)
# ---------------------------------------------------------------------------


@dataclass
class Hgnc(_AnnotatedDataFrame):
    """HGNC (Human Gene Nomenclature Committee) data."""


@dataclass
class Mgi(_AnnotatedDataFrame):
    """MGI (Mouse Genome Informatics) data."""


@dataclass
class NcbiGeneInfo(_AnnotatedDataFrame):
    """NCBI gene information."""


@dataclass
class NcbiGeneHistory(_AnnotatedDataFrame):
    """NCBI gene history."""


@dataclass
class EnsemblToNcbi(_AnnotatedDataFrame):
    """Ensembl-to-NCBI gene identifier mapping."""


@dataclass
class NcbiToEnsembl(_AnnotatedDataFrame):
    """NCBI-to-Ensembl gene identifier mapping."""


@dataclass
class GeneToSymbol(_AnnotatedDataFrame):
    """Gene identifier-to-symbol mapping."""


@dataclass
class TxToGene(_AnnotatedDataFrame):
    """Transcript-to-gene identifier mapping."""


@dataclass
class JaxHumanToMouse(_AnnotatedDataFrame):
    """JAX human-to-mouse ortholog mapping."""


@dataclass
class ProteinToGene(_AnnotatedDataFrame):
    """Protein-to-gene identifier mapping."""
