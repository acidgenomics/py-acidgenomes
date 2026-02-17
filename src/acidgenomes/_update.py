"""Update gene symbols to current nomenclature.

Ported from ``updateGeneSymbols.R``.
"""

from __future__ import annotations

import logging

import pandas as pd

from acidgenomes._constructors import make_hgnc, make_ncbi_gene_info

logger = logging.getLogger(__name__)


def update_gene_symbols(
    gene_names: list[str],
    organism: str,
) -> dict[str, str]:
    """Update gene symbols to current nomenclature.

    For Homo sapiens, checks against HGNC.
    For other organisms, checks against NCBI.

    Parameters
    ----------
    gene_names : list[str]
    organism : str

    Returns
    -------
    dict[str, str]
        Mapping of input gene names to current official gene names.

    Raises
    ------
    ValueError
        If any gene name cannot be matched.
    """
    if not gene_names:
        return {}
    if len(gene_names) != len(set(gene_names)):
        raise ValueError("gene_names must not contain duplicates.")
    if any(pd.isna(g) for g in gene_names):
        raise ValueError("gene_names must not contain NA values.")
    if organism == "Homo sapiens":
        ref = make_hgnc()
        df = ref.data.copy()
        search_cols = ["gene_name", "alias_symbol", "prev_symbol"]
    else:
        ref = make_ncbi_gene_info(organism=organism)
        df = ref.data.copy()
        search_cols = ["gene_name", "gene_synonyms"]
    search_cols = [c for c in search_cols if c in df.columns]
    if "gene_name" not in search_cols:
        raise ValueError("Reference data missing 'gene_name' column.")
    result: dict[str, str] = {}
    failures: list[str] = []
    for g in gene_names:
        matched_name = _match_nested(g, df, search_cols)
        if matched_name is not None:
            result[g] = matched_name
        else:
            failures.append(g)
    if failures:
        raise ValueError(
            f"Failed to match {len(failures)} gene symbol(s): " + ", ".join(failures[:20])
        )
    return result


def _match_nested(
    query: str,
    df: pd.DataFrame,
    search_cols: list[str],
) -> str | None:
    """Match a gene name against multiple columns."""
    for col in search_cols:
        if col not in df.columns:
            continue
        for idx, val in df[col].items():
            if pd.isna(val):
                continue
            val_str = str(val)
            if val_str == query:
                return str(df.at[idx, "gene_name"])
            if "|" in val_str:
                parts = [p.strip() for p in val_str.split("|")]
                if query in parts:
                    return str(df.at[idx, "gene_name"])
    return None
