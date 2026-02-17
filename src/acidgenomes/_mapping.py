"""Gene name / ID mapping and import functions.

Ported from ``mapGeneNamesToEnsembl.R``, ``mapGeneNamesToHgnc.R``,
``mapGeneNamesToNcbi.R``, ``mapGencodeToEnsembl.R``,
``mapEnsemblReleaseToUrl.R``, ``importTxToGene.R``,
``mapHumanOrthologs.R``.
"""

from __future__ import annotations

import logging
import re

import pandas as pd
import requests

from acidgenomes._cache import fetch_text
from acidgenomes._classes import Hgnc, NcbiGeneInfo, TxToGene
from acidgenomes._constructors import (
    make_hgnc,
    make_ncbi_gene_info,
    make_tx_to_gene,
)
from acidgenomes._detect import detect_organism
from acidgenomes._strip_versions import (
    strip_gene_versions,
    strip_transcript_versions,
)

logger = logging.getLogger(__name__)


# -------------------------------------------------------------------------
# mapGeneNamesToHgnc
# -------------------------------------------------------------------------


def map_gene_names_to_hgnc(
    genes: list[str],
    *,
    ignore_case: bool = False,
    hgnc: Hgnc | None = None,
) -> list[int]:
    """Map human gene names (symbols) to HGNC identifiers.

    Parameters
    ----------
    genes : list[str]
    ignore_case : bool
    hgnc : Hgnc or None

    Returns
    -------
    list[int]

    Raises
    ------
    ValueError
        If any gene name cannot be matched.
    """
    if hgnc is None:
        hgnc = make_hgnc()
    df = hgnc.data.copy()
    required = ["gene_name", "hgnc_id"]
    for col in required:
        if col not in df.columns:
            raise ValueError(f"HGNC data missing column '{col}'.")
    lookup = _build_name_lookup(df, "hgnc_id", ignore_case)
    for col in ("alias_symbol", "prev_symbol"):
        if col in df.columns:
            _extend_lookup_from_pipe_col(lookup, df, col, "hgnc_id", ignore_case)
    result, failures = _resolve_genes(genes, lookup, ignore_case)
    if failures:
        raise ValueError(f"{len(failures)} mapping failure(s): " + ", ".join(failures[:20]))
    return result


# -------------------------------------------------------------------------
# mapGeneNamesToNcbi
# -------------------------------------------------------------------------


def map_gene_names_to_ncbi(
    genes: list[str],
    organism: str,
    *,
    taxonomic_group: str | None = None,
    ignore_case: bool = False,
    ncbi: NcbiGeneInfo | None = None,
) -> list[int]:
    """Map gene names to NCBI (Entrez) gene identifiers.

    Parameters
    ----------
    genes : list[str]
    organism : str
    taxonomic_group : str or None
    ignore_case : bool
    ncbi : NcbiGeneInfo or None

    Returns
    -------
    list[int]
    """
    if ncbi is None:
        ncbi = make_ncbi_gene_info(organism=organism, taxonomic_group=taxonomic_group)
    df = ncbi.data.copy()
    if "gene_name" not in df.columns:
        raise ValueError("NCBI data missing 'gene_name' column.")
    id_col = "gene_id" if "gene_id" in df.columns else df.columns[0]
    lookup = _build_name_lookup(df, id_col, ignore_case)
    if "gene_synonyms" in df.columns:
        _extend_lookup_from_pipe_col(lookup, df, "gene_synonyms", id_col, ignore_case)
    result, failures = _resolve_genes(genes, lookup, ignore_case)
    if failures:
        raise ValueError(f"{len(failures)} mapping failure(s): " + ", ".join(failures[:20]))
    return result


# -------------------------------------------------------------------------
# mapGeneNamesToEnsembl
# -------------------------------------------------------------------------


def map_gene_names_to_ensembl(
    genes: list[str],
    organism: str,
    *,
    ignore_case: bool = False,
    hgnc: Hgnc | None = None,
    ncbi: NcbiGeneInfo | None = None,
) -> list[str]:
    """Map gene names to Ensembl gene identifiers.

    For Homo sapiens, defaults to HGNC; otherwise uses NCBI.

    Parameters
    ----------
    genes : list[str]
    organism : str
    ignore_case : bool
    hgnc : Hgnc or None
    ncbi : NcbiGeneInfo or None

    Returns
    -------
    list[str]
    """
    if organism == "Homo sapiens" and ncbi is None:
        return _map_genes_to_ensembl_via_hgnc(genes, ignore_case=ignore_case, hgnc=hgnc)
    return _map_genes_to_ensembl_via_ncbi(
        genes, organism=organism, ignore_case=ignore_case, ncbi=ncbi
    )


def _map_genes_to_ensembl_via_hgnc(
    genes: list[str],
    *,
    ignore_case: bool = False,
    hgnc: Hgnc | None = None,
) -> list[str]:
    """Map genes to Ensembl IDs via HGNC (Homo sapiens)."""
    if hgnc is None:
        hgnc = make_hgnc()
    hgnc_df = hgnc.data
    if "ensembl_gene_id" not in hgnc_df.columns:
        raise ValueError("HGNC data missing 'ensembl_gene_id' column.")
    hids = map_gene_names_to_hgnc(genes, ignore_case=ignore_case, hgnc=hgnc)
    id_map = dict(
        zip(
            hgnc_df["hgnc_id"].dropna().astype(int),
            hgnc_df["ensembl_gene_id"],
            strict=False,
        )
    )
    result: list[str] = []
    failures: list[str] = []
    for g, hid in zip(genes, hids, strict=True):
        ens = id_map.get(hid)
        if ens and pd.notna(ens):
            result.append(str(ens))
        else:
            failures.append(g)
            result.append("")
    if failures:
        raise ValueError(f"{len(failures)} mapping failure(s): " + ", ".join(failures[:20]))
    return result


def _map_genes_to_ensembl_via_ncbi(
    genes: list[str],
    *,
    organism: str,
    ignore_case: bool = False,
    ncbi: NcbiGeneInfo | None = None,
) -> list[str]:
    """Map genes to Ensembl IDs via NCBI db_xrefs."""
    if ncbi is None:
        ncbi = make_ncbi_gene_info(organism=organism)
    ncbi_df = ncbi.data
    if "db_xrefs" not in ncbi_df.columns:
        raise ValueError("NCBI data missing 'db_xrefs' column needed for Ensembl mapping.")
    id_col = "gene_id" if "gene_id" in ncbi_df.columns else ncbi_df.columns[0]
    ens_map = _build_ensembl_xref_map(ncbi_df, id_col)
    ncbi_ids = map_gene_names_to_ncbi(
        genes,
        organism=organism,
        ignore_case=ignore_case,
        ncbi=ncbi,
    )
    result: list[str] = []
    failures: list[str] = []
    for g, nid in zip(genes, ncbi_ids, strict=True):
        ens = ens_map.get(nid)
        if ens:
            result.append(ens)
        else:
            failures.append(g)
            result.append("")
    if failures:
        raise ValueError(f"{len(failures)} mapping failure(s): " + ", ".join(failures[:20]))
    return result


def _build_ensembl_xref_map(df: pd.DataFrame, id_col: str) -> dict[int, str]:
    """Extract Ensembl gene IDs from NCBI db_xrefs column."""
    ens_map: dict[int, str] = {}
    for _, row in df.iterrows():
        xrefs = str(row.get("db_xrefs", ""))
        for ref in xrefs.split("|"):
            if ref.startswith("Ensembl:"):
                ens_map[int(row[id_col])] = ref.replace("Ensembl:", "")
                break
    return ens_map


# -------------------------------------------------------------------------
# mapGencodeToEnsembl
# -------------------------------------------------------------------------


def map_gencode_to_ensembl(release: int | str) -> int:
    """Map a GENCODE release to its corresponding Ensembl release.

    Parameters
    ----------
    release : int or str

    Returns
    -------
    int
    """
    release_str = str(release)
    short = "mouse" if release_str.startswith("M") else "human"
    url = f"https://www.gencodegenes.org/{short}/releases.html"
    text = fetch_text(url)
    for line in text.splitlines():
        if release_str in line and "Ensembl" in line:
            m = re.search(r"Ensembl\s+(\d+)", line)
            if m:
                return int(m.group(1))
    rows = re.findall(
        r"<td[^>]*>\s*" + re.escape(release_str) + r"\s*</td>"
        r".*?<td[^>]*>\s*(\d+)\s*</td>",
        text,
        re.DOTALL,
    )
    if rows:
        return int(rows[0])
    raise ValueError(f"Failed to match GENCODE release '{release}' to Ensembl.")


# -------------------------------------------------------------------------
# mapEnsemblReleaseToUrl
# -------------------------------------------------------------------------


def map_ensembl_release_to_url(release: int | None = None) -> str:
    """Map an Ensembl release version to its archive URL.

    Parameters
    ----------
    release : int or None

    Returns
    -------
    str
    """
    current = "https://useast.ensembl.org"
    if release is None:
        return current
    url = "https://useast.ensembl.org/info/website/archives/index.html"
    text = fetch_text(url)
    pattern = r"Ensembl\s+(\d+)\D+?([A-Za-z]{3})\s+(\d{4})"
    for m in re.finditer(pattern, text):
        ver, month, year = m.group(1), m.group(2), m.group(3)
        if int(ver) == release:
            slug = f"{month.lower()}{year}"
            return f"https://{slug}.archive.ensembl.org"
    raise ValueError(f"Unsupported Ensembl release: {release}.")


# -------------------------------------------------------------------------
# importTxToGene
# -------------------------------------------------------------------------


def import_tx_to_gene(
    file: str,
    *,
    organism: str | None = None,
    genome_build: str | None = None,
    release: int | str | None = None,
    ignore_tx_version: bool = False,
    ignore_gene_version: bool = False,
) -> TxToGene:
    """Import transcript-to-gene annotations from a two-column file.

    Parameters
    ----------
    file : str
    organism : str or None
    genome_build : str or None
    release : int or str or None
    ignore_tx_version : bool
    ignore_gene_version : bool

    Returns
    -------
    TxToGene
    """
    df = pd.read_csv(file, header=None, names=["tx_id", "gene_id"])
    if ignore_tx_version:
        df["tx_id"] = strip_transcript_versions(df["tx_id"].tolist())
    if ignore_gene_version:
        df["gene_id"] = strip_gene_versions(df["gene_id"].tolist())
    return make_tx_to_gene(
        df,
        organism=organism,
        genome_build=genome_build,
        release=release,
    )


# -------------------------------------------------------------------------
# mapHumanOrthologs
# -------------------------------------------------------------------------


def map_human_orthologs(
    genes: list[str],
    organism: str | None = None,
    ensembl_release: int | None = None,
) -> pd.DataFrame:
    """Map gene identifiers to human orthologs via Ensembl REST API.

    Parameters
    ----------
    genes : list[str]
    organism : str or None
    ensembl_release : int or None

    Returns
    -------
    pd.DataFrame
    """
    if organism is None:
        organism = detect_organism(genes)
    if organism == "Homo sapiens":
        raise ValueError("Input genes must not be Homo sapiens.")
    base = map_ensembl_release_to_url(ensembl_release)
    results: list[dict[str, str]] = []
    batch_size = 50
    for i in range(0, len(genes), batch_size):
        batch = genes[i : i + batch_size]
        for gene_id in batch:
            _fetch_ortholog(base, gene_id, results)
    if not results:
        raise RuntimeError("Failed to map any genes to human orthologs.")
    df = pd.DataFrame(results).drop_duplicates()
    df = df.drop_duplicates(subset=["gene_id"], keep="first")
    df = df.drop_duplicates(subset=["human_gene_id"], keep="first")
    df = df.sort_values(["gene_id"]).reset_index(drop=True)
    return df


def _fetch_ortholog(
    base: str,
    gene_id: str,
    results: list[dict[str, str]],
) -> None:
    """Fetch a single ortholog from Ensembl REST API."""
    url = (
        f"{base}/homology/id/{gene_id}"
        f"?type=orthologues;target_species=homo_sapiens"
        f";content-type=application/json"
    )
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        homologies = data.get("data", [{}])[0].get("homologies", [])
        for h in homologies:
            target = h.get("target", {})
            human_id = target.get("id")
            if human_id:
                results.append({"gene_id": gene_id, "human_gene_id": human_id})
                break
    except Exception:
        logger.warning("Failed to fetch orthologs for %s", gene_id)


# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------


def _resolve_genes(
    genes: list[str],
    lookup: dict[str, int],
    ignore_case: bool,
) -> tuple[list[int], list[str]]:
    """Resolve a list of gene names against a lookup dict.

    Returns
    -------
    tuple[list[int], list[str]]
        A tuple of (results, failures).
    """
    result: list[int] = []
    failures: list[str] = []
    for g in genes:
        key = g.upper() if ignore_case else g
        gid = lookup.get(key)
        if gid is not None:
            result.append(int(gid))
        else:
            failures.append(g)
            result.append(-1)
    return result, failures


def _build_name_lookup(
    df: pd.DataFrame,
    id_col: str,
    ignore_case: bool,
) -> dict[str, int]:
    """Build gene_name -> id lookup dict."""
    out: dict[str, int] = {}
    for _, row in df.iterrows():
        name = row.get("gene_name")
        gid = row.get(id_col)
        if pd.isna(name) or pd.isna(gid):
            continue
        key = str(name).upper() if ignore_case else str(name)
        if key not in out:
            out[key] = gid
    return out


def _extend_lookup_from_pipe_col(
    lookup: dict[str, int],
    df: pd.DataFrame,
    col: str,
    id_col: str,
    ignore_case: bool,
) -> None:
    """Extend lookup dict from a pipe-delimited column."""
    for _, row in df.iterrows():
        val = row.get(col)
        gid = row.get(id_col)
        if pd.isna(val) or pd.isna(gid):
            continue
        for raw_alias in str(val).split("|"):
            cleaned = raw_alias.strip()
            if not cleaned:
                continue
            key = cleaned.upper() if ignore_case else cleaned
            if key not in lookup:
                lookup[key] = gid
