"""Gene name / ID mapping and import functions.

Ported from ``mapGeneNamesToEnsembl.R``, ``mapGeneNamesToHgnc.R``,
``mapGeneNamesToNcbi.R``, ``mapGencodeToEnsembl.R``,
``mapEnsemblReleaseToUrl.R``, ``importTxToGene.R``,
``mapHumanOrthologs.R``, ``classifyCuratedGeneGroups.R``.
"""

from __future__ import annotations

import logging
import re

import pandas as pd
import requests

from acidgenomes._cache import fetch_text
from acidgenomes._classes import Hgnc, JaxHumanToMouse, Mgi, NcbiGeneInfo, TxToGene
from acidgenomes._constructors import (
    make_hgnc,
    make_jax_human_to_mouse,
    make_mgi,
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
        Gene names (symbols) to map.
    ignore_case : bool
        Match case-insensitively.
    hgnc : Hgnc or None
        HGNC reference dataset. Downloaded via :func:`make_hgnc` if
        ``None``.

    Returns
    -------
    list[int]
        HGNC identifiers, in the same order as ``genes``.

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
        Gene names (symbols) to map.
    organism : str
        Latin organism name.
    taxonomic_group : str or None
        NCBI FTP taxonomic group. Auto-detected if ``None``.
    ignore_case : bool
        Match case-insensitively.
    ncbi : NcbiGeneInfo or None
        NCBI gene info reference dataset. Downloaded via
        :func:`make_ncbi_gene_info` if ``None``.

    Returns
    -------
    list[int]
        NCBI (Entrez) gene identifiers, in the same order as ``genes``.
    """
    if ncbi is None:
        ncbi = make_ncbi_gene_info(organism=organism, taxonomic_group=taxonomic_group)
    df = ncbi.data.copy()
    if "gene_name" not in df.columns:
        raise ValueError("NCBI data missing 'gene_name' column.")
    id_col = "gene_id" if "gene_id" in df.columns else str(df.columns[0])
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
        Gene names (symbols) to map.
    organism : str
        Latin organism name.
    ignore_case : bool
        Match case-insensitively.
    hgnc : Hgnc or None
        HGNC reference dataset (used for Homo sapiens). Downloaded via
        :func:`make_hgnc` if ``None``.
    ncbi : NcbiGeneInfo or None
        NCBI gene info reference dataset (used for other organisms, or
        when explicitly provided). Downloaded via
        :func:`make_ncbi_gene_info` if ``None`` and required.

    Returns
    -------
    list[str]
        Ensembl gene identifiers, in the same order as ``genes``.
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
    id_col = "gene_id" if "gene_id" in ncbi_df.columns else str(ncbi_df.columns[0])
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
                ens_map[int(str(row[id_col]))] = ref.replace("Ensembl:", "")
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
        GENCODE release (e.g. ``46`` or ``"M35"`` for mouse).

    Returns
    -------
    int
        Corresponding Ensembl release version.
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
        Ensembl release version. Returns the current (non-archived) site
        URL if ``None``.

    Returns
    -------
    str
        Archive site URL for the requested release, or the current site
        URL if ``release`` is ``None``.
    """
    current = "https://rest.ensembl.org"
    if release is None:
        return current
    url = "https://www.ensembl.org/info/website/archives/index.html"
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
        Path to a headerless, two-column (transcript ID, gene ID) file.
    organism : str or None
        Latin organism name, recorded in the returned object's metadata.
    genome_build : str or None
        Genome build, recorded in the returned object's metadata.
    release : int or str or None
        Annotation release version, recorded in the returned object's
        metadata.
    ignore_tx_version : bool
        Strip version suffixes from transcript identifiers.
    ignore_gene_version : bool
        Strip version suffixes from gene identifiers.

    Returns
    -------
    TxToGene
        Transcript ID to gene ID mapping.
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
        Ensembl gene identifiers from a non-human organism.
    organism : str or None
        Latin organism name of ``genes``. Auto-detected if ``None``.
    ensembl_release : int or None
        Ensembl release version to query. Uses the current release if
        ``None``.

    Returns
    -------
    pd.DataFrame
        Deduplicated mapping with ``gene_id`` and ``human_gene_id``
        columns.
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
# classifyCuratedGeneGroups
# -------------------------------------------------------------------------

# HGNC gene_group_id -> curated tag. Verified live against genenames.org
# 2026-09-24 (hgnc_complete_set.txt): 728 "S ribosomal proteins", 729
# "L ribosomal proteins" (cytoplasmic large/small subunit), 646
# "Mitochondrial ribosomal proteins", 940 "Hemoglobin subunits". The
# "Ribosomal protein S6 kinase family" (1156/1691/3524) is a distinct HGNC
# family (RPS6KA*/RPS6KB*) and is deliberately absent from this map -- do
# not add it. A sibling R package carries the identical constant and
# function as a hand-ported twin; a change here must land there too, in
# the same release.
_CURATED_GENE_GROUP_IDS: dict[int, str] = {
    728: "ribo_cyto",
    729: "ribo_cyto",
    646: "ribo_mito",
    940: "hemoglobin",
}


def classify_curated_gene_groups(
    ensembl_gene_ids: list[str],
    organism: str,
    *,
    hgnc: Hgnc | None = None,
    jax: JaxHumanToMouse | None = None,
    mgi: Mgi | None = None,
) -> dict[str, list[str]]:
    """Classify genes into curated HGNC gene-group tags.

    Tags each gene ``"ribo_cyto"`` (cytoplasmic ribosomal protein),
    ``"ribo_mito"`` (mitochondrial ribosomal protein), ``"hemoglobin"``
    (hemoglobin subunit), or none of these. Tags are sourced from HGNC's own
    curated ``gene_group`` assignments, never a symbol regex; a gene not in
    any curated group gets an empty list, not an omitted key.

    For Mus musculus, human HGNC groups are propagated via a fully
    identifier-based chain (HGNC ``hgnc_id`` -> JAX ortholog
    ``mouse_mgi_id`` -> MGI ``ensembl_gene_id``), with no gene-symbol
    matching at any step.

    This is deliberately independent of ``broad_class``
    (:func:`._add_broad_class`): ``broad_class`` is single-valued and every
    ribosomal/hemoglobin gene already has a value from it (``"coding"``,
    ``"pseudo"``, etc). Do not fold these tags into ``broad_class``.

    Parameters
    ----------
    ensembl_gene_ids : list[str]
        Ensembl gene identifiers to classify.
    organism : str
        Latin organism name. Only ``"Homo sapiens"`` and ``"Mus musculus"``
        are supported.
    hgnc : Hgnc or None
        HGNC reference dataset. Downloaded via :func:`make_hgnc` if
        ``None``.
    jax : JaxHumanToMouse or None
        JAX human-to-mouse ortholog dataset. Downloaded via
        :func:`make_jax_human_to_mouse` if ``None``. Ignored for
        ``"Homo sapiens"``.
    mgi : Mgi or None
        MGI reference dataset. Downloaded via :func:`make_mgi` if ``None``.
        Ignored for ``"Homo sapiens"``.

    Returns
    -------
    dict[str, list[str]]
        Mapping of each input Ensembl gene ID to its curated tag list
        (``[]`` if none).

    Raises
    ------
    ValueError
        If ``organism`` is not ``"Homo sapiens"`` or ``"Mus musculus"``.

    Examples
    --------
    >>> tags = classify_curated_gene_groups(["ENSG00000244734"], "Homo sapiens")
    >>> tags["ENSG00000244734"]  # HBB
    ['hemoglobin']
    """
    if organism == "Homo sapiens":
        tags_by_ensembl_id = _hgnc_curated_tags_by_ensembl_id(hgnc)
    elif organism == "Mus musculus":
        tags_by_ensembl_id = _mouse_curated_tags_by_ensembl_id(hgnc, jax, mgi)
    else:
        raise ValueError(
            f"Unsupported organism '{organism}'; "
            "only 'Homo sapiens' and 'Mus musculus' are supported."
        )
    return {gid: tags_by_ensembl_id.get(gid, []) for gid in ensembl_gene_ids}


def _hgnc_curated_tags_by_ensembl_id(hgnc: Hgnc | None) -> dict[str, list[str]]:
    """Build human Ensembl gene ID -> curated tag list from HGNC gene groups."""
    if hgnc is None:
        hgnc = make_hgnc()
    df = hgnc.data
    for col in ("ensembl_gene_id", "gene_group_id"):
        if col not in df.columns:
            raise ValueError(f"HGNC data missing column '{col}'.")
    out: dict[str, list[str]] = {}
    for _, row in df.iterrows():
        ensembl_id = row.get("ensembl_gene_id")
        group_ids = row.get("gene_group_id")
        if ensembl_id is None or group_ids is None:
            continue
        if bool(pd.isna(ensembl_id)) or bool(pd.isna(group_ids)):
            continue
        tags = _tags_for_group_id_string(str(group_ids))
        if tags:
            out[str(ensembl_id)] = tags
    return out


def _tags_for_group_id_string(group_ids: str) -> list[str]:
    """Resolve a pipe-delimited HGNC ``gene_group_id`` string to curated tags."""
    tags: list[str] = []
    for entry in group_ids.split("|"):
        cleaned = entry.strip()
        if not cleaned:
            continue
        try:
            gid = int(cleaned)
        except ValueError:
            continue
        tag = _CURATED_GENE_GROUP_IDS.get(gid)
        if tag is not None and tag not in tags:
            tags.append(tag)
    return tags


def _tags_by_hgnc_id(hgnc_df: pd.DataFrame) -> dict[int, list[str]]:
    """Build human HGNC ID -> curated tag list from HGNC gene groups."""
    for col in ("hgnc_id", "gene_group_id"):
        if col not in hgnc_df.columns:
            raise ValueError(f"HGNC data missing column '{col}'.")
    out: dict[int, list[str]] = {}
    for _, row in hgnc_df.iterrows():
        hgnc_id = row.get("hgnc_id")
        group_ids = row.get("gene_group_id")
        if hgnc_id is None or group_ids is None:
            continue
        if bool(pd.isna(hgnc_id)) or bool(pd.isna(group_ids)):
            continue
        tags = _tags_for_group_id_string(str(group_ids))
        if tags:
            out[int(hgnc_id)] = tags
    return out


def _mouse_mgi_ids_by_human_hgnc_id(jax_df: pd.DataFrame) -> dict[int, set[int]]:
    """Build human HGNC ID -> mouse MGI ID set from the JAX ortholog table.

    A human gene can have more than one mouse paralog (e.g. hemoglobin's
    Hbb-bh0/bh1/bh2/bh3/y expansion) -- every mouse MGI ID reachable from a
    given human HGNC ID is collected, not just the first/last one seen.
    """
    for col in ("human_hgnc_id", "mouse_mgi_id"):
        if col not in jax_df.columns:
            raise ValueError(f"JAX human-to-mouse data missing column '{col}'.")
    out: dict[int, set[int]] = {}
    for _, row in jax_df.iterrows():
        human_hgnc_id = row.get("human_hgnc_id")
        mouse_mgi_id = row.get("mouse_mgi_id")
        if human_hgnc_id is None or mouse_mgi_id is None:
            continue
        if bool(pd.isna(human_hgnc_id)) or bool(pd.isna(mouse_mgi_id)):
            continue
        out.setdefault(int(human_hgnc_id), set()).add(int(mouse_mgi_id))
    return out


def _mouse_ensembl_id_by_mgi_id(mgi_df: pd.DataFrame) -> dict[int, str]:
    """Build mouse MGI ID -> Ensembl gene ID from the MGI gene-model report."""
    for col in ("mgi_accession_id", "ensembl_gene_id"):
        if col not in mgi_df.columns:
            raise ValueError(f"MGI data missing column '{col}'.")
    out: dict[int, str] = {}
    for _, row in mgi_df.iterrows():
        mgi_id = row.get("mgi_accession_id")
        ensembl_id = row.get("ensembl_gene_id")
        if mgi_id is None or ensembl_id is None:
            continue
        if bool(pd.isna(mgi_id)) or bool(pd.isna(ensembl_id)):
            continue
        out[int(mgi_id)] = str(ensembl_id)
    return out


def _mouse_curated_tags_by_ensembl_id(
    hgnc: Hgnc | None,
    jax: JaxHumanToMouse | None,
    mgi: Mgi | None,
) -> dict[str, list[str]]:
    """Propagate human HGNC curated tags to mouse Ensembl gene IDs.

    Fully identifier-based: HGNC ``hgnc_id`` -> JAX ``human_hgnc_id`` /
    ``mouse_mgi_id`` -> MGI ``mgi_accession_id`` / ``ensembl_gene_id``. No
    gene-symbol matching at any step, unlike a naive mouse-symbol regex
    (which silently returns zero matches -- hemoglobin symbols are
    hyphenated in mouse, e.g. ``Hba-a1``, and share no substring with any
    human-derived pattern).
    """
    if hgnc is None:
        hgnc = make_hgnc()
    if jax is None:
        jax = make_jax_human_to_mouse(unique=False)
    if mgi is None:
        mgi = make_mgi()

    tags_by_hgnc_id = _tags_by_hgnc_id(hgnc.data)
    mgi_ids_by_hgnc_id = _mouse_mgi_ids_by_human_hgnc_id(jax.data)
    ensembl_id_by_mgi_id = _mouse_ensembl_id_by_mgi_id(mgi.data)

    out: dict[str, list[str]] = {}
    for hgnc_id, tags in tags_by_hgnc_id.items():
        for mgi_id in mgi_ids_by_hgnc_id.get(hgnc_id, ()):
            ensembl_id = ensembl_id_by_mgi_id.get(mgi_id)
            if ensembl_id is not None:
                out[ensembl_id] = tags
    return out


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
        if gid is None or bool(pd.isna(name)) or bool(pd.isna(gid)):
            continue
        key = str(name).upper() if ignore_case else str(name)
        if key not in out:
            out[key] = int(gid)
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
        if gid is None or bool(pd.isna(val)) or bool(pd.isna(gid)):
            continue
        for raw_alias in str(val).split("|"):
            cleaned = raw_alias.strip()
            if not cleaned:
                continue
            key = cleaned.upper() if ignore_case else cleaned
            if key not in lookup:
                lookup[key] = int(gid)
