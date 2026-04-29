"""Constructor functions that fetch data and return validated class instances.

Ported from ``Hgnc.R``, ``Mgi.R``, ``NcbiGeneInfo.R``, ``NcbiGeneHistory.R``,
``JaxHumanToMouse.R``, ``EnsemblGenes.R``, ``EnsemblToNcbi-methods.R``,
``NcbiToEnsembl-methods.R``, ``GeneToSymbol-methods.R``,
``TxToGene-methods.R``.
"""

from __future__ import annotations

import gzip
import io
import logging
import re
from datetime import date
from typing import Any

import pandas as pd

from acidgenomes._cache import cache_url
from acidgenomes._classes import (
    EnsemblGenes,
    EnsemblToNcbi,
    GeneToSymbol,
    Hgnc,
    JaxHumanToMouse,
    Mgi,
    NcbiGeneHistory,
    NcbiGeneInfo,
    NcbiToEnsembl,
    TxToGene,
)
from acidgenomes._data import NCBI_TAX_IDS, NCBI_TAXONOMIC_GROUPS
from acidgenomes._detect import detect_organism
from acidgenomes._genome_build import current_ensembl_genome_build
from acidgenomes._genome_version import current_ensembl_version

logger = logging.getLogger(__name__)


# -------------------------------------------------------------------------
# Helper: camelCase -> snake_case
# -------------------------------------------------------------------------


def _to_snake(name: str) -> str:
    """Convert camelCase or PascalCase to snake_case."""
    s1 = re.sub(r"([A-Z]+)([A-Z][a-z])", r"\1_\2", name)
    s2 = re.sub(r"([a-z\d])([A-Z])", r"\1_\2", s1)
    return s2.lower().replace(".", "_").replace(" ", "_").replace("#", "")


def _rename_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Rename columns from camelCase/Title to snake_case."""
    df = df.copy()
    df.columns = [_to_snake(c) for c in df.columns]
    return df


# -------------------------------------------------------------------------
# HGNC
# -------------------------------------------------------------------------


def make_hgnc() -> Hgnc:
    """Import HGNC (Human Gene Nomenclature Committee) complete set.

    Downloads the current complete HGNC dataset from genenames.org and
    returns a validated Hgnc object.

    Returns
    -------
    Hgnc

    See Also
    --------
    https://www.genenames.org/download/statistics-and-files/
    """
    logger.info("Importing HGNC complete set.")
    url = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
    path = cache_url(url)
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df = _rename_columns(df)
    renames: dict[str, str] = {
        "entrez_id": "ncbi_gene_id",
        "symbol": "gene_name",
        "name": "description",
    }
    df = df.rename(columns={k: v for k, v in renames.items() if k in df.columns})
    if "hgnc_id" in df.columns:
        df["hgnc_id"] = df["hgnc_id"].astype(str).str.replace("HGNC:", "", regex=False)
        df["hgnc_id"] = pd.to_numeric(df["hgnc_id"], errors="coerce")
        df = df.sort_values("hgnc_id").reset_index(drop=True)
    for col in (
        "date_approved_reserved",
        "date_modified",
        "date_name_changed",
        "date_symbol_changed",
    ):
        if col in df.columns:
            df[col] = pd.to_datetime(df[col], errors="coerce")
    if "ncbi_gene_id" in df.columns:
        df["ncbi_gene_id"] = pd.to_numeric(df["ncbi_gene_id"], errors="coerce")
    if "ensembl_gene_id" not in df.columns and "ensembl_id" in df.columns:
        df = df.rename(columns={"ensembl_id": "ensembl_gene_id"})
    meta: dict[str, Any] = {
        "date": date.today(),
        "organism": "Homo sapiens",
        "url": url,
    }
    return Hgnc(data=df, metadata=meta)


# -------------------------------------------------------------------------
# MGI
# -------------------------------------------------------------------------


def make_mgi() -> Mgi:
    """Import Mouse Genome Informatics (MGI) metadata.

    Returns
    -------
    Mgi

    See Also
    --------
    https://www.informatics.jax.org/
    """
    logger.info("Importing MGI metadata.")
    url = "https://www.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt"
    path = cache_url(url)
    with open(path) as fh:
        lines = fh.readlines()
    header_parts = lines[0].strip().split("\t")
    col_names = [re.sub(r"^\d+\.\s*", "", p) for p in header_parts]
    col_names = [_to_snake(c) for c in col_names]
    col_names = ["ncbi_gene_id" if c == "entrez_gene_id" else c for c in col_names]
    data_text = "".join(lines[1:])
    ncols = len(col_names)
    df = pd.read_csv(
        io.StringIO(data_text),
        sep="\t",
        header=None,
        names=col_names + ["_extra"],
        na_values=["null", "NA"],
        dtype=str,
    )
    if "_extra" in df.columns:
        df = df.drop(columns=["_extra"])
    df = df.iloc[:, :ncols]
    if "mgi_accession_id" in df.columns:
        df["mgi_accession_id"] = (
            df["mgi_accession_id"].astype(str).str.replace("MGI:", "", regex=False)
        )
        df["mgi_accession_id"] = pd.to_numeric(df["mgi_accession_id"], errors="coerce")
        df = df.sort_values("mgi_accession_id").reset_index(drop=True)
    meta: dict[str, Any] = {
        "date": date.today(),
        "organism": "Mus musculus",
        "url": url,
    }
    return Mgi(data=df, metadata=meta)


# -------------------------------------------------------------------------
# NCBI Gene Info
# -------------------------------------------------------------------------


def make_ncbi_gene_info(
    organism: str,
    taxonomic_group: str | None = None,
) -> NcbiGeneInfo:
    """Import NCBI gene information for an organism.

    Parameters
    ----------
    organism : str
        Latin organism name (e.g. 'Homo sapiens').
    taxonomic_group : str or None
        NCBI FTP taxonomic group. Auto-detected if None.

    Returns
    -------
    NcbiGeneInfo
    """
    if taxonomic_group is None:
        groups = NCBI_TAXONOMIC_GROUPS.get(organism)
        if groups is None:
            raise ValueError(
                f"Unknown taxonomic group for '{organism}'; set taxonomic_group explicitly."
            )
        taxonomic_group = groups["gene_info"]
    org_file = organism.replace(" ", "_")
    url = (
        f"https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/"
        f"{taxonomic_group}/{org_file}.gene_info.gz"
    )
    logger.info("Downloading %s gene info from NCBI.", organism)
    path = cache_url(url)
    df = pd.read_csv(path, sep="\t", na_values=["-"], low_memory=False)
    df = _rename_columns(df)
    renames = {"symbol": "gene_name", "synonyms": "gene_synonyms"}
    df = df.rename(columns={k: v for k, v in renames.items() if k in df.columns})
    if "gene_id" in df.columns:
        df["gene_id"] = pd.to_numeric(df["gene_id"], errors="coerce")
        df = df.sort_values("gene_id").reset_index(drop=True)
    meta: dict[str, Any] = {
        "date": date.today(),
        "organism": organism,
        "taxonomic_group": taxonomic_group,
        "url": url,
    }
    return NcbiGeneInfo(data=df, metadata=meta)


# -------------------------------------------------------------------------
# NCBI Gene History
# -------------------------------------------------------------------------


def make_ncbi_gene_history(organism: str) -> NcbiGeneHistory:
    """Import NCBI gene history for an organism.

    Parameters
    ----------
    organism : str

    Returns
    -------
    NcbiGeneHistory
    """
    tax_id = NCBI_TAX_IDS.get(organism)
    if tax_id is None:
        raise ValueError(f"No NCBI taxonomy ID known for '{organism}'.")
    url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_history.gz"
    logger.info("Downloading NCBI gene history.")
    path = cache_url(url)
    df = pd.read_csv(path, sep="\t", na_values=["-"], comment="#", low_memory=False)
    df = _rename_columns(df)
    tax_col = _find_tax_column(df)
    df[tax_col] = pd.to_numeric(df[tax_col], errors="coerce")
    df = df[df[tax_col] == tax_id].copy()
    df = df.drop(columns=[tax_col], errors="ignore")
    df = _rename_columns(df)
    _coerce_date_columns(df)
    meta: dict[str, Any] = {
        "date": date.today(),
        "organism": organism,
        "taxonomy_id": tax_id,
        "url": url,
    }
    return NcbiGeneHistory(data=df, metadata=meta)


def _find_tax_column(df: pd.DataFrame) -> str:
    """Find the taxonomy column in a DataFrame."""
    for candidate in ("x_tax_id", "#tax_id", "tax_id"):
        c_snake = _to_snake(candidate)
        if c_snake in df.columns:
            return c_snake
        if candidate in df.columns:
            return candidate
    return df.columns[0]


def _coerce_date_columns(df: pd.DataFrame) -> None:
    """Coerce date-like columns to datetime in place."""
    for col in df.columns:
        if "date" in col.lower():
            df[col] = pd.to_datetime(
                df[col]
                .astype(str)
                .str.replace(
                    r"^(\d{4})(\d{2})(\d{2})$",
                    r"\1-\2-\3",
                    regex=True,
                ),
                errors="coerce",
            )


# -------------------------------------------------------------------------
# JaxHumanToMouse
# -------------------------------------------------------------------------


def make_jax_human_to_mouse(
    unique: bool = True,
) -> JaxHumanToMouse:
    """Import JAX human-to-mouse ortholog mapping.

    Parameters
    ----------
    unique : bool

    Returns
    -------
    JaxHumanToMouse
    """
    logger.info("Importing JAX human-to-mouse orthologs.")
    url = "https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt"
    path = cache_url(url)
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df = _rename_columns(df)
    df = _clean_jax_ids(df)
    renames = {
        "entrez_gene_id": "ncbi_gene_id",
        "symbol": "gene_name",
    }
    df = df.rename(columns={k: v for k, v in renames.items() if k in df.columns})
    out = _merge_jax_species(df)
    if unique:
        out = out.drop_duplicates(subset=["human_gene_name"], keep="first")
        out = out.drop_duplicates(subset=["mouse_gene_name"], keep="first")
    out = out.sort_values(["human_gene_name", "mouse_gene_name"]).reset_index(drop=True)
    meta: dict[str, Any] = {
        "date": date.today(),
        "unique": unique,
        "url": url,
    }
    return JaxHumanToMouse(data=out, metadata=meta)


def _clean_jax_ids(df: pd.DataFrame) -> pd.DataFrame:
    """Strip prefixes from HGNC/MGI/OMIM ID columns."""
    for col, prefix in [
        ("hgnc_id", "HGNC:"),
        ("mouse_mgi_id", "MGI:"),
        ("omim_gene_id", "OMIM:"),
    ]:
        if col in df.columns:
            df[col] = df[col].astype(str).str.replace(prefix, "", regex=False)
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def _merge_jax_species(df: pd.DataFrame) -> pd.DataFrame:
    """Split JAX data by species and merge on db_class_key."""
    tax_col = None
    for c in df.columns:
        if "tax" in c:
            tax_col = c
            break
    if tax_col is None:
        raise RuntimeError("No taxonomy column found in JAX data.")
    hs = df[df[tax_col] == 9606].copy()
    mm = df[df[tax_col] == 10090].copy()
    hs_renames = {
        "ncbi_gene_id": "human_ncbi_gene_id",
        "gene_name": "human_gene_name",
        "hgnc_id": "human_hgnc_id",
        "omim_gene_id": "human_omim_gene_id",
    }
    hs = hs.rename(columns={k: v for k, v in hs_renames.items() if k in hs.columns})
    mm_renames = {
        "ncbi_gene_id": "mouse_ncbi_gene_id",
        "gene_name": "mouse_gene_name",
    }
    mm = mm.rename(columns={k: v for k, v in mm_renames.items() if k in mm.columns})
    merge_col = _find_merge_column(hs, mm)
    mm_cols = [merge_col]
    for c in [
        "mouse_gene_name",
        "mouse_mgi_id",
        "mouse_ncbi_gene_id",
    ]:
        if c in mm.columns:
            mm_cols.append(c)
    out = hs.merge(mm[mm_cols], on=merge_col, how="inner")
    return out.dropna(subset=["human_gene_name", "mouse_gene_name"])


def _find_merge_column(hs: pd.DataFrame, mm: pd.DataFrame) -> str:
    """Find the db_class_key merge column."""
    for c in ("db_class_key", "d_b_class_key"):
        if c in hs.columns and c in mm.columns:
            return c
    raise RuntimeError("No db_class_key column found in JAX data.")


# -------------------------------------------------------------------------
# EnsemblGenes (from a DataFrame)
# -------------------------------------------------------------------------


_NONCODING_BIOTYPES = frozenset({
    "known_ncrna",
    "lincRNA",
    "lncRNA",
    "non_coding",
})

_SMALL_BIOTYPES = frozenset({
    "miRNA",
    "misc_RNA",
    "ribozyme",
    "rRNA",
    "scaRNA",
    "scRNA",
    "snoRNA",
    "snRNA",
    "sRNA",
})

_DECAYING_BIOTYPES = frozenset({
    "non_stop_decay",
    "nonsense_mediated_decay",
})


def _apply_broad_class(
    biotype: str | None,
    chromosome: str | None,
    gene_name: str | None,
) -> str:
    """Classify a gene into a broad semantic class based on its biotype.

    Ported from R `.applyBroadClass`.
    """
    chr_str = str(chromosome) if chromosome is not None else ""
    name_str = str(gene_name) if gene_name is not None else ""
    if re.match(r"(?i)^MT", chr_str) or re.match(r"(?i)^mt[:\-]", name_str):
        return "mito"
    if biotype is None:
        return "other"
    if biotype == "protein_coding":
        return "coding"
    if biotype in _NONCODING_BIOTYPES:
        return "noncoding"
    if "pseudo" in biotype.lower():
        return "pseudo"
    if biotype in _SMALL_BIOTYPES:
        return "small"
    if biotype in _DECAYING_BIOTYPES:
        return "decaying"
    if re.match(r"(?i)^ig_", biotype):
        return "ig"
    if re.match(r"(?i)^tr_", biotype):
        return "tcr"
    return "other"


def _add_broad_class(df: pd.DataFrame) -> pd.DataFrame:
    """Add a ``broad_class`` column derived from biotype, chromosome, and gene name."""
    biotype_col = None
    for c in ("gene_biotype", "biotype"):
        if c in df.columns:
            biotype_col = c
            break
    chr_col = None
    for c in ("seqnames", "chromosome", "seq_name", "chr"):
        if c in df.columns:
            chr_col = c
            break
    name_col = "gene_name" if "gene_name" in df.columns else None
    df["broad_class"] = [
        _apply_broad_class(
            biotype=row.get(biotype_col) if biotype_col else None,
            chromosome=row.get(chr_col) if chr_col else None,
            gene_name=row.get(name_col) if name_col else None,
        )
        for _, row in df.iterrows()
    ]
    return df


def _ensembl_mysql_subdir(
    organism: str,
    release: int,
    genome_build: str,
) -> str:
    """Construct the Ensembl MySQL FTP subdirectory name."""
    slug = organism.lower().replace(" ", "_")
    suffix = re.sub(r"^[A-Za-z]+", "", re.sub(r"\.p\d+$", "", genome_build))
    return f"{slug}_core_{release}_{suffix}"


def _ensembl_ftp_gene_metadata(
    organism: str,
    genome_build: str,
    release: int,
) -> pd.DataFrame | None:
    """Fetch extra gene-level metadata from the Ensembl FTP server.

    Downloads gene descriptions, gene synonyms, and NCBI gene ID
    cross-references from the Ensembl MySQL dump files and the
    Ensembl-to-Entrez TSV mapping.

    Returns None if any required download fails.
    """
    genome_build = re.sub(r"\.p\d+$", "", genome_build)
    ftp_base = f"https://ftp.ensembl.org/pub/release-{release}"
    mysql_subdir = _ensembl_mysql_subdir(organism, release, genome_build)
    # -- gene.txt.gz ----------------------------------------------------------
    gene_url = f"{ftp_base}/mysql/{mysql_subdir}/gene.txt.gz"
    try:
        gene_path = cache_url(gene_url)
        gene = pd.read_csv(
            gene_path, sep="\t", header=None, na_values=["\\N"], quoting=3,
        )
    except Exception:
        logger.warning("Failed to download gene metadata: %s", gene_url)
        return None
    ncol = gene.shape[1]
    if ncol == 16:
        gene_col_map = {"mysql_id": 7, "gene_id": 12, "description": 9}
    elif ncol == 17:
        gene_col_map = {"mysql_id": 7, "gene_id": 13, "description": 10}
    else:
        logger.warning("Unsupported gene.txt column count: %d", ncol)
        return None
    df_gene = gene.iloc[:, list(gene_col_map.values())].copy()
    df_gene.columns = list(gene_col_map.keys())
    df_gene = df_gene.dropna().drop_duplicates()
    # -- external_synonym.txt.gz ----------------------------------------------
    syn_url = f"{ftp_base}/mysql/{mysql_subdir}/external_synonym.txt.gz"
    try:
        syn_path = cache_url(syn_url)
        synonym = pd.read_csv(
            syn_path, sep="\t", header=None, na_values=["\\N"], quoting=3,
        )
    except Exception:
        logger.warning("Failed to download synonyms: %s", syn_url)
        return None
    synonym = synonym.dropna().drop_duplicates()
    df_syn = (
        synonym.groupby(synonym.columns[0])[synonym.columns[1]]
        .apply(lambda x: "|".join(sorted(set(x.astype(str)))))
        .reset_index()
    )
    df_syn.columns = ["mysql_id", "gene_synonyms"]
    # -- Entrez TSV -----------------------------------------------------------
    org_underscore = organism.replace(" ", "_")
    entrez_url = (
        f"{ftp_base}/tsv/{org_underscore.lower()}/"
        f"{org_underscore}.{genome_build}.{release}.entrez.tsv.gz"
    )
    try:
        entrez_path = cache_url(entrez_url)
        entrez = pd.read_csv(
            entrez_path, sep="\t", na_values=["\\N"], quoting=3,
        )
    except Exception:
        logger.warning("Failed to download Entrez mapping: %s", entrez_url)
        return None
    entrez = entrez[["gene_stable_id", "xref"]].dropna().drop_duplicates()
    entrez = entrez[entrez["xref"].astype(str).str.fullmatch(r"\d+")]
    entrez["xref"] = entrez["xref"].astype(int)
    df_entrez = (
        entrez.groupby("gene_stable_id")["xref"]
        .apply(lambda x: "|".join(sorted(set(x.astype(str)))))
        .reset_index()
    )
    df_entrez.columns = ["gene_id", "ncbi_gene_id"]
    # -- Join all three -------------------------------------------------------
    out = df_gene.merge(df_syn, on="mysql_id", how="left")
    out = out.merge(df_entrez, on="gene_id", how="left")
    out = out.drop(columns=["mysql_id"])
    out = out[sorted(out.columns)]
    return out


def _ensembl_gtf_genes(
    organism: str,
    genome_build: str,
    release: int,
) -> pd.DataFrame:
    """Download an Ensembl GTF and return gene-level rows as a DataFrame."""
    genome_build = re.sub(r"\.p\d+$", "", genome_build)
    slug = organism.replace(" ", "_")
    url = (
        f"https://ftp.ensembl.org/pub/release-{release}/gtf/{slug.lower()}/"
        f"{slug}.{genome_build}.{release}.gtf.gz"
    )
    path = cache_url(url)
    rows: list[dict] = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] != "gene":
                continue
            attrs: dict[str, str] = {}
            for pair in fields[8].split(";"):
                pair = pair.strip()
                if not pair:
                    continue
                m = re.match(r'(\w+)\s+"(.+?)"', pair)
                if m:
                    attrs[m.group(1)] = m.group(2)
            rows.append(
                {
                    "end": int(fields[4]),
                    "gene_biotype": attrs.get("gene_biotype", ""),
                    "gene_id": attrs.get("gene_id", ""),
                    "gene_name": attrs.get("gene_name", ""),
                    "seqnames": fields[0],
                    "start": int(fields[3]),
                    "strand": fields[6],
                }
            )
    return pd.DataFrame(rows)


def make_ensembl_genes_from_gtf(
    organism: str,
    *,
    genome_build: str | None = None,
    release: int | None = None,
    ignore_version: bool = True,
) -> EnsemblGenes:
    """Create an EnsemblGenes object from the Ensembl GTF.

    Downloads the current Ensembl GTF for the given organism, parses
    gene-level annotations, and enriches with metadata from the Ensembl
    FTP server. This is the Python equivalent of R ``makeGRangesFromEnsembl``.

    Parameters
    ----------
    organism : str
        Latin organism name (e.g. ``'Homo sapiens'``, ``'Mus musculus'``).
    genome_build : str or None
        Ensembl genome build (e.g. ``'GRCh38'``).
        Auto-detected if ``None``.
    release : int or None
        Ensembl release version (e.g. ``115``).
        Auto-detected if ``None``.
    ignore_version : bool
        Strip version suffixes from gene identifiers.

    Returns
    -------
    EnsemblGenes

    Examples
    --------
    >>> genes = make_ensembl_genes_from_gtf("Homo sapiens")
    >>> genes = make_ensembl_genes_from_gtf("Mus musculus")
    """
    if genome_build is None:
        genome_build = current_ensembl_genome_build(organism)
    if release is None:
        release = current_ensembl_version()
    df = _ensembl_gtf_genes(
        organism=organism,
        genome_build=genome_build,
        release=release,
    )
    return make_ensembl_genes(
        df,
        organism=organism,
        genome_build=genome_build,
        release=release,
        ignore_version=ignore_version,
    )


def make_ensembl_genes(
    df: pd.DataFrame,
    organism: str | None = None,
    genome_build: str | None = None,
    release: int | None = None,
    ignore_version: bool = True,
) -> EnsemblGenes:
    """Create an EnsemblGenes object from a DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
    organism : str or None
    genome_build : str or None
    release : int or None
    ignore_version : bool

    Returns
    -------
    EnsemblGenes
    """
    df = df.copy()
    df.columns = [_to_snake(c) for c in df.columns]
    if "gene_id" not in df.columns:
        raise ValueError("DataFrame must contain a 'gene_id' column.")
    if organism is None:
        organism = detect_organism(df["gene_id"].dropna().tolist())
    if genome_build is None:
        try:
            genome_build = current_ensembl_genome_build(organism)
        except Exception:
            logger.warning("Failed to detect genome build for '%s'.", organism)
    if release is None:
        try:
            release = current_ensembl_version()
        except Exception:
            logger.warning("Failed to detect current Ensembl release.")
    needs_enrichment = not {"description", "gene_synonyms", "ncbi_gene_id"}.issubset(
        df.columns,
    )
    if needs_enrichment and genome_build is not None and release is not None:
        logger.info("Downloading extra gene-level metadata from Ensembl.")
        extra = _ensembl_ftp_gene_metadata(
            organism=organism,
            genome_build=genome_build,
            release=release,
        )
        if extra is not None:
            gene_id_col = "gene_id"
            if not ignore_version and "gene_id_no_version" in df.columns:
                gene_id_col = "gene_id_no_version"
            for col in ("description", "gene_synonyms", "ncbi_gene_id"):
                if col in df.columns:
                    extra = extra.drop(columns=[col], errors="ignore")
            df = df.merge(extra, left_on=gene_id_col, right_on="gene_id", how="left", suffixes=("", "_ftp"))
            if "gene_id_ftp" in df.columns:
                df = df.drop(columns=["gene_id_ftp"])
    if "broad_class" not in df.columns and (
        "gene_biotype" in df.columns or "biotype" in df.columns
    ):
        df = _add_broad_class(df)
    meta: dict[str, Any] = {
        "date": date.today(),
        "ignore_version": ignore_version,
        "level": "genes",
        "organism": organism,
        "provider": "Ensembl",
    }
    if genome_build is not None:
        meta["genome_build"] = genome_build
    if release is not None:
        meta["release"] = release
    return EnsemblGenes(data=df, metadata=meta)


# -------------------------------------------------------------------------
# EnsemblToNcbi / NcbiToEnsembl
# -------------------------------------------------------------------------


def make_ensembl_to_ncbi(
    df: pd.DataFrame,
    organism: str | None = None,
) -> EnsemblToNcbi:
    """Create an EnsemblToNcbi mapping from a DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
    organism : str or None

    Returns
    -------
    EnsemblToNcbi
    """
    df = df.copy()
    df.columns = [_to_snake(c) for c in df.columns]
    cols = ["ensembl_gene_id", "ncbi_gene_id"]
    for c in cols:
        if c not in df.columns:
            raise ValueError(f"Missing required column: '{c}'.")
    df = df[cols].dropna().drop_duplicates()
    df = df[~df["ensembl_gene_id"].duplicated(keep="first")]
    df = df[~df["ncbi_gene_id"].duplicated(keep="first")]
    df = df.sort_values(cols).reset_index(drop=True)
    if organism is None:
        organism = detect_organism(df["ensembl_gene_id"].tolist())
    meta: dict[str, Any] = {
        "date": date.today(),
        "organism": organism,
    }
    return EnsemblToNcbi(data=df, metadata=meta)


def make_ncbi_to_ensembl(
    df: pd.DataFrame,
    organism: str | None = None,
) -> NcbiToEnsembl:
    """Create a NcbiToEnsembl mapping from a DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
    organism : str or None

    Returns
    -------
    NcbiToEnsembl
    """
    df = df.copy()
    df.columns = [_to_snake(c) for c in df.columns]
    cols = ["ncbi_gene_id", "ensembl_gene_id"]
    for c in cols:
        if c not in df.columns:
            raise ValueError(f"Missing required column: '{c}'.")
    df = df[cols].dropna().drop_duplicates()
    df = df[~df["ncbi_gene_id"].duplicated(keep="first")]
    df = df[~df["ensembl_gene_id"].duplicated(keep="first")]
    df = df.sort_values(cols).reset_index(drop=True)
    if organism is None:
        organism = detect_organism(df["ensembl_gene_id"].tolist())
    meta: dict[str, Any] = {
        "date": date.today(),
        "organism": organism,
    }
    return NcbiToEnsembl(data=df, metadata=meta)


# -------------------------------------------------------------------------
# GeneToSymbol
# -------------------------------------------------------------------------


def make_gene_to_symbol(
    df: pd.DataFrame,
    format: str = "make_unique",
    quiet: bool = False,
) -> GeneToSymbol:
    """Create a GeneToSymbol mapping from a DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
    format : str
    quiet : bool

    Returns
    -------
    GeneToSymbol
    """
    valid_formats = ("make_unique", "1:1", "unmodified")
    if format not in valid_formats:
        raise ValueError(f"format must be one of {valid_formats}.")
    df = df.copy()
    df.columns = [_to_snake(c) for c in df.columns]
    cols = ["gene_id", "gene_name"]
    for c in cols:
        if c not in df.columns:
            raise ValueError(f"Missing required column: '{c}'.")
    df = df[cols].copy()
    df = df.drop_duplicates()
    df = _drop_incomplete_symbols(df, quiet=quiet)
    if df["gene_id"].duplicated().any():
        df = df.drop_duplicates(subset=["gene_id"], keep="first")
    df = df.sort_values(cols).reset_index(drop=True)
    df = _apply_symbol_format(df, format=format)
    df = df.sort_values(cols).reset_index(drop=True)
    meta: dict[str, Any] = {
        "date": date.today(),
        "format": format,
    }
    return GeneToSymbol(data=df, metadata=meta)


def _drop_incomplete_symbols(df: pd.DataFrame, *, quiet: bool = False) -> pd.DataFrame:
    """Drop rows with missing gene_name values."""
    complete = df["gene_name"].notna()
    if complete.all():
        return df
    n_drop = (~complete).sum()
    if not quiet:
        logger.warning(
            "Dropping %d identifier(s) without defined gene symbol.",
            n_drop,
        )
    return df[complete].copy()


def _apply_symbol_format(df: pd.DataFrame, *, format: str) -> pd.DataFrame:
    """Apply the requested symbol format."""
    if format == "make_unique":
        df["gene_name"] = _make_names_unique(df["gene_name"].fillna("unannotated").tolist())
    elif format == "1:1":
        df = df.sort_values(["gene_name", "gene_id"])
        df = df.drop_duplicates(subset=["gene_name"], keep="first")
    return df


def _make_names_unique(names: list[str]) -> list[str]:
    """Make a list of names unique by appending `.N` suffixes."""
    seen: dict[str, int] = {}
    unique: list[str] = []
    for name in names:
        if name in seen:
            seen[name] += 1
            unique.append(f"{name}.{seen[name]}")
        else:
            seen[name] = 0
            unique.append(name)
    return unique


# -------------------------------------------------------------------------
# TxToGene
# -------------------------------------------------------------------------


def make_tx_to_gene(
    df: pd.DataFrame,
    organism: str | None = None,
    genome_build: str | None = None,
    release: int | str | None = None,
    quiet: bool = False,
) -> TxToGene:
    """Create a TxToGene mapping from a DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
    organism : str or None
    genome_build : str or None
    release : int or str or None
    quiet : bool

    Returns
    -------
    TxToGene
    """
    df = df.copy()
    df.columns = [_to_snake(c) for c in df.columns]
    renames = {
        "transcript_id": "tx_id",
        "transcriptid": "tx_id",
    }
    df = df.rename(columns={k: v for k, v in renames.items() if k in df.columns})
    cols = ["tx_id", "gene_id"]
    for c in cols:
        if c not in df.columns:
            raise ValueError(f"Missing required column: '{c}'.")
    df = df[cols].copy()
    df = df.drop_duplicates()
    complete = df.notna().all(axis=1)
    if not complete.all():
        n_drop = (~complete).sum()
        if not quiet:
            logger.warning(
                "Dropping %d element(s) without transcript-to-gene mapping.",
                n_drop,
            )
        df = df[complete].copy()
    if df["tx_id"].duplicated().any():
        df = df.drop_duplicates(subset=["tx_id"], keep="first")
    df = df.sort_values(cols).reset_index(drop=True)
    meta: dict[str, Any] = {
        "date": date.today(),
    }
    if organism is not None:
        meta["organism"] = organism
    if genome_build is not None:
        meta["genome_build"] = genome_build
    if release is not None:
        meta["release"] = release
    return TxToGene(data=df, metadata=meta)
