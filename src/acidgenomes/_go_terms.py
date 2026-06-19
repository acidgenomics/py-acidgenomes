"""Gene Ontology (GO) term fetching and gene-to-GO mapping."""

import pandas as pd

from acidgenomes._cache import fetch_json, fetch_text
from acidgenomes._mapping import map_gene_names_to_ensembl

_GO_OBO_URL = "https://purl.obolibrary.org/obo/go/go-basic.obo"


def map_go_terms() -> pd.DataFrame:
    """Download the GO OBO file and return a DataFrame of GO term ID and name.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ``go_id`` and ``go_name``, excluding obsolete terms.
    """
    text = fetch_text(_GO_OBO_URL)
    records: list[dict[str, str]] = []
    current: dict[str, str] = {}
    in_term = False

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if line == "[Term]":
            in_term = True
            current = {}
        elif line == "" and in_term:
            if current.get("id") and not current.get("is_obsolete"):
                records.append({"go_id": current["id"], "go_name": current.get("name", "")})
            in_term = False
        elif in_term:
            if line.startswith("id: "):
                current["id"] = line[4:]
            elif line.startswith("name: "):
                current["name"] = line[6:]
            elif line.startswith("is_obsolete: true"):
                current["is_obsolete"] = "true"

    if in_term and current.get("id") and not current.get("is_obsolete"):
        records.append({"go_id": current["id"], "go_name": current.get("name", "")})

    return pd.DataFrame(records)


def go_terms_per_gene_name(
    gene_names: list[str],
    organism: str,
) -> pd.DataFrame:
    """Retrieve GO terms associated with a list of gene names via Ensembl REST.

    Resolves gene names to Ensembl IDs using the Ensembl lookup endpoint,
    then fetches GO cross-references for each gene.

    Parameters
    ----------
    gene_names : list[str]
        Gene symbol names (e.g. ``["TP53", "BRCA1"]``).
    organism : str
        Latin organism name (e.g. ``"Homo sapiens"``).

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ``gene_name``, ``go_id``, ``go_aspect``, ``go_term``.
    """
    ensembl_ids = map_gene_names_to_ensembl(gene_names, organism)
    gene_id_to_name = dict(zip(ensembl_ids, gene_names, strict=False))

    records: list[dict] = []
    for gene_id, gene_name in gene_id_to_name.items():
        if not gene_id:
            continue
        url = f"https://rest.ensembl.org/xrefs/id/{gene_id}?content-type=application/json&external_db=GO"
        try:
            data = fetch_json(url)
            if isinstance(data, list):
                for item in data:
                    records.append(
                        {
                            "gene_name": gene_name,
                            "go_id": item.get("primary_id"),
                            "go_aspect": item.get("info_type"),
                            "go_term": item.get("description"),
                        }
                    )
        except Exception:
            continue

    return pd.DataFrame(records)
