"""Protein-to-gene identifier mapping via Ensembl REST API."""

import pandas as pd
import requests

from acidgenomes._classes import ProteinToGene

_ENSEMBL_LOOKUP_URL = "https://rest.ensembl.org/lookup/id"


def make_protein_to_gene_from_ensembl(
    protein_ids: list[str],
    *,
    batch_size: int = 1000,
) -> ProteinToGene:
    """Map Ensembl protein IDs (ENSP) to gene IDs and gene names.

    Uses the Ensembl REST API ``/lookup/id`` endpoint with POST requests
    to batch-resolve protein identifiers.

    Parameters
    ----------
    protein_ids : list[str]
        Ensembl protein identifiers, e.g. ``["ENSP00000000233", ...]``.
    batch_size : int
        Number of IDs per POST request. Maximum allowed by Ensembl is 1000.

    Returns
    -------
    ProteinToGene
        DataFrame with columns ``protein_id``, ``gene_id``, ``gene_name``.
    """
    records: list[dict] = []
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    for i in range(0, len(protein_ids), batch_size):
        batch = protein_ids[i : i + batch_size]
        payload = {"ids": batch}
        resp = requests.post(
            f"{_ENSEMBL_LOOKUP_URL}?content-type=application/json",
            json=payload,
            headers=headers,
            timeout=60,
        )
        resp.raise_for_status()
        data = resp.json()
        for protein_id, info in data.items():
            if info is None:
                continue
            gene_id = info.get("Parent") or info.get("parent")
            gene_name = info.get("display_name") or info.get("gene_name") or gene_id
            records.append({
                "protein_id": protein_id,
                "gene_id": gene_id,
                "gene_name": gene_name,
            })

    df = pd.DataFrame(records, columns=["protein_id", "gene_id", "gene_name"])
    return ProteinToGene(data=df)
