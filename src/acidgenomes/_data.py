"""Internal data constants for AcidGenomes.

These are ported from the R package's ``data-raw/`` directory:
``detectOrganism.csv``, ``mapNcbiTaxId.csv``, and ``sysdata.R``.
"""

from __future__ import annotations

import re
from typing import TypedDict


class OrganismEntry(TypedDict):
    organism: str
    pattern: re.Pattern[str]


# Ported from data-raw/detectOrganism.csv.
DETECT_ORGANISM_DATA: list[OrganismEntry] = [
    {"organism": "Homo sapiens", "pattern": re.compile(r"^ENSG\d{11}")},
    {"organism": "Mus musculus", "pattern": re.compile(r"^ENSMUSG\d{11}")},
    {"organism": "Rattus norvegicus", "pattern": re.compile(r"^ENSRNOG\d{11}")},
    {"organism": "Danio rerio", "pattern": re.compile(r"^ENSDARG\d{11}")},
    {"organism": "Drosophila melanogaster", "pattern": re.compile(r"^FBgn\d{7}")},
    {"organism": "Caenorhabditis elegans", "pattern": re.compile(r"^WBGene\d{8}")},
    {"organism": "Gallus gallus", "pattern": re.compile(r"^ENSGALG\d{11}")},
    {"organism": "Ovis aries", "pattern": re.compile(r"^ENSOARG\d{11}")},
    {"organism": "Sus scrofa", "pattern": re.compile(r"^ENSSSCG\d{11}")},
]

# Ported from data-raw/mapNcbiTaxId.csv.
NCBI_TAX_IDS: dict[str, int] = {
    "Homo sapiens": 9606,
    "Mus musculus": 10090,
    "Rattus norvegicus": 10116,
    "Danio rerio": 7955,
    "Drosophila melanogaster": 7227,
    "Caenorhabditis elegans": 6239,
    "Gallus gallus": 9031,
    "Ovis aries": 9940,
    "Sus scrofa": 9823,
}

# Genome build → UCSC genome identifier mapping.
GENOME_BUILD_TO_UCSC: dict[str, str] = {
    "GRCh38": "hg38",
    "GRCh37": "hg19",
    "GRCm39": "mm39",
    "GRCm38": "mm10",
}

# NCBI FTP taxonomic group mapping.
NCBI_TAXONOMIC_GROUPS: dict[str, dict[str, str]] = {
    "Homo sapiens": {"gene_info": "Mammalia", "refseq": "vertebrate_mammalian"},
    "Mus musculus": {"gene_info": "Mammalia", "refseq": "vertebrate_mammalian"},
    "Rattus norvegicus": {"gene_info": "Mammalia", "refseq": "vertebrate_mammalian"},
    "Danio rerio": {
        "gene_info": "Non-mammalian_vertebrates",
        "refseq": "vertebrate_other",
    },
    "Drosophila melanogaster": {"gene_info": "Invertebrates", "refseq": "invertebrate"},
    "Caenorhabditis elegans": {"gene_info": "Invertebrates", "refseq": "invertebrate"},
    "Gallus gallus": {
        "gene_info": "Non-mammalian_vertebrates",
        "refseq": "vertebrate_other",
    },
    "Ovis aries": {"gene_info": "Mammalia", "refseq": "vertebrate_mammalian"},
    "Sus scrofa": {"gene_info": "Mammalia", "refseq": "vertebrate_mammalian"},
}

# GRanges level ordering (for consistent column ordering).
GRANGES_LEVELS: list[str] = [
    "genes",
    "transcripts",
    "exons",
    "cds",
    "5utr",
    "3utr",
]

# Base URL for AcidGenomes test data.
ACIDGENOMES_TESTS_URL: str = "https://r.acidgenomics.com/testdata/acidgenomes"
