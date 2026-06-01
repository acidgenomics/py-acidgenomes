"""Genome download functions."""

from acidgenomes.download._ensembl import download_ensembl_genome
from acidgenomes.download._gencode import download_gencode_genome
from acidgenomes.download._refseq import download_refseq_genome
from acidgenomes.download._ucsc import download_ucsc_genome

__all__ = [
    "download_ensembl_genome",
    "download_gencode_genome",
    "download_refseq_genome",
    "download_ucsc_genome",
]
