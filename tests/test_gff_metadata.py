"""Tests for GFF metadata extraction and provider detection."""

from pathlib import Path

from acidgenomes.gff._metadata import get_gff_metadata
from acidgenomes.gff._providers import GffFormat, Provider

FIXTURES = Path(__file__).parent / "fixtures"


def test_ensembl_gtf_provider():
    meta = get_gff_metadata(FIXTURES / "ensembl_genes.gtf")
    assert meta.provider == Provider.ENSEMBL


def test_ensembl_gtf_format():
    meta = get_gff_metadata(FIXTURES / "ensembl_genes.gtf")
    assert meta.format == GffFormat.GTF


def test_ensembl_gtf_genome_build():
    meta = get_gff_metadata(FIXTURES / "ensembl_genes.gtf")
    assert meta.genome_build == "GRCh38"


def test_gencode_gtf_provider():
    meta = get_gff_metadata(FIXTURES / "gencode_genes.gtf")
    assert meta.provider == Provider.GENCODE


def test_flybase_gtf_provider():
    meta = get_gff_metadata(FIXTURES / "flybase_genes.gtf")
    assert meta.provider == Provider.FLYBASE


def test_wormbase_gtf_provider():
    meta = get_gff_metadata(FIXTURES / "wormbase_genes.gtf")
    assert meta.provider == Provider.WORMBASE


def test_wormbase_release():
    meta = get_gff_metadata(FIXTURES / "wormbase_genes.gtf")
    assert meta.release == "WS289"
