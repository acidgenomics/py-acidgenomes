"""Unit tests for AcidGenomes.

Tests that require no network access.
"""

from __future__ import annotations

import pandas as pd
import pytest

from acidgenomes._detect import _detect_single, detect_organism


class TestDetectOrganism:
    def test_human_ensg(self) -> None:
        assert _detect_single("ENSG00000000003") == "Homo sapiens"

    def test_mouse_ensmusg(self) -> None:
        assert _detect_single("ENSMUSG00000000001") == "Mus musculus"

    def test_fly_fbgn(self) -> None:
        assert _detect_single("FBgn0000003") == "Drosophila melanogaster"

    def test_worm_wbgene(self) -> None:
        assert _detect_single("WBGene00000001") == "Caenorhabditis elegans"

    def test_zebrafish_ensdarg(self) -> None:
        assert _detect_single("ENSDARG00000000001") == "Danio rerio"

    def test_unknown_returns_none(self) -> None:
        assert _detect_single("FOO12345") is None

    def test_detect_organism_list(self) -> None:
        ids = ["ENSG00000000003", "ENSG00000000005"]
        assert detect_organism(ids) == "Homo sapiens"

    def test_detect_organism_empty_raises(self) -> None:
        with pytest.raises(ValueError):
            detect_organism([])

    def test_detect_organism_mixed_raises(self) -> None:
        with pytest.raises(ValueError):
            detect_organism(["ENSG00000000003", "ENSMUSG00000000001"])


from acidgenomes._strip_versions import (
    strip_exon_versions,
    strip_gene_versions,
    strip_transcript_versions,
)


class TestStripVersions:
    def test_strip_gene_versions_ensembl(self) -> None:
        result = strip_gene_versions(["ENSG00000000003.14", "ENSG00000000005"])
        assert result == ["ENSG00000000003", "ENSG00000000005"]

    def test_strip_transcript_versions(self) -> None:
        result = strip_transcript_versions(["ENST00000000233.10", "ENST00000000412"])
        assert result == ["ENST00000000233", "ENST00000000412"]

    def test_strip_exon_versions(self) -> None:
        result = strip_exon_versions(["ENSE00000000001.2", "ENSE00000000002"])
        assert result == ["ENSE00000000001", "ENSE00000000002"]

    def test_strip_non_versioned_unchanged(self) -> None:
        ids = ["GENE1", "GENE2"]
        assert strip_gene_versions(ids) == ids


from acidgenomes._classes import (
    EnsemblGenes,
    EnsemblToNcbi,
    GeneToSymbol,
    Hgnc,
    NcbiGeneInfo,
    TxToGene,
)


class TestAnnotatedDataFrame:
    def test_hgnc_creation(self) -> None:
        df = pd.DataFrame({"hgnc_id": [1, 2], "gene_name": ["A", "B"]})
        obj = Hgnc(data=df, metadata={"organism": "Homo sapiens"})
        assert len(obj) == 2
        assert obj.metadata["organism"] == "Homo sapiens"

    def test_ncbi_gene_info(self) -> None:
        df = pd.DataFrame({"gene_id": [1], "gene_name": ["X"]})
        obj = NcbiGeneInfo(data=df, metadata={})
        assert len(obj) == 1

    def test_ensembl_genes(self) -> None:
        df = pd.DataFrame({"gene_id": ["ENSG001"]})
        obj = EnsemblGenes(data=df, metadata={"organism": "Homo sapiens"})
        assert obj.organism == "Homo sapiens"

    def test_ensembl_to_ncbi(self) -> None:
        df = pd.DataFrame({"ensembl_gene_id": ["ENSG001"], "ncbi_gene_id": [1]})
        obj = EnsemblToNcbi(data=df, metadata={})
        assert len(obj) == 1

    def test_gene_to_symbol(self) -> None:
        df = pd.DataFrame({"gene_id": ["G1"], "gene_name": ["TP53"]})
        obj = GeneToSymbol(data=df, metadata={"format": "make_unique"})
        assert len(obj) == 1

    def test_tx_to_gene(self) -> None:
        df = pd.DataFrame({"tx_id": ["T1"], "gene_id": ["G1"]})
        obj = TxToGene(data=df, metadata={})
        assert len(obj) == 1


from acidgenomes._constructors import (
    make_ensembl_genes,
    make_ensembl_to_ncbi,
    make_gene_to_symbol,
    make_ncbi_to_ensembl,
    make_tx_to_gene,
)


class TestConstructors:
    def test_make_ensembl_genes(self) -> None:
        df = pd.DataFrame({"gene_id": ["ENSG00000000003", "ENSG00000000005"]})
        obj = make_ensembl_genes(df, organism="Homo sapiens")
        assert isinstance(obj, EnsemblGenes)
        assert len(obj) == 2
        assert obj.organism == "Homo sapiens"

    def test_make_ensembl_to_ncbi(self) -> None:
        df = pd.DataFrame({"ensembl_gene_id": ["ENSG001"], "ncbi_gene_id": [1]})
        obj = make_ensembl_to_ncbi(df, organism="Homo sapiens")
        assert isinstance(obj, EnsemblToNcbi)

    def test_make_ncbi_to_ensembl(self) -> None:
        df = pd.DataFrame({"ncbi_gene_id": [1], "ensembl_gene_id": ["ENSG001"]})
        from acidgenomes._classes import NcbiToEnsembl

        obj = make_ncbi_to_ensembl(df, organism="Homo sapiens")
        assert isinstance(obj, NcbiToEnsembl)

    def test_make_gene_to_symbol_unique(self) -> None:
        df = pd.DataFrame(
            {
                "gene_id": ["G1", "G2", "G3"],
                "gene_name": ["TP53", "TP53", "BRCA1"],
            }
        )
        obj = make_gene_to_symbol(df, format="make_unique")
        assert isinstance(obj, GeneToSymbol)
        names = obj.data["gene_name"].tolist()
        assert len(names) == len(set(names)), "Names should be unique"

    def test_make_gene_to_symbol_1to1(self) -> None:
        df = pd.DataFrame(
            {
                "gene_id": ["G1", "G2", "G3"],
                "gene_name": ["TP53", "TP53", "BRCA1"],
            }
        )
        obj = make_gene_to_symbol(df, format="1:1")
        assert isinstance(obj, GeneToSymbol)
        assert obj.data["gene_name"].value_counts().max() == 1

    def test_make_tx_to_gene(self) -> None:
        df = pd.DataFrame({"tx_id": ["T1", "T2"], "gene_id": ["G1", "G2"]})
        obj = make_tx_to_gene(df)
        assert isinstance(obj, TxToGene)
        assert len(obj) == 2

    def test_make_tx_to_gene_deduplicates(self) -> None:
        df = pd.DataFrame({"tx_id": ["T1", "T1"], "gene_id": ["G1", "G1"]})
        obj = make_tx_to_gene(df)
        assert len(obj) == 1

    def test_make_tx_to_gene_missing_col(self) -> None:
        df = pd.DataFrame({"foo": [1]})
        with pytest.raises(ValueError, match="Missing required column"):
            make_tx_to_gene(df)


from acidgenomes._data import (
    DETECT_ORGANISM_DATA,
    NCBI_TAX_IDS,
    NCBI_TAXONOMIC_GROUPS,
)


class TestData:
    def test_detect_organism_data_has_entries(self) -> None:
        assert len(DETECT_ORGANISM_DATA) >= 9

    def test_ncbi_tax_ids_homo(self) -> None:
        assert NCBI_TAX_IDS["Homo sapiens"] == 9606

    def test_ncbi_tax_ids_mouse(self) -> None:
        assert NCBI_TAX_IDS["Mus musculus"] == 10090

    def test_taxonomic_groups_homo(self) -> None:
        assert NCBI_TAXONOMIC_GROUPS["Homo sapiens"]["gene_info"] == "Mammalia"


from acidgenomes._cache import get_cache_dir


class TestCache:
    def test_get_cache_dir_returns_path(self) -> None:
        p = get_cache_dir()
        assert "acidgenomes" in str(p)
