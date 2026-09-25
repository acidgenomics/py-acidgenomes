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
        with pytest.raises(ValueError, match="must not be empty"):
            detect_organism([])

    def test_detect_organism_mixed_raises(self) -> None:
        with pytest.raises(ValueError, match="Multiple organisms"):
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
    _apply_broad_class,
    _merge_jax_species,
    make_ensembl_genes,
    make_ensembl_to_ncbi,
    make_gene_to_symbol,
    make_ncbi_to_ensembl,
    make_tx_to_gene,
)


class TestConstructors:
    def test_merge_jax_species_no_column_collision_suffix(self) -> None:
        # The raw JAX report stacks one column per field across both
        # species' rows (e.g. "mouse_mgi_id" appears, always NaN, on human
        # rows too). Before the fix, merging on db_class_key without first
        # dropping that meaningless human-side copy let pandas silently
        # rename the real value to "mouse_mgi_id_y" instead of the plain
        # name every downstream consumer expects.
        df = pd.DataFrame(
            {
                "db_class_key": [1, 1],
                "tax_id": [9606, 10090],
                "gene_name": ["HBB", "Hbb-bs"],
                "ncbi_gene_id": [3043, 15129],
                "hgnc_id": ["HGNC:4827", pd.NA],
                "omim_gene_id": [141900, pd.NA],
                "mouse_mgi_id": [pd.NA, "MGI:96021"],
            }
        )
        out = _merge_jax_species(df)
        assert list(out.columns) == [c for c in out.columns if not c.endswith(("_x", "_y"))]
        assert out.loc[0, "mouse_mgi_id"] == "MGI:96021"
        assert out.loc[0, "human_gene_name"] == "HBB"
        assert out.loc[0, "mouse_gene_name"] == "Hbb-bs"

    def test_apply_broad_class_nan_biotype(self) -> None:
        # Ensembl release 116 GTF files can yield a float NaN biotype instead
        # of None or a str, which previously raised AttributeError on .lower().
        result = _apply_broad_class(biotype=float("nan"), chromosome="1", gene_name="GENE1")
        assert result == "other"

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


from acidgenomes._classes import JaxHumanToMouse, Mgi
from acidgenomes._mapping import (
    _mouse_ensembl_id_by_mgi_id,
    _mouse_mgi_ids_by_human_hgnc_id,
    _tags_by_hgnc_id,
    _tags_for_group_id_string,
    classify_curated_gene_groups,
)

# Real HGNC hgnc_complete_set.txt values, verified live 2026-09-24 (see
# classify_curated_gene_groups's own docstring for the group-ID map this
# exercises): HBB (hemoglobin, protein-coding), HBAP1 (hemoglobin,
# pseudogene -- HGNC's own group includes it), RPS4X (cytoplasmic
# ribosomal), MRPL1 (mitochondrial ribosomal, dual-tagged 646|1994), and
# RPS6KA1 (S6 kinase family -- must NOT be tagged ribosomal).
_REAL_HGNC_ROWS = pd.DataFrame(
    {
        "hgnc_id": [4827, 4825, 10424, 14275, 10430],
        "gene_name": ["HBB", "HBAP1", "RPS4X", "MRPL1", "RPS6KA1"],
        "gene_group_id": ["940", "940", "728", "646|1994", "1156|1691|3524"],
        "ensembl_gene_id": [
            "ENSG00000244734",
            "ENSG00000225323",
            "ENSG00000198034",
            "ENSG00000169288",
            "ENSG00000117676",
        ],
    }
)


class TestTagsForGroupIdString:
    def test_ribosomal_cyto(self) -> None:
        assert _tags_for_group_id_string("728") == ["ribo_cyto"]

    def test_ribosomal_mito_dedup(self) -> None:
        # 646 and 1994 both map to "ribo_mito" -- must not appear twice.
        assert _tags_for_group_id_string("646|1994") == ["ribo_mito"]

    def test_hemoglobin(self) -> None:
        assert _tags_for_group_id_string("940") == ["hemoglobin"]

    def test_s6_kinase_family_excluded(self) -> None:
        assert _tags_for_group_id_string("1156|1691|3524") == []

    def test_unknown_group_id_ignored(self) -> None:
        assert _tags_for_group_id_string("940|99999") == ["hemoglobin"]

    def test_empty_string(self) -> None:
        assert _tags_for_group_id_string("") == []


class TestTagsByHgncId:
    def test_named_gene_memberships(self) -> None:
        result = _tags_by_hgnc_id(_REAL_HGNC_ROWS)
        assert result[4827] == ["hemoglobin"]  # HBB
        assert result[4825] == ["hemoglobin"]  # HBAP1, a pseudogene
        assert result[10424] == ["ribo_cyto"]  # RPS4X
        assert result[14275] == ["ribo_mito"]  # MRPL1
        assert 10430 not in result  # RPS6KA1: no curated tag, not stored

    def test_missing_column_raises(self) -> None:
        df = pd.DataFrame({"hgnc_id": [1]})
        with pytest.raises(ValueError, match="gene_group_id"):
            _tags_by_hgnc_id(df)


class TestMouseOrthologPropagation:
    def test_multiple_mouse_paralogs_all_collected(self) -> None:
        # A human gene (HBB) with three mouse paralogs (the real Hbb-bs/
        # Hbb-bt/Hbb-y expansion has more; three is enough to prove no
        # first/last-wins collapse).
        jax_df = pd.DataFrame(
            {
                "human_hgnc_id": [4827, 4827, 4827],
                "mouse_mgi_id": [96040, 96041, 96043],
            }
        )
        result = _mouse_mgi_ids_by_human_hgnc_id(jax_df)
        assert result[4827] == {96040, 96041, 96043}

    def test_mouse_ensembl_id_by_mgi_id(self) -> None:
        mgi_df = pd.DataFrame(
            {
                "mgi_accession_id": [96040, 96041],
                "ensembl_gene_id": ["ENSMUSG00000052305", "ENSMUSG00000045791"],
            }
        )
        result = _mouse_ensembl_id_by_mgi_id(mgi_df)
        assert result[96040] == "ENSMUSG00000052305"
        assert result[96041] == "ENSMUSG00000045791"

    def test_missing_column_raises(self) -> None:
        with pytest.raises(ValueError, match="mouse_mgi_id"):
            _mouse_mgi_ids_by_human_hgnc_id(pd.DataFrame({"human_hgnc_id": [1]}))


class TestClassifyCuratedGeneGroups:
    def test_human_named_assertions(self) -> None:
        hgnc = Hgnc(data=_REAL_HGNC_ROWS, metadata={"organism": "Homo sapiens"})
        result = classify_curated_gene_groups(
            [
                "ENSG00000244734",  # HBB
                "ENSG00000225323",  # HBAP1 (pseudogene)
                "ENSG00000198034",  # RPS4X
                "ENSG00000169288",  # MRPL1
                "ENSG00000117676",  # RPS6KA1
                "ENSG00000000000",  # not in the table at all
            ],
            "Homo sapiens",
            hgnc=hgnc,
        )
        assert result["ENSG00000244734"] == ["hemoglobin"]
        assert result["ENSG00000225323"] == ["hemoglobin"]
        assert result["ENSG00000198034"] == ["ribo_cyto"]
        assert result["ENSG00000169288"] == ["ribo_mito"]
        assert result["ENSG00000117676"] == []
        assert result["ENSG00000000000"] == []

    def test_mouse_propagation_end_to_end(self) -> None:
        hgnc = Hgnc(
            data=pd.DataFrame(
                {
                    "hgnc_id": [4827],
                    "gene_group_id": ["940"],
                    "ensembl_gene_id": ["ENSG00000244734"],
                }
            ),
            metadata={"organism": "Homo sapiens"},
        )
        jax = JaxHumanToMouse(
            data=pd.DataFrame(
                {
                    "human_hgnc_id": [4827, 4827],
                    "mouse_mgi_id": [96040, 96041],
                }
            ),
            metadata={},
        )
        mgi = Mgi(
            data=pd.DataFrame(
                {
                    "mgi_accession_id": [96040, 96041, 999999],
                    "ensembl_gene_id": [
                        "ENSMUSG00000052305",
                        "ENSMUSG00000045791",
                        "ENSMUSG00000000000",
                    ],
                }
            ),
            metadata={"organism": "Mus musculus"},
        )
        result = classify_curated_gene_groups(
            ["ENSMUSG00000052305", "ENSMUSG00000045791", "ENSMUSG00000000000"],
            "Mus musculus",
            hgnc=hgnc,
            jax=jax,
            mgi=mgi,
        )
        assert result["ENSMUSG00000052305"] == ["hemoglobin"]
        assert result["ENSMUSG00000045791"] == ["hemoglobin"]
        # 999999 was never reached by any ortholog edge -- must stay untagged.
        assert result["ENSMUSG00000000000"] == []

    def test_unsupported_organism_raises(self) -> None:
        with pytest.raises(ValueError, match="Unsupported organism"):
            classify_curated_gene_groups(["ENSG001"], "Danio rerio")
