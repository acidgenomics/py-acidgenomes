# acidgenomes

Toolkit for downloading and processing genome annotations.

The package downloads and normalizes gene/transcript/exon annotations from Ensembl,
GENCODE, RefSeq, UCSC, FlyBase, and WormBase into consistent `GenomicRanges` objects,
plus reference datasets (HGNC, MGI, NCBI gene info) and gene-name mapping utilities.

## Installation

This package is hosted at [python.acidgenomics.com](https://python.acidgenomics.com/).
We recommend using [uv](https://docs.astral.sh/uv/) to install.

```sh
uv pip install \
    --index-url 'https://python.acidgenomics.com/simple/' \
    acidgenomes
```

Or add the index to your project's `pyproject.toml`:

```toml
[[tool.uv.index]]
url = "https://python.acidgenomics.com/simple/"
```

Then install:

```sh
uv add acidgenomes
```

## Detecting organisms and versions

```pycon
>>> import acidgenomes as ag
>>> ag.detect_organism(["ENSG00000000003", "ENSG00000000005"])
'Homo sapiens'
>>> ag.current_ensembl_genome_build("Homo sapiens")  # doctest: +SKIP
'GRCh38'
```

`current_ensembl_version`, `current_gencode_version`, `current_refseq_version`, and
similar functions query the corresponding provider for the current release number.

## Building genome annotations

`make_ensembl_genes_from_gtf`/`make_granges_from_ensembl` download the current
Ensembl GTF for an organism and return an `EnsemblGenes` object wrapping a
`GenomicRanges`. `download_ensembl_genome`, `download_gencode_genome`,
`download_refseq_genome`, and `download_ucsc_genome` fetch genome FASTA/GTF files
directly, and `make_granges_from_gff` parses any GFF3/GTF file with automatic
provider detection.

```pycon
>>> genes = ag.make_ensembl_genes_from_gtf("Homo sapiens")  # doctest: +SKIP
```

## ID and version utilities

```pycon
>>> ag.strip_gene_versions(["ENSG00000000003.14"])
['ENSG00000000003']
```

`strip_transcript_versions` and `strip_exon_versions` work the same way for their
respective identifier types.

## Transcript-to-gene and gene-to-symbol mappings

```pycon
>>> import pandas as pd
>>> df = pd.DataFrame({"tx_id": ["ENST001", "ENST002"], "gene_id": ["ENSG001", "ENSG002"]})
>>> t2g = ag.make_tx_to_gene(df)
>>> t2g.data
     tx_id  gene_id
0  ENST001  ENSG001
1  ENST002  ENSG002
```

```pycon
>>> df = pd.DataFrame({"gene_id": ["ENSG001"], "gene_name": ["TP53"]})
>>> g2s = ag.make_gene_to_symbol(df, format="make_unique")
>>> g2s.data
   gene_id gene_name
0  ENSG001      TP53
```

`import_tx_to_gene` builds the same mapping directly from a two-column file.

## Reference datasets

`make_hgnc`, `make_mgi`, `make_ncbi_gene_info`, `make_ncbi_gene_history`, and
`make_jax_human_to_mouse` download and normalize the corresponding reference
dataset:

```pycon
>>> hgnc = ag.make_hgnc()  # doctest: +SKIP
>>> mgi = ag.make_mgi()  # doctest: +SKIP
```

## Gene name mapping

`map_gene_names_to_hgnc`, `map_gene_names_to_ncbi`, and `map_gene_names_to_ensembl`
resolve gene symbols to the corresponding database identifier; `map_human_orthologs`
maps non-human gene IDs to their human orthologs via the Ensembl REST API;
`update_gene_symbols` resolves outdated symbols to current nomenclature.

```pycon
>>> ag.map_gene_names_to_hgnc(["TP53", "BRCA1"])  # doctest: +SKIP
[11998, 1100]
```

```{toctree}
:maxdepth: 1
:caption: Contents
:hidden:

reference/index
changelog
```
