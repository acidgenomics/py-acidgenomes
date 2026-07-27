# acidgenomes

Toolkit for downloading and processing genome annotations.

## Installation

This is a [Python][] package hosted at [python.acidgenomics.com][].
We recommend using [uv][] to install.

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

## Quick start

```python
import acidgenomes as ag

# Detect organism from Ensembl gene IDs
ag.detect_organism(["ENSG00000000003", "ENSG00000000005"])
# => "Homo sapiens"

# Current Ensembl genome build
ag.current_ensembl_genome_build("Homo sapiens")
# => "GRCh38"

# Strip version suffixes
ag.strip_gene_versions(["ENSG00000000003.14"])
# => ["ENSG00000000003"]

# Build a transcript-to-gene mapping
import pandas as pd
df = pd.DataFrame({"tx_id": ["ENST001", "ENST002"], "gene_id": ["ENSG001", "ENSG002"]})
t2g = ag.make_tx_to_gene(df)

# Build a gene-to-symbol mapping
df = pd.DataFrame({"gene_id": ["ENSG001"], "gene_name": ["TP53"]})
g2s = ag.make_gene_to_symbol(df, format="make_unique")
```

### Downloading reference data

```python
# HGNC (human gene nomenclature)
hgnc = ag.make_hgnc()

# MGI (mouse)
mgi = ag.make_mgi()

# NCBI gene info
ncbi = ag.make_ncbi_gene_info("Homo sapiens")

# JAX human-to-mouse orthologs
jax = ag.make_jax_human_to_mouse()

# Map gene names to HGNC IDs
ag.map_gene_names_to_hgnc(["TP53", "BRCA1"])

# Update outdated gene symbols
ag.update_gene_symbols(["ZCCHC11"], organism="Homo sapiens")
```

[python]: https://www.python.org/
[python.acidgenomics.com]: https://python.acidgenomics.com
[uv]: https://docs.astral.sh/uv/


## License

Apache-2.0 — Copyright 2026 Acid Genomics LLC — see [LICENSE](LICENSE).
