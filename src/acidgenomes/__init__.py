"""AcidGenomes -- Toolkit for downloading and processing genome annotations.

This is a Python port of the R/Bioconductor AcidGenomes package.
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("acidgenomes")
except PackageNotFoundError:
    __version__ = "0.0.0.dev0"

# ---- Classes ----------------------------------------------------------------
from acidgenomes._classes import (
    EnsemblExons,
    EnsemblGenes,
    EnsemblToNcbi,
    EnsemblTranscripts,
    FlybaseExons,
    FlybaseGenes,
    FlybaseTranscripts,
    GencodeExons,
    GencodeGenes,
    GencodeTranscripts,
    GeneToSymbol,
    Hgnc,
    JaxHumanToMouse,
    Mgi,
    NcbiGeneHistory,
    NcbiGeneInfo,
    NcbiToEnsembl,
    ProteinToGene,
    RefseqExons,
    RefseqGenes,
    RefseqTranscripts,
    TxToGene,
    UcscExons,
    UcscGenes,
    UcscTranscripts,
    WormbaseExons,
    WormbaseGenes,
    WormbaseTranscripts,
)

# ---- Constructors -----------------------------------------------------------
from acidgenomes._constructors import (
    make_ensembl_genes,
    make_ensembl_to_ncbi,
    make_gene_to_symbol,
    make_hgnc,
    make_jax_human_to_mouse,
    make_mgi,
    make_ncbi_gene_history,
    make_ncbi_gene_info,
    make_ncbi_to_ensembl,
    make_tx_to_gene,
)

# ---- Detection and utilities ------------------------------------------------
from acidgenomes._detect import detect_organism
from acidgenomes._genome_build import (
    current_ensembl_genome_build,
    current_gencode_genome_build,
    current_refseq_genome_build,
    current_ucsc_genome_build,
)
from acidgenomes._genome_version import (
    current_ensembl_version,
    current_flybase_version,
    current_gencode_version,
    current_refseq_version,
    current_wormbase_version,
)

# ---- Mapping ----------------------------------------------------------------
from acidgenomes._mapping import (
    import_tx_to_gene,
    map_ensembl_release_to_url,
    map_gencode_to_ensembl,
    map_gene_names_to_ensembl,
    map_gene_names_to_hgnc,
    map_gene_names_to_ncbi,
    map_human_orthologs,
)

# ---- Strip versions --------------------------------------------------------
from acidgenomes._strip_versions import (
    strip_exon_versions,
    strip_gene_versions,
    strip_transcript_versions,
)

# ---- Update -----------------------------------------------------------------
from acidgenomes._update import update_gene_symbols

__all__ = [
    "EnsemblExons",
    "EnsemblGenes",
    "EnsemblToNcbi",
    "EnsemblTranscripts",
    "FlybaseExons",
    "FlybaseGenes",
    "FlybaseTranscripts",
    "GencodeExons",
    "GencodeGenes",
    "GencodeTranscripts",
    "GeneToSymbol",
    "Hgnc",
    "JaxHumanToMouse",
    "Mgi",
    "NcbiGeneHistory",
    "NcbiGeneInfo",
    "NcbiToEnsembl",
    "ProteinToGene",
    "RefseqExons",
    "RefseqGenes",
    "RefseqTranscripts",
    "TxToGene",
    "UcscExons",
    "UcscGenes",
    "UcscTranscripts",
    "WormbaseExons",
    "WormbaseGenes",
    "WormbaseTranscripts",
    "make_ensembl_genes",
    "make_ensembl_to_ncbi",
    "make_gene_to_symbol",
    "make_hgnc",
    "make_jax_human_to_mouse",
    "make_mgi",
    "make_ncbi_gene_history",
    "make_ncbi_gene_info",
    "make_ncbi_to_ensembl",
    "make_tx_to_gene",
    "current_ensembl_genome_build",
    "current_ensembl_version",
    "current_flybase_version",
    "current_gencode_genome_build",
    "current_gencode_version",
    "current_refseq_genome_build",
    "current_refseq_version",
    "current_ucsc_genome_build",
    "current_wormbase_version",
    "detect_organism",
    "import_tx_to_gene",
    "map_ensembl_release_to_url",
    "map_gene_names_to_ensembl",
    "map_gene_names_to_hgnc",
    "map_gene_names_to_ncbi",
    "map_gencode_to_ensembl",
    "map_human_orthologs",
    "strip_exon_versions",
    "strip_gene_versions",
    "strip_transcript_versions",
    "update_gene_symbols",
]
