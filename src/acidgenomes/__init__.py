"""AcidGenomes -- Toolkit for downloading and processing genome annotations."""

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
    make_ensembl_genes_from_gtf,
    make_ensembl_to_ncbi,
    make_gene_to_symbol,
    make_granges_from_ensembl,
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

# ---- Additional mapping -----------------------------------------------------
from acidgenomes._gencode_history import gencode_release_history, map_ensembl_to_gencode
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
from acidgenomes._go_terms import go_terms_per_gene_name, map_go_terms

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
from acidgenomes._protein_to_gene import make_protein_to_gene_from_ensembl

# ---- Strip versions --------------------------------------------------------
from acidgenomes._strip_versions import (
    strip_exon_versions,
    strip_gene_versions,
    strip_transcript_versions,
)
from acidgenomes._tx_to_gene_fasta import make_tx_to_gene_from_fasta

# ---- Update -----------------------------------------------------------------
from acidgenomes._update import update_gene_symbols

# ---- Genome download --------------------------------------------------------
from acidgenomes.download import (
    download_ensembl_genome,
    download_gencode_genome,
    download_refseq_genome,
    download_ucsc_genome,
)

# ---- GFF/GTF parsing --------------------------------------------------------
from acidgenomes.gff import make_granges_from_gff

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
    "download_ensembl_genome",
    "download_gencode_genome",
    "download_refseq_genome",
    "download_ucsc_genome",
    "gencode_release_history",
    "go_terms_per_gene_name",
    "import_tx_to_gene",
    "make_ensembl_genes",
    "make_ensembl_genes_from_gtf",
    "make_ensembl_to_ncbi",
    "make_gene_to_symbol",
    "make_granges_from_ensembl",
    "make_granges_from_gff",
    "make_hgnc",
    "make_jax_human_to_mouse",
    "make_mgi",
    "make_ncbi_gene_history",
    "make_ncbi_gene_info",
    "make_ncbi_to_ensembl",
    "make_protein_to_gene_from_ensembl",
    "make_tx_to_gene",
    "make_tx_to_gene_from_fasta",
    "map_ensembl_release_to_url",
    "map_ensembl_to_gencode",
    "map_gencode_to_ensembl",
    "map_gene_names_to_ensembl",
    "map_gene_names_to_hgnc",
    "map_gene_names_to_ncbi",
    "map_go_terms",
    "map_human_orthologs",
    "strip_exon_versions",
    "strip_gene_versions",
    "strip_transcript_versions",
    "update_gene_symbols",
]
