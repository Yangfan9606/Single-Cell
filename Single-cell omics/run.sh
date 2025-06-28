#!/bin/bash

# =============================================================================
# Comprehensive snRNA-seq and snATAC-seq Analysis Pipeline
# Based on methods described in: 
#   Morabito, S., Miyoshi, E., Michael, N. et al. 
#   Single-nucleus chromatin accessibility and transcriptomic characterization of Alzheimer’s disease. 
#   Nat Genet 53, 1143–1155 (2021). https://doi.org/10.1038/s41588-021-00894-z
# =============================================================================

Genome_dir=./10x_reference
Transcriptome=refdata-gex-GRCh38-2024-A
scRNA_dir=./scRNAseq_fastq_files

### scRNA-seq
cellranger count --id=sample_rna \
        --transcriptome=$Genome_dir/$Transcriptome \
        --fastqs=$scRNA_dir \
        --sample=RNA1 \
        --chemistry=SC3Pv3 \
        --localcores=8 \
        --localmem=28 \
        --create-bam=true

Genome=refdata-cellranger-arc-GRCh38-2024-A
scATAC_dir=./scATACseq_fastq_files
### scATAC-seq
cellranger-atac count --id=sample_atac \
        --reference=$Genome_dir/$Genome \
        --fastqs=$scATAC_dir \
        --sample=ATAC1 \
        --localcores=8 \
        --localmem=27
