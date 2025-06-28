# Single-Cell Multi-Omics Analysis Pipeline

This repository contains a comprehensive pipeline for integrated analysis of single-nucleus RNA-seq (snRNA-seq) and single-nucleus ATAC-seq (snATAC-seq) data, based on methods described in [Morabito et al., Nat Genet 53, 1143–1155 (2021)](https://doi.org/10.1038/s41588-021-00894-z).

## Overview
The pipeline consists of multiple scripts:

`run.sh` - **Data preprocessing** with Cell Ranger
`1.snRNA.R` - **snRNA-seq data processing** and cell annotation
`2.snATAC.R` - **snATAC-seq data processing**
`3.Integraged_snRNA_snATAC.R` - **Integrated analysis** of both modalities
`4.Integraged_Cicero_CoAccessibility.R` - **Chromatin co-accessibility analysis** with Cicero

## Prerequisites
- 10x Genomics Cell Ranger
- R (v4.4.3 or higher)
- Required R packages: Seurat, Signac, Cicero, chromVAR, and others listed in the scripts

## Step-by-Step Guide

### Step 0: Data Preprocessing (Cell Ranger)
Convert raw FASTQ files to expression/accessibility matrices:

1. Download reference genomes from [10x Genomics](https://www.10xgenomics.com/support/software/downloads)
2. Edit `run.sh` with your paths.

**Outputs:**
- `sample_rna/outs/` - RNA expression data
- `sample_atac/outs/` - ATAC accessibility data

### Step 1: snRNA-seq Data Analysis
Process snRNA-seq data and perform cell annotation:
Key operations:
*Quality control filtering
*Normalization and dimensionality reduction
*Cell clustering using Seurat
*Cell type annotation with canonical markers

**Outputs:**
- `sample.snRNA_processed.rds` - Annotated Seurat object
- `plots/snRNA_QC_before_filtering.pdf` - Pre-filtering QC metrics
- `plots/snRNA_elbow_plot.pdf` - PCA elbow plot
- `plots/snRNA_initial_clusters.pdf` - Initial clustering results
- `plots/snRNA_marker_expression.pdf` - Marker gene expression
- `plots/snRNA_annotated_celltypes.pdf` - Final annotated cell types
