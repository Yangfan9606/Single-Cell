#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-06-28
# version:   1.0
# license:   MIT
# brief:     scATAC-seq analysis pipeline.
#------------------------------------------#

# =============================================================================
# Comprehensive snRNA-seq and snATAC-seq Analysis Pipeline
# Based on methods described in: 
#   Morabito, S., Miyoshi, E., Michael, N. et al. 
#   Single-nucleus chromatin accessibility and transcriptomic characterization of Alzheimer’s disease. 
#   Nat Genet 53, 1143–1155 (2021). https://doi.org/10.1038/s41588-021-00894-z
# =============================================================================

# Load required libraries
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(stringr)

# =============================================================================
# SingleCell_ATAC_seq: PARAMETER SETUP
# =============================================================================
# Sample name
output_sample_name='CTL'

# File paths - MODIFY THESE TO YOUR DATA PATHS
atac_h5_path <- "./sample_atac/outs/filtered_peak_bc_matrix.h5"
atac_meta_path <- "./sample_atac/outs/singlecell.csv"
## 	zcat fragments.tsv.gz|cut -f1-5 | bgzip > corrected_fragments.tsv.gz
## 	tabix -p vcf corrected_fragments.tsv.gz
atac_frag_path <- "./sample_atac/outs/corrected_fragments.tsv.gz"

# QC Parameters
ATAC_MIN_FRAGMENTS <- 1000; ATAC_MAX_FRAGMENTS <- 100000; ATAC_MIN_PEAK_PCT <- 15; ATAC_MAX_NUCSIG <- 2; ATAC_MIN_TSS <- 2

# Analysis parameters
LSI_DIMS <- 2:30; CLUSTER_RESOLUTION <- 0.5

# =============================================================================
# SingleCell_ATAC_seq: snATAC-seq DATA PROCESSING
# =============================================================================

# 1. Load and create Seurat object
atac_counts <- Read10X_h5(filename = atac_h5_path)
atac_metadata <- read.csv(file = atac_meta_path, header = TRUE, row.names = 1)

# Get gene annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
chrom_assay <- CreateChromatinAssay(counts = atac_counts, sep = c(":", "-"), genome = "hg38", fragments = atac_frag_path, min.cells = 10, min.features = 200, annotation = annotations)
atac_obj <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC", meta.data = atac_metadata)

# 2. Quality control metrics
atac_obj <- NucleosomeSignal(object = atac_obj)
atac_obj <- TSSEnrichment(object = atac_obj, fast = FALSE)
atac_obj$pct_reads_in_peaks <- atac_obj$peak_region_fragments / atac_obj$passed_filters * 100

# SAVE/PLOT: ATAC QC metrics
p_atac_qc <- VlnPlot(atac_obj, features = c("peak_region_fragments", "TSS.enrichment", "nucleosome_signal", "pct_reads_in_peaks"), ncol = 4, pt.size = 0)
ggsave("plots/snATAC_QC_metrics.pdf", p_atac_qc, width = 16, height = 4)

atac_obj <- subset(x = atac_obj, subset = peak_region_fragments > ATAC_MIN_FRAGMENTS & peak_region_fragments < ATAC_MAX_FRAGMENTS & pct_reads_in_peaks > ATAC_MIN_PEAK_PCT & nucleosome_signal < ATAC_MAX_NUCSIG & TSS.enrichment > ATAC_MIN_TSS)

# 3. Normalization and dimensionality reduction
atac_obj <- RunTFIDF(atac_obj)
atac_obj <- FindTopFeatures(atac_obj, min.cutoff = "q0")
atac_obj <- RunSVD(atac_obj)

# SAVE/PLOT: LSI depth correlation
p_depth_cor <- DepthCor(atac_obj, n = 10)
ggsave("plots/snATAC_LSI_depth_correlation.pdf", p_depth_cor, width = 10, height = 6)

atac_obj <- RunUMAP(object = atac_obj, reduction = "lsi", dims = LSI_DIMS)
atac_obj <- FindNeighbors(object = atac_obj, reduction = "lsi", dims = LSI_DIMS)
atac_obj <- FindClusters(object = atac_obj, verbose = FALSE, algorithm = 3, resolution = CLUSTER_RESOLUTION)

# SAVE/PLOT: ATAC clustering
p_atac_clusters <- DimPlot(object = atac_obj, label = TRUE) + NoLegend() +
                   ggtitle("snATAC-seq Initial Clusters")
ggsave("plots/snATAC_initial_clusters.pdf", p_atac_clusters, width = 8, height = 6)

# SAVE: snATAC-seq object
saveRDS(atac_obj, str_c(output_sample_name,".snATAC_processed.rds"))
