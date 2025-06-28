#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-06-28
# version:   1.0
# license:   MIT
# brief:     scRNA-seq analysis pipeline.
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
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)

# =============================================================================
# SingleCell_RNA_seq: PARAMETER SETUP
# =============================================================================
# Sample name
output_sample_name='CTL'
# File paths - MODIFY THESE TO YOUR DATA PATHS
rna_dir <- "./sample_rna/outs/filtered_feature_bc_matrix/"

# QC Parameters
RNA_MIN_FEATURES <- 200; RNA_MAX_FEATURES <- 7500; RNA_MIN_COUNTS <- 1000; RNA_MAX_COUNTS <- 50000; RNA_MAX_MT_PERCENT <- 20

# Analysis parameters
PCA_DIMS <- 1:30; CLUSTER_RESOLUTION <- 0.5

# =============================================================================
# SingleCell_RNA_seq: snRNA-seq DATA PROCESSING
# =============================================================================

# 1. Load and create Seurat object
rna_counts <- Read10X(data.dir = rna_dir)
rna_obj <- CreateSeuratObject(counts = rna_counts, project = "snRNA", min.cells = 3, min.features = 200)

# 2. Quality control
rna_obj[["percent.mt"]] <- PercentageFeatureSet(rna_obj, pattern = "^MT-")
rna_obj[["percent.rb"]] <- PercentageFeatureSet(rna_obj, pattern = "^RP[SL]")

# SAVE/PLOT: QC visualization before filtering
p_qc_before <- VlnPlot(rna_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) +
               plot_annotation(title = "snRNA-seq QC Before Filtering")
ggsave("plots/snRNA_QC_before_filtering.pdf", p_qc_before, width = 8, height = 6)

rna_obj <- subset(rna_obj, subset = nFeature_RNA > RNA_MIN_FEATURES & nFeature_RNA < RNA_MAX_FEATURES & nCount_RNA > RNA_MIN_COUNTS & nCount_RNA < RNA_MAX_COUNTS & percent.mt < RNA_MAX_MT_PERCENT)

# 3. Normalization, scaling, and dimensionality reduction
rna_obj <- NormalizeData(rna_obj)
rna_obj <- FindVariableFeatures(rna_obj, selection.method = "vst", nfeatures = 3000)
rna_obj <- ScaleData(rna_obj)
rna_obj <- RunPCA(rna_obj, features = VariableFeatures(rna_obj), verbose = FALSE)

# SAVE/PLOT: Elbow plot for PC selection
p_elbow <- ElbowPlot(rna_obj, ndims = 50)
ggsave("plots/snRNA_elbow_plot.pdf", p_elbow, width = 8, height = 6)

rna_obj <- FindNeighbors(rna_obj, dims = PCA_DIMS)
rna_obj <- FindClusters(rna_obj, resolution = CLUSTER_RESOLUTION)
rna_obj <- RunUMAP(rna_obj, dims = PCA_DIMS)

# SAVE/PLOT: Initial clustering
p_clusters <- DimPlot(rna_obj, reduction = "umap", label = TRUE, group.by = "seurat_clusters")
ggsave("plots/snRNA_initial_clusters.pdf", p_clusters, width = 8, height = 6)

# 4. Cell type annotation using canonical markers
# https://www.nature.com/articles/s41586-019-1195-2#Sec33
cell_type_markers <- list(
    "Excitatory_Neurons" = c("NRGN", "SLC17A7", "CAMK2A", "RBFOX3"),
    "Inhibitory_Neurons" = c("GAD1", "GAD2", "GPS2", "DHFR"),
    "Astrocytes" = c("GFAP", "AQP4", "ALDH1L1"),
    "Oligodendrocytes" = c("MBP", "MOG", "PLP1"),
    "OPC" = c("VCAN", "PDGFRA", "SOX10", "CSPG4"),
    "Microglia" = c("CSF1R", "CD74", "AIF1", "P2RY12", "TMEM119"),
    "Pericytes_Endothelial" = c("CD34","PECAM1","PDGFRB")
)

# SAVE/PLOT: Marker gene expression
p_markers <- DotPlot(rna_obj, features = cell_type_markers, group.by = "seurat_clusters") +
             RotatedAxis() + ggtitle("Canonical Marker Expression")
ggsave("plots/snRNA_marker_expression.pdf", p_markers, width = 12, height = 8)

# Manual annotation (MODIFY based on 'snRNA_initial_clusters.pdf and snRNA_marker_expression.pdf')
new_cluster_ids <- case_when(
  rna_obj$seurat_clusters %in% c("9", "11") ~ "Excitatory_Neurons",
  rna_obj$seurat_clusters %in% c("5", "8", "13") ~ "Inhibitory_Neurons",
  rna_obj$seurat_clusters %in% c("4", "10", "12") ~ "Astrocytes",
  rna_obj$seurat_clusters %in% c("0", "1", "2", "3") ~ "Oligodendrocytes",
  rna_obj$seurat_clusters == "6" ~ "OPC",
  rna_obj$seurat_clusters == "7" ~ "Microglia",
  TRUE ~ as.character(rna_obj$seurat_clusters)
)
rna_obj$cell_type <- new_cluster_ids
Idents(rna_obj) <- rna_obj$cell_type

# SAVE/PLOT: Annotated cell types
p_annotated <- DimPlot(rna_obj, reduction = "umap", label = TRUE, group.by = "cell_type")
ggsave("plots/snRNA_annotated_celltypes.pdf", p_annotated, width = 10, height = 8)

# SAVE: snRNA-seq object
saveRDS(rna_obj, str_c(output_sample_name,".snRNA_processed.rds"))
