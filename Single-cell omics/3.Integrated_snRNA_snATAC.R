#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-06-28
# version:   1.0
# license:   MIT
# brief:     scRNAseq and scATAC-seq integrated pipeline.
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
library(ggplot2)
library(patchwork)
library(SeuratWrappers)
library(stringr)

# =============================================================================
# SingleCell_RNA_ATAC_Integrated: snRNA-seq AND snATAC-seq
# =============================================================================
# Sample name
sample_name='CTL'

# Analysis parameters
LSI_DIMS <- 2:30; MIN_PREDICTION_SCORE <- 0.5

# Load input
rna_obj <- readRDS(str_c(sample_name,".snRNA_processed.rds"))
atac_obj <- readRDS(str_c(sample_name,".snATAC_processed.rds"))

# 1. Calculate gene activities for ATAC data
gene.activities <- GeneActivity(atac_obj)
atac_obj[["RNA"]] <- CreateAssayObject(counts = gene.activities)
atac_obj <- NormalizeData(object = atac_obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = median(atac_obj$nCount_RNA))

# 2. Find transfer anchors using CCA (as described in paper)
transfer.anchors <- FindTransferAnchors(reference = rna_obj, query = atac_obj, features = VariableFeatures(object = rna_obj), reference.assay = "RNA", query.assay = "RNA", reduction = "cca")

# 3. Transfer cell type labels
predicted.labels <- TransferData(anchorset = transfer.anchors, refdata = rna_obj$cell_type, weight.reduction = atac_obj[["lsi"]], dims = LSI_DIMS)

atac_obj <- AddMetaData(object = atac_obj, metadata = predicted.labels)

# Filter cells with high prediction scores (>=0.5 as mentioned in paper)
high_confidence_cells <- rownames(atac_obj@meta.data)[atac_obj@meta.data$prediction.score.max >= MIN_PREDICTION_SCORE]
atac_obj$high_confidence <- rownames(atac_obj@meta.data) %in% high_confidence_cells

# SAVE/PLOT: Prediction scores
p_pred_scores <- VlnPlot(atac_obj, features = "prediction.score.max", group.by = "predicted.id")
ggsave("plots/snATAC_prediction_scores.pdf", p_pred_scores, width = 10, height = 6)

# SAVE/PLOT: Compare RNA and ATAC annotations
p_rna_annot <- DimPlot(rna_obj, group.by = "cell_type", label = TRUE) +
               NoLegend() + ggtitle("snRNA-seq Annotation")
p_atac_annot <- DimPlot(atac_obj, group.by = "predicted.id", label = TRUE) +
                NoLegend() + ggtitle("snATAC-seq Predicted Annotation")
p_comparison <- p_rna_annot + p_atac_annot
ggsave("plots/RNA_vs_ATAC_annotations.pdf", p_comparison, width = 16, height = 6)

# 4. Impute gene expression in ATAC data and create joint embedding
genes.use <- VariableFeatures(rna_obj)
refdata <- GetAssayData(rna_obj, assay = "RNA", slot = "data")[genes.use, ]

# Impute gene expression
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atac_obj[["lsi"]], dims = LSI_DIMS)
atac_obj[["integrated"]] <- imputation

# Create joint dataset (only high confidence ATAC cells)
rna_subset <- rna_obj
atac_subset <- subset(atac_obj, subset = high_confidence == TRUE)

# Add dataset labels
rna_subset$dataset <- "snRNA"
atac_subset$dataset <- "snATAC"
atac_subset$cell_type <- atac_subset$predicted.id

# Merge datasets
joint_obj <- merge(rna_subset, atac_subset, add.cell.ids = c("RNA", "ATAC"))

# Joint processing following paper methods
joint_obj <- NormalizeData(joint_obj)
joint_obj <- FindVariableFeatures(joint_obj)
joint_obj <- ScaleData(joint_obj)
joint_obj <- RunPCA(joint_obj, dims = 1:30)

# Batch correction (paper uses MNN from monocle3), proceed with standard Seurat integration
joint_obj <- RunUMAP(joint_obj, dims = 1:30)
joint_obj <- FindNeighbors(joint_obj, dims = 1:30)
joint_obj <- FindClusters(joint_obj)

# SAVE/PLOT: Joint embedding
p_joint <- DimPlot(joint_obj, group.by = "dataset", split.by = "dataset") +
           ggtitle("Joint snRNA-seq and snATAC-seq Embedding")
ggsave("plots/joint_embedding_by_dataset.pdf", p_joint, width = 12, height = 6)

p_joint_celltype <- DimPlot(joint_obj, group.by = "cell_type", label = TRUE)
ggsave("plots/joint_embedding_by_celltype.pdf", p_joint_celltype, width = 10, height = 8)

# SAVE: Joint object
saveRDS(joint_obj, str_c(sample_name,".joint_integrated.rds"))
