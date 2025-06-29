#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-06-29
# version:   1.0
# license:   MIT
# brief:     TF motif and ChromVAR analysis.
#------------------------------------------#

# =============================================================================
# Comprehensive snRNA-seq and snATAC-seq Analysis Pipeline
# Based on methods described in: 
#   Morabito, S., Miyoshi, E., Michael, N. et al. 
#   Single-nucleus chromatin accessibility and transcriptomic characterization of Alzheimer’s disease. 
#   Nat Genet 53, 1143–1155 (2021). https://doi.org/10.1038/s41588-021-00894-z
# =============================================================================

# Load required libraries
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(GenomicRanges)
library(JASPAR2024)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVAR)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(colorRamp2)
library(RColorBrewer)
library(viridis)

# Output dir
output_dirs <- c("results/motif", "plots/motif")
sapply(output_dirs, function(x) if(!dir.exists(x)) dir.create(x, recursive = TRUE))

# Parameters
input_file = "CTL.joint_integrated.rds"
min_motif_ic = 8  # motif min IC value
chromvar_iterations = 5000  # ChromVAR iteration
target_tfs = c("SOX10", "SOX9", "OLIG2", "MEF2C", "GFAP", "NeuN", "FOXO1")  # Target TFs
target_cell_types = c("Oligodendrocytes","Inhibitory_Neurons","Microglia","Astrocytes","OPC","Excitatory_Neurons")  # Target cell types
seed = 12345  # Random seed
set.seed(seed)

start_time <- Sys.time()
cat("Time:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")

# ==============================================================================
# 1. Loading and preprocessing
# ==============================================================================
atac_obj <- readRDS(input_file)
cat("   - Loaded ATAC object with", ncol(atac_obj), "cells and", nrow(atac_obj), "peaks\n")
# Check cell type annotation
if (!"cell_type" %in% colnames(atac_obj@meta.data)) {
  stop("NO cell_type in ATAC obj...")
}
cell_type_counts <- table(atac_obj$cell_type)
cat("   - Cell type distribution:\n")
for (ct in names(cell_type_counts)) {
  cat("     *", ct, ":", cell_type_counts[ct], "cells\n")
}
# Set default assay
DefaultAssay(atac_obj) <- "ATAC"
seqlevels(atac_obj)
# Only keep standed Chr, avoid ChromVAR ERROR
gr <- granges(atac_obj)
seq_keep <- seqnames(gr) %in% seqnames(BSgenome.Hsapiens.UCSC.hg38) 
seq_keep <- as.vector(seq_keep)
feat.keep <- GRangesToString(grange = gr[seq_keep])
atac_obj[['ATAC']] <- subset(atac_obj[["ATAC"]], features = feat.keep)

# ==============================================================================
# 2. Retrieving and filtering JASPAR motifs
# ==============================================================================
# JASPAR2024
jaspar <- JASPAR2024::JASPAR2024()
tryCatch({
  pfm_all <- getMatrixSet(
    jaspar@db,
    opts = list(
      collection = "CORE",
#      tax_group = "vertebrates",
      species = "Homo sapiens",
      all_versions = FALSE
    )
  )
  cat("   - Retrieved", length(pfm_all), "motifs from JASPAR2024\n")
}, error = function(e) {
  stop("Retriev JASPAR motifs failed: ", e$message)
})

# Calculate PFM Information Coefficient(IC)
compute_motif_ic <- function(pfm) {
  mat <- Matrix(pfm)  # Data matric
  mat <- mat + 0.8  # Add pseudo to avoid log(0)
  col_sums <- colSums(mat)  # Sum of cols
  prob_mat <- t(t(mat) / col_sums)  # Probability matrix
  # IC for each col
  ic_cols <- apply(prob_mat, 2, function(col_probs) {
    entropy <- -sum(col_probs * log2(col_probs), na.rm = TRUE)  # entropy = -sum(p * log2(p))
    2 - entropy  # IC = 2 - entropy
  })
  sum(ic_cols)
}

cat("Calculating motif information content (IC)...\n")
motif_ic <- sapply(pfm_all, function(motif) {
  tryCatch({
    compute_motif_ic(motif)
  }, error = function(e) {
    warning("IC calculated failed: ", name(motif), " (", ID(motif), "): ", e$message)
    return(0)
  })
})
# show IC summary
ic_summary <- summary(motif_ic)
print(ic_summary)
# filter IC
high_quality_mask <- motif_ic >= min_motif_ic
pfm_filtered <- pfm_all[high_quality_mask]

cat("   - Filtered to", length(pfm_filtered), "high-quality motifs (IC >=", min_motif_ic, ")\n")

motif_ids <- TFBSTools::ID(pfm_all)
motif_names <- TFBSTools::name(pfm_all)
motif_mapping <- data.frame(
  motif_id = motif_ids[high_quality_mask],
  motif_name = motif_names[high_quality_mask],
  ic_score = motif_ic[high_quality_mask],
  stringsAsFactors = FALSE
)
# SAVE/PLOT:
write.csv(motif_mapping, "results/motif/motif_mapping.csv", row.names = FALSE)

# ==============================================================================
# 3. Adding motif information to ATAC object
# ==============================================================================
genome <- BSgenome.Hsapiens.UCSC.hg38
tryCatch({
  # add motifs to Seurat obj
  atac_obj <- AddMotifs(
    object = atac_obj,
    genome = genome,
    pfm = pfm_filtered,
    assay = "ATAC"
  )
  cat("   - Successfully added", length(pfm_filtered), "motifs to ATAC object\n")
  motif_assay <- atac_obj[["ATAC"]]@motifs
  cat("   - Motif matches found in", nrow(motif_assay@data), "peaks\n")
  
}, error = function(e) {
  stop("Add motifs failed: ", e$message)
})

# ==============================================================================
# 4. ChromVAR analysis
# ==============================================================================
tryCatch({
  atac_obj <- RunChromVAR(
    object = atac_obj,
    genome = genome,
    assay = "ATAC",
    new.assay.name = "chromvar",
    verbose = TRUE
  )
  cat("   - ChromVAR analysis completed successfully\n")
  cat("   - Created chromvar assay with", nrow(atac_obj[["chromvar"]]), "motifs\n")
}, error = function(e) {
  stop("ChromVAR analysis failed: ", e$message)
})

# ==============================================================================
# 5. Cell type-specific TF activities
# ==============================================================================
DefaultAssay(atac_obj) <- "chromvar"
Idents(atac_obj) <- atac_obj$cell_type

available_types <- intersect(target_cell_types, unique(atac_obj$cell_type))
atac_subset <- atac_obj

avg_activity <- AverageExpression(
  atac_subset,
  assays = "chromvar",
  group.by = "cell_type"
)$chromvar

# SAVE/PLOT: average tf activity
write.csv(avg_activity, "results/motif/average_tf_activity.csv")

# Extract chromvar matrix (row = feature,col = cell)
data_matrix <- GetAssayData(atac_subset, assay = "chromvar", slot = "data")
feature_variances <- apply(data_matrix, 1, var)  # var for each feature
var_features <- names(sort(feature_variances, decreasing = TRUE))[1:200]  # top 200 var
VariableFeatures(atac_subset, assay = "chromvar") <- var_features  # set to variable feature

# ==============================================================================
# 6. Differential TF activity analysis
# ==============================================================================
# Paired differential
diff_results <- list()

for (i in 1:(length(available_types)-1)) {
  for (j in (i+1):length(available_types)) {
    ct1 <- available_types[i]
    ct2 <- available_types[j]
    cat("   - Comparing", ct1, "vs", ct2, "...\n")
    tryCatch({
      diff_tf <- FindMarkers(
        atac_subset,
        ident.1 = ct1,
        ident.2 = ct2,
        assay = "chromvar",
        test.use = "wilcox",
        min.pct = 0.1,
        logfc.threshold = 0.1
      )
      # Add motif anno
      diff_tf$motif_id <- rownames(diff_tf)
      diff_tf <- merge(diff_tf, motif_mapping, by = "motif_id", all.x = TRUE)
      diff_tf <- diff_tf[order(diff_tf$p_val_adj), ]
      
      comparison_name <- paste0(ct1, "_vs_", ct2)
      diff_results[[comparison_name]] <- diff_tf
      # SAVE/PLOT:
      write.csv(
        diff_tf, 
        paste0("results/motif/differential_TF_", comparison_name, ".csv"),
        row.names = FALSE
      )
      cat("     * Found", sum(diff_tf$p_val_adj < 0.05), "significantly different TFs\n")
      
    }, error = function(e) {
      cat("     * Error in comparison:", e$message, "\n")
    })
  }
}

# ==============================================================================
# 7. TF activity visualizations
# ==============================================================================
# 7.1 Heatmap for top variable TFs
top_vars <- head(VariableFeatures(atac_subset, assay = "chromvar"), 50)
# row correspond to top_vars
sub_mat <- avg_activity[rownames(avg_activity) %in% top_vars, ]

common_vars <- intersect(top_vars, rownames(sub_mat))
sub_mat <- sub_mat[common_vars, ]  # keep common rows
sub_mat <- sub_mat[order(match(rownames(sub_mat), top_vars)), ]  # rank by top_vars
# Convert to normal matrix ( if is dgeMatrix)
heatmap_data <- as.matrix(sub_mat)
# Normalize (Z-score)
heatmap_scaled <- t(scale(t(heatmap_data)))
# Add motif name
motif_names_map <- setNames(motif_mapping$motif_name, motif_mapping$motif_id)
rownames(heatmap_scaled) <- motif_names_map[rownames(heatmap_scaled)]
rownames(heatmap_scaled)[is.na(rownames(heatmap_scaled))] <- rownames(heatmap_data)[is.na(rownames(heatmap_scaled))]
# SAVE/PLOT: 
pdf("plots/motif/tf_activity_heatmap.pdf", width = 10, height = 12)
pheatmap(
  heatmap_scaled,
  name = "Z-score",
  scale = "row",          # row normalize (Z-score)
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  color = colorRampPalette(c("navy", "white", "red"))(100)
)
dev.off()

# 7.2 violin plot for target TFs
target_motif_ids <- c()
for (tf in target_tfs) {
  matching_ids <- motif_mapping$motif_id[grep(tf, motif_mapping$motif_name, ignore.case = TRUE)]
  if (length(matching_ids) > 0) {
    target_motif_ids <- c(target_motif_ids, matching_ids[1])  # 取第一个匹配
  }
}

if (length(target_motif_ids) > 0) {
  # SAVE/PLOT: 
  p_violin <- VlnPlot(
    atac_subset,
    features = target_motif_ids,
    ncol = min(3, length(target_motif_ids)),
    pt.size = 0.1,
    assay = "chromvar"
  )
  
  ggsave("plots/motif/target_tf_violin_plots.pdf", p_violin, width = 15, height = 10)
}

# 7.3 UMAP of selected TFs
if ("umap" %in% names(atac_subset@reductions) && length(target_motif_ids) > 0) {
  p_umap <- FeaturePlot(
    atac_subset,
    features = target_motif_ids,
    reduction = "umap",
    ncol = 3,
    min.cutoff = NA,
    max.cutoff = NA,
    pt.size = 0.3,
    order = TRUE
  )
  # SAVE/PLOT: 
  ggsave("plots/motif/target_tf_umap_plots.pdf", p_umap, width = 18, height = 12)
}

# 7.4 Differential TFs barplot
if (length(diff_results) > 0) {
  for (comparison in names(diff_results)) {
    df <- diff_results[[comparison]]
    # Top and Bottom TFs
    top_up <- head(df[df$avg_log2FC > 0 & df$p_val_adj < 0.05, ], 10)
    top_down <- head(df[df$avg_log2FC < 0 & df$p_val_adj < 0.05, ], 10)
    
    if (nrow(top_up) > 0 || nrow(top_down) > 0) {
      plot_data <- rbind(top_up, top_down)
      plot_data$direction <- ifelse(plot_data$avg_log2FC > 0, "Up", "Down")
      plot_data$neg_log10_p <- -log10(plot_data$p_val_adj)
      
      p_bar <- ggplot(plot_data, aes(x = reorder(motif_name, avg_log2FC), y = avg_log2FC)) +
        geom_col(aes(fill = direction), alpha = 0.8) +
        scale_fill_manual(values = c("Up" = "red", "Down" = "blue")) +
        coord_flip() +
        labs(
          title = paste("Top Differential TFs:", comparison),
          x = "Transcription Factor",
          y = "Average Log2 Fold Change",
          fill = "Direction"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14),
          axis.text.y = element_text(size = 10)
        )
      # SAVE/PLOT: 
      ggsave(
        paste0("plots/motif/differential_tf_barplot_", comparison, ".pdf"),
        p_bar,
        width = 10,
        height = 8
      )
    }
  }
}

# ==============================================================================
# 8. Detailed analysis of target transcription factors
# ==============================================================================
target_tf_results <- list()

for (tf_name in target_tfs) {
  cat("   - Analyzing", tf_name, "...\n")
  # match motifs
  matching_motifs <- motif_mapping[grep(tf_name, motif_mapping$motif_name, ignore.case = TRUE), ]
  
  if (nrow(matching_motifs) == 0) {
    cat("     * No motifs found for", tf_name, "\n")
    next
  }
  # best match
  best_motif <- matching_motifs[which.max(matching_motifs$ic_score), ]
  motif_id <- best_motif$motif_id
  
  cat("     * Using motif:", best_motif$motif_name, "(", motif_id, ")\n")
  
  # TF activity in each cell type
  tf_activity <- FetchData(
    atac_subset,
    vars = c(motif_id, "cell_type"),
    slot = "data"
  )
  
  tf_stats <- tf_activity %>%
    group_by(cell_type) %>%
    summarise(
      mean_activity = mean(get(motif_id)),
      median_activity = median(get(motif_id)),
      sd_activity = sd(get(motif_id)),
      .groups = "drop"
    )
  
  target_tf_results[[tf_name]] <- list(
    motif_info = best_motif,
    activity_stats = tf_stats,
    raw_activity = tf_activity
  )
  # SAVE/PLOT: TF activity boxplot
  p_box <- ggplot(tf_activity, aes(x = cell_type, y = get(motif_id), fill = cell_type)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.1) +
    labs(
      title = paste("TF Activity:", tf_name),
      subtitle = paste("Motif:", best_motif$motif_name),
      x = "Cell Type",
      y = "ChromVAR Activity Score",
      fill = "Cell Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    scale_fill_viridis_d()
  
  ggsave(
    paste0("plots/motif/tf_activity_boxplot_", tf_name, ".pdf"),
    p_box,
    width = 8,
    height = 6
  )
}
# SAVE/PLOT: TF analysis results
saveRDS(target_tf_results, "results/motif/target_tf_analysis.rds")

# ==============================================================================
# 9. TF co-activity network
# ==============================================================================
# TF activity
tf_activity_matrix <- FetchData(
  atac_subset,
  vars = target_motif_ids,
  slot = "data"
)
# correlation matrix
cor_matrix <- cor(tf_activity_matrix, use = "complete.obs")
tf_names_for_ids <- sapply(target_motif_ids, function(id) {
  name <- motif_mapping$motif_name[motif_mapping$motif_id == id]
  if (length(name) == 0) return(id)
  return(name[1])
})
rownames(cor_matrix) <- tf_names_for_ids
colnames(cor_matrix) <- tf_names_for_ids

# SAVE/PLOT:
pdf("plots/motif/tf_correlation_heatmap.pdf", width = 10, height = 10)
Heatmap(
  cor_matrix,
  name = "Correlation",
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", cor_matrix[i, j]), x, y, gp = gpar(fontsize = 8))
  }
)
dev.off()
# SAVE/PLOT: TF correlation matrix
write.csv(cor_matrix, "results/motif/tf_correlation_matrix.csv")

# ==============================================================================
# 10. Quality control and validation
# ==============================================================================
# 10.1 ChromVAR quality 
chromvar_qc <- data.frame(
  metric = c(
    "Total_motifs",
    "Variable_motifs",
    "Mean_activity_range",
    "Cells_analyzed"
  ),
  value = c(
    nrow(atac_subset[["chromvar"]]),
    length(VariableFeatures(atac_subset, assay = "chromvar")),
    round(mean(apply(avg_activity, 1, function(x) max(x) - min(x))), 3),
    ncol(atac_subset)
  )
)
# SAVE/PLOT: ChromVAR quality
write.csv(chromvar_qc, "results/motif/chromvar_qc_metrics.csv", row.names = FALSE)

# 10.2 Check TFs markers
known_markers <- list(
  "Oligodendrocytes" = c("SOX10", "OLIG2", "MYF"),
  "Astrocytes" = c("SOX9", "NFIA", "NFIB"),
  "Neurons" = c("NEUROD", "MEF2", "FOXO")
)

marker_validation <- data.frame()

for (cell_type in names(known_markers)) {
  if (cell_type %in% available_types) {
    for (marker in known_markers[[cell_type]]) {
      matching_motifs <- motif_mapping[grep(marker, motif_mapping$motif_name, ignore.case = TRUE), ]
      
      if (nrow(matching_motifs) > 0) {
        motif_id <- matching_motifs$motif_id[1]
        
        # Related activity of matched TF
        cell_type_activity <- mean(FetchData(
          subset(atac_subset, idents = cell_type),
          vars = motif_id,
          slot = "data"
        )[, 1])
        
        other_activity <- mean(FetchData(
          subset(atac_subset, idents = available_types[available_types != cell_type]),
          vars = motif_id,
          slot = "data"
        )[, 1])
        
        fold_change <- cell_type_activity / (other_activity + 0.001)
        
        marker_validation <- rbind(marker_validation, data.frame(
          cell_type = cell_type,
          tf_marker = marker,
          motif_id = motif_id,
          motif_name = matching_motifs$motif_name[1],
          target_activity = round(cell_type_activity, 3),
          other_activity = round(other_activity, 3),
          fold_change = round(fold_change, 3)
        ))
      }
    }
  }
}

if (nrow(marker_validation) > 0) {
  # SAVE/PLOT: TF marker validation
  write.csv(marker_validation, "results/motif/marker_tf_validation.csv", row.names = FALSE)
  
  p_validation <- ggplot(marker_validation, aes(x = tf_marker, y = fold_change, fill = cell_type)) +
    geom_col(position = "dodge", alpha = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    facet_wrap(~cell_type, scales = "free_x") +
    labs(
      title = "Known TF Marker Validation",
      x = "TF Marker",
      y = "Fold Change (Target vs Others)",
      fill = "Cell Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 12)
    )
  # SAVE/PLOT: TF marker
  ggsave("plots/motif/tf_marker_validation.pdf", p_validation, width = 12, height = 8)
}

# ==============================================================================
# 11. Saving final results
# ==============================================================================
# SAVE/PLOT: chromvar ATAC obj
saveRDS(atac_obj, "results/motif/atac_with_chromvar.rds")

# Motif summary
motif_summary <- list(
  parameters = params,
  motif_mapping = motif_mapping,
  differential_results = diff_results,
  target_tf_results = target_tf_results,
  average_activity = avg_activity,
  qc_metrics = chromvar_qc,
  marker_validation = if(exists("marker_validation")) marker_validation else NULL
)
# SAVE/PLOT: motif summary
saveRDS(motif_summary, "results/motif/motif_analysis_summary.rds")

# ==============================================================================
# 12. Analysis report
# ==============================================================================
report <- paste0(
  "=== Motif and ChromVAR Analysis Report ===\n",
  "Analysis Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
  "Input File: ", input_file, "\n\n",
  "Data Summary:\n",
  "- Total cells analyzed: ", ncol(atac_subset), "\n",
  "- Cell types: ", paste(available_types, collapse = ", "), "\n",
  "- Total motifs from JASPAR: ", length(pfm_all), "\n",
  "- High-quality motifs used: ", length(pfm_filtered), "\n\n",
  "ChromVAR Results:\n",
  "- Motif activity scores calculated for ", nrow(atac_obj[["chromvar"]]), " motifs\n",
  "- Variable motifs identified: ", length(VariableFeatures(atac_subset, assay = "chromvar")), "\n\n",
  "Differential Analysis:\n"
)

# Differential comparisons results
if (length(diff_results) > 0) {
  for (comparison in names(diff_results)) {
    n_sig <- sum(diff_results[[comparison]]$p_val_adj < 0.05)
    report <- paste0(report, "- ", comparison, ": ", n_sig, " significantly different TFs\n")
  }
} else {
  report <- paste0(report, "- No differential comparisons performed\n")
}

report <- paste0(report, "\nTarget TF Analysis:\n")
for (tf in names(target_tf_results)) {
  motif_name <- target_tf_results[[tf]]$motif_info$motif_name
  report <- paste0(report, "- ", tf, " (", motif_name, "): analyzed\n")
}

report <- paste0(report, "\nOutput Files:\n",
                 "- Motif mapping: results/motif/motif_mapping.csv\n",
                 "- Average TF activity: results/motif/average_tf_activity.csv\n",
                 "- Differential results: results/motif/differential_TF_*.csv\n",
                 "- Target TF analysis: results/motif/target_tf_analysis.rds\n",
                 "- ChromVAR object: data/motif/atac_with_chromvar.rds\n",
                 "- Visualization plots: plots/motif/\n"
)
# SAVE/PLOT: analysis report
writeLines(report, "results/motif/analysis_report.txt")

# Run time
end_time <- Sys.time()
runtime <- round(difftime(end_time, start_time, units = "mins"), 2)

cat("\n=== Motif and ChromVAR Analysis Completed ===\n")
cat("Total runtime:", runtime, "minutes\n")
cat("Results saved in: results/motif/\n")
cat("Plots saved in: plots/motif/\n")

# Key Results
cat("\nKey Results Summary:\n")
cat("- Total motifs analyzed:", nrow(atac_obj[["chromvar"]]), "\n")
cat("- Variable motifs:", length(VariableFeatures(atac_subset, assay = "chromvar")), "\n")
if (length(diff_results) > 0) {
  total_sig <- sum(sapply(diff_results, function(x) sum(x$p_val_adj < 0.05)))
  cat("- Total significant differential TFs:", total_sig, "\n")
}
cat("- Target TFs analyzed:", length(target_tf_results), "\n")

# Top variable TFs
top_variable <- head(VariableFeatures(atac_subset, assay = "chromvar"), 5)
top_variable_names <- sapply(top_variable, function(x) {
  name <- motif_mapping$motif_name[motif_mapping$motif_id == x]
  if (length(name) > 0) return(name[1]) else return(x)
})
cat("- Top variable TFs:", paste(top_variable_names, collapse = ", "), "\n")
