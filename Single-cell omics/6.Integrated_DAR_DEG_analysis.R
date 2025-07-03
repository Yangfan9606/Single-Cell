#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-07-03
# version:   1.2
# license:   MIT
# brief:     Differentially accessible region (DAR) and differential gene expression (DGE) analysis.
#------------------------------------------#

# Load required libraries
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(VennDiagram)
library(RColorBrewer)
library(viridis)
library(patchwork)
library(gridExtra)
library(cowplot)
# Load Peak annotation libraries
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(annotatr)

# Output dir
output_dirs <- c("results/dar_deg", "plots/dar_deg")
sapply(output_dirs, function(x) if(!dir.exists(x)) dir.create(x, recursive = TRUE))

# Parameters
joint_file = "CTL.joint_integrated.rds"

# comparison group - based on normal neuronal cell types
comparisons = list(
  c("Astrocytes", "Oligodendrocytes"), c("Astrocytes", "Excitatory_Neurons"), c("Oligodendrocytes", "Excitatory_Neurons"), 
  c("Inhibitory_Neurons", "Astrocytes"), c("Inhibitory_Neurons", "Excitatory_Neurons"), c("Microglia", "Astrocytes")
)

# DAR
dar_params = list(
  min_pct = 0.05,
  logfc_threshold = 0.25,
  test_use = "LR",
  latent_vars = c("nCount_ATAC", "nFeature_ATAC")
)
# DEG
deg_params = list(
  min_pct = 0.1,
  logfc_threshold = 0.25,
  test_use = "wilcox"
)

fdr_threshold = 0.05; min_cells_per_group = 20; seed = 12345
set.seed(seed) # Random seed

start_time <- Sys.time()
cat("Time:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")

# ==============================================================================
# 1. Loading and validating data
# ==============================================================================
joint_obj <- readRDS(joint_file)
cat("   - Loaded joint object:", ncol(joint_obj), "cells\n")
# RNA and ATAC data in joint_obj
rna_cells <- joint_obj$dataset == "snRNA"
atac_cells <- joint_obj$dataset == "snATAC"
cat("   - RNA cells:", sum(rna_cells), "\n")
cat("   - ATAC cells:", sum(atac_cells), "\n")
# Extract RNA and ATAC data
rna_obj <- subset(joint_obj, cells = colnames(joint_obj)[rna_cells])
DefaultAssay(rna_obj) <- "RNA"

atac_obj <- subset(joint_obj, cells = colnames(joint_obj)[atac_cells])
DefaultAssay(atac_obj) <- "ATAC"

# Check cell type annotation
check_cell_types <- function(obj, obj_name) {
  if (is.null(obj)) return(NULL)
  possible_names <- c("cell_type", "predicted.id", "seurat_clusters", "celltype")
  cell_type_col <- NULL
  for (name in possible_names) {
    if (name %in% colnames(obj@meta.data)) {
      cell_type_col <- name
      break
    }
  }
  if (is.null(cell_type_col)) {
    cat("   - Warning:", obj_name, "object lacks cell type annotation\n")
    cat("   - Available columns:", paste(colnames(obj@meta.data), collapse = ", "), "\n")
    return(NULL)
  }
  # Standard cell type col names
  if (cell_type_col != "cell_type") {
    obj$cell_type <- obj@meta.data[[cell_type_col]]
  }
  cell_types <- unique(obj$cell_type)
  cell_counts <- table(obj$cell_type)
  cat("   -", obj_name, "cell types:\n")
  for (ct in names(cell_counts)) {
    cat("     *", ct, ":", cell_counts[ct], "cells\n")
  }
  return(cell_types)
}

atac_cell_types <- check_cell_types(atac_obj, "ATAC")
rna_cell_types <- check_cell_types(rna_obj, "RNA")

# Check comparison groups
valid_comparisons <- list()
for (comp in comparisons) {
  if (all(comp %in% atac_cell_types)) {
    n1 <- sum(atac_obj$cell_type == comp[1])
    n2 <- sum(atac_obj$cell_type == comp[2])
    if (n1 >= min_cells_per_group && n2 >= min_cells_per_group) {
      valid_comparisons <- append(valid_comparisons, list(comp))
    } else {
      cat("   - Warning: Comparison", paste(comp, collapse = " vs "), 
          "skipped - insufficient cells (", n1, "vs", n2, ")\n")
    }
  } else {
    cat("   - Warning: Comparison", paste(comp, collapse = " vs "), 
        "skipped - cell types not available\n")
  }
}
if (length(valid_comparisons) == 0) {  stop("NO Comparison Avaliable !")  }
cat("   - Valid comparisons:", length(valid_comparisons), "\n")

# ==============================================================================
# 2. Differential Accessibility Region (DAR) analysis
# ==============================================================================
Idents(atac_obj) <- atac_obj$cell_type

dar_results <- list()
dar_summary <- data.frame()

for (i in seq_along(valid_comparisons)) {
  comp <- valid_comparisons[[i]]
  comp_name <- paste(comp[1], "vs", comp[2])
  cat("   - Analyzing", comp_name, "...\n")
  n1 <- sum(atac_obj$cell_type == comp[1])
  n2 <- sum(atac_obj$cell_type == comp[2])
  cat("     *", comp[1], ":", n1, "cells\n")
  cat("     *", comp[2], ":", n2, "cells\n")
  tryCatch({
    # DAR analysis
    dar <- FindMarkers(
      object = atac_obj,
      ident.1 = comp[1],
      ident.2 = comp[2],
      min.pct = dar_params$min_pct,
      logfc.threshold = dar_params$logfc_threshold,
      test.use = dar_params$test_use,
      latent.vars = dar_params$latent_vars,
      verbose = FALSE
    )
    # Add peak info and anno
    dar$peak <- rownames(dar)
    dar$comparison <- comp_name
    dar$group1 <- comp[1]
    dar$group2 <- comp[2]
    dar$direction <- ifelse(dar$avg_log2FC > 0, 
                            paste("Higher in", comp[1]), 
                            paste("Higher in", comp[2]))
    dar_sig <- dar[dar$p_val_adj < fdr_threshold, ]  # Padj filter
    n_total <- nrow(dar)
    n_sig <- nrow(dar_sig)
    n_up_g1 <- sum(dar_sig$avg_log2FC > 0)
    n_up_g2 <- sum(dar_sig$avg_log2FC < 0)
    cat("     * Total peaks tested:", n_total, "\n")
    cat("     * Significant DARs:", n_sig, "(", round(n_sig/n_total*100, 1), "%)\n")
    cat("     * Higher in", comp[1], ":", n_up_g1, "\n")
    cat("     * Higher in", comp[2], ":", n_up_g2, "\n")
    # summary
    dar_summary <- rbind(dar_summary, data.frame(
      comparison = comp_name,
      group1 = comp[1],
      group2 = comp[2],
      n_cells_g1 = n1,
      n_cells_g2 = n2,
      n_peaks_tested = n_total,
      n_dar_significant = n_sig,
      n_higher_g1 = n_up_g1,
      n_higher_g2 = n_up_g2,
      pct_significant = round(n_sig/n_total*100, 2)
    ))
    
    dar_results[[comp_name]] <- list(
      all_results = dar,
      significant = dar_sig,
      comparison = comp,
      stats = list(
        n_total = n_total,
        n_sig = n_sig,
        n_up_g1 = n_up_g1,
        n_up_g2 = n_up_g2
      )
    )

    comp_name_clean <- gsub(" ", "_", comp_name)
    write.csv(dar, paste0("results/dar_deg/DAR_all_", comp_name_clean, ".csv"), row.names = TRUE)
    write.csv(dar_sig, paste0("results/dar_deg/DAR_significant_", comp_name_clean, ".csv"), row.names = TRUE)
    
  }, error = function(e) {
    cat("     * Error in DAR analysis:", e$message, "\n")
  })
}

# SAVE:
write.csv(dar_summary, "results/dar_deg/DAR_summary.csv", row.names = FALSE)

# ==============================================================================
# 3. Differential Expression Gene (DEG) analysis
# ==============================================================================

deg_results <- list()
deg_summary <- data.frame()

if (!is.null(rna_obj)) {
  Idents(rna_obj) <- rna_obj$cell_type
  for (i in seq_along(valid_comparisons)) {
    comp <- valid_comparisons[[i]]
    comp_name <- paste(comp[1], "vs", comp[2])
    # Check comparsion group
    if (!all(comp %in% unique(rna_obj$cell_type))) {
      cat("   - Skipping", comp_name, "- not available in RNA data\n")
      next
    }
    cat("   - Analyzing", comp_name, "...\n")
    n1 <- sum(rna_obj$cell_type == comp[1])
    n2 <- sum(rna_obj$cell_type == comp[2])
    cat("     *", comp[1], ":", n1, "cells\n")
    cat("     *", comp[2], ":", n2, "cells\n")
    if (n1 < min_cells_per_group || n2 < min_cells_per_group) {
      cat("     * Warning: Too few cells for reliable analysis\n")
      next
    }
    tryCatch({
      # DEG analysis
      deg <- FindMarkers(
        object = rna_obj,
        ident.1 = comp[1],
        ident.2 = comp[2],
        min.pct = deg_params$min_pct,
        logfc.threshold = deg_params$logfc_threshold,
        test.use = deg_params$test_use,
        verbose = FALSE
      )
      
      deg$gene <- rownames(deg)
      deg$comparison <- comp_name
      deg$group1 <- comp[1]
      deg$group2 <- comp[2]
      deg$direction <- ifelse(deg$avg_log2FC > 0, 
                              paste("Higher in", comp[1]), 
                              paste("Higher in", comp[2]))
      deg_sig <- deg[deg$p_val_adj < fdr_threshold, ]  #  Padj filter
      n_total <- nrow(deg)
      n_sig <- nrow(deg_sig)
      n_up_g1 <- sum(deg_sig$avg_log2FC > 0)
      n_up_g2 <- sum(deg_sig$avg_log2FC < 0)
      cat("     * Total genes tested:", n_total, "\n")
      cat("     * Significant DEGs:", n_sig, "(", round(n_sig/n_total*100, 1), "%)\n")
      cat("     * Higher in", comp[1], ":", n_up_g1, "\n")
      cat("     * Higher in", comp[2], ":", n_up_g2, "\n")
      # summary
      deg_summary <- rbind(deg_summary, data.frame(
        comparison = comp_name,
        group1 = comp[1],
        group2 = comp[2],
        n_cells_g1 = n1,
        n_cells_g2 = n2,
        n_genes_tested = n_total,
        n_deg_significant = n_sig,
        n_higher_g1 = n_up_g1,
        n_higher_g2 = n_up_g2,
        pct_significant = round(n_sig/n_total*100, 2)
      ))
      
      deg_results[[comp_name]] <- list(
        all_results = deg,
        significant = deg_sig,
        comparison = comp,
        stats = list(
          n_total = n_total,
          n_sig = n_sig,
          n_up_g1 = n_up_g1,
          n_up_g2 = n_up_g2
        )
      )
      comp_name_clean <- gsub(" ", "_", comp_name)
      write.csv(deg, paste0("results/dar_deg/DEG_all_", comp_name_clean, ".csv"), row.names = TRUE)
      write.csv(deg_sig, paste0("results/dar_deg/DEG_significant_", comp_name_clean, ".csv"), row.names = TRUE)
      
    }, error = function(e) {
      cat("     * Error in DEG analysis:", e$message, "\n")
    })
  }
  
  # SAVE:
  if (nrow(deg_summary) > 0) {
    write.csv(deg_summary, "results/dar_deg/DEG_summary.csv", row.names = FALSE)
  }
} else {  cat("\n3. DEG analysis skipped - no RNA data available\n")  }

# ==============================================================================
# 4. DAR visualizations
# ==============================================================================
create_dar_plots <- function(comp_name, dar_data, dar_sig) {
  # DAR data
  dar_data$significant <- ifelse(dar_data$p_val_adj < fdr_threshold, 
                                 "Significant", "Not Significant")
  dar_data$neg_log10_padj <- -log10(dar_data$p_val_adj)
  dar_data$neg_log10_padj[is.infinite(dar_data$neg_log10_padj)] <- 
    max(dar_data$neg_log10_padj[!is.infinite(dar_data$neg_log10_padj)]) + 10
  plots <- list()
  # 4.1 volcano
  plots$volcano <- ggplot(dar_data, aes(x = avg_log2FC, y = neg_log10_padj)) +
    geom_point(aes(color = significant), alpha = 0.6, size = 0.8) +
    scale_color_manual(values = c("Not Significant" = "grey70", "Significant" = "#E31A1C")) +
    geom_vline(xintercept = c(-dar_params$logfc_threshold, dar_params$logfc_threshold), 
               linetype = "dashed", color = "#1F78B4", alpha = 0.7) +
    geom_hline(yintercept = -log10(fdr_threshold), 
               linetype = "dashed", color = "#1F78B4", alpha = 0.7) +
    labs(
      title = paste("DAR Analysis:", comp_name),
      subtitle = paste("Significant DARs:", nrow(dar_sig)),
      x = "Average Log2 Fold Change",
      y = "-Log10 Adjusted P-value",
      color = "Significance"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "bottom"
    )
  # 4.2 MA plot
  dar_data$avg_expr <- (dar_data$pct.1 + dar_data$pct.2) / 2
  plots$ma <- ggplot(dar_data, aes(x = avg_expr, y = avg_log2FC)) +
    geom_point(aes(color = significant), alpha = 0.6, size = 0.8) +
    scale_color_manual(values = c("Not Significant" = "grey70", "Significant" = "#E31A1C")) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", alpha = 0.5) +
    geom_hline(yintercept = c(-dar_params$logfc_threshold, dar_params$logfc_threshold), 
               linetype = "dashed", color = "#1F78B4", alpha = 0.7) +
    labs(
      title = paste("DAR MA Plot:", comp_name),
      x = "Average Expression (%)",
      y = "Average Log2 Fold Change",
      color = "Significance"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom"
    )
  # 4.3 Sig DAR distribution
  if (nrow(dar_sig) > 0) {
    plots$distribution <- ggplot(dar_sig, aes(x = avg_log2FC, fill = direction)) +
      geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
      scale_fill_brewer(type = "qual", palette = "Set1") +
      labs(
        title = paste("Distribution of Significant DARs:", comp_name),
        x = "Average Log2 Fold Change",
        y = "Count",
        fill = "Direction"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom"
      )
  }
  return(plots)
}

# Create DAR plot for each comparsion group
dar_plot_list <- list()
for (comp_name in names(dar_results)) {
  cat("   - Creating DAR plots for", comp_name, "...\n")
  dar_data <- dar_results[[comp_name]]$all_results
  dar_sig <- dar_results[[comp_name]]$significant
  plots <- create_dar_plots(comp_name, dar_data, dar_sig)
  dar_plot_list[[comp_name]] <- plots
  # SAVE：
  comp_name_clean <- gsub(" ", "_", comp_name)
  ggsave(paste0("plots/dar_deg/DAR_volcano_", comp_name_clean, ".pdf"), 
         plots$volcano, width = 10, height = 8)
  ggsave(paste0("plots/dar_deg/DAR_MA_", comp_name_clean, ".pdf"), 
         plots$ma, width = 10, height = 8)
  if ("distribution" %in% names(plots)) {
    ggsave(paste0("plots/dar_deg/DAR_distribution_", comp_name_clean, ".pdf"), 
           plots$distribution, width = 10, height = 6)
  }
}

# DAR summary plot
if (nrow(dar_summary) > 0) {
  # Sig DAR number Bar plot
  p_dar_counts <- ggplot(dar_summary, aes(x = reorder(comparison, n_dar_significant))) +
    geom_col(aes(y = n_dar_significant), fill = "#E31A1C", alpha = 0.7) +
    geom_text(aes(y = n_dar_significant, label = n_dar_significant), 
              hjust = -0.1, size = 3) +
    coord_flip() +
    labs(
      title = "Number of Significant DARs by Comparison",
      x = "Comparison",
      y = "Number of Significant DARs"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"))
  ggsave("plots/dar_deg/DAR_summary_counts.pdf", p_dar_counts, width = 10, height = 8)
  
  # Sig percent
  p_dar_pct <- ggplot(dar_summary, aes(x = reorder(comparison, pct_significant))) +
    geom_col(aes(y = pct_significant), fill = "#1F78B4", alpha = 0.7) +
    geom_text(aes(y = pct_significant, label = paste0(pct_significant, "%")), 
              hjust = -0.1, size = 3) +
    coord_flip() +
    labs(
      title = "Percentage of Significant DARs by Comparison",
      x = "Comparison",
      y = "Percentage of Significant DARs (%)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"))
  ggsave("plots/dar_deg/DAR_summary_percentage.pdf", p_dar_pct, width = 10, height = 8)
}

# ==============================================================================
# 5. DEG visualizations
# ==============================================================================
if (length(deg_results) > 0) {
  create_deg_plots <- function(comp_name, deg_data, deg_sig) {
    # DEG data
    deg_data$significant <- ifelse(deg_data$p_val_adj < fdr_threshold, "Significant", "Not Significant")
    max_neg_log10 <- max(
      -log10(deg_data$p_val_adj[!is.infinite(-log10(deg_data$p_val_adj))]),
      na.rm = TRUE) + 10
    deg_data$neg_log10_padj <- -log10(deg_data$p_val_adj)
    deg_data$neg_log10_padj[is.infinite(deg_data$neg_log10_padj)] <- max_neg_log10
    
    deg_sig$neg_log10_padj <- -log10(deg_sig$p_val_adj)
    deg_sig$neg_log10_padj[is.infinite(deg_sig$neg_log10_padj)] <- max_neg_log10
    
    plots <- list()
    # 5.1 Volcano plot
    plots$volcano <- ggplot(deg_data, aes(x = avg_log2FC, y = neg_log10_padj)) +
      geom_point(aes(color = significant), alpha = 0.6, size = 0.8) +
      scale_color_manual(values = c("Not Significant" = "grey70", "Significant" = "#E31A1C")) +
      geom_vline(xintercept = c(-deg_params$logfc_threshold, deg_params$logfc_threshold), 
                 linetype = "dashed", color = "#1F78B4", alpha = 0.7) +
      geom_hline(yintercept = -log10(fdr_threshold), 
                 linetype = "dashed", color = "#1F78B4", alpha = 0.7) +
      labs(
        title = paste("DEG Analysis:", comp_name),
        subtitle = paste("Significant DEGs:", nrow(deg_sig)),
        x = "Average Log2 Fold Change",
        y = "-Log10 Adjusted P-value",
        color = "Significance"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        legend.position = "bottom"
      )
    
    # gene signature
    if (nrow(deg_sig) > 0) {
      top_genes <- deg_sig[order(deg_sig$p_val_adj)[1:min(10, nrow(deg_sig))], ]
      plots$volcano <- plots$volcano +
        ggrepel::geom_text_repel(
          data = top_genes,
          aes(label = gene),
          size = 3,
          max.overlaps = 10
        )
    }
    return(plots)
  }
  
  # DEG for each comparsion group
  for (comp_name in names(deg_results)) {
    cat("   - Creating DEG plots for", comp_name, "...\n")
    deg_data <- deg_results[[comp_name]]$all_results
    deg_sig <- deg_results[[comp_name]]$significant
    plots <- create_deg_plots(comp_name, deg_data, deg_sig)
    comp_name_clean <- gsub(" ", "_", comp_name)
    ggsave(paste0("plots/dar_deg/DEG_volcano_", comp_name_clean, ".pdf"), 
           plots$volcano, width = 12, height = 10)
  }

  # DEG summary
  if (nrow(deg_summary) > 0) {
    # DEG sig number bar plot
    p_deg_counts <- ggplot(deg_summary, aes(x = reorder(comparison, n_deg_significant))) +
      geom_col(aes(y = n_deg_significant), fill = "#E31A1C", alpha = 0.7) +
      geom_text(aes(y = n_deg_significant, label = n_deg_significant), 
                hjust = -0.1, size = 3) +
      coord_flip() +
      labs(
        title = "Number of Significant DEGs by Comparison",
        x = "Comparison",
        y = "Number of Significant DEGs"
      ) +
      theme_minimal() +
      theme(plot.title = element_text(size = 14, face = "bold"))
    ggsave("plots/dar_deg/DEG_summary_counts.pdf", p_deg_counts, width = 10, height = 8)
  }
}

# ==============================================================================
# 6. Functional enrichment
# ==============================================================================
if (length(deg_results) > 0) {
  enrichment_results <- list()
  perform_enrichment <- function(genes, comp_name, direction) {
    if (length(genes) < 5) {
      cat("     * Too few genes for enrichment:", length(genes), "\n")
      return(NULL)
    }
    tryCatch({
      # GO
      ego <- enrichGO(gene = genes,
                      universe = rownames(rna_obj),
                      OrgDb = org.Hs.eg.db,
                      keyType = 'SYMBOL',
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)
      # KEGG
      gene_entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", 
                          OrgDb = org.Hs.eg.db)$ENTREZID
      if (length(gene_entrez) >= 5) {
        kegg <- enrichKEGG(gene = gene_entrez,
                           organism = 'hsa',
                           keyType = 'kegg',
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH",
                           qvalueCutoff = 0.2)
      } else {
        kegg <- NULL
      }
      return(list(GO = ego, KEGG = kegg))
    }, error = function(e) {
      cat("     * Error in enrichment analysis:", e$message, "\n")
      return(NULL)
    })
  }
  # For each comparsion group
  for (comp_name in names(deg_results)) {
    cat("   - Enrichment analysis for", comp_name, "...\n")
    deg_sig <- deg_results[[comp_name]]$significant
    comp <- deg_results[[comp_name]]$comparison
    if (nrow(deg_sig) == 0) {
      cat("     * No significant DEGs for enrichment\n")
      next
    }
    up_genes <- deg_sig$gene[deg_sig$avg_log2FC > 0]
    down_genes <- deg_sig$gene[deg_sig$avg_log2FC < 0]
    enrichment_results[[comp_name]] <- list()
    if (length(up_genes) > 0) {
      cat("     * Analyzing", length(up_genes), "upregulated genes in", comp[1], "...\n")
      enrich_up <- perform_enrichment(up_genes, comp_name, paste("Up in", comp[1]))
      enrichment_results[[comp_name]][["up"]] <- enrich_up
    }
    if (length(down_genes) > 0) {
      cat("     * Analyzing", length(down_genes), "upregulated genes in", comp[2], "...\n")
      enrich_down <- perform_enrichment(down_genes, comp_name, paste("Up in", comp[2]))
      enrichment_results[[comp_name]][["down"]] <- enrich_down
    }
  }
  # SAVE:
  for (comp_name in names(enrichment_results)) {
    comp_name_clean <- gsub(" ", "_", comp_name)
    for (direction in names(enrichment_results[[comp_name]])) {
      enrich_data <- enrichment_results[[comp_name]][[direction]]
      if (is.null(enrich_data)) next
      # GO results
      if (!is.null(enrich_data$GO) && nrow(enrich_data$GO@result) > 0) {
        write.csv(enrich_data$GO@result, 
                  paste0("results/dar_deg/GO_enrichment_", comp_name_clean, "_", direction, ".csv"))
        # GO plot
        if (nrow(enrich_data$GO@result) >= 5) {
          p_go <- dotplot(enrich_data$GO, showCategory = 20) +
            ggtitle(paste("GO Enrichment:", comp_name, "-", direction))
          ggsave(paste0("plots/dar_deg/GO_enrichment_", comp_name_clean, "_", direction, ".pdf"),
                 p_go, width = 12, height = 10)
        }
      }
      # KEGG results
      if (!is.null(enrich_data$KEGG) && nrow(enrich_data$KEGG@result) > 0) {
        write.csv(enrich_data$KEGG@result, 
                  paste0("results/dar_deg/KEGG_enrichment_", comp_name_clean, "_", direction, ".csv"))
        # KEGG plot
        if (nrow(enrich_data$KEGG@result) >= 5) {
          p_kegg <- dotplot(enrich_data$KEGG, showCategory = 15) +
            ggtitle(paste("KEGG Enrichment:", comp_name, "-", direction))
          ggsave(paste0("plots/dar_deg/KEGG_enrichment_", comp_name_clean, "_", direction, ".pdf"),
                 p_kegg, width = 12, height = 8)
        }
      }
    }
  }
}

# ==============================================================================
# 7. DAR-DEG association analysis
# ==============================================================================
# Function to parse peak coordinates
parse_peak_coordinates <- function(peak_names) {
  # Convert peak names like "chr8-27611120-27612035" to GRanges object
  peak_parts <- strsplit(peak_names, "-")
  chr <- sapply(peak_parts, function(x) x[1])
  start <- as.numeric(sapply(peak_parts, function(x) x[2]))
  end <- as.numeric(sapply(peak_parts, function(x) x[3]))
  # Create GRanges object
  peak_gr <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end), peak_id = peak_names)
  return(peak_gr)
}
# Function to annotate peaks with nearest genes
annotate_peaks <- function(peak_gr, txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                           distance_threshold = 50000) {  # Peak to TSS distance, default 50k
  # Annotate peaks using ChIPseeker
  peak_anno <- annotatePeak(peak_gr, 
                            TxDb = txdb,
                            annoDb = "org.Hs.eg.db",
                            assignGenomicAnnotation = TRUE,
                            genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", 
                                                          "Exon", "Intron", "Downstream", 
                                                          "Intergenic"))
  # Convert to data frame
  peak_anno_df <- as.data.frame(peak_anno)
  # Add distance filtering (optional - keep peaks within distance_threshold)
  peak_anno_df$within_threshold <- abs(peak_anno_df$distanceToTSS) <= distance_threshold
  return(peak_anno_df)
}

# Main function for DAR-DEG association analysis
perform_dar_deg_association <- function(dar_results, deg_results, 
                                        txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                                        distance_threshold = 50000,
                                        output_dir = "results/dar_deg/",
                                        plot_dir = "plots/dar_deg/") {
  # Create output directories
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  if (length(dar_results) == 0 || length(deg_results) == 0) {
    cat("No DAR or DEG results provided\n")
    return(NULL)
  }
  association_results <- list()
  for (comp_name in names(dar_results)) {
    if (!comp_name %in% names(deg_results)) {
      cat("   - Skipping", comp_name, "- no matching DEG data\n")
      next
    }
    cat("   - Analyzing associations for", comp_name, "...\n")
    
    dar_sig <- dar_results[[comp_name]]$significant
    deg_sig <- deg_results[[comp_name]]$significant
    
    if (nrow(dar_sig) == 0 || nrow(deg_sig) == 0) {
      cat("     * Insufficient significant results for association\n")
      next
    }
    tryCatch({
      # Step 1: Parse peak coordinates
      cat("     * Parsing peak coordinates...\n")
      peak_gr <- parse_peak_coordinates(dar_sig$peak)
      # Step 2: Annotate peaks
      cat("     * Annotating peaks with nearest genes...\n")
      peak_anno_df <- annotate_peaks(peak_gr, txdb, distance_threshold)
      # Step 3: Merge with DAR results
      dar_annotated <- merge(dar_sig, peak_anno_df, 
                             by.x = "peak", by.y = "peak_id", all.x = TRUE)
      # Step 4: Find overlapping genes
      dar_genes <- unique(dar_annotated$SYMBOL[!is.na(dar_annotated$SYMBOL)])
      deg_genes <- unique(deg_sig$gene)
      common_genes <- intersect(dar_genes, deg_genes)
      if (length(common_genes) > 0) {
        cat("     * Found", length(common_genes), "genes with both DAR and DEG signals\n")
        # Step 5: Create association dataframe
        association_df <- data.frame(
          gene = character(),
          dar_logfc = numeric(),
          deg_logfc = numeric(),
          dar_padj = numeric(),
          deg_padj = numeric(),
          peak_annotation = character(),
          distance_to_tss = numeric(),
          stringsAsFactors = FALSE
        )
        for (gene in common_genes) {
          # Get DAR info for this gene (may have multiple peaks)
          dar_gene_data <- dar_annotated[dar_annotated$SYMBOL == gene & 
                                           !is.na(dar_annotated$SYMBOL), ]
          # Get DEG info for this gene
          deg_gene_data <- deg_sig[deg_sig$gene == gene, ]
          if (nrow(dar_gene_data) > 0 && nrow(deg_gene_data) > 0) {
            # If multiple peaks for same gene, take the one with strongest signal
            if (nrow(dar_gene_data) > 1) {
              dar_gene_data <- dar_gene_data[which.min(dar_gene_data$p_val_adj), ]
            }
            association_df <- rbind(association_df, data.frame(
              gene = gene,
              dar_logfc = dar_gene_data$avg_log2FC[1],
              deg_logfc = deg_gene_data$avg_log2FC[1],
              dar_padj = dar_gene_data$p_val_adj[1],
              deg_padj = deg_gene_data$p_val_adj[1],
              peak_annotation = dar_gene_data$annotation[1],
              distance_to_tss = dar_gene_data$distanceToTSS[1],
              stringsAsFactors = FALSE
            ))
          }
        }
        association_df <- association_df[complete.cases(association_df), ]
        if (nrow(association_df) >= 3) {
          # Step 6: Calculate correlation
          correlation <- cor(association_df$dar_logfc, association_df$deg_logfc, 
                             method = "pearson")
          cor_pvalue <- cor.test(association_df$dar_logfc, association_df$deg_logfc)$p.value
          cat("     * Correlation between DAR and DEG log2FC:", round(correlation, 3), 
              "(p =", round(cor_pvalue, 4), ")\n")
          # Step 7: Create visualizations
          # Scatter plot with correlation
          p_association <- ggplot(association_df, aes(x = dar_logfc, y = deg_logfc)) +
            geom_point(aes(color = -log10(dar_padj + deg_padj)), 
                       alpha = 0.7, size = 3) +
            geom_smooth(method = "lm", se = TRUE, color = "red", alpha = 0.3) +
            scale_color_viridis_c(name = "-log10(combined p)") +
            labs(
              title = paste("DAR-DEG Association:", comp_name),
              subtitle = paste("Correlation =", round(correlation, 3), 
                               "(p =", round(cor_pvalue, 4), ")"),
              x = "DAR Log2 Fold Change (Accessibility)",
              y = "DEG Log2 Fold Change (Expression)"
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(size = 14, face = "bold"),
              plot.subtitle = element_text(size = 12),
              legend.position = "right"
            )
          # Add gene labels for top genes
          if (nrow(association_df) <= 20) {
            p_association <- p_association + 
              ggrepel::geom_text_repel(aes(label = gene), size = 3, max.overlaps = 10)
          }
          # Peak annotation distribution plot
          p_annotation <- ggplot(association_df, aes(x = peak_annotation)) +
            geom_bar(fill = "steelblue", alpha = 0.7) +
            labs(
              title = "Peak Annotation Distribution",
              x = "Genomic Annotation",
              y = "Number of Genes"
            ) +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
          # Distance to TSS distribution
          p_distance <- ggplot(association_df, aes(x = abs(distance_to_tss)/1000)) +
            geom_histogram(bins = 20, fill = "orange", alpha = 0.7) +
            labs(
              title = "Distance to TSS Distribution",
              x = "Distance to TSS (kb)",
              y = "Number of Genes"
            ) +
            theme_minimal()
          # Combined plot
          p_combined <- plot_grid(
            p_association,
            plot_grid(p_annotation, p_distance, ncol = 1),
            ncol = 2, rel_widths = c(2, 1)
          )
          # Step 8: Save results
          comp_name_clean <- gsub("[^A-Za-z0-9_]", "_", comp_name)
          # SAVE plot:
          ggsave(paste0(plot_dir, "DAR_DEG_association_", comp_name_clean, ".pdf"),
                 p_combined, width = 16, height = 10)
          ggsave(paste0(plot_dir, "DAR_DEG_scatter_", comp_name_clean, ".pdf"),
                 p_association, width = 10, height = 8)
          # Save association data
          write.csv(association_df, 
                    paste0(output_dir, "DAR_DEG_association_", comp_name_clean, ".csv"),
                    row.names = FALSE)
          # Save annotated DAR data
          write.csv(dar_annotated,
                    paste0(output_dir, "DAR_annotated_", comp_name_clean, ".csv"),
                    row.names = FALSE)
          # Step 9: Functional enrichment analysis on common genes
          if (length(common_genes) >= 10) {
            cat("     * Performing functional enrichment analysis...\n")
            # GO enrichment
            ego <- enrichGO(gene = common_genes,
                            OrgDb = org.Hs.eg.db,
                            keyType = 'SYMBOL',
                            ont = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.2)
            if (!is.null(ego) && nrow(ego@result) > 0) {
              # Save GO results
              write.csv(ego@result, 
                        paste0(output_dir, "GO_enrichment_", comp_name_clean, ".csv"),
                        row.names = FALSE)
              # Plot GO results
              p_go <- dotplot(ego, showCategory = 10) + 
                ggtitle(paste("GO Enrichment:", comp_name))
              
              ggsave(paste0(plot_dir, "GO_enrichment_", comp_name_clean, ".pdf"),
                     p_go, width = 12, height = 8)
            }
          }
          # Store results
          association_results[[comp_name]] <- list(
            data = association_df,
            annotated_peaks = dar_annotated,
            correlation = correlation,
            correlation_pvalue = cor_pvalue,
            n_genes = nrow(association_df),
            common_genes = common_genes
          )
        } else {
          cat("     * Too few genes for meaningful correlation analysis\n")
        }
      } else {
        cat("     * No common genes found between DARs and DEGs\n")
      }
    }, error = function(e) {
      cat("     * Error in association analysis:", e$message, "\n")
    })
  }
  return(association_results)
}

# Complete association analysis
association_results <- perform_dar_deg_association(dar_results, deg_results)

# ==============================================================================
# 8. Cell type specificity analysis
# ==============================================================================
# comprehensive gene frequency analysis across multiple comparisons
perform_gene_frequency_analysis <- function(dar_results, deg_results, 
                                            association_results = NULL,
                                            txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                                            distance_threshold = 50000,
                                            output_dir = "results/dar_deg/",
                                            plot_dir = "plots/dar_deg/") {
  # Create output directories
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  cat("Performing comprehensive gene frequency analysis...\n")
  # Initialize storage lists
  all_dar_genes <- list()
  all_deg_genes <- list()
  dar_gene_details <- data.frame()
  deg_gene_details <- data.frame()
  # Process DAR results
  if (length(dar_results) > 0) {
    cat("Processing DAR results...\n")
    for (comp_name in names(dar_results)) {
      dar_sig <- dar_results[[comp_name]]$significant
      if (nrow(dar_sig) > 0) {
        cat(paste("  - Processing", comp_name, ":", nrow(dar_sig), "significant peaks\n"))
        # Annotate peaks to get associated genes
        tryCatch({
          # Parse peak coordinates
          peak_gr <- parse_peak_coordinates(dar_sig$peak)
          # Annotate peaks
          peak_anno_df <- annotate_peaks(peak_gr, txdb, distance_threshold)
          # Extract genes (filter out NA values)
          annotated_genes <- peak_anno_df$SYMBOL[!is.na(peak_anno_df$SYMBOL)]
          if (length(annotated_genes) > 0) {
            all_dar_genes[[comp_name]] <- annotated_genes
            # Store detailed information
            comp_details <- data.frame(
              gene = annotated_genes,
              comparison = comp_name,
              peak_id = peak_anno_df$peak_id[!is.na(peak_anno_df$SYMBOL)],
              annotation = peak_anno_df$annotation[!is.na(peak_anno_df$SYMBOL)],
              distance_to_tss = peak_anno_df$distanceToTSS[!is.na(peak_anno_df$SYMBOL)],
              stringsAsFactors = FALSE
            )
            # Merge with DAR significance data
            dar_merged <- merge(comp_details, dar_sig, by.x = "peak_id", by.y = "peak", all.x = TRUE)
            dar_merged$type <- "DAR"
            dar_merged$comparison <- comp_name
            dar_gene_details <- rbind(dar_gene_details, dar_merged)
          }
        }, error = function(e) {
          cat(paste("    Warning: Failed to annotate peaks for", comp_name, ":", e$message, "\n"))
          # Fallback: use simple extraction (though not recommended)
          simple_genes <- gsub(".*-", "", dar_sig$peak)
          simple_genes <- simple_genes[simple_genes != "" & !is.na(simple_genes)]
          if (length(simple_genes) > 0) {
            all_dar_genes[[comp_name]] <- simple_genes
          }
        })
      }
    }
  }
  # Process DEG results
  if (length(deg_results) > 0) {
    cat("Processing DEG results...\n")
    for (comp_name in names(deg_results)) {
      deg_sig <- deg_results[[comp_name]]$significant
      if (nrow(deg_sig) > 0) {
        cat(paste("  - Processing", comp_name, ":", nrow(deg_sig), "significant genes\n"))
        genes <- deg_sig$gene
        all_deg_genes[[comp_name]] <- genes
        # Store detailed information
        comp_details <- deg_sig
        comp_details$comparison <- comp_name
        comp_details$type <- "DEG"
        deg_gene_details <- rbind(deg_gene_details, comp_details)
      }
    }
  }
  # Calculate gene frequencies
  dar_gene_freq <- calculate_gene_frequencies(all_dar_genes, "DAR")
  deg_gene_freq <- calculate_gene_frequencies(all_deg_genes, "DEG")
  # Combine frequency data
  gene_freq_combined <- rbind(dar_gene_freq, deg_gene_freq)
  if (nrow(gene_freq_combined) > 0) {
    cat("Creating visualizations and summary statistics...\n")
    # 1. Basic frequency distribution plot
    p_freq_dist <- create_frequency_distribution_plot(gene_freq_combined)
    # 2. Top genes across comparisons
    p_top_genes <- create_top_genes_plot(gene_freq_combined, top_n = 20)
    # 3. Comparison matrix heatmap
    p_comparison_matrix <- create_comparison_matrix_plot(all_dar_genes, all_deg_genes)
    # 4. Summary statistics table
    summary_stats <- calculate_summary_statistics(dar_gene_freq, deg_gene_freq, 
                                                  all_dar_genes, all_deg_genes)
    # 5. Overlap analysis
    overlap_analysis <- perform_overlap_analysis(all_dar_genes, all_deg_genes)
    # 6. Consistency analysis (genes that appear in multiple comparisons)
    consistency_analysis <- analyze_gene_consistency(gene_freq_combined)
    # Create combined visualization
    p_combined <- create_combined_summary_plot(p_freq_dist, p_top_genes, 
                                               p_comparison_matrix, summary_stats)
    # Save plots
    ggsave(paste0(plot_dir, "gene_frequency_comprehensive_analysis.pdf"), 
           p_combined, width = 20, height = 16)
    ggsave(paste0(plot_dir, "gene_frequency_distribution.pdf"), 
           p_freq_dist, width = 12, height = 8)
    ggsave(paste0(plot_dir, "top_genes_across_comparisons.pdf"), 
           p_top_genes, width = 14, height = 10)
    ggsave(paste0(plot_dir, "comparison_matrix_heatmap.pdf"), 
           p_comparison_matrix, width = 12, height = 10)
    # Save data
    write.csv(gene_freq_combined, 
              paste0(output_dir, "gene_frequency_summary.csv"), row.names = FALSE)
    write.csv(dar_gene_details, 
              paste0(output_dir, "dar_gene_detailed_annotations.csv"), row.names = FALSE)
    write.csv(deg_gene_details, 
              paste0(output_dir, "deg_gene_detailed_results.csv"), row.names = FALSE)
    write.csv(summary_stats, 
              paste0(output_dir, "frequency_analysis_summary_statistics.csv"), row.names = FALSE)
    write.csv(overlap_analysis$overlap_summary, 
              paste0(output_dir, "dar_deg_overlap_analysis.csv"), row.names = FALSE)
    write.csv(consistency_analysis, 
              paste0(output_dir, "gene_consistency_analysis.csv"), row.names = FALSE)
    # Print summary
    print_analysis_summary(summary_stats, overlap_analysis, consistency_analysis)
    return(list(
      gene_frequencies = gene_freq_combined,
      dar_details = dar_gene_details,
      deg_details = deg_gene_details,
      summary_stats = summary_stats,
      overlap_analysis = overlap_analysis,
      consistency_analysis = consistency_analysis
    ))
  } else {
    cat("No significant genes found in any comparison.\n")
    return(NULL)
  }
}
# Helper functions
calculate_gene_frequencies <- function(gene_lists, analysis_type) {
  if (length(gene_lists) == 0) return(data.frame())
  all_genes <- unlist(gene_lists)
  gene_freq <- table(all_genes)
  gene_counts <- data.frame(
    gene = names(gene_freq),
    frequency = as.numeric(gene_freq),
    type = analysis_type,
    stringsAsFactors = FALSE
  )
  # Add percentage
  gene_counts$percentage <- round((gene_counts$frequency / length(gene_lists)) * 100, 2)
  return(gene_counts)
}

create_frequency_distribution_plot <- function(gene_freq_combined) {
  # Create a more sophisticated frequency distribution plot
  max_freq <- max(gene_freq_combined$frequency)
  p <- ggplot(gene_freq_combined, aes(x = frequency, fill = type)) +
    geom_histogram(bins = max_freq, alpha = 0.7, position = "dodge") +
    geom_density(aes(y = ..scaled.. * max(..count..)), alpha = 0.3) +
    scale_fill_brewer(type = "qual", palette = "Set1", name = "Analysis Type") +
    scale_x_continuous(breaks = 1:max_freq) +
    labs(
      title = "Distribution of Gene Frequencies Across Comparisons",
      subtitle = paste("Total genes analyzed: DAR =", 
                       sum(gene_freq_combined$type == "DAR"),
                       ", DEG =", sum(gene_freq_combined$type == "DEG")),
      x = "Number of Comparisons (Frequency)",
      y = "Number of Genes"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "top"
    ) +
    facet_wrap(~type, scales = "free_y")
  return(p)
}

create_top_genes_plot <- function(gene_freq_combined, top_n = 20) {
  # Get top genes by frequency
  top_genes <- gene_freq_combined %>%
    group_by(type) %>%
    top_n(top_n, frequency) %>%
    arrange(type, desc(frequency))
  p <- ggplot(top_genes, aes(x = reorder(gene, frequency), y = frequency, fill = type)) +
    geom_col(alpha = 0.8) +
    coord_flip() +
    scale_fill_brewer(type = "qual", palette = "Set1", name = "Analysis Type") +
    labs(
      title = paste("Top", top_n, "Most Frequent Genes Across Comparisons"),
      x = "Gene",
      y = "Frequency (Number of Comparisons)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 10)
    ) +
    facet_wrap(~type, scales = "free")
  return(p)
}

create_comparison_matrix_plot <- function(all_dar_genes, all_deg_genes) {
  # Create a matrix showing which genes appear in which comparisons
  all_comparisons <- unique(c(names(all_dar_genes), names(all_deg_genes)))
  all_genes_unique <- unique(c(unlist(all_dar_genes), unlist(all_deg_genes)))
  # Create presence/absence matrix
  gene_matrix <- matrix(0, nrow = length(all_genes_unique), ncol = length(all_comparisons))
  rownames(gene_matrix) <- all_genes_unique
  colnames(gene_matrix) <- all_comparisons
  # Fill matrix for DAR genes
  for (comp in names(all_dar_genes)) {
    genes <- all_dar_genes[[comp]]
    gene_matrix[genes[genes %in% all_genes_unique], comp] <- 1
  }
  # Fill matrix for DEG genes (use 2 to distinguish from DAR)
  for (comp in names(all_deg_genes)) {
    genes <- all_deg_genes[[comp]]
    existing_values <- gene_matrix[genes[genes %in% all_genes_unique], comp]
    gene_matrix[genes[genes %in% all_genes_unique], comp] <- 
      ifelse(existing_values == 1, 3, 2)  # 3 = both DAR and DEG
  }
  # Convert to data frame for ggplot
  gene_matrix_df <- expand.grid(Gene = rownames(gene_matrix), 
                                Comparison = colnames(gene_matrix))
  gene_matrix_df$Value <- as.vector(gene_matrix)
  gene_matrix_df$Status <- factor(gene_matrix_df$Value, 
                                  levels = c(0, 1, 2, 3),
                                  labels = c("Not significant", "DAR only", 
                                             "DEG only", "Both DAR & DEG"))
  # Filter to show only genes that appear in at least 2 comparisons
  genes_to_show <- names(which(rowSums(gene_matrix > 0) >= 2))
  gene_matrix_df <- gene_matrix_df[gene_matrix_df$Gene %in% genes_to_show, ]
  
  if (nrow(gene_matrix_df) > 0) {
    p <- ggplot(gene_matrix_df, aes(x = Comparison, y = Gene, fill = Status)) +
      geom_tile(color = "white", size = 0.5) +
      scale_fill_manual(values = c("white", "#FF6B6B", "#4ECDC4", "#45B7D1")) +
      labs(
        title = "Gene-Comparison Matrix",
        subtitle = "Genes appearing in multiple comparisons",
        x = "Comparison",
        y = "Gene"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 14, face = "bold")
      )
  } else {
    p <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, 
               label = "No genes appear in multiple comparisons", 
               size = 6) +
      theme_void()
  }
  return(p)
}

calculate_summary_statistics <- function(dar_gene_freq, deg_gene_freq, 
                                         all_dar_genes, all_deg_genes) {
  stats <- data.frame(
    Metric = character(),
    DAR = numeric(),
    DEG = numeric(),
    stringsAsFactors = FALSE
  )
  # Basic counts
  stats <- rbind(stats, data.frame(
    Metric = "Total unique genes",
    DAR = ifelse(nrow(dar_gene_freq) > 0, nrow(dar_gene_freq), 0),
    DEG = ifelse(nrow(deg_gene_freq) > 0, nrow(deg_gene_freq), 0)
  ))
  stats <- rbind(stats, data.frame(
    Metric = "Total comparisons",
    DAR = length(all_dar_genes),
    DEG = length(all_deg_genes)
  ))
  # Frequency statistics
  if (nrow(dar_gene_freq) > 0) {
    stats <- rbind(stats, data.frame(
      Metric = "Mean frequency per gene",
      DAR = round(mean(dar_gene_freq$frequency), 2),
      DEG = ifelse(nrow(deg_gene_freq) > 0, round(mean(deg_gene_freq$frequency), 2), 0)
    ))
    stats <- rbind(stats, data.frame(
      Metric = "Max frequency",
      DAR = max(dar_gene_freq$frequency),
      DEG = ifelse(nrow(deg_gene_freq) > 0, max(deg_gene_freq$frequency), 0)
    ))
    stats <- rbind(stats, data.frame(
      Metric = "Genes in >1 comparison",
      DAR = sum(dar_gene_freq$frequency > 1),
      DEG = ifelse(nrow(deg_gene_freq) > 0, sum(deg_gene_freq$frequency > 1), 0)
    ))
  }
  return(stats)
}

perform_overlap_analysis <- function(all_dar_genes, all_deg_genes) {
  # Calculate overlaps between DAR and DEG for each comparison
  overlap_results <- data.frame()
  common_comparisons <- intersect(names(all_dar_genes), names(all_deg_genes))
  for (comp in common_comparisons) {
    dar_genes <- all_dar_genes[[comp]]
    deg_genes <- all_deg_genes[[comp]]
    overlap <- intersect(dar_genes, deg_genes)
    overlap_results <- rbind(overlap_results, data.frame(
      comparison = comp,
      dar_genes = length(dar_genes),
      deg_genes = length(deg_genes),
      overlap_genes = length(overlap),
      overlap_percentage = round(length(overlap) / min(length(dar_genes), length(deg_genes)) * 100, 2),
      stringsAsFactors = FALSE
    ))
  }
  return(list(
    overlap_summary = overlap_results,
    common_comparisons = common_comparisons
  ))
}

analyze_gene_consistency <- function(gene_freq_combined) {
  # Analyze genes that appear consistently across multiple comparisons
  consistency_analysis <- gene_freq_combined %>%
    group_by(gene) %>%
    summarise(
      total_appearances = sum(frequency),
      analysis_types = paste(sort(unique(type)), collapse = ", "),
      max_frequency = max(frequency),
      .groups = 'drop'
    ) %>%
    arrange(desc(total_appearances), desc(max_frequency))
  return(consistency_analysis)
}

create_combined_summary_plot <- function(p_freq_dist, p_top_genes, 
                                         p_comparison_matrix, summary_stats) {
  # Create a summary statistics plot
  p_stats <- ggplot(summary_stats, aes(x = Metric)) +
    geom_col(aes(y = DAR, fill = "DAR"), alpha = 0.7, position = "dodge") +
    geom_col(aes(y = DEG, fill = "DEG"), alpha = 0.7, position = "dodge") +
    scale_fill_brewer(type = "qual", palette = "Set1", name = "Type") +
    labs(title = "Summary Statistics", x = "", y = "Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # Combine plots
  p_combined <- plot_grid(
    p_freq_dist,
    plot_grid(p_stats, p_comparison_matrix, ncol = 1),
    ncol = 2, rel_widths = c(1.5, 1)
  )
  return(p_combined)
}

print_analysis_summary <- function(summary_stats, overlap_analysis, consistency_analysis) {
  cat(paste0("\n", paste(rep("=", 60), collapse = ""),'\n'))
  cat("GENE FREQUENCY ANALYSIS SUMMARY\n")
  cat(paste0(paste(rep("=", 60), collapse = ""), "\n"))
  cat("\nBasic Statistics:\n")
  print(summary_stats)
  if (nrow(overlap_analysis$overlap_summary) > 0) {
    cat("\nDAR-DEG Overlap Analysis:\n")
    print(overlap_analysis$overlap_summary)
  }
  cat("\nTop 10 Most Consistent Genes:\n")
  print(head(consistency_analysis, 10))
  cat(paste0(paste(rep("=", 60), collapse = ""), "\n"))
}

frequency_results <- perform_gene_frequency_analysis(
  dar_results = dar_results,
  deg_results = deg_results,
  distance_threshold = 50000,  # 50kb 阈值
  output_dir = "results/dar_deg/",
  plot_dir = "plots/dar_deg/"
)
