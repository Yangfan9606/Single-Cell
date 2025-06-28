#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-06-28
# version:   1.0
# license:   MIT
# brief:     Analyze chromatin co-accessibility and identify interacting regulatory elements.
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
library(cicero)
library(BPCells)
library(monocle3)
library(SeuratWrappers)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(Matrix)
library(stringr)
### INSTALL SeuratWrappers : remotes::install_github('satijalab/seurat-wrappers') 
### INSTALL monocle3 : devtools::install_github('cole-trapnell-lab/monocle3', ref="develop")
### INSTALL Cicero for moncle3 : devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")

# =============================================================================
# SingleCell_RNA_ATAC_Integrated: Chromatin Co-Accessibility analysis (cicero)
# =============================================================================
# Sample name
sample_name='CTL'

# Parameters
input_file=str_c(sample_name,'.joint_integrated.rds')
genome_size_file = "reference/hg38.chrom.sizes"
sample_num = 200  # Increase the sampling quantity to improve accuracy.
coaccess_threshold = 0.15  # Co-accessibility threshold
target_genes = c("NRGN", "SLC17A7", "CAMK2A", "RBFOX3", "GAD1", "GAD2", "GPS2", "DHFR", "GFAP", "AQP4", "ALDH1L1", "MBP", "MOG", "PLP1", "VCAN", "PDGFRA", "SOX10", "CSPG4", "CSF1R", "CD74", "AIF1", "P2RY12", "TMEM119", "CD34","PECAM1","PDGFRB")  # Target genes
seed = 12345  # Random seed
set.seed(seed)
# Output dir
output_dirs <- c("results/cicero", "plots/cicero")
sapply(output_dirs, function(x) if(!dir.exists(x)) dir.create(x, recursive = TRUE))

# ==============================================================================
# 1. Data loading and preprocessing.
# ==============================================================================
# Cicero requires a cell_data_set (cds) object, which can be obtained by converting from a Seurat object.
atac_obj <- readRDS(input_file)
# Check if there are UMAP coordinates
if (!"umap" %in% names(atac_obj@reductions)) {
  stop("The ATAC object is missing UMAP coordinates; please run dimensionality reduction analysis first.")
}

# ==============================================================================
# 2. Create Cicero input object.
# ==============================================================================
atac_counts <- GetAssayData(atac_obj, assay = "ATAC", layer = "counts")
atac_meta <- atac_obj@meta.data

# Fix: Compatibility of rowSums calculation with different versions of Matrix.
safe_rowSums <- function(x) {
  tryCatch({    Matrix::rowSums(x > 0)  }, error = function(e1) {
    tryCatch({      base::rowSums(as.matrix(x > 0))    }, error = function(e2) {
      apply(x, 1, function(row) sum(row > 0))    })  })
}

# Filter low-quality cells and peaks. 
min_cells <- ceiling(ncol(atac_counts) * 0.05)  # Retain peaks with signals in at least 5% of cells.
peak_counts <- safe_rowSums(atac_counts)
keep_peaks <- peak_counts >= min_cells
cat("   - Filtering peaks: keeping", sum(keep_peaks), "out of", length(keep_peaks), "peaks\n")
atac_counts <- atac_counts[keep_peaks, ]

# Create peak information data frame.
peak_info <- data.frame(
  site_name = rownames(atac_counts),
  chr = sapply(strsplit(rownames(atac_counts), "-"), `[`, 1),
  bp1 = as.numeric(sapply(strsplit(rownames(atac_counts), "-"), `[`, 2)),
  bp2 = as.numeric(sapply(strsplit(rownames(atac_counts), "-"), `[`, 3))
)
rownames(peak_info) <- peak_info$site_name

# Validate chromosome information.
valid_chrs <- paste0("chr", c(1:22, "X", "Y"))
peak_info <- peak_info[peak_info$chr %in% valid_chrs, ]
atac_counts <- atac_counts[rownames(peak_info), ]
cat("   - Kept", nrow(peak_info), "peaks on standard chromosomes\n")

# Ensure matrix format compatibility.
if (ncol(atac_counts) != nrow(atac_meta)) {
  cat("   - Warning: Cell metadata dimension mismatch. Subsetting to common cells...\n")
  common_cells <- intersect(colnames(atac_counts), rownames(atac_meta))
  atac_counts <- atac_counts[, common_cells]
  atac_meta <- atac_meta[common_cells, ]
}
cat("   - Final dimensions:", nrow(atac_counts), "peaks and", ncol(atac_counts), "cells\n")

# Create cell_data_set (cds) object.
tryCatch({
  input_cds <- new_cell_data_set(
    expression_data = as.matrix(atac_counts),
    cell_metadata = atac_meta,
    gene_metadata = peak_info
  )
  umap_coords <- Embeddings(atac_obj, reduction = "umap")[colnames(input_cds), ]
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
  reducedDims(input_cds)$UMAP <- umap_coords
  cat("   - Successfully created cell_data_set object\n")
}, error = function(e) {  stop("Failed to create cell_data_set object: ", e$message)  })

# ==============================================================================
# 3. Run Cicero analysis.
# ==============================================================================
genome_size <- read.table(genome_size_file, header = FALSE)
colnames(genome_size) <- c("chromosome", "size")

# Create Cicero CDS - Use the correct parameter names.
cat("   - Creating Cicero CDS...\n")
tryCatch({
  cicero_cds <- tryCatch({
    cicero::make_cicero_cds(
      cds = input_cds,
      reduced_coordinates = reducedDims(input_cds)$UMAP,
      k = 50,
      size_factor_normalize = FALSE
    )
  }, error = function(e) {
    cicero::make_cicero_cds(
      input_cds = input_cds,
      reduced_coordinates = reducedDims(input_cds)$UMAP,
      k = 50,
      size_factor_normalize = FALSE
    )
  })
  cat("   - Cicero CDS created successfully\n")
}, error = function(e) {  stop("创建Cicero CDS失败: ", e$message)  })

# Run Cicero
tryCatch({
  conns <- run_cicero(
    cds = cicero_cds, 
    genomic_coords = genome_size_file, 
    sample_num = sample_num,
    window = 500000,  # Sliding window size (default 500kb)
  )
  cat("   - Cicero analysis completed\n")
  cat("   - Found", nrow(conns), "potential connections\n")
}, error = function(e) {  stop("Cicero analysis failed: ", e$message)  })

# SAVE/PLOT: Original connection results
write.csv(conns, str_c("results/cicero/",sample_name,".all_connections.csv"), row.names = FALSE)

# ==============================================================================
# 4. Filtering and annotating significant connections.
# ==============================================================================
significant_conns <- conns %>%
  filter(coaccess >= coaccess_threshold) %>%
  arrange(desc(coaccess))
cat("   - Significant connections (coaccess >=", coaccess_threshold, "):", nrow(significant_conns), "\n")

# Function for calculating distance
calculate_distance <- function(peak1, peak2) {
  p1 <- strsplit(peak1, "-")[[1]]
  p2 <- strsplit(peak2, "-")[[1]]
  if (p1[1] != p2[1]) return(NA_real_)   # Check if the chromosomes are the same.
  start1 <- as.numeric(p1[2]); end1 <- as.numeric(p1[3])
  start2 <- as.numeric(p2[2]); end2 <- as.numeric(p2[3])
  if (end1 < start2) {
    return(start2 - end1)   # Peak1 is on to the left side of Peak2.
  } else if (end2 < start1) {
    return(start1 - end2)   # Peak1 is on the right side of Peak2.
  } else {
    return(0)               # Peak overlap
  }
}
add_peak_info <- function(df) {
  df %>%
    separate(Peak1, into = c("Peak1_chr", "Peak1_start", "Peak1_end"), 
             sep = "-", remove = FALSE, convert = TRUE) %>%
    separate(Peak2, into = c("Peak2_chr", "Peak2_start", "Peak2_end"), 
             sep = "-", remove = FALSE, convert = TRUE) %>%
    mutate(
      Peak1_mid = (Peak1_start + Peak1_end) / 2,
      Peak2_mid = (Peak2_start + Peak2_end) / 2
    ) %>%
    select(Peak1, Peak2, coaccess, distance, 
           Peak1_chr, Peak1_mid, Peak2_chr, Peak2_mid,
           everything())
}

# Main processing flow
significant_conns <- significant_conns %>%
  mutate(Peak1 = as.character(Peak1),
         Peak2 = as.character(Peak2)) %>%
  rowwise() %>%
  mutate(distance = calculate_distance(Peak1, Peak2)) %>%
  ungroup() %>%
  add_peak_info()
# Statistical information
cat("   - Distance statistics:\n")
cat("     * Mean distance:", round(mean(significant_conns$distance)), "bp\n")
cat("     * Median distance:", round(median(significant_conns$distance)), "bp\n")
cat("     * Max distance:", round(max(significant_conns$distance)), "bp\n")

# Plotting the distance distribution map.
p_distance <- ggplot(significant_conns, aes(x = distance/1000)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  labs(
    title = "Distance Distribution of Co-accessible Peaks",
    x = "Distance (kb)",
    y = "Number of Connections"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )
# SAVE/PLOT: Connection distance distribution
ggsave(str_c("plots/cicero/",sample_name,".connection_distance_distribution.pdf"), p_distance, width = 8, height = 6)

# Plotting the distribution of co-accessibility scores.
p_coaccess <- ggplot(significant_conns, aes(x = coaccess)) +
  geom_histogram(bins = 30, fill = "coral", alpha = 0.7) +
  geom_vline(xintercept = coaccess_threshold, linetype = "dashed", color = "red") +
  labs(
    title = "Co-accessibility Score Distribution",
    x = "Co-accessibility Score",
    y = "Number of Connections"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )
# SAVE/PLOT: Co-Access score distribution
ggsave(str_c("plots/cicero/",sample_name,".coaccess_score_distribution.pdf"), p_coaccess, width = 8, height = 6)

# ==============================================================================
# 5. Gene-specific connectivity analysis.
# ==============================================================================
gene_anno <- genes(EnsDb.Hsapiens.v86, filter = AnnotationFilterList(GeneIdFilter("ENSG", "startsWith")))
gene_anno <- keepStandardChromosomes(gene_anno, pruning.mode = "coarse")

# Analyze the connections for each target gene.
gene_connections <- list()

for (gene in target_genes) {
  cat("   - Analyzing connections for", gene, "...\n")
  gene_info <- gene_anno[gene_anno$gene_name == gene]
  if (length(gene_info) == 0) {
    cat("     * Warning: Gene", gene, "not found in annotation\n")
    next
  }
  gene_extended <- resize(gene_info, width = width(gene_info) + 1000000, fix = "center")   # Expanded gene region (±500kb)
  # Find connections that overlap with this gene region
  gene_conns <- significant_conns %>%
    filter(
      (Peak1_chr == as.character(seqnames(gene_extended)) & 
         Peak1_mid >= start(gene_extended) & Peak1_mid <= end(gene_extended)) |
        (Peak2_chr == as.character(seqnames(gene_extended)) & 
           Peak2_mid >= start(gene_extended) & Peak2_mid <= end(gene_extended))
    ) %>%
    arrange(desc(coaccess))
  gene_connections[[gene]] <- gene_conns
  cat("     * Found", nrow(gene_conns), "connections\n")
}

# ==============================================================================
# 6. Visualizing connectivity of specific gene regions.
# ==============================================================================
# Create a connectivity map for each gene
for (gene in names(gene_connections)) {
  if (nrow(gene_connections[[gene]]) == 0) next
  cat("   - Creating connection plot for", gene, "...\n")
  gene_info <- gene_anno[gene_anno$gene_name == gene]
  chr <- as.character(seqnames(gene_info))[1]
  start_pos <- max(1, start(gene_info) - 500000)
  end_pos <- end(gene_info) + 500000
  
  tryCatch({
    p_conn <- plot_connections(
      conns = significant_conns,
      chr = chr,
      minbp = start_pos,
      maxbp = end_pos,
      gene_model = gene_anno,
      coaccess_cutoff = coaccess_threshold,
      connection_width = 0.5,
      collapseTranscripts = TRUE,
      alpha_by_coaccess = TRUE
    ) +
      labs(title = paste("Co-accessibility Connections around", gene)) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.text = element_text(size = 10)
      )
    #SAVE/PLOT:
    ggsave(
      paste0("plots/cicero/",sample_name,".connections_", gene, ".pdf"), 
      p_conn, 
      width = 12, 
      height = 8,
      limitsize = FALSE
    )
    
  }, error = function(e) {
    cat("     * Warning: Failed to create plot for", gene, ":", e$message, "\n")
  })
}

# ==============================================================================
# 7. Save the results.
# ==============================================================================
# SAVE/PLOT: Significant connections
write.csv(significant_conns, str_c("results/cicero/",sample_name,".significant_connections.csv"), row.names = FALSE)
# SAVE/PLOT: Gene-specific connections
saveRDS(gene_connections, str_c("results/cicero/",sample_name,".gene_specific_connections.rds"))

# Statistics summary of connections
connection_summary <- data.frame(
  Gene = names(gene_connections),
  Connections = sapply(gene_connections, nrow),
  Max_Coaccess = sapply(gene_connections, function(x) if(nrow(x) > 0) max(x$coaccess) else 0),
  Mean_Distance = sapply(gene_connections, function(x) if(nrow(x) > 0) round(mean(x$distance)) else 0)
)
# SAVE/PLOT: Connection_summary
write.csv(connection_summary, str_c("results/cicero/",sample_name,".connection_summary.csv"), row.names = FALSE)

# Save the processed Cicero object
cicero_results <- list(
  cds = input_cds,
  all_connections = conns,
  significant_connections = significant_conns,
  gene_connections = gene_connections,
  parameters = params
)
saveRDS(cicero_results, "results/cicero/cicero_results.rds")

# ==============================================================================
# 8. Generate analysis report.
# ==============================================================================
report <- paste0(
  "=== Cicero Co-accessibility Analysis Report ===\n",
  "Analysis Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
  "Input File: ", input_file, "\n\n",
  "Data Summary:\n",
  "- Total cells: ", ncol(atac_counts), "\n",
  "- Total peaks analyzed: ", nrow(atac_counts), "\n", 
  "- Sampling number: ", sample_num, "\n\n",
  "Connection Results:\n",
  "- Total connections found: ", nrow(conns), "\n",
  "- Significant connections (>= ", coaccess_threshold, "): ", nrow(significant_conns), "\n",
  "- Mean connection distance: ", round(mean(significant_conns$distance)), " bp\n",
  "- Median connection distance: ", round(median(significant_conns$distance)), " bp\n\n",
  "Gene-specific Analysis:\n"
)

for (i in 1:nrow(connection_summary)) {
  report <- paste0(report,
                   "- ", connection_summary$Gene[i], ": ", connection_summary$Connections[i], 
                   " connections (max score: ", round(connection_summary$Max_Coaccess[i], 3), ")\n"
  )
}

report <- paste0(report, "\nOutput Files:\n",
                 "- All connections: results/cicero/all_connections.csv\n",
                 "- Significant connections: results/cicero/significant_connections.csv\n",
                 "- Gene-specific connections: results/cicero/gene_specific_connections.rds\n",
                 "- Connection plots: plots/cicero/\n"
)

writeLines(report, "results/cicero/analysis_report.txt")

cat("\n=== Cicero Co-accessibility Analysis Completed ===\n")
cat("Results saved in: results/cicero/\n")
cat("Plots saved in: plots/cicero/\n")

print(connection_summary)
