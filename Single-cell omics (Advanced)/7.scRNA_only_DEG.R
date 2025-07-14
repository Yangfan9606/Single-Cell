#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-07-14
# version:   1.0
# license:   MIT
# brief:     scRNA only Differential Analysis Script.
#------------------------------------------#

setwd('D:\\Data\\0_yangfan\\0_SingleCell_MultiOmics\\0_Optimized')
# 加载所需库 / Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(writexl)
library(RColorBrewer)
library(viridis)
library(ComplexHeatmap)
library(gridExtra)
library(tidyverse)

# 设置参数 / Set parameters
sample_names <- c('ad_sample_rna',"ctl_sample_rna")  # 样品名称 / Sample names
group_names <- c("AD", "CTL")  # 组别 / Group names
data_dir <- "D:/Data/0_yangfan/0_SingleCell_MultiOmics/0_Optimized"  # 数据路径 / Data directory
output_dir <- "snRNA_analysis_output"  # 输出路径 / Output directory
project_name <- "AD_vs_CTL_snRNA"  # 项目名称 / Project name
Log2FC_threshold <- 1
# 创建输出目录 / Create output directories
dir.create(file.path(output_dir, "differential_analysis"), recursive = TRUE, showWarnings = F)
dir.create(file.path(output_dir, "differential_analysis", "figures"), recursive = TRUE, showWarnings = F)
dir.create(file.path(output_dir, "differential_analysis", "tables"), recursive = TRUE, showWarnings = F)
dir.create(file.path(output_dir, "differential_analysis", "source_data"), recursive = TRUE, showWarnings = F)

# 加载已注释的Seurat对象 / Load annotated Seurat object
cat("Loading annotated Seurat object...\n")
seurat_obj <- readRDS(file.path(output_dir, "cell_annotation", paste0(project_name, "_annotated.rds")))

# 检查对象信息 / Check object information
cat("Seurat object dimensions:", dim(seurat_obj), "\n")
cat("Cell types:", unique(seurat_obj@meta.data$cell_type), "\n")
cat("Groups:", unique(seurat_obj@meta.data$group), "\n")

# =============================================================================
# 1. 数据预处理和质量检查 / Data preprocessing and quality check
# =============================================================================
# 添加细胞类型-组别复合标识 / Add cell type-group composite identifier
seurat_obj$celltype_group <- paste0(seurat_obj$cell_type_major, "_", seurat_obj$group)

# 统计每个细胞类型在各组中的细胞数 / Count cells per cell type and group
cell_counts <- table(seurat_obj$cell_type_major, seurat_obj$group)
cat("Cell counts per cell type and group:\n")
print(cell_counts)

# 保存细胞计数表 / Save cell count table
write.table(cell_counts, 
            file = file.path(output_dir, "differential_analysis", "source_data", "cell_counts_per_group.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE)

# 可视化细胞计数 / Visualize cell counts
cell_counts_df <- as.data.frame(cell_counts)
colnames(cell_counts_df) <- c("CellType", "Group", "Count")

p_cell_counts <- ggplot(cell_counts_df, aes(x = CellType, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cell Counts per Cell Type and Group",
       x = "Cell Type", y = "Cell Count") +
  scale_fill_manual(values = c("AD" = "#E31A1C", "CTL" = "#1F78B4"))
ggsave(file.path(output_dir, "differential_analysis", "figures", "cell_counts_barplot.pdf"), 
       p_cell_counts, width = 10, height = 6, dpi = 300)

# =============================================================================
# 2. 差异基因分析 / Differential gene expression analysis
# =============================================================================
cat("Performing differential gene expression analysis...\n")

# 筛选有足够细胞数的细胞类型 / Filter cell types with sufficient cells
min_cells_per_group <- 10  # 最小细胞数阈值 / Minimum cells threshold
valid_celltypes <- names(which(apply(cell_counts >= min_cells_per_group, 1, all)))

cat("Valid cell types for DEG analysis:", valid_celltypes, "\n")

# 存储差异基因结果 / Store DEG results
deg_results <- list()
deg_summary <- data.frame()

# 对每个细胞类型进行差异分析 / Perform differential analysis for each cell type
for (celltype in valid_celltypes) {
  cat(paste("Analyzing cell type:", celltype, "\n"))
  # 设置细胞标识 / Set cell identities
  Idents(seurat_obj) <- "celltype_group"
  # 执行差异基因分析 / Perform differential gene expression analysis
  tryCatch({
    deg_markers <- FindMarkers(
      object = seurat_obj,  # Seurat对象 / Seurat object
      ident.1 = paste0(celltype, "_AD"),  # 实验组标识 / Treatment group identifier
      ident.2 = paste0(celltype, "_CTL"),  # 对照组标识 / Control group identifier
      test.use = "wilcox",  # 统计检验方法 / Statistical test method
      # 可选参数 / Optional parameters: "wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"
      logfc.threshold = Log2FC_threshold,  # log2倍数变化阈值 / Log2 fold change threshold
      # 可选参数 / Optional parameters: 任何正数 / any positive number
      min.pct = 0.1,  # 最小表达比例 / Minimum expression percentage
      # 可选参数 / Optional parameters: 0到1之间的数值 / value between 0 and 1
      min.diff.pct = 0.1,  # 最小差异比例 / Minimum difference percentage
      # 可选参数 / Optional parameters: 0到1之间的数值 / value between 0 and 1
      only.pos = FALSE,  # 是否仅返回上调基因 / Whether to return only upregulated genes
      # 可选参数 / Optional parameters: TRUE, FALSE
      max.cells.per.ident = Inf,  # 每个组别的最大细胞数 / Maximum cells per group
      # 可选参数 / Optional parameters: 任何正整数或Inf / any positive integer or Inf
      random.seed = 42,  # 随机种子 / Random seed
      # 可选参数 / Optional parameters: 任何整数 / any integer
      latent.vars = NULL,  # 潜在变量 / Latent variables
      # 可选参数 / Optional parameters: NULL, 变量名向量 / vector of variable names
      min.cells.feature = 3,  # 基因表达的最小细胞数 / Minimum cells expressing a feature
      # 可选参数 / Optional parameters: 任何正整数 / any positive integer
      min.cells.group = 3,  # 每组的最小细胞数 / Minimum cells per group
      # 可选参数 / Optional parameters: 任何正整数 / any positive integer
      pseudocount.use = 1,  # 伪计数 / Pseudocount
      # 可选参数 / Optional parameters: 任何正数 / any positive number
      fc.name = NULL,  # 倍数变化列名 / Fold change column name
      # 可选参数 / Optional parameters: NULL, 字符串 / string
      base = 2,  # 对数底数 / Logarithm base
      # 可选参数 / Optional parameters: 任何正数 / any positive number
      verbose = TRUE  # 是否显示进度信息 / Whether to show progress
      # 可选参数 / Optional parameters: TRUE, FALSE
    )
    # 添加额外信息 / Add additional information
    deg_markers$gene <- rownames(deg_markers)
    deg_markers$cell_type <- celltype
    deg_markers$direction <- ifelse(deg_markers$avg_log2FC > 0, "Up", "Down")
    # 保存结果 / Save results
    deg_results[[celltype]] <- deg_markers
    # 保存差异基因表 / Save DEG table
    write.table(deg_markers, 
                file = file.path(output_dir, "differential_analysis", "tables", paste0("DEG_", celltype, ".txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    # 统计摘要 / Summary statistics
    n_up <- sum(deg_markers$avg_log2FC > 0 & deg_markers$p_val_adj < 0.05)
    n_down <- sum(deg_markers$avg_log2FC < 0 & deg_markers$p_val_adj < 0.05)
    n_total <- nrow(deg_markers)
    
    deg_summary <- rbind(deg_summary, data.frame(
      CellType = celltype,
      Total_DEGs = n_total,
      Upregulated = n_up,
      Downregulated = n_down,
      Significant_DEGs = n_up + n_down
    ))
    cat(paste("  - Total DEGs:", n_total, "\n"))
    cat(paste("  - Upregulated (FDR<0.05):", n_up, "\n"))
    cat(paste("  - Downregulated (FDR<0.05):", n_down, "\n"))
  }, error = function(e) {
    cat(paste("Error in", celltype, ":", e$message, "\n"))
  })
}

# 保存差异基因摘要 / Save DEG summary
write.table(deg_summary, 
            file = file.path(output_dir, "differential_analysis", "source_data", "DEG_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# =============================================================================
# 3. 可视化差异基因结果 / Visualize DEG results
# =============================================================================
cat("Creating visualizations...\n")
# 3.1 差异基因数量条形图 / DEG count bar plot
p_deg_summary <- ggplot(deg_summary, aes(x = CellType)) +
  geom_bar(aes(y = Upregulated), stat = "identity", fill = "#E31A1C", alpha = 0.7) +
  geom_bar(aes(y = -Downregulated), stat = "identity", fill = "#1F78B4", alpha = 0.7) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Differential Gene Expression Summary",
       x = "Cell Type", y = "Number of DEGs") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5)

ggsave(file.path(output_dir, "differential_analysis", "figures", "DEG_summary_barplot.pdf"), 
       p_deg_summary, width = 6, height = 4, dpi = 300)

# 3.2 为每个细胞类型创建火山图 / Create volcano plots for each cell type
for (celltype in names(deg_results)) {
  deg_data <- deg_results[[celltype]]
  # 火山图 / Volcano plot
  p_volcano <- EnhancedVolcano(
    toptable = deg_data,  # 差异基因表 / DEG table
    lab = deg_data$gene,  # 基因标签 / Gene labels
    x = 'avg_log2FC',  # X轴变量 / X-axis variable
    y = 'p_val_adj',  # Y轴变量 / Y-axis variable
    title = paste('Volcano Plot:', celltype),  # 图标题 / Plot title
    pCutoff = 0.05,  # p值阈值 / P-value threshold
    FCcutoff = Log2FC_threshold,  # 倍数变化阈值 / Fold change threshold
    pointSize = 2.0,  # 点大小 / Point size
    labSize = 3.0,  # 标签大小 / Label size
    colAlpha = 0.7,  # 颜色透明度 / Color transparency
    # 可选参数 / Optional parameters: 0到1之间的数值 / value between 0 and 1
    legendPosition = 'right',  # 图例位置 / Legend position
    # 可选参数 / Optional parameters: 'top', 'bottom', 'left', 'right', 'none'
    legendLabSize = 12,  # 图例标签大小 / Legend label size
    legendIconSize = 4.0,  # 图例图标大小 / Legend icon size
    drawConnectors = TRUE,  # 是否绘制连接线 / Whether to draw connectors
    widthConnectors = 0.75,  # 连接线宽度 / Connector width
    max.overlaps = 10  # 最大重叠数 / Maximum overlaps
  )
  ggsave(file.path(output_dir, "differential_analysis", "figures", 
                   paste0("volcano_plot_", celltype, ".png")), 
         p_volcano, width = 10, height = 6, dpi = 300)
  # 保存火山图数据 / Save volcano plot data
  volcano_data <- data.frame(
    gene = deg_data$gene,
    log2FC = deg_data$avg_log2FC,
    pvalue = deg_data$p_val,
    adj_pvalue = deg_data$p_val_adj,
    significant = deg_data$p_val_adj < 0.05 & abs(deg_data$avg_log2FC) > Log2FC_threshold
  )
  write.table(volcano_data, 
              file = file.path(output_dir, "differential_analysis", "source_data", 
                               paste0("volcano_data_", celltype, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# =============================================================================
# 4. 功能富集分析 / Functional enrichment analysis
# =============================================================================
cat("Performing functional enrichment analysis...\n")
# 存储富集结果 / Store enrichment results
enrichment_results <- list()
for (celltype in names(deg_results)) {
  deg_data <- deg_results[[celltype]]
  # 获取显著差异基因 / Get significant DEGs
  sig_genes <- deg_data[deg_data$p_val_adj < 0.05 & abs(deg_data$avg_log2FC) > 0.25, "gene"]
  if (length(sig_genes) > 10) {  # 至少10个基因才进行富集分析 / At least 10 genes for enrichment
    tryCatch({
      # GO富集分析 / GO enrichment analysis
      ego <- enrichGO(
        gene = sig_genes,  # 基因列表 / Gene list
        OrgDb = org.Hs.eg.db,  # 物种数据库 / Species database
        keyType = 'SYMBOL',  # 基因ID类型 / Gene ID type
        # 可选参数 / Optional parameters: 'ENTREZID', 'ENSEMBL', 'SYMBOL', 'GENENAME'
        ont = "BP",  # GO本体类型 / GO ontology type
        # 可选参数 / Optional parameters: "BP", "MF", "CC", "ALL"
        pAdjustMethod = "BH",  # 多重检验校正方法 / Multiple testing correction method
        # 可选参数 / Optional parameters: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
        pvalueCutoff = 0.05,  # p值阈值 / P-value threshold
        qvalueCutoff = 0.2,  # q值阈值 / Q-value threshold
        # 可选参数: q=0.2(控制假阳性<20%), q=0.3(手动筛选后续结果), q=0.05(确保高置信度)
        readable = TRUE,  # 是否转换基因ID / Whether to convert gene IDs
        minGSSize = 10,  # 最小基因集大小 / Minimum gene set size
        maxGSSize = 500  # 最大基因集大小 / Maximum gene set size
      )
      if (nrow(ego@result) > 0) {
        enrichment_results[[celltype]] <- ego
        # 保存富集结果 / Save enrichment results
        write.table(ego@result, 
        file = file.path(output_dir, "differential_analysis", "tables", 
                         paste0("GO_enrichment_", celltype, ".txt")), sep = "\t", quote = FALSE, row.names = FALSE)
        # 绘制富集结果 / Plot enrichment results
        p_enrich <- dotplot(ego, showCategory = 20) + 
          ggtitle(paste("GO Enrichment:", celltype))
        
        ggsave(file.path(output_dir, "differential_analysis", "figures", 
                         paste0("GO_enrichment_", celltype, ".png")), 
               p_enrich, width = 8, height = 6, dpi = 300)
      }
      
    }, error = function(e) {
      cat(paste("Error in enrichment analysis for", celltype, ":", e$message, "\n"))
    })
  }
}

# =============================================================================
# 5. 热图可视化 / Heatmap visualization
# =============================================================================
cat("Creating heatmap visualization...\n")
# 获取所有显著差异基因 / Get all significant DEGs
all_sig_genes <- c()
for (celltype in names(deg_results)) {
  deg_data <- deg_results[[celltype]]
  sig_genes <- deg_data[deg_data$p_val_adj < 0.05 & abs(deg_data$avg_log2FC) > Log2FC_threshold, "gene"]
  all_sig_genes <- c(all_sig_genes, sig_genes)
}

# 去重并选择前50个基因 / Remove duplicates and select top 50 genes
all_sig_genes <- unique(all_sig_genes)
if (length(all_sig_genes) > 50) {
  all_sig_genes <- all_sig_genes[1:50]
}

if (length(all_sig_genes) > 0) {
  # 计算平均表达 / Calculate average expression
  avg_exp <- AverageExpression(
    object = seurat_obj,  # Seurat对象 / Seurat object
    features = all_sig_genes,  # 特征基因 / Feature genes
    group.by = "celltype_group",  # 分组变量 / Grouping variable
    slot = "data",  # 数据槽 / Data slot
    # 可选参数 / Optional parameters: "counts", "data", "scale.data"
    assays = NULL,  # 分析类型 / Assays
    # 可选参数 / Optional parameters: NULL, 字符串向量 / character vector
    return.seurat = FALSE,  # 是否返回Seurat对象 / Whether to return Seurat object
    add.ident = NULL,  # 添加标识 / Add identifier
    # 可选参数 / Optional parameters: NULL, 字符串向量 / character vector
    verbose = TRUE  # 是否显示进度 / Whether to show progress
  )
  # 准备热图数据 / Prepare heatmap data
  heatmap_data <- avg_exp$integrated
  # 绘制热图 / Create heatmap
  p_heatmap <- pheatmap(
    mat = as.matrix(heatmap_data),
    scale = "row",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    main = "Top 50 Differential Genes Heatmap",
    fontsize = 16,
    fontsize_row = 10,
    fontsize_col = 12
  )
  pdf(file=file.path(output_dir, "differential_analysis", "figures", "DEG_heatmap.pdf"),wi=12,he=10)
  p_heatmap
  dev.off()
  # 保存热图数据 / Save heatmap data
  write.table(heatmap_data, 
              file = file.path(output_dir, "differential_analysis", "source_data", "heatmap_data.txt"),
              sep = "\t", quote = FALSE, row.names = TRUE)
  
  ###  补充热图
  # 数据预处理：行标准化 (Z-score)
  scaled_data <- heatmap_data %>% 
    t() %>%  # 转置以便按行计算
    scale() %>%  # scale(center = TRUE, scale = TRUE)
    t() %>%  # 转置回原始方向
    as.data.frame()
  # 层次聚类（行）
  row_clust <- hclust(
    d = dist(scaled_data, method = "euclidean"), # 距离计算方法
    method = "complete"  # 聚类方法（可选：single, complete, average, centroid等）
  )
  # 层次聚类（列）
  col_clust <- hclust(
    d = dist(t(scaled_data), method = "euclidean"),
    method = "complete"
  )
  # 重组数据为ggplot长格式
  heatmap_long <- scaled_data %>%
    rownames_to_column(var = "Gene") %>%  # 添加基因名列
    pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression"
    ) %>%
    mutate(
      Gene = factor(Gene, levels = row_clust$labels[row_clust$order]), # 按聚类排序基因
      Sample = factor(Sample, levels = colnames(heatmap_data))  # 保持样品原始顺序（不聚类）
      #Sample = factor(Sample, levels = col_clust$labels[col_clust$order]) # 按聚类排序样本
    )
  # 创建颜色梯度
  color_palette <- colorRampPalette(
    colors = c("#2166AC", "white", "#B2182B"),  # 颜色渐变范围
    space = "rgb"  # 颜色空间（可选：rgb, Lab）
  )(100)

  #传统蓝-白-红优化版
  color_palette <- colorRampPalette( 
    c("#2166AC", "#4575b4", "#74add1", "#abd9e9", "#e0f3f8", 
      "#ffffbf", "#fee090", "#fdae61", "#f46d43", "#d73027", "#a50026")
  )(100)
  
  # 绘制热图
  heatmap_plot <- ggplot(data = heatmap_long, mapping = aes(x = Sample, y = Gene, fill = Expression)
  ) +
    geom_tile(
      color = 'grey90',  # Cell边框颜色（设为NA表示无边框）
      width = 1, height = 1
    ) +
    scale_fill_gradientn(colours = color_palette, 
      name = "Average\nExpression",  # 图例标题
      na.value = "grey50"  # NA值颜色
    ) +
    labs(title = "Top 50 Differential Genes Heatmap", x = NULL, y = NULL) +
    theme_minimal(base_size = 14) +  # 基础字体大小
    theme(
      plot.title = element_text(size = 16,face = "bold", hjust = 0.5  # 标题居中
      ),
      axis.text.x = element_text(size = 12, angle = 45, face = 'bold',color='black', 
        hjust = 1,   # 水平对齐（0左1右）
        vjust = 1    # 垂直对齐（0上1下）
      ),
      axis.text.y = element_text(size = 10,face = 'bold',color='black'),
      legend.position = 'right',  # 图例位置（可选：none, left, top, bottom）
      legend.key.height = unit(2, "cm"),  # 图例高度
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 16,face = 'bold',color='black')
    )
  heatmap_plot
  ggsave(filename = file.path(output_dir, "differential_analysis", "figures", "top50_DEG_heatmap.png"),
         plot = heatmap_plot, width = 12, height = 10, dpi=300)
}

# =============================================================================
# 6. 生成综合报告 / Generate comprehensive report
# =============================================================================
cat("Generating comprehensive report...\n")
# 参数设置报告 / Parameter settings report
parameters_report <- data.frame(
  Parameter = c("Minimum cells per group", "Log2FC threshold", "Minimum expression percentage", 
                "P-value adjustment method", "Statistical test", "Minimum gene set size"),
  Value = c(min_cells_per_group, Log2FC_threshold, 0.1, "BH", "wilcox", 10),
  Description = c("Minimum cells required per group for analysis", 
                  "Minimum log2 fold change for DEG calling",
                  "Minimum percentage of cells expressing a gene",
                  "Multiple testing correction method",
                  "Statistical test used for differential expression",
                  "Minimum gene set size for enrichment analysis")
)

# 分析摘要 / Analysis summary
analysis_summary <- data.frame(
  Metric = c("Total cell types analyzed", "Total DEGs identified", "Cell types with >50 DEGs", 
             "Average DEGs per cell type", "Enrichment analyses performed"),
  Value = c(length(valid_celltypes), 
            sum(deg_summary$Total_DEGs), 
            sum(deg_summary$Significant_DEGs > 50),
            round(mean(deg_summary$Significant_DEGs), 2),
            length(enrichment_results))
)

# 保存Excel报告 / Save Excel report
report_list <- list(
  "Parameters" = parameters_report,
  "Analysis_Summary" = analysis_summary,
  "DEG_Summary" = deg_summary,
  "Cell_Counts" = as.data.frame(cell_counts)
)

write_xlsx(report_list, 
           path = file.path(output_dir, "differential_analysis", "Differential_Analysis_Report.xlsx"))

# 保存分析参数到txt文件 / Save analysis parameters to txt file
writeLines(c(
  paste("Analysis Date:", Sys.Date()),
  paste("Project Name:", project_name),
  paste("Input File:", file.path(output_dir, "cell_annotation", paste0(project_name, "_annotated.rds"))),
  paste("Minimum cells per group:", min_cells_per_group),
  paste("Log2FC threshold:", 0.25),
  paste("Statistical test:", "wilcox"),
  paste("P-value adjustment:", "BH"),
  paste("Valid cell types:", paste(valid_celltypes, collapse = ", ")),
  paste("Total DEGs identified:", sum(deg_summary$Total_DEGs)),
  paste("Significant DEGs (FDR<0.05):", sum(deg_summary$Significant_DEGs))
), file.path(output_dir, "differential_analysis", "analysis_parameters.txt"))

cat("Differential analysis completed successfully!\n")
cat("Results saved to:", file.path(output_dir, "differential_analysis"), "\n")
cat("- Tables: DEG results for each cell type\n")
cat("- Figures: Volcano plots, summary plots, heatmaps\n")
cat("- Source data: Raw data for all plots\n")
cat("- Report: Comprehensive Excel report\n")
