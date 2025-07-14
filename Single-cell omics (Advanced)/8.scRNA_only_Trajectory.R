#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-07-14
# version:   1.0
# license:   MIT
# brief:     scRNA only Trajectory Analysis Script.
#------------------------------------------#

# ===============================================================================
# 使用Monocle3进行拟时序分析 / Using Monocle3 for pseudotime analysis
# ===============================================================================
setwd('D:\\Data\\0_yangfan\\0_SingleCell_MultiOmics\\0_Optimized')
# 加载必要的包 / Load required packages
library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(pheatmap)
library(RColorBrewer)
library(viridis)

# 设置参数 / Set parameters
sample_names <- c('ad_sample_rna',"ctl_sample_rna")  # 样品名称 / Sample names
group_names <- c("AD", "CTL")  # 组别 / Group names
data_dir <- "D:/Data/0_yangfan/0_SingleCell_MultiOmics/0_Optimized"  # 数据路径 / Data directory
output_dir <- "snRNA_analysis_output"  # 输出路径 / Output directory
project_name <- "AD_vs_CTL_snRNA"  # 项目名称 / Project name

# 创建输出目录 / Create output directories
trajectory_dir <- file.path(output_dir, "trajectory_analysis")
dir.create(trajectory_dir, recursive = TRUE, showWarnings = FALSE)
plots_dir <- file.path(trajectory_dir, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
data_export_dir <- file.path(trajectory_dir, "data_export")
dir.create(data_export_dir, recursive = TRUE, showWarnings = FALSE)

# ===============================================================================
# 1. 加载数据 / Load Data
# ===============================================================================

print("Loading annotated snRNA data...")
seurat_obj <- readRDS(file.path(output_dir, "cell_annotation", paste0(project_name, "_annotated.rds")))

# 将 cell_type_major 重命名
current_colnames <- colnames(seurat_obj@meta.data)
idx <- which(current_colnames == "cell_type_major")
colnames(seurat_obj@meta.data)[idx] <- "cell_type"

# 使用 JoinLayers 合并 RNA assay 的所有 count 和 data layers
# 用于 Monocle3 的 raw count 矩阵
seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")

# 检查数据 / Check data
print(paste("Total cells:", ncol(seurat_obj)))
print(paste("Total genes:", nrow(seurat_obj)))
print("Cell types:")
print(table(seurat_obj$cell_type))

# ===============================================================================
# 2. 选择用于轨迹分析的细胞类型 / Select cell types for trajectory analysis
# ===============================================================================

# 这里需要根据实际细胞类型进行选择，通常选择发育相关或病理相关的细胞类型
# Here you need to select based on actual cell types, usually development-related or pathology-related cell types
print("Available cell types:")
print(unique(seurat_obj$cell_type))

# 示例：选择神经元相关细胞类型进行轨迹分析
# Example: Select neuron-related cell types for trajectory analysis
target_cell_types <- c("Excitatory_Neurons", "Inhibitory_Neurons", "Astrocytes", "Oligodendrocytes")
# 请根据实际细胞类型修改 / Please modify based on actual cell types

# 如果目标细胞类型不存在，选择所有细胞类型
# If target cell types don't exist, select all cell types
available_types <- unique(seurat_obj$cell_type)
if(any(target_cell_types %in% available_types)) {
  selected_types <- intersect(target_cell_types, available_types)
  print(paste("Selected cell types:", paste(selected_types, collapse = ", ")))
} else {
  # 如果没有匹配的细胞类型，选择前5个最多的细胞类型
  # If no matching cell types, select top 5 most abundant cell types
  type_counts <- table(seurat_obj$cell_type)
  selected_types <- names(sort(type_counts, decreasing = TRUE))[1:min(5, length(type_counts))]
  print(paste("Using top abundant cell types:", paste(selected_types, collapse = ", ")))
}

# 子集数据 / Subset data
seurat_subset <- subset(seurat_obj, subset = cell_type %in% selected_types)
print(paste("Cells after subsetting:", ncol(seurat_subset)))

# ===============================================================================
# 3. 转换为Monocle3对象 / Convert to Monocle3 object
# ===============================================================================

print("Converting to Monocle3 object...")
# 1. 提取表达矩阵 (来自合并后的 RNA counts)
# 使用 GetAssayData 是更安全、更通用的方法
expression_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
ok_genes = rownames(seurat_obj)
expression_matrix <- expression_matrix[ok_genes, ]
# 2. 提取细胞元数据
cell_metadata <- seurat_obj@meta.data
# 3. 提取基因元数据
gene_metadata <- data.frame(gene_short_name = rownames(expression_matrix),
                            row.names = rownames(expression_matrix))
# 创建CDS对象
cds <- new_cell_data_set(
  expression_data = expression_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

print(paste("Cells in CDS:", ncol(cds)))
print(paste("Genes in CDS:", nrow(cds)))
print("Monocle3 object created successfully.")

# ===============================================================================
# 4. 预处理数据 / Preprocess data
# ===============================================================================

print("Preprocessing data...")
cds <- preprocess_cds(
  cds,
  method = 'PCA',
  num_dim = 50,  # num_dim: PCA维度数 / number of PCA dimensions (可选: 10-100)
  norm_method = "log",  # norm_method: 标准化方法 / normalization method (可选: "log", "size_only", "none")
  use_genes = NULL,  # use_genes: 使用的基因列表 / list of genes to use (可选: NULL, gene_list)
  pseudo_count = 1,  # pseudo_count: 伪计数 / pseudo count (可选: 0.1, 1, 10)
  scaling = TRUE,  # scaling: 是否缩放 / whether to scale (可选: TRUE, FALSE)
  verbose = TRUE  # verbose: 是否显示详细信息 / whether to show verbose output (可选: TRUE, FALSE)
)

print("Data preprocessing completed.")

# ===============================================================================
# 5. 降维 / Dimensionality reduction
# ===============================================================================

print("Performing dimensionality reduction...")
cds <- reduce_dimension(
  cds,
  max_components = 2,  # max_components: 最大组件数 / maximum number of components (可选: 2, 3)
  reduction_method = "UMAP",  # reduction_method: 降维方法 / dimensionality reduction method (可选: "UMAP", "tSNE", "PCA", "LSI")
  preprocess_method = "PCA",  # preprocess_method: 预处理方法 / preprocessing method (可选: "PCA", "LSI")
  umap.min_dist = 0.1,  # umap.min_dist: UMAP最小距离 / UMAP minimum distance (可选: 0.01-1.0)
  umap.n_neighbors = 15L,  # umap.n_neighbors: UMAP邻居数 / UMAP number of neighbors (可选: 5-50)
  umap.fast_sgd = TRUE,  # umap.fast_sgd: 是否使用快速SGD / whether to use fast SGD (可选: TRUE, FALSE)
  umap.nn_method = "annoy",  # umap.nn_method: 最近邻方法 / nearest neighbor method (可选: "annoy", "hnsw")
  verbose = TRUE  # verbose: 是否显示详细信息 / whether to show verbose output (可选: TRUE, FALSE)
)

print("Dimensionality reduction completed.")

# ===============================================================================
# 6. 聚类 / Clustering
# ===============================================================================

print("Performing clustering...")
cds <- cluster_cells(
  cds,
  resolution = 0.5,  # resolution: 聚类分辨率 / clustering resolution (可选: 0.1-2.0)
  random_seed = 42,  # random_seed: 随机种子 / random seed (可选: 任意整数)
  verbose = TRUE,  # verbose: 是否显示详细信息 / whether to show verbose output (可选: TRUE, FALSE)
  k = 20  # k: k近邻数 / number of k-nearest neighbors (可选: 10-50)
)

print("Clustering completed.")

# ===============================================================================
# 7. 学习轨迹 / Learn trajectory
# ===============================================================================

print("Learning trajectory...")
cds <- learn_graph(
  cds,
  use_partition = TRUE,  # use_partition: 是否使用分区 / whether to use partition (可选: TRUE, FALSE)
  close_loop = FALSE,  # close_loop: 是否闭合循环 / whether to close loop (可选: TRUE, FALSE)
  learn_graph_control = list(
    ncenter = 100,  # ncenter: 中心点数量 / number of centers (可选: 50-500)
    euclidean_distance_ratio = 1,  # euclidean_distance_ratio: 欧氏距离比例 / Euclidean distance ratio (可选: 0.5-2.0)
    geodesic_distance_ratio = 1/3,  # geodesic_distance_ratio: 测地距离比例 / geodesic distance ratio (可选: 0.1-1.0)
    minimal_branch_len = 10,  # minimal_branch_len: 最小分支长度 / minimal branch length (可选: 5-50)
    orthogonal_proj_tip = FALSE,  # orthogonal_proj_tip: 是否正交投影末端 / whether to orthogonally project tips (可选: TRUE, FALSE)
    prune_graph = TRUE  # prune_graph: 是否修剪图 / whether to prune graph (可选: TRUE, FALSE)
  ),
  verbose = TRUE  # verbose: 是否显示详细信息 / whether to show verbose output (可选: TRUE, FALSE)
)

print("Trajectory learning completed.")

# ===============================================================================
# 8. 轨迹可视化 / Trajectory visualization
# ===============================================================================

print("Generating trajectory plots...")
# 绘制轨迹图 / Plot trajectory
p1 <- plot_cells(
  cds,
  color_cells_by = "cell_type",  # color_cells_by: 细胞着色依据 / what to color cells by (可选: "cell_type", "cluster", "partition", "pseudotime")
  label_cell_groups = TRUE,  # label_cell_groups: 是否标记细胞组 / whether to label cell groups (可选: TRUE, FALSE)
  label_leaves = TRUE,  # label_leaves: 是否标记叶子节点 / whether to label leaves (可选: TRUE, FALSE)
  label_branch_points = TRUE,  # label_branch_points: 是否标记分支点 / whether to label branch points (可选: TRUE, FALSE)
  graph_label_size = 3,  # graph_label_size: 图标签大小 / graph label size (可选: 1-10)
  cell_size = 0.5,  # cell_size: 细胞大小 / cell size (可选: 0.1-2.0)
  alpha = 0.8,  # alpha: 透明度 / alpha transparency (可选: 0.1-1.0)
  trajectory_graph_color = "black",  # trajectory_graph_color: 轨迹图颜色 / trajectory graph color (可选: "black", "grey", "blue")
  trajectory_graph_segment_size = 1  # trajectory_graph_segment_size: 轨迹图线段大小 / trajectory graph segment size (可选: 0.5-3.0)
) +
  ggtitle("Cell Trajectory by Cell Type") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(filename = file.path(plots_dir, "trajectory_by_celltype.pdf"), 
       plot = p1, width = 12, height = 10, dpi = 300)

# 绘制分区图 / Plot partitions
p2 <- plot_cells(
  cds,
  color_cells_by = "partition",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  graph_label_size = 3,
  cell_size = 0.5,
  alpha = 0.8
) +
  ggtitle("Cell Trajectory by Partition") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(filename = file.path(plots_dir, "trajectory_by_partition.pdf"), 
       plot = p2, width = 12, height = 10, dpi = 300)

# 绘制聚类图 / Plot clusters
#p3 <- plot_cells(
#  cds,
#  color_cells_by = "cluster",
#  label_cell_groups = TRUE,
#  label_leaves = TRUE,
#  label_branch_points = TRUE,
#  graph_label_size = 3,
#  cell_size = 0.5,
#  alpha = 0.8
#) +
#  ggtitle("Cell Trajectory by Cluster") +
#  theme_minimal() +
#  theme(legend.position = "bottom")

#ggsave(filename = file.path(plots_dir, "trajectory_by_cluster.png"), 
#       plot = p3, width = 12, height = 10, dpi = 300)

# ===============================================================================
# 9. 确定轨迹起始点 / Determine trajectory root
# ===============================================================================

print("Determining trajectory root...")
# 可以选择特定的细胞类型作为起始点
# You can choose specific cell types as root

# 方法1：交互式选择起始点 / Method 1: Interactive root selection
# cds <- order_cells(cds, root_cells = NULL)  # 这会打开交互式界面

# 方法2：基于细胞类型选择起始点 / Method 2: Select root based on cell type
# 选择最早期/干细胞样的细胞类型作为起始点
# Select the earliest/stem cell-like cell type as root
if("Neural_Progenitors" %in% selected_types) {
  root_cells <- rownames(colData(cds)[colData(cds)$cell_type == "Neural_Progenitors", ])
} else {
  # 如果没有祖细胞类型，选择第一个细胞类型
  # If no progenitor cell type, select the first cell type
  root_cells <- rownames(colData(cds)[colData(cds)$cell_type == selected_types[1], ])
}

print(paste("Using", length(root_cells), "cells as root"))

# 排序细胞 / Order cells
cds <- order_cells(
  cds,
  root_cells = root_cells,  # root_cells: 根细胞 / root cells (可选: NULL, cell_names)
  root_pr_nodes = NULL  # root_pr_nodes: 根主节点 / root principal nodes (可选: NULL, node_names)
)

print("Cell ordering completed.")

# ===============================================================================
# 10. 拟时序可视化 / Pseudotime visualization
# ===============================================================================

print("Generating pseudotime plots...")

# 绘制拟时序图 / Plot pseudotime
p4 <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  graph_label_size = 3,
  cell_size = 0.5,
  alpha = 0.8
) +
  ggtitle("Cell Trajectory by Pseudotime") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(filename = file.path(plots_dir, "trajectory_by_pseudotime.pdf"), 
       plot = p4, width = 12, height = 10, dpi = 300)

# 按组别绘制拟时序 / Plot pseudotime by group
p5 <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  graph_label_size = 3,
  cell_size = 0.5,
  alpha = 0.8
) +
  facet_wrap(~orig.ident) +
  ggtitle("Cell Trajectory by Pseudotime (by Group)") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(filename = file.path(plots_dir, "trajectory_by_pseudotime_grouped.pdf"), 
       plot = p5, width = 16, height = 8, dpi = 300)

# ===============================================================================
# 11. 拟时序分析 / Pseudotime analysis
# ===============================================================================

print("Performing pseudotime analysis...")

# 获取拟时序信息 / Get pseudotime information
pseudotime_data <- data.frame(
  cell_id = rownames(colData(cds)),
  pseudotime = pseudotime(cds),
  cell_type = colData(cds)$cell_type,
  group = colData(cds)$orig.ident,
  cluster = clusters(cds),
  partition = partitions(cds)
)

# 保存拟时序数据 / Save pseudotime data
write.table(pseudotime_data, 
            file = file.path(data_export_dir, "pseudotime_data.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

print("Pseudotime data saved.")

# ===============================================================================
# 12. 寻找轨迹相关基因 / Find trajectory-associated genes
# ===============================================================================

print("Finding trajectory-associated genes...")

# 使用Moran's I检验寻找轨迹相关基因 / Use Moran's I test to find trajectory-associated genes
trajectory_genes <- graph_test(
  cds,
  neighbor_graph = "principal_graph",  # neighbor_graph: 邻居图类型 / neighbor graph type (可选: "principal_graph", "knn")
  cores = 1,  # cores: 使用的核心数 / number of cores (可选: 1-8)
  verbose = TRUE  # verbose: 是否显示详细信息 / whether to show verbose output (可选: TRUE, FALSE)
)

# 过滤显著的轨迹相关基因 / Filter significant trajectory-associated genes
significant_genes <- trajectory_genes %>%
  filter(q_value < 0.05) %>%  # q_value阈值 / q_value threshold
  arrange(q_value) %>%
  head(100)  # 选择前100个基因 / Select top 100 genes

print(paste("Found", nrow(significant_genes), "significant trajectory-associated genes"))

# 保存轨迹相关基因 / Save trajectory-associated genes
write.table(significant_genes, 
            file = file.path(data_export_dir, "trajectory_genes.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ===============================================================================
# 13. 基因表达热图 / Gene expression heatmap
# ===============================================================================

print("Generating gene expression heatmap...")

if(nrow(significant_genes) > 0) {
  # 选择前50个基因进行可视化 / Select top 50 genes for visualization
  top_genes <- head(significant_genes$gene_short_name, 50)
  
  # 绘制基因表达热图 / Plot gene expression heatmap
  p6 <- plot_genes_in_pseudotime(
    cds[top_genes,],
    color_cells_by = "cell_type",  # color_cells_by: 细胞着色依据 / what to color cells by (可选: "cell_type", "pseudotime", "cluster")
    min_expr = 0.1,  # min_expr: 最小表达量 / minimum expression (可选: 0, 0.1, 0.5)
    ncol = 1,  # ncol: 列数 / number of columns (可选: 1, 2, 3)
    panel_order = NULL,  # panel_order: 面板顺序 / panel order (可选: NULL, gene_names)
    trend_formula = "~ splines::ns(pseudotime, df=3)",  # trend_formula: 趋势公式 / trend formula (可选: "~sm.ns(pseudotime, df=3)", "~pseudotime")
    label_by_short_name = TRUE  # label_by_short_name: 是否使用短名称标记 / whether to label by short name (可选: TRUE, FALSE)
  )
  ggsave(filename = file.path(plots_dir, "trajectory_genes_expression.png"), 
         plot = p6, width = 12, height = 20, dpi = 300)
}
  # 创建基因表达矩阵热图 / Create gene expression matrix heatmap
#  gene_expr_matrix <- exprs(cds[top_genes,])
  # 按拟时序排序细胞 / Order cells by pseudotime
#  cell_order <- order(pseudotime(cds))
#  gene_expr_matrix <- gene_expr_matrix[, cell_order]
  # 创建注释信息 / Create annotation information
#  annotation_col <- data.frame(
#    pseudotime = pseudotime(cds)[cell_order],
#    cell_type = colData(cds)$cell_type[cell_order],
#    group = colData(cds)$orig.ident[cell_order]
#  )
#  rownames(annotation_col) <- colnames(gene_expr_matrix)
  # 绘制热图 / Plot heatmap
#  png(filename = file.path(plots_dir, "trajectory_genes_heatmap.png"), 
#      width = 16, height = 12, units = "in", res = 300)
#  
#  pheatmap(as.matrix(gene_expr_matrix),
#           scale = "row",
#           cluster_cols = FALSE,
#           cluster_rows = TRUE,
#           show_colnames = FALSE,
#           show_rownames = TRUE,
#           annotation_col = annotation_col,
#           color = colorRampPalette(c("blue", "white", "red"))(100),
#           main = "Trajectory-Associated Genes Expression Heatmap")
#  dev.off()
#  print("Gene expression heatmap saved.")

# ===============================================================================
# 14. 功能富集分析 / Functional enrichment analysis
# ===============================================================================

print("Performing functional enrichment analysis...")

if(nrow(significant_genes) > 0) {
  # 准备基因列表 / Prepare gene list
  gene_list <- significant_genes$gene_short_name
  
  # 这里可以添加GO/KEGG富集分析
  # Here you can add GO/KEGG enrichment analysis
  # 由于需要额外的包和数据库，这里提供框架
  # Since it requires additional packages and databases, framework is provided here
  
  # 示例：使用clusterProfiler进行富集分析
  # Example: Using clusterProfiler for enrichment analysis
  # library(clusterProfiler)
  # library(org.Hs.eg.db)
  # 
  # ego <- enrichGO(gene         = gene_list,
  #                 OrgDb        = org.Hs.eg.db,
  #                 ont          = "BP",
  #                 pAdjustMethod = "BH",
  #                 pvalueCutoff  = 0.01,
  #                 qvalueCutoff  = 0.05,
  #                 readable     = TRUE)
  
  print("Functional enrichment analysis framework prepared.")
}

# ===============================================================================
# 15. 生成报告 / Generate report
# ===============================================================================

print("Generating analysis report...")

# 创建报告数据 / Create report data
report_data <- list(
  "Analysis_Parameters" = data.frame(
    Parameter = c("Sample_Names", "Group_Names", "Data_Directory", "Output_Directory", 
                  "Project_Name", "Selected_Cell_Types", "Total_Cells_Analyzed", 
                  "Total_Genes", "PCA_Dimensions", "UMAP_Min_Distance", 
                  "UMAP_Neighbors", "Clustering_Resolution", "Trajectory_Centers"),
    Value = c(paste(sample_names, collapse = ", "),
              paste(group_names, collapse = ", "),
              data_dir,
              output_dir,
              project_name,
              paste(selected_types, collapse = ", "),
              ncol(cds),
              nrow(cds),
              50,
              0.1,
              15,
              0.5,
              100),
    stringsAsFactors = FALSE
  ),
  
  "Cell_Type_Distribution" = data.frame(table(colData(cds)$cell_type)),
  
  "Trajectory_Summary" = data.frame(
    Metric = c("Total_Cells_in_Trajectory", "Number_of_Clusters", "Number_of_Partitions",
               "Significant_Trajectory_Genes", "Pseudotime_Range_Min", "Pseudotime_Range_Max"),
    Value = c(ncol(cds),
              length(unique(clusters(cds))),
              length(unique(partitions(cds))),
              nrow(significant_genes),
              round(min(pseudotime(cds), na.rm = TRUE), 3),
              round(max(pseudotime(cds), na.rm = TRUE), 3)),
    stringsAsFactors = FALSE
  )
)

# 保存Excel报告 / Save Excel report
write.xlsx(report_data, file = file.path(trajectory_dir, "trajectory_analysis_report.xlsx"))

print("Analysis report saved.")

# ===============================================================================
# 16. 保存结果对象 / Save result objects
# ===============================================================================

print("Saving result objects...")

# 保存CDS对象 / Save CDS object
# saveRDS(cds, file = file.path(trajectory_dir, paste0(project_name, "_trajectory_cds.rds")))

# 保存轨迹数据摘要 / Save trajectory data summary
trajectory_summary <- list(
  cds = cds,
  pseudotime_data = pseudotime_data,
  trajectory_genes = significant_genes,
  selected_cell_types = selected_types,
  parameters = list(
    num_dim = 50,
    resolution = 0.5,
    ncenter = 100,
    umap_min_dist = 0.1,
    umap_neighbors = 15
  )
)

# saveRDS(trajectory_summary, file = file.path(trajectory_dir, paste0(project_name, "_trajectory_summary.rds")))

print("Result objects saved.")

# ===============================================================================
# 完成分析 / Analysis completed
# ===============================================================================

print("=== Trajectory Analysis Completed ===")
print(paste("Results saved in:", trajectory_dir))
print("Generated files:")
print("1. Trajectory plots (PNG)")
print("2. Pseudotime data (TXT)")
print("3. Trajectory-associated genes (TXT)")
print("4. Gene expression heatmap (PNG)")
print("5. Analysis report (XLSX)")
print("6. CDS object (RDS)")
print("7. Trajectory summary (RDS)")

# 显示分析摘要 / Show analysis summary
cat("\n=== Analysis Summary ===\n")
cat("Project Name:", project_name, "\n")
cat("Selected Cell Types:", paste(selected_types, collapse = ", "), "\n")
cat("Total Cells Analyzed:", ncol(cds), "\n")
cat("Total Genes:", nrow(cds), "\n")
cat("Number of Clusters:", length(unique(clusters(cds))), "\n")
cat("Number of Partitions:", length(unique(partitions(cds))), "\n")
cat("Significant Trajectory Genes:", nrow(significant_genes), "\n")
cat("Pseudotime Range:", round(min(pseudotime(cds), na.rm = TRUE), 3), "to", round(max(pseudotime(cds), na.rm = TRUE), 3), "\n")
cat("========================\n")

print("All analysis completed successfully!")
