#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-07-14
# version:   1.0
# license:   MIT
# brief:     scRNA only cell Type Annotation Script.
#------------------------------------------#

setwd('D:\\Data\\0_yangfan\\0_SingleCell_MultiOmics\\0_Optimized')
# 加载必要的库 / Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(SingleR)
library(celldex)
library(writexl)
library(pheatmap)
library(tidyr)

# 设置参数 / Set parameters
sample_names <- c('ad_sample_rna',"ctl_sample_rna")  # 样品名称 / Sample names
group_names <- c("AD", "CTL")  # 组别 / Group names
data_dir <- "D:/Data/0_yangfan/0_SingleCell_MultiOmics/0_Optimized"  # 数据路径 / Data directory
output_dir <- "snRNA_analysis_output"  # 输出路径 / Output directory
project_name <- "AD_vs_CTL_snRNA"  # 项目名称 / Project name

# 创建输出目录 / Create output directory
dir.create(file.path(data_dir, output_dir, "cell_annotation"), showWarnings = FALSE, recursive = TRUE)
annotation_dir <- file.path(data_dir, output_dir, "cell_annotation")

# 加载整合后的Seurat对象 / Load integrated Seurat object
seurat_obj <- readRDS(file.path(data_dir, output_dir, paste0(project_name, "_integrated.rds")))

cat("=== 细胞类型注释开始 / Cell Type Annotation Started ===\n")
cat("细胞总数 / Total cells:", ncol(seurat_obj), "\n")
cat("基因总数 / Total genes:", nrow(seurat_obj), "\n")

# 1. 基于已知marker基因的初步注释 / Initial annotation based on known marker genes
cat("\n1. 基于marker基因的初步注释 / Initial annotation based on marker genes\n")
# 定义主要细胞类型的marker基因 / Define marker genes for major cell types
marker_genes <- list(
  "Excitatory_Neurons" = c("SLC17A7", 'SLC17A6', 'SLC17A8', "CAMK2A", "NRGN", "RBFOX3", "NEUROD6", "SATB2"),
  "Inhibitory_Neurons" = c("GAD1", "GAD2", "SLC32A1", "PVALB", "SST", "VIP", "LAMP5"),
  "Astrocytes" = c('GFAP','SLC1A2','GLUL','S100B','ALDH1L1','AQP4','IGFBP3','ATP13A4','CBS','SOX9','CD80','CD86','C5AR1'),
  "Oligodendrocytes" = c('OLIG1','OLIG2','OLIG3','CLDN11','MBP','MOG','MAG','GALC','CNP','SOX10','FA2H','UGT8'),
  "Oligodendrocyte_precursor" = c('LHFPL3','MEGF11','PCDH15','PDGFRA','CSPG4','RNF5'),
  "Microglia" = c('P2RY12','ITGAM','CD68','AIF1','CX3CR1','TMEM119','ADGRE1','C1QA','NOS2','TNF','ISYNA1','CCL4','ADORA3','ADRB2','BHLHE41','BIN1','KLF2','NAV3','RHOB','SALL1','SIGLEC8','SPRY1','TAL1'),
  "Endothelial" = c('PECAM1','VWF','A2M','APOLD1','FLT1','TM4SF1','CD34','MCAM','ENG','VCAM1','CDH5')
)

# 检查marker基因在数据中的存在情况 / Check marker gene presence in data
available_markers <- list()
for(cell_type in names(marker_genes)) {
  available_markers[[cell_type]] <- marker_genes[[cell_type]][marker_genes[[cell_type]] %in% rownames(seurat_obj)]
  cat(paste0(cell_type, ": ", length(available_markers[[cell_type]]), "/", length(marker_genes[[cell_type]]), " markers available\n"))
}

# 添加marker基因模块得分 / Add marker gene module scores
for(cell_type in names(available_markers)) {
  if(length(available_markers[[cell_type]]) > 0) {
    seurat_obj <- AddModuleScore(
      object = seurat_obj,
      features = list(available_markers[[cell_type]]),
      name = paste0(cell_type, "_score"),
      ctrl = 5,  # ctrl: 控制基因数量 / Number of control genes [1-100]
      seed = 42  # seed: 随机种子 / Random seed [any integer]
    )
  }
}

# 2. 使用SingleR进行自动注释 / Automatic annotation using SingleR
cat("\n2. 使用SingleR进行自动注释 / Automatic annotation using SingleR\n")

# 获取参考数据集 / Get reference datasets
ref_hpca <- HumanPrimaryCellAtlasData()  # 人类初级细胞图谱 / Human Primary Cell Atlas
ref_blueprint <- BlueprintEncodeData()    # Blueprint/ENCODE数据 / Blueprint/ENCODE data

# 提取表达矩阵 / Extract expression matrix
expr_matrix <- GetAssayData(seurat_obj, assay = "integrated", layer = "data")

# 运行SingleR注释 / Run SingleR annotation
singler_hpca <- SingleR(
  test = expr_matrix,
  ref = ref_hpca,
  labels = ref_hpca$label.main,
  method = "classic",  # method: 分类方法 / Classification method ["classic", "single", "wilcox"]
  clusters = NULL,     # clusters: 聚类信息 / Cluster information [NULL, vector]
  genes = "de",        # genes: 基因选择策略 / Gene selection strategy ["de", "sd", "all"]
  quantile = 0.8,      # quantile: 分位数阈值 / Quantile threshold [0-1]
  fine.tune = TRUE,    # fine.tune: 是否精细调优 / Whether to fine-tune [TRUE, FALSE]
  tune.thresh = 0.05,  # tune.thresh: 调优阈值 / Tuning threshold [0-1]
  prune = TRUE         # prune: 是否修剪低质量注释 / Whether to prune low-quality annotations [TRUE, FALSE]
)

singler_blueprint <- SingleR(
  test = expr_matrix,
  ref = ref_blueprint,
  labels = ref_blueprint$label.main,
  method = "classic",
  clusters = NULL,
  genes = "de",
  quantile = 0.8,
  fine.tune = TRUE,
  tune.thresh = 0.05,
  prune = TRUE
)

# 添加SingleR注释结果到Seurat对象 / Add SingleR results to Seurat object
seurat_obj$SingleR_HPCA <- singler_hpca$labels
seurat_obj$SingleR_Blueprint <- singler_blueprint$labels
seurat_obj$SingleR_HPCA_score <- singler_hpca$scores[cbind(1:nrow(singler_hpca$scores), 
                                                           match(singler_hpca$labels, colnames(singler_hpca$scores)))]
seurat_obj$SingleR_Blueprint_score <- singler_blueprint$scores[cbind(1:nrow(singler_blueprint$scores), 
                                                                     match(singler_blueprint$labels, colnames(singler_blueprint$scores)))]

# 3. 基于聚类结果的手动注释 / Manual annotation based on clustering results
cat("\n3. 基于聚类结果的手动注释 / Manual annotation based on clustering results\n")

# 计算每个聚类的marker基因平均表达 / Calculate average marker gene expression per cluster
cluster_scores <- data.frame(cluster = levels(seurat_obj$seurat_clusters))
for(cell_type in names(available_markers)) {
  score_col <- paste0(cell_type, "_score1")
  if(score_col %in% colnames(seurat_obj@meta.data)) {
    avg_scores <- aggregate(seurat_obj@meta.data[[score_col]], 
                            by = list(seurat_obj$seurat_clusters), 
                            FUN = mean)
    cluster_scores[[cell_type]] <- avg_scores$x
  }
}
# 为每个聚类分配细胞类型
# 创建命名向量，索引为聚类编号
cluster_annotation <- setNames(
  as.character(cluster_scores$cluster),  # 初始值设为聚类编号（字符型）
  nm = cluster_scores$cluster            # 命名向量使用聚类编号作为名字
)

# 基于最高得分进行分配
for(i in 1:nrow(cluster_scores)) {
  # 获取当前行（排除cluster列）并移除NA值
  scores <- unlist(cluster_scores[i, -1])
  scores <- scores[!is.na(scores)]
  if(length(scores) > 0) {
    # 确定最高得分对应的细胞类型
    best_match <- names(which.max(scores))
    # 获取当前行的实际聚类编号
    current_cluster <- cluster_scores$cluster[i]
    # 通过聚类编号（而非行索引）更新注释
    cluster_annotation[as.character(current_cluster)] <- best_match
  }
}

# 创建细胞类型注释向量（确保名称正确）
cluster_annotation <- setNames(rep(NA_character_, nrow(cluster_scores)),nm = as.character(cluster_scores$cluster))

# 基于最高得分分配细胞类型
for(i in 1:nrow(cluster_scores)) {
  scores <- unlist(cluster_scores[i, -1])
  valid_scores <- scores[!is.na(scores) & is.finite(scores)]
  if(length(valid_scores) > 0) {
    best_type <- names(which.max(valid_scores))
    cluster_annotation[as.character(cluster_scores$cluster[i])] <- best_type
  }
}
# 确保注释向量有正确的名称
names(cluster_annotation) <- as.character(cluster_scores$cluster)

# 创建与Seurat对象细胞顺序匹配的注释向量
cell_type_vector <- cluster_annotation[as.character(seurat_obj$seurat_clusters)]
# 添加细胞名称标识符
names(cell_type_vector) <- colnames(seurat_obj)
# 安全添加到Seurat对象
seurat_obj <- AddMetaData(object = seurat_obj,metadata = cell_type_vector,col.name = "cell_type")

# 4. 精细化注释 / Fine-tuning annotation
cat("\n4. 精细化注释 / Fine-tuning annotation\n")
# 结合多种注释方法进行精细化 / Combine multiple annotation methods for refinement
# 创建投票系统 / Create voting system
annotation_votes <- data.frame(
  cell_barcode = colnames(seurat_obj),
  manual = seurat_obj$cell_type,
  singler_hpca = seurat_obj$SingleR_HPCA,
  singler_blueprint = seurat_obj$SingleR_Blueprint,
  stringsAsFactors = FALSE
)

# 简化注释类别 / Simplify annotation categories
simplify_annotation <- function(annotation) {
  annotation <- toupper(annotation)
  case_when(
    grepl("ExNero", annotation) ~ "Excitatory_Neurons",
    grepl("InNero", annotation) ~ "Inhibitory_Neurons",
    grepl("ASTRO", annotation) ~ "Astrocytes",
    grepl("OLIGO", annotation) ~ "Oligodendrocytes",
    grepl("OPC", annotation) ~ "Oligodendrocyte_precursor",
    grepl("MICRO", annotation) ~ "Microglia",
    grepl("ENDOTH", annotation) ~ "Endothelial",
    TRUE ~ "Other"
  )
}

# 应用简化规则 / Apply simplification rules
annotation_votes$manual_simplified <- simplify_annotation(annotation_votes$manual)
annotation_votes$singler_hpca_simplified <- simplify_annotation(annotation_votes$singler_hpca)
annotation_votes$singler_blueprint_simplified <- simplify_annotation(annotation_votes$singler_blueprint)

### 最终注释（基于一致性投票） / Final annotation based on consensus voting
### final_annotation <- apply(annotation_votes[, c("manual_simplified", "singler_hpca_simplified", "singler_blueprint_simplified")], 1, function(x) {
###  x <- x[!is.na(x)]
###  if(length(x) == 0) return("Unknown")
  # 返回出现频率最高的注释 / Return most frequent annotation
###    tbl <- table(x)
###    names(tbl)[which.max(tbl)]
###  })
###  seurat_obj$cell_type_final <- final_annotation

# 最终注释（优先手动注释） / Final annotation prioritizing manual annotation
final_annotation <- apply(annotation_votes[, c("manual_simplified", "singler_hpca_simplified", "singler_blueprint_simplified")], 1, function(x) {
  # 优先使用手动注释（如果存在） / Prioritize manual annotation if available
  if(!is.na(x["manual_simplified"]) && x["manual_simplified"] != "") {
    return(x["manual_simplified"])
  }
  # 移除NA值 / Remove NA values
  valid_x <- x[!is.na(x)]
  # 如果无有效注释 / If no valid annotations
  if(length(valid_x) == 0) {  return("Unknown")  }
  # 返回出现频率最高的注释 / Return most frequent annotation
  tbl <- table(valid_x)
  top_types <- names(tbl)[tbl == max(tbl)]
  # 如果有多个最高频类型，取第一个 / If tie, take first
  return(top_types[1])
})
# 添加到Seurat对象 / Add to Seurat object
seurat_obj$cell_type_final <- final_annotation

#
#
#    由于加入了手动注释，这里先用 手动注释的 cell_type
#
#

# 5. 生成可视化图表 / Generate visualization plots
cat("\n5. 生成可视化图表 / Generate visualization plots\n")

# 5.1 UMAP图显示细胞类型 / UMAP plot showing cell types
p1 <- DimPlot(seurat_obj, group.by = "cell_type", label = TRUE, label.size = 3) +
  ggtitle("UMAP of Cell Type Annotation") +
  theme(legend.position = "bottom")
ggsave(file.path(annotation_dir, "cell_type_umap.pdf"), p1, width = 10, height = 8, dpi = 300)

# 保存UMAP坐标和细胞类型信息 / Save UMAP coordinates and cell type info
umap_data <- data.frame(
  Cell_Barcode = colnames(seurat_obj),
  UMAP_1 = seurat_obj@reductions$umap@cell.embeddings[, 1],
  UMAP_2 = seurat_obj@reductions$umap@cell.embeddings[, 2],
  Cell_Type = seurat_obj$cell_type,
  Sample = seurat_obj$orig.ident,
  Group = seurat_obj$group
)
write.table(umap_data, file.path(annotation_dir, "umap_cell_types.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# 5.2 样本和组别的细胞类型分布 / Cell type distribution by sample and group
p2 <- ggplot(umap_data, aes(x = Cell_Type, fill = Group)) +
  geom_bar(position = "dodge") +
  theme_minimal() +
  labs(title = "Cell Type Distribution by Group",
       x = "Cell Type", y = "Cell Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(annotation_dir, "cell_type_distribution.pdf"), p2, width = 12, height = 6, dpi = 300)
# 保存细胞类型分布数据 / Save cell type distribution data
cell_type_dist <- table(umap_data$Cell_Type, umap_data$Group)
write.table(cell_type_dist, file.path(annotation_dir, "cell_type_distribution.txt"), sep = "\t", quote = FALSE)

# 5.3 Marker基因热图 / Marker gene heatmap
if(length(available_markers) > 0) {
  # 随机选择每个细胞类型的代表细胞 / Randomly select representative cells for each cell type
  set.seed(42)
  cells_to_plot <- c()
  for(ct in unique(seurat_obj$cell_type)) {
    # 获取当前细胞类型的所有细胞标识符
    ct_cells <- colnames(seurat_obj)[which(seurat_obj$cell_type == ct)]
    if(length(ct_cells) > 100) {
      cells_to_plot <- c(cells_to_plot, sample(ct_cells, 100))
    } else {
      cells_to_plot <- c(cells_to_plot, ct_cells)
    }
  }
  # 创建所有marker基因的列表 / Create list of all marker genes
  all_markers <- unique(unlist(available_markers))
  all_markers <- all_markers[all_markers %in% rownames(seurat_obj)]
  if(length(all_markers) > 0) {
    # 检查数据层可用性并选择合适的数据
    assay <- DefaultAssay(seurat_obj)
    available_layers <- Layers(seurat_obj, assay = assay)
    if("data" %in% available_layers) {
      data_layer <- "data"
    } else if("scale.data" %in% available_layers) {
      data_layer <- "scale.data"
    } else if("counts" %in% available_layers) {
      data_layer <- "counts"
    } else {      stop("没有可用的数据层: counts, data, 或 scale.data")    }
    # 创建热图
    p3 <- DoHeatmap(seurat_obj, 
                    features = all_markers,
                    cells = cells_to_plot,
                    group.by = "cell_type",
                    size = 3,
                    angle = 90) +
      ggtitle("Marker基因表达热图")
    ggsave(file.path(annotation_dir, "marker_genes_heatmap.png"), p3, width = 15, height = 10, dpi = 300)
    
    # 保存热图数据
    heatmap_data <- GetAssayData(seurat_obj, assay = assay, layer = data_layer)[all_markers, cells_to_plot]
    
    write.table(as.data.frame(as.matrix(heatmap_data)), 
      file.path(annotation_dir, "marker_genes_heatmap_data.txt"), 
      sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
  }
}

# 5.4 SingleR注释质量评估 / SingleR annotation quality assessment
# 提取 SingleR 评分数据
scores_df <- as.data.frame(singler_hpca$scores)
#colnames(scores_df) <- singler_hpca$labels
scores_df$Cell <- rownames(scores_df)
# 转换为长格式
scores_long <- pivot_longer(scores_df, cols = -Cell, names_to = "CellType", values_to = "Score")
# 创建 ggplot 热图
p4_ggplot <- ggplot(scores_long, aes(x = Cell, y = CellType, fill = Score)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1)) +
  labs(title = "SingleR HPCA Annotation Quality",
       x = "Cells", y = "Cell Types") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())

# 保存 ggplot 版本
ggsave(file.path(annotation_dir, "SingleR_HPCA_quality_ggplot.pdf"), 
       p4_ggplot, width = 12, height = 8, dpi = 300)

# 保存SingleR得分数据 / Save SingleR score data
singler_scores <- as.data.frame(singler_hpca$scores)
singler_scores$predicted_label <- singler_hpca$labels
singler_scores$cell_barcode <- rownames(singler_scores)
write.table(singler_scores, file.path(annotation_dir, "singler_hpca_scores.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# 6. 统计摘要 / Statistical summary
cat("\n6. 生成统计摘要 / Generate statistical summary\n")

# 计算基本统计信息 / Calculate basic statistics
total_cells <- ncol(seurat_obj)
cell_type_counts <- table(seurat_obj$cell_type)
cell_type_percentages <- round(prop.table(cell_type_counts) * 100, 2)

# 按组别统计 / Statistics by group
group_stats <- table(seurat_obj$cell_type, seurat_obj$group)
group_percentages <- round(prop.table(group_stats, margin = 2) * 100, 2)

# 创建统计摘要表 / Create statistical summary table
summary_stats <- data.frame(
  Cell_Type = names(cell_type_counts),
  Total_Cells = as.numeric(cell_type_counts),
  Percentage = as.numeric(cell_type_percentages),
  AD_Cells = as.numeric(group_stats[, "AD"]),
  AD_Percentage = as.numeric(group_percentages[, "AD"]),
  CTL_Cells = as.numeric(group_stats[, "CTL"]),
  CTL_Percentage = as.numeric(group_percentages[, "CTL"]),
  stringsAsFactors = FALSE
)

# 7. 创建Excel报告 / Create Excel report
cat("\n7. 创建Excel报告 / Create Excel report\n")

# 参数设置摘要 / Parameter settings summary
parameter_summary <- data.frame(
  Parameter = c("Total_Cells", "Total_Genes", "Reference_Datasets", "Marker_Gene_Categories", 
                "SingleR_Method", "SingleR_Quantile", "SingleR_Fine_Tune", "Random_Seed"),
  Value = c(total_cells, nrow(seurat_obj), "HPCA + Blueprint", length(available_markers),
            "classic", 0.8, TRUE, 42),
  Description = c("分析的细胞总数", "检测到的基因总数", "SingleR参考数据集", "使用的marker基因类别数",
                  "SingleR分类方法", "SingleR分位数阈值", "是否进行精细调优", "随机种子"),
  stringsAsFactors = FALSE
)

# 创建Excel文件 / Create Excel file
excel_data <- list(
  "Parameter_Summary" = parameter_summary,
  "Cell_Type_Statistics" = summary_stats,
  "Detailed_Annotation" = data.frame(
    Cell_Barcode = colnames(seurat_obj),
    Sample = seurat_obj$orig.ident,
    Group = seurat_obj$group,
    Cluster = seurat_obj$seurat_clusters,
    Cell_Type_Final = seurat_obj$cell_type,
    SingleR_HPCA = seurat_obj$SingleR_HPCA,
    SingleR_Blueprint = seurat_obj$SingleR_Blueprint,
    SingleR_HPCA_Score = round(seurat_obj$SingleR_HPCA_score, 3),
    SingleR_Blueprint_Score = round(seurat_obj$SingleR_Blueprint_score, 3),
    stringsAsFactors = FALSE
  )
)

write_xlsx(excel_data, file.path(annotation_dir, "cell_annotation_report.xlsx"))

# 8. 保存注释后的Seurat对象 / Save annotated Seurat object
cat("\n8. 保存注释后的Seurat对象 / Save annotated Seurat object\n")
saveRDS(seurat_obj, file.path(annotation_dir, paste0(project_name, "_annotated.rds")))

# 输出最终统计信息 / Output final statistics
cat("\n=== 细胞类型注释完成 / Cell Type Annotation Completed ===\n")
cat("总细胞数 / Total cells:", total_cells, "\n")
cat("识别的细胞类型 / Identified cell types:", length(unique(seurat_obj$cell_type)), "\n")
cat("细胞类型分布 / Cell type distribution:\n")
print(summary_stats)

cat("\n输出文件 / Output files:\n")
cat("- 注释后的Seurat对象 / Annotated Seurat object:", file.path(annotation_dir, paste0(project_name, "_annotated.rds")), "\n")
cat("- Excel报告 / Excel report:", file.path(annotation_dir, "cell_annotation_report.xlsx"), "\n")
cat("- 可视化图表 / Visualization plots:", annotation_dir, "\n")
cat("- 源数据文件 / Source data files:", annotation_dir, "\n")

cat("\n脚本执行完成 / Script execution completed!\n")
