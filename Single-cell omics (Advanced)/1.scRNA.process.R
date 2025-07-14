#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-07-14
# version:   1.0
# license:   MIT
# brief:     scRNA data quality control and Seurat v5 integration.
#------------------------------------------#

# Script 1: snRNA data quality control and Seurat v5 integration

# 加载必要的包 / Load required packages
library(Seurat)
library(SeuratWrappers)
library(scDblFinder)  # 多细胞查找 / Find multi cells
library(ggplot2)
library(dplyr)
library(openxlsx)
library(patchwork)

# 设置参数 / Set parameters
sample_names <- c('ad_sample_rna',"ctl_sample_rna")  # 样品名称 / Sample names
group_names <- c("AD", "CTL")  # 组别 / Group names
data_dir <- "/home/yangfan/Data/Data/yangfan/multi_omics/optimized_pipeline/"  # 数据路径 / Data directory
output_dir <- "snRNA_analysis_output"  # 输出路径 / Output directory
project_name <- "AD_vs_CTL_snRNA"  # 项目名称 / Project name

# 创建输出目录 / Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = F)

# 初始化结果记录 / Initialize result tracking
sample_summary <- data.frame(
  Sample = character(),
  Group = character(),
  Initial_Cells = numeric(),
  Initial_Genes = numeric(),
  After_QC_Cells = numeric(),
  After_QC_Genes = numeric(),
  Doublets_Removed = numeric(),
  Median_UMI = numeric(),
  Median_Genes = numeric(),
  Median_Mito_Pct = numeric(),
  stringsAsFactors = FALSE
)

# 存储所有样品对象 / Store all sample objects
seurat_list <- list()

# 逐个处理样品 / Process each sample
for (i in 1:length(sample_names)) {
  cat("Processing sample:", sample_names[i], "\n")
  # 读取数据 / Read data
  data_path <- file.path(data_dir, sample_names[i], "/outs/filtered_feature_bc_matrix")
  counts <- Read10X(data.dir = data_path)
  # 创建Seurat对象 / Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = counts, project = sample_names[i],
    min.cells = 3,  # 至少在3个细胞中表达的基因 / genes expressed in at least 3 cells
    min.features = 200  # 至少表达200个基因的细胞 / cells expressing at least 200 genes
  )
  # 添加样品信息 / Add sample information
  seurat_obj$sample <- sample_names[i]; seurat_obj$group <- group_names[i]
  # 记录初始统计 / Record initial statistics
  initial_cells <- ncol(seurat_obj); initial_genes <- nrow(seurat_obj)
 
  # 计算质控指标 / Calculate QC metrics
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj,
    pattern = "^MT-"  # 线粒体基因模式 / mitochondrial gene pattern
  )
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj,
    pattern = "^RP[SL]"  # 核糖体基因模式 / ribosomal gene pattern
  )
 
  # 双细胞检测 / Doublet detection
  sce <- as.SingleCellExperiment(seurat_obj)
  sce <- scDblFinder(sce,
    samples = "sample",  # 样品列名 / sample column name
    BPPARAM = BiocParallel::SerialParam()  # 并行参数 / parallel parameters
  )
  seurat_obj$doublet_score <- sce$scDblFinder.score
  seurat_obj$is_doublet <- sce$scDblFinder.class == "doublet"
 
  # 可视化质控指标 / Visualize QC metrics
  qc_plots <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
    ncol = 4, pt.size = 0.1
  )
  # 保存质控图 / Save QC plots
  ggsave(filename = file.path(output_dir, paste0(sample_names[i], "_QC_violin.pdf")),
    plot = qc_plots, width = 12, height = 6)
 
  # 特征散点图 / Feature scatter plot
#  scatter_plot <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "is_doublet")
#  ggsave(filename = file.path(output_dir, paste0(sample_names[i], "_scatter.pdf")), plot = scatter_plot, width = 5, height = 4)
  Splot_data=FetchData(seurat_obj, vars = c("nCount_RNA", "percent.mt", "is_doublet"))
  plot_data <- data.frame(feature1 = Splot_data$nCount_RNA, feature2 = Splot_data$percent.mt, group = Splot_data$is_doublet)
  true_count <- sum(plot_data$group == "TRUE", na.rm = TRUE); total_count <- nrow(plot_data); true_percent <- round(true_count / total_count * 100, 2)
  pdf(file = file.path(output_dir, paste0(sample_names[i], "_scatter.pdf")),wi=5,he=4)
  scatter_plot = ggplot(plot_data, aes(x = feature1, y = feature2, color = group, shape = group))+
	  geom_point() +
	  scale_shape_manual(values = c(1, 16), name = "is_doublet") +
	  scale_color_manual(values = c("#63BDC2", "#E27C71"), name = "is_doublet") +
	  theme_classic()+
	  guides(color = guide_legend(override.aes = list(shape = c(1, 16), size = 3))) +
	  annotate("text", # 添加doublet数量和比例标注
           x = max(plot_data$feature1, na.rm = TRUE) * 0.8,
           y = max(plot_data$feature2, na.rm = TRUE) * 0.95,
           label = paste0("Doublets: ", true_count, " (", true_percent, "%)"),
           size = 4, color = "#E27C71", fontface = "bold")
  print(scatter_plot)
  dev.off()
 
  # 设置质控阈值 / Set QC thresholds
  # 可根据数据调整 / Adjust according to data
  min_features <- 200
  max_features <- 5000
  max_mito_pct <- 20
  max_ribo_pct <- 50
 
  # 过滤细胞 / Filter cells
  seurat_obj <- subset(seurat_obj,
    subset = nFeature_RNA > min_features & nFeature_RNA < max_features &
             percent.mt < max_mito_pct & percent.ribo < max_ribo_pct &
             !is_doublet)
 
  # 记录过滤后统计 / Record post-filtering statistics
  after_qc_cells <- ncol(seurat_obj); after_qc_genes <- nrow(seurat_obj); doublets_removed <- sum(sce$scDblFinder.class == "doublet")
  # 计算中位数统计 / Calculate median statistics
  median_umi <- median(seurat_obj$nCount_RNA); median_genes <- median(seurat_obj$nFeature_RNA); median_mito <- median(seurat_obj$percent.mt)
 
  # 更新汇总表 / Update summary table
  sample_summary <- rbind(sample_summary, data.frame(
    Sample = sample_names[i],
    Group = group_names[i],
    Initial_Cells = initial_cells,
    Initial_Genes = initial_genes,
    After_QC_Cells = after_qc_cells,
    After_QC_Genes = after_qc_genes,
    Doublets_Removed = doublets_removed,
    Median_UMI = median_umi,
    Median_Genes = median_genes,
    Median_Mito_Pct = median_mito
  ))
 
  # 保存单个样品的质控数据 / Save individual sample QC data
  qc_data <- seurat_obj@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "is_doublet")]
  write.table(qc_data, sep = "\t", quote = FALSE, row.names = TRUE,
    file = file.path(output_dir, paste0(sample_names[i], "_QC_data.txt"))
  )
  # 保存处理后的对象 / Save processed object
  ### saveRDS(seurat_obj, file = file.path(output_dir, paste0(sample_names[i], "_processed.rds")))
 
  # 添加到列表 / Add to list
  seurat_list[[sample_names[i]]] <- seurat_obj
  cat("Sample", sample_names[i], "processed:", after_qc_cells, "cells retained\n")
}

# 整合分析 / Integration analysis
cat("Starting integration analysis...\n")

# SCTransform标准化 / SCTransform normalization
seurat_list <- lapply(seurat_list, function(x) {
  SCTransform(x, verbose = FALSE,
    method = "glmGamPoi",  # 使用glmGamPoi方法 / use glmGamPoi method
    vars.to.regress = c("percent.mt", "nCount_RNA")  # 回归变量 / variables to regress
  )
})

# 选择整合特征 / Select integration features
features <- SelectIntegrationFeatures(object.list = seurat_list,
  nfeatures = 3000  # 选择特征数量 / number of features to select. Default: 2000. (建议2k-5k). Increase when 数据集大且细胞类型复杂/批次效应强/关注低丰度或特定基因-通路/整合后批次混合不佳（相同细胞类型仍按批次分离）/FindVariableFeatures 的肘部图显示高变基因范围较广
)

# 准备SCT整合 / Prepare SCT integration
seurat_list <- PrepSCTIntegration(object.list = seurat_list,
  anchor.features = features  # 锚点特征 / anchor features
)

# 使用FindIntegrationAnchors + IntegrateData
anchors <- FindIntegrationAnchors(object.list = seurat_list, verbose = FALSE,
  normalization.method = "SCT",
  anchor.features = features
)

# 整合数据 / Integrate data
merged_seurat <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)

# 设置默认assay为integrated
DefaultAssay(merged_seurat) <- "integrated"
# 降维分析 / Dimensionality reduction
merged_seurat <- RunPCA(merged_seurat, verbose = FALSE,
  npcs = 50,  # 主成分数量 / number of principal components
)

# 确定PC数量 / Determine number of PCs
elbow_plot <- ElbowPlot(merged_seurat, ndims = 50)
ggsave(filename = file.path(output_dir, paste0(project_name, "_elbow_plot.pdf")),
  plot = elbow_plot, width = 8, height = 6
)

# 邻域图和聚类 / Neighborhood graph and clustering
PCA_DIMS=1:30
UMAP_DIMS=PCA_DIMS
merged_seurat <- FindNeighbors(merged_seurat, dims = PCA_DIMS,
  reduction = "pca"  # 使用的降维方法 / reduction method to use
)

merged_seurat <- FindClusters(merged_seurat,
  resolution = 0.5,  # 聚类分辨率 / clustering resolution
  algorithm = 1  # 聚类算法 / clustering algorithm: 1= Louvain; 2= Louvain refined; 3= SLM; 4= Leiden
)

# UMAP降维 / UMAP dimensionality reduction
merged_seurat <- RunUMAP(merged_seurat, verbose = FALSE, dims = UMAP_DIMS,
  reduction = "pca"  # 输入降维方法 / input reduction method
)

# 可视化结果 / Visualize results
# 按样品分组 / Group by sample
sample_plot <- DimPlot(merged_seurat,
  reduction = "umap",  # 使用的降维方法 / reduction method to use
  group.by = "sample",  # 分组变量 / grouping variable
  label = TRUE,  # 显示标签 / show labels
  label.size = 3  # 标签大小 / label size
)

# 按组别分组 / Group by condition
group_plot <- DimPlot(merged_seurat,
  reduction = "umap",
  group.by = "group",
  label = TRUE,
  label.size = 3
)

# 按聚类分组 / Group by cluster
cluster_plot <- DimPlot(merged_seurat,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  label.size = 3
)

# 合并图形 / Combine plots
combined_plot <- (sample_plot | group_plot) / cluster_plot
ggsave(filename = file.path(output_dir, paste0(project_name, "_integration_UMAP.pdf")),
  plot = combined_plot, width = 12, height = 10
)

# 保存UMAP坐标 / Save UMAP coordinates
umap_coords <- merged_seurat@reductions$umap@cell.embeddings
umap_data <- data.frame(
  Cell = rownames(umap_coords),
  UMAP_1 = umap_coords[, 1],
  UMAP_2 = umap_coords[, 2],
  Sample = merged_seurat$sample,
  Group = merged_seurat$group,
  Cluster = merged_seurat$seurat_clusters
)
write.table(umap_data,
  file = file.path(output_dir, paste0(project_name, "_UMAP_coordinates.txt")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# 质控特征可视化 / QC feature visualization
qc_feature_plot <- FeaturePlot(merged_seurat,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  ncol = 3, reduction = "umap"
)
ggsave(filename = file.path(output_dir, paste0(project_name, "_QC_features_UMAP.pdf")),
  plot = qc_feature_plot,
  width = 15, height = 5
)

# 保存最终对象 / Save final object
saveRDS(merged_seurat, file = file.path(output_dir, paste0(project_name, "_integrated.rds")))

# 生成报告 / Generate report
# 整合参数 / Integration parameters
dims=paste0(PCA_DIMS[1],":",PCA_DIMS[length(PCA_DIMS)])
integration_params <- data.frame(
  Parameter = c("SCTransform_method", "vars_to_regress", "integration_features",
                "integration_method", "PCA_dims", "clustering_resolution",
                "clustering_algorithm", "UMAP_dims"),
  Value = c("glmGamPoi", "percent.mt,nCount_RNA", "3000", "RPCAIntegration",
            dims, "0.5", "1", dims),
  Description = c("SCTransform标准化方法", "回归变量", "整合特征数量",
                 "整合方法", "PCA维度", "聚类分辨率", "聚类算法", "UMAP维度")
)

# 聚类统计 / Cluster statistics
cluster_stats <- merged_seurat@meta.data %>%
  group_by(seurat_clusters, group) %>%
  summarise(
    n_cells = n(),
    mean_UMI = mean(nCount_RNA),
    mean_genes = mean(nFeature_RNA),
    mean_mito = mean(percent.mt),
    .groups = 'drop'
  )

# 总体统计 / Overall statistics
overall_stats <- data.frame(
  Metric = c("Total_cells_after_integration", "Total_genes", "Number_of_clusters",
             "Number_of_samples", "Number_of_groups"),
  Value = c(ncol(merged_seurat), nrow(merged_seurat),
            length(unique(merged_seurat$seurat_clusters)),
            length(unique(merged_seurat$sample)),
            length(unique(merged_seurat$group)))
)

# 写入Excel文件 / Write to Excel file
wb <- createWorkbook()
addWorksheet(wb, "Sample_Summary")
addWorksheet(wb, "Integration_Parameters")
addWorksheet(wb, "Cluster_Statistics")
addWorksheet(wb, "Overall_Statistics")

writeData(wb, "Sample_Summary", sample_summary)
writeData(wb, "Integration_Parameters", integration_params)
writeData(wb, "Cluster_Statistics", cluster_stats)
writeData(wb, "Overall_Statistics", overall_stats)

saveWorkbook(wb, file = file.path(output_dir, paste0(project_name, "_integration_report.xlsx")), overwrite = TRUE)

cat("snRNA integration analysis completed!\n")
cat("Output files saved in:", output_dir, "\n")
cat("Final integrated object:", paste0(project_name, "_integrated.rds"), "\n")
cat("Excel report:", paste0(project_name, "_integration_report.xlsx"), "\n")

# 输出兼容性文件用于其他软件分析 / Output compatibility files for other software
# 导出表达矩阵用于其他分析 / Export expression matrix for other analysis
# 可用于: SCENIC (基因调控网络), CellChat (细胞通讯), 等
