#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-07-14
# version:   2.0
# license:   MIT
# brief:     Advanced scATAC-seq Integration.
#------------------------------------------#

setwd('D:\\Data\\0_yangfan\\0_SingleCell_MultiOmics\\0_Optimized')
# Required libraries
library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(EnsDb.Hsapiens.v86)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(circlize)
library(gridExtra)
library(viridis)
library(scales)
library(IRanges)
library(rtracklayer)
library(ChIPseeker)
library(ReactomePA)
library(enrichplot)
library(DOSE)
library(org.Hs.eg.db)

# Create output directory structure
# 创建输出目录结构
create_output_directories <- function(base_dir = "scATAC_output") {
  dirs <- c(
    file.path(base_dir, "plots"),
    file.path(base_dir, "source_data"),
    file.path(base_dir, "integration_results")
  )
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      cat("Created directory:", dir, "\n")}}
  return(base_dir)
}

# Define data parameters / 定义数据参数
data_dir <- "D:\\Data\\0_yangfan\\0_SingleCell_MultiOmics\\0_Optimized"
samples <- c("ad_sample_atac", "ctl_sample_atac")
conditions <- c("AD", "CTL")
output_dir <- "snATAC_analysis_output"

data_dir = data_dir
samples = samples
conditions = conditions
output_dir = output_dir
min_cells = 10
min_features = 200
ncount_min = 200
ncount_max = 50000
pct_reads_in_peaks_min = 10
blacklist_ratio_max = 0.025
nucleosome_signal_max = 2
tss_enrichment_min = 1
harmony_theta = 2
resolution = 0.5

# # run_advanced_snATAC_integration
dims_harmony = 2:30
dims_neighbors = 2:30

  cat("Starting advanced snATAC-seq integration pipeline...\n")
  cat("开始高级snATAC-seq整合流水线...\n")
  # Create output directories
  create_output_directories(output_dir)
  # Load annotation
  cat("Loading genome annotations...\n")
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- "UCSC"
  # Initialize list to store sample objects / 初始化列表以存储样本对象
  sample_objects <- list()
  
  # Process each sample
  for (i in seq_along(samples)) {
    sample <- samples[i]
    condition <- conditions[i]
    cat("Processing sample:", sample, "with condition:", condition, "\n")
    # Load data using the specified input format
    # 使用指定的输入格式加载数据
    counts <- Read10X_h5(filename = file.path(data_dir, sample, "outs/filtered_peak_bc_matrix.h5"))
    metadata <- read.csv(
      file.path(data_dir, sample, "outs/singlecell.csv"),
      header = TRUE, 
      row.names = 1
    )
    # Create ChromatinAssay object / 创建ChromatinAssay对象
    chrom_assay <- CreateChromatinAssay(
      counts = counts,  # count matrix / count矩阵
      sep = c(":", "-"),  # separators ["-", ":", "_"] / 分隔符
      fragments = file.path(data_dir, sample, "outs/corrected_fragments.tsv.gz"),  # fragments file / 片段文件
      annotation = annotations  # genome annotations / 基因组注释
    )
    # Create Seurat object / 创建Seurat对象
    sample_obj <- CreateSeuratObject(
      counts = chrom_assay,
      assay = "peaks",
      meta.data = metadata
    )
    # Add condition information / 添加条件信息
    sample_obj$condition <- condition
    sample_obj$sample <- sample
    # Store in list
    sample_objects[[sample]] <- sample_obj
  }
## After each sample done ...
  # Merge all samples / 合并所有样本
  cat("Merging samples...\n")
  if (length(sample_objects) == 1) {
    combined_obj <- sample_objects[[1]]
  } else {
    combined_obj <- merge(
      x = sample_objects[[1]],
      y = sample_objects[-1],
      add.cell.ids = names(sample_objects)
    )
  }
  
  # Calculate QC metrics
  cat("Calculating QC metrics...\n")
  combined_obj <- NucleosomeSignal(combined_obj)
  combined_obj <- TSSEnrichment(combined_obj,fast=FALSE)  # fast=FALSE Computes enrichment matrix
  # Add additional QC metrics / 添加额外的质量控制指标
  combined_obj$pct_reads_in_peaks <- combined_obj$peak_region_fragments / combined_obj$passed_filters * 100
  
  if ("blacklist_region_fragments" %in% colnames(combined_obj@meta.data)) {
    combined_obj$blacklist_ratio <- combined_obj$blacklist_region_fragments / 
      combined_obj$peak_region_fragments
  } else {
    combined_obj$blacklist_ratio <- 0
  }
  
  # Save pre-filtering QC metrics / 保存过滤前的质量控制指标
  write.csv(combined_obj@meta.data, file.path(output_dir, "source_data", "pre_filter_metadata.csv"), row.names = TRUE)
  # Create pre-filtering QC dashboard / 创建过滤前的质量控制仪表板
  cat("Creating pre-filtering QC dashboard...\n")
#####  create_comprehensive_qc_dashboard 
    object = combined_obj
    group.by = "condition"
    output_dir = output_dir

    cat("  Creating comprehensive QC dashboard...\n")
    cat("  创建综合质量控制仪表板...\n")
    # 1. Fragment length distribution / 片段长度分布
    cat("  Creating fragment length distribution plot...\n")
    p1 <- FragmentHistogram(object = object, group.by = group.by) +
      ggtitle("Fragment Length Distribution\n片段长度分布") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    # 2. TSS enrichment plot / TSS富集图
    cat("  Creating TSS enrichment plot...\n")
#    DefaultAssay(object) <- "ATAC"
    p2 <- TSSPlot(object = object, group.by = group.by) +
      ggtitle("TSS Enrichment Profile\nTSS富集谱") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    # 3. Quality metrics violin plots / 质量指标小提琴图
    cat("  Creating quality metrics violin plots...\n")
    qc_metrics <- c("nCount_peaks", "nFeature_peaks", "nucleosome_signal", "TSS.enrichment")
    p3 <- VlnPlot(
      object = object,
      features = qc_metrics,
      group.by = group.by,
      ncol = 2,
      pt.size = 0.1
    ) + 
      plot_annotation(title = "Quality Metrics by Condition\n按条件分组的质量指标")
    # 4. Peak region fragments / Peak区域片段数
    cat("  Creating peak region fragments plot...\n")
    p4 <- VlnPlot(
      object = object,
      features = "peak_region_fragments",
      group.by = group.by,
      pt.size = 0.1
    ) + 
      ggtitle("Peak Region Fragments\nPeak区域片段数") +
      theme(plot.title = element_text(hjust = 0.5))
    # 5. Blacklist ratio / 黑名单比例
    if ("blacklist_ratio" %in% colnames(object@meta.data)) {
      cat("    Creating blacklist ratio plot...\n")
      p5 <- VlnPlot(
        object = object,
        features = "blacklist_ratio",
        group.by = group.by,
        pt.size = 0.1
      ) + 
        ggtitle("Blacklist Ratio\n黑名单比例") +
        theme(plot.title = element_text(hjust = 0.5))
    } else {
      p5 <- ggplot() + theme_void()
    }
    # Combine plots
    combined_qc <- (p1 | p2) / (p3) / (p4 | p5)
    # Save QC dashboard
    ggsave(file.path(output_dir, "plots", "qc_dashboard.pdf"), combined_qc, width = 14, height = 16, dpi = 300)
    # Save QC metrics source data
    qc_data <- object@meta.data[, c(group.by, qc_metrics, "peak_region_fragments")]
    if ("blacklist_ratio" %in% colnames(object@meta.data)) {
      qc_data$blacklist_ratio <- object@meta.data$blacklist_ratio
    }
    write.csv(qc_data, file.path(output_dir, "source_data", "qc_metrics_data.csv"), row.names = TRUE)
    cat("  QC dashboard created successfully!\n")
  
#################################################################################
#
#################################################################################
  # Filter cells based on QC metrics / 基于质量控制指标过滤细胞
  cat("Filtering cells...\n")
  cat("Before filtering:", ncol(combined_obj), "cells\n")
  ncol(object)
  combined_obj <- subset(
    object,
    subset = nCount_peaks > ncount_min & 
      nCount_peaks < ncount_max & 
      pct_reads_in_peaks > pct_reads_in_peaks_min & 
      blacklist_ratio < blacklist_ratio_max & 
      nucleosome_signal < nucleosome_signal_max & 
      TSS.enrichment > tss_enrichment_min
  )
  cat("After filtering:", ncol(combined_obj), "cells\n")
  # Save post-filtering metadata / 保存过滤后的元数据
  write.csv(combined_obj@meta.data, file.path(output_dir, "source_data", "post_filter_metadata.csv"), row.names = TRUE)
  
  # Normalization and feature selection / 标准化和特征选择
  cat("Performing normalization and feature selection...\n")
  combined_obj <- RunTFIDF(combined_obj)
  combined_obj <- FindTopFeatures(combined_obj, min.cutoff = 50)
  # Dimensional reduction
  cat("Performing dimensional reduction...\n")
  combined_obj <- RunSVD(combined_obj)
  # Integration: https://stuartlab.org/signac/articles/integrate_atac
  cat("Running ATAC integration...\n")
  # 分割对象
  obj_list <- SplitObject(combined_obj, split.by = "condition")
  # find integration anchors
  integration.anchors <- FindIntegrationAnchors(
    object.list = obj_list,
    anchor.features = rownames(combined_obj),
    reduction = "rlsi",
    dims = 2:30
  )
  # integrate LSI embeddings
  combined_obj <- IntegrateEmbeddings(
    anchorset = integration.anchors,
    reductions = combined_obj[["lsi"]],
    new.reduction.name = "integrated_lsi",
    dims.to.integrate = 1:30
  )
  
  # Clustering and UMAP / 聚类和UMAP
  cat("Performing clustering and UMAP...\n")
  combined_obj <- FindNeighbors(
    combined_obj, 
    reduction = "integrated_lsi", 
    dims = dims_neighbors
  )
  combined_obj <- FindClusters(
    combined_obj, 
    resolution = resolution,
    algorithm = 1  # Louvain algorithm / Louvain算法
  )
  combined_obj <- RunUMAP(
    combined_obj, 
    reduction = "integrated_lsi", 
    dims = dims_neighbors
  )
  # Calculate gene activity scores / 计算基因活性得分
  cat("Calculating gene activity scores...\n")
  gene_activities <- GeneActivity(combined_obj)
  combined_obj[["RNA"]] <- CreateAssayObject(counts = gene_activities)
  combined_obj <- NormalizeData(
    combined_obj, 
    assay = "RNA", 
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
  # Save final integrated object / 保存最终整合对象
  cat("Saving integrated object...\n")
  saveRDS(combined_obj, file.path(output_dir, "integration_results", "integrated_snatac_object.rds"))
  
  # Create integration assessment / 创建整合评估
  # ============================================================================
  # Integration Analysis Functions
  # 整合分析函数
  # ============================================================================
  
  #' Create integration quality assessment plots
  #' 创建整合质量评估图
  #' @param object Seurat object / Seurat对象
  #' @param group.by Grouping variable / 分组变量
  #' @param reduction Reduction to use / 使用的降维
  #'   Options: "umap", "tsne", "pca", "lsi", "harmony"
  #' @param output_dir Output directory / 输出目录
  create_integration_assessment <- function(object, 
                                            group.by = "condition",
                                            reduction = "umap",
                                            output_dir = "output") {
    cat("Creating integration quality assessment...\n")
    cat("创建整合质量评估...\n")
    # 1. UMAP plot colored by condition / 按条件着色的UMAP图
    p1 <- DimPlot(
      object = object,
      reduction = reduction,
      group.by = group.by,
      pt.size = 0.1
    ) + 
      ggtitle("Integration by Condition\n按条件显示的整合") +
      theme(plot.title = element_text(hjust = 0.5))
    # 2. UMAP plot colored by clusters / 按聚类着色的UMAP图
    p2 <- DimPlot(
      object = object,
      reduction = reduction,
      group.by = "seurat_clusters",
      label = TRUE,
      pt.size = 0.1
    ) + 
      ggtitle("Integration by Clusters\n按聚类显示的整合") +
      theme(plot.title = element_text(hjust = 0.5))
    # 3. Split plot by condition / 按条件分割图
    p3 <- DimPlot(
      object = object,
      reduction = reduction,
      group.by = group.by,
      split.by = group.by,
      pt.size = 0.1
    ) + 
      ggtitle("Integration Split by Condition\n按条件分割的整合") +
      theme(plot.title = element_text(hjust = 0.5))
    # 4. Feature plots for key QC metrics / 关键质量控制指标的特征图
    p4 <- FeaturePlot(
      object = object,
      features = c("nCount_peaks", "nFeature_peaks", "nucleosome_signal", "TSS.enrichment"),
      reduction = reduction,
      ncol = 2,
      pt.size = 0.1
    ) + 
      plot_annotation(title = "QC Metrics on Integration\n整合上的质量控制指标")
    # Combine integration plots
    combined_integration <- (p1 | p2) / p3 / p4
    # Save integration assessment
    ggsave(file.path(output_dir, "plots", "integration_assessment.pdf"),
           combined_integration, width = 12, height = 16, dpi = 300)
    
    cat("Integration quality assessment completed!\n")
    cat("整合质量评估完成！\n")
    return(combined_integration)
  }
  cat("Creating integration assessment...\n")
  integration_assessment <- create_integration_assessment(
    combined_obj,
    group.by = "condition",
    reduction = "umap",
    output_dir = output_dir
  )
  # Save final metadata / 保存最终元数据
  write.csv(combined_obj@meta.data, file.path(output_dir, "source_data", "final_integrated_metadata.csv"), row.names = TRUE)
  # Create summary report
  integration_summary <- data.frame(
    Metric = c(
      "Total Samples",
      "Total Cells (final)",
      "Total Peaks",
      "Number of Clusters",
      "Conditions",
      "Harmony Theta",
      "Clustering Resolution"
    ),
    Value = c(
      length(samples),
      ncol(combined_obj),
      nrow(combined_obj),
      length(unique(combined_obj$seurat_clusters)),
      paste(unique(combined_obj$condition), collapse = ", "),
      harmony_theta,
      resolution
    )
  )
  
  write.csv(integration_summary, file.path(output_dir, "source_data", "integration_summary.csv"), row.names = FALSE)
  
  cat("Advanced snATAC-seq integration completed successfully!\n")
####################################################################################
  # ============================================================================
  # Advanced peak annotation with detailed genomic features
  # 使用详细基因组特征进行高级peak注释
  # ============================================================================
  
  #' @param object Seurat object / Seurat对象
  #' @param txdb Transcript database / 转录本数据库
  #'   Options: TxDb.Hsapiens.UCSC.hg38.knownGene, TxDb.Mmusculus.UCSC.mm10.knownGene
  #' @param annoDb Annotation database / 注释数据库
  #'   Options: "org.Hs.eg.db" (human), "org.Mm.eg.db" (mouse)
  #' @param tssRegion TSS region definition / TSS区域定义
  #'   Options: c(-3000, 3000), c(-1000, 1000), c(-5000, 5000)
  #' @param output_dir Output directory / 输出目录
  # Peak annotation
  
    output_dir = output_dir
    object= combined_obj
    txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
    annoDb = "org.Hs.eg.db"
    tssRegion = c(-3000, 3000)
  
    cat("Starting advanced peak annotation...\n")
    # Get peak ranges / 获取peak范围
    peaks_gr <- granges(object)
    # Annotate peaks with ChIPseeker / 使用ChIPseeker注释peaks
    cat("Annotating peaks with ChIPseeker...\n")
    peakAnno <- annotatePeak(
      peaks_gr,
      tssRegion = tssRegion,
      TxDb = txdb,
      annoDb = annoDb,
      verbose = FALSE
    )
    # Store entire annotation data frame in misc
    object@misc$peak_annotation <- as.data.frame(peakAnno)
    # Add annotation to metadata
    peak_df <- as.data.frame(peakAnno)
    list_cols <- sapply(peak_df, is.list)
    peak_df[list_cols] <- lapply(peak_df[list_cols], function(x) sapply(x, paste, collapse=";"))
    # Save annotation data / 保存注释数据
    write.csv(peak_df, file.path(output_dir, "source_data", "peak_annotations.csv"), row.names = FALSE)
    # Create visualization plots
    # 1. Genomic annotation pie chart / 基因组注释饼图
    cat("Creating genomic annotation pie chart...\n")
    # 提取注释统计数据
    anno_df <- as.data.frame(peakAnno@annoStat)
    # 用 ggplot2 绘制饼图
    p1 <- ggplot(anno_df, aes(x = "", y = Frequency, fill = Feature)) +
      geom_bar(width = 1, stat = "identity") +
      coord_polar("y", start = 0) +
      labs(title = "Peak Genomic Annotation Distribution") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5))
    # 2. Distance to TSS distribution / 到TSS距离分布
    cat("Creating distance to TSS distribution plot...\n")
    p2 <- plotDistToTSS(peakAnno) + 
      ggtitle("Distance to TSS Distribution\n") +
      theme(plot.title = element_text(hjust = 0.5))
    # 3. Annotation bar plot / 注释柱状图
    cat("Creating annotation bar plot...\n")
    p3 <- plotAnnoBar(peakAnno) + 
      ggtitle("Peak Annotation Categories\nPeak注释类别") +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1))
    # Combine plots
    combined_plot <- (p1 | p2) / p3
    # Save plots and source data
    ggsave(file.path(output_dir, "plots", "peak_annotation_summary.pdf"),
           combined_plot, width = 12, height = 8, dpi = 300)
    # Save annotation summary statistics
    anno_summary <- as.data.frame(peakAnno@anno)
    anno_counts <- table(anno_summary$annotation)
    write.csv(as.data.frame(anno_counts), file.path(output_dir, "source_data", "annotation_summary_counts.csv"), row.names = FALSE)
    cat("Peak annotation completed successfully!\n")
#####################  ################################
