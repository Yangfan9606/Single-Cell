#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-07-29
# version:   1.0
# license:   MIT
# brief:     Cell type prediction use output of RandomForest_classifier.R.
#------------------------------------------#

## ====================================================
##  Cell Type Predictor (Random Forest)
## ====================================================

# 加载必需的库
library(dplyr)          # 数据处理
library(Seurat)         # 单细胞RNA分析
library(randomForest)   # 随机森林算法

# 设置输入和输出路径
input_dir <- ('D:/Data/0_yangfan/0_SingleCell_RNA/scRNA_analysis_results/data')
my_obj = readRDS(paste0(input_dir,'/seurat_clustered_object.rds'))
output_dir <- 'D:/Data/0_yangfan/0_SingleCell_Scripts_Paper_2024/celltype_classifier_output/randomforest_classifier_output'

seurat_obj = my_obj

## 加载训练好的随机森林分类器
# Random Forest模型通常保存为包含模型对象的RDS文件
classifier_path = file.path(output_dir, "RandomForest_classifier_rf_nt500.rds")
model_data <- readRDS(classifier_path)
best_model = model_data$model   # 模型对象

# 从训练好的模型中提取所需的基因特征
# Random Forest模型的特征变量存储在模型对象中
required_genes <- rownames(best_model$importance)  # 获取特征重要性排序的基因

## 关键参数设置
confidence_threshold = 0.7    # 置信度阈值：Random Forest通常使用稍低的阈值（0.7 vs 0.9）
# 因为RF输出的概率分布相对更平滑
entropy_threshold = 1.0      # 熵阈值：RF允许稍高的熵值，因为其预测更稳定
min_vote_ratio = 0.5        # 最小投票比例：RF特有参数，确保预测类别得到足够的树投票

# 输出列名定义
output_col = "predicted_celltype"        # 预测细胞类型列名
confidence_col = "prediction_confidence" # 预测置信度列名
entropy_col = "prediction_entropy"       # 预测熵列名
vote_col = "prediction_votes"           # RF特有：投票数列名

## 准备特征矩阵
# 找到模型所需基因与当前数据集的交集
feature.space <- intersect(required_genes, rownames(seurat_obj))
missing.features <- setdiff(required_genes, rownames(seurat_obj))

if(length(missing.features) > 0) {
  message("Warning: ", length(missing.features), " required genes are missing from the dataset")
  message("Missing genes: ", paste(head(missing.features, 10), collapse = ", "))
}

## 数据预处理 - 与训练时保持一致
# 设置默认检测方法为RNA
DefaultAssay(seurat_obj) = "RNA"

# 标准化和缩放数据，确保与训练时的预处理步骤一致
temp <- NormalizeData(seurat_obj, verbose = F) %>% 
  ScaleData(features = feature.space, verbose = F)

## 提取缩放后的数据矩阵
# Random Forest需要样本为行，特征为列的矩阵
x_new <- t(GetAssayData(temp, assay = 'RNA', layer = "scale.data"))
rm(temp)  # 释放内存

## 处理基因特征匹配问题
missing_genes <- setdiff(required_genes, colnames(x_new))
extra_genes <- setdiff(colnames(x_new), required_genes)

## 处理缺失基因 - 用0填充
# Random Forest对缺失特征相对鲁棒，但仍需保持特征维度一致
if (length(missing_genes) > 0) {
  message("Filling ", length(missing_genes), " missing genes with zeros")
  fill_matrix <- matrix(
    0,
    nrow = nrow(x_new),
    ncol = length(missing_genes),
    dimnames = list(rownames(x_new), missing_genes)
  )
  x_new <- cbind(x_new, fill_matrix)
}

## 按训练时的基因顺序重新排列特征
x_new <- x_new[, required_genes, drop = FALSE]

## 数据完整性检查
stopifnot(ncol(x_new) == length(required_genes))
stopifnot(all(colnames(x_new) == required_genes))

message("Feature matrix prepared: ", nrow(x_new), " cells × ", ncol(x_new), " genes")

## 使用Random Forest进行预测
message("Making predictions with Random Forest...")

# Random Forest预测：同时获取类别预测和概率预测
pred_class_raw <- predict(best_model, newdata = x_new, type = "class")
pred_prob_raw <- predict(best_model, newdata = x_new, type = "prob")

# 获取投票信息（Random Forest特有）
pred_votes <- predict(best_model, newdata = x_new, type = "vote")

## 计算预测指标
# 1. 最大概率作为置信度
max_prob <- apply(pred_prob_raw, 1, max)

# 2. 计算Shannon熵来衡量预测不确定性
entropy <- apply(pred_prob_raw, 1, function(p) {
  p <- p[p > 0]  # 移除概率为0的类别
  if(length(p) == 0) return(0)
  -sum(p * log2(p))  # Shannon熵公式
})

# 3. Random Forest特有：计算最高投票比例
max_votes <- apply(pred_votes, 1, max)
total_trees <- rowSums(pred_votes)
vote_ratio <- max_votes / total_trees

## 应用质量控制阈值识别不确定预测
# 多重条件判断：置信度低 OR 熵值高 OR 投票比例低
uncertain_pred <- (max_prob < confidence_threshold) | 
  (entropy > entropy_threshold) | 
  (vote_ratio < min_vote_ratio)

# 将不确定的预测标记为"Uncertain"
final_pred <- as.character(pred_class_raw)
final_pred[uncertain_pred] <- "Uncertain"

## 将预测结果添加到Seurat对象的元数据中
seurat_obj@meta.data[[output_col]] <- final_pred
seurat_obj@meta.data[[confidence_col]] <- max_prob
seurat_obj@meta.data[[entropy_col]] <- entropy
seurat_obj@meta.data[[vote_col]] <- vote_ratio

## 输出预测结果摘要
message("Random Forest predictions completed!")
message("Total cells: ", length(final_pred))
message("Uncertain predictions: ", sum(uncertain_pred), " (", 
        round(100 * sum(uncertain_pred) / length(uncertain_pred), 2), "%)")

# 显示预测类别分布
cat("\nPredicted cell type distribution:\n")
print(table(final_pred))

# 显示前几个细胞的预测结果
cat("\nFirst few predictions:\n")
print(head(seurat_obj@meta.data[, c(output_col, confidence_col, entropy_col, vote_col)]))

## 可视化预测结果
# 在UMAP/t-SNE图上显示预测的细胞类型
print(DimPlot(seurat_obj, group.by = "predicted_celltype", label = T, pt.size = 0.5) +
        ggtitle("Random Forest Cell Type Predictions"))

## 输出关键指标用于后续分析
cat("\nPrediction quality metrics:\n")
cat("Mean confidence:", round(mean(max_prob), 3), "\n")
cat("Mean entropy:", round(mean(entropy), 3), "\n")
cat("Mean vote ratio:", round(mean(vote_ratio), 3), "\n")
