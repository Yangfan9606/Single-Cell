#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-07-29
# version:   1.0
# license:   MIT
# brief:     Cell type prediction use output of SVM_classifier.R.
#------------------------------------------#

## ====================================================
##  Cell Type Predictor (Support Vector Machine)
## ====================================================

# 加载必需的库
library(dplyr)    # 数据处理
library(Seurat)   # 单细胞RNA分析
library(e1071)    # SVM算法实现

# 设置输入和输出路径
input_dir <- ('D:/Data/0_yangfan/0_SingleCell_RNA/scRNA_analysis_results/data')
my_obj = readRDS(paste0(input_dir,'/seurat_clustered_object.rds'))
output_dir <- 'D:/Data/0_yangfan/0_SingleCell_Scripts_Paper_2024/celltype_classifier_output/svm_classifier_output'

seurat_obj = my_obj

## 加载训练好的SVM分类器
# SVM模型文件通常包含模型对象和相关的预处理信息
classifier_path = file.path(output_dir, "SVM_classifier_kernelradial_cost1_gamma0.001429.rds")
model_data <- readRDS(classifier_path)

# 提取SVM模型和所需基因列表
best_model <- model_data$model          # SVM模型对象
required_genes <- colnames(best_model$SV)   # 训练时使用的基因特征
# 如果模型文件结构不同，可能需要调整：
# required_genes <- attr(best_model$terms, "term.labels")  # 备选方法

## 关键参数设置
confidence_threshold = 0.7    # 置信度阈值：SVM通常需要较高的阈值（0.8）
# 因为SVM的决策边界相对严格
entropy_threshold = 0.6      # 熵阈值：SVM的概率输出相对稳定，使用中等阈值
distance_threshold = 0.3     # 决策边界距离阈值：SVM特有参数，衡量样本到决策边界的距离

# 输出列名定义
output_col = "predicted_celltype"        # 预测细胞类型列名
confidence_col = "prediction_confidence" # 预测置信度列名
entropy_col = "prediction_entropy"       # 预测熵列名
distance_col = "decision_distance"       # SVM特有：决策距离列名

## 准备特征矩阵
# 找到模型所需基因与当前数据集的交集
feature.space <- intersect(required_genes, rownames(seurat_obj))
missing.features <- setdiff(required_genes, rownames(seurat_obj))

if(length(missing.features) > 0) {
  message("Warning: ", length(missing.features), " required genes are missing from the dataset")
  message("Missing genes: ", paste(head(missing.features, 10), collapse = ", "))
}

## 数据预处理 - 与训练时保持严格一致
# SVM对数据缩放非常敏感，必须确保预处理步骤完全一致
DefaultAssay(seurat_obj) = "RNA"

# 标准化和缩放数据
temp <- NormalizeData(seurat_obj, verbose = F) %>% 
  ScaleData(features = feature.space, verbose = F)

## 提取缩放后的数据矩阵
# SVM需要样本为行，特征为列的数据格式
x_new <- t(GetAssayData(temp, assay = 'RNA', layer = "scale.data"))
rm(temp)  # 释放内存

## 处理基因特征匹配问题
missing_genes <- setdiff(required_genes, colnames(x_new))
extra_genes <- setdiff(colnames(x_new), required_genes)

## 处理缺失基因 - 用0填充
# SVM对缺失特征敏感，需要谨慎处理
if (length(missing_genes) > 0) {
  message("Filling ", length(missing_genes), " missing genes with zeros")
  message("Note: SVM performance may be affected by missing features")
  
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

## 使用SVM进行预测
message("Making predictions with Support Vector Machine...")

# 转换为数据框格式（e1071包要求）
x_new_df <- as.data.frame(x_new)

# SVM类别预测
pred_class_raw <- predict(best_model, newdata = x_new_df, type = "class")

# SVM概率预测（需要模型训练时设置probability=TRUE）
pred_prob_raw <- predict(best_model, newdata = x_new_df, probability = TRUE)
pred_prob_matrix <- attr(pred_prob_raw, "probabilities")

# 获取决策距离（SVM特有的置信度指标）
decision_values <- predict(best_model, newdata = x_new_df, decision.values = TRUE)
decision_distances <- attr(decision_values, "decision.values")

# 对于多类问题，计算到最近决策边界的距离
if(ncol(decision_distances) > 1) {
  # 多类SVM：使用最大决策值作为距离度量
  min_decision_distance <- apply(abs(decision_distances), 1, min)
} else {
  # 二类SVM：直接使用决策值的绝对值
  min_decision_distance <- abs(decision_distances[, 1])
}

## 计算预测指标
# 1. 最大概率作为置信度
max_prob <- apply(pred_prob_matrix, 1, max)

# 2. 计算Shannon熵来衡量预测不确定性
entropy <- apply(pred_prob_matrix, 1, function(p) {
  p <- p[p > 0]  # 移除概率为0的类别
  if(length(p) == 0) return(0)
  -sum(p * log2(p))  # Shannon熵公式
})

# 3. SVM特有：标准化决策距离作为额外的置信度指标
normalized_distance <- min_decision_distance / max(min_decision_distance)

## 应用质量控制阈值识别不确定预测
# 多重条件判断：置信度低 OR 熵值高 OR 距离决策边界太近
uncertain_pred <- (max_prob < confidence_threshold) | 
  (entropy > entropy_threshold) | 
  (normalized_distance < distance_threshold)

# 将不确定的预测标记为"Uncertain"
final_pred <- as.character(pred_class_raw)
final_pred[uncertain_pred] <- "Uncertain"

## 将预测结果添加到Seurat对象的元数据中
seurat_obj@meta.data[[output_col]] <- final_pred
seurat_obj@meta.data[[confidence_col]] <- max_prob
seurat_obj@meta.data[[entropy_col]] <- entropy
seurat_obj@meta.data[[distance_col]] <- normalized_distance

## 输出预测结果摘要
message("SVM predictions completed!")
message("Total cells: ", length(final_pred))
message("Uncertain predictions: ", sum(uncertain_pred), " (", 
        round(100 * sum(uncertain_pred) / length(uncertain_pred), 2), "%)")

# 显示预测类别分布
cat("\nPredicted cell type distribution:\n")
print(table(final_pred))

# 显示前几个细胞的预测结果
cat("\nFirst few predictions:\n")
print(head(seurat_obj@meta.data[, c(output_col, confidence_col, entropy_col, distance_col)]))

## 可视化预测结果
# 在UMAP/t-SNE图上显示预测的细胞类型
print(DimPlot(seurat_obj, group.by = "predicted_celltype", label = T, pt.size = 0.5) +
        ggtitle("SVM Cell Type Predictions"))

## 输出关键指标用于后续分析
cat("\nPrediction quality metrics:\n")
cat("Mean confidence:", round(mean(max_prob), 3), "\n")
cat("Mean entropy:", round(mean(entropy), 3), "\n")
cat("Mean decision distance:", round(mean(normalized_distance), 3), "\n")

# SVM特有的模型信息
cat("\nSVM model information:\n")
cat("Kernel type:", best_model$kernel, "\n")
cat("Number of support vectors:", best_model$tot.nSV, "\n")
cat("Number of classes:", length(best_model$levels), "\n")
