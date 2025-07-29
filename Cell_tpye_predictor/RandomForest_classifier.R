#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-07-29
# version:   1.0
# license:   MIT
# brief:     Model training based on Random Forest (RF).
#------------------------------------------#

## ====================================================
##  Cell Type Classifier (randomforest)
## ====================================================

library(randomForest)
library(doMC)
library(caret)  # 用于混淆矩阵和评估指标
library(MLmetrics)  # 用于对数损失计算

input_dir = 'D:/Data/0_yangfan/0_SingleCell_Scripts_Paper_2024/celltype_classifier_output'
output_dir <- 'D:/Data/0_yangfan/0_SingleCell_Scripts_Paper_2024/celltype_classifier_output/randomforest_classifier_output'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

intput_data = file.path(input_dir,'datasets.rds')
datasets = readRDS(intput_data)

output_path = output_dir
algorithm = 'Random_Forest'
model_name = "RandomForest_classifier"
## Random Forest 
ncores = 8
ntree = 500

## 1. Model Training -----------------------------------------------------------
registerDoMC(cores = ncores)
start_time <- Sys.time() 
fit <- randomForest(
  x = datasets$train$x,
  y = as.factor(datasets$train$y),
  ntree = ntree,
  importance = TRUE
)
end_time <- Sys.time()
training_duration <- as.numeric(difftime(end_time, start_time, units = "secs")) # Time record
  
## Wrap model with metadata
model_obj <- list(
  model = fit,
  algorithm = algorithm,
  name = model_name,
  training_time = Sys.time(),
  parameters = list(
    ntree = ntree,
    mtry = fit$mtry,  # 每次分割时考虑的变量数
    n_features = ncol(datasets$train$x),  # 特征数量
    n_samples = nrow(datasets$train$x),   # 训练样本数
    n_classes = length(levels(datasets$train$y))  # 类别数量
  )
)

param_str <- paste0("rf_nt",ntree)
filename <- paste0(model_name, "_", param_str, ".rds")
file_path=file.path(output_path, filename)
saveRDS(model_obj, file = file_path)
cat(sprintf("模型已保存: %s (训练耗时: %.2f秒)\n", filename, training_duration))

## 2. Model Evaluation ---------------------------------------------------------
model = model_obj
fit <- model_obj$model

pred_class <- as.character(predict(fit, datasets$test$x))
pred_prob <- predict(fit, datasets$test$x, type = "prob")
importance_scores <- importance(fit)

## 确保预测结果和真实标签格式一致
true_labels <- as.factor(datasets$test$y)
pred_type <- as.factor(pred_class)
## 确保预测类别的levels与真实标签一致
levels(pred_type) <- levels(true_labels)
## 2.1 总体准确率 (Overall Accuracy)
accuracy <- mean(pred_type == true_labels)
## 2.2 标准化准确率 (Normalized Accuracy) - 每类的准确率的平均
class_accuracies <- sapply(levels(true_labels), function(class) {
  class_indices <- which(true_labels == class)
  if(length(class_indices) > 0) {
    mean(pred_type[class_indices] == true_labels[class_indices])
  } else {
    0
  }
})
normalized_accuracy <- round(mean(class_accuracies, na.rm = TRUE), 4)
## 2.3 平衡准确率 (Balanced Accuracy) | 各类别召回率的平均值
cm <- confusionMatrix(pred_type, true_labels)  # 生成混淆矩阵，计算分类性能指标
balanced_accuracy <- cm$byClass[,"Balanced Accuracy"]  # 提取每个类别的平衡准确率
if (is.vector(balanced_accuracy) && length(balanced_accuracy) > 1) {
  balanced_accuracy <- mean(balanced_accuracy, na.rm = TRUE)
}
## 2.4 特征重要性统计 (Feature Importance Statistics) | 统计重要特征数量（基于Mean Decrease Accuracy > 阈值）
importance_threshold <- quantile(importance_scores[,"MeanDecreaseAccuracy"], 0.75)  # 75分位数作为阈值
n_important_features <- sum(importance_scores[,"MeanDecreaseAccuracy"] > importance_threshold)
## 2.5 对数损失 (Log Loss) | 将概率矩阵转换为向量形式 | 三维数组情况：转换为矩阵
if (length(dim(pred_prob)) == 3) { pred_prob_matrix <- pred_prob[,,1]
  } else {  pred_prob_matrix <- pred_prob  }
## 确保概率矩阵列名与类别标签一致
if (!is.null(colnames(pred_prob_matrix))) {
  pred_prob_matrix <- pred_prob_matrix[, levels(true_labels), drop = FALSE] # 重新排序概率矩阵以匹配真实标签的level顺序
}
## 计算对数损失
log_loss <- tryCatch({
  MLmetrics::MultiLogLoss(y_true = as.numeric(true_labels), y_pred = pred_prob_matrix)
}, error = function(e) {
  ## 如果计算失败，使用手动计算
  eps <- 1e-15  # 避免log(0)
  pred_prob_clipped <- pmax(pmin(pred_prob_matrix, 1-eps), eps)
  true_labels_numeric <- as.numeric(true_labels)
  -mean(log(pred_prob_clipped[cbind(seq_along(true_labels_numeric), true_labels_numeric)]))
})
## 2.6 Random Forest 特有指标
oob_error <- fit$err.rate[ntree, "OOB"]  # Out-of-Bag错误率
mean_decrease_accuracy <- mean(importance_scores[,"MeanDecreaseAccuracy"])
mean_decrease_gini <- mean(importance_scores[,"MeanDecreaseGini"])
## 2.7 详细的混淆矩阵指标
cm_detailed <- confusionMatrix(pred_type, true_labels)  # 为每个类别计算精确率、召回率、F1分数
## 提取每个类别的指标
if (length(levels(true_labels)) > 2) {
  ## 多分类情况
  precision_per_class <- cm_detailed$byClass[,"Pos Pred Value"]
  recall_per_class <- cm_detailed$byClass[,"Sensitivity"]
  f1_per_class <- cm_detailed$byClass[,"F1"]
  specificity_per_class <- cm_detailed$byClass[,"Specificity"]
  ## 宏平均指标
  macro_precision <- mean(precision_per_class, na.rm = TRUE)
  macro_recall <- mean(recall_per_class, na.rm = TRUE)
  macro_f1 <- mean(f1_per_class, na.rm = TRUE)
  macro_specificity <- mean(specificity_per_class, na.rm = TRUE)
} else {  stop("二分类情况, 请检查")}
  
## 收集评估结果
evaluation_results <- list(
  ntree = ntree,  # 树的数量
  mtry = fit$mtry,  # 每次分割考虑的变量数
  accuracy = accuracy, # 总体准确率，表示正确预测的样本比例
  normalized_accuracy = normalized_accuracy,  # 归一化准确率，参考paper计算
  balanced_accuracy = balanced_accuracy,  # 平衡准确率，是每个类别召回率的平均值
  macro_precision = macro_precision,  # 宏平均精确率（Macro Precision），是每个类别的精确率（Precision）的平均
  macro_recall = macro_recall,  # 宏平均召回率（Macro Recall），是每个类别召回率的平均值
  macro_f1 = macro_f1,  # 宏平均 F1 分数（Macro F1），是每个类别 F1 分数的平均值，综合精确率和召回率
  macro_specificity = macro_specificity, # 宏平均特异性
  log_loss = log_loss,  # 对数损失（Log Loss），衡量预测概率与真实标签之间的差异
  n_important_features = n_important_features,  # 重要特征数量
  oob_error = oob_error,  # Out-of-Bag错误率
  mean_decrease_accuracy = mean_decrease_accuracy,  # 平均准确度下降
  mean_decrease_gini = mean_decrease_gini,  # 平均基尼不纯度下降
  training_duration = training_duration  # 训练时间
)

# 保存评估结果
evaluation_file <- file.path(output_path, "RF_evaluation_results.rds")
saveRDS(evaluation_results, file = evaluation_file)

## 创建评估结果汇总表
summary_df <- data.frame(
  Metric = c("Accuracy", "Normalized_Accuracy", "Balanced_Accuracy", 
             "Macro_Precision", "Macro_Recall", "Macro_F1", "Macro_Specificity",
             "Log_Loss", "N_Important_Features", "OOB_Error", 
             "Mean_Decrease_Accuracy", "Mean_Decrease_Gini", "Training_Duration_Sec"),
  Value = c(evaluation_results$accuracy, evaluation_results$normalized_accuracy,
            evaluation_results$balanced_accuracy, evaluation_results$macro_precision,
            evaluation_results$macro_recall, evaluation_results$macro_f1,
            evaluation_results$macro_specificity, evaluation_results$log_loss,
            evaluation_results$n_important_features, evaluation_results$oob_error,
            evaluation_results$mean_decrease_accuracy, evaluation_results$mean_decrease_gini,
            evaluation_results$training_duration)
)

## 保存汇总表
summary_file <- file.path(output_path, "RF_model_summary.csv")
write.csv(summary_df, file = summary_file, row.names = FALSE)

## 保存特征重要性
importance_df <- data.frame(
  Feature = rownames(importance_scores),
  MeanDecreaseAccuracy = importance_scores[,"MeanDecreaseAccuracy"],
  MeanDecreaseGini = importance_scores[,"MeanDecreaseGini"]
)
importance_file <- file.path(output_path, "RF_feature_importance.csv")
write.csv(importance_df, file = importance_file, row.names = FALSE)

## 输出结果汇总
cat("\n=== Random Forest 分类器评估结果 ===\n")
cat(sprintf("总体准确率: %.4f\n", evaluation_results$accuracy))
cat(sprintf("归一化准确率: %.4f\n", evaluation_results$normalized_accuracy))
cat(sprintf("平衡准确率: %.4f\n", evaluation_results$balanced_accuracy))
cat(sprintf("宏平均F1分数: %.4f\n", evaluation_results$macro_f1))
cat(sprintf("对数损失: %.4f\n", evaluation_results$log_loss))
cat(sprintf("OOB错误率: %.4f\n", evaluation_results$oob_error))
cat(sprintf("重要特征数量: %d\n", evaluation_results$n_important_features))
cat(sprintf("训练时间: %.2f秒\n", evaluation_results$training_duration))
cat(sprintf("模型参数: ntree=%d, mtry=%d\n", evaluation_results$ntree, evaluation_results$mtry))

print(summary_df)
