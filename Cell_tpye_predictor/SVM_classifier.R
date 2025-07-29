#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-07-29
# version:   1.0
# license:   MIT
# brief:     Model training based on Support Vector Machine (SVM).
#------------------------------------------#

## ====================================================
##  Cell Type Classifier (SVM)
## ====================================================

library(e1071)
library(doMC)
library(caret)  # 用于混淆矩阵和评估指标
library(MLmetrics)  # 用于对数损失计算

input_dir = 'D:/Data/0_yangfan/0_SingleCell_Scripts_Paper_2024/celltype_classifier_output'
output_dir <- 'D:/Data/0_yangfan/0_SingleCell_Scripts_Paper_2024/celltype_classifier_output/svm_classifier_output'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_data = file.path(input_dir,'datasets.rds')
datasets = readRDS(input_data)

output_path = output_dir
algorithm = 'Support_Vector_Machine'
model_name = "SVM_classifier"

## SVM 参数
ncores = 8
kernel_type = "radial"  # 核函数类型: "linear", "polynomial", "radial", "sigmoid"
cost_value = 1  # 正则化参数
gamma_value = 1 / ncol(datasets$train$x)  # 径向基核函数参数，设为1/特征数
cross_validation = 10  # 交叉验证折数

## 1. Model Training -----------------------------------------------------------
registerDoMC(cores = ncores)
start_time <- Sys.time() 

# 数据预处理：确保标签为因子类型
train_labels <- as.factor(datasets$train$y)
test_labels <- as.factor(datasets$test$y)

# 数据清理：处理NA/NaN/Inf值
cat("检查训练数据...\n")
cat(sprintf("训练数据维度: %d x %d\n", nrow(datasets$train$x), ncol(datasets$train$x)))
cat(sprintf("NA值数量: %d\n", sum(is.na(datasets$train$x))))
cat(sprintf("Inf值数量: %d\n", sum(is.infinite(datasets$train$x))))

# 移除或替换异常值
if(sum(is.na(datasets$train$x)) > 0 || sum(is.infinite(datasets$train$x)) > 0) {
  cat("发现异常值，进行处理...\n")
  # 用列均值替换NA值
  train_x_clean <- datasets$train$x
  test_x_clean <- datasets$test$x
  
  for(i in 1:ncol(train_x_clean)) {
    # 处理训练集
    col_mean <- mean(train_x_clean[,i], na.rm = TRUE)
    train_x_clean[is.na(train_x_clean[,i]) | is.infinite(train_x_clean[,i]), i] <- col_mean
    
    # 处理测试集
    test_x_clean[is.na(test_x_clean[,i]) | is.infinite(test_x_clean[,i]), i] <- col_mean
  }
  
  datasets$train$x <- train_x_clean
  datasets$test$x <- test_x_clean
  
  cat("异常值处理完成\n")
}

# 重新检查gamma值
if(gamma_value <= 0 || is.na(gamma_value) || is.infinite(gamma_value)) {
  gamma_value <- 1 / ncol(datasets$train$x)
  cat(sprintf("重置gamma值为: %f\n", gamma_value))
}

# 训练SVM模型
fit <- svm(
  x = datasets$train$x,
  y = train_labels,
  kernel = kernel_type,
  cost = cost_value,
  gamma = gamma_value,
  probability = TRUE,  # 启用概率预测
  cross = cross_validation,  # 交叉验证
  scale = TRUE  # 数据标准化
)

end_time <- Sys.time()
training_duration <- as.numeric(difftime(end_time, start_time, units = "secs")) # Time record
print(training_duration)
## Wrap model with metadata
model_obj <- list(
  model = fit,
  algorithm = algorithm,
  name = model_name,
  training_time = Sys.time(),
  parameters = list(
    kernel = kernel_type,
    cost = cost_value,
    gamma = gamma_value,
    cross_validation = cross_validation,
    scale = TRUE,
    n_features = ncol(datasets$train$x),  # 特征数量
    n_samples = nrow(datasets$train$x),   # 训练样本数
    n_classes = length(levels(train_labels)),  # 类别数量
    n_support_vectors = fit$tot.nSV  # 支持向量总数
  )
)

param_str <- paste0("kernel", kernel_type, "_cost", cost_value, "_gamma", 
                    round(gamma_value, 6))
filename <- paste0(model_name, "_", param_str, ".rds")
file_path <- file.path(output_path, filename)
saveRDS(model_obj, file = file_path)
cat(sprintf("模型已保存: %s (训练耗时: %.2f秒)\n", filename, training_duration))

## 2. Model Evaluation ---------------------------------------------------------
model <- model_obj
fit <- model_obj$model

# 预测
pred_class <- as.character(predict(fit, datasets$test$x))
pred_prob <- predict(fit, datasets$test$x, probability = TRUE)
pred_prob_matrix <- attr(pred_prob, "probabilities")

# 确保预测结果和真实标签格式一致
true_labels <- test_labels
pred_type <- as.factor(pred_class)

# 确保预测类别的levels与真实标签一致
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

## 2.4 SVM特有指标
n_support_vectors <- fit$tot.nSV  # 支持向量总数
sv_per_class <- fit$nSV  # 每个类别的支持向量数量
cv_accuracy <- fit$tot.accuracy / 100  # 交叉验证准确率（转换为0-1范围）

## 2.5 对数损失 (Log Loss)
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

## 2.6 详细的混淆矩阵指标
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
} else {
  ## 二分类情况
  precision_per_class <- cm_detailed$byClass["Pos Pred Value"]
  recall_per_class <- cm_detailed$byClass["Sensitivity"]
  f1_per_class <- cm_detailed$byClass["F1"]
  specificity_per_class <- cm_detailed$byClass["Specificity"]
  
  macro_precision <- precision_per_class
  macro_recall <- recall_per_class
  macro_f1 <- f1_per_class
  macro_specificity <- specificity_per_class
}

## 收集评估结果
evaluation_results <- list(
  kernel = kernel_type,  # 核函数类型
  cost = cost_value,  # 正则化参数
  gamma = gamma_value,  # 核函数参数
  accuracy = accuracy, # 总体准确率，表示正确预测的样本比例
  normalized_accuracy = normalized_accuracy,  # 归一化准确率，参考paper计算
  balanced_accuracy = balanced_accuracy,  # 平衡准确率，是每个类别召回率的平均值
  macro_precision = macro_precision,  # 宏平均精确率（Macro Precision），是每个类别的精确率（Precision）的平均
  macro_recall = macro_recall,  # 宏平均召回率（Macro Recall），是每个类别召回率的平均值
  macro_f1 = macro_f1,  # 宏平均 F1 分数（Macro F1），是每个类别 F1 分数的平均值，综合精确率和召回率
  macro_specificity = macro_specificity, # 宏平均特异性
  log_loss = log_loss,  # 对数损失（Log Loss），衡量预测概率与真实标签之间的差异
  n_support_vectors = n_support_vectors,  # 支持向量总数
  cv_accuracy = cv_accuracy,  # 交叉验证准确率
  training_duration = training_duration  # 训练时间
)

# 保存评估结果
evaluation_file <- file.path(output_path, "evaluation_results.rds")
saveRDS(evaluation_results, file = evaluation_file)

## 创建评估结果汇总表
summary_df <- data.frame(
  Metric = c("Accuracy", "Normalized_Accuracy", "Balanced_Accuracy", 
             "Macro_Precision", "Macro_Recall", "Macro_F1", "Macro_Specificity",
             "Log_Loss", "N_Support_Vectors", "CV_Accuracy", 
             "Training_Duration_Sec"),
  Value = c(evaluation_results$accuracy, evaluation_results$normalized_accuracy,
            evaluation_results$balanced_accuracy, evaluation_results$macro_precision,
            evaluation_results$macro_recall, evaluation_results$macro_f1,
            evaluation_results$macro_specificity, evaluation_results$log_loss,
            evaluation_results$n_support_vectors, evaluation_results$cv_accuracy,
            evaluation_results$training_duration)
)

## 保存汇总表
summary_file <- file.path(output_path, "model_summary.csv")
write.csv(summary_df, file = summary_file, row.names = FALSE)

## 保存支持向量信息
# 确保sv_per_class有正确的类别名称
if(is.null(names(sv_per_class))) {
  names(sv_per_class) <- levels(train_labels)
}

sv_info_df <- data.frame(
  Class = names(sv_per_class),
  N_Support_Vectors = as.numeric(sv_per_class),
  Percentage = round(sv_per_class / sum(sv_per_class) * 100, 2)
)
sv_info_file <- file.path(output_path, "support_vectors_info.csv")
write.csv(sv_info_df, file = sv_info_file, row.names = FALSE)

## 保存每个类别的详细指标
class_metrics_df <- data.frame(
  Class = levels(true_labels),
  Precision = precision_per_class,
  Recall = recall_per_class,
  F1_Score = f1_per_class,
  Specificity = specificity_per_class
)
class_metrics_file <- file.path(output_path, "class_metrics.csv")
write.csv(class_metrics_df, file = class_metrics_file, row.names = FALSE)

## 输出结果汇总
cat("\n=== SVM 分类器评估结果 ===\n")
cat(sprintf("总体准确率: %.4f\n", evaluation_results$accuracy))
cat(sprintf("归一化准确率: %.4f\n", evaluation_results$normalized_accuracy))
cat(sprintf("平衡准确率: %.4f\n", evaluation_results$balanced_accuracy))
cat(sprintf("宏平均F1分数: %.4f\n", evaluation_results$macro_f1))
cat(sprintf("对数损失: %.4f\n", evaluation_results$log_loss))
cat(sprintf("交叉验证准确率: %.4f\n", evaluation_results$cv_accuracy))
cat(sprintf("支持向量总数: %d\n", evaluation_results$n_support_vectors))
cat(sprintf("训练时间: %.2f秒\n", evaluation_results$training_duration))
cat(sprintf("模型参数: kernel=%s, cost=%.2f, gamma=%.6f\n", 
            evaluation_results$kernel, evaluation_results$cost, evaluation_results$gamma))

print(summary_df)

cat("\n=== 各类别支持向量分布 ===\n")
print(sv_info_df)
