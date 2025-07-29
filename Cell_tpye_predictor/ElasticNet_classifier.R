#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-07-29
# version:   1.0
# license:   MIT
# brief:     Model training based on Elastic Net (EN).
#------------------------------------------#

## ====================================================
##  Cell Type Classifier (elastic net)
## ====================================================
library(glmnet)
library(doMC)
library(caret)  # 用于混淆矩阵和评估指标
library(MLmetrics)  # 用于对数损失计算

input_dir = 'D:/Data/0_yangfan/0_SingleCell_Scripts_Paper_2024/celltype_classifier_output'
output_dir <- 'D:/Data/0_yangfan/0_SingleCell_Scripts_Paper_2024/celltype_classifier_output/elasticnet_classifier_output'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

intput_data = file.path(input_dir,'datasets.rds')
datasets = readRDS(intput_data)

output_path = output_dir
algorithm = 'Elastic_Net'
model_name = "ElasticNet_classifier"
## Elastic Net 
alphas = c (0, 0.25, 0.5, 0.75, 1)  # "Ridge"=0, "Lasso"=1
lambda_seq <- 10^seq(-8, 3, length.out = 200)
cv_folds = 10
cv_repeats = 1
ncores = 8

## 1. Model Training -----------------------------------------------------------
registerDoMC(cores = ncores)
for (i in seq_along(alphas)) {
  alpha <- alphas[i]
  cat(sprintf("训练 alpha=%.2f 的模型 (%d/%d)...\n", alpha, i, length(alphas)))
  start_time <- Sys.time() 
  fit <- cv.glmnet(
    x = datasets$train$x,  # x 训练集特征矩阵（细胞x基因）
    y = datasets$train$y,  # y 训练集标签（细胞类型因子向量）
    weights = datasets$train$w,  # y 训练集标签（细胞类型因子向量）
    family = "multinomial",  # family 多分类问题指定（二分类用"binomial"）
    type.measure = "deviance",  # type.measure 评估指标：偏差（=2*对数损失，概率校准敏感）
    alpha = alpha,  # alpha 弹性网混合参数（0=Ridge, 1=Lasso）
    lambda = lambda_seq,  # lambdas 创建从10^-8到10^3的对数均匀分布序列（覆盖典型范围）
    nfolds = cv_folds,  # nfolds 交叉验证折数（通常5-10）default: 10
    parallel = TRUE  # parallel 启用并行计算（加速交叉验证）
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
      alpha = alpha,
      cv_folds = cv_folds,
      family = "multinomial",
      type.measure = "deviance",
      lambda_length = length(lambda_seq),  # lambda序列长度
      n_features = ncol(datasets$train$x),  # 特征数量
      n_samples = nrow(datasets$train$x),   # 训练样本数
      n_classes = length(levels(datasets$train$y))  # 类别数量
    )
  )
  param_str <- paste0("a", format(alpha, nsmall=2), "_l", length(lambda_seq))
  filename <- paste0(model_name, "_", param_str, ".rds")
  file_path=file.path(output_path, filename)
  saveRDS(model_obj, file = file_path)
  cat(sprintf("模型已保存: %s (训练耗时: %.2f秒)\n", filename, training_duration))
}

## 2. Model Evaluation ---------------------------------------------------------
evaluation_results <- list()
for (i in seq_along(alphas)) {
  alpha <- alphas[i]
  cat(sprintf("评估 alpha=%.2f 的模型 (%d/%d)...\n", alpha, i, length(alphas)))
  ## 加载训练好的模型
  param_str <- paste0("a", format(alpha, nsmall=2), "_l", length(lambda_seq))
  filename <- paste0(model_name, "_", param_str, ".rds")
  file_path=file.path(output_path, filename)
  model = readRDS(file_path)
  fit <- model_obj$model
  ## 获取最优lambda对应的预测结果
  pred_prob <- predict(fit, newx = datasets$test$x, 
                       s = "lambda.1se",  # 误差在最小值一倍标准差内的最简模型, 相对不精确，但需要数据量小，且更适合防止过拟合 
                       type = "response")
  pred_class <- predict(fit, newx = datasets$test$x, s = "lambda.1se", type = "class")
  pred_type <- as.factor(pred_class[,1])  # 提取预测类别, glmnet 内部会自动对概率矩阵（type="response" 的结果）执行类似 which.max 的操作
####  pred.prob <- predict(model, test$x, s="lambda.min", type="response")
####  pred.type <- colnames(pred.prob)[apply(pred.prob, 1, which.max)]  # 手动提取 
  ## 获取真实标签
  true_labels <- as.factor(datasets$test$y)
  
  ## 2.1 总体准确率 (Overall Accuracy)
  accuracy <- mean(pred_type == true_labels)
  ## 2.2 标准化准确率 (Normalized Accuracy)  每类的准确率 的 平均
  normalized_accuracy=round(1 * sum((pred_type == true_labels)/table(true_labels)[true_labels])/length(unique(true_labels)), 4)
  ## 2.3 平衡准确率 (Balanced Accuracy) | 各类别召回率的平均值
  cm <- confusionMatrix(pred_type, true_labels)  # 生成混淆矩阵，计算分类性能指标
  balanced_accuracy <- cm$byClass[,"Balanced Accuracy"]  # 提取每个类别的平衡准确率
  if (is.vector(balanced_accuracy) && length(balanced_accuracy) > 1) {
    balanced_accuracy <- mean(balanced_accuracy, na.rm = TRUE)  }
  ## 2.4 非零系数数量 (Number of Non-zero Coefficients)
  coef_matrix <- coef(fit, s = "lambda.1se")  # 获取最优lambda下的系数矩阵
  if (is.list(coef_matrix)) {
    ## 多分类情况：统计所有类别的非零系数
    n_coef <- sum(sapply(coef_matrix, function(x) sum(x != 0))) - length(coef_matrix)  # 减去截距项
  } else {    n_coef <- sum(coef_matrix != 0) - 1  # 二分类情况 | 减去截距项
  }
  ## 2.5 对数损失 (Log Loss)
    ## 将概率矩阵转换为向量形式
  if (length(dim(pred_prob)) == 3) {
    ## 三维数组情况：转换为矩阵
    pred_prob_matrix <- pred_prob[,,1]
  } else {    pred_prob_matrix <- pred_prob  }
    ## 确保概率矩阵列名与类别标签一致
  if (!is.null(colnames(pred_prob_matrix))) {pred_prob_matrix <- pred_prob_matrix[, levels(true_labels), drop = FALSE] # 重新排序概率矩阵以匹配真实标签的level顺序
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
  ## 2.6 交叉验证偏差 (Cross-validation Deviance)
  deviance_min <- min(fit$cvm)           # 最小交叉验证偏差
  deviance_1se <- fit$cvm[fit$lambda == fit$lambda.1se]  # lambda.1se对应的偏差
  lambda_min <- fit$lambda.min
  lambda_1se <- fit$lambda.1se
  
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
  } else {    stop("二分类情况, 请检查")  }
  
  ## 收集当前模型结果
  evaluation_results[[i]] <- list(
    alpha = alpha,  # Elastic Net 模型中控制 L1（Lasso）和 L2（Ridge）正则化权重的混合参数
    accuracy = accuracy, # 总体准确率，表示正确预测的样本比例
    normalized_accuracy = normalized_accuracy,  # 归一化准确率，参考paper计算
    balanced_accuracy = balanced_accuracy,  # 平衡准确率，是每个类别召回率的平均值
    macro_precision = macro_precision,  # 宏平均精确率（Macro Precision），是每个类别的精确率（Precision）的平均
    macro_recall = macro_recall,  # 宏平均召回率（Macro Recall），是每个类别召回率的平均值
    macro_f1 = macro_f1,  # 宏平均 F1 分数（Macro F1），是每个类别 F1 分数的平均值，综合精确率和召回率
    macro_specificity = macro_specificity, #
    log_loss = log_loss,  # 对数损失（Log Loss），衡量预测概率与真实标签之间的差异
    n_coef = n_coef,  # 非零特征数（Non-Zero Features），表示选定 lambda 值下模型中非零系数的特征数量
    deviance_min = deviance_min,  # 交叉验证中最小偏差（Deviance）
    deviance_1se = deviance_1se,  # 交叉验证中误差在一个标准误范围内最简洁模型的偏差
    lambda_min = lambda_min,  # 交叉验证中使偏差最小的正则化参数 lambda
    lambda_1se = lambda_1se  # 交叉验证中误差在一个标准误范围内最简单模型的 lambda（最强正则化）
  )
}
all_results_file <- file.path(output_path, "all_evaluations_summary.rds")
saveRDS(evaluation_results, file = all_results_file)

## 统计信息：
stats <- do.call(rbind, lapply(evaluation_results, function(m) data.frame(
  m$accuracy,
  m$normalized_accuracy,
  m$balanced_accuracy,
  m$macro_precision,
  m$macro_recall,
  m$macro_f1,
  m$log_loss,
  m$n_coef,
  m$deviance_min,
  m$deviance_1se,
  m$lambda_min,
  m$lambda_1se
  )))
rownames(stats) <- lapply(evaluation_results, function(m) paste0('Aplha_',m$alpha))
cat("Models Comparison:")
print(stats)

## 创建评估结果比较表
comparison_df <- data.frame(
  Alpha = sapply(evaluation_results, function(x) x$alpha),
  Accuracy = sapply(evaluation_results, function(x) x$accuracy),
  Normalized_Accuracy = sapply(evaluation_results, function(x) x$normalized_accuracy),
  Balanced_Accuracy = sapply(evaluation_results, function(x) x$balanced_accuracy),
  Macro_Precision = sapply(evaluation_results, function(x) x$macro_precision),
  Macro_Recall = sapply(evaluation_results, function(x) x$macro_recall),
  Macro_F1 = sapply(evaluation_results, function(x) x$macro_f1),
  Log_Loss = sapply(evaluation_results, function(x) x$log_loss),
  N_Coefficients = sapply(evaluation_results, function(x) x$n_coef),
  CV_Deviance_Min = sapply(evaluation_results, function(x) x$deviance_min),
  CV_Deviance_1se = sapply(evaluation_results, function(x) x$deviance_1se),
  Lambda_Min = sapply(evaluation_results, function(x) x$lambda_min),
  Lambda_1SE = sapply(evaluation_results, function(x) x$lambda_1se)
)
## 保存比较表
comparison_file <- file.path(output_path, "model_comparison.csv")
write.csv(comparison_df, file = comparison_file, row.names = FALSE)

# 找到最佳模型（基于平衡准确率）
best_model_idx <- which.max(comparison_df$Balanced_Accuracy)
best_alpha <- comparison_df$Alpha[best_model_idx]
cat(sprintf("\n最佳模型: Alpha = %.2f (平衡准确率: %.4f)\n", 
            best_alpha, comparison_df$Balanced_Accuracy[best_model_idx]))
