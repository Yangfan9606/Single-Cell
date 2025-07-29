#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-07-29
# version:   1.0
# license:   MIT
# brief:     Cell type prediction use output of ElasticNet_classifier.R.
#------------------------------------------#

## ====================================================
##  Cell Type Predictor (elastic net)
## ====================================================
library(dplyr)
library(Seurat)
library(glmnet)

input_dir <- ('D:/Data/0_yangfan/0_SingleCell_RNA/scRNA_analysis_results/data')
my_obj = readRDS(paste0(input_dir,'/seurat_clustered_object.rds'))
output_dir <- 'D:/Data/0_yangfan/0_SingleCell_Scripts_Paper_2024/celltype_classifier_output/elasticnet_classifier_output'

seurat_obj = my_obj

## Load the trained classifier
classifier_path = file.path(output_dir, "ElasticNet_classifier_a0.00_l200.rds")
best_model <- readRDS(classifier_path)
required_genes <- rownames(best_model$model$glmnet.fit$beta[[1]])

## Parameters
confidence_threshold = 0.9
entropy_threshold = 0.5
output_col = "predicted_celltype"
confidence_col = "prediction_confidence"
entropy_col = "prediction_entropy"


## Arrange design matrix
feature.space <- intersect(required_genes, rownames(seurat_obj))
missing.features <- setdiff(required_genes, rownames(seurat_obj))

## Ensure the data is processed the same way as training data
DefaultAssay(seurat_obj) = "RNA"
temp <- NormalizeData(seurat_obj, verbose = F) %>% ScaleData(features=feature.space, verbose = F)
## Extract scaled data
x_new <- t(GetAssayData(temp,assay = 'RNA', layer = "scale.data"))
rm(temp)
## 检查并修正基因特征匹配问题
missing_genes <- setdiff(required_genes, colnames(x_new))
extra_genes <- setdiff(colnames(x_new), required_genes)

## 处理缺失基因（填充0）
if (length(missing_genes) > 0) {
  fill_matrix <- matrix(
    0,
    nrow = nrow(x_new),
    ncol = length(missing_genes),
    dimnames = list(rownames(x_new), missing_genes)
  )
  x_new <- cbind(x_new, fill_matrix)
}

## 移除多余基因并按训练顺序排序
x_new <- x_new[, required_genes, drop = FALSE]

## 确认维度匹配
stopifnot(ncol(x_new) == length(required_genes))
stopifnot(all(colnames(x_new) == required_genes))

## 进行预测
pred_prob_raw <- predict(best_model$model, newx = x_new, s = "lambda.1se", type = "response")

## Make predictions based on algorithm
pred_prob <- pred_prob_raw[,,1]
pred_class <- colnames(pred_prob)[apply(pred_prob, 1, which.max)]

## Calculate confidence and entropy
max_prob <- apply(pred_prob, 1, max)
entropy <- apply(pred_prob, 1, function(p) {
  p <- p[p > 0]
  if(length(p) == 0) return(0)
  -sum(p * log2(p))
})
## Apply thresholds for uncertain predictions
uncertain_pred <- max_prob < confidence_threshold | entropy > entropy_threshold
pred_class[uncertain_pred] <- "Uncertain"  
## Add predictions to metadata
seurat_obj@meta.data[[output_col]] <- pred_class
seurat_obj@meta.data[[confidence_col]] <- max_prob
seurat_obj@meta.data[[entropy_col]] <- entropy  
message("Predictions completed. ", sum(uncertain_pred), " cells marked as 'Uncertain'")  

head(seurat_obj$predicted_celltype)
DimPlot(seurat_obj, group.by = "predicted_celltype", label=T)
seurat_obj$prediction_entropy
seurat_obj$prediction_confidence
