#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-07-29
# version:   1.0
# license:   MIT
# brief:     Dataset preparation for traning.
#------------------------------------------#

## ====================================================
##  Cell Type prediction Data Preparation 
## ====================================================

library(Seurat)
library(dplyr)
library(caret)

setwd('D:/Data/0_yangfan/0_SingleCell_Scripts_Paper_2024')
output_path = 'D:/Data/0_yangfan/0_SingleCell_Scripts_Paper_2024/celltype_classifier_output'
seurat_obj = readRDS("Reference_MTG_RNAseq_all-nuclei.2022-06-07.rds")
table(seurat_obj$subclass_label)

#### L5/6 NP、L5 IT、L6 CT、L4 IT、L2/3 IT、L6 IT Car3、L6 IT、L6b、L5 ET：这些是皮层不同层的兴奋性神经元，分别位于5/6层、5层、6层、4层、2/3层等，负责在皮层内部或与其他脑区（如纹状体、丘脑）通信。例如，L5 IT神经元连接纹状体，参与运动控制。
#### Vip、Sst Chodl、Sncg、Pvalb、Lamp5 Lhx6、Sst、Lamp5、Chandelier：这些是抑制性中间神经元，调节兴奋性神经元的活动。例如，Pvalb神经元快速抑制信号，维持大脑平衡
Excitatory_neuron = c("L5/6 NP", "L5 IT", "L6 CT", "L4 IT", "L2/3 IT", "L6 IT Car3", "L6 IT", "L6b", "L5 ET")
Inhibitory_neuron = c("Vip", "Sst Chodl", "Sncg", "Pvalb", "Lamp5 Lhx6", "Sst", "Lamp5", "Chandelier")
## 归为一类 | 提取元数据并修改
new_metadata <- seurat_obj@meta.data %>%
  mutate(cell_type = case_when(
    subclass_label %in% Excitatory_neuron ~ "Excitatory_neuron",
    subclass_label %in% Inhibitory_neuron ~ "Inhibitory_neuron",
    TRUE ~ as.character(subclass_label)  # 保留其他标签不变
  )) 
# 将修改后的元数据添加回 Seurat 对象
seurat_obj$cell_type <- new_metadata$cell_type

## 选取 target cells
target_cells <- c("Oligodendrocyte",  # 少突胶质细胞：为轴突提供髓鞘，加速信号传导
                  "Astrocyte",  # 星形胶质细胞：支持神经元，维护血脑屏障，帮助神经递质回收
                  "Microglia-PVM",  # 小胶质细胞-血管周围巨噬细胞：大脑的免疫细胞，清除损伤或感染
                  "OPC",  # 少突胶质前体细胞：可分化为少突胶质细胞，参与髓鞘修复
                  "Endothelial", # 内皮细胞：构成血管内壁，维护血脑屏障
                  "Pax6",  # 一种转录因子，更多用于标记神经发育中的祖细胞，而非具体细胞类型
                  "VLMC", # 血管脑膜细胞：与血管系统相关，支持大脑血流
                  "Excitatory_neuron", "Inhibitory_neuron")  

seurat_subset <- subset(seurat_obj, cell_type %in% target_cells)
rm(seurat_obj)
gc()
table(seurat_subset$cell_type)
## 随机选择 500 个细胞的 ID
sampled_cells <- sample(colnames(seurat_subset), size = 40000, replace = FALSE)
seurat_downsampled <- subset(seurat_subset, cells = sampled_cells)
table(seurat_downsampled$cell_type)
##  prepare_celltype_data 

seurat_obj = seurat_subset
seurat_obj = seurat_downsampled
filter_genes = "original" # c("original", "hvg", "hvg_filter")
n_features = 700
remove_prefixes = c("^(AL|AC|LINC)\\d+", "^MT-", "^RPS", "^RPL")  # Alu element, Actin, Long Intergenic Non-coding | 技术噪音 (重复序列干扰) | 管家基因（无特异性，高表达）| 低信息量（低表达、高稀疏性） | Added mitochondrial and ribosomal genes
train_ratio = 0.75
stratify = TRUE
weight_method = "balanced" # c("inverse_freq", "balanced", "none")
cell_type_col = "cell_type"
min_cells_per_type = 10  # Filter rare cell types
output_path = NULL
seed = 42

## Clear scale data to save memory
seurat_obj@assays$RNA@scale.data <- matrix(0)
# Filter rare cell types
cell_type_counts <- table(seurat_obj@meta.data[[cell_type_col]])
valid_types <- names(cell_type_counts)[cell_type_counts >= min_cells_per_type]
seurat_obj <- subset(seurat_obj, cells = which(seurat_obj@meta.data[[cell_type_col]] %in% valid_types))

message("Filtered ", length(cell_type_counts) - length(valid_types), " rare cell types (< ", min_cells_per_type, " cells)")
message("  Process ", length(valid_types), " cell types")

## Gene filtering
if (filter_genes %in% c("original", "hvg_filter")) {
  keep_genes <- !grepl(paste(remove_prefixes, collapse = "|"), rownames(seurat_obj))
  seurat_obj <- subset(seurat_obj, features = rownames(seurat_obj)[keep_genes])
}
## Data processing with improved parameters
if (filter_genes == "original") {
  seurat_obj <- seurat_obj %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = n_features, verbose = FALSE) %>%
    ScaleData(verbose = FALSE)  # ORIGINAL: Process all genes
} else {
  seurat_obj <- seurat_obj %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = n_features, verbose = FALSE) %>%
    ScaleData(features = VariableFeatures(seurat_obj), verbose = FALSE) # IMPROVED: Use highly variable genes
}
  
## Extract data
x <- t(seurat_obj@assays$RNA@scale.data)
y <- as.character(seurat_obj@meta.data[[cell_type_col]])


## Remove NA labels
x <- x[!is.na(y), ]
y <- y[!is.na(y)]  
  
## Data splitting with improved stratification
if (stratify) {
  ## IMPROVED: Stratified sampling | 分层抽样 (确保训练集/测试集中各类别比例与原始数据一致) | 分类问题（尤其不均衡数据）必须使用分层抽样；回归问题或均衡数据可随机抽样。
  ## library(caret)
  train_idx <- createDataPartition(y, p = train_ratio, list = FALSE)[,1]
} else {
  train_idx <- sample(seq_len(nrow(x)), size = floor(train_ratio * nrow(x)))  # ORIGINAL: Random sampling | 简单随机抽样
}

train_x <- x[train_idx, ]
train_y <- y[train_idx]
test_x <- x[-train_idx, ]
test_y <- y[-train_idx]

## Calculate weights | 防止模型偏向多数类 | 提升少数类的识别能力 | 使损失函数更关注难分类样本
if (weight_method == "inverse_freq") {
  ## ORIGINAL: Inverse frequency weighting | 逆频率加权, 权重与类别频率成反比
  train_w <- 1 / (table(train_y)[train_y] * length(unique(train_y)))
  test_w <- 1 / (table(test_y)[test_y] * length(unique(test_y)))
} else if (weight_method == "balanced") {
  ## IMPROVED: Sklearn-style balanced weights | 平衡加权, 使每个类别的总权重相等（非样本权重相等） | 当某类样本极少时(如<5%)，可额外增加权重
  class_counts <- table(train_y)
  train_w <- length(train_y) / (length(class_counts) * class_counts[train_y])
  test_w <- length(test_y) / (length(unique(test_y)) * table(test_y)[test_y])
} else {
  train_w <- rep(1, length(train_y)) # No weighting
  test_w <- rep(1, length(test_y))
}

datasets <- list(
  train = list(x = train_x, y = train_y, w = as.numeric(train_w)),
  test = list(x = test_x, y = test_y, w = as.numeric(test_w)),
  metadata = list(
    n_features = ncol(x),
    n_cells_train = nrow(train_x),
    n_cells_test = nrow(test_x),
    n_classes = length(unique(y)),
    class_distribution = table(y),
    train_distribution = table(train_y),
    test_distribution = table(test_y),
    filter_method = filter_genes,
    weight_method = weight_method,
    seed = seed
  )
)  

saveRDS(datasets, file = file.path(output_path,'datasets.rds'))
