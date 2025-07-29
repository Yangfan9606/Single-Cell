### Single-Cell Cell-type prediction based on Machine Learning (ML)

This repository contains R sciprts (classifier/predictor) for different algorithm.
Annotation data were download from: https://portal.brain-map.org/

## Current ML Algorithms: Elastic Net, Random Forest, and Support Vector Machine

```
├── Cell_type_prediction_data_preparation.R  # Data preprocessing and splitting
├── ElasticNet_classifier.R                  # Elastic Net model training
├── ElasticNet_predictor.R                   # Elastic Net prediction
├── RandomForest_classifier.R               # Random Forest model training
├── RandomForest_predictor.R                # Random Forest prediction
├── SVM_classifier.R                        # SVM model training
└── SVM_predictor.R                         # SVM prediction
```

## Dependencies
### Required R Packages
```r
# Core packages
library(Seurat)      # Single-cell analysis
library(dplyr)       # Data manipulation
library(caret)       # Machine learning utilities

# Algorithm-specific packages
library(glmnet)      # Elastic Net
library(randomForest) # Random Forest
library(e1071)       # SVM

# Evaluation and utilities
library(MLmetrics)   # Model evaluation metrics
library(doMC)        # Parallel processing
```

## Usage

### 1. Data Preparation

```r
source("Cell_type_prediction_data_preparation.R")
```

**Key Parameters:**
- `target_cells`: Cell types to include in classification
- `n_features`: Number of highly variable genes (default: 700)
- `train_ratio`: Training/test split ratio (default: 0.75)
- `weight_method`: Class balancing method ("balanced", "inverse_freq", "none")
- `min_cells_per_type`: Minimum cells required per cell type (default: 10)

**Features:**
- Automatic cell type grouping (excitatory/inhibitory neurons)
- Gene filtering (removes ribosomal, mitochondrial, and repetitive elements)
- Stratified sampling for balanced train/test splits
- Class weight calculation for imbalanced datasets

### 2. Model Training

#### Elastic Net
```r
source("ElasticNet_classifier.R")
```
- **Alpha values**: 0 (Ridge), 0.25, 0.5, 0.75, 1 (Lasso)
- **Cross-validation**: 10-fold CV with lambda optimization
- **Regularization**: Automatic lambda sequence generation

#### Random Forest
```r
source("RandomForest_classifier.R")
```
- **Trees**: 500 trees (configurable)
- **Features**: Automatic mtry optimization
- **Importance**: Built-in feature importance scoring

#### Support Vector Machine
```r
source("SVM_classifier.R")
```
- **Kernel**: Radial basis function (configurable)
- **Parameters**: Automatic gamma calculation
- **Cross-validation**: Built-in accuracy estimation

### 3. Model Prediction

#### Making Predictions
```r
# Load your single-cell data
seurat_obj <- readRDS("your_seurat_object.rds")

# Run prediction (example with Elastic Net)
source("ElasticNet_predictor.R")
```

**Quality Control Features:**
- **Confidence thresholds**: Filter low-confidence predictions
- **Entropy measures**: Identify uncertain classifications
- **Algorithm-specific metrics**: Decision distances (SVM), voting ratios (RF)

## Model Comparison

### Evaluation Metrics

| Metric | Description | Use Case |
|--------|-------------|----------|
| **Overall Accuracy** | Proportion of correct predictions | General performance |
| **Normalized Accuracy** | Average per-class accuracy | Balanced evaluation |
| **Macro F1-Score** | Harmonic mean of precision/recall | Imbalanced datasets |
| **Log Loss** | Probabilistic prediction quality | Confidence calibration |
| **Balanced Accuracy** | Average recall across classes | Class imbalance |

### Algorithm-Specific Metrics

#### Elastic Net
- **Non-zero coefficients**: Feature selection effectiveness
- **Lambda optimization**: Regularization strength
- **Cross-validation deviance**: Model complexity vs. performance

#### Random Forest
- **OOB Error**: Out-of-bag prediction accuracy
- **Feature Importance**: Variable significance ranking
- **Vote ratios**: Consensus strength across trees

#### Support Vector Machine
- **Support vectors**: Model complexity indicator
- **Decision distances**: Confidence near decision boundaries
- **Cross-validation accuracy**: Generalization estimate
