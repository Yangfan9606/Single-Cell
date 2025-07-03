### Single-Cell Multi-Omics Analysis Pipeline

This repository contains a comprehensive pipeline for integrated analysis of single-nucleus RNA-seq (snRNA-seq) and single-nucleus ATAC-seq (snATAC-seq) data, based on methods described in [Morabito et al., Nat Genet 53, 1143–1155 (2021)](https://doi.org/10.1038/s41588-021-00894-z).

## Overview
The pipeline consists of multiple scripts:
- `run.sh` : **Data preprocessing** with Cell Ranger
- `1.snRNA.R` : **snRNA-seq data processing** and cell annotation
- `2.snATAC.R` : **snATAC-seq data processing**
- `3.Integrated_snRNA_snATAC.R` : **Integrated analysis** of both modalities
- `4.Integrated_Cicero_CoAccessibility.R` : **Chromatin co-accessibility analysis** with Cicero
- `5.Integrated_chromVAR_Motif_activity.R` : **TF motif and ChromVAR analysis**

## Prerequisites
- 10x Genomics Cell Ranger
- R (v4.4.3 or higher)
- Required R packages: Seurat, Signac, Cicero, chromVAR, and others listed in the scripts

## Step-by-Step Guide

### Step 0: Data Preprocessing (Cell Ranger)
Convert raw FASTQ files to expression/accessibility matrices:

1. Download reference genomes from [10x Genomics](https://www.10xgenomics.com/software)
2. Edit `run.sh` with your paths.

**Outputs:**
| File Path | Description |
|-----------|-------------|
| `sample_rna/outs/` | RNA expression data |
| `sample_atac/outs/` | ATAC accessibility data |

### Step 1: snRNA-seq Data Analysis
Process snRNA-seq data and perform cell annotation:

Edit `1.snRNA.R` with your paths and sample name.

**Key operations:**
- Quality control filtering
- Normalization and dimensionality reduction
- Cell clustering using Seurat
- Cell type annotation with canonical markers

**Outputs:**
| File Path | Description |
|-----------|-------------|
| `CTL.snRNA_processed.rds` | Annotated Seurat object |
| `plots/snRNA_QC_before_filtering.pdf` | Pre-filtering QC metrics |
| `plots/snRNA_elbow_plot.pdf` | PCA elbow plot |
| `plots/snRNA_initial_clusters.pdf` | Initial clustering results |
| `plots/snRNA_marker_expression.pdf` | Marker gene expression |
| `plots/snRNA_annotated_celltypes.pdf` | Final annotated cell types |
> **Note**: The output sample name "CTL" can be modified by changing input file names and output object name in the script.

### Step 2: snATAC-seq Data Analysis
Process snATAC-seq data

Edit `2.snATAC.R` with your paths and sample name.

**Key operations:**
- Quality control filtering
- Dimensionality reduction
- Clustering and visualization

**Outputs:**
| File Path | Description |
|-----------|-------------|
| `CTL.snATAC_processed.rds` | Annotated ATAC Seurat object |
| `plots/snATAC_QC_metrics.pdf` | Violin plots of QC metrics |
| `plots/snATAC_LSI_depth_correlation.pdf` | Correlation plot between LSI components and sequencing depth |
| `plots/snATAC_initial_clusters.pdf` | UMAP visualization of initial cell clusters |
> **Note**: The output sample name "CTL" can be modified by changing input file names and output object name in the script.

### Step 3: Integrated snRNA-seq and snATAC-seq Analysis
Integrate both modalities:

Edit `3.Integrated_snRNA_snATAC.R` with your paths and sample name.

**Key operations:**
- Gene activity calculation
- Cross-modal integration
- High-confidence cell filtering
- Joint embedding creation

**Outputs:**
| File Path | Description |
|-----------|-------------|
| `CTL.joint_integrated.rds` | Integrated Seurat object containing both modalities |
| `plots/snATAC_prediction_scores.pdf` | Violin plots of prediction scores by cell type |
| `plots/RNA_vs_ATAC_annotations.pdf` | Side-by-side UMAPs comparing RNA and ATAC annotations |
| `plots/joint_embedding_by_dataset.pdf` | Joint UMAP split by dataset (RNA/ATAC) |
| `plots/joint_embedding_by_celltype.pdf` | Joint UMAP colored by cell type |
> **Note**: The output sample name "CTL" can be modified by changing input file names and output object name in the script.

### Step 4: Chromatin Co-Accessibility Analysis (Cicero)
Identify long-range chromatin interactions:

Edit `4.Integrated_Cicero_CoAccessibility.R` with your paths and sample name.

**Key operations:**
- Cicero Analysis：
  - Create Cicero cds using UMAP coordinates
  - Find co-accessible peak pairs within a 500kb window
- Filtering and annotation:
  - Filter connections by co-accessibility threshold (default: 0.15)
  - Calculate genomic distance between connected peaks
- Gene-specific analysis:
  - For each target gene (predefined list), find connections within a 1Mb window centered on the gene

**Outputs:**
| File Path | Description |
|-----------|-------------|
|`results/cicero/all_connections.csv`|	All detected peak-peak connections (raw) |
|`results/cicero/significant_connections.csv`|	Filtered connections (coaccess ≥ 0.15) with genomic coordinates |
|`results/cicero/gene_specific_connections.rds`|	Gene-targeted connection sets |
|`results/cicero/connection_summary.csv`|	Statistics per target gene: Connection count + Max coaccess score + Mean distance |
|`results/cicero/cicero_results.rds`|	Complete analysis objects |
|`plots/cicero/connection_distance_distribution.pdf`|	Histogram of connection distances |
|`plots/cicero/coaccess_score_distribution.pdf`|	Histogram of coaccess scores |
|`plots/cicero/connections_[GENE].pdf`|	Co-accessibility plots around target genes |
|`results/cicero/analysis_report.txt`|	Summary report with key statistics |

### Step 5: TF Motif Activity and ChromVAR Analysis
Identify cell type-specific transcription factor activity:

Edit `5.Integrated_chromVAR_Motif_activity.R` with your paths and sample name.

**Key operations:**
- Motif database processing:
  - Retrieves JASPAR2024 CORE human motifs
  - Filters motifs by information coefficient (IC ≥ 8):
- Chromatin accessibility analysis
- Cell-type specific TF activity
- Differential TF activity
- Target TF validation

**Outputs:**
| Category | File Path | Description |
|----------|-----------|-------------|
| Results | `results/motif/motif_mapping.csv` |	Filtered motifs with IC scores |
| | `results/motif/average_tf_activity.csv` |	Mean TF activity per cell type |
| | `results/motif/differential_TF_*.csv` |	Pairwise differential activity results |
| | `results/motif/target_tf_analysis.rds` |	Detailed target TF analysis |
| | `results/motif/tf_correlation_matrix.csv` |	TF co-activity network |
| | `results/motif/atac_with_chromvar.rds` |	ChromVAR-enhanced Seurat object |
| Visualizations | `plots/motif/tf_activity_heatmap.pdf` |	Top 50 variable TFs (z-scores) |
| | `plots/motif/target_tf_violin_plots.pdf` |	Activity distribution for target TFs |
| | `plots/motif/target_tf_umap_plots.pdf` |	Spatial activity patterns |
| | `plots/motif/differential_tf_barplot_*.pdf` |	Top differential TFs per comparison |
| | `plots/motif/tf_activity_boxplot_*.pdf` |	Cell-type specific activity for target TFs |
| | `plots/motif/tf_correlation_heatmap.pdf` |	TF co-activity network |
| | `plots/motif/tf_marker_validation.pdf` |	Known marker validation |

### Step 6: DAR and DEG Analysis
Perform differential accessibility (DAR) and differential expression (DEG) analysis:

Edit `6.Integrated_DAR_DEG_analysis.R` with your paths and sample name.
> **Note**: Replace <comparison> with comparison names (e.g., Astrocytes_vs_Oligodendrocytes) and <direction> with up or down for enrichment results.

**Key operations:**
- DAR analysis:
  - Identifies differential accessibility regions (peaks) between cell types
  - Applies LR test with latent variable correction (`nCount_ATAC`, `nFeature_ATAC`)
- DEG analysis:
  - Identifies differential expression genes between cell types
  - Uses Wilcoxon rank-sum test
- Visualization:
  - Volcano plots, MA plots, and histograms for DARs/DEGs
  - Summary bar plots of significant features per comparison
- GO/KEGG analysis for DEGs (up/down-regulated)
- DAR-DEG integration:
  - Annotates DAR peaks to nearest genes
  - Correlates chromatin accessibility with gene expression
  - Identifies genes with both accessibility and expression changes
- Cell-type specificity:
  - Gene frequency analysis across comparisons
  - Overlap statistics between DAR/DEG results
 
**Outputs:**
| Category | File Path | Description |
|----------|-----------|-------------|
| Results | `results/dar_deg/DAR_all_<comparison>.csv` |	Full DAR results per comparison |
| | `results/dar_deg/DAR_significant_<comparison>.csv` |	Significant DARs (FDR < 0.05) |
| | `results/dar_deg/DAR_summary.csv` |	Summary statistics for DAR analysis |
| | `results/dar_deg/DEG_all_<comparison>.csv` |	Full DEG results per comparison |
| | `results/dar_deg/DEG_significant_<comparison>.csv` |	Significant DEGs (FDR < 0.05) |
| | `results/dar_deg/DAR_DEG_association_<comparison>.csv` |	Genes with correlated DAR/DEG changes |
| | `results/dar_deg/gene_frequency_summary.csv` |	Gene occurrence frequency across comparisons |
| | `results/dar_deg/GO_enrichment_<comparison>_<direction>.csv` |	GO enrichment for DEGs (up/down) |
| | `results/dar_deg/KEGG_enrichment_<comparison>_<direction>.csv` |	KEGG enrichment for DEGs (up/down) |
| Visualizations | `plots/dar_deg/DAR_volcano_<comparison>.pdf` |	Volcano plot of DAR results |
| | `plots/dar_deg/DAR_MA_<comparison>.pdf` |	MA plot of accessibility vs. fold-change |
| | `plots/dar_deg/DAR_distribution_<comparison>.pdf` |	Histogram of significant DAR logFC |
| | `plots/dar_deg/DEG_volcano_<comparison>.pdf` |	Volcano plot of DEG results (gene labels) |
| | `plots/dar_deg/DAR_summary_counts.pdf` |	Bar plot of significant DAR counts per comparison |
| | `plots/dar_deg/DAR_DEG_association_<comparison>.pdf` |	Scatter plot: DAR vs. DEG logFC + peak annotations |
| | `plots/dar_deg/GO_enrichment_<comparison>_<direction>.pdf` |	Dot plot of GO enrichment |
| | `plots/dar_deg/gene_frequency_comprehensive_analysis.pdf` |	Combined gene frequency analysis plots |
