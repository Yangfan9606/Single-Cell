### Advanced Single-Cell Multi-Omics Analysis Pipeline

This repository contains an advanced single-cell multi-omics pipeline for integrated analysis of single-nucleus RNA-seq (snRNA-seq) and single-nucleus ATAC-seq (snATAC-seq) data, based on methods described in [Morabito et al., Nat Genet 53, 1143–1155 (2021)](https://doi.org/10.1038/s41588-021-00894-z).

## Overview
The pipeline consists of multiple scripts:
- `1.scRNA.process.R` : **snRNA-seq data processing** and multiple samples integration
- `2.scATAC.process.R` : **snATAC-seq data processing** and multiple samples integration
- `3.Integrate_Unmatched_snRNA_snATAC.R` : [UnFinished] **Integrated analysis** for **Unmatched** datasets
- `4.Integrate_Matched_snRNA_snATAC.R` : [UnFinished] **Integrated analysis** for **Matched** datasets
- `5.scRNA_only_Cell_Annotation.R` : **Cell type annotation** use scRNA data only
- `6.Integrate_Cell_Annotation.R` : [UnFinished] **Cell type annotation** use integrated data
- `7.scRNA_only_DEG.R` : **Differential gene expression analysis** use scRNA data
- `8.scRNA_only_Trajectory.R` : **Trajectory analysis** use scRNA data only
- `9.Integrate_Trajectory.R` : [UnFinished] **Trajectory analysis** use integrated data
- `10.scRNA_only_GRN.py` : [UnFinished] **Gene regulatory network Analysis** use scRNA data only
- `11.Integrate_GRN` : [UnFinished] **Gene regulatory network Analysis** use integrated data

## Step-by-Step Guide
### Step 1: scRNA-seq data processing
Process scRNA-seq data:

Edit `1.scRNA.process.R` with your paths and sample name.

**Key operations:**
- Data preprocessing: Load data from 10X Genomics outputs
- QC indicator calculation:
    - Mitochondrial gene ratio (percent.mt)
    - Ribosomal gene ratio (percent.ribo)
    - Dual cell detection: Use scDblFinder to identify and remove double cells
- Data integration:
    - SCTransform standardization: using glmGamPoi method
    - Select 5000 highly variable genes as integration features
    - RPCA integration method based on SCT
- PCA dimensionality reduction (1-30 principal components)
- Louvain clustering (resolution 0.5)

**Outputs:** `**snRNA_analysis_output/**`
| File Path | Description |
|-----------|-------------|
| `AD_vs_CTL_snRNA_integrated.rds` | Integrated Seurat object |
| `ad_sample_rna_scatter.pdf` | Dual cell detection scatter plot |
| `AD_vs_CTL_snRNA_elbow_plot.pdf` | PCA elbow plot |
| `AD_vs_CTL_snRNA_integration_UMAP.pdf` | Integrated UMAP plot |
| `AD_vs_CTL_snRNA_QC_features_UMAP.pdf` | QC indicator UMAP plot  |
| `ad_sample_rna_QC_data.txt` | Sample QC data |
| `AD_vs_CTL_snRNA_UMAP_coordinates.txt` | UMAP coordinates |
| `AD_vs_CTL_snRNA_integration_report.xlsx` | Analysis report |

### Step 2: scATAC-seq data processing
Process scATAC-seq data and perform peak annotation:

Edit `2.scATAC.process.R` with your paths and sample name.

**Key operations:**
- Data loading and processing:
    - Load ATAC-seq data from 10X Genomics outputs
    - Create ChromatinAssay objects for chromatin analysis
    - Sample metadata processing and integration
- Quality control:
    - Calculate multiple QC indicators (fragment length, TSS enrichment, nucleosome signal, etc.)
    - Generate a comprehensive QC dashboard
    - Cell filtering based on QC indicators
- Data integration and analysis:
    - Use Signac's integration method to merge multi-sample data
    - Dimensionality reduction processing (LSI)
    - Cluster analysis (Louvain algorithm)
    - UMAP visualization
- Advanced annotation:
    - Peak genomic position annotation
    - Gene activity score calculation
    - Detailed annotation visualization (genomic distribution, TSS distance, etc.)

**Outputs:** `**snATAC_analysis_output/**`
| File Path | Description |
|-----------|-------------|
| `integration_results/integrated_snatac_object.rds` | Integrated ATAC Seurat object |
| `plots/qc_dashboard.pdf` | Quality control dashboard |
| `plots/integration_assessment.pdf` | Integration assessment |
| `plots/peak_annotation_summary.pdf` | Peak annotation summary |
| `source_data/pre_filter_metadata.csv` | Pre-filter metadata |
| `source_data/post_filter_metadata.csv` | Post-filter metadata |
| `source_data/qc_metrics_data.csv` |QC metrics data |
| `source_data/peak_annotations.csv` | Peak annotation details |

### Step 3: Integrated unmatched scRNA and scATAC data

### Step 4: Integrated matched scRNA and scATAC data

### Step 5: Cell type annotation (scRNA only)

Edit `5.scRNA_only_Cell_Annotation.R` with your paths and sample name.

**Key operations:**
- Multi-strategy cell annotation:
  - Marker gene annotation: using predefined cell type-specific gene sets
  - SingleR automatic annotation: based on HumanPrimaryCellAtlas and Blueprint reference datasets
- Annotation optimization and integration:
  - Annotation simplification: unifying annotation naming from different sources
  - Voting system: integrating the results of multiple annotation methods
  - Priority strategy: manual annotation takes precedence over automatic annotation

**Outputs:** `snRNA_analysis_output/├──cell_annotation/`
| File Path | Description |
|-----------|-------------|
| `AD_vs_CTL_snRNA_annotated.rds`| Annotated Seurat object |
| `cell_type_umap.pdf` | UMAP distribution of cell types |
| `cell_type_distribution.pdf` | Distribution of cell types between groups |
| `marker_genes_heatmap.png`| Heatmap of marker gene expression |
| `SingleR_HPCA_quality_ggplot.pdf`| SingleR annotation quality |
| `umap_cell_types.txt`| UMAP coordinates + cell types |
| `cell_type_distribution.txt`| Cell type distribution data |
| `marker_genes_heatmap_data.txt`| Heatmap source data |
| `singler_hpca_scores.txt`| SingleR score |
| `cell_annotation_report.xlsx`| Excel comprehensive report |

### Step 6: Cell type annotation (scRNA + scATAC)

### Step 7: DEG (scRNA only)
Differential gene expression (DGE) analysis on scRNA-seq data using the Seurat package:

Edit `7.scRNA_only_DEG.R` with your paths and sample name.

**Key operations:**
- Data Preprocessing:
  - Loads and processes annotated Seurat objects
  - Checks cell type and group distributions
  - filters cell types with sufficient cell counts
- Differential Gene Expression (DGE): using the Wilcoxon test
- Visualization
- Functional Enrichment

**Outputs:** `snRNA_analysis_output/├──differential_analysis/`
| File Path | Description |
|-----------|-------------|
| **figures/** | |
| `cell_counts_barplot.pdf` | Bar plot of cell counts by cell type and group |
| `DEG_summary_barplot.pdf` | Bar plot of upregulated and downregulated DEGs per cell type |
| `volcano_plot_[celltype].png` | Volcano plot for each cell type showing log2 fold change vs. adjusted p-value |
| `GO_enrichment_[celltype].png` | Dot plot of GO enrichment results for each cell type |
| `DEG_heatmap.pdf` and `top50_DEG_heatmap.png` | Heatmaps of average expression for the top 50 DEGs |
| **tables/** | |
| `DEG_[celltype].txt` | DGE results for each cell type |
| `GO_enrichment_[celltype].txt` | GO enrichment results for each cell |
| **source_data/** | |
| `cell_counts_per_group.txt` | Cell counts by cell type and group |
| `DEG_summary.txt` | Summary of DEGs (total, upregulated, downregulated, significant) |
| `volcano_data_[celltype].txt` | Data used for volcano plots |
| `heatmap_data.txt` | Data used for heatmaps |
| **Reports:** | |
| `Differential_Analysis_Report.xlsx` | Excel file with parameters, analysis summary, DEG summary, and cell counts |
| `analysis_parameters.txt` | Text file summarizing key analysis parameters and results |


### Step 8: Trajectory (scRNA only)
Pseudotime trajectory analysis on scRNA-seq data using the Monocle3 package:

Edit `8.scRNA_only_Trajectory.R` with your paths and sample name.

**Key operations:**
- Data Loading and Preprocessing:
  - Loads an annotated Seurat object, subsets relevant cell types
  - Converts data to a Monocle3 cell data set (CDS) object
- Trajectory Analysis
- Visualization
- Gene Identification: Identifies genes significantly associated with the trajectory using Moran's I test

**Outputs:** `snRNA_analysis_output/├──trajectory_analysis/`
| File Path | Description |
|-----------|-------------|
| ├──**Saved Objects** | |
| `[project_name]_trajectory_cds.rds` | Monocle3 CDS object (commented out) |
| `[project_name]_trajectory_summary.rds` | Summary of trajectory analysis (commented out) |
| ├──**plots/** | |
| `trajectory_by_celltype.pdf` | Trajectory plot colored by cell type |
| `trajectory_by_partition.pdf` | Trajectory plot colored by partition |
| `trajectory_by_pseudotime.pdf` | Trajectory plot colored by pseudotime |
| `trajectory_by_pseudotime_grouped.pdf` | Trajectory plot by pseudotime, faceted by group |
| `trajectory_genes_expression.png` | Heatmap of expression for top 50 trajectory-associated genes |
| ├──**data_export/** | |
| `pseudotime_data.txt` | Cell-level pseudotime, cell type, group, cluster, and partition information |
| `trajectory_genes.txt` | Significant trajectory-associated genes with Moran's I test results |
| ├──**Reports:** | |
| `trajectory_analysis_report.xlsx` | Excel file with analysis parameters, cell type distribution, and trajectory summary |

### Step 9: Trajectory (scRNA + scATAC)

### Step 10: GRN (scRNA only)

### Step 11: GRN (scRNA + scATAC)
