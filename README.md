# T47D Breast Cancer Fate Manuscript — scRNA-seq Analysis Pipeline

Single-cell RNA-seq analysis pipeline for characterising clonal fate in T47D breast cancer cells across a treatment time course (day 0, 7, 15). Integrates transcriptomics, CNV inference, lineage barcoding, and external cohort validation to identify transcriptional predictors of hybrid/senescence fate.

## Quick Start

1. Edit `config.R` and set every path to match your system.
2. Set your R working directory to the project root (`setwd("<path>/fate_manuscript")`).
3. Run scripts sequentially in the order given below, in the **same R session** (objects produced early are consumed later).

## Pipeline Overview

| Script | Step | What it does | Key outputs |
|--------|------|-------------|-------------|
| `R/01_preprocessing.R` | 1 | Load 10x data, QC, doublet removal, Seurat workflow, NMF metaprograms, cell-state clustering, pathway enrichment | `T47D` Seurat object, Figure 1 panels |
| `R/02_infercnv.R` | 2 | InferCNV on naive T47D cells vs. HBCA normal epithelial reference | `INFERCNV_DIR/run.final.infercnv_obj`, `CNV_table.csv`, Figure S1E |
| `R/03_cnv_subclonal.R` | 3 | Integrate WES CNV + inferCNV to resolve subclonal architecture | `Fig1F_clonal_expression`, `FigS1D_WES_clonal` |
| `R/04_trajectory.R` | 4 | Hybrid cell identification, DE, TF activity (decoupleR), Monocle3 pseudotime, signature scores | Figure 2 panels |
| `R/05_lineage_fate.R` | 5 | Lineage-barcode integration, CNV identity, XGBoost confidence filtering, clonal diversity, fate assignment | Figure 3 panels |
| `R/06_fate_priming.R` | 6 | Clonal compactness permutation test, fate-predictive markers, volcano, Euler diagram, violin plots | Figure 4 panels |
| `R/07_coraleen_signature.R` | 7 | Elastic-net signature predicting surgery ROR-S from baseline; bootstrap stability, LOOCV | Figure 5 panels, `pos_genes` / `neg_genes` objects |
| `R/08_neopalana_validation.R` | 8 | External validation in GSE93204 (NeoPalAna cohort) using signature from step 7 | Figure 5E panels |

### Additional Analysis

`R/additional_analysis/04_hybrid_meta_analysis.R` — Standalone meta-analysis across seven published scRNA-seq datasets (GSE131984, GSE134836, GSE146544, GSE168668, GSE234181, GSE280381, GSE297393). Applies the hybrid-detection logic to each cohort and produces **Figure 2F** (bubble plot of hybrid proportions across treatment conditions). **Does not require `T47D` in session**; runs independently of the main pipeline. Set `HYBRID_DATASETS_DIR` in `config.R` to the parent directory containing one subdirectory per dataset.

### Shared utilities

`R/utils/functions.R` — sourced automatically by each script. Contains:

- `filter_de_hits()` — filter a `FindMarkers` result table by adj. p-value, log2FC, and expression ratio
- `nmf_cluster_cells()` — cross-replicate robust NMF metaprogram discovery
- `save_figure()` / `save_pheatmap()` / `save_euler()` — save plot objects as PDF + PNG
- `save_table()` — save a data frame as CSV (+ XLSX if `writexl` is installed)
- `run_qc_and_doublet_removal()` — MAD-based QC + DoubletFinder on a list of Seurat objects

## Session Dependencies

Scripts 04–06 expect `T47D` (produced by 01) to already exist in the session.  
Script 08 expects `pos_genes` and `neg_genes` (produced by 07) to be in the session.  
Script 05 reads `CNV_table.csv` written to `INFERCNV_DIR` by script 02.  
Script 06 calls `filter_de_hits()` (defined in `R/utils/functions.R`, loaded via `source`).  
CNV clone labels (`cnv_clone`) are added to `T47D` in script 05 after inferCNV.

## Configuration (`config.R`)

All file paths and compute settings live in a single `config.R` at the project root. Edit every variable before running. The key variables:

| Variable | Description |
|----------|-------------|
| `TENX_DIR` | Filtered 10x feature-barcode matrix directory |
| `HBCA_H5AD` | HBCA normal epithelial reference `.h5ad` |
| `LINEAGE_DB_PATH` | Lineage barcode database CSV |
| `TF_NETWORK` | `human_prior_tri.csv` for decoupleR ULM |
| `WES_DIR` | WES pipeline output root (CNVkit + PyClone) |
| `GTF_FILE` | Gencode v38 annotation GTF |
| `SEURAT_RDS` | Path where the processed Seurat object is saved/loaded |
| `INFERCNV_DIR` | InferCNV output directory (also source of `CNV_table.csv`) |
| `OUTPUT_DIR` | Root results folder for all figures and tables |
| `CORALEEN_DIR` | Directory containing `merged_metadata.csv` and `expression_matrix_counts.csv` for the CORALEEN bulk RNA-seq cohort (used by script 07) |
| `HYBRID_DATASETS_DIR` | Parent directory containing one subdirectory per external scRNA-seq dataset (used by `additional_analysis/04_hybrid_meta_analysis.R`) |
| `NUM_THREADS` | Threads for inferCNV |

## Output Structure

```
OUTPUT_DIR/
├── figures/
│   ├── Figure_1/          ← scripts 01, 02 (FigS1E), 03
│   ├── Figure_2/          ← scripts 04, additional_analysis/04
│   ├── Figure_3/          ← script 05
│   ├── Figure_4/          ← script 06
│   └── Figure_5/
│       └── supp/
└── supplementary_tables/
```

WES/CNV plots from script 03 also write to `WES_DIR/variants/`.

## Dependencies

R version: **4.3.3**

### CRAN / R packages

| Package | Version |
|---------|---------|
| Seurat | 5.3.0 |
| patchwork | 1.3.1 |
| ggplot2 | 3.5.2 |
| ggpubr | 0.6.1 |
| ggrepel | 0.9.6 |
| tidyverse | 2.0.0 |
| dplyr | 1.1.4 |
| tidyr | 1.3.1 |
| tibble | 3.3.0 |
| purrr | 1.1.0 |
| readr | 2.1.5 |
| stringr | 1.5.1 |
| scales | 1.4.0 |
| plyr | 1.8.9 |
| data.table | 1.17.8 |
| Matrix | 1.6-5 |
| RcppML | 0.3.7 |
| igraph | 2.1.4 |
| pheatmap | 1.0.13 |
| RColorBrewer | 1.1-3 |
| fastcluster | 1.3.0 |
| parallelDist | 0.2.7 |
| vegan | 2.7-1 |
| xgboost | 3.0.5 |
| eulerr | 7.0.2 |
| glmnet | 4.1-10 |
| DoubletFinder | 2.0.4 |
| monocle3 | 1.3.1 |
| msigdbr | 26.1.0 |
| presto | 1.0.0 |
| decoupleR | 2.14.0 |
| SeuratWrappers | 0.4.0 |
| ggpmisc | 0.7.0 |
| tidytext | 0.4.3 |
| writexl | 1.5.4 |
| mclust | 6.1.2 |
| jsonlite | 2.0.0 |
| mclust | 6.1.2 |
| mclust | 6.1.2 |
| mclust | 6.1.2 |

### Bioconductor packages

| Package | Version |
|---------|---------|
| fgsea | 1.28.0 |
| limma | 3.58.1 |
| destiny | 3.16.0 |
| infercnv | 1.20.0 |
| GEOquery | 2.70.0 |
| genefu | 2.34.0 |
| biomaRt | 2.58.0 |
| fishplot | 0.5.3 |
| SingleCellExperiment | 1.24.0 |
| DESeq2 | — |
