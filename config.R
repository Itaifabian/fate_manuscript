# =============================================================================
# config.R — T47D Fate Manuscript: All user-configurable paths and parameters
# =============================================================================
# Edit this file before running any pipeline script.
# All scripts source this file from the project root.

# -----------------------------------------------------------------------------
# INPUT DATA
# -----------------------------------------------------------------------------
TENX_DIR        <- "/yoren64/users/shared_files/Itai_for_marina/filtered_feature_bc_matrix"
HBCA_H5AD       <- "/yoren64/users/itaifabian/cell_lines/T47D_analysis/cnv/normal_reference.h5ad"
LINEAGE_DB_PATH <- "/yoren64/users/shared_files/Itai_for_marina/lineage.db.csv"
TF_NETWORK      <- "/yoren64/users/itaifabian/T47D_analysis/decouple/human_prior_tri.csv"
WES_DIR         <- "/yoren64/users/itaifabian/cell_lines/T47D_analysis/WES_new/analysis/T47D_Orig"
GTF_FILE        <- "/yoren64/users/itaifabian/OSI Bulk RNA/gencode.v38.annotation.gtf"

# -----------------------------------------------------------------------------
# INTERMEDIATE DATA (written by early steps, read by later ones)
# -----------------------------------------------------------------------------
INFERCNV_DIR <- "/yoren64/users/itaifabian/infercnv"

# -----------------------------------------------------------------------------
# OUTPUT
# -----------------------------------------------------------------------------
OUTPUT_DIR <- "/yoren64/users/itaifabian/fate_manuscript/results"
FIGURES_DIR   <- file.path(OUTPUT_DIR, "figures")
table_dir <- file.path(OUTPUT_DIR, "supplementary_tables")

# -----------------------------------------------------------------------------
# INTERMEDIATE DATA (saved Seurat object — written by 01, read by 04–06)
# -----------------------------------------------------------------------------
SEURAT_RDS <- "/path/to/your/T47D_processed.rds"           # edit this

# -----------------------------------------------------------------------------
# ADDITIONAL ANALYSIS DATA
# -----------------------------------------------------------------------------
CORALEEN_DIR        <- "/path/to/your/CORALEEN"            # edit this
HYBRID_DATASETS_DIR <- "/path/to/your/hybrid_datasets"     # edit this

# -----------------------------------------------------------------------------
# COMPUTE
# -----------------------------------------------------------------------------
NUM_THREADS <- 16
