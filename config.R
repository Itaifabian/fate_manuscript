# =============================================================================
# config.R — T47D Fate Manuscript: All user-configurable paths and parameters
# =============================================================================
# Edit this file before running any pipeline script.
# All scripts source this file from the project root.

# -----------------------------------------------------------------------------
# INPUT DATA
# -----------------------------------------------------------------------------
TENX_DIR        <- "/path/to/your/filtered_feature_bc_matrix"
DATA_DIR        <- "data"
HBCA_H5AD       <- "/path/to/your/normal_reference.h5ad"
LINEAGE_DB_PATH <- "data/lineage.db.csv"
WES_DIR         <- "/path/to/your/WES_DIR"
GTF_FILE        <- "/path/to/your/gencode.v38.annotation.gtf"

# -----------------------------------------------------------------------------
# INTERMEDIATE DATA (written by early steps, read by later ones)
# -----------------------------------------------------------------------------
INFERCNV_DIR <- "/path/to/your/infercnv"

# -----------------------------------------------------------------------------
# OUTPUT
# -----------------------------------------------------------------------------
OUTPUT_DIR <- '/path/to/your/output_dir'
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
