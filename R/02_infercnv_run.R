source("config.R")
source("R/utils/functions.R")

library(Seurat)
library(infercnv)
library(anndata)
library(Matrix)
library(biomaRt)
library(readr)
library(dplyr)
library(ggplot2)
library(ggnewscale)

# ─── CONFIG ───────────────────────────────────────────────────────────────────
HBCA_H5AD     <- HBCA_H5AD
fig_dir <- file.path(FIGURES_DIR, "Figure_1")

MIN_CELLS <- 50
CHR_ORDER <- paste0("chr", c(1:22, "X", "Y"))
# How many normal cells to subsample per cell type (500 is plenty for reference)
NORMAL_CELLS_PER_TYPE <- 500

# Column in HBCA obs that holds cell type labels
# From previous work: "level1_global" (LASP / LHS / BMYO)
HBCA_CELLTYPE_COL <- "level1_global"

# Ensembl mart (initialised here so it is available for both HBCA ID mapping
# in Part 2 and the gene-coordinate query in Part 4)
mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl",
                   mirror = "www")

# ──────────────────────────────────────────────────────────────────────────────

dir.create(INFERCNV_DIR, showWarnings = FALSE, recursive = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# PART 1 — Load T47D tumor counts
# ══════════════════════════════════════════════════════════════════════════════

cat("Loading T47D Seurat object...\n")
naive_T47D <- subset(T47D, time == 0)
naive_T47D <- FindVariableFeatures(naive_T47D)
naive_T47D <- NormalizeData(naive_T47D)
naive_T47D <- ScaleData(naive_T47D)
naive_T47D <- RunPCA(naive_T47D)
ElbowPlot(naive_T47D)
naive_T47D <- RunUMAP(naive_T47D, dims = 1:20)

cat(sprintf("  %d cells, %d genes\n", ncol(naive_T47D), nrow(naive_T47D)))

tumor_counts <- GetAssayData(naive_T47D, assay = "RNA", layer = "counts")
tumor_counts  <- as(tumor_counts, "dgCMatrix")

annotations_tumor <- data.frame(
  cell_type = rep("Tumor", ncol(tumor_counts)),
  row.names  = colnames(tumor_counts)
)

cat(sprintf("  Tumor cells: %d\n", ncol(tumor_counts)))

# ══════════════════════════════════════════════════════════════════════════════
# PART 2 — Load HBCA normal epithelial reference
# ══════════════════════════════════════════════════════════════════════════════

cat("\nLoading HBCA h5ad normal reference...\n")
ad <- read_h5ad(HBCA_H5AD)

cat("  HBCA obs columns:", paste(colnames(ad$obs), collapse = ", "), "\n")
cat("  HBCA cell type distribution:\n")
print(table(ad$obs[[HBCA_CELLTYPE_COL]]))

# Extract raw counts — HBCA stores normalized values in X, raw in raw slot
if (!is.null(ad$raw)) {
  cat("  Using ad$raw$X for raw counts\n")
  counts_normal_full <- t(ad$raw$X)
} else {
  cat("  WARNING: No raw slot found, using ad$X (verify these are integers)\n")
  counts_normal_full <- t(ad$X)
}

# Set gene names from var
if (!is.null(ad$raw)) {
  rownames(counts_normal_full) <- rownames(ad$raw$var)
} else {
  rownames(counts_normal_full) <- rownames(ad$var)
}
colnames(counts_normal_full) <- rownames(ad$obs)

cat(sprintf("  Raw check — value range: %.3f to %.3f (should be integers)\n",
            min(counts_normal_full@x), max(counts_normal_full@x)))

# Subsample normal cells per type
cat(sprintf("  Subsampling up to %d cells per cell type...\n", NORMAL_CELLS_PER_TYPE))
set.seed(42)
cell_types <- unique(ad$obs[[HBCA_CELLTYPE_COL]])
keep_cells <- unlist(lapply(cell_types, function(ct) {
  cells <- rownames(ad$obs)[ad$obs[[HBCA_CELLTYPE_COL]] == ct]
  sample(cells, min(NORMAL_CELLS_PER_TYPE, length(cells)))
}))

counts_normal <- counts_normal_full[, keep_cells]

# Label normals with prefix so they're clearly reference groups
annotations_normal <- data.frame(
  cell_type = paste0("normal_", ad$obs[keep_cells, HBCA_CELLTYPE_COL]),
  row.names  = keep_cells
)
cat("  Normal annotation distribution after subsampling:\n")
print(table(annotations_normal$cell_type))

cat("Converting HBCA Ensembl IDs to gene symbols...\n")

ensembl_ids <- rownames(counts_normal)

# Strip version numbers if present (ENSG00000243485.5 → ENSG00000243485)
ensembl_ids_clean <- gsub("\\..*$", "", ensembl_ids)

symbol_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = ensembl_ids_clean,
  mart       = mart
)

# Remove entries with empty symbols
symbol_map <- symbol_map[symbol_map$hgnc_symbol != "", ]
symbol_map <- symbol_map[!duplicated(symbol_map$ensembl_gene_id), ]

cat(sprintf("  Mapped %d / %d Ensembl IDs to symbols\n",
            nrow(symbol_map), length(ensembl_ids)))

# Filter and rename rows of counts_normal
keep_idx <- match(ensembl_ids_clean, symbol_map$ensembl_gene_id)
keep_mask <- !is.na(keep_idx)

counts_normal <- counts_normal[keep_mask, ]
rownames(counts_normal) <- symbol_map$hgnc_symbol[keep_idx[keep_mask]]

cat(sprintf("  counts_normal after renaming: %d genes x %d cells\n",
            nrow(counts_normal), ncol(counts_normal)))
head(rownames(counts_normal))

# ══════════════════════════════════════════════════════════════════════════════
# PART 3 — Combine tumor + normal on shared genes
# ══════════════════════════════════════════════════════════════════════════════

cat("\nMerging tumor and normal matrices...\n")
shared_genes <- intersect(rownames(tumor_counts), rownames(counts_normal))
cat(sprintf("  Shared genes: %d\n", length(shared_genes)))

counts_combined <- cbind(
  tumor_counts[shared_genes, ],
  counts_normal[shared_genes, ]
)
cat(sprintf("  Combined matrix: %d genes x %d cells\n",
            nrow(counts_combined), ncol(counts_combined)))

annotations_combined <- rbind(annotations_tumor, annotations_normal)

# Sanity check
stopifnot(all(colnames(counts_combined) == rownames(annotations_combined)))

annot_file <- file.path(INFERCNV_DIR, "cell_annotations.txt")
write.table(annotations_combined, annot_file,
            sep = "\t", quote = FALSE, col.names = FALSE)
cat(sprintf("  Annotations written to %s\n", annot_file))

# ══════════════════════════════════════════════════════════════════════════════
# PART 4 — Gene order file from Ensembl (hg38)
# ══════════════════════════════════════════════════════════════════════════════

cat("\nFetching gene coordinates from Ensembl (hg38)...\n")

gene_coords <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name",
                 "start_position", "end_position"),
  filters    = "hgnc_symbol",
  values     = rownames(counts_combined),
  mart       = mart
)

# Keep standard chromosomes only, one entry per gene
gene_coords <- gene_coords[
  gene_coords$chromosome_name %in% c(as.character(1:22), "X", "Y"), ]
gene_coords <- gene_coords[!duplicated(gene_coords$hgnc_symbol), ]
gene_coords <- gene_coords[order(
  match(gene_coords$chromosome_name, c(as.character(1:22), "X", "Y")),
  gene_coords$start_position), ]

gene_order_file <- file.path(INFERCNV_DIR, "gene_order.txt")
write.table(gene_coords, gene_order_file,
            sep = "\t", quote = FALSE,
            col.names = FALSE, row.names = FALSE)
cat(sprintf("  %d genes with coordinates written to %s\n",
            nrow(gene_coords), gene_order_file))

# Subset combined matrix to genes with coordinates
counts_combined <- counts_combined[
  rownames(counts_combined) %in% gene_coords$hgnc_symbol, ]
cat(sprintf("  Final matrix after gene filter: %d genes x %d cells\n",
            nrow(counts_combined), ncol(counts_combined)))

# ══════════════════════════════════════════════════════════════════════════════
# PART 5 — Create inferCNV object
# ══════════════════════════════════════════════════════════════════════════════

cat("\nCreating inferCNV object...\n")

# Reference group names must match the prefixed labels above
reference_groups <- unique(annotations_normal$cell_type)
cat("  Reference groups:", paste(reference_groups, collapse = ", "), "\n")

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts_combined,
  annotations_file  = annot_file,
  delim             = "\t",
  gene_order_file   = gene_order_file,
  ref_group_names   = reference_groups
)

# ══════════════════════════════════════════════════════════════════════════════
# PART 6 — Run inferCNV
# ══════════════════════════════════════════════════════════════════════════════

cat("\nRunning inferCNV...\n")

infercnv_obj <- infercnv::run(
  infercnv_obj,
  
  cutoff             = 0.1,
  out_dir            = INFERCNV_DIR,
  save_rds           = TRUE,
  no_plot            = TRUE,
  
  cluster_by_groups  = TRUE,
  cluster_references = FALSE,
  analysis_mode      = "subclusters",
  
  # Leiden but constrained — much faster than random_trees
  tumor_subcluster_partition_method = "leiden",
  leiden_method                     = "PCA",
  leiden_function                   = "CPM",
  leiden_resolution                 = 5e-5,   
  tumor_subcluster_pval             = 0.05,
  k_nn                              = 20,
  z_score_filter                    = 1.5,    # stricter gene filtering
  
  HMM                = TRUE,
  HMM_type           = "i6",
  BayesMaxPNormal    = 0.5,
  
  smooth_method      = "pyramidinal",
  window_length      = 151,
  denoise            = TRUE,
  sd_amplifier       = 1.5,
  
  num_threads        = NUM_THREADS,
  resume_mode        = FALSE
)

cat("Done.\n")



# ══════════════════════════════════════════════════════════════════════════════
# PART 7 — assign infercnv to the seurat object
# ══════════════════════════════════════════════════════════════════════════════


# First assign ALL subclusters including singletons
all_subclusters <- infercnv_obj@tumor_subclusters$subclusters[["Tumor"]]

all_flat <- unlist(lapply(names(all_subclusters), function(s) {
  setNames(rep(s, length(all_subclusters[[s]])),
           names(all_subclusters[[s]]))
}))

naive_T47D$infercnv_subcluster <- as.character(all_flat[colnames(naive_T47D)])

cat("All subcluster sizes:\n")
print(table(naive_T47D$infercnv_subcluster, useNA="ifany"))

# Now remove cells in subclusters below MIN_CELLS
subcluster_sizes <- table(naive_T47D$infercnv_subcluster)
valid_names      <- names(subcluster_sizes[subcluster_sizes >= MIN_CELLS])

naive_T47D <- naive_T47D[, naive_T47D$infercnv_subcluster %in% valid_names]

cat(sprintf("\nAfter removing subclusters < %d cells: %d cells remaining\n",
            MIN_CELLS, ncol(naive_T47D)))

table(naive_T47D$infercnv_subcluster)

DimPlot(naive_T47D, group.by = "infercnv_subcluster",
        cols = c('Tumor_s1' = "purple", 'Tumor_s2'= "darkgreen"))

df <- data.frame(
  cell = names(all_flat),
  label = unname(all_flat)
)

write_csv(df, paste0(INFERCNV_DIR, "/CNV_table.csv"))
# ══════════════════════════════════════════════════════════════════════════════
# PART 8 — Build gene table for plotting
# ══════════════════════════════════════════════════════════════════════════════
cat("Building gene table...\n")

cnv_mat        <- infercnv_obj@expr.data
infercnv_genes <- rownames(cnv_mat)   # already in chromosomal order

# Load gene coordinates from GTF for just these genes
# (uses the gene_order file already written above)
gene_order_raw <- read.table(gene_order_file, sep="\t", header=FALSE,
                             col.names=c("gene","chr","start","end"))

gene_table <- data.frame(
  gene       = infercnv_genes,
  stringsAsFactors = FALSE
)
gene_table$chr   <- paste0("chr", gene_order_raw$chr[
  match(gene_table$gene, gene_order_raw$gene)])
gene_table$start <- gene_order_raw$start[
  match(gene_table$gene, gene_order_raw$gene)]
gene_table$chr   <- factor(gene_table$chr, levels=CHR_ORDER)
gene_table$genome_idx <- seq_len(nrow(gene_table))

# Get the subclusters for the Tumor group
valid_subclusters <- infercnv_obj@tumor_subclusters$subclusters[["Tumor"]]

# Filter to subclusters with enough cells (uses MIN_CELLS from config)
subcluster_sizes <- sapply(valid_subclusters, length)
valid_subclusters <- valid_subclusters[subcluster_sizes >= MIN_CELLS]

subcluster_cols <- names(valid_subclusters)

# Per-subcluster mean expression (inferCNV scale: neutral = 1)
# Do NOT multiply by 2 — the 0.7–1.4 scale is what inferCNV uses internally
# and what makes the color mapping meaningful
for (s in subcluster_cols) {
  cells <- intersect(names(valid_subclusters[[s]]), colnames(cnv_mat))
  gene_table[[paste0(s, "_raw")]] <- rowMeans(cnv_mat[, cells, drop=FALSE])
}

# Convenience aliases for exactly 2 subclusters
s1_col <- paste0(subcluster_cols[1], "_raw")
s2_col <- paste0(subcluster_cols[2], "_raw")

gene_table$delta <- gene_table[[s2_col]] - gene_table[[s1_col]]

cat(sprintf("  s1 range: %.3f – %.3f\n",
            min(gene_table[[s1_col]]), max(gene_table[[s1_col]])))
cat(sprintf("  s2 range: %.3f – %.3f\n",
            min(gene_table[[s2_col]]), max(gene_table[[s2_col]])))
cat(sprintf("  delta range: %.4f – %.4f\n",
            min(gene_table$delta), max(gene_table$delta)))

# ══════════════════════════════════════════════════════════════════════════════
# PART 9 — Categorise subclonal events
# Logic: direction is determined by delta; gain/loss label by expression vs 1
#   s2 > s1 AND s2 > 1  → s2_gain
#   s2 > s1 AND s1 < 1  → s1_loss  (s2 is higher but s1 is in loss territory)
#   s1 > s2 AND s1 > 1  → s1_gain
#   s1 > s2 AND s2 < 1  → s2_loss
#   |delta| < threshold → neutral
# ══════════════════════════════════════════════════════════════════════════════
DELTA_THRESHOLD <- quantile(abs(gene_table$delta), 0.9, na.rm=TRUE)
cat("Auto delta threshold (75th percentile):", DELTA_THRESHOLD, "\n")

gene_table$category <- case_when(
  abs(gene_table$delta) < DELTA_THRESHOLD                                          ~ "neutral",
  gene_table$delta >  DELTA_THRESHOLD & gene_table[[s2_col]] > 1                  ~ "s2_gain",
  gene_table$delta >  DELTA_THRESHOLD & gene_table[[s1_col]] < 1                  ~ "s1_loss",
  gene_table$delta < -DELTA_THRESHOLD & gene_table[[s1_col]] > 1                  ~ "s1_gain",
  gene_table$delta < -DELTA_THRESHOLD & gene_table[[s2_col]] < 1                  ~ "s2_loss",
  TRUE                                                                              ~ "neutral"
)

cat("  Category distribution:\n")
print(table(gene_table$category))

# ══════════════════════════════════════════════════════════════════════════════
# PART 10 — Chromosome axis stats
# ══════════════════════════════════════════════════════════════════════════════
chr_stats <- gene_table %>%
  group_by(chr) %>%
  summarise(x_mid=mean(genome_idx), x_min=min(genome_idx), .groups="drop") %>%
  arrange(chr)

chr_bounds_x <- chr_stats$x_min[-1] - 0.5

# ══════════════════════════════════════════════════════════════════════════════
# PART 11 — Build plot data
# ══════════════════════════════════════════════════════════════════════════════
expr_tracks <- rbind(
  data.frame(genome_idx = gene_table$genome_idx,
             chr        = gene_table$chr,
             track      = subcluster_cols[2],   # s2 on top
             value      = gene_table[[s2_col]]),
  data.frame(genome_idx = gene_table$genome_idx,
             chr        = gene_table$chr,
             track      = subcluster_cols[1],
             value      = gene_table[[s1_col]])
)
expr_tracks$track <- factor(expr_tracks$track,
                            levels=c(subcluster_cols[2], subcluster_cols[1]))

delta_cat <- data.frame(
  genome_idx = gene_table$genome_idx,
  chr        = gene_table$chr,
  track      = "Subclonal Events",
  category   = gene_table$category
)
delta_cat$track <- factor(delta_cat$track, levels="Subclonal Events")

cat_colors <- c(
  "neutral" = "white",
  "s1_gain" = "#1a7a1a",   # dark green  — s1 gain
  "s1_loss" = "#a8d8a8",   # light green — s1 loss
  "s2_gain" = "#6a0dad",   # dark purple — s2 gain
  "s2_loss" = "#c9a8e8"    # light purple — s2 loss
)

# ══════════════════════════════════════════════════════════════════════════════
# PART 12 — Plot
# ══════════════════════════════════════════════════════════════════════════════
cat("Plotting...\n")

p <- ggplot() +
  # Expression tracks (s2 and s1)
  geom_tile(
    data = expr_tracks,
    aes(x=genome_idx, y=0.5, fill=value), height=1
  ) +
  scale_fill_gradientn(
    colours = c("#2166ac","#92c5de","white","#f4a582","#d6604d"),
    values  = scales::rescale(c(0.7, 0.9, 1.0, 1.1, 1.3)),
    limits  = c(0.7, 1.3),
    oob     = scales::squish,
    name    = "Smoothed\nExpression\n(1=neutral)",
    breaks  = c(0.7, 1.0, 1.3),
    labels  = c("0.7\n(loss)", "1.0\nneutral", "1.3\n(gain)")
  ) +
  new_scale_fill() +
  # Categorical subclonal event track
  geom_tile(
    data = delta_cat,
    aes(x=genome_idx, y=0.5, fill=category), height=1
  ) +
  scale_fill_manual(
    values = cat_colors,
    name   = "Subclonal Event",
    breaks = c("s1_gain","s1_loss","s2_gain","s2_loss","neutral")
    ) +
  facet_grid(track ~ chr, scales="free_x", space="free_x", switch="y") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  labs(
    x = NULL, y = NULL
  ) +
  theme_classic(base_size=11) +
  theme(
    strip.text.x     = element_text(face="bold", size=7),
    strip.text.y     = element_text(face="bold", size=9, angle=0),
    strip.background = element_blank(),
    strip.placement  = "outside",
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    axis.line        = element_blank(),
    panel.border     = element_rect(color="grey50", fill=NA, linewidth=0.5),
    panel.spacing.x  = unit(0.08, "cm"),
    panel.spacing.y  = unit(0.4,  "cm"),
    legend.position  = "right"
  )

save_figure(p, "FigS1E_cnv_heatmap", width = 8, height = 5)

cat("\nDone!\n")
cat("Plot saved to:", file.path(fig_dir, "FigS1E_cnv_heatmap.pdf"), "\n")
