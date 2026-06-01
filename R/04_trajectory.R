# ============================================================================
# HYBRID TRAJECTORY ANALYSIS PIPELINE
# Step 4 of 8: trajectory, TF activity, Monocle3 pseudotime (Figure 2)
# Required objects in session:
#   T47D : Seurat object with metadata: time, state_cluster, Phase, S.Score, G2M.Score
# ============================================================================

source("config.R")
source("R/utils/functions.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(viridis)
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
  library(ggpubr)
  library(monocle3)
  library(SeuratWrappers)
  library(purrr)
  library(msigdbr)
  library(patchwork)
  library(decoupleR)
  library(limma)
  library(ggrepel)
})

# -----------------------------------------------------------------------------
# 0. Parameters & output setup
# -----------------------------------------------------------------------------
if (!exists("T47D")) stop("Object `T47D` not found in the session.")

fig_dir   <- file.path(FIGURES_DIR, "Figure_2")
table_dir <- file.path(OUTPUT_DIR, "supplementary_tables")
dir.create(fig_dir,   showWarnings = FALSE, recursive = TRUE)

cat("Figures →", fig_dir,   "\n")
cat("Tables  →", table_dir, "\n")

# Analysis parameters
hybrid_cluster_id  <- "c4"
de_latent_vars     <- c("S.Score", "G2M.Score")
de_test_use        <- "MAST"
hybrid_umap_dims   <- 1:15
monocle_num_dim    <- 15
trajectory_cores   <- 10
n_heatmap_bins     <- 100

# -----------------------------------------------------------------------------
# 1. Validation
# -----------------------------------------------------------------------------
required_meta <- c("time", "state_cluster", "Phase", "S.Score", "G2M.Score")
missing_meta  <- setdiff(required_meta, colnames(T47D@meta.data))
if (length(missing_meta) > 0)
  stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "))

DefaultAssay(T47D) <- "RNA"
message("Validation complete.")

# -----------------------------------------------------------------------------
# 2. Helper functions
# -----------------------------------------------------------------------------

extract_time_numeric <- function(x) {
  out <- suppressWarnings(as.numeric(as.character(x)))
  missing_idx <- is.na(out)
  if (any(missing_idx)) {
    parsed <- stringr::str_extract(as.character(x[missing_idx]), "[0-9]+\\.?[0-9]*")
    out[missing_idx] <- suppressWarnings(as.numeric(parsed))
  }
  out
}

# -----------------------------------------------------------------------------
# 3. Hybrid labels
# -----------------------------------------------------------------------------
T47D$hybrid <- dplyr::case_when(
  as.character(T47D$time) == "0"                         ~ "naive",
  as.character(T47D$state_cluster) == hybrid_cluster_id  ~ "hybrid",
  TRUE                                                   ~ "non-hybrid"
)
T47D$hybrid <- factor(T47D$hybrid, levels = c("naive", "non-hybrid", "hybrid"))
Idents(T47D) <- "hybrid"

hybrid_counts <- as.data.frame(table(T47D$hybrid, T47D$time))
colnames(hybrid_counts) <- c("hybrid", "time", "n_cells")
print(hybrid_counts)
save_table(hybrid_counts, "TableS6_hybrid_cell_counts_by_time")

# -----------------------------------------------------------------------------
# 4. Differential expression: hybrid vs naive / non-hybrid
# -----------------------------------------------------------------------------
de_naive <- FindMarkers(
  T47D, ident.1 = "hybrid", ident.2 = "naive",
  latent.vars = de_latent_vars, test.use = de_test_use, verbose = FALSE
)
de_non_hybrid <- FindMarkers(
  T47D, ident.1 = "hybrid", ident.2 = "non-hybrid",
  latent.vars = de_latent_vars, test.use = de_test_use, verbose = FALSE
)

# Table S7/S8 — full DE results
save_table(rownames_to_column(de_naive,      "gene"), "TableS7_DE_hybrid_vs_naive")
save_table(rownames_to_column(de_non_hybrid, "gene"), "TableS8_DE_hybrid_vs_non_hybrid")

de_naive_sig      <- filter_de_hits(de_naive)
de_non_hybrid_sig <- filter_de_hits(de_non_hybrid)

consistent_up <- intersect(
  rownames(de_naive_sig[de_naive_sig$avg_log2FC > 0, , drop = FALSE]),
  rownames(de_non_hybrid_sig[de_non_hybrid_sig$avg_log2FC > 0, , drop = FALSE])
)
consistent_down <- intersect(
  rownames(de_naive_sig[de_naive_sig$avg_log2FC < 0, , drop = FALSE]),
  rownames(de_non_hybrid_sig[de_non_hybrid_sig$avg_log2FC < 0, , drop = FALSE])
)

message("Consistently upregulated:   ", length(consistent_up))
message("Consistently downregulated: ", length(consistent_down))

# Table S9 — consistent DE gene sets
save_table(
  data.frame(
    gene      = c(consistent_up, consistent_down),
    direction = c(rep("up", length(consistent_up)), rep("down", length(consistent_down)))
  ),
  "TableS9_consistent_hybrid_DE_genes"
)

# -----------------------------------------------------------------------------
# 5. Figure 2A — hybrid signature heatmap
# -----------------------------------------------------------------------------
all_genes <- unique(c(consistent_up, consistent_down))
if (length(all_genes) == 0) stop("No consistent DE genes detected.")

hybrid_avg <- AverageExpression(
  T47D, features = all_genes, group.by = "hybrid", slot = "data", verbose = FALSE
)
hybrid_mat        <- hybrid_avg$RNA
hybrid_mat_scaled <- t(scale(t(hybrid_mat)))
hybrid_mat_scaled <- hybrid_mat_scaled[complete.cases(hybrid_mat_scaled), , drop = FALSE]

ph_2a <- pheatmap(
  hybrid_mat_scaled,
  color          = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  breaks         = seq(-1.5, 1.5, length.out = 101),
  cluster_rows   = TRUE,
  cluster_cols   = FALSE,
  main           = "",
  border_color   = NA,
  fontsize_row   = 8,
  cellwidth      = 50,
  treeheight_row = 0,
  silent         = TRUE
)
save_pheatmap(ph_2a, "Fig2A_hybrid_signature_heatmap", width = 6, height = 10)

# -----------------------------------------------------------------------------
# 5b. Figure 2B — violin plots for cell-cycle arrest genes
# -----------------------------------------------------------------------------
p_fig2b <- VlnPlot(
  T47D,
  features = c("CDKN1A", "CDKN2B"),
  pt.size  = 0,
  group.by = "hybrid",
  cols     = c("grey", "#bdbdbd", "#636363")
)
save_figure(p_fig2b, "Fig2B_CDKN_vlnplot", width = 6, height = 5)

# -----------------------------------------------------------------------------
# 6. TF activity (decoupleR)
# -----------------------------------------------------------------------------
# Step 1 — Pseudobulk RNA (use counts + sum, then CPM normalize)
T47D$hybrid_replicate <- paste0(T47D$time, "_", T47D$hybrid, "_", T47D$replicate)
DefaultAssay(T47D) <- "RNA"

pb_rna <- AggregateExpression(
  T47D,
  group.by = "hybrid_replicate",
  assays = "RNA",
  slot = "counts",          # raw counts, not log-normalized
  return.seurat = FALSE,
  fun = "sum"               # sum counts, normalize after
)$RNA %>% as.matrix()

# CPM normalize
pb_rna_cpm <- sweep(pb_rna, 2, colSums(pb_rna) / 1e6, "/")

# Step 2 — Run ULM on pseudobulk CPM profiles
net <- read.csv(TF_NETWORK) %>%
  mutate(source = toupper(source), target = toupper(target))

tf_pb <- decoupleR::run_ulm(
  mat     = pb_rna_cpm,
  net     = net,
  .source = "source",
  .target = "target",
  .mor    = "mor",
  minsize = 15
)

# Step 3 — Pivot to matrix (TFs x 10 groups)
tf_pb_mat <- tf_pb %>%
  dplyr::filter(statistic == "ulm") %>%
  select(source, condition, score) %>%
  pivot_wider(names_from = condition, values_from = score) %>%
  tibble::column_to_rownames("source") %>%
  as.matrix()


meta <- data.frame(sample = colnames(tf_pb_mat)) %>%
  mutate(
    is_hybrid = factor(
      ifelse(grepl("hybrid", sample) & !grepl("non", sample), "hybrid", "other"),
      levels = c("other", "hybrid")  # "other" as reference → positive logFC = higher in hybrid
    )
  )

stopifnot(all(meta$sample == colnames(tf_pb_mat)))

design <- model.matrix(~ is_hybrid, data = meta)
fit    <- lmFit(tf_pb_mat, design)
fit    <- eBayes(fit)

# Get full results (no filter)
all_results <- topTable(fit, coef = "is_hybridhybrid", n = Inf) %>%  # fixed coef name
  rownames_to_column("TF") %>%
  mutate(
    sig   = adj.P.Val < 0.05,
    label = ifelse(sig, TF, NA)
  )

p_fig2d <- ggplot(all_results, aes(x = logFC, y = -log10(adj.P.Val), color = sig, label = label)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_text_repel(size = 3, max.overlaps = 20, na.rm = TRUE, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick")) +
  labs(
    x        = "logFC (hybrid - other)",
    y        = "-log10(p-value)",
    color    = "FDR < 0.05"
  ) +
  theme_bw() +
  theme_classic()

save_figure(p_fig2d, "Fig2D_TF_vlocanoplot", width = 6, height = 5)

# -----------------------------------------------------------------------------
# 7. Hybrid-only analysis
# Isolating the hybrid cells in order to study the transcriptional dynamic of hybrid cells
# -----------------------------------------------------------------------------
Hybrid_T47D <- subset(T47D, subset = hybrid == "hybrid")
Hybrid_T47D <- FindVariableFeatures(Hybrid_T47D)
Hybrid_T47D <- ScaleData(Hybrid_T47D, vars.to.regress = c("Phase", "time"))
Hybrid_T47D <- RunPCA(Hybrid_T47D)
Hybrid_T47D <- RunUMAP(Hybrid_T47D, dims = 1:10)


# -----------------------------------------------------------------------------
# 7a. Monocle3 pseudotime analysis
# -----------------------------------------------------------------------------
cds <- as.cell_data_set(Hybrid_T47D)
cds <- detect_genes(cds)
cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds, method = "PCA",
                      use_genes = VariableFeatures(Hybrid_T47D),
                      num_dim   = monocle_num_dim)
cds <- reduce_dimension(cds)

# Use Seurat UMAP (cell-cycle + time regressed) as the trajectory space
reducedDims(cds)$UMAP <- Embeddings(Hybrid_T47D, "umap")

cds <- cluster_cells(cds, k = 10)
cds <- learn_graph(cds)
cds <- order_cells(cds)   # interactive root selection

p_fig2d <- plot_cells(
  cds,
  color_cells_by          = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_leaves            = FALSE,
  label_branch_points     = FALSE,
  show_trajectory_graph   = FALSE,
  cell_size               = 1
)
save_figure(p_fig2d, "Fig2D_pseudotime_UMAP", width = 6, height = 5)

# -----------------------------------------------------------------------------
# 8. Trajectory-correlated genes (Moran's I)
# -----------------------------------------------------------------------------
deg_graph <- graph_test(cds, neighbor_graph = "principal_graph", cores = trajectory_cores)
deg_graph$gene <- rownames(deg_graph)

traj_graph_genes <- deg_graph %>%
  dplyr::filter(q_value < 0.05, morans_I >= 0.2) %>%
 dplyr::arrange(desc(morans_I))

save_table(traj_graph_genes, "TableS11_trajectory_Morans_I_genes")

genes_to_model <- traj_graph_genes$gene
genes_to_model <- genes_to_model[!grepl("^RPL|^RPS|^MT", genes_to_model)]

test_cds          <- cds[genes_to_model, ]
valid_cells       <- which(is.finite(pseudotime(test_cds)))
test_cds_filtered <- test_cds[, valid_cells]

test_model <- fit_models(
  test_cds_filtered,
  model_formula_str = "~ splines::ns(pseudotime, df=3)",
  clean_model       = TRUE
)

t_stat_df <- coefficient_table(test_model) %>%
  dplyr::filter(term %in% c(
    "splines::ns(pseudotime, df = 3)1",
    "splines::ns(pseudotime, df = 3)2",
    "splines::ns(pseudotime, df = 3)3"
  )) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(
    net_t_value   = sum(test_val),
    model_q_value = unique(q_value)[1],
    .groups       = "drop"
  ) %>%
  dplyr::filter(model_q_value < 0.05 & abs(net_t_value) > 10)

save_table(t_stat_df, "TableS12_pseudotime_spline_gene_stats")

top_up   <- dplyr::filter(t_stat_df, net_t_value > 0) %>% dplyr::arrange(desc(net_t_value)) %>% head(20) %>% pull(gene_id)
top_down <- dplyr::filter(t_stat_df, net_t_value < 0) %>% dplyr::arrange(net_t_value)       %>% head(20) %>% pull(gene_id)
genes_to_plot <- c(top_down, top_up)

# -----------------------------------------------------------------------------
# 9. Figure 2E — pseudotime heatmap
# -----------------------------------------------------------------------------
pt.matrix  <- pseudotime(cds)
cell_order <- names(pt.matrix)[order(pt.matrix)]
pt.matrix  <- pt.matrix[cell_order]
exprs_mat  <- normalized_counts(cds)[, cell_order]

heatmap_matrix <- matrix(NA, nrow = length(genes_to_plot), ncol = n_heatmap_bins,
                         dimnames = list(genes_to_plot, NULL))

for (i in seq_along(genes_to_plot)) {
  gene       <- genes_to_plot[i]
  smooth_fit <- smooth.spline(pt.matrix, exprs_mat[gene, ], df = 3)
  heatmap_matrix[i, ] <- predict(
    smooth_fit,
    seq(min(pt.matrix), max(pt.matrix), length.out = n_heatmap_bins)
  )$y
}

heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))
heatmap_matrix_scaled[heatmap_matrix_scaled >  2] <-  2
heatmap_matrix_scaled[heatmap_matrix_scaled < -2] <- -2

pseudotime_range <- seq(min(pt.matrix), max(pt.matrix), length.out = n_heatmap_bins)
colnames(heatmap_matrix_scaled) <- pseudotime_range

col_annotation    <- data.frame(Pseudotime = pseudotime_range)
rownames(col_annotation) <- colnames(heatmap_matrix_scaled)
annotation_colors <- list(Pseudotime = colorRampPalette(c("purple", "yellow"))(100))

ph_2e <- pheatmap(
  heatmap_matrix_scaled,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  show_rownames     = TRUE,
  show_colnames     = FALSE,
  annotation_col    = col_annotation,
  annotation_colors = annotation_colors,
  main              = "",
  color             = colorRampPalette(c("blue", "lightblue", "yellow", "red"))(100),
  treeheight_row    = 0,
  silent            = TRUE
)
save_pheatmap(ph_2e, "Fig2E_pseudotime_heatmap", width = 10, height = 8)

# -----------------------------------------------------------------------------
# 10. Figure 2F — pseudotime vs signature score correlations
# -----------------------------------------------------------------------------
Hybrid_T47D$pseudotime <- pt.matrix[colnames(Hybrid_T47D)]
  
m_df <- msigdbr::msigdbr()
  
  
Hybrid_T47D <- AddModuleScore(
    Hybrid_T47D,
    features = list(
      m_df[m_df$gs_name == "FRIDMAN_SENESCENCE_UP",              ]$gene_symbol,
      m_df[m_df$gs_name == "FISCHER_DREAM_TARGETS",              ]$gene_symbol,
      m_df[m_df$gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION", ]$gene_symbol
    )
  )
Hybrid_T47D$senescence_score <- Hybrid_T47D$Cluster1
Hybrid_T47D$prolif_score     <- Hybrid_T47D$Cluster2
Hybrid_T47D$OXPHOS_score     <- Hybrid_T47D$Cluster3
  
plot_data <- Hybrid_T47D@meta.data
  
# Bin before correlation to avoid p-value inflation from cell-level non-independence
binned <- plot_data %>%
    mutate(bin = cut(pseudotime, breaks = 50)) %>%
    group_by(bin) %>%
    dplyr::summarise(
      pseudotime       = mean(pseudotime),
      prolif_score     = mean(prolif_score),
      senescence_score = mean(senescence_score),
      oxphos_score     = mean(OXPHOS_score),
      n_cells          = n(),
      .groups          = "drop"
    )
  
res_pro    <- cor.test(binned$pseudotime, binned$prolif_score)
res_sen <- cor.test(binned$pseudotime, binned$senescence_score, weight = binned$n_cells)
res_oxphos <- cor.test(binned$pseudotime, binned$oxphos_score)

message(sprintf("Proliferation vs pseudotime:  r = %.3f, p = %.3e", res_pro$estimate,    res_pro$p.value))
message(sprintf("Senescence vs pseudotime:     r = %.3f, p = %.3e", res_sen$estimate,    res_sen$p.value))
message(sprintf("OXPHOS vs pseudotime:         r = %.3f, p = %.3e", res_oxphos$estimate, res_oxphos$p.value))

# Table S13 — correlation statistics
save_table(
  data.frame(
    signature   = c("Proliferation", "Senescence", "OXPHOS"),
    r           = c(res_pro$estimate, res_sen$estimate, res_oxphos$estimate),
    p_value     = c(res_pro$p.value,  res_sen$p.value,  res_oxphos$p.value),
    n_bins      = nrow(binned)
  ),
  "TableS13_pseudotime_signature_correlations"
)

plot_score_pseudotime <- function(df, score_column, title) {
  ggplot(df, aes_string(x = "pseudotime", y = score_column, color = "pseudotime")) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_smooth(method = "lm", color = "black", linewidth = 1, se = TRUE) +
    scale_color_viridis_c(option = "plasma") +
    theme_classic() +
    labs(title = title, y = "Score", x = "Pseudotime")
}

p_fig2g <- (plot_score_pseudotime(plot_data, "prolif_score",     "Proliferation Score") /
              plot_score_pseudotime(plot_data, "senescence_score", "Senescence Score")    /
              plot_score_pseudotime(plot_data, "OXPHOS_score",     "OXPHOS Score"))

save_figure(p_fig2g, "Fig2G_pseudotime_signature_scores", width = 5, height = 12)

cat("\n=== Figure 2 analysis complete ===\n")
cat("Figures →", fig_dir,   "\n")
cat("Tables  →", table_dir, "\n")

