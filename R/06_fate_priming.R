# =============================================================================
# Fate Map & T47D Clone Analysis Pipeline  —  FIGURE 4
# Step 6 of 8
# =============================================================================

source("config.R")
source("R/utils/functions.R")
# Purpose:
#   Characterises transcriptional fate priming in naive T47D breast cancer cells
#   by linking clone-level hybrid fate outcomes to pre-existing gene expression
#   signatures.  Identifies markers that predict whether a cell will enter a
#   senescent or hybrid (EMT-like) state after treatment.
#
# Figures produced & auto-saved to output_dir/figures/Figure_4/:
#   Fig 4B  - Bar chart:  fraction of hybrid-fated clones per CNV cluster
#   Fig 4C  - Euler diagram: overlap between naive-state and fate markers
#   Fig 4D  - Volcano plot: differential expression (senescence vs. hybrid primed)
#   Fig 4E  - Violin plots: fate marker expression split by fate category
#   Fig 4F  - Violin plots: LRP10 expression by hybrid fraction bin (CNV-1 & CNV-2)
#
# Dependencies: Seurat, ggplot2, ggrepel, dplyr, eulerr, scales, patchwork
# =============================================================================


# ── 0. Libraries ──────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(eulerr)
  library(scales)
  library(patchwork)  
})


# ── 1. Parameters & Output Setup ──────────────────────────────────────────────

if (!exists("T47D")) stop("Object `T47D` not found in the session.")

# Figure-specific sub-folder
fig_dir   <- file.path(FIGURES_DIR, "Figure_4")
table_dir <- file.path(OUTPUT_DIR, "supplementary_tables")
dir.create(fig_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

cat("Figures →", fig_dir,   "\n")
cat("Tables  →", table_dir, "\n")



safe_mean_dist <- function(labels, embedding_matrix) {
  temp_df <- tibble(cell = names(labels), lineage = as.character(labels))
  lineage_means <- c()
  
  for (lineage_id in unique(temp_df$lineage)) {
    lineage_cells <- temp_df %>%
      filter(lineage == lineage_id) %>%
      pull(cell)
    
    if (length(lineage_cells) < 2) next
    
    d <- dist(embedding_matrix[lineage_cells, , drop = FALSE])
    lineage_means <- c(lineage_means, mean(as.vector(d)))
  }
  
  if (length(lineage_means) == 0) return(NA_real_)
  
  mean(lineage_means)
}

# =============================================================================
# STAGE 1 — CLONAL COMPACTNESS PERMUTATION TEST
# =============================================================================
clonal_similarity_pipeline <- function(seurat_obj, timepoint, n_permutations = 1000) {
  
  seu_sub <- subset(seurat_obj, time == timepoint & !is.na(fraction_hybrid))
  seu_sub <- FindVariableFeatures(seu_sub, verbose = FALSE)
  seu_sub <- ScaleData(seu_sub, verbose = FALSE)
  seu_sub <- RunPCA(seu_sub, verbose = FALSE)
  
  # BUG FIX 1: was [, 20, drop = FALSE] — only PC20; now takes first 20 PCs
  embeddings <- Embeddings(seu_sub, "pca")[, 1:20, drop = FALSE]
  
  lineage_labels <- seu_sub$lineage_barcode
  names(lineage_labels) <- colnames(seu_sub)
  
  real_mean <- safe_mean_dist(lineage_labels, embeddings)
  
  pseudo_means <- replicate(n_permutations, {
    shuffled        <- sample(lineage_labels)
    names(shuffled) <- names(lineage_labels)
    safe_mean_dist(shuffled, embeddings)
  })
  
  # Lower tail: real clones are compact when real_mean is unusually small
  threshold <- quantile(pseudo_means, 0.05, na.rm = TRUE)
  
  # BUG FIX 2: p-value = fraction of null means as small or smaller than real
  # (tests whether real clones are MORE compact than chance)
  p_val <- mean(pseudo_means <= real_mean, na.rm = TRUE)
  
  perm_summary <- tibble(
    timepoint           = timepoint,
    real_mean           = real_mean,
    null_threshold_5pct = threshold,
    empirical_p_value   = p_val,
    n_permutations      = n_permutations
  )
  
  plot_df_perm <- tibble(mean_dist = pseudo_means)
  
  # BUG FIX 3: xintercept moved outside aes() to avoid NA-rendering issues
  # and ensure legend entries are drawn via dummy aes colour mapping
  p <- ggplot(plot_df_perm, aes(x = mean_dist)) +
    xlim(
      min(real_mean, min(pseudo_means, na.rm = TRUE), na.rm = TRUE) - 0.1,
      max(pseudo_means, na.rm = TRUE) + 0.1
    ) +
    geom_histogram(aes(fill = "Permuted data"), bins = 40,
                   color = "darkgray", alpha = 0.4) +
    geom_vline(aes(color = "Mean clonal distance"),
               xintercept = real_mean,
               linetype = "dashed", linewidth = 1) +
    geom_vline(aes(color = "Threshold of significance (5%)"),
               xintercept = threshold,
               linetype = "dashed", linewidth = 1) +
    # BUG FIX 4: segment y positions computed from data, not Inf, to avoid clipping
    annotate("segment",
             x    = real_mean, xend = threshold,
             y    = max(table(cut(pseudo_means, 40))) * 0.95,
             yend = max(table(cut(pseudo_means, 40))) * 0.95,
             arrow = arrow(ends = "both", length = unit(0.1, "cm"))) +
    annotate("text",
             x     = (real_mean + threshold) / 2,
             y     = Inf,
             label = paste0("p = ", signif(p_val, 3)),
             vjust = 1.5) +
    scale_fill_manual(name   = "", values = c("Permuted data" = "grey70")) +
    scale_color_manual(
      name   = "",
      values = c("Mean clonal distance"        = "blue",
                 "Threshold of significance (5%)" = "red")
    ) +
    labs(
      title = paste0("Clonal compactness — Day ", timepoint),
      x     = "Mean within-lineage distance",
      y     = "Frequency"
    ) +
    theme_classic() +
    theme(legend.position = "right")
  
  list(plot = p, summary = perm_summary)
}

# =============================================================================
# FIGURE 4B — Compactness permutation test
# =============================================================================
res_d0  <- clonal_similarity_pipeline(T47D, "0")
res_d7  <- clonal_similarity_pipeline(T47D, "7")
res_d15 <- clonal_similarity_pipeline(T47D, "15")

p_fig4b    <- res_d0$plot
p_sup_d7   <- res_d7$plot
p_sup_d15  <- res_d15$plot

save_figure(p_fig4b,   "Fig4B_clonal_compactness_day0",  width = 7, height = 5)
save_figure(p_sup_d7,  "FigS3_clonal_compactness_day7",  width = 7, height = 5)
save_figure(p_sup_d15, "FigS3_clonal_compactness_day15", width = 7, height = 5)

perm_summary_all <- bind_rows(
  res_d0$summary,
  res_d7$summary,
  res_d15$summary
)
save_table(perm_summary_all, "Fig4B_permutation_test_summary")


# ── 3. Subset Naive Cells with Known Fate ─────────────────────────────────────
# Restrict to time == 0 (pre-treatment / naive state) so we capture the
# transcriptional landscape *before* any perturbation.
# Cells without a fate assignment are excluded (fate is NA).

naive_cells_T47D <- subset(T47D, time == 0 & !is.na(T47D$fate))


# ── 5. Fate-Predictive Marker Discovery (Volcano Genes) ───────────────────────
# Differential expression: senescence_primed vs. hybrid_primed naive cells.
# ident.2 is specified explicitly; with only two fate categories the result is
# identical numerically, but the explicit call documents the comparison clearly
# and guards against unexpected behaviour if a third fate category is added.
#
# Filters:  min.pct = 0.3    → gene expressed in ≥ 30 % of either group
#           logfc.threshold = 0.1 → loose pre-filter; tighter cut applied below

Idents(naive_cells_T47D) <- "fate"

volcano_plot_genes <- FindMarkers(
  naive_cells_T47D,
  ident.1         = "senescence_primed",
  ident.2         = "hybrid_primed",      
  min.pct         = 0.3,
  logfc.threshold = 0.1
)

fate_markers_all <- filter_de_hits(volcano_plot_genes, min_log2fc = 0.5)

save_table(fate_markers_all, "fate_markers_all", row.names = TRUE)


# ── 6. Naive-State Marker Discovery ───────────────────────────────────────────
# Identifies genes specifically upregulated in the naive (time == 0) state
# relative to all other time points — captures the pre-existing transcriptional
# context that may underlie fate priming.

Idents(T47D) <- "time"

naive_state_markers <- FindMarkers(T47D, ident.1 = 0, min.pct = 0.3) %>%
  filter_de_hits(min_log2fc = 0.5)

save_table(naive_state_markers, "naive_state_markers", row.names = TRUE)


# ── 7. Figure 4C — Euler Diagram: Naive-State ∩ Fate Markers ─────────────────
# Visualises the overlap between genes that define the naive state broadly
# and genes predictive of senescence fate specifically.
# Overlapping genes are candidate "pre-loaded" fate determinants.

fit <- euler(list(
  "Naive state markers" = rownames(naive_state_markers),
  "Fate markers"        = rownames(fate_markers_all)
))

fig4c <- plot(fit,
              quantities = TRUE,
              fills      = c("yellow", "red"),
              alpha      = 0.5)


# ── 8. Figure 4D — Volcano Plot: Senescence vs. Hybrid Primed ─────────────────
# x-axis: avg_log2FC  (positive = higher in senescence_primed)
# y-axis: -log10 adjusted p-value
# Red points cross both the FC and FDR thresholds.
# A curated gene list is labelled for biological interpretation.

volcano_plot_genes$log10_pval  <- -log10(volcano_plot_genes$p_val_adj)
volcano_plot_genes$significant <- volcano_plot_genes$p_val_adj < 0.05 &
  abs(volcano_plot_genes$avg_log2FC) > 0.5

# Genes selected for manual labelling (biologically relevant candidates)
fate_markers_label <- c("MAGEA3", "LRRN3", "ZNF503", "LRP10",
                        "GAS6",   "PREX1", "PTGES",  "CLDN1")

volcano_plot_genes$label <- ifelse(
  rownames(volcano_plot_genes) %in% fate_markers_label,
  rownames(volcano_plot_genes),
  NA
)

fig4d <- ggplot(volcano_plot_genes,
                aes(x = avg_log2FC, y = log10_pval, color = significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("grey", "red")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05),  linetype = "dashed", color = "black") +
  geom_text_repel(aes(label = label),
                  color        = "black",
                  max.overlaps = 1000,
                  size         = 3) +
  labs(x     = "Average log2 Fold Change",
       y     = "-log10 Adjusted p-value",
       color = "Significant") +
  theme_classic()


# ── 9. Bin Naive Cells by Hybrid Fate Propensity ──────────────────────────────
# mean_hybrid is the background hybrid rate at day 15 across all cells —
# used as a reference threshold to split clones into three propensity tiers:
#   Hybrid_low    : clone produced zero hybrid daughters
#   Hybrid_medium : clone hybrid fraction ≤ population mean  (≥ 0, ≤ mean)
#   Hybrid_high   : clone hybrid fraction > population mean
#
# Using <= for the lower boundary and > for the upper boundary ensures
# every non-zero value is captured; fraction_hybrid == mean_hybrid is assigned
# to Hybrid_medium rather than falling through to NA.

mean_hybrid <- sum(T47D$time == 15 & T47D$hybrid == "hybrid") /
  sum(T47D$time == 15)

naive_cells_T47D$fate_bins <- case_when(
  naive_cells_T47D$fraction_hybrid == 0                           ~ "Hybrid_low",
  naive_cells_T47D$fraction_hybrid <= mean_hybrid                 ~ "Hybrid_medium",
  naive_cells_T47D$fraction_hybrid >  mean_hybrid                 ~ "Hybrid_high",
  TRUE                                                            ~ NA_character_
)

naive_cells_T47D$fate_bins <- factor(
  naive_cells_T47D$fate_bins,
  levels = c("Hybrid_high", "Hybrid_medium", "Hybrid_low")
)


# ── 10. Figure 4E — Violin Plots: Key Fate Markers by Fate Category ───────────
# Confirms that the DE markers identified above show consistent expression
# differences when cells are grouped by their experimentally-determined fate.

fig4e <- VlnPlot(naive_cells_T47D,
                 features  = c("LRP10", "LRRN3", "ZNF503", "MAGEA3"),
                 group.by  = "fate",
                 pt.size   = 0,
                 ncol      = 2) &
  scale_fill_manual(values = c("hybrid_primed"      = "purple",
                               "senescence_primed"  = "orange"))


# ── 11. Figure 4F — LRP10 Expression by Hybrid Bin, Stratified by CNV Clone ───
# Splitting by CNV clone (cnv_1 vs cnv_2) controls for copy-number-driven
# expression differences, isolating the fate-priming signal in LRP10.

naive_cnv1_cells <- subset(naive_cells_T47D, cnv_clone == "cnv1")
naive_cnv2_cells <- subset(naive_cells_T47D, cnv_clone == "cnv2")

bin_colors <- c("Hybrid_low"    = "#fee8c8",
                "Hybrid_medium" = "#fdbb84",
                "Hybrid_high"   = "#e34a33")

p_cnv1 <- VlnPlot(naive_cnv1_cells, features = "RIMS1",
                  group.by = "fate_bins", pt.size = 0) +
  scale_fill_manual(values = bin_colors) +
  ggtitle("CNV Clone 1")

p_cnv2 <- VlnPlot(naive_cnv2_cells, features = "RIMS1",
                  group.by = "fate_bins", pt.size = 0) +
  scale_fill_manual(values = bin_colors) +
  ggtitle("CNV Clone 2")

fig4f <- p_cnv1 + p_cnv2   # patchwork composition


# ── 12. Save All Figures ───────────────────────────────────────────────────────

save_euler(fig4c, "fig4c_euler_naive_vs_fate",   width = 5,  height = 5)
save_figure(fig4d, "fig4d_volcano",                 width = 7,  height = 6)
save_figure(fig4e, "fig4e_violin_fate_markers",     width = 8,  height = 8)
save_figure(fig4f, "fig4f_cnv_bins",           width = 10, height = 5)

cat("\nAll Figure 4 panels saved to:", fig_dir, "\n")
cat("Pipeline complete.\n")
