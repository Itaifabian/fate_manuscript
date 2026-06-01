# =============================================================================
# LINEAGE / CLONAL FATE ANALYSIS PIPELINE  —  FIGURE 3
# Step 5 of 8
# =============================================================================

source("config.R")
source("R/utils/functions.R")
# Purpose
#   Integrates lineage-barcode information with a Seurat object, filters
#   low-confidence CNV assignments, quantifies clonal composition over time,
#   measures diversity, evaluates replicate concordance, assigns lineage fate,
#   and tests whether future lineage fate is prefigured in naive-state structure.
#
# Figure mapping
#   - Fig 3A : Stacked bar — clonal proportions by replicate and time
#   - Fig 3B : Bubble plot — top-20 lineage trajectory
#   - Fig 3C : Simpson diversity by replicate and time
#   - Fig 3D : Replicate concordance — lineage size
#   - Fig 3E : Replicate concordance — hybrid fraction
#   - Fig 3F : lineage fate composition at day 15
#   - Fig 3G : Permutation test — naive-state clonal compactness
#   - Sup D7 : Compactness test at day 7
#   - Sup D15: Compactness test at day 15
#
# Required input objects in session
#   - T47D  : Seurat object
#
# Required metadata in T47D
#   - cnv_cluster · hybrid · repeat · sample · time
#
# External input
#   - lineage.db.csv
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(vegan)
  library(tidytext)
  library(RColorBrewer)
  library(Matrix)
  library(xgboost)
  library(ggpmisc)
  library(patchwork)
  library(grid)
  library(stringr)
  library(tibble)
  library(purrr)
  library(dplyr)
})

# =============================================================================
# STAGE 0 — PATHS AND PARAMETERS
# =============================================================================

# ---------- PATHS ------------------------------------------------------------
cnv_file        <- file.path(INFERCNV_DIR, "CNV_table.csv")

# Figure-specific sub-folders
fig_dir   <- file.path(FIGURES_DIR, "Figure_3")
table_dir <- file.path(OUTPUT_DIR, "supplementary_tables")
dir.create(fig_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

cat("Figures →", fig_dir,   "\n")
cat("Tables  →", table_dir, "\n")

# ---------- ANALYSIS PARAMETERS ----------------------------------------------
valid_sample_types    <- c("T47D-lag", "T47D-naive", "T47D-late")
min_lineage_umi_count <- 2
max_barcodes_per_cell <- 1

# XGBoost confidence-filtering
xgb_conf_threshold <- 0.95
xgb_opp_threshold  <- 0.05
xgb_n_boot         <- 50
xgb_sample_frac    <- 0.8
xgb_max_cells_boot <- 5000
xgb_nrounds        <- 150
xgb_max_depth      <- 6
xgb_eta            <- 0.1
xgb_ncores         <- 4

# lineage prevalence
min_lineage_cells <- 2

# Permutation test
n_permutations      <- 1000
pca_dims_for_distance <- 1:20

set.seed(123)


# =============================================================================
# STAGE 1 — VALIDATE INPUTS
# =============================================================================

if (!exists("T47D")) {
  stop("Object `T47D` was not found in the current R session.")
}

required_meta <- c("hybrid", "replicate", "sample", "time")
missing_meta  <- setdiff(required_meta, colnames(T47D@meta.data))
if (length(missing_meta) > 0) {
  stop("T47D is missing required metadata columns: ",
       paste(missing_meta, collapse = ", "))
}

if (!file.exists(LINEAGE_DB_PATH)) {
  stop("Lineage database file not found at: ", LINEAGE_DB_PATH)
}

DefaultAssay(T47D) <- "RNA"
message("Input validation complete.")

# =============================================================================
# STAGE 2 — HELPER FUNCTIONS
# =============================================================================

standardize_cell_name <- function(sample_name, cell_barcode) {
  suffix <- case_when(
    sample_name == "T47D-naive-1" ~ "-1",
    sample_name == "T47D-naive-2" ~ "-2",
    sample_name == "T47D-lag-1"   ~ "-3",
    sample_name == "T47D-lag-2"   ~ "-4",
    sample_name == "T47D-late-1"  ~ "-5",
    sample_name == "T47D-late-2"  ~ "-6",
    TRUE ~ NA_character_
  )
  paste0(cell_barcode, suffix)
}

sample_type_to_time <- function(sample_type) {
  case_when(
    sample_type == "T47D-naive" ~ 0,
    sample_type == "T47D-lag"   ~ 7,
    sample_type == "T47D-late"  ~ 15,
    TRUE ~ NA_real_
  )
}

sample_to_replicate <- function(sample_name) {
  case_when(
    grepl("R1",  sample_name) ~ "Replicate 1",
    grepl("R2",  sample_name) ~ "Replicate 2",
    grepl("-1$", sample_name) ~ "Replicate 1",
    grepl("-2$", sample_name) ~ "Replicate 2",
    TRUE ~ NA_character_
  )
}

# =============================================================================
# STAGE 3 — LOAD AND CLEAN LINEAGE DATABASE
# =============================================================================

lineage_raw <- read_csv(LINEAGE_DB_PATH, show_col_types = FALSE)

required_lineage_cols <- c(
  "sample_type", "sample", "cell_10x_barcode", "standard_barcode", "count"
)
missing_lineage_cols <- setdiff(required_lineage_cols, colnames(lineage_raw))
if (length(missing_lineage_cols) > 0) {
  stop("lineage.db.csv is missing required columns: ",
       paste(missing_lineage_cols, collapse = ", "))
}

lineage_db <- lineage_raw %>%
  group_by(sample_type, sample, cell_10x_barcode, standard_barcode) %>%
  mutate(lineage_barcode = standard_barcode) %>%
  ungroup() %>%
  filter(sample_type %in% valid_sample_types, count > min_lineage_umi_count) %>%
  mutate(cell = standardize_cell_name(sample, cell_10x_barcode)) %>%
  filter(!is.na(cell)) %>%
  group_by(cell) %>%
  mutate(barcode_counts = n()) %>%
  ungroup() %>%
  filter(barcode_counts == max_barcodes_per_cell)

lineage_lookup <- lineage_db %>%
  distinct(cell, .keep_all = TRUE) %>%
  { setNames(.$lineage_barcode, .$cell) }

message("Clean lineage table built: ", nrow(lineage_db), " barcode assignments retained.")

T47D$cell_name <- colnames(T47D)
T47D <- AddMetaData(T47D, metadata = lineage_lookup, col.name = "lineage_barcode")

# =============================================================================
# STAGE 4 — DEFINE lineage-LEVEL CNV IDENTITY
# =============================================================================

cnv_cluster <- read_csv(cnv_file)
lineage_clusters <- lineage_db %>%
  inner_join(
    cnv_cluster)

# Identify misassigned cells and filter by XGBoost confidence
mixed_cells <- lineage_clusters %>%
  group_by(lineage_barcode) %>%
  mutate(fraction_c2 = mean(label == 'Tumor_s2'),
         cnv_cluster = ifelse(fraction_c2 < 0.5, "cnv1", "cnv2"),
    mixed = case_when(
      fraction_c2 < 0.5 & label == "Tumor_s2" ~ TRUE,
      fraction_c2 > 0.5 & label == "Tumor_s1" ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  dplyr::filter(!is.na(fraction_c2))

misassigned_cells <- mixed_cells[mixed_cells$mixed, ]$cell
message("Total misassigned cells: ", length(misassigned_cells))

cnv_lookup <- mixed_cells %>%
  group_by(lineage_barcode) %>%
  summarise(cnv_cluster = names(which.max(table(cnv_cluster))),
            .groups = "drop")

# Set lineage_barcode as rowname for direct lookup
cnv_lookup_df <- cnv_lookup %>%
  column_to_rownames("lineage_barcode")

T47D$cnv_clone <- cnv_lookup_df[T47D$lineage_barcode, "cnv_cluster"]

# ---------- XGBoost confidence filter ----------------------------------------
xgb_copykat_confidence_filtered <- function(
    seu,
    label_col        = "cnv_clone",
    assay            = "RNA",
    layer            = "counts",
    use_hvg          = TRUE,
    ncores           = 4,
    max_cells_boot   = 5000,
    n_boot           = 50,
    sample_frac      = 0.8,
    nrounds          = 150,
    max_depth        = 6,
    eta              = 0.1,
    seed             = 123,
    conf_threshold   = 0.99,
    opp_threshold    = 0.1
) {
  set.seed(seed)
  if (!inherits(seu, "Seurat")) stop("`seu` must be a Seurat object.")
  
  labels_raw      <- seu@meta.data[[label_col]]
  keep_cells_idx  <- !is.na(labels_raw)
  seu_sub         <- subset(seu, cells = colnames(seu)[keep_cells_idx])
  labels_raw      <- seu_sub@meta.data[[label_col]]
  
  label_levels <- sort(unique(as.character(labels_raw)))
  if (length(label_levels) != 2) stop("Confidence filter expects exactly 2 CNV classes.")
  
  y <- ifelse(labels_raw == label_levels[2], 1L, 0L)
  names(y) <- colnames(seu_sub)
  
  mat <- tryCatch(
    GetAssayData(seu_sub, assay = assay, layer = layer),
    error = function(e) GetAssayData(seu_sub, assay = assay, slot = layer)
  )
  
  if (use_hvg) {
    hvg        <- VariableFeatures(seu_sub)
    keep_genes <- if (length(hvg) < 100) Matrix::rowSums(mat) > 0 else rownames(mat) %in% hvg
  } else {
    keep_genes <- Matrix::rowSums(mat) > 0
  }
  
  mat      <- mat[keep_genes, , drop = FALSE]
  X        <- Matrix::t(log1p(mat))
  n_cells  <- nrow(X)
  cell_names <- rownames(X)
  
  prob_mat <- matrix(NA_real_,    nrow = n_cells, ncol = n_boot)
  pred_mat <- matrix(NA_integer_, nrow = n_cells, ncol = n_boot)
  dX_all   <- xgb.DMatrix(data = X)
  
  for (b in seq_len(n_boot)) {
    n_sample <- min(ceiling(sample_frac * n_cells), max_cells_boot)
    idx      <- sample(seq_len(n_cells), size = n_sample, replace = TRUE)
    dtrain   <- xgb.DMatrix(data = X[idx, , drop = FALSE], label = y[idx])
    model    <- xgb.train(
      params  = list(booster = "gbtree", objective = "binary:logistic",
                     max_depth = max_depth, eta = eta, nthread = ncores),
      data    = dtrain, nrounds = nrounds, verbose = 0
    )
    p1               <- predict(model, newdata = dX_all)
    prob_mat[, b]    <- p1
    pred_mat[, b]    <- ifelse(p1 >= 0.5, 1L, 0L)
  }
  
  assigned_prob_mat             <- prob_mat
  assigned_prob_mat[y == 0, ]   <- 1 - assigned_prob_mat[y == 0, ]
  opposite_rate                 <- rowMeans(sweep(pred_mat, 1, y, "!="))
  mean_assigned_conf            <- rowMeans(assigned_prob_mat)
  aggregated_confidence         <- mean_assigned_conf * (1 - opposite_rate)
  
  summary_df <- tibble(
    cell                  = cell_names,
    mean_assigned_conf    = mean_assigned_conf,
    opposite_rate         = opposite_rate,
    aggregated_confidence = aggregated_confidence,
    pass                  = aggregated_confidence > conf_threshold & opposite_rate < opp_threshold
  )
  
  list(
    passing_cells = summary_df %>% filter(pass) %>% pull(cell),
    summary       = summary_df
  )
}

xgb_result <- xgb_copykat_confidence_filtered(
  T47D,
  conf_threshold = xgb_conf_threshold,
  opp_threshold  = xgb_opp_threshold,
  n_boot         = xgb_n_boot,
  sample_frac    = xgb_sample_frac,
  max_cells_boot = xgb_max_cells_boot,
  nrounds        = xgb_nrounds,
  max_depth      = xgb_max_depth,
  eta            = xgb_eta,
  ncores         = xgb_ncores
)

# Join XGBoost summary with misassignment status
conf_df <- xgb_result$summary %>%
  left_join(
    data.frame(
      cell       = colnames(T47D),
      misassigned = colnames(T47D) %in% misassigned_cells
    ),
    by = "cell"
  ) %>%
  mutate(assignment = ifelse(misassigned,
                             "Misassigned\n(cross-subclone)",
                             "Well-assigned\n(same subclone)"))

# Violin plot
p_confidence <- ggplot(conf_df,
       aes(x=assignment, y=mean_assigned_conf, fill=assignment)) +
  geom_boxplot(outlier.size=0.5, width=0.5) +
  scale_fill_manual(values=c(
    "Misassigned\n(cross-subclone)"  = "#d62728",
    "Well-assigned\n(same subclone)" = "#1f77b4")) +
  stat_compare_means(method="wilcox", label="p.signif",
                     comparisons=list(c("Misassigned\n(cross-subclone)",
                                        "Well-assigned\n(same subclone)"))) +
  labs(x=NULL, y="Mean assigned confidence",
       title="Subclone assignment confidence") +
  theme_classic(base_size=12) +
  theme(legend.position="none",
        plot.title=element_text(face="bold"))

save_figure(p_confidence, "Sup_XGBoost_confidence_violin", width=6, height=6)
p_confidence

# =============================================================================
# STAGE 6 — NORMALIZE ABUNDANCE
# =============================================================================

# Prevalent lineages present at both day 0 and day 15
prevalent_lineages <- as.data.frame.matrix(
  table(T47D$lineage_barcode, T47D$time)
) %>%
  rownames_to_column("lineage_barcode") %>%
  filter(`0` > min_lineage_cells, `15` > min_lineage_cells)


sample_totals <- T47D@meta.data %>%
  filter(!is.na(lineage_barcode)) %>%
  group_by(sample) %>%
  summarise(total_barcoded_cells = n(), .groups = "drop")

normalized_lineages <- T47D@meta.data %>%
  filter(!is.na(lineage_barcode)) %>%
  group_by(sample, time, lineage_barcode, lineage_repeat) %>%
  summarise(raw_count = n(), .groups = "drop") %>%
  mutate(
    replicate  = sub(".*_", "", sample),
    time       = factor(as.character(time), levels = c("0", "7", "15")),
    base_lineage = lineage_barcode
  ) %>%
  group_by(lineage_repeat) %>%
  filter(raw_count > min_lineage_cells, n_distinct(time) == 3) %>%
  ungroup() %>%
  left_join(sample_totals, by = "sample") %>%
  mutate(
    freq      = (raw_count / total_barcoded_cells) * 10000,
    replicate = sample_to_replicate(sample)
  )

lineage_order <- normalized_lineages %>%
  filter(time == "0") %>%
  group_by(base_lineage) %>%
  summarise(total_0 = sum(raw_count), .groups = "drop") %>%
  arrange(total_0) %>%
  pull(base_lineage)

all_lineages        <- unique(normalized_lineages$base_lineage)
final_lineage_order <- c(lineage_order, setdiff(all_lineages, lineage_order))

normalized_lineages <- normalized_lineages %>%
  mutate(base_lineage = factor(base_lineage, levels = final_lineage_order))

# lineage colour palette
num_lineages   <- length(final_lineage_order)
base_colors  <- c("#264653", "#2a9d8f", "#8ab17d", "#b8c1ec",
                  "#457b9d", "#a8dadc", "#f1faee", "#6d6875")
pool_size    <- num_lineages * 3
color_pool   <- colorRampPalette(base_colors)(pool_size)
stride_idx   <- seq(1, pool_size, by = 3)[seq_len(num_lineages)]
shuffled_pal <- sample(stride_idx)
better_colors <- setNames(color_pool[shuffled_pal], final_lineage_order)

# =============================================================================
# FIGURE 3A — Stacked bar: clonal proportions across time & replicates
# =============================================================================

p_fig3a <- ggplot(normalized_lineages,
                  aes(x = replicate, y = freq, fill = base_lineage)) +
  geom_col(position = "fill", width = 0.8) +
  scale_fill_manual(values = better_colors) +
  facet_wrap(~time) +
  labs(x = "", y = "Clonal proportion") +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 14),
    text       = element_text(family = "sans", face = "bold", size = 12)
  )

save_figure(p_fig3a, "Fig3A_clonal_proportions", width = 8, height = 5)

# =============================================================================
# FIGURE 3B — Bubble plot: top-20 lineages across time
# =============================================================================

top20_lineages <- normalized_lineages %>%
  filter(time == 0) %>%
  group_by(base_lineage) %>%
  summarise(overall_mean = mean(freq), .groups = "drop") %>%
  slice_max(overall_mean, n = 20) %>%
  pull(base_lineage)

bubble_data <- normalized_lineages %>%
  filter(base_lineage %in% top20_lineages) %>%
  group_by(base_lineage, time) %>%
  summarise(mean_freq = mean(freq), sd_freq = sd(freq), .groups = "drop") %>%
  mutate(base_lineage = factor(base_lineage, levels = rev(top20_lineages)))

p_fig3b <- ggplot(bubble_data,
                  aes(x = time, y = base_lineage, size = mean_freq, colour = base_lineage)) +
  geom_point(alpha = 0.85) +
  scale_size_continuous(name = "Mean freq\n(CP10K)", range = c(2, 16)) +
  scale_colour_manual(values = better_colors[as.character(top20_lineages)], guide = "none") +
  labs(x = "Time (days)", y = "Lineage barcode") +
  theme_classic() +
  theme(
    text           = element_text(family = "sans", face = "bold", size = 12),
    axis.text.y    = element_blank(),
    legend.position = "right"
  )

save_figure(p_fig3b, "Fig3B_top20_bubble", width = 6, height = 6)

# =============================================================================
# FIGURE 3C — Simpson diversity per sample
# =============================================================================

lineage_mat <- T47D@meta.data %>%
  filter(!is.na(lineage_barcode)) %>%
  group_by(sample, time, lineage_barcode) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = lineage_barcode, values_from = n, values_fill = 0)

sample_meta <- lineage_mat %>% dplyr::select(sample, time)

x <- lineage_mat %>% dplyr::select(-sample, -time) %>% as.data.frame()
rownames(x) <- sample_meta$sample
sample_meta$simpson <- vegan::diversity(x, index = "simpson")

sample_meta <- sample_meta %>%
  mutate(
    time      = as.numeric(as.character(time)),
    replicate = sample_to_replicate(sample)
  )

p_fig3c <- ggplot(sample_meta,
                  aes(x = replicate, y = simpson, fill = factor(time))) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  labs(x = "", y = "Simpson diversity", fill = "Time (days)") +
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("0" = "#b2e2e2", "7" = "#2ca25f", "15" = "#006d2c"))

save_figure(p_fig3c, "Fig3C_simpson_diversity", width = 6, height = 5)
save_table(sample_meta, "Fig3C_simpson_diversity")

# =============================================================================
# FIGURES 3D / 3E — Replicate concordance (lineage size & hybrid fraction)
# =============================================================================

hybrid_size_df <- T47D@meta.data %>%
  filter(!is.na(lineage_barcode) & time == 15) %>%
  group_by(sample, lineage_barcode) %>%
  summarise(
    lineage_size       = n(),
    hybrid_fraction  = mean(hybrid == "hybrid", na.rm = TRUE),
    .groups          = "drop"
  ) %>%
  mutate(
    time      = sub("_.*", "", sample),
    replicate = sub(".*_", "", sample)
  ) %>%
  filter(lineage_size > 3)

shared_barcodes <- hybrid_size_df %>%
  filter(time == "15") %>%
  group_by(lineage_barcode) %>%
  summarise(replicates = n_distinct(replicate), .groups = "drop") %>%
  filter(replicates == 2)

wide_df <- hybrid_size_df %>%
  filter(lineage_barcode %in% shared_barcodes$lineage_barcode, time != "0") %>%
  dplyr::select(time, replicate, lineage_barcode, hybrid_fraction, lineage_size) %>%
  pivot_wider(
    names_from  = replicate,
    values_from = c(hybrid_fraction, lineage_size)
  ) %>%
  filter(complete.cases(.))

p_fig3d <- ggplot(wide_df, aes(x = lineage_size_R1, y = lineage_size_R2)) +
  geom_point() +
  geom_smooth(method = "lm", col = "black") +
  stat_fit_glance(
    method = "lm", geom = "text",
    aes(label = paste("R^2 =", signif(..r.squared.., 3),
                      "p =", signif(..p.value.., 3))),
    vjust = 0.1, hjust = 0.1, size = 3, parse = FALSE
  ) +
  theme_classic() +
  labs(x = "lineage size (R1)", y = "lineage size (R2)")

p_fig3e <- ggplot(wide_df, aes(x = hybrid_fraction_R1, y = hybrid_fraction_R2)) +
  geom_point() +
  geom_smooth(method = "lm", col = "black") +
  stat_fit_glance(
    method = "lm", geom = "text",
    aes(label = paste("R^2 =", signif(..r.squared.., 3),
                      "p =", signif(..p.value.., 3))),
    vjust = 0.1, hjust = 0.1, size = 3, parse = FALSE
  ) +
  theme_classic() +
  labs(x = "Hybrid fraction (R1)", y = "Hybrid fraction (R2)")

save_figure(p_fig3d, "Fig3D_replicate_concordance_lineagesize",     width = 5, height = 5)
save_figure(p_fig3e, "Fig3E_replicate_concordance_hybridfraction", width = 5, height = 5)
save_table(wide_df, "Fig3D_3E_replicate_concordance")

# =============================================================================
# STAGE 9 — lineage FATE AT DAY 15
# =============================================================================

fate_map <- T47D@meta.data %>%
  filter(lineage_barcode %in% prevalent_lineages$lineage_barcode, time == 15) %>%
  group_by(lineage_repeat) %>%
  summarise(
    lineage_size       = n(),
    sum_hybrid       = sum(hybrid == "hybrid", na.rm = TRUE),
    fraction_hybrid  = mean(hybrid == "hybrid", na.rm = TRUE),
    .groups          = "drop"
  ) %>%
  filter(lineage_size > min_lineage_cells) %>%
  arrange(fraction_hybrid) %>%
  mutate(
    bins = case_when(
      fraction_hybrid == 0 ~ "senescence_primed",
      fraction_hybrid >  0 ~ "hybrid_primed",
      TRUE ~ NA_character_
    )
  )

plot_df_fate <- fate_map %>%
  mutate(
    fraction_hybrid_order = fraction_hybrid,
    non_hybrid            = 1 - fraction_hybrid
  ) %>%
  pivot_longer(
    cols      = c(fraction_hybrid, non_hybrid),
    names_to  = "state_raw",
    values_to = "fraction"
  ) %>%
  mutate(
    state = recode(state_raw,
                   fraction_hybrid = "Hybrid",
                   non_hybrid      = "Non-hybrid"),
    state = factor(state, levels = c("Hybrid", "Non-hybrid"))
  )

# =============================================================================
# FIGURE 3F — lineage fate composition at day 15
# =============================================================================

p_fig3f <- ggplot(
  plot_df_fate,
  aes(
    x    = reorder_within(lineage_repeat, fraction_hybrid_order, bins),
    y    = fraction,
    fill = state
  )
) +
  geom_col() +
  scale_x_reordered() +
  labs(x = "lineage (lineage barcode)", y = "Fraction of cells", fill = "State") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x        = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  scale_fill_manual(
    values = c("Non-hybrid" = "lightgrey", "Hybrid" = "blue"),
    breaks = c("Hybrid", "Non-hybrid")
  )

save_figure(p_fig3f, "Fig3F_lineage_fate_day15", width = 9, height = 5)
save_table(fate_map, "Fig3F_lineage_fate_map")

# Transfer fate labels back onto cells
indices              <- match(T47D$lineage_repeat, fate_map$lineage_repeat)
T47D$fraction_hybrid <- fate_map$fraction_hybrid[indices]
T47D$fate <- case_when(
  T47D$fraction_hybrid == 0 ~ "senescence_primed",
  T47D$fraction_hybrid >  0 ~ "hybrid_primed",
  TRUE ~ NA_character_
)

# Pull the cluster info using those indices
indices <- match(fate_map$lineage_repeat, T47D$lineage_repeat)
fate_map$clone <- T47D@meta.data$cnv_clone[indices]
library(ggplot2)
fate_map$lineage_count <- 1


p_fig3g <- ggplot(fate_map, aes(x = clone, y = lineage_count, fill = fraction_hybrid)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_gradient(
    low  = "#fee8c8",
    high = "#e34a33",
    name = "Hybrid fraction"
  ) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "",
    y = "Fraction of Clones"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_figure(p_fig3g, "Fig3G_CNV_hybrid_Barplot", width = 9, height = 5)
