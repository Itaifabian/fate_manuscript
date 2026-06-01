# ==============================================================================
# MASTER HYBRID DETECTION PIPELINE
# Datasets: GSE131984, GSE134836, GSE146544, GSE168668,
#           GSE234181, GSE280381, GSE297393
# Output:   all_hybrid_proportions.csv + Fig2F bubble plot (PDF + PNG)
#
# Convention: run_qc_and_doublet_removal() accepts either
#   (a) a single Seurat object  → split internally by split_by
#   (b) a named list of Seurat objects already split per sample
#   percent.mt must be computed on each object BEFORE passing to the wrapper.
#
# Data paths: set HYBRID_DATASETS_DIR in config.R.
#   Each dataset lives in a subdirectory of that root; subfolder names are
#   listed per-section below and may need adjustment to match your layout.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(Matrix)
  library(stringr)
  library(mclust)
  library(msigdbr)
  library(purrr)
  library(tibble)
})

# ==============================================================================
# 0.  Source shared utilities
#     Provides: run_qc_and_doublet_removal(), save_figure(), detect_hybrids()
# ==============================================================================

source("R/utils/functions.R")
source("config.R")

fig_dir <- file.path(FIGURES_DIR, "Figure_2")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(42)

# ==============================================================================
# 1.  Gene sets  (shared by all datasets)
# ==============================================================================

message("Loading MSigDB gene sets...")
m_df         <- msigdbr::msigdbr()
prolif_genes <- m_df[m_df$gs_name == "HALLMARK_G2M_CHECKPOINT", ]$gene_symbol
senesc_genes <- m_df[m_df$gs_name %in% c("FRIDMAN_SENESCENCE_UP",
                                         "KEGG_LYSOSOME"),      ]$gene_symbol

# ==============================================================================
# 2.  Helper: tidy proportion table for one Seurat object
# ==============================================================================

detect_hybrids <- function(seu,
                           prolif_genes,
                           senesc_genes,
                           condition_col = "condition",
                           ctr_label     = "ctr",
                           tis_label     = "tis",
                           seed          = 42) {

  library(dplyr); library(ggplot2); library(patchwork); library(mclust)

  set.seed(seed)
  seu <- AddModuleScore(seu,
                        features = list(prolif_genes, senesc_genes),
                        name     = "State",
                        seed     = seed)

  ctr_idx <- seu[[condition_col, drop = TRUE]] == ctr_label
  tis_idx <- seu[[condition_col, drop = TRUE]] == tis_label

  # ── Tier 1: GMM (Mclust, 2-component) ────────────────────────────────────────
  find_gmm_boundary <- function(scores, label = "") {

    fit <- tryCatch(
      Mclust(scores, G = 2, verbose = FALSE),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      cat(sprintf("GMM (%s): Mclust failed to converge.\n", label))
      return(NULL)
    }

    means <- sort(fit$parameters$mean)
    sds   <- sqrt(fit$parameters$variance$sigmasq)
    if (length(sds) == 1) sds <- rep(sds, 2)
    sds   <- sds[order(fit$parameters$mean)]

    overlap <- exp(
      -0.25 * log(0.25 * (sds[1]^2/sds[2]^2 + sds[2]^2/sds[1]^2 + 2)) +
        -0.25 * (means[2] - means[1])^2 / (sds[1]^2 + sds[2]^2)
    )

    cat(sprintf("GMM (%s): means=[%.3f, %.3f]  sds=[%.3f, %.3f]  Bhattacharyya overlap=%.3f\n",
                label, means[1], means[2], sds[1], sds[2], overlap))

    if (overlap > 0.3) {
      cat(sprintf("GMM (%s): components overlap too much (%.3f > 0.3) — rejecting.\n",
                  label, overlap))
      return(NULL)
    }

    x_seq    <- seq(min(scores), max(scores), length.out = 5000)
    post     <- predict(fit, newdata = x_seq)$z
    boundary <- x_seq[which.min(abs(post[, 1] - post[, 2]))]

    cat(sprintf("GMM (%s): boundary = %.4f\n", label, boundary))
    list(boundary  = boundary,
         fit       = fit,
         method    = "GMM",
         label     = label,
         means     = means,
         sds       = sds,
         scores    = scores)
  }

  # ── Tier 2: KDE valley ────────────────────────────────────────────────────────
  find_kde_valley <- function(scores, label = "") {

    kde    <- density(scores, bw = "SJ", adjust = 0.5, n = 4000)
    window <- max(10L, min(100L, as.integer(0.025 * length(kde$x))))

    is_max <- sapply(seq(window + 1, length(kde$y) - window), function(i) {
      kde$y[i] == max(kde$y[(i - window):(i + window)])
    })
    is_max    <- c(rep(FALSE, window), is_max, rep(FALSE, window))
    sig_peaks <- which(is_max & kde$y > 0.10 * max(kde$y))

    cat(sprintf("KDE (%s): %d peak(s) at x = %s\n",
                label, length(sig_peaks),
                paste(sprintf("%.3f", kde$x[sig_peaks]), collapse = ", ")))

    if (length(sig_peaks) < 2) return(NULL)

    peak_y     <- kde$y[sig_peaks]
    top2       <- sig_peaks[order(peak_y, decreasing = TRUE)[1:2]]
    top2_x     <- sort(kde$x[top2])
    valley_idx <- which(kde$x > top2_x[1] & kde$x < top2_x[2])
    boundary   <- kde$x[valley_idx[which.min(kde$y[valley_idx])]]

    cat(sprintf("KDE (%s): boundary = %.4f\n", label, boundary))
    list(boundary = boundary,
         kde      = kde,
         method   = "KDE",
         label    = label,
         peak_x   = top2_x,
         peak_y   = kde$y[top2[order(kde$x[top2])]],
         scores   = scores)
  }

  # ── Dispatch: Tier 1 → Tier 2a → Tier 2b ─────────────────────────────────────
  cat("\n── Proliferation threshold ─────────────────────────────────────────────\n")

  p_out <- find_gmm_boundary(seu$State1[ctr_idx], "CTR-GMM")

  if (is.null(p_out)) {
    cat("Tier 1 failed → trying KDE on CTR only...\n")
    p_out <- find_kde_valley(seu$State1[ctr_idx], "CTR-KDE")
  }

  if (is.null(p_out)) {
    cat("Tier 2 failed → trying KDE on all cells combined...\n")
    p_out <- find_kde_valley(seu$State1, "All-KDE")
  }

  if (is.null(p_out)) {
    stop(paste(
      "All tiers failed — no bimodal proliferation signal detected",
      "(CTR-GMM, CTR-KDE, All-GMM).",
      "Inspect: ggplot(seu@meta.data, aes(x=State1, colour=condition)) + geom_density()",
      "Check gene coverage: table(prolif_genes %in% rownames(seu))"
    ))
  }

  prolif_boundary <- p_out$boundary
  method_used     <- p_out$method
  tier_used       <- p_out$label

  cat(sprintf("\n✓ Method: %s (%s)  |  Boundary: %.4f\n",
              method_used, tier_used, prolif_boundary))
  cat(sprintf("  CTR prolif: %.1f%%  Tis prolif: %.1f%%\n",
              mean(seu$State1[ctr_idx] >= prolif_boundary) * 100,
              mean(seu$State1[tis_idx] >= prolif_boundary) * 100))

  # ── Senescence threshold — CTR 95th percentile ───────────────────────────────
  cat("\n── Senescence threshold ────────────────────────────────────────────────\n")

  senesc_boundary <- as.numeric(quantile(seu$State2[ctr_idx], 0.95))

  cat(sprintf("CTR 95th pct: %.4f\n", senesc_boundary))
  cat(sprintf("  CTR senesc: %.1f%%  Tis senesc: %.1f%%\n",
              mean(seu$State2[ctr_idx] >= senesc_boundary) * 100,
              mean(seu$State2[tis_idx] >= senesc_boundary) * 100))

  # ── Classify ──────────────────────────────────────────────────────────────────
  seu$state_cluster <- case_when(
    seu$State1 >= prolif_boundary & seu$State2 >= senesc_boundary ~ "Hybrid",
    seu$State1 >= prolif_boundary & seu$State2 <  senesc_boundary ~ "Proliferative",
    seu$State1 <  prolif_boundary & seu$State2 >= senesc_boundary ~ "Senescent",
    TRUE                                                            ~ "Low/Quiescent"
  )
  seu$threshold_method <- tier_used

  cat("\nState distribution:\n")
  print(table(seu$state_cluster, seu[[condition_col, drop = TRUE]]))

  # ── Diagnostic plots ──────────────────────────────────────────────────────────
  if (method_used == "GMM") {

    x_seq  <- seq(min(p_out$scores), max(p_out$scores), length.out = 500)
    props  <- p_out$fit$parameters$pro
    props  <- props[order(p_out$fit$parameters$mean)]
    comp1  <- dnorm(x_seq, p_out$means[1], p_out$sds[1]) * props[1]
    comp2  <- dnorm(x_seq, p_out$means[2], p_out$sds[2]) * props[2]
    gmm_df <- data.frame(
      x         = rep(x_seq, 2),
      y         = c(comp1, comp2),
      Component = rep(c("Low", "High"), each = 500)
    )

    p_prolif <- ggplot() +
      geom_histogram(data = data.frame(x = p_out$scores),
                     aes(x = x, y = after_stat(density)),
                     bins = 60, fill = "grey80", color = "grey60") +
      geom_line(data = gmm_df, aes(x = x, y = y, color = Component),
                linewidth = 1.2) +
      geom_vline(xintercept = prolif_boundary,
                 linetype = "dashed", linewidth = 1, color = "#D85A30") +
      annotate("text", x = prolif_boundary + 0.01, y = Inf,
               label = sprintf("GMM boundary = %.3f", prolif_boundary),
               hjust = 0, vjust = 1.5, size = 3.2) +
      scale_color_manual(values = c("Low" = "#4E9BCA", "High" = "#1D9E75")) +
      labs(title    = sprintf("Proliferation — GMM (%s)", tier_used),
           subtitle = sprintf("CTR=%.1f%%  Tis=%.1f%%",
                              mean(seu$State1[ctr_idx] >= prolif_boundary) * 100,
                              mean(seu$State1[tis_idx] >= prolif_boundary) * 100),
           x = "Proliferation score", y = "Density") +
      theme_classic(base_size = 12)

  } else {

    kde_df   <- data.frame(x = p_out$kde$x, y = p_out$kde$y)
    peaks_df <- data.frame(x = p_out$peak_x, y = p_out$peak_y)

    p_prolif <- ggplot() +
      geom_histogram(data = data.frame(x = p_out$scores),
                     aes(x = x, y = after_stat(density)),
                     bins = 60, fill = "grey80", color = "grey60") +
      geom_line(data = kde_df, aes(x = x, y = y),
                color = "#4E9BCA", linewidth = 1.2) +
      geom_point(data = peaks_df, aes(x = x, y = y),
                 shape = 8, size = 4, color = "#1D9E75", stroke = 1.5) +
      geom_vline(xintercept = prolif_boundary,
                 linetype = "dashed", linewidth = 1, color = "#D85A30") +
      annotate("text", x = prolif_boundary + 0.01, y = Inf,
               label = sprintf("KDE valley = %.3f", prolif_boundary),
               hjust = 0, vjust = 1.5, size = 3.2) +
      labs(title    = sprintf("Proliferation — KDE (%s)", tier_used),
           subtitle = sprintf("CTR=%.1f%%  Tis=%.1f%%",
                              mean(seu$State1[ctr_idx] >= prolif_boundary) * 100,
                              mean(seu$State1[tis_idx] >= prolif_boundary) * 100),
           x = "Proliferation score", y = "Density") +
      theme_classic(base_size = 12)

  }

  p_senesc <- ggplot(data.frame(x = seu$State2[ctr_idx]), aes(x = x)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 60, fill = "grey80", color = "grey60") +
    geom_vline(xintercept = senesc_boundary,
               linetype = "dashed", linewidth = 1, color = "#D85A30") +
    annotate("rect",
             xmin = senesc_boundary, xmax = Inf,
             ymin = 0,              ymax = Inf,
             fill = "#D85A30", alpha = 0.08) +
    annotate("text", x = senesc_boundary + 0.002, y = Inf,
             label = sprintf("CTR 95th pct = %.3f", senesc_boundary),
             hjust = 0, vjust = 1.5, size = 3.2) +
    labs(title    = "Senescence — CTR 95th percentile",
         subtitle = sprintf("CTR=%.1f%%  Tis=%.1f%%",
                            mean(seu$State2[ctr_idx] >= senesc_boundary) * 100,
                            mean(seu$State2[tis_idx] >= senesc_boundary) * 100),
         x = "Senescence score (CTR only)", y = "Density") +
    theme_classic(base_size = 12)

  p_scatter <- ggplot(seu@meta.data,
                      aes(x = State1, y = State2,
                          color = .data[[condition_col]])) +
    geom_point(size = 0.3, alpha = 0.4) +
    geom_vline(xintercept = prolif_boundary,
               linetype = "dashed", linewidth = 0.8) +
    geom_hline(yintercept = senesc_boundary,
               linetype = "dashed", linewidth = 0.8) +
    scale_color_manual(values = setNames(c("#4E9BCA", "#D85A30"),
                                         c(ctr_label, tis_label))) +
    labs(title    = sprintf("Hybrid detection  [%s]", tier_used),
         subtitle = sprintf("Prolif = %.3f (%s) | Senesc = %.3f (CTR 95th pct)",
                            prolif_boundary, method_used, senesc_boundary),
         x = "Proliferation score", y = "Senescence score") +
    theme_classic(base_size = 12)

  print((p_prolif | p_senesc) / p_scatter)

  return(seu)
}


compute_proportions <- function(seurat_obj,
                                dataset,
                                cell_line,
                                condition_col = "condition",
                                ctr_label     = "ctr",
                                tis_label     = "tis") {

  meta <- seurat_obj@meta.data

  if (!condition_col %in% colnames(meta)) {
    warning(sprintf("[%s] Column '%s' not found — skipping.", dataset, condition_col))
    return(NULL)
  }
  if (!"state_cluster" %in% colnames(meta)) {
    warning(sprintf("[%s] 'state_cluster' missing — was detect_hybrids() run?", dataset))
    return(NULL)
  }

  meta %>%
    dplyr::filter(.data[[condition_col]] %in% c(ctr_label, tis_label)) %>%
    dplyr::group_by(condition_val = .data[[condition_col]]) %>%
    dplyr::summarise(
      n_cells    = n(),
      n_hybrid   = sum(state_cluster == "Hybrid"),
      pct_hybrid = 100 * mean(state_cluster == "Hybrid"),
      .groups    = "drop"
    ) %>%
    dplyr::mutate(
      dataset   = dataset,
      cell_line = cell_line,
      group     = if_else(condition_val == ctr_label, "Control", "Treated"),
      y_pos     = if_else(group == "Control", 1L, 2L)
    )
}

# Accumulator
all_proportions <- list()

# ==============================================================================
# 3.  GSE131984 — SUM159 TNBC  (DMSO vs JQ1+Pal)
# ==============================================================================
{
  message("\n\n══════════════════════════════════════════════════════════")
  message("  GSE131984 — SUM159 TNBC")
  message("══════════════════════════════════════════════════════════")

  data_dir <- file.path(HYBRID_DATASETS_DIR, "GSE131984")

  message("  Loading metadata...")
  meta_df           <- fread(file.path(data_dir, "GSE131984_scRNA_metadata.txt.gz"),
                             data.table = FALSE)
  rownames(meta_df) <- meta_df[, 1]
  if (!"sample" %in% colnames(meta_df) && "treatment" %in% colnames(meta_df))
    meta_df$sample <- meta_df$treatment

  message("  Loading count matrix (may take a minute)...")
  dense_counts           <- fread(file.path(data_dir, "GSE131984_scRNA_counts.txt.gz"),
                                  data.table = FALSE)
  rownames(dense_counts) <- dense_counts[, 1]
  dense_counts           <- dense_counts[, -1]
  sparse_counts          <- Matrix(as.matrix(dense_counts), sparse = TRUE)
  rm(dense_counts); gc()

  common_cells  <- intersect(colnames(sparse_counts), rownames(meta_df))
  sparse_counts <- sparse_counts[, common_cells]
  meta_df       <- meta_df[common_cells, ]

  seurat_list <- list()
  for (batch in unique(meta_df$orig.ident)) {
    message("  Building object for batch: ", batch)
    bc   <- rownames(meta_df)[meta_df$orig.ident == batch]
    sobj <- CreateSeuratObject(counts    = sparse_counts[, bc],
                               project   = batch,
                               meta.data = meta_df[bc, ])
    sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
    seurat_list[[batch]] <- sobj
    rm(sobj); gc()
  }
  rm(sparse_counts, meta_df); gc()

  seurat_list <- run_qc_and_doublet_removal(seurat_list)

  tnbc_seurat <- merge(seurat_list[[1]], y = seurat_list[-1],
                       add.cell.ids = names(seurat_list))
  tnbc_seurat <- JoinLayers(tnbc_seurat)
  rm(seurat_list); gc()

  tnbc_seurat           <- subset(tnbc_seurat, orig.ident %in% c("DMSO", "JQ1+Pal"))
  tnbc_seurat$condition <- ifelse(tnbc_seurat$orig.ident == "DMSO", "ctr", "tis")

  tnbc_seurat <- NormalizeData(tnbc_seurat, verbose = FALSE)
  tnbc_seurat <- detect_hybrids(tnbc_seurat, prolif_genes, senesc_genes,
                                ctr_label = "ctr", tis_label = "tis")

  all_proportions[["GSE131984"]] <- compute_proportions(
    tnbc_seurat, dataset = "GSE131984", cell_line = "SUM159",
    ctr_label = "ctr", tis_label = "tis"
  )
  rm(tnbc_seurat); gc()
}

# ==============================================================================
# 4.  GSE134836 — PC9 NSCLC  (Day 0 vs Day 3 Erlotinib)
# ==============================================================================
{
  message("\n\n══════════════════════════════════════════════════════════")
  message("  GSE134836 — PC9 NSCLC")
  message("══════════════════════════════════════════════════════════")

  data_dir    <- file.path(HYBRID_DATASETS_DIR, "GSE134836_RAW")
  tar_files   <- list.files(data_dir, pattern = "\\.tar\\.gz$", full.names = TRUE)
  seurat_list <- list()

  for (tar_path in tar_files) {
    gsm_id    <- str_extract(basename(tar_path), "GSM[0-9]+")
    day_token <- ifelse(grepl("D0", tar_path), "PC9_D0", "PC9_D3_Erl")
    message("  Unpacking: ", day_token, " (", gsm_id, ")")

    tmp_dir <- file.path(data_dir, paste0("tmp_", gsm_id))
    dir.create(tmp_dir, showWarnings = FALSE)
    untar(tar_path, exdir = tmp_dir)

    all_dirs <- list.dirs(tmp_dir, recursive = TRUE)
    mtx_dir  <- all_dirs[file.exists(file.path(all_dirs, "matrix.mtx.gz"))][1]
    if (is.na(mtx_dir)) stop("matrix.mtx.gz not found in: ", tar_path)

    counts <- Read10X(data.dir = mtx_dir)
    if (is.list(counts)) counts <- counts[["Gene Expression"]]
    colnames(counts) <- paste0(colnames(counts), "_", day_token)
    unlink(tmp_dir, recursive = TRUE)

    sobj              <- CreateSeuratObject(counts = counts, project = day_token)
    sobj$timepoint    <- day_token
    sobj$condition    <- ifelse(day_token == "PC9_D0", "ctr", "tis")
    sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
    seurat_list[[day_token]] <- sobj
    rm(counts, sobj); gc()
  }

  seurat_list <- run_qc_and_doublet_removal(seurat_list)

  pc9_seurat <- merge(seurat_list[[1]], y = seurat_list[-1],
                      add.cell.ids = names(seurat_list))
  pc9_seurat <- JoinLayers(pc9_seurat)
  rm(seurat_list); gc()

  pc9_seurat <- NormalizeData(pc9_seurat, verbose = FALSE)
  pc9_seurat <- detect_hybrids(pc9_seurat, prolif_genes, senesc_genes,
                               condition_col = "condition",
                               ctr_label = "ctr", tis_label = "tis")

  all_proportions[["GSE134836"]] <- compute_proportions(
    pc9_seurat, dataset = "GSE134836", cell_line = "PC9",
    ctr_label = "ctr", tis_label = "tis"
  )
}

# ==============================================================================
# 5.  GSE146544 — AML1566  (Control vs AraC)
# ==============================================================================
{
  message("\n\n══════════════════════════════════════════════════════════")
  message("  GSE146544 — AML1566")
  message("══════════════════════════════════════════════════════════")

  data_dir <- file.path(HYBRID_DATASETS_DIR, "GSE146544")

  raw_ctrl <- ReadMtx(mtx      = file.path(data_dir, "GSM4396313_Ctrl_matrix.mtx.gz"),
                      cells    = file.path(data_dir, "GSM4396313_Ctrl_barcodes.tsv.gz"),
                      features = file.path(data_dir, "GSM4396313_Ctrl_features.tsv.gz"))
  raw_arac <- ReadMtx(mtx      = file.path(data_dir, "GSM4396314_AraC_matrix.mtx.gz"),
                      cells    = file.path(data_dir, "GSM4396314_AraC_barcodes.tsv.gz"),
                      features = file.path(data_dir, "GSM4396314_AraC_features.tsv.gz"))

  seu_ctrl <- CreateSeuratObject(counts = raw_ctrl, project = "AML1566_Ctrl")
  seu_ctrl$condition   <- "Control"
  seu_ctrl$patient     <- "AML1566"
  seu_ctrl[["percent.mt"]] <- PercentageFeatureSet(seu_ctrl, pattern = "^MT-")

  seu_arac <- CreateSeuratObject(counts = raw_arac, project = "AML1566_AraC")
  seu_arac$condition   <- "AraC"
  seu_arac$patient     <- "AML1566"
  seu_arac[["percent.mt"]] <- PercentageFeatureSet(seu_arac, pattern = "^MT-")
  rm(raw_ctrl, raw_arac); gc()

  seurat_list <- list(Control = seu_ctrl, AraC = seu_arac)
  rm(seu_ctrl, seu_arac)
  seurat_list <- run_qc_and_doublet_removal(seurat_list)

  aml_seurat <- merge(seurat_list$Control, y = seurat_list$AraC,
                      add.cell.ids = c("AML1566_Ctrl", "AML1566_AraC"))
  aml_seurat <- JoinLayers(aml_seurat)
  rm(seurat_list); gc()

  aml_seurat <- NormalizeData(aml_seurat, verbose = FALSE)
  aml_seurat <- detect_hybrids(aml_seurat, prolif_genes, senesc_genes,
                               ctr_label = "Control", tis_label = "AraC")

  all_proportions[["GSE146544"]] <- compute_proportions(
    aml_seurat, dataset = "GSE146544", cell_line = "AML1566",
    ctr_label = "Control", tis_label = "AraC"
  )
  rm(aml_seurat); gc()
}

# ==============================================================================
# 6.  GSE168668 — LNCaP prostate  (DMSO vs ENZ48)
# ==============================================================================
{
  message("\n\n══════════════════════════════════════════════════════════")
  message("  GSE168668 — LNCaP prostate")
  message("══════════════════════════════════════════════════════════")

  data_dir  <- file.path(HYBRID_DATASETS_DIR, "GSE168668_RAW")
  prefixes  <- gsub("_matrix\\.mtx\\.gz$", "",
                    list.files(data_dir, pattern = "_matrix\\.mtx\\.gz$"))
  seurat_list <- list()

  for (prefix in prefixes) {
    message("  Loading: ", prefix)
    counts <- ReadMtx(
      mtx            = file.path(data_dir, paste0(prefix, "_matrix.mtx.gz")),
      cells          = file.path(data_dir, paste0(prefix, "_barcodes.tsv.gz")),
      features       = file.path(data_dir, paste0(prefix, "_features.tsv.gz")),
      feature.column = 2
    )
    token <- str_extract(prefix, "LNCaP-[A-Za-z0-9]+")
    sobj  <- CreateSeuratObject(counts = counts, project = token)
    sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
    seurat_list[[prefix]] <- sobj
    rm(counts, sobj); gc()
  }

  seurat_list <- run_qc_and_doublet_removal(seurat_list)

  prostate_seurat <- merge(seurat_list[[1]], y = seurat_list[-1],
                           add.cell.ids = names(seurat_list))
  prostate_seurat <- JoinLayers(prostate_seurat)
  rm(seurat_list); gc()

  prostate_seurat$condition <- case_when(
    prostate_seurat$orig.ident == "LNCaP-DMSO"                          ~ "untreated",
    prostate_seurat$orig.ident == "LNCaP-ENZ48"                         ~ "idling_persistent",
    prostate_seurat$orig.ident %in% c("LNCaP-RESA", "LNCaP-RESB")      ~ "resistant_escaped",
    TRUE                                                                  ~ "Unknown"
  )

  prostate_seurat <- NormalizeData(prostate_seurat, verbose = FALSE)
  prostate_seurat <- detect_hybrids(prostate_seurat, prolif_genes, senesc_genes,
                                    condition_col = "condition",
                                    ctr_label     = "untreated",
                                    tis_label     = "idling_persistent")

  all_proportions[["GSE168668"]] <- compute_proportions(
    prostate_seurat, dataset = "GSE168668", cell_line = "LNCaP",
    condition_col = "condition",
    ctr_label = "untreated", tis_label = "resistant_escaped"
  )
}

# ==============================================================================
# 7.  GSE234181 — MCF7 & T47D breast dormancy  (T0 vs Dormant)
# ==============================================================================
{
  message("\n\n══════════════════════════════════════════════════════════")
  message("  GSE234181 — MCF7 & T47D dormancy")
  message("══════════════════════════════════════════════════════════")

  DATA_DIR <- file.path(HYBRID_DATASETS_DIR, "GSE234181")

  samples <- list(
    MCF7_T0_1   = list(cell_line="MCF7", condition="T0",      treatment="None", timepoint="T0",
                       folder="filtered_feature_bc_matrix_T0_1",   zip="filtered_feature_bc_matrix_T0_1.zip"),
    MCF7_T0_2   = list(cell_line="MCF7", condition="T0",      treatment="None", timepoint="T0",
                       folder="filtered_feature_bc_matrix_T0_2",   zip="filtered_feature_bc_matrix_T0_2.zip"),
    MCF7_AI1m1  = list(cell_line="MCF7", condition="Dormant", treatment="AI",   timepoint="1m",
                       folder="filtered_feature_bc_matrix_AI1m1",  zip="filtered_feature_bc_matrix_AI1m1.zip"),
    MCF7_AI1m2  = list(cell_line="MCF7", condition="Dormant", treatment="AI",   timepoint="1m",
                       folder="filtered_feature_bc_matrix_AI1m2",  zip="filtered_feature_bc_matrix_AI1m2.zip"),
    MCF7_AI2m1  = list(cell_line="MCF7", condition="Dormant", treatment="AI",   timepoint="2m",
                       folder="filtered_feature_bc_matrix_AI2m1",  zip="filtered_feature_bc_matrix_AI2m1.zip"),
    MCF7_AI2m2  = list(cell_line="MCF7", condition="Dormant", treatment="AI",   timepoint="2m",
                       folder="filtered_feature_bc_matrix_AI2m2",  zip="filtered_feature_bc_matrix_AI2m2.zip"),
    MCF7_TAM1m1 = list(cell_line="MCF7", condition="Dormant", treatment="TAM",  timepoint="1m",
                       folder="filtered_feature_bc_matrix_TAM1m1", zip="filtered_feature_bc_matrix_TAM1m1.zip"),
    MCF7_TAM1m2 = list(cell_line="MCF7", condition="Dormant", treatment="TAM",  timepoint="1m",
                       folder="filtered_feature_bc_matrix_TAM1m2", zip="filtered_feature_bc_matrix_TAM1m2.zip"),
    MCF7_TAM2m1 = list(cell_line="MCF7", condition="Dormant", treatment="TAM",  timepoint="2m",
                       folder="filtered_feature_bc_matrix_TAM2m1", zip="filtered_feature_bc_matrix_TAM2m1.zip"),
    MCF7_TAM2m2 = list(cell_line="MCF7", condition="Dormant", treatment="TAM",  timepoint="2m",
                       folder="filtered_feature_bc_matrix_TAM2m2", zip="filtered_feature_bc_matrix_TAM2m2.zip"),
    T47D_T0_1   = list(cell_line="T47D", condition="T0",      treatment="None", timepoint="T0",
                       folder="filtered_feature_bc_matrix_T47D_T0_1", zip="filtered_feature_bc_matrix_T47D_T0_1.zip"),
    T47D_T0_2   = list(cell_line="T47D", condition="T0",      treatment="None", timepoint="T0",
                       folder="filtered_feature_bc_matrix_T47D_T0_2", zip="filtered_feature_bc_matrix_T47D_T0_2.zip"),
    T47D_1m1    = list(cell_line="T47D", condition="Dormant", treatment="AI",   timepoint="1m",
                       folder="filtered_feature_bc_matrix_T47D_1m1",  zip="filtered_feature_bc_matrix_T47D_1m1.zip"),
    T47D_2m1    = list(cell_line="T47D", condition="Dormant", treatment="AI",   timepoint="1m",
                       folder="filtered_feature_bc_matrix_T47D_2m1",  zip="filtered_feature_bc_matrix_T47D_2m1.zip")
  )

  find_matrix_dir <- function(root) {
    dirs <- c(root, list.dirs(root, recursive = TRUE, full.names = TRUE))
    for (d in dirs) {
      f <- list.files(d)
      if (any(grepl("barcodes", f)) && any(grepl("features|genes", f)) && any(grepl("matrix", f)))
        return(d)
    }
    stop("No 10x triplet found under: ", root)
  }

  resolve_sample_dir <- function(meta, base) {
    fp <- file.path(base, meta$folder)
    if (dir.exists(fp)) return(fp)
    zp <- file.path(base, meta$zip)
    if (!file.exists(zp)) stop("Missing folder and zip for: ", meta$folder)
    unzip(zp, exdir = base)
    if (dir.exists(fp)) return(fp)
    stop("Unzipped but folder still missing: ", meta$folder)
  }

  load_10x <- function(mdir) {
    counts <- Read10X(data.dir = mdir)
    if (is.list(counts)) {
      key <- names(counts)[grepl("Gene Expression|RNA|gex", names(counts), ignore.case = TRUE)][1]
      counts <- counts[[key]]
    }
    counts
  }

  seurat_list <- list()
  for (sid in names(samples)) {
    meta <- samples[[sid]]
    message("  Loading: ", sid)
    root <- tryCatch(resolve_sample_dir(meta, DATA_DIR), error = function(e) { warning(e$message); NULL })
    if (is.null(root)) next
    mdir <- tryCatch(find_matrix_dir(root),              error = function(e) { warning(e$message); NULL })
    if (is.null(mdir)) next
    cnts <- tryCatch(load_10x(mdir),                     error = function(e) { warning(e$message); NULL })
    if (is.null(cnts)) next

    seu <- CreateSeuratObject(counts = cnts, project = sid,
                              min.cells = 3, min.features = 200)
    seu$sample_id <- sid
    seu$cell_line <- meta$cell_line
    seu$condition <- meta$condition
    seu$treatment <- meta$treatment
    seu$timepoint <- meta$timepoint
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
    seurat_list[[sid]] <- seu
    rm(cnts, seu); gc()
  }

  seurat_list <- run_qc_and_doublet_removal(seurat_list)

  for (cl in c("MCF7", "T47D")) {
    ids  <- names(seurat_list)[sapply(seurat_list, function(s) s$cell_line[1] == cl)]
    objs <- seurat_list[ids]
    if (length(objs) == 0) { warning("No samples for ", cl); next }

    merged <- if (length(objs) == 1) objs[[1]] else
      merge(objs[[1]], y = objs[-1], add.cell.ids = ids, project = cl)
    merged <- JoinLayers(merged)
    merged <- NormalizeData(merged, verbose = FALSE)
    merged <- detect_hybrids(merged, prolif_genes, senesc_genes,
                             condition_col = "condition",
                             ctr_label = "T0", tis_label = "Dormant")

    all_proportions[[paste0("GSE234181_", cl)]] <- compute_proportions(
      merged, dataset = paste0("GSE234181_", cl), cell_line = cl,
      condition_col = "condition", ctr_label = "T0", tis_label = "Dormant"
    )
  }
}

# ==============================================================================
# 8.  GSE280381 — T47D & MCF7 doxorubicin  (ctr vs tis)
# ==============================================================================
{
  message("\n\n══════════════════════════════════════════════════════════")
  message("  GSE280381 — T47D & MCF7 doxorubicin")
  message("══════════════════════════════════════════════════════════")

  data_dir     <- file.path(HYBRID_DATASETS_DIR, "GSE280381")
  matrix_files <- list.files(data_dir, pattern = "matrix\\.mtx\\.gz$", full.names = TRUE)
  if (length(matrix_files) == 0) stop("No matrix files found in: ", data_dir)

  sample_table <- map_dfr(matrix_files, function(mtx_path) {
    prefix <- sub("_matrix\\.mtx\\.gz$", "", basename(mtx_path))
    parts1 <- strsplit(prefix, "_")[[1]]; parts2 <- strsplit(parts1[2], "-")[[1]]
    tibble(
      sample_id = prefix,
      gsm       = parts1[1],
      cell_line = parts2[1],
      condition = parts2[2],
      replicate = parts2[3],
      mtx       = mtx_path,
      barcodes  = file.path(data_dir, paste0(prefix, "_barcodes.tsv.gz")),
      features  = file.path(data_dir, paste0(prefix, "_features.tsv.gz"))
    )
  })

  seu_list <- setNames(vector("list", nrow(sample_table)), sample_table$sample_id)
  for (i in seq_len(nrow(sample_table))) {
    row <- sample_table[i, ]
    message("  Loading: ", row$sample_id)
    counts <- ReadMtx(mtx = row$mtx, features = row$features, cells = row$barcodes)
    seu    <- CreateSeuratObject(counts = counts, project = row$sample_id,
                                 min.cells = 3, min.features = 200)
    seu$sample_id <- row$sample_id
    seu$cell_line <- row$cell_line
    seu$condition <- row$condition
    seu$replicate <- row$replicate
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
    seu_list[[i]] <- seu
    rm(counts, seu); gc()
  }

  seu_list <- run_qc_and_doublet_removal(seu_list)

  for (cl in c("T47D", "MCF7")) {
    ids <- sample_table$sample_id[sample_table$cell_line == cl]
    ids <- ids[ids %in% names(seu_list)]
    if (length(ids) == 0) { warning("No samples for ", cl); next }

    merged <- merge(seu_list[[ids[1]]], y = seu_list[ids[-1]],
                    add.cell.ids = ids, project = cl)
    merged <- JoinLayers(merged)
    merged <- NormalizeData(merged, verbose = FALSE)
    merged <- detect_hybrids(merged, prolif_genes, senesc_genes,
                             condition_col = "condition",
                             ctr_label = "ctr", tis_label = "tis")

    all_proportions[[paste0("GSE280381_", cl)]] <- compute_proportions(
      merged, dataset = paste0("GSE280381_", cl), cell_line = cl,
      condition_col = "condition", ctr_label = "ctr", tis_label = "tis"
    )
  }
}

# ==============================================================================
# 9.  GSE297393 — MCF7 etoposide  (UT vs 25)
# ==============================================================================
{
  message("\n\n══════════════════════════════════════════════════════════")
  message("  GSE297393 — MCF7 etoposide")
  message("══════════════════════════════════════════════════════════")

  MCF_eto <- readRDS(file.path(HYBRID_DATASETS_DIR, "GSE297393",
                               "GSE297393_merged_filtered_etop.rds"))
  MCF_eto <- subset(MCF_eto, Dose %in% c("UT", "25"))
  MCF_eto$condition <- ifelse(MCF_eto$Dose == "UT", "ctr", "tis")

  MCF_eto[["percent.mt"]] <- PercentageFeatureSet(MCF_eto, pattern = "^MT-")

  seurat_list <- SplitObject(MCF_eto, split.by = "condition")
  rm(MCF_eto); gc()
  seurat_list <- run_qc_and_doublet_removal(seurat_list)

  MCF_eto <- merge(seurat_list[[1]], y = seurat_list[-1],
                   add.cell.ids = names(seurat_list))
  rm(seurat_list); gc()

  MCF_eto <- NormalizeData(MCF_eto, verbose = FALSE)
  MCF_eto <- detect_hybrids(MCF_eto, prolif_genes, senesc_genes,
                            ctr_label = "ctr", tis_label = "tis")

  all_proportions[["GSE297393"]] <- compute_proportions(
    MCF_eto, dataset = "GSE297393", cell_line = "MCF7-eto",
    ctr_label = "ctr", tis_label = "tis"
  )
}

# ==============================================================================
# 10.  Combine all proportions & save
# ==============================================================================

message("\n\n══════════════════════════════════════════════════════════")
message("  Assembling master proportions table")
message("══════════════════════════════════════════════════════════")

prop_df <- bind_rows(Filter(Negate(is.null), all_proportions))
stopifnot("No proportions computed — check upstream errors." = nrow(prop_df) > 0)
print(prop_df)

# ==============================================================================
# 11.  Bubble plot  —  Fig 2F
#       x    = dataset × cell_line
#       y    = 1 (Control)  /  2 (Treated)
#       size = pct_hybrid
# ==============================================================================

message("Generating bubble plot...")

prop_df <- prop_df %>%
  mutate(
    x_label = paste0(dataset, "\n(", cell_line, ")"),
    y_pos   = as.numeric(y_pos)
  )

x_order <- prop_df %>%
  filter(y_pos == 2) %>%
  dplyr::arrange(desc(pct_hybrid)) %>%
  pull(x_label)
x_order <- c(x_order, sort(setdiff(unique(prop_df$x_label), x_order)))

bubble_plot <- ggplot(prop_df,
                      aes(x    = factor(x_label, levels = x_order),
                          y    = y_pos,
                          size = pct_hybrid,
                          fill = group)) +

  geom_point(shape = 21, stroke = 0.5, color = "grey30", alpha = 0.85) +

  scale_size(range = c(0.5, 15), name = "% Hybrid cells") +
  scale_fill_manual(values = c("Control" = "grey", "Treated" = "darkgrey"),
                    name   = "Condition") +
  scale_y_continuous(breaks = c(1, 2),
                     labels = c("Control", "Treated"),
                     limits = c(0.3, 2.7)) +

  labs(x = NULL, y = NULL) +

  theme_classic(base_size = 13) +
  theme(
    axis.text.x        = element_text(angle = 35, hjust = 1, size = 9),
    axis.text.y        = element_text(size = 10, face = "bold"),
    axis.ticks.y       = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4),
    plot.title         = element_text(face = "bold", size = 14),
    legend.position    = "right"
  )

plot_width <- max(10, length(x_order) * 1.2)

save_figure(bubble_plot, "Fig2F_hybrid_meta_analysis_bubble",
            width = plot_width, height = 6)

message("Fig2F saved to: ", file.path(fig_dir, "Fig2F_hybrid_meta_analysis_bubble.pdf"))
message("\n✓  Pipeline complete.")
