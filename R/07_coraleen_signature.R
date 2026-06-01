# =============================================================================
# CORALEEN — Baseline Transcriptomic Signature Predicting Surgery ROR-S
# Step 7 of 8
#
# Pipeline overview:
#   1.  Load and process metadata + count matrix
#   2.  Compute ROR-S score on surgery samples (endpoint)
#   3.  Compute proliferation score across all timepoints
#   4.  Baseline-only VST (no leakage from later timepoints)
#   5.  Bootstrap stability selection — defines the candidate gene set
#   6.  Elastic net on stable genes — fits the signature
#   7.  Stacked LOOCV — unbiased performance estimate
#   8.  Directional score on all baseline samples
#   9.  Plots + auto-save (PDF + PNG) → results/figures/Figure_5/
#   10. Supplementary tables (CSV + XLSX) → results/supplementary_tables/
# =============================================================================

# =============================================================================
# 0. CONSTANTS — all tunable parameters in one place
# =============================================================================
source("config.R")
source("R/utils/functions.R")

MIN_COUNT_FILTER <- 10
R_THRESHOLD      <- 0.33
STABILITY_THRESH <- 0.50
B_BOOT           <- 500
ALPHAS           <- seq(0, 1, by = 0.1)
BASE_DIR      <- dirname(OUTPUT_DIR)
FIGURE_WIDTH  <- 7
FIGURE_HEIGHT <- 6
FIGURE_DPI    <- 300

# =============================================================================
# 0b. OUTPUT DIRECTORY SETUP
#     Mirrors the existing results structure:
#       results/figures/Figure_5/          ← main figure (p2)
#       results/figures/Figure_5/supp/     ← supplementary figures
#       results/supplementary_tables/      ← all CSV/XLSX tables (flat)
# =============================================================================
RESULTS_DIR  <- OUTPUT_DIR
FIG5_DIR     <- file.path(RESULTS_DIR, "figures", "Figure_5")
SUPP_FIG_DIR <- file.path(RESULTS_DIR, "figures", "Figure_5", "supp")
SUPP_DIR     <- file.path(RESULTS_DIR, "supplementary_tables")
fig_dir      <- FIG5_DIR      # for save_figure() default lookup
table_dir    <- SUPP_DIR      # for save_table() default lookup

dir.create(FIG5_DIR,     recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP_DIR,     recursive = TRUE, showWarnings = FALSE)

message("\n=== Output directory: ", BASE_DIR, " ===\n")

# =============================================================================
# 1. LIBRARIES
# =============================================================================
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(jsonlite)
library(DESeq2)
library(biomaRt)
library(genefu)
library(glmnet)
library(ggplot2)
library(ggpubr)
library(scales)

# =============================================================================
# 2. METADATA
# =============================================================================
meta <- read_csv(file.path(CORALEEN_DIR, 'merged_metadata.csv'))

parse_extra_attributes <- function(x) {
  if (is.na(x) || !nzchar(x)) return(tibble(tag = character(), value = character()))
  out <- tryCatch(fromJSON(x), error = function(e) NULL)
  if (is.null(out) || nrow(out) == 0) return(tibble(tag = character(), value = character()))
  tibble(tag = out$tag, value = as.character(out$value))
}

# Build wide data manually
extra_wide_list <- lapply(seq_len(nrow(meta)), function(i) {
  tbl <- parse_extra_attributes(meta$extra_attributes[i])
  if (nrow(tbl) == 0) return(data.frame(row_id = i))
  tbl <- tbl[!duplicated(tbl$tag), ]
  wide <- setNames(as.data.frame(t(tbl$value)), paste0("extra_wide", tbl$tag))
  wide$row_id <- i
  wide
})

extra_wide_df <- dplyr::bind_rows(extra_wide_list)

meta2 <- meta %>%
  dplyr::mutate(row_id = row_number()) %>%
  left_join(extra_wide_df, by = "row_id") %>%
  dplyr::select(-row_id)

setdiff(colnames(meta2), colnames(meta))

meta_fil <- dplyr::filter(meta2,
                          extra_widedescription %in%
                     c('Letro+Ribo_SCR', 'Letro+Ribo_W03', 'Letro+Ribo_SUR'))

# =============================================================================
# 3. COUNT MATRIX
# =============================================================================
counts    <- read_csv(file.path(CORALEEN_DIR, 'expression_matrix_counts.csv'))
ens_clean <- sub("\\..*$", "", counts[[1]])

mart <- useMart("ensembl", 
                dataset = "hsapiens_gene_ensembl",
                host = "https://useast.ensembl.org")  # or uswest.ensembl.org

gene_map <- getBM(
  attributes = c('ensembl_transcript_id', 'external_gene_name'),
  filters    = 'ensembl_transcript_id',
  values     = ens_clean,
  mart       = mart
)

counts$gene_name <- gene_map$external_gene_name[
  match(ens_clean, gene_map$ensembl_transcript_id)]
counts$entrez_id <- gene_map$entrezgene_id[
  match(ens_clean, gene_map$ensembl_transcript_id)]

gene_level_counts <- counts %>%
  filter(!is.na(gene_name) & gene_name != "") %>%
  group_by(gene_name) %>%
  summarise(across(where(is.numeric), sum)) %>%
  ungroup()

final_matrix           <- as.data.frame(gene_level_counts)
rownames(final_matrix) <- final_matrix$gene_name
final_matrix$gene_name <- NULL
colnames(final_matrix) <- sub("^prod_ega-box-\\d+_", "", colnames(final_matrix))

common_samples <- intersect(colnames(final_matrix), meta_fil$standardized)
stopifnot(length(common_samples) > 0)
counts_fil <- final_matrix[, common_samples, drop = FALSE]
meta_fil   <- meta_fil %>% filter(standardized %in% common_samples)

# =============================================================================
# 4. ALL-SAMPLE VST  (ROR-S + proliferation trajectory only)
# =============================================================================
dds_all <- DESeqDataSetFromMatrix(
  countData = round(counts_fil),
  colData   = meta_fil,
  design    = ~ 1
)
dds_all      <- dds_all[rowSums(counts(dds_all)) >= MIN_COUNT_FILTER, ]
vsd_all      <- varianceStabilizingTransformation(dds_all, blind = TRUE)
expr_vst_all <- assay(vsd_all)

# =============================================================================
# 5. ROR-S SCORE ON SURGERY SAMPLES
# =============================================================================
data(pam50.robust)

surgery_ids  <- meta_fil %>% filter(timepoint == "surgery") %>% pull(standardized)
common_pam50 <- intersect(rownames(expr_vst_all), rownames(pam50.robust$centroids))
exp_input_surgery <- t(expr_vst_all[common_pam50, surgery_ids, drop = FALSE])

annot_df <- data.frame(
  probe          = colnames(exp_input_surgery),
  EntrezGene.ID  = colnames(exp_input_surgery),
  row.names      = colnames(exp_input_surgery),
  stringsAsFactors = FALSE
)

ror_surgery <- rorS(data = exp_input_surgery, annot = annot_df, do.mapping = FALSE)

ror_df <- tibble(
  standardized      = rownames(exp_input_surgery),
  Surgery_ROR_Score = as.numeric(ror_surgery$score),
  Surgery_ROR_Risk  = as.vector(ror_surgery$risk)
) %>%
  left_join(meta_fil %>% dplyr::select(standardized, subject_id),
            by = "standardized") %>%
  dplyr::select(subject_id, Surgery_ROR_Score, Surgery_ROR_Risk)

meta_fil <- meta_fil %>% left_join(ror_df, by = "subject_id")

# =============================================================================
# 6. PROLIFERATION SCORE (all timepoints)
# =============================================================================
prolif_genes  <- c("BIRC5","CCNB1","CDC20","NUF2","CEP55",
                   "NDC80","MKI67","PTTG1","RRM2","TYMS","UBE2C")
prolif_probes <- intersect(rownames(expr_vst_all), prolif_genes)
expr_prolif   <- t(scale(t(expr_vst_all[prolif_probes, , drop = FALSE])))
prolif_score  <- colMeans(expr_prolif, na.rm = TRUE)

meta_fil <- meta_fil %>%
  left_join(tibble(standardized = names(prolif_score),
                   prolif_score  = prolif_score),
            by = "standardized")

meta_fil$timepoint <- factor(meta_fil$timepoint,
                             levels = c('pre treatment', 'on-treatment', 'surgery'))

# =============================================================================
# 7. BASELINE-ONLY VST
# =============================================================================
pre_ids <- meta_fil %>%
  filter(timepoint == "pre treatment", !is.na(Surgery_ROR_Score)) %>%
  pull(standardized)

meta_pre <- meta_fil %>%
  filter(standardized %in% pre_ids,
         !is.na(Surgery_ROR_Score),
         !is.na(prolif_score))

dds_pre <- DESeqDataSetFromMatrix(
  countData = round(counts_fil[, pre_ids, drop = FALSE]),
  colData   = meta_fil %>% filter(standardized %in% pre_ids),
  design    = ~ 1
)
dds_pre      <- dds_pre[rowSums(counts(dds_pre)) >= MIN_COUNT_FILTER, ]
vsd_pre      <- varianceStabilizingTransformation(dds_pre, blind = TRUE)
expr_vst_pre <- assay(vsd_pre)

# =============================================================================
# 8. FEATURE MATRIX
# =============================================================================
stopifnot(exists("fate_markers_all"))
genes_test <- intersect(rownames(fate_markers_all), rownames(expr_vst_pre))
stopifnot(length(genes_test) > 0)
message(sprintf("Candidate genes in VST matrix: %d", length(genes_test)))

X_raw     <- t(expr_vst_pre[genes_test, meta_pre$standardized, drop = FALSE])
y_ror     <- meta_pre$Surgery_ROR_Score
n_samples <- nrow(X_raw)
message(sprintf("Baseline samples for modelling: %d", n_samples))

# =============================================================================
# 9. BOOTSTRAP STABILITY SELECTION
# =============================================================================
set.seed(321)

boot_genes <- replicate(B_BOOT, simplify = FALSE, {
  idx <- sample(seq_len(n_samples), replace = TRUE)
  Xb  <- X_raw[idx, , drop = FALSE]
  yb  <- y_ror[idx]
  rb  <- apply(Xb, 2, function(x)
    tryCatch(cor(x, yb, method = "spearman", use = "complete.obs"),
             error = function(e) NA_real_))
  rb  <- rb[!is.na(rb)]
  names(rb)[abs(rb) >= R_THRESHOLD]
})

gene_freq    <- sort(table(unlist(boot_genes)) / B_BOOT, decreasing = TRUE)
stable_genes <- names(gene_freq[gene_freq >= STABILITY_THRESH])

stability_df <- data.frame(
  gene           = names(gene_freq),
  selection_freq = as.numeric(gene_freq),
  stable         = names(gene_freq) %in% stable_genes
) %>% arrange(desc(selection_freq))

message(sprintf("Stable genes: %s", paste(stable_genes, collapse = ", ")))
print(head(stability_df, 20))
stopifnot(length(stable_genes) >= 2)

# --- Save Supplementary Table S5A: Bootstrap stability scores ----------------
message("\n[Saving] Supplementary tables ...")
save_table(stability_df, "TableS5A_bootstrap_stability",
           sheet_name = "Bootstrap Stability")
# =============================================================================
# 10. ELASTIC NET ON STABLE GENES
# =============================================================================
X_stable     <- scale(X_raw[, stable_genes, drop = FALSE])
scale_center <- attr(X_stable, "scaled:center")
scale_sd     <- attr(X_stable, "scaled:scale")

message(sprintf("\nFitting elastic net on %d stable genes × %d samples ...",
                length(stable_genes), n_samples))

set.seed(42)
foldid <- sample(rep(seq(5), length.out = n_samples))

cv_results <- lapply(ALPHAS, function(a) {
  cv.glmnet(X_stable, y_ror,
            alpha        = a,
            foldid       = foldid,      # same folds for every alpha
            standardize  = FALSE,
            type.measure = "mse")
})

min_mse    <- sapply(cv_results, function(m) min(m$cvm))
best_idx   <- which.min(min_mse)
best_alpha <- ALPHAS[best_idx]
best_model <- cv_results[[best_idx]]

message(sprintf("Best alpha: %.2f  |  lambda.min: %.4f  |  5-fold CV MSE: %.3f",
                best_alpha, best_model$lambda.min, min_mse[best_idx]))

final_coeffs <- as.matrix(coef(best_model, s = "lambda.min"))
signature_df <- data.frame(
  gene        = rownames(final_coeffs),
  coefficient = as.numeric(final_coeffs)
) %>%
  filter(gene != "(Intercept)", coefficient != 0) %>%
  arrange(desc(abs(coefficient)))

message(sprintf("\nFinal signature genes (n=%d):", nrow(signature_df)))
print(signature_df)

pos_genes <- signature_df$gene[signature_df$coefficient > 0]
neg_genes  <- signature_df$gene[signature_df$coefficient < 0]

# --- Save Supplementary Table S5B: Elastic net signature coefficients --------
save_table(signature_df, "TableS5B_signature_coefficients",
           sheet_name = "Signature Coefficients")

# =============================================================================
# 11. STACKED LOOCV
# =============================================================================
loocv_preds <- numeric(n_samples)

for (i in seq_len(n_samples)) {
  X_tr  <- X_raw[-i, stable_genes, drop = FALSE]
  X_te  <- X_raw[ i, stable_genes, drop = FALSE]
  y_tr  <- y_ror[-i]
  ctr   <- colMeans(X_tr)
  sds   <- apply(X_tr, 2, sd); sds[sds == 0] <- 1
  Xtr_s <- scale(X_tr, center = ctr, scale = sds)
  Xte_s <- scale(X_te, center = ctr, scale = sds)
  
  m_i <- glmnet(Xtr_s, y_tr,
                alpha       = best_alpha,
                lambda      = best_model$lambda.min,
                standardize = FALSE)
  loocv_preds[i] <- as.numeric(predict(m_i, newx = Xte_s,
                                       s = best_model$lambda.min))
}

loocv_r  <- cor(loocv_preds, y_ror, method = "spearman")
loocv_r2 <- loocv_r^2

message(sprintf("\nLOOCV Spearman r  : %.3f", loocv_r))
message(sprintf("LOOCV Spearman R2 : %.3f", loocv_r2))

# --- Save Supplementary Table S5C: LOOCV predictions vs true ROR-S ----------
loocv_df <- tibble(
  standardized    = meta_pre$standardized,
  subject_id      = meta_pre$subject_id,
  Surgery_ROR_Score     = y_ror,
  LOOCV_predicted_score = loocv_preds
)
save_table(loocv_df, "TableS5C_LOOCV_predictions",
           sheet_name = "LOOCV Predictions")

# =============================================================================
# 12. DIRECTIONAL SCORE ON ALL BASELINE SAMPLES
# =============================================================================
X_z <- scale(X_raw[, stable_genes, drop = FALSE],
             center = scale_center, scale = scale_sd)

compute_dir_score <- function(mat, pos, neg) {
  pos_ok <- intersect(pos, colnames(mat))
  neg_ok <- intersect(neg, colnames(mat))
  pos_s  <- if (length(pos_ok) > 0) rowMeans(mat[, pos_ok, drop = FALSE]) else 0
  neg_s  <- if (length(neg_ok) > 0) rowMeans(mat[, neg_ok, drop = FALSE]) else 0
  pos_s - neg_s
}

mean_sig_score <- compute_dir_score(X_z, pos_genes, neg_genes)
names(mean_sig_score) <- rownames(X_z)

meta_pre <- meta_pre %>%
  left_join(tibble(standardized   = names(mean_sig_score),
                   mean_sig_score = as.numeric(mean_sig_score)),
            by = "standardized") %>%
  mutate(Surgery_ROR_Risk = factor(Surgery_ROR_Risk,
                                   levels = c("Low", "Intermediate", "High")))

# =============================================================================
# 13. VALIDATION METRICS
# =============================================================================
cor_ror    <- cor(meta_pre$mean_sig_score, meta_pre$Surgery_ROR_Score,
                  method = "spearman", use = "complete.obs")
cor_prolif <- cor(meta_pre$mean_sig_score, meta_pre$prolif_score,
                  method = "spearman", use = "complete.obs")
kw         <- kruskal.test(mean_sig_score ~ Surgery_ROR_Risk, data = meta_pre)

message(sprintf("\nIn-sample Spearman r (score vs ROR)   : %.3f", cor_ror))
message(sprintf("In-sample Spearman r (score vs prolif): %.3f  [exploratory]", cor_prolif))
message(sprintf("Kruskal-Wallis across Low/Int/High    : p = %.4f", kw$p.value))
message(sprintf("LOOCV Spearman r                      : %.3f", loocv_r))
message(sprintf("LOOCV Spearman R2                     : %.3f", loocv_r2))

# --- Save Supplementary Table S5D: per-sample scores + metadata --------------
scores_out <- meta_pre %>%
  dplyr::select(standardized, subject_id, timepoint,
                Surgery_ROR_Score, Surgery_ROR_Risk,
                mean_sig_score, prolif_score)
save_table(scores_out, "TableS5D_per_sample_scores",
           sheet_name = "Per-sample Scores")

# --- Save Supplementary Table S5E: summary statistics -----------------------
summary_stats <- data.frame(
  metric  = c("In-sample Spearman r (score vs ROR-S)",
              "In-sample Spearman r (score vs prolif)",
              "Kruskal-Wallis p (Low/Int/High)",
              "LOOCV Spearman r",
              "LOOCV Spearman R2",
              "n (baseline samples)",
              "n (stable genes)",
              "n (signature genes)",
              "Bootstrap resamples",
              "Stability threshold",
              "Spearman r threshold per bootstrap",
              "Best elastic net alpha",
              "lambda.min"),
  value   = c(round(cor_ror, 4),
              round(cor_prolif, 4),
              signif(kw$p.value, 4),
              round(loocv_r, 4),
              round(loocv_r2, 4),
              n_samples,
              length(stable_genes),
              nrow(signature_df),
              B_BOOT,
              STABILITY_THRESH,
              R_THRESHOLD,
              round(best_alpha, 2),
              round(best_model$lambda.min, 6))
)
save_table(summary_stats, "TableS5E_summary_statistics",
           sheet_name = "Summary Statistics")

# =============================================================================
# 14. PLOTS
#     Main figure : p2 only → Figure_5/Figure5C_*
#     All others  → Figure_5/supplementary/figures/FigureS5*
# =============================================================================
message("\n[Saving] Figures ...")

# -- MAIN Figure 5C: Signature score across ROR risk groups -------------------
p2 <- ggplot(meta_pre, aes(x = Surgery_ROR_Risk, y = mean_sig_score)) +
  geom_boxplot(aes(fill = Surgery_ROR_Risk), alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
  scale_fill_manual(values = c("Low"          = "lightgrey",
                               "Intermediate" = "grey",
                               "High"         = "darkgrey")) +
  stat_compare_means(
    comparisons = list(c("Low","Intermediate"),
                       c("Intermediate","High"),
                       c("Low","High")),
    label = "p.signif", method = "wilcox.test") +
  labs(x = "Surgery ROR Risk",
       y = "Signature score (Pre-treatment)") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

print(p2)
save_figure(p2, "Figure5C_signature_by_ROR_risk_group",
            dir = BASE_DIR)   # <-- main figure directory

# =============================================================================
# 15. SUPPLEMENTARY FIGURES → Figure_5/supplementary/figures/
# =============================================================================

# -- Supplementary Figure S5A: Signature score vs continuous ROR-S -----------
p1 <- ggplot(meta_pre, aes(x = mean_sig_score, y = Surgery_ROR_Score,
                           colour = Surgery_ROR_Risk)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey40", linetype = "dashed") +
  scale_colour_manual(values = c("Low" = "#2166ac",
                                 "Intermediate" = "#f4a460",
                                 "High" = "#d6604d")) +
  labs(title    = "Baseline signature score vs Surgery ROR Score",
       subtitle = sprintf("LOOCV R\u00b2 = %.2f  |  n = %d", loocv_r2, nrow(meta_pre)),
       x        = "Signature score (Pre-treatment)",
       y        = "Surgery ROR Score",
       colour   = "ROR Risk") +
  theme_classic(base_size = 13)

print(p1)
save_figure(p1, "FigureS5A_signature_vs_ROR_scatter",
            dir = SUPP_FIG_DIR)

# -- Supplementary Figure S5B: Proliferation score across timepoints ----------
plot_data_prolif <- meta_fil %>%
  dplyr::select(timepoint, prolif_score, Surgery_ROR_Risk) %>%
  filter(!is.na(Surgery_ROR_Risk)) %>%
  mutate(Surgery_ROR_Risk = factor(Surgery_ROR_Risk,
                                   levels = c("Low", "Intermediate", "High"))) %>%
  na.omit()

p3 <- ggplot(plot_data_prolif,
             aes(x = timepoint, y = prolif_score, fill = Surgery_ROR_Risk)) +
  geom_boxplot(alpha = 0.7) +
  stat_compare_means(label = "p.signif") +
  scale_fill_manual(values = c("Low"          = "#2166ac",
                               "Intermediate" = "#f4a460",
                               "High"         = "#d6604d")) +
  labs(title = "Proliferation score across timepoints",
       x     = "Timepoint",
       y     = "Proliferation Score",
       fill  = "Surgery ROR Risk") +
  theme_classic(base_size = 13)

print(p3)
save_figure(p3, "FigureS5B_proliferation_by_timepoint",
            dir = SUPP_FIG_DIR)

# -- Supplementary Figure S5C: Bootstrap stability bar chart ------------------
top_n_genes  <- min(30, nrow(stability_df))
stab_plot_df <- head(stability_df, top_n_genes) %>%
  mutate(gene = factor(gene, levels = rev(gene)))

p4 <- ggplot(stab_plot_df, aes(x = gene, y = selection_freq, fill = stable)) +
  geom_col() +
  geom_hline(yintercept = STABILITY_THRESH, linetype = "dashed",
             colour = "red", linewidth = 0.7) +
  scale_fill_manual(values = c("TRUE" = "#4d9221", "FALSE" = "#bababa"),
                    labels = c("TRUE" = "Stable", "FALSE" = "Not stable")) +
  coord_flip() +
  labs(title    = "Bootstrap stability selection",
       subtitle = sprintf("Top %d genes  |  threshold = %.0f%%  |  B = %d resamples",
                          top_n_genes, STABILITY_THRESH * 100, B_BOOT),
       x        = "Gene",
       y        = "Selection frequency",
       fill     = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")

print(p4)
save_figure(p4, "FigureS5C_bootstrap_stability",
            width = 6, height = max(5, top_n_genes * 0.22 + 2),
            dir = SUPP_FIG_DIR)

# -- Supplementary Figure S5D: Elastic net coefficients ----------------------
coef_plot_df <- signature_df %>%
  mutate(direction = ifelse(coefficient > 0, "Positive", "Negative"),
         gene      = factor(gene, levels = gene[order(coefficient)]))

p5 <- ggplot(coef_plot_df, aes(x = gene, y = coefficient, fill = direction)) +
  geom_col() +
  scale_fill_manual(values = c("Positive" = "#d6604d", "Negative" = "#2166ac")) +
  coord_flip() +
  labs(title    = "Elastic net signature coefficients",
       subtitle = sprintf("alpha = %.2f  |  lambda.min = %.4f  |  n genes = %d",
                          best_alpha, best_model$lambda.min, nrow(signature_df)),
       x        = "Gene",
       y        = "Coefficient",
       fill     = "Direction") +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")

print(p5)
save_figure(p5, "FigureS5D_elastic_net_coefficients",
            height = max(4, nrow(signature_df) * 0.35 + 2),
            dir = SUPP_FIG_DIR)

# -- Supplementary Figure S5E: LOOCV predicted vs true ROR-S -----------------
p6 <- ggplot(loocv_df, aes(x = LOOCV_predicted_score, y = Surgery_ROR_Score)) +
  geom_point(size = 3, alpha = 0.8, colour = "#4d9221") +
  geom_smooth(method = "lm", se = TRUE, colour = "grey40", linetype = "dashed") +
  labs(title    = "LOOCV: Predicted vs Observed ROR-S",
       subtitle = sprintf("Spearman r = %.3f  |  R\u00b2 = %.3f  |  n = %d",
                          loocv_r, loocv_r2, n_samples),
       x        = "LOOCV predicted score",
       y        = "Surgery ROR Score") +
  theme_classic(base_size = 13)

print(p6)
save_figure(p6, "FigureS5E_LOOCV_predicted_vs_observed",
            dir = SUPP_FIG_DIR)

# =============================================================================
# 16. SESSION INFO — reproducibility record
# =============================================================================
session_file <- file.path(FIG5_DIR, "session_info.txt")
sink(session_file)
cat("CORALEEN Figure 5 pipeline — session info\n")
cat("Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")
print(sessionInfo())
sink()
message(sprintf("  [saved] session_info.txt"))

# =============================================================================
# 17. FINAL SUMMARY
# =============================================================================
message("\n", strrep("=", 60))
message("Figure 5 outputs written to:")
message("  Main figure : ", FIG5_DIR)
message("  Supp figs   : ", SUPP_FIG_DIR)
message("  Supp tables : ", SUPP_DIR)
message(strrep("=", 60))
message(sprintf("  Stable genes (n=%d) : %s",
                length(stable_genes), paste(stable_genes, collapse = ", ")))
message(sprintf("  Signature genes     : %d", nrow(signature_df)))
message(sprintf("  LOOCV Spearman r    : %.3f  (R2 = %.3f)", loocv_r, loocv_r2))
message(strrep("=", 60), "\n")

