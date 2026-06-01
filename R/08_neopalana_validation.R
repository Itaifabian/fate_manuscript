# =============================================================================
# NeoPalAna (GSE93204) — External Validation of Baseline Hybrid Signature
# Step 8 of 8
#
# Strategy: apply the signature (pos_genes / neg_genes) derived from CORALEEN
# to baseline samples in GSE93204. Split patients into hybrid_high / hybrid_low
# by median score. Validate that hybrid_high patients show a different
# proliferation trajectory across treatment (Baseline → C1D15 → Surgery)
# compared to hybrid_low.
#
# Endpoint: PAM50 proliferation score (independent of ROR-S, avoiding the
# known collinearity between ROR and proliferation at baseline).
#
# Prerequisites: pos_genes and neg_genes must be loaded from the CORALEEN
# analysis before running this script.
#
# Outputs:
#   Main figure : results/figures/Figure_5/Figure5E_*
#   Supp figs   : results/figures/Figure_5/supp/FigureS5*
#   Supp tables : results/supplementary_tables/TableS5*
# =============================================================================

# =============================================================================
# 0. OUTPUT DIRECTORY SETUP
# =============================================================================

source("config.R")
source("R/utils/functions.R")

FIGURE_WIDTH  <- 7
FIGURE_HEIGHT <- 6
FIGURE_DPI    <- 300

RESULTS_DIR  <- OUTPUT_DIR
FIG5_DIR     <- file.path(RESULTS_DIR, "figures", "Figure_5")
SUPP_FIG_DIR <- file.path(FIG5_DIR, "supp")
SUPP_DIR     <- file.path(RESULTS_DIR, "supplementary_tables")
fig_dir      <- FIG5_DIR   # for save_figure() default lookup
table_dir    <- SUPP_DIR   # for save_table() default lookup

dir.create(FIG5_DIR,     recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP_DIR,     recursive = TRUE, showWarnings = FALSE)

message("\n=== NeoPalAna Figure 5E — output root: ", RESULTS_DIR, " ===\n")

# =============================================================================
# 1. LIBRARIES
# =============================================================================
library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(genefu)

stopifnot(exists("pos_genes"), exists("neg_genes"))

# =============================================================================
# 2. DOWNLOAD AND PARSE GSE93204
# =============================================================================
gset <- getGEO("GSE93204", GSEMatrix = TRUE, Annot = TRUE, getGPL = TRUE)
idx  <- if (length(gset) > 1) grep("GPL21185", attr(gset, "names")) else 1
gset <- gset[[idx]]

expression_data <- gset@assayData$exprs
pData           <- gset@phenoData@data

# =============================================================================
# 3. SAMPLE METADATA
# =============================================================================
colData <- pData %>%
  dplyr::select(geo_accession, characteristics_ch2.2, characteristics_ch2.3) %>%
  mutate(
    time       = gsub("^treatment: ", "", characteristics_ch2.2),
    time       = factor(time, levels = c("Baseline", "C1D1", "C1D15", "Surgery")),
    patient_id = gsub("^patient id: ", "", characteristics_ch2.3)
  )

# =============================================================================
# 4. PROBE → GENE COLLAPSE
# =============================================================================
annotation_data <- gset@featureData@data

dat <- as.data.frame(expression_data) %>%
  tibble::rownames_to_column("probe") %>%
  inner_join(annotation_data, by = c("probe" = "ID")) %>%
  filter(!is.na(`Gene symbol`), `Gene symbol` != "")

if ("GI" %in% colnames(dat)) dat$GI <- NULL

expr_gene <- dat %>%
  dplyr::select(-probe) %>%
  group_by(`Gene symbol`) %>%
  summarise(across(where(is.numeric), median), .groups = "drop")

expr_mat           <- as.matrix(expr_gene[, -1])
rownames(expr_mat) <- expr_gene$`Gene symbol`

# =============================================================================
# 5. ALIGN METADATA TO EXPRESSION COLUMNS
# =============================================================================
shared_samples <- intersect(colnames(expr_mat), colData$geo_accession)
message(sprintf("Shared samples between expression and metadata: %d", length(shared_samples)))

expr_mat <- expr_mat[, shared_samples, drop = FALSE]
colData  <- colData[match(shared_samples, colData$geo_accession), , drop = FALSE]
rownames(colData) <- colData$geo_accession

stopifnot(identical(rownames(colData), colnames(expr_mat)))

# =============================================================================
# 6. ROR-S ON SURGERY SAMPLES
# =============================================================================
data(pam50.robust)

pre_treatment_samples <- colData %>%
  filter(time == "Baseline") %>%
  pull(geo_accession)

endpoint_samples <- colData %>%
  filter(time == "Surgery")

common_genes <- intersect(rownames(pam50.robust$centroids), rownames(expr_mat))
message(sprintf("PAM50 genes found in expression matrix: %d / %d",
                length(common_genes), nrow(pam50.robust$centroids)))

baseline_medians <- rowMedians(expr_mat[common_genes, pre_treatment_samples], na.rm = TRUE)
names(baseline_medians) <- common_genes
expr_mat_centered <- expr_mat[common_genes, ] - baseline_medians

exp_input <- t(expr_mat_centered[common_genes, endpoint_samples$geo_accession])

annot_df <- data.frame(
  probe         = common_genes,
  EntrezGene.ID = common_genes,
  gene.symbol   = common_genes,
  row.names     = common_genes,
  stringsAsFactors = FALSE
)

ror_res <- rorS(data = exp_input, annot = annot_df, do.mapping = FALSE)

ror_endpoint <- data.frame(
  patient_id = endpoint_samples$patient_id,
  ROR_Score  = as.numeric(ror_res$score),
  ROR_Risk   = as.vector(ror_res$risk),
  stringsAsFactors = FALSE
)

colData <- left_join(colData, ror_endpoint, by = "patient_id") %>%
  dplyr::rename(Endpoint_ROR_Score = ROR_Score,
                Endpoint_ROR_Risk  = ROR_Risk)

# =============================================================================
# 7. PROLIFERATION SCORE (all timepoints)
# =============================================================================
prolif_genes  <- c("BIRC5","CCNB1","CDC20","NUF2","CEP55",
                   "NDC80","MKI67","PTTG1","RRM2","TYMS","UBE2C")
prolif_probes <- intersect(rownames(expr_mat), prolif_genes)
message(sprintf("Proliferation genes found: %d / %d",
                length(prolif_probes), length(prolif_genes)))

expr_prolif  <- t(scale(t(expr_mat[prolif_probes, , drop = FALSE])))
prolif_score <- colMeans(expr_prolif, na.rm = TRUE)

colData <- colData %>%
  left_join(tibble(geo_accession = names(prolif_score),
                   prolif_score  = as.numeric(prolif_score)),
            by = "geo_accession")

# =============================================================================
# 8. BASELINE SIGNATURE SCORE
# =============================================================================
pos_present <- intersect(pos_genes, rownames(expr_mat))
neg_present <- intersect(neg_genes, rownames(expr_mat))

message(sprintf("Signature genes found — positive: %d / %d  |  negative: %d / %d",
                length(pos_present), length(pos_genes),
                length(neg_present), length(neg_genes)))

stopifnot(length(pos_present) + length(neg_present) > 0)

all_sig_genes <- c(pos_present, neg_present)

baseline_expr  <- expr_mat[all_sig_genes, pre_treatment_samples, drop = FALSE]
gene_means     <- rowMeans(baseline_expr, na.rm = TRUE)
gene_sds       <- apply(baseline_expr, 1, sd, na.rm = TRUE)
gene_sds[gene_sds == 0] <- 1

sig_expr_all <- expr_mat[all_sig_genes, , drop = FALSE]
sig_z        <- (sig_expr_all - gene_means) / gene_sds

# =============================================================================
# 9. BASELINE GROUP ASSIGNMENT (hybrid_high / hybrid_low)
# =============================================================================
colData_naive <- colData %>%
  filter(time %in% c("C1D1", "Baseline")) %>%
  mutate(time = factor(time, levels = c("C1D1", "Baseline"))) %>%
  arrange(patient_id, time) %>%
  distinct(patient_id, .keep_all = TRUE)

colData_naive$escape_sig <- colMeans(
  sig_z[pos_present, colData_naive$geo_accession, drop = FALSE], na.rm = TRUE)

colData_naive$sen_sig <- colMeans(
  sig_z[neg_present, colData_naive$geo_accession, drop = FALSE], na.rm = TRUE)

colData_naive$combined_score <- colData_naive$escape_sig - colData_naive$sen_sig

score_median <- median(colData_naive$combined_score, na.rm = TRUE)
message(sprintf("Signature score median (used for group split): %.4f", score_median))

colData_naive$group <- ifelse(
  colData_naive$combined_score > score_median, "hybrid_high", "hybrid_low"
)

message("Group distribution:")
print(table(colData_naive$group))

# =============================================================================
# 10. COMBINE AND FILTER FOR TRAJECTORY PLOT
# =============================================================================
colData_combined <- left_join(
  colData,
  colData_naive %>% dplyr::select(patient_id, group, combined_score),
  by = "patient_id"
) %>%
  filter(!is.na(group))

message(sprintf("Patients with group assignment: %d", n_distinct(colData_combined$patient_id)))
message(sprintf("Samples in trajectory analysis: %d", nrow(colData_combined)))

# =============================================================================
# 11. SUPPLEMENTARY TABLES
# =============================================================================
message("\n[Saving] Supplementary tables ...")

# TableS5F: per-sample metadata, group assignment, scores
save_table(
  colData_combined %>%
    dplyr::select(geo_accession, patient_id, time,
                  group, combined_score,
                  prolif_score, Endpoint_ROR_Score, Endpoint_ROR_Risk),
  "TableS5F_NeoPalAna_per_sample_scores",
  sheet_name = "NeoPalAna Per-sample Scores"
)

# TableS5G: baseline group assignment summary
save_table(
  colData_naive %>%
    dplyr::select(geo_accession, patient_id, time,
                  escape_sig, sen_sig, combined_score, group),
  "TableS5G_NeoPalAna_group_assignment",
  sheet_name = "NeoPalAna Group Assignment"
)

# TableS5H: signature gene coverage in GSE93204
coverage_df <- data.frame(
  direction       = c("positive", "negative"),
  genes_in_sig    = c(length(pos_genes), length(neg_genes)),
  genes_found     = c(length(pos_present), length(neg_present)),
  pct_coverage    = round(
    c(length(pos_present) / length(pos_genes),
      length(neg_present) / length(neg_genes)) * 100, 1)
)
save_table(coverage_df, "TableS5H_NeoPalAna_signature_coverage",
           sheet_name = "Signature Coverage")

# =============================================================================
# 12. PLOTS
# =============================================================================
message("\n[Saving] Figures ...")

# -- MAIN Figure 5E: Proliferation trajectory by hybrid group -----------------
p1 <- ggplot(
  colData_combined %>% filter(time %in% c("Baseline", "C1D15", "Surgery")),
  aes(x = time, y = prolif_score, fill = group)
) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
  stat_compare_means(label = "p.signif") +
  scale_fill_manual(values = c("hybrid_high" = "darkgrey",
                               "hybrid_low"  = "grey")) +
  labs(x      = "Timepoint",
       y      = "Proliferation score",
       fill   = "Hybrid group") +
  theme_classic(base_size = 13)

print(p1)
save_figure(p1, "Figure5E_NeoPalAna_proliferation_trajectory",
            dir = FIG5_DIR)   # <-- main figure directory

# =============================================================================
# 13. SUPPLEMENTARY FIGURES → Figure_5/supp/
# =============================================================================

# -- Supplementary FigureS5F: Baseline score distribution --------------------
p2 <- ggplot(colData_naive, aes(x = combined_score, fill = group)) +
  geom_histogram(bins = 20, colour = "white", alpha = 0.8) +
  geom_vline(xintercept = score_median, linetype = "dashed", colour = "grey30") +
  scale_fill_manual(values = c("hybrid_high" = "darkgrey",
                               "hybrid_low"  = "grey")) +
  labs(title = "Baseline hybrid signature score distribution",
       x     = "Combined score (pos - neg)",
       y     = "Count",
       fill  = "Group") +
  theme_classic(base_size = 13)

print(p2)
save_figure(p2, "FigureS5F_NeoPalAna_score_distribution",
            dir = SUPP_FIG_DIR)

# -- Supplementary FigureS5G: Proliferation by ROR Risk across timepoints ----
p3 <- colData %>%
  filter(!is.na(Endpoint_ROR_Risk), !is.na(prolif_score)) %>%
  mutate(Endpoint_ROR_Risk = factor(Endpoint_ROR_Risk,
                                    levels = c("Low", "Intermediate", "High"))) %>%
  ggplot(aes(x = time, y = prolif_score, fill = Endpoint_ROR_Risk)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  scale_fill_manual(values = c("Low"          = "lightgrey",
                               "Intermediate" = "grey",
                               "High"         = "darkgrey")) +
  stat_compare_means(comparisons = list(c("Low","Intermediate"),
                                        c("Intermediate","High"),
                                        c("Low","High"))) +
  labs(title    = "Proliferation score across timepoints by Surgery ROR Risk",
       subtitle = "Baseline collinearity motivates trajectory-based validation",
       x        = "Timepoint",
       y        = "Proliferation score",
       fill     = "Surgery ROR Risk") +
  theme_classic(base_size = 13)

print(p3)
save_figure(p3, "FigureS5G_NeoPalAna_proliferation_by_ROR_risk",
            dir = SUPP_FIG_DIR)

# =============================================================================
# 14. SESSION INFO
# =============================================================================
session_file <- file.path(FIG5_DIR, "session_info_NeoPalAna.txt")
sink(session_file)
cat("NeoPalAna Figure 5E pipeline — session info\n")
cat("Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")
print(sessionInfo())
sink()
message("  [saved] session_info_NeoPalAna.txt")

# =============================================================================
# 15. FINAL SUMMARY
# =============================================================================
message("\n", strrep("=", 60))
message("NeoPalAna Figure 5E outputs written to:")
message("  Main figure : ", FIG5_DIR)
message("  Supp figs   : ", SUPP_FIG_DIR)
message("  Supp tables : ", SUPP_DIR)
message(strrep("=", 60))
message(sprintf("  Patients assigned to groups : %d",
                n_distinct(colData_combined$patient_id)))
message(sprintf("  hybrid_high / hybrid_low    : %d / %d",
                sum(colData_naive$group == "hybrid_high"),
                sum(colData_naive$group == "hybrid_low")))
message(sprintf("  Score median (split point)  : %.4f", score_median))
message(sprintf("  Pos genes found : %d / %d  |  Neg genes found : %d / %d",
                length(pos_present), length(pos_genes),
                length(neg_present), length(neg_genes)))
message(strrep("=", 60), "\n")

