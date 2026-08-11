# ============================================================
# T47D scRNA-seq time-course analysis
# Step 1 of 8: QC, doublet removal, Seurat workflow, NMF metaprograms
# ============================================================

source("config.R")
source("R/utils/functions.R")

# -----------------------------
# 0. Libraries
# -----------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(ggpubr)
  library(tidyr)
  library(tidyverse)
  library(plyr)
  library(presto)
  library(destiny)
  library(readr)
  library(data.table)
  library(msigdbr)
  library(fgsea)
  library(dplyr)
  library(DoubletFinder)
  library(RcppML)
  library(Matrix)
  library(purrr)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(fastcluster)
  library(parallelDist)
  library(igraph)
})

# -----------------------------
# 1. Paths & output setup
# -----------------------------
tenx_dir   <- TENX_DIR
output_dir <- OUTPUT_DIR

fig_dir <- file.path(FIGURES_DIR, "Figure_1")

# Sub-folders
dir.create(fig_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

cat("Figures   →", fig_dir,   "\n")
cat("Tables    →", table_dir, "\n")

# Analysis parameters
downsample_n <- 500
n_pcs        <- 20
set.seed(123)

# -----------------------------
# 2. Load 10x data
# -----------------------------
counts <- Read10X(data.dir = tenx_dir)
T47D   <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200)
rm(counts)

# -----------------------------
# 3. QC metrics & metadata
# -----------------------------
T47D[["percent.mt"]] <- PercentageFeatureSet(T47D, pattern = "^MT-")

sample_index <- as.numeric(gsub(".*-(.*)", "\\1", Cells(T47D), perl = TRUE))

sample_labels <- mapvalues(
  sample_index,
  from = 1:6,
  to   = c("0_R1", "0_R2", "7_R1", "7_R2", "15_R1", "15_R2")
)
names(sample_labels) <- Cells(T47D)
T47D <- AddMetaData(T47D, metadata = sample_labels, col.name = "sample")

sample_split <- strsplit(T47D$sample, "_")
sample_mat   <- do.call(rbind, sample_split)
colnames(sample_mat) <- c("time", "replicate")
sample_df <- as.data.frame(sample_mat, stringsAsFactors = FALSE)
rownames(sample_df) <- Cells(T47D)
T47D <- AddMetaData(T47D, metadata = sample_df)
T47D$time <- factor(T47D$time, levels = c("0", "7", "15"))
Idents(T47D) <- "sample"

# Figure S1A — QC before filtering
p_qc_before <- VlnPlot(
  T47D,
  layer    = "counts",
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol     = 1,
  pt.size  = 0.3
)
save_figure(p_qc_before, "FigS1A_QC_before_filtering", width = 8, height = 10)

# -----------------------------
# 4. QC filtering + doublet removal
# -----------------------------
sample.list <- SplitObject(T47D, split.by = "sample")

sample.list <- run_qc_and_doublet_removal(sample.list)

T47D <- merge(sample.list[[1]], y = sample.list[2:length(sample.list)])
T47D <- JoinLayers(T47D)

# Figure S1B — QC after filtering
p_qc_after <- VlnPlot(
  T47D,
  group.by = "sample",
  layer    = "counts",
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol     = 1,
  pt.size  = 0.3
)
save_figure(p_qc_after, "FigS1B_QC_after_filtering", width = 8, height = 10)

# -----------------------------
# 5. Standard Seurat workflow
# -----------------------------
T47D <- NormalizeData(T47D)
T47D <- FindVariableFeatures(T47D)
T47D <- ScaleData(T47D)
T47D <- CellCycleScoring(T47D, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
T47D <- RunPCA(T47D)

p_elbow <- ElbowPlot(T47D, ndims = 50)
save_figure(p_elbow, "FigS1C_PCA_elbow", width = 6, height = 4)

T47D <- RunUMAP(T47D, dims = 1:n_pcs)
T47D$time <- factor(T47D$time, levels = c("0", "7", "15"))

# Figure 1E — UMAP by time
p1E <- DimPlot(T47D, group.by = "time", cols = c("#b2e2e2", "#2ca25f", "#006d2c"))
save_figure(p1E, "Fig1E_UMAP_by_time")

# -----------------------------
# 6. CNV cluster labels
# NOTE: CNV cluster assignments are added to T47D after running
# 02_infercnv.R and 03_cnv_subclonal.R. The cnv_clone metadata
# column is written back to the object in 05_lineage_fate.R.
# -----------------------------

# -----------------------------
# 7. NMF metaprograms
# -----------------------------
T47D_MPs <- nmf_cluster_cells(
  seu            = T47D,
  split_by       = "sample",
  group_col      = "time",
  n_metaprograms = 4
)

curated_map <- T47D_MPs$curated_map

# Table S2 — full MP gene assignments (before manual curation)
save_table(curated_map, "TableS2_MP_gene_assignments_raw")

# -----------------------------
# 8. Heatmap
# -----------------------------
cat("\n=== Preparing heatmap ===\n")
Idents(T47D)  <- "sample"
seurat_small  <- subset(T47D, downsample = downsample_n)
seurat_small  <- ScaleData(seurat_small, features = curated_map$Gene)
mat           <- GetAssayData(seurat_small, layer = "scale.data")

genes_to_plot  <- intersect(curated_map$Gene, rownames(mat))
curated_map_plot <- curated_map[curated_map$Gene %in% genes_to_plot, ] %>% arrange(MP)

mat <- mat[curated_map_plot$Gene, ]
mat[mat >  2] <-  2
mat[mat < -2] <- -2

row_ann <- data.frame(Program = factor(curated_map_plot$MP))
rownames(row_ann) <- curated_map_plot$Gene

col_ann <- data.frame(Day = as.factor(seurat_small$time))
rownames(col_ann) <- colnames(mat)

n_mp <- length(unique(curated_map_plot$MP))
ann_colors <- list(
  Day = c("0" = "#b2e2e2", "7" = "#2ca25f", "15" = "#006d2c"),
  Program = setNames(
    brewer.pal(max(3, min(n_mp, 8)), "Dark2")[seq_len(n_mp)],
    sort(unique(as.character(curated_map_plot$MP)))
  )
)

cat("Plotting", nrow(curated_map_plot), "genes across", ncol(mat), "cells\n")

ph <- ComplexHeatmap::pheatmap(
  mat,
  clustering_method = "ward.D2",
  annotation_col    = col_ann,
  annotation_row    = row_ann,
  annotation_colors = ann_colors,
  treeheight_col    = 2,
  cluster_rows      = FALSE,
  cluster_cols      = TRUE,
  show_colnames     = FALSE,
  show_rownames     = FALSE,
  color             = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
  cutree_cols       = 4,
  fontsize          = 10,
  silent            = TRUE,
  
  # --- the part that matters ---
  use_raster        = TRUE,
  raster_device     = "agg_png",   # or "agg_png" if you have {ragg}
  raster_quality    = 8,            # oversampling factor
  raster_by_magick  = FALSE
)

save_pheatmap(ph, "Fig1G_MP_heatmap", width = 12, height = 8)

# -----------------------------
# 9. Cell-state clustering
# -----------------------------
expr_mat <- as.matrix(GetAssayData(T47D, layer = "data")[curated_map_plot$Gene, ])
expr_mat <- scale(t(expr_mat))

dist_matrix <- parallelDist::parDist(expr_mat, method = "euclidean")
hc_ward     <- fastcluster::hclust(dist_matrix, method = "ward.D2")
clusters    <- cutree(hc_ward, k = 4)

T47D$state_cluster <- paste0("c", clusters[colnames(T47D)])

cat("\nState cluster × time:\n")
print(table(T47D$state_cluster, T47D$time))

# Relabel to match biological ordering — adjust case_when to match your data
T47D$state_cluster <- case_when(
  T47D$state_cluster == "c1" ~ "c1",
  T47D$state_cluster == "c2" ~ "c4",
  T47D$state_cluster == "c3" ~ "c2",
  TRUE                       ~ "c3"
)

# Figure 1G — UMAP by state cluster
p_state <- DimPlot(T47D, group.by = "state_cluster")
save_figure(p_state, "FigS1D_UMAP_state_clusters")

# Table S4 — cell state cluster assignments
cluster_df <- data.frame(cell = colnames(T47D), cluster = T47D$state_cluster)
save_table(cluster_df, "TableS4_cell_state_cluster_assignments")

# -----------------------------
# 10. Enrichment function
# -----------------------------
Enrichment_by_cluster <- function(
    seurat_object,
    min_pct     = 1,
    cluster     = "seurat_clusters",
    pathways,
    collection  = "H",
    species     = "human",
    padj_cutoff = 0.05
) {
  genes_de <- wilcoxauc(
    X = seurat_object[["RNA"]]$data,
    y = as.vector(seurat_object@meta.data[[cluster]])
  ) %>%
    filter(pct_in > min_pct | pct_out > min_pct)
  
  msig_species <- if (species == "human") "Homo sapiens" else "Mus musculus"
  m_all        <- msigdbr(species = msig_species)
  all_pathways <- unique(c(pathways, m_all[m_all$gs_collection == collection, ]$gs_name))
  fgsea_sets   <- dplyr::filter(m_all, gs_name %in% all_pathways) %>%
    split(x = .$gene_symbol, f = .$gs_name)
  
  all_results <- lapply(unique(genes_de$group), function(cl) {
    ranks <- genes_de %>%
      filter(group == cl) %>%
      mutate(rank = auc - mean(auc)) %>%
      arrange(desc(rank)) %>%
      { setNames(.$rank, .$feature) }
    
    fgsea(pathways = fgsea_sets, stats = ranks) %>%
      as_tibble() %>%
      mutate(Cluster = cl) %>%
      arrange(desc(NES))
  })
  
  gsea_df <- bind_rows(all_results)
  
  plot_data <- gsea_df %>%
    filter(pathway %in% pathways) %>%
    mutate(
      significant = !is.na(padj) & padj <= padj_cutoff,
      log_pval    = pmin(-log10(padj), 10),
      pathway     = factor(pathway, levels = rev(pathways))
    )
  
  p <- ggplot(
    filter(plot_data, significant),
    aes(x = Cluster, y = pathway)
  ) +
    geom_point(aes(size = log_pval, fill = NES), shape = 21,
               color = "black", stroke = 0.5) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    scale_size_continuous(range = c(3, 9), name = "-log10(adj. P)") +
    theme_bw() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y  = element_text(color = "black", size = 10),
      axis.title   = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(fill = "NES")
  
  list(plot = p, gsea_df = gsea_df)
}

# -----------------------------
# 11. Figure 1H — pathway enrichment
# -----------------------------
selected_pathways <- c(
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_TGF_BETA_SIGNALING",
  "FRIDMAN_SENESCENCE_UP",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE"
)

enrich_out <- Enrichment_by_cluster(
  seurat_object = T47D,
  cluster       = "state_cluster",
  pathways      = selected_pathways
)

save_figure(enrich_out$plot, "Fig1H_pathway_enrichment_by_state_cluster", width = 8, height = 5)

# Table S5 — full GSEA results
save_table(
  enrich_out$gsea_df %>% mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")),
  "TableS5_GSEA_results_by_state_cluster"
)

cat("\n=== Analysis complete ===\n")
cat("All figures saved to:", fig_dir, "\n")
cat("All tables saved to: ", table_dir, "\n")

