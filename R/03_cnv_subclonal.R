source("config.R")
source("R/utils/functions.R")

library(ggplot2)
library(patchwork)
library(dplyr)
library(data.table)
library(infercnv)
library(fishplot)
library(png)
library(scales)

# ============================================================================
# PATHS & CONFIG
# ============================================================================
WES_OUT <- WES_DIR
GTF     <- GTF_FILE
fig_dir <- file.path(FIGURES_DIR, "Figure_1")

MIN_CELLS <- 50
CHR_ORDER <- paste0("chr", c(1:22, "X"))

REGIONS <- list(
  list(chrom="chr1",  start=143129293, end=248918363, label="chr1q"),
  list(chrom="chr17", start=18241885,  end=83094575,  label="chr17"),
  list(chrom="chr20", start=1619869,   end=64273600,  label="chr20")
)

# ============================================================================
# 1. LOAD WES DATA
# ============================================================================
cnr <- read.table(paste0(WES_OUT, "/cnv/T47D_Orig.cnr"), sep="\t", header=TRUE)
cnr$midpos <- (as.numeric(cnr$start) + as.numeric(cnr$end)) / 2

seg <- read.table(paste0(WES_OUT, "/cnv/T47D_Orig.seg"), sep="\t", header=TRUE)

vcf <- fread(
  cmd=paste0("zcat ", WES_OUT, "/variants/T47D_Orig_PASS_SNPs.vcf.gz | grep -v '^#'"),
  sep="\t", header=FALSE,
  col.names=c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE"))
vcf$AF <- as.numeric(sapply(strsplit(as.character(vcf$SAMPLE), ":"), function(x) {
  fmt <- strsplit(as.character(vcf$FORMAT[1]), ":")[[1]]
  x[which(fmt == "AF")]
}))

# ============================================================================
# 2. LOAD INFERCNV
# ============================================================================
infercnv_obj <- readRDS(file.path(INFERCNV_DIR, "run.final.infercnv_obj"))
cnv_mat      <- infercnv_obj@expr.data
raw_counts   <- infercnv_obj@count.data

valid_subclusters <- Filter(
  function(x) length(x) >= MIN_CELLS,
  infercnv_obj@tumor_subclusters$subclusters[["Tumor"]]
)
subcluster_cols <- names(valid_subclusters)
s1 <- subcluster_cols[1]
s2 <- subcluster_cols[2]

subcluster_means <- lapply(subcluster_cols, function(s) {
  cells <- intersect(names(valid_subclusters[[s]]), colnames(cnv_mat))
  rowMeans(cnv_mat[, cells, drop=FALSE])
})
names(subcluster_means) <- subcluster_cols
s1_col <- paste0(s1, "_raw")
s2_col <- paste0(s2, "_raw")

# ============================================================================
# 3. GENE COORDINATES
# ============================================================================
gene_coords <- fread(
  cmd=paste0("grep -v '^#' '", GTF, "' | awk '$3==\"gene\"' | ",
             "awk '{match($0, /gene_name \"([^\"]+)\"/, a); ",
             "print $1\"\\t\"$4\"\\t\"$5\"\\t\"a[1]}'"),
  sep="\t", header=FALSE, col.names=c("chr","start","end","gene"))
gene_coords <- gene_coords[gene_coords$chr %in% paste0("chr", c(1:22,"X","Y")), ]
gene_coords <- gene_coords[!duplicated(gene_coords$gene), ]

# ============================================================================
# 4. GENE TABLE
# ============================================================================
gene_table           <- data.frame(gene=rownames(cnv_mat), stringsAsFactors=FALSE)
gene_table$chr       <- gene_coords$chr[match(gene_table$gene, gene_coords$gene)]
gene_table$start     <- gene_coords$start[match(gene_table$gene, gene_coords$gene)]
gene_table[[s1_col]] <- subcluster_means[[s1]][gene_table$gene]
gene_table[[s2_col]] <- subcluster_means[[s2]][gene_table$gene]

# ============================================================================
# 5. COLOR SCHEME  —  CNV1=green, CNV2=purple
# ============================================================================
subclone_colors <- c(setNames("#1a7a1a", s1), setNames("#6a0dad", s2))
subclone_labels <- c(setNames("CNV1", s1),    setNames("CNV2", s2))

# ============================================================================
# 6. EVIDENCE TABLE (for highlight colors)
# ============================================================================
evidence <- do.call(rbind, lapply(REGIONS, function(r) {
  genes_r <- gene_table[!is.na(gene_table$chr) &
                          gene_table$chr   == r$chrom &
                          gene_table$start >= r$start &
                          gene_table$start <= r$end, ]
  s1_m <- mean(genes_r[[s1_col]], na.rm=TRUE)
  s2_m <- mean(genes_r[[s2_col]], na.rm=TRUE)
  data.frame(region=r$label,
             higher_subclone=ifelse(s1_m > s2_m, s1, s2),
             stringsAsFactors=FALSE)
}))

# ============================================================================
# 7. FISHPLOT  (shared between both figures)
# ============================================================================
pyclone_path <- file.path(WES_OUT, "variants/t47d_strict_subclones.tsv")
mut_data     <- read.delim(pyclone_path)

fish_meta <- mut_data %>%
  group_by(cluster_id) %>%
  dplyr::summarise(mean_ccf=mean(cellular_prevalence), mut_count=n(), .groups="drop") %>%
  dplyr::filter(mut_count > 50) %>%
  dplyr::arrange(desc(mean_ccf)) %>%
  dplyr::mutate(
    rank     = row_number(),
    is_trunk = rank == 1,
    ccf_pct  = round(mean_ccf * 100, 1),
    display  = case_when(
      is_trunk ~ "Trunk",
      rank == 2 ~ "CNV1",
      rank == 3 ~ "CNV2",
      TRUE      ~ paste0("Clone ", cluster_id)
    )
  )

n_clones  <- nrow(fish_meta)
trunk_row <- which(fish_meta$is_trunk)
subclone_rows <- which(!fish_meta$is_trunk)
subclone_sum  <- sum(fish_meta$ccf_pct[subclone_rows])
if (subclone_sum > fish_meta$ccf_pct[trunk_row]) {
  scale_f <- (fish_meta$ccf_pct[trunk_row] - 0.1) / subclone_sum
  fish_meta$ccf_pct[subclone_rows] <- round(fish_meta$ccf_pct[subclone_rows] * scale_f, 1)
}

timepoints <- c(0, 50, 100)
frac.table <- matrix(0, nrow=n_clones, ncol=3)
frac.table[trunk_row, 1] <- 1
frac.table[, 2] <- fish_meta$ccf_pct
frac.table[, 3] <- fish_meta$ccf_pct
parents            <- rep(trunk_row, n_clones)
parents[trunk_row] <- 0

extra_colors <- c("#d46b08","#0066cc","#c41d7f","#389e0d")
clone_colors <- character(n_clones)
clone_colors[trunk_row] <- "#aaaaaa"
ei <- 1
for (i in subclone_rows) {
  clone_colors[i] <- switch(fish_meta$display[i],
                            "CNV1" = "#1a7a1a", "CNV2" = "#6a0dad",
                            { col <- extra_colors[min(ei, length(extra_colors))]; ei <- ei+1; col })
}

fish <- createFishObject(frac.table, parents,
                         timepoints=timepoints, clone.labels=fish_meta$display)
fish <- layoutClones(fish)
fish <- setCol(fish, col=clone_colors)

# Render to PNG → ggplot grob (base R device workaround)
make_fish_grob <- function(width=2400, height=900) {
  fish_png <- tempfile(fileext=".png")
  while (dev.cur() > 1) dev.off()
  png(fish_png, width=width, height=height, res=200)
  fishPlot(fish, shape="spline", vlines=50, vlab="T47D",
           bg.type="solid", bg.col="white", cex.title=0)
  drawLegend(fish, nrow=1, cex=0.8)
  dev.off()
  while (dev.cur() > 1) dev.off()
  img <- png::readPNG(fish_png)
  ggplot() +
    annotation_raster(img, xmin=0, xmax=1, ymin=0, ymax=1) +
    coord_cartesian(xlim=c(0,1), ylim=c(0,1), expand=FALSE) +
    theme_void()
}

# ============================================================================
# 8. VAF PANELS  (used in both figures)
# ============================================================================
make_vaf_panel <- function(r) {
  col     <- subclone_colors[evidence$higher_subclone[evidence$region == r$label]]
  vcf_sub <- vcf[vcf$CHROM == r$chrom & vcf$POS >= r$start & vcf$POS <= r$end, ]
  ggplot(vcf_sub, aes(x=POS, y=AF)) +
    geom_point(size=1.5, alpha=0.6, color=col, shape=16) +
    geom_hline(yintercept=0.5,  linetype="dashed", color="black",  linewidth=0.5) +
    geom_hline(yintercept=0.33, linetype="dotted", color="gray60", linewidth=0.4) +
    geom_hline(yintercept=0.67, linetype="dotted", color="gray60", linewidth=0.4) +
    coord_cartesian(ylim=c(0,1)) +
    labs(x=NULL, y="VAF", title=r$label) +
    scale_x_continuous(labels=function(x) paste0(round(x/1e6,0),"Mb")) +
    theme_classic(base_size=9) +
    theme(plot.title=element_text(face="bold", size=9))
}

vaf_panels <- lapply(REGIONS, make_vaf_panel)

# ============================================================================
# 9. SCALED EXPRESSION HEATMAP  (main figure bottom row)
# ============================================================================
cell_totals <- colSums(raw_counts)
norm_counts <- sweep(raw_counts, 2, cell_totals / 1e4, FUN="/")
log_counts  <- log1p(norm_counts)
gene_means  <- rowMeans(log_counts)
gene_sds    <- apply(log_counts, 1, sd)
gene_sds[gene_sds == 0] <- 1
scaled_mat  <- (log_counts - gene_means) / gene_sds

expr_long <- do.call(rbind, lapply(REGIONS, function(r) {
  genes_r <- gene_table[!is.na(gene_table$chr) &
                          gene_table$chr   == r$chrom &
                          gene_table$start >= r$start &
                          gene_table$start <= r$end, ]
  if (nrow(genes_r) == 0) return(NULL)
  genes_r  <- genes_r[order(genes_r$start), ]
  genes_r$gene_rank <- seq_len(nrow(genes_r))
  genes_use <- intersect(genes_r$gene, rownames(scaled_mat))
  genes_r   <- genes_r[genes_r$gene %in% genes_use, ]
  
  do.call(rbind, lapply(c(s1, s2), function(s) {
    cells <- intersect(names(valid_subclusters[[s]]), colnames(scaled_mat))
    mat_s <- scaled_mat[genes_r$gene, cells, drop=FALSE]
    data.frame(region=r$label, subclone=s,
               gene=genes_r$gene, gene_rank=genes_r$gene_rank, pos=genes_r$start,
               mean_expr=rowMeans(mat_s), stringsAsFactors=FALSE)
  }))
})) %>%
  mutate(
    subclone       = factor(subclone, levels=c(s1, s2)),
    region         = factor(region, levels=sapply(REGIONS, `[[`, "label")),
    subclone_label = factor(subclone_labels[as.character(subclone)],
                            levels=rev(unname(subclone_labels))),
    expr_clipped   = pmax(pmin(mean_expr, 1), -1)
  )

make_expr_panel <- function(r) {
  df_r <- expr_long[expr_long$region == r$label, ]
  ranks    <- sort(unique(df_r$gene_rank))
  tick_idx <- ranks[round(seq(1, length(ranks), length.out=5))]
  tick_pos <- df_r$pos[match(tick_idx, df_r$gene_rank)]
  tick_labs <- paste0(round(tick_pos/1e6, 0), "Mb")
  
  ggplot(df_r, aes(x=gene_rank, y=subclone_label, fill=expr_clipped)) +
    geom_tile(height=0.9, width=1) +
    scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b",
                         midpoint=0, limits=c(-1,1), oob=squish,
                         name="Scaled expr") +
    scale_x_continuous(breaks=tick_idx, labels=tick_labs, expand=c(0,0)) +
    labs(x=NULL, y=NULL) +
    theme_classic(base_size=9) +
    theme(axis.text.y=element_text(size=8, face="bold"),
          axis.ticks.y=element_blank(),
          legend.position="bottom",
          legend.key.width=unit(1.2,"cm"),
          legend.key.height=unit(0.3,"cm"),
          legend.text=element_text(size=7),
          legend.title=element_text(size=7))
}

expr_panels <- lapply(REGIONS, make_expr_panel)

# ============================================================================
# 10. LOGRR PANELS  (supplementary only)
# ============================================================================
make_logr_panel <- function(r) {
  col     <- subclone_colors[evidence$higher_subclone[evidence$region == r$label]]
  cnr_sub <- cnr[cnr$chromosome == r$chrom &
                   cnr$start >= r$start & cnr$end <= r$end, ]
  ggplot(cnr_sub, aes(x=midpos, y=log2)) +
    geom_point(size=0.2, alpha=0.3, color="gray50") +
    geom_hline(yintercept=0, linetype="dashed", linewidth=0.4) +
    geom_smooth(method="loess", span=0.15, color=col, linewidth=0.8, se=FALSE) +
    annotate("rect", xmin=r$start, xmax=r$end, ymin=-3, ymax=7,
             fill=col, alpha=0.05) +
    coord_cartesian(ylim=c(-3,5)) +
    labs(x=NULL, y="Log2R", title=r$label) +
    scale_x_continuous(labels=function(x) paste0(round(x/1e6,0),"Mb")) +
    theme_classic(base_size=9) +
    theme(plot.title=element_text(face="bold", size=9))
}

logr_panels <- lapply(REGIONS, make_logr_panel)

# ============================================================================
# 11. GENOME-WIDE WES PLOT  (supplementary only)
# ============================================================================
chr_lengths <- cnr %>%
  dplyr::filter(chromosome %in% CHR_ORDER) %>%
  dplyr::group_by(chromosome) %>%
  dplyr::summarise(max_pos=max(end), .groups="drop") %>%
  dplyr::arrange(factor(chromosome, levels=CHR_ORDER)) %>%
  as.data.frame()
chr_lengths$offset  <- c(0, cumsum(as.numeric(chr_lengths$max_pos[-nrow(chr_lengths)])))
chr_lengths$mid_pos <- chr_lengths$offset + chr_lengths$max_pos/2

cnr2 <- merge(cnr, chr_lengths[,c("chromosome","offset")], by="chromosome")
cnr2$genome_pos <- cnr2$midpos + as.numeric(cnr2$offset)

highlight <- do.call(rbind, lapply(REGIONS, function(r) {
  ev  <- evidence[evidence$region == r$label, ]
  off <- chr_lengths$offset[chr_lengths$chromosome == r$chrom]
  data.frame(genome_start=r$start+off, genome_end=r$end+off,
             label=r$label, higher_subclone=ev$higher_subclone,
             stringsAsFactors=FALSE)
}))
highlight$color <- subclone_colors[highlight$higher_subclone]

p_genome <- ggplot() +
  geom_point(data=cnr2, aes(x=genome_pos, y=log2, color=chromosome),
             size=0.1, alpha=0.25, show.legend=FALSE) +
  geom_rect(data=highlight,
            aes(xmin=genome_start, xmax=genome_end, ymin=-4, ymax=4),
            fill=highlight$color, alpha=0.2, inherit.aes=FALSE) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=0.5) +
  geom_text(data=highlight,
            aes(x=(genome_start+genome_end)/2, y=3.5, label=label),
            color=highlight$color, size=2.8, fontface="bold") +
  scale_x_continuous(breaks=chr_lengths$mid_pos,
                     labels=gsub("chr","",chr_lengths$chromosome)) +
  coord_cartesian(ylim=c(-4,4)) +
  labs(x=NULL, y="Log2R") +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(size=7))

# ============================================================================
# 12. ASSEMBLE & SAVE
# ============================================================================

# --- Fig 1F: fishplot / VAF / scaled expression ---
p_fish_1f <- make_fish_grob(width=2400, height=900)

fig1F <- p_fish_1f /
  wrap_plots(vaf_panels,  ncol=3) /
  wrap_plots(expr_panels, ncol=3) +
  plot_layout(heights=c(1.2, 1.2, 1.0))

save_figure(fig1F, 'Fig1F_clonal_expression')

# --- Fig S1D: genome-wide logR / per-region logR+VAF / fishplot ---
p_fish_s1d <- make_fish_grob(width=2400, height=900)

figS1D <- p_genome /
  wrap_plots(logr_panels, ncol=3) /
  wrap_plots(vaf_panels,  ncol=3) /
  p_fish_s1d +
  plot_layout(heights=c(1.0, 1.2, 1.2, 1.0))

save_figure(figS1D, 'FigS1D_WES_clonal')
