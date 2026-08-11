# =============================================================================
# R/utils/functions.R — shared helpers for the T47D fate manuscript pipeline
# =============================================================================
# Source this file from each analysis script (after sourcing config.R).
# The save_* helpers look for `fig_dir` / `table_dir` in the calling environment
# when `dir` is not passed explicitly, so scripts only need to set those two
# variables and can call save_figure(p, name) without repeating the path.

# -----------------------------------------------------------------------------
# DE helpers
# -----------------------------------------------------------------------------

filter_de_hits <- function(de_tbl,
                           padj_max        = 0.05,
                           min_log2fc      = 1,
                           min.pct         = 0.3,
                           ratio_threshold = 2,
                           eps             = 1e-3) {
  de_tbl %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(
      p_val_adj < padj_max,
      pct.1 > min.pct | pct.2 > min.pct,
      abs(avg_log2FC) > min_log2fc,
      (pct.1 / (pct.2 + eps) > ratio_threshold) |
        (pct.2 / (pct.1 + eps) > ratio_threshold)
    ) %>%
    tibble::column_to_rownames("gene")
}

# -----------------------------------------------------------------------------
# Figure helpers
# -----------------------------------------------------------------------------

save_figure <- function(p, name,
                        dir    = get0("fig_dir", envir = parent.frame(), inherits = TRUE),
                        width  = 7, height = 5, dpi = 300) {
  if (is.null(dir)) stop("set `fig_dir` or pass `dir` explicitly to save_figure().")
  ggsave(file.path(dir, paste0(name, ".pdf")), plot = p,
         width = width, height = height, device = "pdf")
  ggsave(file.path(dir, paste0(name, ".png")), plot = p,
         width = width, height = height, dpi = dpi)
  message("Saved: ", file.path(dir, paste0(name, ".pdf")))
  invisible(p)
}


save_pheatmap <- function(ph, name,
                          dir    = get0("fig_dir", envir = parent.frame(), inherits = TRUE),
                          width  = 10, height = 8) {
  if (is.null(dir)) stop("set `fig_dir` or pass `dir` explicitly to save_pheatmap().")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, paste0(name, ".pdf"))
  
  pdf(path, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  
  grid::grid.newpage()
  if (inherits(ph, c("Heatmap", "HeatmapList"))) {
    ComplexHeatmap::draw(ph, newpage = FALSE)
  } else if (!is.null(ph$gtable)) {
    grid::grid.draw(ph$gtable)
  } else {
    stop("unrecognised heatmap object: ", paste(class(ph), collapse = "/"))
  }
  
  message("Saved: ", path)
  invisible(ph)
}

# eulerr / base-graphics objects (cannot be passed to ggsave)
save_euler <- function(p, name,
                       dir    = get0("fig_dir", envir = parent.frame(), inherits = TRUE),
                       width  = 5, height = 5, dpi = 300) {
  if (is.null(dir)) stop("set `fig_dir` or pass `dir` explicitly to save_euler().")
  pdf_path <- file.path(dir, paste0(name, ".pdf"))
  png_path <- file.path(dir, paste0(name, ".png"))
  pdf(pdf_path, width = width, height = height); print(p); dev.off()
  png(png_path, width = width, height = height, units = "in", res = dpi); print(p); dev.off()
  message("Saved: ", pdf_path)
  invisible(p)
}

# -----------------------------------------------------------------------------
# Table helper
# -----------------------------------------------------------------------------

save_table <- function(df, name,
                       dir        = get0("table_dir", envir = parent.frame(), inherits = TRUE),
                       sheet_name = "Sheet1",
                       row.names  = FALSE) {
  if (is.null(dir)) stop("set `table_dir` or pass `dir` explicitly to save_table().")
  csv_path <- file.path(dir, paste0(name, ".csv"))
  write.csv(df, csv_path, row.names = row.names, quote = TRUE)
  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(setNames(list(df), sheet_name),
                        file.path(dir, paste0(name, ".xlsx")))
  }
  message("Saved: ", csv_path)
  invisible(df)
}


# -----------------------------------------------------------------------------
# run_qc_and_doublet_removal()
#
# Accepts a named list of Seurat objects (one per sample/batch).
# Each object must have percent.mt already computed before calling this.
#
# Per-sample steps:
#   1. MAD-based QC  (linear scale, median ± 3 MAD; counts/features floored at 0)
#   2. Normalise → HVG → Scale → PCA → UMAP → Neighbours → Clusters
#   3. DoubletFinder parameter sweep → classify → keep Singlets only
#
# Returns the cleaned list; caller is responsible for merging.
#
# Args:
#   sample.list        Named list of Seurat objects
#   pcs                PC dims for reductions + DoubletFinder  (default 1:20)
#   cluster_resolution FindClusters resolution                 (default 0.5)
#   doublet_rate       Expected doublet fraction               (default 0.008)
#   pN                 DoubletFinder pN                        (default 0.25)
# -----------------------------------------------------------------------------
run_qc_and_doublet_removal <- function(sample.list,
                                       pcs                = 1:20,
                                       cluster_resolution = 0.5,
                                       doublet_rate       = 0.008,
                                       pN                 = 0.25) {
  
  sample.list <- lapply(sample.list, function(x) {
    
    # ── 1. MAD QC ─────────────────────────────────────────────────────────────
    max.feat  <- median(x$nFeature_RNA) + 3 * mad(x$nFeature_RNA)
    min.feat  <- max(median(x$nFeature_RNA) - 3 * mad(x$nFeature_RNA), 0)
    max.count <- median(x$nCount_RNA)   + 3 * mad(x$nCount_RNA)
    min.count <- max(median(x$nCount_RNA) - 3 * mad(x$nCount_RNA), 0)
    max.mito  <- max(median(x$percent.mt)   + 3 * mad(x$percent.mt), 25)
    min_mito <- min(median(x$percent.mt) - 3 * mad(x$percent.mt), 5)
    
    n_before <- ncol(x)
    x <- subset(x, subset =
                  nFeature_RNA > min.feat  & nFeature_RNA < max.feat  &
                  nCount_RNA   > min.count & nCount_RNA   < max.count &
                  percent.mt   < max.mito)
    
    cat(sprintf(
      "  [QC] %-30s | Kept %d / %d cells | nCount (%.0f, %.0f)  nFeat (%.0f, %.0f)  MT < %.1f%%\n",
      Project(x), ncol(x), n_before,
      min.count, max.count, min.feat, max.feat, max.mito
    ))
    
    # Guard: DoubletFinder needs enough cells to be meaningful
    if (ncol(x) < 100) {
      warning(sprintf("  [SKIP DoubletFinder] %s: only %d cells after QC.",
                      Project(x), ncol(x)))
      return(x)
    }
    
    # ── 2. Preprocessing ───────────────────────────────────────────────────────
    x <- NormalizeData(x,        verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
    x <- ScaleData(x,            verbose = FALSE)
    x <- RunPCA(x,               verbose = FALSE)
    x <- RunUMAP(x,    dims = pcs, verbose = FALSE)
    x <- FindNeighbors(x, dims = pcs, verbose = FALSE)
    x <- FindClusters(x,  resolution = cluster_resolution, verbose = FALSE)
    
    # ── 3. DoubletFinder ──────────────────────────────────────────────────────
    x <- tryCatch({
      sweep.res  <- paramSweep(x, PCs = pcs, sct = FALSE)
      sweep.stat <- summarizeSweep(sweep.res, GT = FALSE)
      bcmvn      <- find.pK(sweep.stat)
      pK_val     <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
      nExp_poi   <- round(doublet_rate * ncol(x))
      
      x <- doubletFinder(x, PCs = pcs, pN = pN, pK = pK_val,
                         nExp = nExp_poi, reuse.pANN = FALSE)
      df_col <- grep("DF.classifications", colnames(x@meta.data), value = TRUE)
      # Use the last matching column (safest when reuse.pANN columns accumulate)
      df_col <- df_col[length(df_col)]
      
      n_before_df <- ncol(x)
      x <- subset(x, cells = colnames(x)[x@meta.data[[df_col]] == "Singlet"])
      cat(sprintf("  [DoubletFinder] %-30s | Removed %d doublets → %d singlets\n",
                  Project(x), n_before_df - ncol(x), ncol(x)))
      x
    }, error = function(e) {
      warning(sprintf("  [DoubletFinder ERROR] %s: %s — returning QC-filtered object.",
                      Project(x), conditionMessage(e)))
      x
    })
    
    return(x)
  })
  
  return(sample.list)
}


# -----------------------------------------------------------------------------
# nmf_cluster_cells — cross-replicate robust NMF metaprograms
# -----------------------------------------------------------------------------
# ARGUMENTS:
#   seu               Merged Seurat object (RNA assay, "data" layer)
#   split_by          Metadata column to split into per-sample objects
#   group_col         Metadata column defining condition/replicate groups
#   k_values          NMF ranks to sweep (default: 3:10)
#   top_n_genes       Top genes per NMF component (default: 60)
#   nmf_maxit         NMF iteration limit (default: 500)
#   overlap_threshold Minimum gene overlap to call two programs equivalent (default: 40)
#   similarity_thr    Jaccard threshold for redundancy removal (default: 0.8)
#   n_metaprograms    K for cutting the MP dendrogram (default: 4)
#   min_frac          Min fraction of programs a gene must appear in per MP (default: 0.25)
#   gene_cap          Max genes per metaprogram (default: 100)
#   exclude_cc        Exclude cell-cycle genes if cc.genes exists (default: TRUE)
#   seed              Random seed (default: 123)
#
# RETURNS: list($MP_genes, $curated_map, $jaccard_matrix, $hclust)

nmf_cluster_cells <- function(
    seu,
    split_by          = "sample",
    group_col         = "condition",
    k_values          = 3:10,
    top_n_genes       = 60,
    nmf_maxit         = 500,
    overlap_threshold = 40,
    similarity_thr    = 0.8,
    n_metaprograms    = 4,
    min_frac          = 0.25,
    gene_cap          = 100,
    exclude_cc        = TRUE,
    seed              = 123
) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(RcppML)
    library(Matrix)
    library(dplyr)
    library(igraph)
  })

  set.seed(seed)

  DefaultAssay(seu) <- "RNA"
  gene_set <- VariableFeatures(seu)

  if (length(gene_set) == 0)
    stop("VariableFeatures(seu) is empty. Run FindVariableFeatures() first.")

  if (exclude_cc && exists("cc.genes"))
    gene_set <- setdiff(gene_set, c(cc.genes$s.genes, cc.genes$g2m.genes))

  hk_pattern <- "^MT-|^RP[LS]|^GAPDH$|^ACTB$"
  gene_set <- setdiff(gene_set, gene_set[grepl(hk_pattern, gene_set, ignore.case = TRUE)])
  cat("Gene set size:", length(gene_set), "\n")
  if (length(gene_set) < 50) warning("Gene set is very small; NMF stability may be poor.")

  sample_list <- SplitObject(seu, split.by = split_by)

  sample_to_group <- setNames(
    sapply(names(sample_list), function(s) {
      grps <- unique(sample_list[[s]][[group_col, drop = TRUE]])
      if (length(grps) > 1)
        stop("Sample '", s, "' has multiple values in group_col '", group_col, "'.")
      grps
    }),
    names(sample_list)
  )

  cat("Sample → group mapping:\n"); print(sample_to_group)

  .run_nmf_sample <- function(s_obj, sample_id) {
    expr <- GetAssayData(s_obj, layer = "data")
    genes_use <- intersect(gene_set, rownames(expr))
    if (length(genes_use) < 20) {
      warning("Skipping sample ", sample_id, ": too few genes.")
      return(list(programs = list(), weights = list()))
    }
    expr <- as.matrix(expr[genes_use, , drop = FALSE])
    expr[is.na(expr)] <- 0
    programs <- list(); weights <- list()
    for (k in k_values) {
      cat("  NMF:", sample_id, "k =", k, "\n")
      fit <- RcppML::nmf(expr, k = k, maxit = nmf_maxit, tol = 1e-5, verbose = FALSE)
      W <- fit$w; rownames(W) <- genes_use
      for (i in seq_len(k)) {
        pid <- paste0(sample_id, "_K", k, "_P", i)
        programs[[pid]] <- genes_use[order(W[, i], decreasing = TRUE)][seq_len(min(top_n_genes, nrow(W)))]
        wv <- W[, i]; names(wv) <- genes_use
        weights[[pid]] <- wv
      }
    }
    list(programs = programs, weights = weights)
  }

  .robust_within <- function(progs) {
    if (length(progs) == 0) return(list())
    prog_names <- names(progs); keep <- character(0)
    for (i in seq_along(progs)) {
      overlaps <- sapply(progs, function(x) length(intersect(progs[[i]], x)))
      if (sum(overlaps >= overlap_threshold) >= 2) keep <- c(keep, prog_names[i])
    }
    progs[unique(keep)]
  }

  .cross_replicate_overlap <- function(prog, other_samples) {
    if (length(other_samples) == 0) return(0)
    max(sapply(other_samples, function(sp) {
      if (length(sp) == 0) return(0)
      max(sapply(sp, function(x) length(intersect(prog, x))))
    }))
  }

  .reduce_consensus <- function(progs) {
    if (length(progs) <= 1) return(progs)
    progs <- unlist(progs, recursive = FALSE)
    all_el <- unique(unlist(progs))
    mat <- sparseMatrix(
      i    = rep(seq_along(progs), sapply(progs, length)),
      j    = match(unlist(progs), all_el),
      x    = 1,
      dims = c(length(progs), length(all_el)),
      dimnames = list(names(progs), all_el)
    )
    ints <- as.matrix(tcrossprod(mat))
    rs   <- rowSums(mat)
    sim  <- ints / (outer(rs, rs, "+") - ints)
    g    <- igraph::graph_from_adjacency_matrix(sim >= similarity_thr, mode = "undirected", diag = FALSE)
    cl   <- igraph::components(g)$membership
    consensus <- list()
    for (cid in unique(cl)) {
      idx <- which(cl == cid)
      if (length(idx) == 1) {
        consensus[[ names(cl)[idx] ]] <- progs[[ idx ]]
      } else {
        sub_mat     <- mat[idx, , drop = FALSE]
        gene_counts <- colSums(sub_mat)
        min_ct      <- ceiling(length(idx) * 0.5)
        cg          <- names(gene_counts)[gene_counts >= min_ct]
        subg        <- igraph::induced_subgraph(g, vids = idx)
        hub         <- names(which.max(igraph::degree(subg)))[1]
        consensus[[hub]] <- cg
      }
    }
    consensus
  }

  .jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

  .make_mp_genes <- function(mp_progs) {
    n     <- length(mp_progs)
    freq  <- sort(table(unlist(mp_progs)), decreasing = TRUE)
    cands <- names(freq[freq / n >= min_frac])
    if (length(cands) > gene_cap) cands <- cands[seq_len(gene_cap)]
    cat("  Programs:", n, "| genes passing min_frac:", length(cands), "\n")
    cands
  }

  .assign_genes <- function(MPs, MP_genes, weights_flat) {
    all_genes <- unique(unlist(MP_genes))
    mp_names  <- names(MPs)
    assignment <- sapply(all_genes, function(gene) {
      scores <- sapply(mp_names, function(mn) {
        pids <- names(MPs[[mn]])
        mean(sapply(pids, function(p) {
          w <- weights_flat[[p]]
          if (!is.null(w) && gene %in% names(w)) w[[gene]] else 0
        }))
      })
      mp_names[which.max(scores)]
    })
    data.frame(Gene = names(assignment), MP = assignment, stringsAsFactors = FALSE)
  }

  cat("\n=== Running NMF on all samples ===\n")
  nmf_runs     <- lapply(names(sample_list), function(s) .run_nmf_sample(sample_list[[s]], s))
  names(nmf_runs) <- names(sample_list)
  nmf_programs <- lapply(nmf_runs, `[[`, "programs")
  nmf_weights  <- lapply(nmf_runs, `[[`, "weights")

  cat("\n=== Within-sample robustness ===\n")
  robust_within <- lapply(nmf_programs, .robust_within)
  cat("Robust programs per sample:\n"); print(sapply(robust_within, length))

  cat("\n=== Cross-replicate robustness ===\n")
  robust_across <- list()
  for (s in names(robust_within)) {
    grp <- sample_to_group[[s]]
    peer_samples <- names(sample_to_group)[sample_to_group == grp & names(sample_to_group) != s]
    others <- robust_within[peer_samples]
    if (length(others) == 0) {
      message("Sample '", s, "' has no peers — applying single-replicate fallback filter.")
      all_sample_progs <- nmf_programs[[s]]
      n_all            <- length(all_sample_progs)
      robust_across[[s]] <- Filter(function(p) {
        if (n_all == 0) return(FALSE)
        n_overlap <- sum(sapply(all_sample_progs,
                                function(x) length(intersect(p, x)) >= overlap_threshold))
        (n_overlap - 1) / max(n_all - 1, 1) >= 0.4
      }, robust_within[[s]])
    } else {
      robust_across[[s]] <- Filter(
        function(p) .cross_replicate_overlap(p, others) >= overlap_threshold,
        robust_within[[s]]
      )
    }
  }
  cat("Cross-replicate robust programs per sample:\n"); print(sapply(robust_across, length))

  cat("\n=== Removing redundant programs ===\n")
  all_programs <- .reduce_consensus(robust_across)
  cat("Non-redundant programs:", length(all_programs), "\n")
  if (length(all_programs) < 2) stop("Too few programs remain for metaprogram clustering.")

  cat("\n=== Metaprogram clustering ===\n")
  n <- length(all_programs)
  S <- matrix(0, n, n, dimnames = list(names(all_programs), names(all_programs)))
  for (i in seq_len(n)) for (j in i:n) {
    v <- .jaccard(all_programs[[i]], all_programs[[j]])
    S[i, j] <- v; S[j, i] <- v
  }
  hc    <- hclust(as.dist(1 - S), method = "average")
  mp_id <- cutree(hc, k = n_metaprograms)
  cat("\nMetaprogram cluster sizes:\n"); print(table(mp_id))

  cat("\n=== Metaprogram gene sets ===\n")
  weights_flat <- unlist(nmf_weights, recursive = FALSE)
  MPs      <- split(all_programs, mp_id)
  MP_genes <- lapply(MPs, .make_mp_genes)
  curated_map <- .assign_genes(MPs, MP_genes, weights_flat)

  mp_levels      <- sort(unique(curated_map$MP))
  curated_map$MP <- paste0("MP", match(curated_map$MP, mp_levels))
  cat("Total unique genes across metaprograms:", nrow(curated_map), "\n")

  list(MP_genes = MP_genes, curated_map = curated_map, jaccard_matrix = S, hclust = hc)
}
