#' Filter expression data by ancestry and phenotype
#'
#' Combines two ancestry-specific matrices (`X`, `Y`) with their metadata (`MX`, `MY`),
#' filters genes using `edgeR::filterByExpr` with groups defined by
#' `interaction(meta[[a_col]], meta[[g_col]])`, and re-splits the filtered data
#' using `subset_phenotype_ancestry()`. Supports only 2×2 designs
#' (two ancestry levels × two phenotype levels).
#'
#' @param X Matrix of counts (rows = samples, cols = features) for ancestry X.
#' @param Y Matrix of counts for ancestry Y.
#' @param MX Metadata for X (data.frame, rownames = sample IDs).
#' @param MY Metadata for Y (data.frame, rownames = sample IDs).
#' @param g_col Name of phenotype column in metadata (factor with 2 levels).
#' @param a_col Name of ancestry column in metadata (1 level per block).
#' @param any_group Logical; if TRUE (default), filter by expression within groups. If FALSE, filter without grouping.
#' @param plot Return mean–variance trend plot if available (default TRUE).
#' @param verbose Print summary (default TRUE).
#'
#' @export
filter_by_expression <- function(
  X,
  Y,
  MX,
  MY,
  g_col,
  a_col,
  any_group = TRUE,
  plot = TRUE,
  verbose = TRUE
){

  ## --- Input data structure check ---
  assert_input(
    X = X, 
    Y = Y, 
    MX = MX, 
    MY = MY, 
    g_col = g_col, 
    a_col = a_col,
    .fun = "filter_by_expression"
  )


  ## --- Factor levels ----
  a_1 <- unique(MX[[a_col]])
  a_2 <- unique(MY[[a_col]])

  g_levels <- levels(MX[[g_col]])
  a_levels <- c(a_1, a_2)

  if (length(g_levels) != 2 || length(a_levels) != 2) {
    stop("[filter_by_expression] Function supports only 2x2 designs (two levels in g_col × two levels a_col).")
  }


  ## --- Bind data ----
  matr <- rbind(X, Y)
  meta <- rbind(MX, MY)

  if (!identical(rownames(matr), rownames(meta))) {
    stop("[filter_by_expression] Combined counts/meta rownames must match exactly.")
  }


  ## --- edgeR filterByExpr ---
  dge <- edgeR::DGEList(counts = t(matr))

  if (any_group) {
    grp <- interaction(meta[[g_col]], meta[[a_col]], drop = TRUE)

    if (length(grp) != ncol(dge)) {
      stop("[filter_by_expression] Length of group != number of samples in DGEList.")
    }

    keep <- edgeR::filterByExpr(dge, group = grp)

  } else {
    # No grouping — treat all samples as one population
    grp <- factor(rep("All", ncol(dge)))
    keep <- edgeR::filterByExpr(dge, group = grp)
  }

  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  matr_filt <- t(dge$counts)
  n_features  <- sum(keep)


  ## --- MV-trend plot ---
  p_trend <- plot_mean_variance_trend(
    X = matr_filt,
    point_size = 0.5
  )

  p_density <- plot_mean_variance_density(
    X = matr_filt,
    point_size = 0.5
  ) 

  # Patchwork
  p <- p_trend + p_density
  if (plot){
    print(p)
  }


  ## --- Resplit into two ancestries ---
  out <- filter_phenotype_ancestry(
    X = matr_filt,
    M = meta,
    g_col = g_col,
    a_col = a_col,
    g_levels = g_levels,
    a_levels = a_levels,
    plot = FALSE,
    verbose = FALSE
  )

  # Validation 
  if (ncol(out$X$matr) != ncol(out$Y$matr)){
    stop("[filter_by_expression] Number of features do not align between X and Y")
  }

  # Add the plot
  out$plot <- p


  ## --- Verbose message ----
  if (verbose) {

    fmt_counts <- function(M, g_col) {
      if (is.null(M) || !nrow(M)) return("")
      tab <- table(M[[g_col]], dnn = NULL)
      paste(sprintf("%s: %-4d", names(tab), as.integer(tab)), collapse = "  ")
    }

    message("\nFilter by expression summary:")
    if (any_group) {
      message(sprintf("%-13s %s", "Groups:", paste(unique(grp), collapse = "  ")))
    } else {
      message(sprintf("%-13s %s", "Groups:", "1"))
    }

    message(sprintf(
     "Ancestry (X): %-10s N: %-5d  %-18s  features: %-5d",
      out$X$ancestry,
      nrow(out$X$meta),
      fmt_counts(out$X$meta, g_col),
      n_features
    ))

    message(sprintf(
      "Ancestry (Y): %-10s N: %-5d  %-18s  features: %-5d",
      out$Y$ancestry,
      nrow(out$Y$meta),
      fmt_counts(out$Y$meta, g_col),
      n_features
    ))
  }

  # if (verbose) {
  #   fmt_counts <- function(M, g_col) {
  #     if (is.null(M) || !nrow(M)) return("")
  #     tab <- table(M[[g_col]], dnn = NULL)
  #     paste(sprintf("%s: %-4d", names(tab), as.integer(tab)), collapse = " ")
  #   }


  #   message("\nFilter by expression summary:")
  #   if (any_group){
  #     message(sprintf("Groups:     %s", paste(unique(grp), collapse = "  ")))
  #   } else{
  #      message(sprintf("Groups:     %s", "1"))
  #   }
  #   message(sprintf("%s (X):    N: %-4d %s features: %-4d", out$X$ancestry, nrow(out$X$meta), fmt_counts(out$X$meta, g_col), n_features))
  #   message(sprintf("%s (Y):    N: %-4d %s features: %-4d", out$Y$ancestry, nrow(out$Y$meta), fmt_counts(out$Y$meta, g_col), n_features))
  # }


  ## --- Return ---
  return(out)
}