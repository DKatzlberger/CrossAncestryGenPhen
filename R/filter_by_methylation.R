#' Filter methylation data by ancestry and phenotype
#'
#' Combines two ancestry-specific matrices (`X`, `Y`) with their metadata (`MX`, `MY`),
#' filters genes/probes by mean beta-value and variance across samples,
#' and re-splits the filtered data using `subset_phenotype_ancestry()`.
#' Supports only 2×2 designs (two ancestry levels × two phenotype levels).
#'
#' @param X Matrix of beta-values (rows = samples, cols = features) for ancestry X.
#' @param Y Matrix of beta-values for ancestry Y.
#' @param MX Metadata for X (data.frame, rownames = sample IDs).
#' @param MY Metadata for Y (data.frame, rownames = sample IDs).
#' @param g_col Name of phenotype column in metadata (factor with 2 levels).
#' @param a_col Name of ancestry column in metadata (1 level per block).
#' @param mean_q Quantile cutoffs for mean beta-values (default = c(0.2, 0.9)).
#' @param var_q Lower quantile cutoff for variance (default = 0.2).
#' @param verbose Print summary (default TRUE).
#' @param plot Return mean–variance scatter plot (default TRUE).
#'
#' @export
filter_by_methylation <- function(
  X,
  Y,
  MX,
  MY,
  g_col,
  a_col,
  mean_q = c(0.2, 0.9),
  var_q = 0.20,
  verbose = TRUE,
  plot = TRUE
){

  ## --- Input data structure check ---
  assert_input(
    X = X, 
    Y = Y, 
    MX = MX, 
    MY = MY, 
    g_col = g_col, 
    a_col = a_col
  )


  ## --- Factor levels ----
  a_1 <- unique(MX[[a_col]])
  a_2 <- unique(MY[[a_col]])

  g_levels <- levels(MX[[g_col]])
  a_levels <- c(a_1, a_2)

  if (length(g_levels) != 2 || length(a_levels) != 2) {
    stop("[filter_by_methylation] Function supports only 2x2 designs (two levels in g_col × two levels a_col).")
  }


  ## --- Bind data ----
  matr <- rbind(X, Y)
  meta <- rbind(MX, MY)

  if (!identical(rownames(matr), rownames(meta))) {
    stop("[filter_by_methylation] Combined methylation/meta rownames must match exactly.")
  }

  ## --- Mean–variance filtering ----
  grp <- interaction(meta[[g_col]], meta[[a_col]], drop = TRUE)

  ## --- Filtering per group ---
  keep_mat <- sapply(levels(grp), function(g) {
      idx <- grp == g
      m   <- colMeans(matr[idx, , drop = FALSE], na.rm = TRUE)
      v   <- apply(matr[idx, , drop = FALSE], 2, var, na.rm = TRUE)

      m_low  <- quantile(m, probs = mean_q[1], na.rm = TRUE)
      m_high <- quantile(m, probs = mean_q[2], na.rm = TRUE)
      v_thr  <- quantile(v, probs = var_q, na.rm = TRUE)

      (m > m_low & m < m_high) & (v > v_thr)
    })

  if (is.vector(keep_mat)) keep_mat <- matrix(keep_mat, ncol = 1)

  keep <- rowSums(keep_mat) > 0

  matr_filt  <- matr[, keep, drop = FALSE]
  n_features <- sum(keep)

  ## --- Mean–variance plot ----
  m_filt <- colMeans(matr_filt, na.rm = TRUE)
  v_filt <- apply(matr_filt, 2, var, na.rm = TRUE)
  df <- data.frame(mean = m_filt, var = v_filt)


  p <- ggplot(df, aes(x = mean, y = log2(var))) +
    geom_point(size = 0.5) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend() +
    theme(legend.position = "none") +
    labs(
      x = "Beta-value",
      y = "Log2 Variance",
    )

  if (plot){
    print(p)
  }


  ## --- Resplit into two ancestries ----
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
  if (ncol(out$X$matr) != ncol(out$Y$matr)) {
    stop("[filter_by_methylation] Number of features do not align between X and Y")
  }

  # Add the plot
  out$plot <- p


  ## --- Verbose message ----
  if (verbose) {
    fmt_counts <- function(M, g_col) {
      if (is.null(M) || !nrow(M)) return("")
      tab <- table(M[[g_col]], dnn = NULL)
      paste(sprintf("%s: %-4d", names(tab), as.integer(tab)), collapse = " ")
    }


    message("\nFilter by methylation summary:")
    message(sprintf("%s (X):    N: %-4d %s features: %-4d", out$X$ancestry, nrow(out$X$meta), fmt_counts(out$X$meta, g_col), n_features))
    message(sprintf("%s (Y):    N: %-4d %s features: %-4d", out$Y$ancestry, nrow(out$Y$meta), fmt_counts(out$Y$meta, g_col), n_features))
  }


  ## --- Return ---
  return(out)
}