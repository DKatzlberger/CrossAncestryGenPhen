#' Filter methylation data by ancestry and phenotype
#'
#' Combines two ancestry-specific matrices (`X`, `Y`) with their metadata (`MX`, `MY`),
#' filters CpGs **only by variance** across samples, optionally within
#' ancestry × phenotype groups. Never uses mean beta-values for filtering.
#'
#' Supports only 2×2 designs (two levels in `g_col` × two levels in `a_col`).
#'
#' @param X Matrix of beta-values (rows = samples, cols = CpGs) for ancestry X.
#' @param Y Matrix of beta-values for ancestry Y.
#' @param MX Metadata for X (data.frame; rownames must match rownames of X).
#' @param MY Metadata for Y (data.frame; rownames must match rownames of Y).
#' @param g_col Name of phenotype column in metadata (factor with 2 levels).
#' @param a_col Name of ancestry column in metadata (factor/character, 1 level per block).
#' @param any_group Logical; if TRUE, filter variance within ancestry×phenotype groups.
#'        If FALSE, treat all samples as one population (global variance filter).
#' @param var_q Lower quantile cutoff for variance (default = 0.20).
#' @param verbose Print summary (default TRUE).
#' @param plot Print variance plot (default TRUE).
#'
#' @export
filter_by_methylation <- function(
  X,
  Y,
  MX,
  MY,
  g_col,
  a_col,
  any_group = TRUE,
  var_q = 0.20,
  verbose = TRUE,
  plot = TRUE
){

  ## --- Input data checks ---
  assert_input(
    X = X, 
    Y = Y, 
    MX = MX, 
    MY = MY, 
    g_col = g_col, 
    a_col = a_col,
    .fun = "filter_by_methylation"
  )

  ## --- Factor levels ---
  a_1 <- unique(MX[[a_col]])
  a_2 <- unique(MY[[a_col]])

  g_levels <- levels(MX[[g_col]])
  a_levels <- c(a_1, a_2)

  if (length(g_levels) != 2 || length(a_levels) != 2) {
    stop("[filter_by_methylation] Function supports only 2×2 designs.")
  }


  ## --- Bind data ----
  matr <- rbind(X, Y)
  meta <- rbind(MX, MY)

  if (!identical(rownames(matr), rownames(meta))) {
    stop("[filter_by_methylation] Combined methylation/meta rownames must match exactly.")
  }


  ## --- Global variance ---
  global_var <- apply(matr, 2, var, na.rm = TRUE)
  global_var[is.na(global_var)] <- 0


  ## --- Validate grouping ---
  grp <- interaction(meta[[g_col]], meta[[a_col]], drop = TRUE)
  if (any_group) {
    group_sizes <- table(grp)
    if (any(group_sizes < 2)) stop("[filter_by_methylation] any_group = TRUE requires ≥2 samples per ancestry×phenotype group.")
  }
    

  ## --- Filtering ---
  if (any_group) {

    ## group-specific filtering
    keep_mat <- sapply(levels(grp), function(g) {
      idx <- grp == g

      ## group variance (no fallback, guaranteed ≥2 samples)
      v <- apply(matr[idx, , drop = FALSE], 2, var, na.rm = TRUE)
      v[is.na(v)] <- 0  # only for numerical cleanup

      v_thr <- quantile(v, probs = var_q, na.rm = TRUE)
      v > v_thr
    })

    if (is.vector(keep_mat)) keep_mat <- matrix(keep_mat, ncol = 1)

    keep <- rowSums(keep_mat) > 0

  } else {

    ## global filtering
    v <- global_var
    v_thr <- quantile(v, probs = var_q, na.rm = TRUE)
    keep <- v > v_thr
  }

  ## --- Apply filtering ---
  matr_filt <- matr[, keep, drop = FALSE]
  n_features <- sum(keep)


  ## --- Mean–variance plot ----
  m_filt <- colMeans(matr_filt, na.rm = TRUE)
  v_filt <- apply(matr_filt, 2, var, na.rm = TRUE)
  df <- data.frame(mean = m_filt, var = v_filt)


  p <- ggplot(df, aes(x = mean, y = log2(var))) +
    geom_point(size = 0.5) +
    theme_CrossAncestryGenPhen() +
    theme(legend.position = "none") +
    labs(
      x = "Beta-value",
      y = "Log2 Variance",
    )

  if (plot) print(p)


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


  ## --- Return ---
  return(out)
}