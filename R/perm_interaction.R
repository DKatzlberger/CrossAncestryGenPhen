#' Permutation-Based Interaction Test with Mandatory Bootstrapped Inference
#'
#' Performs a permutation test (optional) for interaction effects between conditions across ancestries.
#' Always estimates bootstrapped standard errors and confidence intervals for the observed interaction effects
#' via resampling within ancestries.
#'
#' Permutation testing can be toggled using the `permute` flag. When `permute = FALSE`, only bootstrapped
#' confidence intervals and standard errors are provided, and p-values are returned as NA.
#'
#' The input expression matrices must have \strong{samples as rows} and \strong{genes as columns}.
#'
#' @param X A numeric matrix of expression values for ancestry A (rows = samples, columns = genes).
#' @param Y A numeric matrix of expression values for ancestry B (same format as \code{X}).
#' @param MX A data.frame of metadata for \code{X}, must include columns for group and ancestry.
#' @param MY A data.frame of metadata for \code{Y}, must include columns for group and ancestry.
#' @param g_col Column name in metadata indicating condition/group (factor with 2 levels).
#' @param a_col Column name in metadata indicating ancestry (must be unique per ancestry).
#' @param B Integer. Number of permutations and bootstraps to perform. Set to \code{NULL} for adaptive convergence.
#' @param permute Logical. Whether to perform permutation-based interaction testing.
#' @param seed Optional integer to set random seed for reproducibility.
#' @param min_iter Minimum number of permutations/bootstraps to run before checking convergence.
#' @param max_iter Maximum number of permutations/bootstraps allowed.
#' @param tol Tolerance threshold for convergence (applies to empirical p-values and bootstrap SEs).
#' @param batch_size Number of permutations/bootstraps to add per iteration.
#' @param check_convergence Logical. Whether to check convergence in fixed-\code{B} mode (if \code{permute = TRUE}).
#'
#' @return A list with:
#' \describe{
#'   \item{summary_stats}{A data.frame with gene-wise observed statistics (\code{T_obs}),
#'                        bootstrapped standard errors (\code{SE}),
#'                        95\% confidence intervals (\code{CI_lower}, \code{CI_upper}),
#'                        empirical p-values (\code{p_value}, NA if \code{permute = FALSE}),
#'                        and Benjamini-Hochberg adjusted p-values (\code{p_adj}, NA if \code{permute = FALSE}).}
#'   \item{T_perm}{Matrix of permutation-based test statistics (or \code{NULL} if \code{permute = FALSE}).}
#'   \item{T_boot}{Matrix of bootstrap-based observed statistics (rows = bootstraps, columns = genes).}
#'   \item{B_used_perm}{Number of valid permutations used (or \code{NA} if \code{permute = FALSE}).}
#'   \item{B_used_boot}{Number of valid bootstrap replicates used.}
#'   \item{converged_perm}{Logical indicating if permutation p-values converged (or \code{NA} if \code{permute = FALSE}).}
#'   \item{converged_boot}{Logical indicating if bootstrap SEs converged.}
#' }
#'
#' @importFrom stats p.adjust quantile complete.cases
#' @export

perm_interaction <- function(
  X,
  Y,
  MX,
  MY,
  g_col,
  a_col,
  B = 1000,
  seed = NULL,
  permute = TRUE,
  min_iter = 500,
  max_iter = 10000,
  tol = 1e-3,
  batch_size = 100,
  check_convergence = TRUE
) {
  if (!is.null(seed)) set.seed(seed)

  stopifnot(is.matrix(X), is.matrix(Y))
  stopifnot(ncol(X) == ncol(Y))

  # --- Validate inputs ---
  g_X <- MX[[g_col]]
  g_Y <- MY[[g_col]]
  validate_groups(g_X, g_Y)
  g1 <- levels(g_X)[1]
  g2 <- levels(g_X)[2]

  a_X <- unique(MX[[a_col]])
  a_Y <- unique(MY[[a_col]])
  validate_ancestry(a_X, a_Y)

  # --- Observed statistic ---
  delta_X <- mean_diff_by_group(X, g_X, g1, g2)
  delta_Y <- mean_diff_by_group(Y, g_Y, g1, g2)
  T_obs <- delta_Y - delta_X

  # --- Bootstrapped SE and CI for observed ---
  boot_res <- bootstrap_T_obs(
    X = X, 
    Y = Y, 
    MX = MX, 
    MY = MY, 
    g_col = g_col, 
    a_col = a_col, 
    g1 = g1, 
    g2 = g2, 
    B = B, 
    seed = if (!is.null(seed)) seed + 1000 else NULL
  )

  # --- Initialize default outputs for permutation ---
  T_perm <- NULL
  p_vals <- rep(NA_real_, ncol(X))
  p_adj_vals <- rep(NA_real_, ncol(X))
  converged <- NA
  B_used <- NA_integer_

  if (permute) {
    # --- Combine for permutation ---
    XY <- rbind(X, Y)
    g_XY <- c(as.character(g_X), as.character(g_Y))
    n1 <- nrow(X)
    n_total <- nrow(XY)
    n_feat <- ncol(XY)
    max_B <- ifelse(is.null(B), max_iter, B)

    T_perm <- matrix(NA_real_, nrow = max_B, ncol = n_feat)
    colnames(T_perm) <- colnames(XY)

    p_vals_prev <- rep(NA_real_, n_feat)
    converged <- FALSE
    B_used <- 0

    if (!is.null(B)) {
      if (check_convergence) {
        batch_iters <- ceiling(B / batch_size)
        for (i in seq_len(batch_iters)) {
          for (j in seq_len(batch_size)) {
            b <- (i - 1) * batch_size + j
            if (b > B) break
            idx <- permute_indices(n_total, n1)
            T_perm[b, ] <- compute_T_stat(XY, g_XY, idx$p1, idx$p2, g1, g2)
          }

          current_rows <- 1:min(i * batch_size, B)
          T_valid <- T_perm[current_rows, , drop = FALSE]
          T_valid <- T_valid[complete.cases(T_valid), , drop = FALSE]

          if (nrow(T_valid) > 10) {
            p_vals <- compute_empirical_pvalues(T_valid, T_obs)
            if (!anyNA(p_vals_prev)) {
              delta <- max(abs(p_vals - p_vals_prev))
              if (delta < tol && nrow(T_valid) >= min_iter) {
                converged <- TRUE
                break
              }
            }
            p_vals_prev <- p_vals
          }
        }
      } else {
        for (b in seq_len(B)) {
          idx <- permute_indices(n_total, n1)
          T_perm[b, ] <- compute_T_stat(XY, g_XY, idx$p1, idx$p2, g1, g2)
        }
      }

      T_perm <- T_perm[complete.cases(T_perm), , drop = FALSE]
      B_used <- nrow(T_perm)
    } else {
      while (!converged && B_used < max_iter) {
        for (i in seq_len(batch_size)) {
          B_used <- B_used + 1
          idx <- permute_indices(n_total, n1)
          T_perm[B_used, ] <- compute_T_stat(XY, g_XY, idx$p1, idx$p2, g1, g2)
        }

        T_valid <- T_perm[1:B_used, , drop = FALSE]
        T_valid <- T_valid[complete.cases(T_valid), , drop = FALSE]

        p_vals <- compute_empirical_pvalues(T_valid, T_obs)

        if (!anyNA(p_vals_prev)) {
          delta <- max(abs(p_vals - p_vals_prev))
          converged <- B_used >= min_iter && delta < tol
        }

        p_vals_prev <- p_vals
      }

      T_perm <- T_perm[1:B_used, , drop = FALSE]
      T_perm <- T_perm[complete.cases(T_perm), , drop = FALSE]
    }

    if (B_used == 0) stop("All permutations failed due to empty subgroup combinations.")
    p_vals <- compute_empirical_pvalues(T_perm, T_obs)
    p_adj_vals <- p.adjust(p_vals, method = "BH")
  }

  # --- Assemble output ---
  summary_stats <- data.frame(
    feature = colnames(X),
    T_obs = T_obs,
    SE = boot_res$summary_stats$SE,
    CI_lower = boot_res$summary_stats$CI_lower,
    CI_upper = boot_res$summary_stats$CI_upper,
    p_value = p_vals,
    p_adj = p_adj_vals,
    row.names = NULL
  )

  res <- list(
    summary_stats = summary_stats,
    T_perm = T_perm,
    T_boot = boot_res$T_boot,
    B_used_perm = B_used,
    B_used_boot = boot_res$B_used,
    converged_perm = converged,
    converged_boot = boot_res$converged
  )

  return(res)
}
