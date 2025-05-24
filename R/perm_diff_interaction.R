#' Permutation-Based Interaction Test (Adaptive or Fixed)
#'
#' Performs a permutation test for interaction between conditions across ancestries.
#' Supports both a fixed number of permutations (`B`) and an adaptive mode (`B = NULL`),
#' where permutations are run until empirical p-values converge.
#'
#' The input expression matrices must have **samples as rows** and **genes as columns**.
#' Meta-data must contain consistent group and ancestry annotations for each sample.
#'
#' @param X A numeric matrix of expression values for ancestry A.
#'          Rows = samples, columns = genes.
#' @param Y A numeric matrix of expression values for ancestry B.
#'          Rows = samples, columns = genes.
#' @param MX A data.frame of metadata for X. Must include group and ancestry info.
#' @param MY A data.frame of metadata for Y. Must include group and ancestry info.
#' @param g_col Name of the column in metadata indicating the condition/group.
#' @param a_col Name of the column in metadata indicating ancestry.
#' @param B Integer. Number of permutations (set to NULL to enable adaptive mode).
#' @param seed Optional random seed for reproducibility.
#' @param min_iter Minimum number of permutations to run before checking convergence (adaptive mode).
#' @param max_iter Maximum number of permutations allowed (adaptive mode).
#' @param tol Tolerance threshold for convergence of empirical p-values (adaptive mode).
#' @param batch_size Number of permutations to add per iteration when checking convergence.
#' @param check_convergence Logical. Whether to check convergence in fixed-B mode.
#'
#' @return A list containing:
#' \describe{
#'   \item{summary_stats}{A data.frame with one row per gene, containing: \code{T_obs}, \code{p_value}, \code{ci_lower}, and \code{ci_upper}.}
#'   \item{T_boot}{Matrix of permutation test statistics. One row per permutation.}
#'   \item{B_used}{Number of valid permutations actually used.}
#'   \item{converged}{Logical indicating whether convergence was reached (adaptive or fixed mode).}
#' }
#'
#' @examples
#' set.seed(1)
#' n_per_group <- 25
#' p <- 2000
#'
#' X <- matrix(rnorm(n_per_group * 2 * p), nrow = n_per_group * 2, ncol = p)
#' Y <- matrix(rnorm(n_per_group * 2 * p), nrow = n_per_group * 2, ncol = p)
#' colnames(X) <- colnames(Y) <- paste0("Gene_", seq_len(p))
#'
#' MX <- data.frame(
#'   condition = factor(rep(c("A", "B"), each = n_per_group)),
#'   ancestry = "EUR"
#' )
#' MY <- data.frame(
#'   condition = factor(rep(c("A", "B"), each = n_per_group)),
#'   ancestry = "AFR"
#' )
#'
#' result <- perm_diff_interaction(
#'   X = X,
#'   Y = Y,
#'   MX = MX,
#'   MY = MY,
#'   g_col = "condition",
#'   a_col = "ancestry",
#'   B = 100,
#'   seed = 42
#' )
#'
#' head(result$summary_stats)
#'
#' @importFrom stats colMeans p.adjust quantile complete.cases
#' @importFrom data.table data.table
#' @export

perm_diff_interaction <- function(
  X,
  Y,
  MX,
  MY,
  g_col,
  a_col,
  B = NULL,
  seed = NULL,
  min_iter = 500,
  max_iter = 50000,
  tol = 1e-3,
  batch_size = 100,
  check_convergence = TRUE
) {
  if (!is.null(seed)) set.seed(seed)

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
  T_obs <- delta_X - delta_Y

  # --- Combine data ---
  XY <- rbind(X, Y)
  g_XY <- c(as.character(g_X), as.character(g_Y))

  n1 <- nrow(X)
  n_total <- nrow(XY)
  n_feat <- ncol(XY)
  max_B <- ifelse(is.null(B), max_iter, B)

  T_boot <- matrix(NA_real_, nrow = max_B, ncol = n_feat)
  colnames(T_boot) <- colnames(XY)

  converged <- NA
  B_used <- 0
  p_vals_prev <- rep(NA_real_, n_feat)

  if (!is.null(B)) {
    # --- Fixed B Mode ---
    if (check_convergence) {
      batch_iters <- ceiling(B / batch_size)
      converged <- FALSE

      for (i in seq_len(batch_iters)) {
        for (j in seq_len(batch_size)) {
          b <- (i - 1) * batch_size + j
          if (b > B) break
          idx <- permute_indices(n_total, n1)
          T_boot[b, ] <- compute_T_stat(XY, g_XY, idx$p1, idx$p2, g1, g2)
        }

        current_rows <- 1:min(i * batch_size, B)
        T_valid <- T_boot[current_rows, , drop = FALSE]
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
        T_boot[b, ] <- compute_T_stat(XY, g_XY, idx$p1, idx$p2, g1, g2)
      }
    }

    T_boot <- T_boot[complete.cases(T_boot), , drop = FALSE]
    B_used <- nrow(T_boot)
  } else {
    # --- Adaptive Mode ---
    converged <- FALSE
    while (!converged && B_used < max_iter) {
      for (i in seq_len(batch_size)) {
        B_used <- B_used + 1
        idx <- permute_indices(n_total, n1)
        T_boot[B_used, ] <- compute_T_stat(XY, g_XY, idx$p1, idx$p2, g1, g2)
      }

      T_valid <- T_boot[1:B_used, , drop = FALSE]
      T_valid <- T_valid[complete.cases(T_valid), , drop = FALSE]
      B_valid <- nrow(T_valid)

      p_vals <- compute_empirical_pvalues(T_valid, T_obs)

      if (!anyNA(p_vals_prev)) {
        delta <- max(abs(p_vals - p_vals_prev))
        converged <- B_valid >= min_iter && delta < tol
      }

      p_vals_prev <- p_vals
    }

    T_boot <- T_boot[1:B_used, , drop = FALSE]
    T_boot <- T_boot[complete.cases(T_boot), , drop = FALSE]
    B_used <- nrow(T_boot)
  }

  # --- Check final validity ---
  if (B_used == 0) {
    stop("All permutations failed due to empty subgroup combinations. Consider increasing sample size.")
  }

  # --- Compute p-values and CI ---
  p_vals <- compute_empirical_pvalues(T_boot, T_obs)
  ci_bounds <- compute_CI_bounds(T_boot)

  # --- Return results ---
  return(list(
    summary_stats = data.frame(
      feature = colnames(XY),
      T_obs = T_obs,
      ci_lower = ci_bounds[1, ],
      ci_upper = ci_bounds[2, ],
      p_value = p_vals,
      row.names = NULL 
    ),
    T_boot = T_boot,
    B_used = B_used,
    converged = converged
  ))
}
