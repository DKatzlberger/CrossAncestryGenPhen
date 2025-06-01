#' Validate Group Factors
#'
#' Ensures group labels are factors, share levels, and have exactly 2 levels.
#' @param g_X A factor vector representing group labels for dataset X.
#' @param g_Y A factor vector representing group labels for dataset Y.
#' @keywords internal
validate_groups <- function(
  g_X,
  g_Y
) {
  stopifnot(is.factor(g_X), is.factor(g_Y))
  stopifnot(identical(levels(g_X), levels(g_Y)))
  stopifnot(length(levels(g_X)) == 2)
}

#' Validate Ancestry Constraints
#'
#' Ensures each dataset has exactly one ancestry level.
#' @param a_X An ancestry label for dataset X.
#' @param a_Y An ancestry label for dataset Y.
#' @keywords internal
validate_ancestry <- function(
  a_X,
  a_Y
) {
  stopifnot(length(a_X) == 1, length(a_Y) == 1)
}

#' Mean Difference Between Group Levels
#'
#' Computes column mean difference between two group levels.
#' @param data A matrix or data frame of numeric features.
#' @param group A factor indicating group membership for each row.
#' @param g1 The reference level.
#' @param g2 The comparison level.
#' @return A numeric vector of mean differences.
#' @keywords internal
mean_diff_by_group <- function(
  data,
  group,
  g1,
  g2
) {
  colMeans(data[group == g2, , drop = FALSE]) -
    colMeans(data[group == g1, , drop = FALSE])
}

#' Random Permutation of Indices
#'
#' Returns a random split of indices into two groups.
#' @param n_total Total number of samples.
#' @param n1 Size of the first group.
#' @return A list with two index vectors: `p1` and `p2`.
#' @keywords internal
permute_indices <- function(
  n_total,
  n1
) {
  idx <- sample.int(n_total)
  list(
    p1 = idx[1:n1],
    p2 = idx[(n1 + 1):n_total]
  )
}

#' Compute Test Statistic for Permutation
#'
#' Computes the interaction effect: difference of group-level differences.
#' Returns NA vector if any subgroup is empty.
#' @param XY A matrix of features.
#' @param g_XY A factor of group labels corresponding to rows of `XY`.
#' @param p1 Indices for the first permutation group.
#' @param p2 Indices for the second permutation group.
#' @param g1 First group label.
#' @param g2 Second group label.
#' @return A numeric vector of test statistics or NA values.
#' @keywords internal
compute_T_stat <- function(
  XY,
  g_XY,
  p1,
  p2,
  g1,
  g2
) {
  g1_p1 <- XY[p1[g_XY[p1] == g1], , drop = FALSE]
  g2_p1 <- XY[p1[g_XY[p1] == g2], , drop = FALSE]
  g1_p2 <- XY[p2[g_XY[p2] == g1], , drop = FALSE]
  g2_p2 <- XY[p2[g_XY[p2] == g2], , drop = FALSE]

  if (
    min(
      nrow(g1_p1),
      nrow(g2_p1),
      nrow(g1_p2),
      nrow(g2_p2)
    ) == 0
  ) {
    return(rep(NA_real_, ncol(XY)))
  }

  d1 <- colMeans(g1_p1) - colMeans(g2_p1)
  d2 <- colMeans(g1_p2) - colMeans(g2_p2)

  d1 - d2
}

#' Empirical P-Values
#'
#' Computes two-sided empirical p-values.
#' @param T_boot A matrix of bootstrapped test statistics.
#' @param T_obs A numeric vector of observed test statistics.
#' @return A numeric vector of p-values.
#' @keywords internal
compute_empirical_pvalues <- function(
  T_boot,
  T_obs
) {
  colMeans(
    abs(T_boot) >= abs(
      matrix(
        T_obs,
        nrow = nrow(T_boot),
        ncol = length(T_obs),
        byrow = TRUE
      )
    )
  )
}

#' Internal: Bootstrap Interaction Statistic (Adaptive or Fixed)
#'
#' Computes bootstrapped standard errors and 95% confidence intervals
#' for observed interaction statistics using nonparametric sampling within ancestries.
#'
#' @keywords internal
bootstrap_T_obs <- function(
  X,
  Y,
  MX,
  MY,
  g_col,
  a_col,
  g1,
  g2,
  B = NULL,
  seed = NULL,
  min_iter = 500,
  max_iter = 10000,
  tol = 1e-3,
  batch_size = 100,
  check_convergence = TRUE
) {
  if (!is.null(seed)) set.seed(seed)

  n_feat <- ncol(X)
  SE_prev <- rep(NA_real_, n_feat)
  max_B <- if (is.null(B)) max_iter else B
  T_boot <- matrix(NA_real_, nrow = max_B, ncol = n_feat)
  colnames(T_boot) <- colnames(X)

  converged <- if (is.null(B) || check_convergence) FALSE else NA
  B_used <- 0

  while ((is.null(B) && !isTRUE(converged) && B_used < max_iter) ||
         (!is.null(B) && B_used < B)) {

    for (i in seq_len(batch_size)) {
      if (!is.null(B) && B_used >= B) break

      B_used <- B_used + 1

      idx_X <- sample(nrow(X), replace = TRUE)
      idx_Y <- sample(nrow(Y), replace = TRUE)

      Xb <- X[idx_X, , drop = FALSE]
      Yb <- Y[idx_Y, , drop = FALSE]
      MXb <- MX[idx_X, , drop = FALSE]
      MYb <- MY[idx_Y, , drop = FALSE]

      g_Xb <- MXb[[g_col]]
      g_Yb <- MYb[[g_col]]

      if (any(table(g_Xb)[c(g1, g2)] == 0) || any(table(g_Yb)[c(g1, g2)] == 0)) {
        T_boot[B_used, ] <- NA_real_
        next
      }

      dX <- mean_diff_by_group(Xb, g_Xb, g1, g2)
      dY <- mean_diff_by_group(Yb, g_Yb, g1, g2)
      T_boot[B_used, ] <- dY - dX
    }

    T_valid <- T_boot[1:B_used, , drop = FALSE]
    T_valid <- T_valid[complete.cases(T_valid), , drop = FALSE]

    if ((is.null(B) || check_convergence) && nrow(T_valid) >= min_iter) {
      SE_new <- apply(T_valid, 2, sd)
      if (!anyNA(SE_prev)) {
        delta <- max(abs(SE_new - SE_prev))
        converged <- delta < tol
      } else {
        converged <- FALSE  # start checking from now
      }
      SE_prev <- SE_new
    }
  }

  # Final cleaned-up bootstrap results
  T_valid <- T_boot[1:B_used, , drop = FALSE]
  T_valid <- T_valid[complete.cases(T_valid), , drop = FALSE]

  SE_final <- apply(T_valid, 2, sd)
  CI_lower <- apply(T_valid, 2, quantile, probs = 0.025)
  CI_upper <- apply(T_valid, 2, quantile, probs = 0.975)

  summary_stats <- data.frame(
    feature = colnames(X),
    SE = SE_final,
    CI_lower = CI_lower,
    CI_upper = CI_upper,
    row.names = NULL
  )

  return(list(
    summary_stats = summary_stats,
    T_boot = T_valid,
    B_used = nrow(T_valid),
    converged = converged
  ))
}

