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
  colMeans(data[group == g1, , drop = FALSE]) -
    colMeans(data[group == g2, , drop = FALSE])
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

#' Confidence Interval Bounds
#'
#' Computes 2.5% and 97.5% quantiles for each feature.
#' @param T_boot A matrix of bootstrapped test statistics.
#' @return A matrix of lower and upper confidence bounds.
#' @keywords internal
compute_CI_bounds <- function(
  T_boot
) {
  apply(
    T_boot,
    2,
    quantile,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
}
