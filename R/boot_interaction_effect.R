#' Simple Bootstrap-Based Interaction Estimation with CI Level Control
#'
#' Estimates variability in the condition-by-ancestry interaction effect using non-stratified bootstrap resampling.
#'
#' @param X Expression matrix for ancestry A (samples x genes)
#' @param Y Expression matrix for ancestry B (samples x genes)
#' @param MX Metadata for X (must include group column)
#' @param MY Metadata for Y
#' @param g_col Column in metadata for condition/group (factor with 2 levels)
#' @param B Number of bootstrap samples
#' @param seed Optional seed for reproducibility
#' @param alpha Significance level (e.g., 0.05 for 95\% CI)
#'
#' @return A list with summary stats and bootstrap replicates
#' @export
boot_interaction_effect <- function(
  X,
  Y,
  MX,
  MY,
  g_col,
  B = 1000,
  seed = NULL,
  alpha = 0.05
) {

  if (!is.null(seed)) set.seed(seed)

  # Define groups that are compared
  g_X <- MX[[g_col]]
  g_Y <- MY[[g_col]]

  # Make sure the reference is correct
  g1 <- levels(g_X)[1]
  g2 <- levels(g_X)[2]

  # Difference in groups
  delta_X <- colMeans(X[g_X == g2, , drop = FALSE]) - colMeans(X[g_X == g1, , drop = FALSE])
  delta_Y <- colMeans(Y[g_Y == g2, , drop = FALSE]) - colMeans(Y[g_Y == g1, , drop = FALSE])

  # Interactions (difference of differences)
  T_obs <- delta_Y - delta_X

  # Initialize boot matrix
  T_boot <- matrix(NA_real_, nrow = B, ncol = ncol(X))
  colnames(T_boot) <- colnames(X)
  B_used <- 0

  for (b in seq_len(B)) {

    # Resample X and Y independently (non-stratified)
    idx_X <- sample(nrow(X), replace = TRUE)
    idx_Y <- sample(nrow(Y), replace = TRUE)

    Xb <- X[idx_X, , drop = FALSE]
    Yb <- Y[idx_Y, , drop = FALSE]

    g_Xb <- g_X[idx_X]
    g_Yb <- g_Y[idx_Y]

    # Skip if any group missing in either ancestry
    if (!all(c(g1, g2) %in% g_Xb) || !all(c(g1, g2) %in% g_Yb)) next

    # Group means and interaction contrast
    dX <- colMeans(Xb[g_Xb == g2, , drop = FALSE]) - colMeans(Xb[g_Xb == g1, , drop = FALSE])
    dY <- colMeans(Yb[g_Yb == g2, , drop = FALSE]) - colMeans(Yb[g_Yb == g1, , drop = FALSE])
    T_boot[B_used + 1, ] <- dY - dX
    B_used <- B_used + 1
  }

  # Keep only valid bootstrap replicates
  T_boot <- T_boot[seq_len(B_used), , drop = FALSE]
  if (B_used == 0) stop("All bootstrap replicates failed due to missing groups.")

  # Compute important statistics
  sd <- apply(T_boot, 2, var)
  se <- apply(T_boot, 2, sd)
  ci <- t(apply(T_boot, 2, quantile, probs = c(alpha / 2, 1 - alpha / 2)))

  summary_stats <- data.frame(
    feature = colnames(X),
    T_obs = T_obs,
    sd = sd,
    se = se,
    ci_lower = ci[, 1],
    ci_upper = ci[, 2],
    row.names = NULL
  )

  return(list(
    summary_stats = summary_stats,
    T_boot = T_boot,
    B_used_boot = B_used
  ))
}
