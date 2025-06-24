#' Simple Permutation Test for Correlation Differences
#'
#' Tests whether the correlation of condition effects between ancestries X and Y
#' with a reference ancestry R differs. No bootstrapping or convergence logic.
#'
#' @param X Matrix of expression values (ancestry 1); samples x genes
#' @param Y Matrix of expression values (ancestry 2)
#' @param R Matrix of expression values (reference ancestry)
#' @param MX Metadata for X; must contain group column
#' @param MY Metadata for Y
#' @param MR Metadata for R
#' @param g_col Column in metadata specifying condition/group (factor with 2 levels)
#' @param method Correlation method: "pearson" (default) or "spearman"
#' @param B Number of permutations
#' @param seed Random seed (optional)
#'
#' @return A list with summary statistics and null distribution
#' @export
perm_correlation_difference <- function(
  X,
  Y,
  R,
  MX,
  MY,
  MR,
  g_col,
  method = c("pearson", "spearman"),
  B = 1000,
  seed = NULL
) {
  method <- match.arg(method)
  if (!is.null(seed)) set.seed(seed)

  # Define groups that are compared
  g_X <- MX[[g_col]]
  g_Y <- MY[[g_col]]
  g_R <- MR[[g_col]]

  # Make sure the reference is correct
  g1 <- levels(factor(g_X))[1]
  g2 <- levels(factor(g_X))[2]

  # Reference group difference (fixed)
  delta_R <- colMeans(R[g_R == g2, , drop = FALSE]) - colMeans(R[g_R == g1, , drop = FALSE])

  # Observed differences in X and Y
  delta_X <- colMeans(X[g_X == g2, , drop = FALSE]) - colMeans(X[g_X == g1, , drop = FALSE])
  delta_Y <- colMeans(Y[g_Y == g2, , drop = FALSE]) - colMeans(Y[g_Y == g1, , drop = FALSE])

  # Observed correlations with reference
  cor_XR <- cor(delta_X, delta_R, use = "complete.obs", method = method)
  cor_YR <- cor(delta_Y, delta_R, use = "complete.obs", method = method)
  T_obs <- cor_YR - cor_XR

  # Permutation setup
  XY <- rbind(X, Y)
  MXY <- rbind(MX, MY)
  g_XY <- MXY[[g_col]]

  # Sample sizes
  n1 <- nrow(X)
  n2 <- nrow(Y)

  # Initilaize perm vector
  T_perm <- rep(NA_real_, B)

  for (b in 1:B) {

    # Randomly split into two pseudo ancestries
    perm_idx <- sample(1:nrow(XY))
    p1 <- perm_idx[1:n1]
    p2 <- perm_idx[(n1 + 1):(n1 + n2)]

    # Group labels for each pseudo ancestry
    g1_p1 <- XY[p1[g_XY[p1] == g1], , drop = FALSE]
    g2_p1 <- XY[p1[g_XY[p1] == g2], , drop = FALSE]

    g1_p2 <- XY[p2[g_XY[p2] == g1], , drop = FALSE]
    g2_p2 <- XY[p2[g_XY[p2] == g2], , drop = FALSE]

   # Skip iteration if any group is empty
    if (any(sapply(list(g1_p1, g2_p1, g1_p2, g2_p2), nrow) == 0)) {
      T_perm[b] <- NA
      next
    }

    # Compute group effects 
    dX <- colMeans(g2_p1) - colMeans(g1_p1)
    dY <- colMeans(g2_p2) - colMeans(g1_p2)

    cXR <- cor(dX, delta_R, use = "complete.obs", method = method)
    cYR <- cor(dY, delta_R, use = "complete.obs", method = method)

    T_perm[b] <- cYR - cXR
  }

  # Remove invalid permutations
  T_perm <- T_perm[!is.na(T_perm)]
  B_used <- length(T_perm)

  if (B_used == 0) stop("No valid permutations completed.")
  p_vals_emp <- (sum(abs(T_perm) >= abs(T_obs)) + 1) / (B_used + 1)

  # Parametric p-values assuming normal distribution
  mu <- mean(T_perm)
  sigma <- sd(T_perm)
  z <- (T_obs - mu) / sigma
  p_vals_param <- 2 * (1 - pnorm(abs(z)))


  summary_stats <- data.frame(
    feature = "Global",
    T_obs = T_obs,
    XR = cor_XR,
    YR = cor_YR,
    p_param_value = p_vals_param,
    p_emp_value = p_vals_emp,
    row.names = NULL
  )

  return(list(
    summary_stats = summary_stats,
    T_null = matrix(T_perm, ncol = 1, dimnames = list(NULL, "T_null")),
    B_used = B_used
  ))
}
