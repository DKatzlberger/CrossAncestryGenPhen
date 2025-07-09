#' Simple Permutation-Based Interaction Test (No Bootstrapping, No Convergence)
#'
#' Tests if condition effects (g1 vs g2) are different between ancestries using permutation.
#'
#' @param X Expression matrix for ancestry A (samples x genes)
#' @param Y Expression matrix for ancestry B (samples x genes)
#' @param MX Metadata for X (must include group and ancestry columns)
#' @param MY Metadata for Y
#' @param g_col Column in metadata for condition/group (factor with 2 levels)
#' @param B Number of permutations
#' @param seed Optional seed for reproducibility
#'
#' @return A list with summary stats, permutation matrix, and number of valid permutations
#' 
#' @useDynLib CrossAncestryGenPhen
#' @importFrom Rcpp evalCpp
#' 
#' @export
perm_interaction_effect_cpp <- function(
  X,
  Y,
  MX,
  MY,
  g_col,
  B = 1000,
  seed = NULL
) {

  if (!is.null(seed)) set.seed(seed)

  # Define groups that are compared
  g_X <- MX[[g_col]]
  g_Y <- MY[[g_col]]

  # Check factor
  stopifnot(is.factor(g_X))
  stopifnot(is.factor(g_Y))
  stopifnot(length(levels(g_X)) == 2)
  stopifnot(length(levels(g_Y)) == 2)

  # Make sure the reference is correct
  g1 <- levels(g_X)[1]
  g2 <- levels(g_X)[2]

  # Difference in groups
  delta_X <- colMeans(X[g_X == g2, , drop = FALSE]) - colMeans(X[g_X == g1, , drop = FALSE])
  delta_Y <- colMeans(Y[g_Y == g2, , drop = FALSE]) - colMeans(Y[g_Y == g1, , drop = FALSE])

  # Interactions (difference of differences)
  T_obs <- delta_Y - delta_X

  XY <- rbind(X, Y)
  MXY <- rbind(MX, MY)
  g_XY <- factor(MXY[[g_col]], levels = levels(g_X))

  ave_expr <- colMeans(XY)
  n1 <- nrow(X)
  n2 <- nrow(Y)

  # Labels must be integer for C++
  g1_label <- which(levels(g_XY) == g1)
  g2_label <- which(levels(g_XY) == g2)

  # Rccp compiled version
  T_perm <- permInteractionCpp(
    XY = as.matrix(XY),
    g_XY = as.integer(g_XY),
    n1 = n1,
    n2 = n2,
    B = B,
    g1 = g1_label,
    g2 = g2_label
  )

  valid_rows <- complete.cases(T_perm)
  T_perm <- T_perm[valid_rows, , drop = FALSE]
  B_actual <- nrow(T_perm)

  p_vals_emp <- numeric(ncol(XY))
  p_vals_param <- numeric(ncol(XY))
  z_scores <- numeric(ncol(XY))

  for (j in seq_len(ncol(XY))) {
    T_j <- T_perm[, j]
    obs <- T_obs[j]

    p_vals_emp[j] <- (sum(abs(T_j) >= abs(obs)) + 1) / (length(T_j) + 1)

    mu <- mean(T_j)
    sigma <- sd(T_j)
    z <- (obs - mu) / sigma
    z_scores[j] <- z
    p_vals_param[j] <- 2 * (1 - pnorm(abs(z)))
  }

  p_adj_emp <- p.adjust(p_vals_emp, method = "BH")
  p_adj_param <- p.adjust(p_vals_param, method = "BH")

  list(
    summary_stats = data.frame(
      feature = colnames(XY),
      T_obs = T_obs,
      z_score = z_scores,
      p_value = p_vals_param,
      p_adj = p_adj_param,
      ave_expr = ave_expr,
      row.names = NULL
    ),
    p_emp_value = p_vals_emp,
    T_null = T_perm,
    B_used = B_actual
  )
}
