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
#' @export
perm_interaction_effect <- function(
  X,          
  Y,          
  MX,         
  MY,         
  g_col,      
  B = 1000,   
  seed = NULL        
){

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

  # Combine the data
  XY <- rbind(X, Y)
  MXY <- rbind(MX, MY)
  # g_XY <- MXY[[g_col]]
  g_XY <- factor(MXY[[g_col]], levels = levels(g_X))

  # Calculate mean expression
  ave_expr <- colMeans(XY) 

  # Sample sizes
  n1 <- nrow(X)
  n2 <- nrow(Y)

  # Initilaize perm matrix
  T_perm <- matrix(0, nrow = B, ncol = ncol(XY))
  colnames(T_perm) <- colnames(XY)

  for (b in 1:B){
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
      T_perm[b, ] <- NA
      next
    }

    # Compute group effects and interaction
    dX <- colMeans(g2_p1) - colMeans(g1_p1)
    dY <- colMeans(g2_p2) - colMeans(g1_p2)
    T_perm[b, ] <- dY - dX
  }

  # Remove any NAs caused by bad splits
  valid_rows <- complete.cases(T_perm)
  T_perm <- T_perm[valid_rows, , drop = FALSE]
  B_actual <- nrow(T_perm)

  # Statistics
  # Initialize vectors
  p_vals_emp <- numeric(ncol(XY))
  p_vals_param <- numeric(ncol(XY))
  z_scores <- numeric(ncol(XY))

  # Loop through features
  for (j in seq_len(ncol(XY))) {
    T_j <- T_perm[, j]
    obs <- T_obs[j]
    
    # Empirical p-value
    p_vals_emp[j] <- (sum(abs(T_j) >= abs(obs)) + 1) / (length(T_j) + 1)
    
    # Parametric (z-based) p-value
    mu <- mean(T_j)
    sigma <- sd(T_j)
    z <- (obs - mu) / sigma
    z_scores[j] <- z
    p_vals_param[j] <- 2 * (1 - pnorm(abs(z)))
  }

  # FDR corrections
  p_adj_emp <- p.adjust(p_vals_emp, method = "BH")
  p_adj_param <- p.adjust(p_vals_param, method = "BH")

  # Return structured results
  return(
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
  )
}
