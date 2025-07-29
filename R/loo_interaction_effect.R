#' Leave-One-Out Diagnostic for interaction effects
#'
#' Computes leave-one-out diagnostics for each feature by dropping each 
#' sample one at a time from each ancestry dataset. Shows how much each 
#' sample influences the interaction effect estimate.
#'
#' @param X A numeric matrix or data.frame of expression values for ancestry A. 
#' @param Y A numeric matrix or data.frame of expression values for ancestry B.
#' @param MX A data.frame of metadata for X. Must include a column with condition/group labels.
#' @param MY A data.frame of metadata for Y. Must include a column with condition/group labels.
#' @param g_col String. Name of the column in `MX` and `MY` that defines the condition/group (must be a factor with exactly 2 levels).
#' @param a_col String. Name of the column in `MX` and `MY` that defines ancestry.
#' @param feature Character or numeric vector. Names or indices of features (columns) to test.
#'
#' @return A tidy data.frame with columns:
#'   \describe{
#'     \item{feature}{Feature name or index.}
#'     \item{sample}{Sample name that was left out.}
#'     \item{T_loo}{Interaction effect with that sample left out.}
#'     \item{T_obs}{Original full-sample interaction effect.}
#'     \item{ancestry}{Ancestry of the sample ("eur" or "afr").}
#'   }
#'
#' @export
loo_interaction_effect <- function(
  X, 
  Y, 
  MX, 
  MY, 
  g_col,
  a_col,
  feature
) {
  results_list <- list()
  
  g_X <- MX[[g_col]]
  g_Y <- MY[[g_col]]

  a_X <- unique(MX[[a_col]])
  a_Y <- unique(MY[[a_col]])

  if (length(a_X) != 1 | length(a_Y) != 1) {
    stop("Each ancestry dataset must contain only one unique ancestry.")
  }

  g1 <- levels(g_X)[1]
  g2 <- levels(g_X)[2]

  # Compute full-sample deltas
  delta_X <- apply(X[g_X == g2, , drop = FALSE], 2, mean) -
             apply(X[g_X == g1, , drop = FALSE], 2, mean)
  delta_Y <- apply(Y[g_Y == g2, , drop = FALSE], 2, mean) -
             apply(Y[g_Y == g1, , drop = FALSE], 2, mean)
  T_obs <- delta_Y - delta_X

  for (f in feature) {
    T_loo <- c()
    samples <- c()
    ancestry <- c()

    # Resolve feature index
    if (is.character(f)) {
      f_idx <- which(colnames(X) == f)
      if (length(f_idx) == 0) {
        warning(paste("Feature not found, skipping:", f))
        next
      }
    } else {
      f_idx <- f
      if (f_idx < 1 || f_idx > ncol(X)) {
        warning(paste("Feature index out of bounds, skipping:", f))
        next
      }
    }

    # Loop over X samples
    for (i in seq_len(nrow(X))) {
      X_loo <- X[-i, , drop = FALSE]
      MX_loo <- MX[-i, , drop = FALSE]
      g_X_loo <- MX_loo[[g_col]]
      a_X_out <- as.character(MX[[a_col]][i])

      delta_X_loo <- mean(X_loo[g_X_loo == g2, f_idx, drop = TRUE]) -
                     mean(X_loo[g_X_loo == g1, f_idx, drop = TRUE])

      T_loo <- c(T_loo, delta_Y[f_idx] - delta_X_loo)
      samples <- c(samples, rownames(X)[i])
      ancestry <- c(ancestry, a_X_out)
    }

    # Loop over Y samples
    for (j in seq_len(nrow(Y))) {
      Y_loo <- Y[-j, , drop = FALSE]
      MY_loo <- MY[-j, , drop = FALSE]
      g_Y_loo <- MY_loo[[g_col]]
      a_Y_out <- as.character(MY[[a_col]][j])

      delta_Y_loo <- mean(Y_loo[g_Y_loo == g2, f_idx, drop = TRUE]) -
                     mean(Y_loo[g_Y_loo == g1, f_idx, drop = TRUE])

      T_loo <- c(T_loo, delta_Y_loo - delta_X[f_idx])
      samples <- c(samples, rownames(Y)[j])
      ancestry <- c(ancestry, a_Y_out)
    }

    # Store results for this feature
    results_list[[as.character(f)]] <- data.frame(
      feature = rep(f, length(samples)),
      sample = samples,
      T_loo = T_loo,
      T_obs = rep(T_obs[f_idx], length(samples)),
      ancestry = ancestry,
      stringsAsFactors = FALSE
    )
  }

  # Combine all features
  results <- do.call(rbind, results_list)
  rownames(results) <- NULL  

  return(results)
}
