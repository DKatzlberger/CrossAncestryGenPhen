#' Build T-statistic Matrix (Features x Iterations)
#'
#' Constructs a matrix of T statistics where rows are features and columns are iterations.
#'
#' @keywords internal
build_T_matrix <- function(x, features, iterations) {
  mat <- matrix(NA_real_, nrow = length(features), ncol = length(iterations))
  rownames(mat) <- features
  colnames(mat) <- as.character(iterations)

  for (f in features) {
    sub <- x[x$feature == f, ]
    mat[f, as.character(sub$iteration)] <- sub$T_obs
  }

  return(mat)
}


#' Compute GLS Estimate
#'
#' Computes a weighted average estimate using Generalized Least Squares, accounting for
#' heteroskedasticity and correlation via a supplied covariance structure.
#'
#' @keywords internal
compute_gls_estimate <- function(
  T_vec, 
  SE_vec, 
  R_mat
) {
  # Check for mismatched lengths or missing values
  if (length(T_vec) != length(SE_vec)) {
    message("T_vec and SE_vec have different lengths.")
    return(list(mu_hat = NA, se_hat = NA))
  }
  if (anyNA(T_vec) || anyNA(SE_vec)) {
    message("Missing values detected in T_vec or SE_vec.")
    return(list(mu_hat = NA, se_hat = NA))
  }
  if (anyNA(R_mat)) {
    message("Missing values detected in correlation matrix R_mat.")
    return(list(mu_hat = NA, se_hat = NA))
  }

  # Build full covariance matrix: S = D * R * D
  D <- diag(SE_vec)
  S <- D %*% R_mat %*% D
  kappa_val <- kappa(S)  # Condition number of S

  # Warning for ill-conditioned covariance matrix
  if (kappa_val > 1e10) {
    message("Ill-conditioned covariance matrix (kappa = ", signif(kappa_val, 3), ").")
  }

  # Invert the covariance matrix S
  S_inv <- tryCatch(solve(S), error = function(e) {
    message("Failed to invert covariance matrix: ", e$message)
    return(NULL)
  })

  # Return NA if inversion failed
  if (is.null(S_inv)) {
    return(list(mu_hat = NA, se_hat = NA))
  }

  # GLS estimator for the common mean across correlated iterations
  n <- length(T_vec)
  one_vec <- rep(1, n)

  weight_sum <- as.numeric(t(one_vec) %*% S_inv %*% one_vec)
  weighted_T_sum <- as.numeric(t(one_vec) %*% S_inv %*% T_vec)

  mu_hat <- weighted_T_sum / weight_sum
  se_hat <- sqrt(1 / weight_sum)

  return(list(mu_hat = mu_hat, se_hat = se_hat))
}
