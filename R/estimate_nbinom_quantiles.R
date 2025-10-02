#' Negative Binomial quantiles per gene
#'
#' Computes quantiles of the Negative Binomial distribution for each gene,
#' given mean (\code{mu}) and dispersion (\code{phi}) values.
#'
#' @param mu Numeric vector of means (>0).
#' @param phi Numeric vector of dispersions (>0), same length as \code{mu}.
#' @param probs Probabilities in [0, 1] at which to compute quantiles.
#'
#' @return Matrix of quantiles with rows = \code{probs} and cols = genes.
#' @export
estimate_nbinom_quantiles <- function(
  mu, 
  phi, 
  probs
){

  ## --- Input check ---
  if (length(mu) != length(phi)){
    stop("[estimate_nbinom_quantiles] Same numbers of mu and phi are required.")
  }

  ## --- Initialize matrix ---
  n_genes  <- length(mu)
  n_quant  <- length(probs)
  q_matrix <- matrix(NA, nrow = n_quant, ncol = n_genes)

  for (j in seq_len(n_genes)) {
    mu_j  <- mu[j]
    phi_j <- phi[j]
    
    if (is.na(mu_j) || is.na(phi_j) || mu_j <= 0 || phi_j <= 0) next
    
    size_j <- 1 / phi_j
    q_matrix[, j] <- qnbinom(p = probs, mu = mu_j, size = size_j)
  }
  rownames(q_matrix) <- paste0(round(probs * 100, 5), "%")
  
  ## --- Return ---
  return(q_matrix)
}