#' Simulate Negative Binomial RNA-seq Counts from Real Data
#'
#' Simulates RNA-seq count data using a Negative Binomial model, with parameters
#' estimated from real count data. Useful for creating realistic baseline data for
#' method development or DEG simulation.
#'
#' @param X A gene-by-sample matrix or data frame of real (integer) RNA-seq counts.
#' @param MX A data frame with rownames matching X and columns for ancestry and group.
#' @param g_col The name of the column in MX indicating experimental group (e.g., "group").
#' @param a_col The name of the column in MX indicating ancestry (e.g., "ancestry").
#' @param seed Optional random seed for reproducibility.
#'
#' @return A list containing:
#' \describe{
#'   \item{counts}{A simulated count matrix of the same shape as X (samples × genes).}
#'   \item{sample_info}{Metadata with ancestry, group, and interaction group.}
#'   \item{params}{A data frame of per-gene estimated means and dispersions.}
#' }
#'
#' @export
simulate_NB_counts <- function(
  X,
  MX,
  g_col,
  a_col,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  # Create interaction-based group column
  MX$group <- interaction(MX[[a_col]], MX[[g_col]], sep = "_")

  # Filter lowly expressed genes (e.g., <5 counts in nearly all samples)
  keep_genes <- colSums(X >= 5) >= 2
  X <- X[, keep_genes, drop = FALSE]

  mu_base <- rowMeans(X)
  dispersion <- apply(X, 1, function(x) {
    m <- mean(x)
    v <- var(x)
    if (v <= m || is.na(v)) return(1)
    (m^2) / (v - m)
  })

  # Simulate counts from Negative Binomial
  counts_sim <- matrix(0, nrow = nrow(X), ncol = ncol(X))
  rownames(counts_sim) <- rownames(X)
  colnames(counts_sim) <- colnames(X)

  for (i in seq_len(ncol(X))) {
    counts_sim[, i] <- rnbinom(n = nrow(X), mu = mu_base, size = dispersion)
  }

  list(
    counts = counts_sim,  # genes × samples
    sample_info = MX,
    params = data.frame(
      feature = rownames(X),
      mu = mu_base,
      dispersion = dispersion,
      row.names = NULL
    )
  )
}
