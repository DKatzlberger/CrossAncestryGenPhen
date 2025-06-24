#' Aggregate gene statistics using Rank Product
#'
#' This function aggregates results from multiple iterations or studies
#' by calculating the Rank Product (RP) of raw p-values, the mean of
#' test statistics (`T_obs`), and the proportion of significance based
#' on FDR-adjusted p-values. Optionally applies jitter to p-values to
#' break ties in low-resolution data. If jitter is not applied, the
#' function reports the fraction of tied p-values in each iteration
#' and the mean fraction across iterations to help users assess whether
#' jittering might be appropriate.
#'
#' @param x A data.frame with at least the following columns:
#'   \itemize{
#'     \item \code{feature}: Gene or feature identifier.
#'     \item \code{T_obs}: Observed test statistic per iteration.
#'     \item \code{p_value}: Raw p-value.
#'     \item \code{p_adj}: FDR-adjusted p-value.
#'     \item \code{iteration}: Replicate or study identifier.
#'   }
#' @param fdr_threshold Numeric. FDR-adjusted p-value threshold used
#'   to determine significance (default is 0.05).
#' @param jitter_p Logical. Whether to apply small uniform random jitter to p-values
#'   to break ties and improve rank resolution (default is FALSE).
#' @param jitter_amount Numeric. Maximum amount of uniform noise to add to
#'   p-values during jittering (default is 1e-6).
#'
#' @return A data.frame with the following columns:
#'   \itemize{
#'     \item \code{feature}: Gene or feature identifier.
#'     \item \code{mean_T_obs}: Mean of observed test statistics across iterations.
#'     \item \code{RP}: Rank Product score (lower = stronger consistent signal).
#'     \item \code{prop_sig}: Proportion of iterations with FDR-adjusted p < \code{fdr_threshold}.
#'   }
#'
#' @examples
#' \dontrun{
#' result <- RP_aggregation(
#'   combined_results,
#'   jitter_p = TRUE,
#'   jitter_amount = 1e-6
#' )
#' head(result)
#' }
#' @export
RP_aggregation <- function(
    x, 
    fdr_threshold = 0.05,
    jitter_p = FALSE,
    jitter_amount = 1e-6
) {
  genes <- unique(x$feature)
  iterations <- sort(unique(x$iteration))
  
  # Construct gene x iteration matrix of raw p-values
  pval_matrix <- matrix(
    NA, nrow = length(genes), 
    ncol = length(iterations),
    dimnames = list(genes, as.character(iterations))
  )
  
  for (i in seq_len(nrow(x))) {
    gene <- x$feature[i]
    iter <- as.character(x$iteration[i])
    pval_matrix[gene, iter] <- x$p_value[i]
  }

  # If no jitter, report detailed tie statistics
  if (!jitter_p) {
    message("Tied p-value report per iteration:")
    tie_fractions <- numeric(ncol(pval_matrix))
    
    for (j in seq_len(ncol(pval_matrix))) {
      p <- na.omit(pval_matrix[, j])
      tab <- table(p)
      ties <- tab[tab > 1]
      n_tied <- sum(ties)
      n_total <- length(p)
      frac_tied <- if (n_total > 0) n_tied / n_total else 0
      tie_fractions[j] <- frac_tied
      message(sprintf("  Iteration %s: %d/%d p-values are tied (%.2f%%)",
                      colnames(pval_matrix)[j],
                      n_tied, n_total, 100 * frac_tied))
    }
    
    mean_frac <- mean(tie_fractions)
    message(sprintf("Mean fraction of tied p-values across iterations: %.2f%%", 100 * mean_frac))
  }

  # Apply jitter if requested
  if (jitter_p) {
    message("Applying jitter to p-values.")
    jitter_matrix <- matrix(runif(length(pval_matrix), 0, jitter_amount), nrow = nrow(pval_matrix))
    pval_matrix <- pval_matrix + jitter_matrix
    pval_matrix[pval_matrix > 1] <- 1  # Keep values in [0, 1]
  }

  # Rank p-values per iteration
  ranked_pvals <- apply(pval_matrix, 2, function(col) rank(col, ties.method = "average", na.last = "keep"))
  
  # Compute Rank Product
  rank_product <- apply(ranked_pvals, 1, function(ranks) {
    ranks <- ranks[!is.na(ranks)]
    prod(ranks)^(1 / length(ranks))
  })
  
  # Mean T_obs
  mean_T_obs <- tapply(x$T_obs, x$feature, mean, na.rm = TRUE)
  
  # Proportion significant
  sig_flags <- x$p_adj < fdr_threshold
  prop_signif <- tapply(sig_flags, x$feature, function(x) mean(x, na.rm = TRUE))
  
  # Assemble result
  result <- data.frame(
    feature = names(rank_product),
    mean_T_obs = mean_T_obs[names(rank_product)],
    RP = rank_product,
    prop_sig = prop_signif[names(rank_product)],
    row.names = NULL
  )
  
  result <- result[order(result$RP), ]
  return(result)
}
