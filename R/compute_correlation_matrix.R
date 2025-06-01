#' Compute Correlation Matrix from Long-format Data
#'
#' Transforms long-format data to wide format and computes a correlation matrix across iterations.
#'
#' @param x A data frame or data.table with columns `feature`, `iteration`, and a value column.
#' @param value_col The name of the column to pivot on (e.g., "p_value" or "T_obs").
#' @param method Correlation method: one of "pearson", "spearman", or "kendall".
#'
#' @return A correlation matrix with NA on the diagonal.
#' @import data.table
#' @export
compute_correlation_matrix <- function(
    x, 
    value_col, 
    method = "pearson"
) {

  stopifnot(
    is.data.frame(x),
    all(c("feature", "iteration", value_col) %in% names(x)),
    method %in% c("pearson", "spearman", "kendall")
  )

  # Correlation matrix computation
  dt <- as.data.table(x)
  wide <- dcast(dt, feature ~ iteration, value.var = value_col)
  mat <- as.matrix(wide[, -1, with = FALSE])  
  cor_mat <- cor(mat, use = "pairwise.complete.obs", method = method)

  n_iter <- ncol(cor_mat)
  names_vec <- as.character(seq_len(n_iter))
  colnames(cor_mat) <- rownames(cor_mat) <- as.character(names_vec)
  
  return(cor_mat)
}
