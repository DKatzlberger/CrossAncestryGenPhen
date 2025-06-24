#' Compute Correlation Matrix from Long-format Data
#'
#' Transforms long-format data to wide format and computes a correlation matrix across iterations.
#'
#' @param x A data frame or data.table with columns `feature`, a user-defined iteration column, and a value column.
#' @param value_col The name of the column to pivot on (e.g., "p_value" or "T_obs").
#' @param iter_col The name of the column to treat as "iteration" (e.g., "sample_id", "perm_id").
#' @param method Correlation method: one of "pearson", "spearman", or "kendall".
#'
#' @return A correlation matrix with NA on the diagonal.
#' @import data.table
#' @export
compute_correlation_matrix <- function(
    x, 
    value_col, 
    iter_col,
    method = "pearson"
) {
  stopifnot(
    is.data.frame(x),
    all(c("feature", iter_col, value_col) %in% names(x)),
    method %in% c("pearson", "spearman", "kendall")
  )

  # Correlation matrix computation
  dt <- as.data.table(x)
  formula_str <- sprintf("feature ~ %s", iter_col)
  wide <- dcast(dt, formula = as.formula(formula_str), value.var = value_col)

  mat <- as.matrix(wide[, -1, with = FALSE])
  cor_mat <- cor(mat, use = "pairwise.complete.obs", method = method)

  colnames(cor_mat) <- rownames(cor_mat) <- colnames(mat)

  return(cor_mat)
}
