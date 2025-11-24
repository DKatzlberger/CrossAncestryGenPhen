#' Beta to M-value transformation
#'
#' Converts DNA methylation beta-values (0–1) to M-values using a logit transform:
#' \deqn{M = \log_2(\beta / (1 - \beta))}
#' M-values are unbounded and better suited for linear models (e.g., limma).
#'
#' @param X Numeric matrix/data.frame of beta-values (0–1).
#' @param offset Small value to avoid log(0); default = 1e-6.
#'
#' @return Matrix of M-values with same dimensions as `X`.
#' @export
beta_to_mval <- function(
  X,
  offset = 1e-6
){

  ## --- Input data structure check ---
  assert_input(
    X = X, 
    .fun = "beta_to_mval"
  )

  # NA values present?
  if (anyNA(X)) {
    warning("[beta_to_mval] NA values detected. These will be preserved.")
  }

  # Negative values?
  if (any(X < 0, na.rm = TRUE)) {
    warning("[beta_to_mval] Negative beta-values detected. These are invalid and will be clipped to 'offset'.")
  }


  ## --- Clip to avoid 0 or 1 ---
  X[X <= 0] <- offset
  X[X >= 1] <- 1 - offset

  ## --- Transform ---
  M <- log2(X / (1 - X))

  
  ## --- Return ---
  return(M)
}