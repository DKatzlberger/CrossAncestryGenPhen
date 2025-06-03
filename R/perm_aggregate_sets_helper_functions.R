#' Empirical P-Values
#'
#' Computes two-sided empirical p-values.
#' @param T_obs A numeric vector of observed test statistics.
#' @param T_null A numeric vector of null distribution test statistics.
#' @return A numeric vector of p-values.
#' @keywords internal
compute_empirical_p <- function(
  T_obs,
  T_null,
  alternative = c("two.sided", "greater", "less")
) {
  alternative <- match.arg(alternative)
  
  if (alternative == "two.sided") {
    p <- (sum(abs(T_null) >= abs(T_obs)) + 1) / (length(T_null) + 1)
  } else if (alternative == "greater") {
    p <- (sum(T_null >= T_obs) + 1) / (length(T_null) + 1)
  } else if (alternative == "less") {
    p <- (sum(T_null <= T_obs) + 1) / (length(T_null) + 1)
  }
  
  return(p)
}