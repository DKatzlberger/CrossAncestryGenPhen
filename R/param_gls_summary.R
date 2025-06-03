#' Parametric GLS Feature Summary Across Iterations
#'
#' Aggregates feature-level test statistics across iterations using a Generalized Least Squares (GLS) framework,
#' assuming normality of the GLS estimator. Computes parametric p-values and summary statistics for each feature.
#'
#' @param x A data frame in long format with columns: `feature`, `iteration`, `T_obs`, `SE`, and `p_adj`.
#' @param alpha Numeric. Significance threshold used to compute the proportion of significant iterations.
#'
#' @return A data frame with feature, GLS estimate, SE, z-score, p-value, FDR-adjusted p-value, and prop_signif.
#' @export
param_gls_summary <- function(
  x, 
  alpha = 0.05
) {
  feature_names <- unique(x$feature)
  iter_names <- unique(x$iteration)

  # Build T matrix to estimate correlation across iterations
  T_mat <- build_T_matrix(x, feature_names, iter_names)
  R_full <- cor(T_mat, use = "pairwise.complete.obs", method = "pearson")

  results <- vector("list", length(feature_names))
  for (f in seq_along(feature_names)) {

    feature <- feature_names[f]
    df_f <- x[x$feature == feature, ]

    iter_ids <- as.character(df_f$iteration)
    T_vec <- df_f$T_obs
    SE_vec <- df_f$SE

    mean_T_obs <- mean(T_vec, na.rm = TRUE)
    prop_signif <- mean(df_f$p_adj < alpha, na.rm = TRUE)

    estimate <- NA_real_
    SE <- NA_real_
    p_val <- NA_real_

    if (length(T_vec) >= 2 && !anyNA(T_vec) && !anyNA(SE_vec)) {
      R_sub <- R_full[iter_ids, iter_ids, drop = FALSE]
      obs_result <- compute_gls_estimate(T_vec, SE_vec, R_sub)
      mu_hat <- obs_result$mu_hat
      se_hat <- obs_result$se_hat

      if (!is.na(mu_hat) && !is.na(se_hat)) {
        estimate <- mu_hat
        SE <- se_hat
        p_val <- 2 * pnorm(-abs(mu_hat / se_hat))
      } else {
        message("GLS failed for feature: ", feature)
      }
    } 

    results[[f]] <- data.frame(
      feature = feature,
      mean_T_obs = mean_T_obs,
      estimate = estimate,
      SE = SE,
      p_value = p_val,
      p_adj = NA_real_,
      prop_signif = prop_signif
    )
  }

  summary_stats <- do.call(rbind, results)
  summary_stats$p_adj <- p.adjust(summary_stats$p_value, method = "BH")

  return(summary_stats)
}
