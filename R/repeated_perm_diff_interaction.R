#' Repeated Permutation Interaction Test on Stratified Subsets
#'
#' Performs multiple stratified subsamples of the overrepresented group,
#' runs `perm_interaction()` on each, and aggregates the results.
#'
#' @param X Expression matrix (samples × genes) for the overrepresented ancestry.
#' @param Y Expression matrix for the underrepresented ancestry.
#' @param MX Metadata for X.
#' @param MY Metadata for Y.
#' @param stratify_cols Character vector of metadata column names to stratify by.
#' @param g_col Name of the grouping variable (e.g., "condition").
#' @param a_col Name of the ancestry variable.
#' @param features Optional vector of feature names to include (default: all).
#' @param n_iter Number of stratified subsample iterations.
#' @param B Number of permutations in each `perm_interaction()` call.
#' @param seed Base random seed for reproducibility.
#'
#' @return A list with:
#' \describe{
#'   \item{aggregated}{Data frame with per-feature summary across iterations.}
#'   \item{all_iterations}{Combined data frame of all summary_stats results.}
#'   \item{metadata}{Info about run parameters.}
#' }
#' @export
repeated_perm_diff_interaction <- function(
  X,
  Y,
  MX,
  MY,
  stratify_cols,
  g_col,
  a_col,
  features = NULL,
  n_iter = 50,
  B = 100,
  seed = 42
) {
  all_results <- list()
  set.seed(seed)

  for (i in seq_len(n_iter)) {
    split <- split_stratified_ancestry_sets(
      X = X,
      Y = Y,
      MX = MX,
      MY = MY,
      stratify_cols = stratify_cols,
      seed = seed + i
    )

    if (nrow(split$test$X) == 0) next

    result <- perm_interaction(
      X = split$test$X,
      Y = split$inference$X,
      MX = split$test$M,
      MY = split$inference$M,
      g_col = g_col,
      a_col = a_col,
      B = B,
      seed = seed
    )

    stats <- result$summary_stats
    if (!is.null(features)) {
      stats <- stats[stats$feature %in% features, ]
    }

    stats$iteration <- i
    all_results[[length(all_results) + 1]] <- stats
  }

  combined <- do.call(rbind, all_results)

  # Aggregate: mean, sd, median p, prob_signif
  agg <- aggregate(
    cbind(T_obs, p_value) ~ feature,
    data = combined,
    FUN = function(x) c(
      mean = mean(x),
      sd = sd(x),
      median = median(x),
      prob_signif = mean(x < 0.05)
    )
  )

  # Flatten matrix columns
  agg_df <- data.frame(
    feature = agg$feature,
    mean_T_obs = agg$T_obs[, "mean"],
    sd_T_obs = agg$T_obs[, "sd"],
    median_p = agg$p_value[, "median"],
    prob_signif = agg$p_value[, "prob_signif"],
    row.names = NULL
  )

  list(
    aggregated  = agg_df,
    all_iterations = combined,
    metadata = list(
      n_iter = n_iter,
      B = B,
      seed = seed,
      features = if (is.null(features)) "all" else features
    )
  )
}
