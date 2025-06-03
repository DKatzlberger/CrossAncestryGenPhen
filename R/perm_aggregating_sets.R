
perm_aggregate_sets <- function(
  x,
  T_col,                                            # The test statistic column to aggregate
  feature_col,                                      # Per feature analysis
  agg_fun = mean,                                   # Aggregation function like mean, median, etc.
  alternative = c("two.sided", "greater", "less"),  # What kind of test to perform
  B = 1000,                                         # Number of permutations
  seed = NULL                                       # Random seed for reproducibility
) {
  if (!is.null(seed)) set.seed(seed)

  features <- unique(x[[feature_col]])
  n_features <- length(features)

  # Observed aggregated stat per feature
  obs_vec <- tapply(x[[T_col]], x[[feature_col]], FUN = agg_fun)

  # Prepare matrix to hold permuted stats
  perm_matrix <- matrix(NA_real_, nrow = n_features, ncol = B)
  rownames(perm_matrix) <- names(obs_vec)

  for (f in seq_along(features)) { 
    # Test each feature independently
    feature_name <- features[f]
    T_f <- x[x[[feature_col]] == feature_name, T_col]

    for (b in seq_len(B)) {
      # Shuffle the labels in shuffle_col
      T_perm <- sample(T_f)
      # Aggregate the permuted values using agg_fun
      perm_matrix[f, b] <- agg_fun(T_perm)
    }
  }

  # Empirical p-values (greater than observed)
  p_value <- numeric(n_features)
  for (f in seq_along(features)) {
    T_obs <- obs_vec[f]
    T_null <- perm_matrix[f, ]
    p_value[f] <- compute_empirical_p(T_obs, T_null, alternative = alternative)
  }

  # 5. Adjust for multiple testing (FDR)
  p_adj <- p.adjust(p_value, method = "BH")

  # 6. Return result
  result <- data.frame(
    feature = names(obs_vec),
    T_obs = obs_vec,
    p_value = p_value,
    p_adj = p_adj,
    alternative = alternative,
    row.names = NULL
  )

  return(result)
}