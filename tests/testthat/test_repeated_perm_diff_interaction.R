test_that("repeated_perm_diff_interaction returns expected structure and stats", {
  set.seed(123)

  # Simulate expression data
  p <- 100  # number of genes
  n_EUR <- 120
  n_AFR <- 40

  X <- matrix(rnorm(n_EUR * p), nrow = n_EUR, ncol = p)
  Y <- matrix(rnorm(n_AFR * p), nrow = n_AFR, ncol = p)
  colnames(X) <- colnames(Y) <- paste0("Gene_", seq_len(p))

  # Simulate imbalanced metadata for EUR and balanced for AFR
  MX <- data.frame(
    condition = factor(rep(c("A", "B"), each = n_EUR / 2)),
    sex = factor(c(rep("F", 70), rep("M", 50))),  # fewer A.M and B.M
    ancestry = "EUR"
  )

  MY <- data.frame(
    condition = factor(rep(c("A", "B"), each = n_AFR / 2)),
    sex = factor(rep(c("M", "F"), length.out = n_AFR)),  # balanced
    ancestry = "AFR"
  )

  # Subset of genes to test
  features_to_test <- paste0("Gene_", 1:10)

  # Run the repeated permutation test and allow expected warning
  expect_warning(
    result <- repeated_perm_diff_interaction(
      X = X,
      Y = Y,
      MX = MX,
      MY = MY,
      stratify_cols = c("condition", "sex"),
      g_col = "condition",
      a_col = "ancestry",
      features = features_to_test,
      n_iter = 5,
      B = 10,
      seed = 42
    ),
    regexp = "could not be matched",
    info = "Expected warning for unmatched strata due to imbalance"
  )

  # Check structure of the returned result
  expect_type(result, "list")
  expect_named(result, c("aggregated", "all_iterations", "metadata"))

  # Check structure of aggregated output
  agg <- result$aggregated
  expect_true(is.data.frame(agg))
  expect_equal(ncol(agg), 5)
  expect_true(all(c("feature", "mean_T_obs", "sd_T_obs", "median_p", "prob_signif") %in% colnames(agg)))
  expect_equal(length(unique(agg$feature)), length(features_to_test))

  # Check all_iterations
  all_iter <- result$all_iterations
  expect_true(is.data.frame(all_iter))
  expect_true(all(c("feature", "T_obs", "ci_lower", "ci_upper", "p_value", "iteration") %in% colnames(all_iter)))
  expect_true(all(all_iter$feature %in% features_to_test))

  # Check that prob_signif is between 0 and 1
  expect_true(all(result$aggregated$prob_signif >= 0 & result$aggregated$prob_signif <= 1))
})
