test_that("perm_diff_interaction runs correctly on simulated data", {
  set.seed(1)
  n_per_group <- 25
  p <- 2000
  
  # Simulate expression data
  X <- matrix(rnorm(n_per_group * 2 * p), nrow = n_per_group * 2, ncol = p)
  Y <- matrix(rnorm(n_per_group * 2 * p), nrow = n_per_group * 2, ncol = p)
  colnames(X) <- colnames(Y) <- paste0("Gene_", seq_len(p))
  
  # Simulate metadata
  MX <- data.frame(
    condition = factor(rep(c("A", "B"), each = n_per_group)),
    ancestry = "EUR"
  )
  MY <- data.frame(
    condition = factor(rep(c("A", "B"), each = n_per_group)),
    ancestry = "AFR"
  )
  
  # Run the function
  result <- perm_diff_interaction(
    X = X,
    Y = Y,
    MX = MX,
    MY = MY,
    g_col = "condition",
    a_col = "ancestry",
    B = 100,
    seed = 42
  )
  
  # Check structure
  expect_type(result, "list")
  expect_true(all(c("summary_stats", "T_boot", "B_used", "converged") %in% names(result)))
  
  # Check summary table
  expect_s3_class(result$summary_stats, "data.frame")
  expect_equal(nrow(result$summary_stats), p)
  expect_true(all(result$summary_stats$p_value >= 0 & result$summary_stats$p_value <= 1))
  expect_true(all(is.finite(result$summary_stats$T_obs)))
  
  # Check permutation matrix
  expect_true(is.matrix(result$T_boot))
  expect_equal(ncol(result$T_boot), p)
  expect_gte(nrow(result$T_boot), 1)
})
