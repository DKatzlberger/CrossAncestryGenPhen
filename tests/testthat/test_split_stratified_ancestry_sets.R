test_that("split_stratified_ancestry_sets works correctly on simulated data", {
  set.seed(42)

  # Simulate data
  n_X <- 150  # Overrepresented ancestry
  n_Y <- 50   # Underrepresented ancestry
  p <- 100    # Number of genes

  X <- matrix(rnorm(n_X * p), nrow = n_X, ncol = p)
  Y <- matrix(rnorm(n_Y * p), nrow = n_Y, ncol = p)

  MX <- data.frame(
    condition = sample(c("A", "B"), n_X, replace = TRUE),
    sex = sample(c("M", "F"), n_X, replace = TRUE)
  )

  MY <- data.frame(
    condition = sample(c("A", "B"), n_Y, replace = TRUE),
    sex = sample(c("M", "F"), n_Y, replace = TRUE)
  )

  # Run the function
  result <- split_stratified_ancestry_sets(
    X = X,
    Y = Y,
    MX = MX,
    MY = MY,
    stratify_cols = c("condition", "sex"),
    seed = 42
  )

  # Check structure
  expect_type(result, "list")
  expect_named(result, c("train", "test", "inference", "strata_info"))

  # Check each split
  for (split in c("train", "test", "inference")) {
    X_mat <- result[[split]]$X
    M_df  <- result[[split]]$M

    expect_true(is.matrix(X_mat), info = paste(split, "$X is not a matrix"))
    expect_true(is.data.frame(M_df), info = paste(split, "$M is not a data.frame"))
    expect_equal(nrow(X_mat), nrow(M_df), info = paste(split, "$X/$M row mismatch"))
  }

  # Train size should match or be less than inference size
  expect_lte(nrow(result$train$X), nrow(result$inference$X))

  # Train + test rows should equal total usable X subset
  total_X_used <- nrow(result$train$X) + nrow(result$test$X)
  expect_lte(total_X_used, n_X)
})
