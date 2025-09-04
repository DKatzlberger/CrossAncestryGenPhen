#' Parallel Subset-Based Interaction Effect Analysis
#'
#' Runs stratified random subsets through a permutation-based interaction 
#' effect pipeline in parallel using C++ core logic. Uses {future} and 
#' {furrr} for parallelization with automatic plan handling.
#'
#' @param X Data frame or matrix of predictors.
#' @param Y Outcome variable or matrix.
#' @param MX Covariates for the test split.
#' @param MY Covariates for the inference split.
#' @param g_col Name of the grouping or genotype column.
#' @param n_iter Number of random subsets to run. Default: 1000.
#' @param B Number of permutations per subset. Default: 1000.
#' @param save_null Store full permutation null matrix? Default: FALSE.
#' @param seed Optional seed for reproducibility. Default: NULL.
#' @param workers Number of parallel workers. Default: availableCores()-1.
#'
#' @return A list with aggregated stats, iteration stats, sample log,
#' permutation nulls (optional), and meta information.
#'
#' @details The function sets a parallel plan with `multisession` and 
#' restores the previous plan automatically.
#'
#' @import future
#' @import furrr
#' @export
subset_perm_interaction_effect <- function(
  X,
  Y,
  MX,
  MY,
  g_col,
  n_iter = 1000,
  B = 1000,
  save_null = FALSE,
  seed = NULL,
  workers = future::availableCores() - 1,
  verbose = TRUE
) {
  # Save current plan and set new one
  orig_plan <- future::plan()
  future::plan(multisession, workers = workers)
  on.exit(future::plan(orig_plan), add = TRUE)

  # Precompute seeds for reproducibility
  seeds <- if (!is.null(seed)) seed + seq_len(n_iter) else rep(list(NULL), n_iter)

  # Prepare arguments for parallelize
  args <- data.frame(
    i = seq_len(n_iter),
    seed_iter = seeds
  )

  # Define the function for one iteration
  run_one <- function(
    i,
    seed_iter
  ) {
    
    # Set RNG 
    if (!is.null(seed_iter)) set.seed(seed_iter)

    # Stratified splits
    split <- split_stratified_ancestry_sets(
      X = X,
      Y = Y,
      MX = MX,
      MY = MY,
      g_col = g_col,
      seed = NULL
    )
    # Track the sample ids
    id <- track_sample_ids(split, i)

    # Run C++ core
    perm_res <- perm_interaction_effect_cpp(
      X = split$test$X,
      Y = split$inference$X,
      MX = split$test$M,
      MY = split$inference$M,
      g_col = g_col,
      B = B,
      seed = NULL
    )

    stats <- perm_res$summary_stats
    stats$iteration <- i

    list(
      stats = stats,
      id = id,
      t_null = if (save_null) perm_res$T_null else NULL
    )
  }

  # Run in parallel
  res <- furrr::future_pmap(
    args,
    run_one, 
    .options = furrr::furrr_options(seed = seed),
    .progress = TRUE
  )

  # Extract and combine
  perm_log <- do.call(rbind, lapply(res, `[[`, "stats"))
  id_log   <- do.call(rbind, lapply(res, `[[`, "id"))
  t_null_log <- if (save_null) lapply(res, `[[`, "t_null") else NULL

  aggregated_perm <- summarize_subsets(
    data = perm_log
  )

  list(
    aggregated_stats = aggregated_perm,
    iteration_stats = perm_log,
    sample_log = id_log,
    permutation_nulls = t_null_log,
    meta = list(
      n_iter = n_iter,
      B = B,
      seed = seed
    )
  )
}

