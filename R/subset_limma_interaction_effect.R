#' Subset Limma Interaction Effect with Stratified Cross-Validation
#'
#' Runs stratified train/test splits with ancestry sets and estimates 
#' interaction effects using limma in parallel across multiple iterations.
#' Returns aggregated statistics, iteration-level stats, and sample logs.
#'
#' @param X A data frame or matrix of features (genotypes or other predictors).
#' @param Y A response vector or matrix.
#' @param MX Additional metadata or covariates for `X`.
#' @param MY Additional metadata or covariates for `Y`.
#' @param g_col Name of the genotype column in `X` used for interaction.
#' @param a_col Name of the ancestry column used for stratified splitting.
#' @param n_iter Integer. Number of iterations to run. Default is 1000.
#' @param seed Optional integer for random seed to ensure reproducibility.
#' @param workers Number of parallel workers to use. 
#'   Defaults to `future::availableCores() - 1` (reserving one core).
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{aggregated_stats}{Data frame with aggregated summary statistics across iterations.}
#'   \item{iteration_stats}{Data frame with iteration-level statistics.}
#'   \item{sample_log}{Data frame logging sample IDs used in each split.}
#'   \item{meta}{List with metadata about the run (e.g., `n_iter` and `seed`).}
#' }
#'
#' @export
subset_limma_interaction_effect <- function(
  X,
  Y,
  MX,
  MY,
  g_col,
  a_col,
  n_iter = 1000,
  seed = NULL,
  workers = future::availableCores() - 1  # optional: reserve 1 core
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
    limma_res <- limma_interaction_effect(
      X = split$test$X,
      Y = split$inference$X,
      MX = split$test$M,
      MY = split$inference$M,
      g_col = g_col,
      a_col = a_col
    )

    stats <- limma_res$summary_stats
    stats$iteration <- i

    list(
      stats = stats,
      id = id
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
  limma_log <- do.call(rbind, lapply(res, `[[`, "stats"))
  id_log   <- do.call(rbind, lapply(res, `[[`, "id"))

  aggregated_limma <- summarize_subsets(
    data = limma_log,
    iter_col = "iteration"
  )

  list(
    aggregated_stats = aggregated_limma,
    iteration_stats = limma_log,
    sample_log = id_log,
    meta = list(
      n_iter = n_iter,
      seed = seed
    )
  )
}

