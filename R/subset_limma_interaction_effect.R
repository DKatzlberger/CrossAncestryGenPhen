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
  covariates = NULL,
  use_voom = TRUE,
  n_iter = 1000,
  workers = future::availableCores() - 1,
  seed = NULL,
  verbose = TRUE
) {

  ## --- Input data structure check ---
  assert_input(
    X = X, 
    Y = Y, 
    MX = MX, 
    MY = MY, 
    g_col = g_col, 
    a_col = a_col
  )

  ## --- Parallelization steps ---
  orig_plan <- future::plan()
  future::plan(multisession, workers = workers)
  on.exit(future::plan(orig_plan), add = TRUE)

  ## --- Seeds for reproducibility ---
  seeds <- if (!is.null(seed)) seed + seq_len(n_iter) else rep(list(NULL), n_iter)


  ## --- Prepare arguments for parallel execution ---
  args <- data.frame(
    i = seq_len(n_iter),
    seed_iter = seeds
  )


  ## --- Define the function for one iteration ---
  run_one <- function(
    i,
    seed_iter
  ) {
    
    # Seed 
    if (!is.null(seed_iter)) set.seed(seed_iter)

    # Stratified splits
    split <- split_stratified_ancestry_sets(
      X = X,
      Y = Y,
      MX = MX,
      MY = MY,
      g_col = g_col,
      a_col = a_col,
      seed = seed_iter,
      verbose = verbose
    )
      
    # Store sample ids
    id <- track_sample_ids(split, i)

    # Run limma on subset
    res <- limma_interaction_effect(
      X = split$X$counts,
      Y = split$Y$counts,
      MX = split$X$meta,
      MY = split$Y$meta,
      g_col = g_col,
      a_col = a_col,
      covariates = covariates,
      use_voom = TRUE,
      verbose = verbose
    )

    # Filter interaction coef
    res <- subset(res, coef_type == "interaction")

    # Add the iteration number
    res$iteration <- i

    list(
      res = res,
      ids = id
    )
  }


  ## --- Run in parallel ---
  parallel_res <- furrr::future_pmap(
    args,
    run_one, 
    .options = furrr::furrr_options(seed = seed),
    .progress = TRUE
  )


  ## --- Extract and combine ---
  res_log <- do.call(rbind, lapply(parallel_res, `[[`, "res"))
  ids_log <- do.call(rbind, lapply(parallel_res, `[[`, "ids"))


  ## --- Aggregation of iterations ---
  agg_res <- summarize_subsets(
    data = res_log
  )

  ## --- Return ---
  return(
    list(
      summary_stats   = agg_res,
      iteration_stats = res_log,
      sample_stats    = ids_log
    )
  )
}

