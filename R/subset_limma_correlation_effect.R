#' Run repeated subset-based limma correlation effects
#'
#' @param X Expression matrix for ancestry X. Rows = samples, columns = genes.
#' @param Y Expression matrix for ancestry Y. Rows = samples, columns = genes.
#' @param MX Metadata for X. Must include group and ancestry columns.
#' @param MY Metadata for Y. Must include group and ancestry columns.
#' @param g_col Name of the genotype column in `X` used for interaction.
#' @param a_col Name of the ancestry column used for stratified splitting.
#' @param covariates Optional vector of covariate column names to adjust for.
#' @param use_voom Logical; whether to use limma-voom (default: TRUE).
#' @param n_iter Integer. Number of iterations to run. Default is 1000.
#' @param method Method for aggregating p-values across iterations. Options: "pearson", "spearman".
#' @param seed Optional integer for random seed to ensure reproducibility.
#' @param verbose Logical; whether to print progress messages.
#'
#' @importFrom furrr future_pmap
#' @importFrom future nbrOfWorkers
#' @importFrom data.table as.data.table rbindlist
#' 
#' @export
subset_limma_correlation_effect <- function(
  X,
  Y,
  MX,
  MY,
  g_col,
  a_col,
  covariates = NULL,
  use_voom = TRUE,
  n_iter = 1000,
  method = c("pearson", "spearman"),
  seed = NULL,
  verbose = TRUE
){

  ## --- Match the method ---
  method <- match.arg(method)


  ## --- Input data structure check ---
  assert_input(
    X = X, 
    Y = Y,
    MX = MX, 
    MY = MY,
    g_col = g_col, 
    a_col = a_col,
    .fun = "subset_limma_correlation_effect"
  )


  ## --- Parallelization setup ---
  n_workers  <- future::nbrOfWorkers()
  message(sprintf("\n[subset_limma_correlation_effect] Workers available: %d", n_workers))


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
    res <- limma_correlation_effect(
      R = split$R$matr,
      X = split$X$matr,
      Y = split$Y$matr,
      MR = split$R$meta,
      MX = split$X$meta,
      MY = split$Y$meta,
      g_col = g_col,
      a_col = a_col,
      covariates = covariates,
      use_voom = use_voom,
      verbose = verbose
    )

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
    .progress = FALSE
  )


  ## --- Extract and combine ---
  res_log <- do.call(rbind, lapply(parallel_res, `[[`, "res"))
  ids_log <- do.call(rbind, lapply(parallel_res, `[[`, "ids"))


  ## --- Aggregation of iterations ---
  agg_log <- summarize_limma_correlation_effect_subsets(
    stats  = res_log,
    method = method
  )

  # Function should return both
  sel_method = agg_log$sel_delta_res
  all_method = agg_log$all_delta_res


  ## --- Return ---
  return(
    list(
      summary_stats = sel_method,
      methods_stats = all_method,
      subsets_stats = res_log,
      samples_stats = ids_log
    )
  )
}