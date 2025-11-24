#' Run repeated subset-based logistic prediction effects
#'
#' @param X Matrix or data frame of features for subset X.
#' @param Y Matrix or data frame of features for subset Y.
#' @param MX Metadata for X (must include grouping columns).
#' @param MY Metadata for Y (must include grouping columns).
#' @param g_col Name of the ancestry/group variable in metadata.
#' @param a_col Name of the label/outcome variable.
#' @param n_folds Number of folds for logistic model cross-validation.
#' @param n_models Number of logistic models to average per fold.
#' @param maxit Maximum iterations for logistic regression (optional).
#' @param n_iter Number of repeated stratified resampling iterations.
#' @param method Performance metric to compute; currently only `"auc"`.
#' @param seed Optional integer seed for reproducibility.
#' @param verbose Logical; print progress messages.
#'
#' @importFrom furrr future_pmap
#' @importFrom future nbrOfWorkers
#' @importFrom data.table as.data.table rbindlist
#'
#' @export
subset_logistic_prediction_effect <- function(
  X,
  Y,
  MX,
  MY,
  g_col,
  a_col,
  n_folds,
  n_models,
  maxit = NULL,
  n_iter = 1000,
  method = c("auc", "logloss"),
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
    .fun = "subset_logistic_prediction_effect"
  )


  ## --- Parallelization setup ---
  n_workers  <- future::nbrOfWorkers()
  message(sprintf("\n[subset_logistic_prediction_effect] Workers available: %d", n_workers))


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

    # Run logistic regression on subset
    res <- logistic_prediction_effect(
      R = split$R$matr,
      X = split$X$matr,
      Y = split$Y$matr,
      MR = split$R$meta,
      MX = split$X$meta,
      MY = split$Y$meta,
      g_col = g_col,
      a_col = a_col,
      n_folds = n_folds,
      n_models = n_models,
      maxit = maxit,
      seed = seed_iter,
      verbose = verbose
    )

    # Disect into summary_stats feature_stats
    summary_stats <- res$summary_stats
    feature_stats <- res$feature_stats

    # Add the iteration number
    summary_stats$iteration <- i
    feature_stats$iteration <- i

    list(
      summary_stats = summary_stats,
      feature_stats = feature_stats,
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
  summary_stats_log <- do.call(rbind, lapply(parallel_res, `[[`, "summary_stats"))
  features_stats_log <- do.call(rbind, lapply(parallel_res, `[[`, "feature_stats"))
  ids_log <- do.call(rbind, lapply(parallel_res, `[[`, "ids"))

  ## --- Aggregation of iterations ---
  agg_log <- summarize_logistic_prediction_effect_subsets(
    stats  = summary_stats_log,
    method = method,
    by = NULL
  )

  # Function should return both
  sel_method = agg_log$sel_delta_res
  all_method = agg_log$all_delta_res


  ## --- Return ---
  return(
    list(
      summary_stats = sel_method,
      methods_stats = all_method,
      feature_stats = features_stats_log,
      subsets_stats = summary_stats_log,
      samples_stats = ids_log
    )
  )
}
