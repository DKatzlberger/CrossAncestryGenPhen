#' Leave-One-Out (LOO) Differential Expression Analysis
#' 
#' @param X Expression matrix for ancestry X. Rows = samples, columns = genes.
#' @param Y Expression matrix for ancestry Y. Rows = samples, columns = genes.
#' @param MX Metadata for X. Must include group and ancestry columns.
#' @param MY Metadata for Y. Must include group and ancestry columns.
#' @param g_col Name of the column indicating group (factor 1).
#' @param a_col Name of the column indicating ancestry (factor 2).
#' @param covariates Optional vector of covariate column names to adjust for.
#' @param verbose Logical, whether to print messages.
#' 
#' @export
loo_interaction_effect <- function(
  X, 
  Y, 
  MX, 
  MY, 
  g_col,
  a_col,
  covariates = NULL,
  use_voom = TRUE,
  verbose = TRUE
) {

  ## --- Input data structure check ---
  assert_input(
    X = X, 
    Y = Y, 
    MX = MX, 
    MY = MY, 
    g_col = g_col, 
    a_col = a_col,
    .fun = "loo_interaction_effect"
  )


  ## --- Initialize results ---
  all_samples  <- c(rownames(X), rownames(Y))
  results_list <- list()


  ## --- Verbose ---
  if (verbose) message("\n[loo_interaction_effect] Removing sample:")


  ## --- Loop over each sample ---‚
  for (s in all_samples){
    if (verbose) message(sprintf("%s", s))
    

    # Drop sample s from both X/Y and MX/MY if present
    X_sub  <- X[setdiff(rownames(X), s), , drop = FALSE]
    Y_sub  <- Y[setdiff(rownames(Y), s), , drop = FALSE]
    MX_sub <- MX[setdiff(rownames(MX), s), , drop = FALSE]
    MY_sub <- MY[setdiff(rownames(MY), s), , drop = FALSE]


    # Run limma on the reduced dataset
    res <- limma_interaction_effect(
      X = X_sub,
      Y = Y_sub,
      MX = MX_sub,
      MY = MY_sub,
      g_col = g_col,
      a_col = a_col,
      covariates = covariates,
      use_voom = use_voom,
      verbose = FALSE
    )

    # Annotate which sample was removed
    res$sample_removed <- s
    results_list[[s]] <- res
  }


  ## --- Combine results into one data frame ---
  all_res <- do.call(rbind, results_list)


  ## --- Return ---
  return(all_res)
}