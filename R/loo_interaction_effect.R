loo_interaction_effect <- function(
  X, 
  Y, 
  MX, 
  MY, 
  g_col,
  a_col,
  covariates,
  use_voom,
  verbose
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
  results_list <- list()

  all_samples  <- c(rownames(X), rownames(Y))
  all_meta     <- rbind(MX, MY)
  all_ancestry <- all_meta[[a_col]]
  

  ## --- Drop individual samples ---
  limma_res <- limma_interaction_effect(
    X = ,
    Y = ,
    MX = ,
    MY = ,
    g_col = g_col,
    a_col = a_col,
    covariates = covariates,
    use_voom = ,
    verbose = FALSE
  )
}
