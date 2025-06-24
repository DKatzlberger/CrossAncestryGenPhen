#' Permutation Test for Prediction Performance Differences
#'
#' This function compares model prediction performance between two datasets 
#' (e.g., different ancestries) using a reference training set. It trains 
#' a classification model using tidymodels, evaluates it on test and 
#' inference datasets, and assesses the statistical significance of 
#' the performance difference via permutation testing.
#'
#' @param X Matrix of predictors for the test group (ancestry 1); samples x features
#' @param Y Matrix of predictors for the inference group (ancestry 2)
#' @param R Matrix of predictors for the reference training group
#' @param MX Data frame of metadata for `X`, must include the outcome column
#' @param MY Data frame of metadata for `Y`
#' @param MR Data frame of metadata for `R`
#' @param g_col Character. Name of the column in `MX`, `MY`, `MR` containing 
#'   the binary outcome (must have exactly 2 levels).
#' @param method Character. Which model to use: `"glmnet"` (default) or `"rf"` (random forest).
#' @param metric Character. Performance metric to optimize and test, currently supports `"roc_auc"` only.
#' @param cv_folds Integer. Number of cross-validation folds (default = 5).
#' @param tune_len Integer. Number of levels per hyperparameter in grid search (default = 10).
#' @param max_iter Integer. Not used currently but reserved for future logic (default = 1000).
#' @param B Integer. Number of permutations to run (default = 1000).
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A list containing:
#' \describe{
#'   \item{summary_stats}{A data frame with the observed statistic, group metrics, and p-value}
#'   \item{T_null}{Null distribution of the test statistic from permutations}
#'   \item{B_used}{Number of successful permutations}
#' }
#'
#'
#' @export
#'
#' @importFrom recipes recipe
#' @importFrom parsnip logistic_reg rand_forest set_engine set_mode fit
#' @importFrom workflows workflow add_recipe add_model
#' @importFrom tune tune_grid select_best control_grid finalize_workflow
#' @importFrom dials parameters grid_regular
#' @importFrom rsample vfold_cv
#' @importFrom yardstick roc_auc metric_set
#' @importFrom hardhat extract_parameter_set_dials
#' @importFrom magrittr %>%
perm_prediction_difference <- function(
  X,
  Y,
  R,
  MX,
  MY,
  MR,
  g_col,
  method = c("glmnet", "rf"),
  metric = c("roc_auc"),
  cv_folds = 5,
  tune_len = 10,
  max_iter = 1000,
  B = 1000,
  seed = NULL
){
  # Set the seed for reproducibility
  if (!is.null(seed)) set.seed(seed)

  # Match method and metric arguments
  method <- match.arg(method)
  metric <- match.arg(metric)
  
  # Prepare data to fit to the model
  df_X <- as.data.frame(X)
  df_Y <- as.data.frame(Y)
  df_R <- as.data.frame(R)
  # Get true labels as factor with matching levels
  if (length(levels(MX[[g_col]])) != 2) {
    stop("Response variable must have exactly 2 levels.")
  }
  g1 <- levels(factor(MX[[g_col]]))[1]
  g2 <- levels(factor(MX[[g_col]]))[2]
  # Relevel
  y_X <- factor(MX[[g_col]], levels = c(g1, g2))
  y_Y <- factor(MY[[g_col]], levels = c(g1, g2))
  y_R <- factor(MR[[g_col]], levels = c(g1, g2))

  # Prepare data frames
  df_X$.y <- y_X
  df_Y$.y <- y_Y
  df_R$.y <- y_R

  # Recipe for training the model
  rec <- recipes::recipe(.y ~ ., data = df_R)

  # Model specifications based on the method
  model_spec <- switch(
    method,
    glmnet = {
      parsnip::logistic_reg(
        penalty = tune(),     
        mixture = tune()       
      ) %>%
      parsnip::set_engine("glmnet") %>%
      parsnip::set_mode("classification")
    },
    rf = {
      parsnip::rand_forest(
        mtry = tune(),         
        trees = 500,
        min_n = tune()         
      ) %>%
      parsnip::set_engine("ranger") %>%
      parsnip::set_mode("classification")
    },
    stop("Unsupported method")
  )

  # Define the worflow for model training
  wf <- workflows::workflow() %>%
    workflows::add_recipe(rec) %>%
    workflows::add_model(model_spec)

  # Define the metric to score cv
  predefined_metrics <- yardstick::metric_set(
    yardstick::roc_auc
  )

  # Grid
  grid <- dials::grid_regular(
    hardhat::extract_parameter_set_dials(model_spec),
    levels = tune_len
  )

  # Define the cv for hyperparameter tuning
  folds <- rsample::vfold_cv(
    df_R, 
    v = cv_folds, 
    strata = .y
  )

  # Tune model
  tuned <- tune::tune_grid(
    wf,
    resamples = folds,
    grid = grid,
    metrics = predefined_metrics,
    control = tune::control_grid(verbose = TRUE)
  )

  # Select the best params and refit the model
  best_params <- tune::select_best(tuned, metric = metric)
  final_wf <- tune::finalize_workflow(wf, best_params)
  fitted_model <- parsnip::fit(final_wf, data = df_R)


  # Prediction setup
  metric_fn <- switch(
    metric,
    "roc_auc" = yardstick::roc_auc,
    stop("Unsupported metric. Only 'roc_auc' is supported.")
  )

  # Train prediction (Done once)
  prob_R <- predict(fitted_model, new_data = df_R, type = "prob")
  prob_R$.y <- df_R$.y
  metric_R <- metric_fn(
    data = prob_R, 
    truth = .y, 
    !!rlang::sym(paste0(".pred_", g1))
  )[[".estimate"]] 

  # Test statisitc
  # Observed prediction difference
  prob_X <- predict(fitted_model, new_data = df_X, type = "prob")
  prob_X$.y <- df_X$.y
  metric_XR <- metric_fn(
    data = prob_X, 
    truth = .y,
    !!rlang::sym(paste0(".pred_", g1))
  )[[".estimate"]] 

  prob_Y <- predict(fitted_model, new_data = df_Y, type = "prob")
  prob_Y$.y <- df_Y$.y
  metric_YR <- metric_fn(
    data = prob_Y, 
    truth = .y, 
    !!rlang::sym(paste0(".pred_", g1))
  )[[".estimate"]] 

  # Delta observed
  T_obs <- metric_YR - metric_XR

  # Permutation setup
  XY <- rbind(X, Y)
  MXY <- rbind(MX, MY)

  # Sample sizes
  n1 <- nrow(X)
  n2 <- nrow(Y)

  # Initilaize perm vector
  T_perm <- rep(NA_real_, B)

  message("Starting permutation test with ", B, " iterations...")
  for (b in 1:B) {

    # Randomly split into two pseudo ancestries
    perm_idx <- sample(1:nrow(XY))
    p1_idx <- perm_idx[1:n1]
    p2_idx <- perm_idx[(n1 + 1):(n1 + n2)]

    # Convert the pseudo ancestries to data frames
    df_p1 <- as.data.frame(XY[p1_idx, , drop = FALSE])
    df_p2 <- as.data.frame(XY[p2_idx, , drop = FALSE])

    # Add the true labels as factors
    df_p1$.y <- factor(MXY[[g_col]][p1_idx], levels = c(g1, g2))
    df_p2$.y <- factor(MXY[[g_col]][p2_idx], levels = c(g1, g2))

    # Predict probabilities for pseudo ancestries
    prob_p1 <- predict(fitted_model, new_data = df_p1, type = "prob")
    prob_p1$.y <- df_p1$.y
    metric_p1 <- metric_fn(
      data = prob_p1, 
      truth = .y, 
      !!rlang::sym(paste0(".pred_", g1))
    )[[".estimate"]]

    prob_p2 <- predict(fitted_model, new_data = df_p2, type = "prob")
    prob_p2$.y <- df_p2$.y
    metric_p2 <- metric_fn(
      data = prob_p2, 
      truth = .y, 
      !!rlang::sym(paste0(".pred_", g1))
    )[[".estimate"]]

    if (any(is.na(c(metric_p1, metric_p2)))) {
      T_perm[b] <- NA
      next
    }
    # Calculate the difference for the permutation
    T_perm[b] <- metric_p2 - metric_p1
  }

  # Remove invalid permutations
  T_perm <- T_perm[!is.na(T_perm)]
  B_used <- length(T_perm)

  if (B_used == 0) stop("No valid permutations completed.")
  p_vals_emp <- (sum(abs(T_perm) >= abs(T_obs)) + 1) / (B_used + 1)

  # Parametric p-values assuming normal distribution
  mu <- mean(T_perm)
  sigma <- sd(T_perm)
  z <- (T_obs - mu) / sigma
  p_vals_param <- 2 * (1 - pnorm(abs(z)))

  summary_stats <- data.frame(
    feature = "Global",
    T_obs = T_obs,
    R = metric_R,
    XR = metric_XR,
    YR = metric_YR,
    p_param_value = p_vals_param,
    p_emp_value = p_vals_emp,
    row.names = NULL
  )

  return(
    list(
      summary_stats = summary_stats,
      T_null = matrix(T_perm, ncol = 1, dimnames = list(NULL, "T_null")),
      B_used = B_used
    )
  )
}
