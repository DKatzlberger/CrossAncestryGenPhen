#' Logistic prediction with glmnet
#'
#' Fits elastic-net logistic models on R, tests on X and Y,
#' checks leakage, tunes hyperparameters, and returns
#' predictions and model coefficients.
#'
#' @param R Expression matrix for refeernce ancestry. Rows = samples, columns = genes.
#' @param X Expression matrix for ancestry X. Rows = samples, columns = genes.
#' @param Y Expression matrix for ancestry Y. Rows = samples, columns = genes.
#' @param MR Metadata for R. Must include group and ancestry columns.
#' @param MX Metadata for X. Must include group and ancestry columns.
#' @param MY Metadata for Y. Must include group and ancestry columns.
#' @param g_col Name of the column indicating group (factor 1).
#' @param a_col Name of the column indicating ancestry (factor 2).
#' @param n_folds CV folds
#' @param n_models Grid size
#' @param maxit Max iter
#' @param seed RNG seed
#' @param verbose Print msgs
#'
#' @return A list with:
#' \describe{
#'   \item{summary_stats}{Predictions + labels}
#'   \item{feature_stats}{Non-zero glmnet coefs}
#' }
#' @import recipes
#' @import parsnip
#' @import workflows
#' @import tune
#' @import rsample
#' @import dials
#' @import yardstick
#' @import broom
#' @importFrom generics tidy
#' @importFrom stats cor
#'
#' @export
logistic_prediction_effect <- function(
  R,
  X,
  Y,
  MR,
  MX,
  MY,
  g_col,
  a_col,
  n_folds,
  n_models,
  maxit = NULL,
  seed = NULL,
  verbose = TRUE
){
  ## --- Seed ---
  if(!is.null(seed)) set.seed(seed)


  ## --- Input data structure check ---
  assert_input(
    R = R,
    X = X, 
    Y = Y,
    MR = MR,
    MX = MX, 
    MY = MY,
    g_col = g_col, 
    a_col = a_col,
    .fun = "logistic_prediction_effect"
  )


  ## --- Check data leakage ---
  if (length(intersect(rownames(R), rownames(X))) > 0) stop("Data leakage: R and X share rownames.")
  if (length(intersect(rownames(R), rownames(Y))) > 0) stop("Data leakage: R and Y share rownames.")
  if (length(intersect(rownames(X), rownames(Y))) > 0) stop("Data leakage: X and Y share rownames.")


  ## --- Ancestry validation ---
  a_1   <- unique(MX[[a_col]])
  a_2   <- unique(MY[[a_col]])
  a_1_1 <- unique(MR[[a_col]])
  if (a_1 != a_1_1) stop("[logistic_prediction_effect] Ancestry level must be the same in Reference (R) as in Subset (X).")

  ## --- Validation ---
  expr_list <- list(R = R, X = X, Y = Y)
  meta_list <- list(R = MR, X = MX, Y = MY)

  frames <- lapply(names(expr_list), function(a_name) {

    matr <- expr_list[[a_name]]
    meta <- meta_list[[a_name]]

    if (!identical(rownames(matr), rownames(meta))) {
      stop(sprintf("[logistic_prediction_effect] Matrix and meta rownames must match exactly."))
    }

    ## --- Ensure 2-level group ----
    if (!is.factor(meta[[g_col]])) stop("[logistic_prediction_effect] g_col is not a factor.")
    g_levels <- levels(meta[[g_col]])
    a_levels <- unique(meta[[a_col]])

    if (length(g_levels) != 2) stop(sprintf("[logistic_prediction_effect] Function currently supports only 2 groups (two levels in g_col)."))
    if (length(a_levels) != 1) stop(sprintf("[logistic_prediction_effect] Function currently supports only 1 ancestry (one level in a_col)."))

    g_1 <- g_levels[1]
    g_2 <- g_levels[2]
    a_1 <- a_levels[1]

    ## --- Create group ---
    meta[["groups"]] <- factor(
      paste(
        meta[[g_col]], 
        meta[[a_col]], 
        sep = "."
      ), 
      levels = c(
        paste(g_1, a_1, sep = "."),  
        paste(g_2, a_1, sep = ".")
      )
    )

    ## --- Frames with label ---
    prediction_frame <- cbind(meta[ , "groups", drop = FALSE], matr)


    ## --- Summary header ---
    groups_levels <- levels(meta$groups)
    summary_frame <- data.frame(
      coef_id   = paste0("relationship_", a_name),
      coef_type = "relationship",
      contrast  = paste0(groups_levels[2], " - ", groups_levels[1]),
      g_1       = g_1,
      g_2       = g_2,
      a_1       = a_1,
      a_2       = a_2,
      row.names = NULL
    )

    ## --- Return ---
    return(
      list(
        data = prediction_frame,
        meta = summary_frame
      )
    )
  })
  prediction_frames <- lapply(frames, `[[`, "data")
  summary_frames    <- lapply(frames, `[[`, "meta")

  names(prediction_frames) <- names(expr_list)
  names(summary_frames)    <- names(expr_list)


  ## --- Label leakage ---
  check_label_leakage <- function(df, label_col = "groups") {

    y <- df[[label_col]]
    feature_names <- setdiff(colnames(df), label_col)

    ## Exact duplicate columns
    exact_dupes <- feature_names[sapply(feature_names, function(f)
      identical(df[[f]], y)
    )]

    ## Perfect correlation (for numeric features)
    perfect_corr <- c()
    if (is.factor(y) && length(levels(y)) == 2) {
      y_num <- as.numeric(y) - 1
      perfect_corr <- feature_names[sapply(feature_names, function(f) {
        x <- df[[f]]
        if (is.numeric(x)) {
          val <- suppressWarnings(cor(x, y_num))
          !is.na(val) && abs(val) == 1
        } else FALSE
      })]
    }

    leaks <- unique(c(exact_dupes, perfect_corr))

    ## Return 
    list(
      leak = length(leaks) > 0,
      features = leaks
    )
  }

  # Run label leakage checks
  res_R <- check_label_leakage(prediction_frames$R)
  res_X <- check_label_leakage(prediction_frames$X)
  res_Y <- check_label_leakage(prediction_frames$Y)

  leakage_detected <- res_R$leak || res_X$leak || res_Y$leak
  leaking_features <- unique(c(res_R$features, res_X$features, res_Y$features))
  if (leakage_detected) {
    message(">>> LABEL LEAKAGE DETECTED <<<")
    message("Leaking features: ", paste(leaking_features, collapse = ", "))
    stop()
  }


  ## --- Model specification ---
  features <- colnames(prediction_frames$R)[-1]
  form_str <- paste("groups ~ .")
  recipe   <- recipes::recipe(groups ~ ., data = prediction_frames$R)

  # Elastic-net
  model_spec <- logistic_reg(
    penalty = tune(), 
    mixture = tune()    
  ) %>%
  set_mode("classification")
  
  if (!is.null(maxit)) {
    model_spec <- model_spec %>% set_engine("glmnet", maxit = maxit)
  } else {
    model_spec <- model_spec %>% set_engine("glmnet")
  }


  ## --- Workflow ---
  workflow <- workflow() %>%
    add_recipe(recipe) %>%
    add_model(model_spec)


  ## --- Cross-validation ---
  folds <- vfold_cv(
    prediction_frames$R, 
    v = n_folds, 
    strata = groups
  )

  # Grid 
  grid <- grid_space_filling(
    extract_parameter_set_dials(model_spec),
    size = n_models
  )


  ## --- Hyperparameter tuning (penalty + mixture) ---
  tune_result <- workflow %>% 
    tune_grid(
      resamples = folds,
      grid      = grid,
      metrics   = metric_set(roc_auc),
      control   = control_grid(verbose = FALSE)
    )
  
  # Extract best parameters
  best <- select_best(tune_result, metric = "roc_auc")


  ## --- Verbose message ---
  if (verbose) {
    message("\nHyperparameter optimization:")
    message(sprintf("%-20s  %s %d features", "Formula:", form_str, length(features)))
    message(sprintf("%-20s  %s", "Groups:", paste(levels(prediction_frames$R$groups), collapse = "  ")))
    message(sprintf("%-20s  %s", "Parameters:", paste("Lambda:", signif(best$penalty, 4), " Alpha:", signif(best$mixture, 4))))
  }


  ## --- Training step ---
  final_wf  <- finalize_workflow(workflow, best)
  final_fit <- final_wf %>% fit(data = prediction_frames$R)


  ## --- Predcition: subset X, inference Y ---
  pos_class <- levels(prediction_frames$R$groups)[2]
  pred_col  <- paste0(".pred_", pos_class)

  pred_R <- predict(final_fit, prediction_frames$R, type = "prob")
  pred_X <- predict(final_fit, prediction_frames$X, type = "prob")
  pred_Y <- predict(final_fit, prediction_frames$Y, type = "prob")

  # Attach predictions and sample ids
  pred_R$group <- prediction_frames$R$groups
  pred_X$group <- prediction_frames$X$groups
  pred_Y$group <- prediction_frames$Y$groups

  pred_R$sample_id <- rownames(R)
  pred_X$sample_id <- rownames(X)
  pred_Y$sample_id <- rownames(Y)

  logloss_R <- mn_log_loss_vec(truth = prediction_frames$R$groups, estimate = pred_R[[pred_col]], event_level = "second")
  logloss_X <- mn_log_loss_vec(truth = prediction_frames$X$groups, estimate = pred_X[[pred_col]], event_level = "second")
  logloss_Y <- mn_log_loss_vec(truth = prediction_frames$Y$groups, estimate = pred_Y[[pred_col]], event_level = "second")

  auc_R <- roc_auc_vec(truth = prediction_frames$R$groups, estimate = pred_R[[pred_col]], event_level = "second")
  auc_X <- roc_auc_vec(truth = prediction_frames$X$groups, estimate = pred_X[[pred_col]], event_level = "second")
  auc_Y <- roc_auc_vec(truth = prediction_frames$Y$groups, estimate = pred_Y[[pred_col]], event_level = "second")


  ## --- Verbose message ---
  if (verbose) {
    message("\nPrediction performance:")
    message(sprintf("%-20s  logLoss %.4f   AUC %.4f ", paste0("Reference R (", unique(summary_frames$R$a_1), "):"), logloss_R, auc_R))
    message(sprintf("%-20s  logLoss %.4f   AUC %.4f ", paste0("Subset    X (", unique(summary_frames$R$a_1), "):"), logloss_X, auc_X))
    message(sprintf("%-20s  logLoss %.4f   AUC %.4f ", paste0("Inference Y (", unique(summary_frames$R$a_2), "):"), logloss_Y, auc_Y))
  }


  ## --- Coefficients ---
  coef <- suppressMessages(tidy(final_fit))
  coef <- coef[coef$term != "(Intercept)", ]
  feature_stats <- data.frame(
    coef_id     = summary_frames$R$coef_id,
    coef_type   = summary_frames$R$coef_type,
    contrast    = summary_frames$R$contrast,
    g_1         = summary_frames$R$g_1,
    g_2         = summary_frames$R$g_2,
    a_1         = summary_frames$R$a_1,
    a_2         = summary_frames$R$a_2,
    feature     = coef$term,
    estimate    = coef$estimate,
    row.names   = NULL
  )


  ## --- Summary stats ---
  make_summary <- function(header, pred_df, pred_col) {
    data.frame(
      coef_id    = header$coef_id,
      coef_type  = header$coef_type,
      contrast   = header$contrast,
      g_1        = header$g_1,
      g_2        = header$g_2,
      a_1        = header$a_1,
      a_2        = header$a_2,
      sample_id  = pred_df$sample_id,
      true       = pred_df$group,
      prob       = pred_df[[pred_col]],
      row.names  = NULL
    )
  }

  summary_X <- make_summary(summary_frames$X, pred_X, pred_col)
  summary_Y <- make_summary(summary_frames$Y, pred_Y, pred_col)
  summary_stats <- rbind(summary_X, summary_Y)


  ## --- Return ---
  return(
    list(
      summary_stats = summary_stats,
      feature_stats = feature_stats
    )
  )
}