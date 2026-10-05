forest_prediction_effect <- function(
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
    .fun = "forest_prediction_effect"
  )


  ## --- Check data leakage ---
  if (length(intersect(rownames(R), rownames(X))) > 0) stop("Data leakage: R and X share rownames.")
  if (length(intersect(rownames(R), rownames(Y))) > 0) stop("Data leakage: R and Y share rownames.")
  if (length(intersect(rownames(X), rownames(Y))) > 0) stop("Data leakage: X and Y share rownames.")


  ## --- Ancestry validation ---
  A_1   <- unique(MX[[a_col]])
  A_2   <- unique(MY[[a_col]])
  A_1_1 <- unique(MR[[a_col]])
  if (A_1 != A_1_1) stop("[forest_prediction_effect] Ancestry level must be the same in Reference (R) as in Subset (X).")

  ## --- Validation ---
  expr_list <- list(R = R, X = X, Y = Y)
  meta_list <- list(R = MR, X = MX, Y = MY)

  frames <- lapply(names(expr_list), function(a_name) {

    matr <- expr_list[[a_name]]
    meta <- meta_list[[a_name]]

    if (!identical(rownames(matr), rownames(meta))) {
      stop(sprintf("[forest_prediction_effect] Matrix and meta rownames must match exactly."))
    }

    ## --- Ensure 2-level group ----
    if (!is.factor(meta[[g_col]])) stop("[forest_prediction_effect] g_col is not a factor.")
    g_levels <- levels(meta[[g_col]])
    a_levels <- unique(meta[[a_col]])

    if (length(g_levels) != 2) stop(sprintf("[forest_prediction_effect] Function currently supports only 2 groups (two levels in g_col)."))
    if (length(a_levels) != 1) stop(sprintf("[forest_prediction_effect] Function currently supports only 1 ancestry (one level in a_col)."))

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
      a_1       = A_1,
      a_2       = A_2,
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

}