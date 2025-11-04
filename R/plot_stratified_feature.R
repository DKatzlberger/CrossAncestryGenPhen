#' Plot stratified feature distributions across splits
#'
#' Generates boxplots of scaled feature values across data splits
#' (X, Y, R) stratified by a grouping column, optionally for
#' specified features.
#'
#' @param X Numeric matrix or data.frame of features for cohort X (train).
#' @param Y Numeric matrix or data.frame of features for cohort Y (test).
#' @param R Numeric matrix or data.frame of features for cohort R (inference).
#' @param MX Metadata for `X` (train).
#' @param MY Metadata for `Y` (test).
#' @param MR Metadata for `R` (inference).
#' @param x_var Name of the metadata column for x-axis (string).
#' @param fill_var Name of the metadata column for fill (string).
#' @param features Character vector of features to plot (defaults to first 9).
#' @param title Plot title (optional).
#' @param x_label Label for x-axis (optional).
#' @param y_label Label for y-axis (optional).
#'
#' @return A ggplot2 object.
#' @export
plot_stratified_feature <- function(
  X, 
  Y, 
  R, 
  MX, 
  MY, 
  MR,
  x_var,
  fill_var, 
  features = NULL,
  title = NULL, 
  x_label = NULL,
  y_label = NULL
) {

  ## --- Input data structure check ---
  assert_input(
    X = X,
    Y = Y,
    R = R,
    MX = MX, 
    MY = MY,
    MR = MR,
    g_col = fill_var,
    a_col = x_var
  )

  ## --- Infer features if not provided ---
  if (is.null(features)) {
    features <- colnames(X)[1:9]
  }

  # Subset matrices to selected features
  X <- X[, features, drop = FALSE]
  Y <- Y[, features, drop = FALSE]
  R <- R[, features, drop = FALSE]


  ## --- Scale counts ---
  RXY <- rbind(R, X, Y)
  RXY_scaled <- scale(RXY)

  ## --- Prepare x_var label ---
  MR[[x_var]] <- paste0(MR[[x_var]], "\n(Reference, R)")
  MX[[x_var]] <- paste0(MX[[x_var]], "\n(Subset, X)")
  MY[[x_var]] <- paste0(MY[[x_var]], "\n(Inference, Y)")

  # Ensure factor levels
  M <- rbind(MR, MX, MY)
  M[[x_var]] <- factor(
    M[[x_var]], 
    levels = c(
      unique(MR[[x_var]]),  
      unique(MX[[x_var]]),  
      unique(MY[[x_var]]) 
    )
  )
  stopifnot(identical(rownames(RXY), rownames(M)))

  ## --- Pivot long ---
  df_long <- data.frame(
    feature   = rep(colnames(RXY_scaled), each = nrow(RXY_scaled)),
    value     = as.vector(RXY_scaled),
    condition = rep(M[[fill_var]], times = length(features)),
    x_group   = rep(M[[x_var]],   times = length(features)),
    stringsAsFactors = FALSE
  )
  # Ensure factor levels
  df_long$feature <- factor(df_long$feature, levels = features)

  ## --- Boxplot ---
  p <- ggplot(
    data = df_long, 
    aes(
      x = x_group, 
      y = value, 
      fill = condition
      )
    ) +
    geom_boxplot(
      position = "dodge",
      color = "black",
      linewidth = 0.1,
      outlier.size = 0.25
    ) +
    facet_wrap(
      ~feature, 
      scales = "free_y"
    ) +
    labs(
      title = title,
      x = ifelse(is.null(x_label), x_var, x_label),
      y = ifelse(is.null(y_label), "Z-score", y_label),
      fill = fill_var
    ) +
    theme_CrossAncestryGenPhen()


  ## --- Return ----
  return(p)
}
