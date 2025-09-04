#' Plot feature distributions across two splits
#'
#' Generates boxplots of scaled feature values across two data splits
#' (X and Y), stratified by a grouping column, optionally for specified features.
#'
#' @param X A numeric matrix or data frame for the test split
#' @param Y A numeric matrix or data frame for the inference split
#' @param MX Metadata for X, containing IDs and group information
#' @param MY Metadata for Y, containing IDs and group information
#' @param features Character vector of feature names to plot.
#'   Defaults to first 9 common features if NULL.
#' @param g_col Column name in metadata indicating the group 
#' @param a_col Column name in metadata indicating the ancestry 
#' @param id_col Column name in metadata for sample identifiers (assumes rownames of matrix are IDs)
#' @param title Optional title for the plot
#'
#' @return A ggplot2 object containing the faceted boxplot
#' @export
plot_feature <- function(
  X, 
  Y, 
  MX, 
  MY, 
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
    MX = MX, 
    MY = MY,
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


  ## --- Ancetsry levels ----
  a_1 <- unique(MX[[x_var]]); a_2 <- unique(MY[[x_var]])
  a_levels <- c(a_1, a_2)


  ## --- Bind frames ---
  M <- rbind(MX, MY)
  M[[x_var]] <- factor(M[[x_var]], levels = a_levels)

  ## --- Scale counts ---
  XY <- rbind(X, Y)
  XY_scaled <- scale(XY)
  stopifnot(identical(rownames(XY), rownames(M)))


  ## --- Pivot long ---
  df_long <- data.frame(
    feature   = rep(colnames(XY_scaled), each = nrow(XY_scaled)),
    value     = as.vector(XY_scaled),
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
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend()

  ## --- Return ----
  return(p)
}
