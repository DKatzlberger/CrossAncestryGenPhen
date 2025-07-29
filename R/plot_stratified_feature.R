#' Plot stratified feature distributions across splits
#'
#' Generates boxplots of scaled feature values across data splits
#' (X, Y, R) stratified by a grouping column, optionally for
#' specified features.
#'
#' @param X A numeric matrix or data frame for the test split
#' @param Y A numeric matrix or data frame for the inference split
#' @param R A numeric matrix or data frame for the train split
#' @param MX Metadata for X, containing IDs and group information
#' @param MY Metadata for Y, containing IDs and group information
#' @param MR Metadata for R, containing IDs and group information
#' @param features Character vector of feature names to plot.
#'   Defaults to first 9 common features if NULL.
#' @param g_col Column name in metadata indicating the grouping
#' @param title Optional title for the plot
#' @param point_size Numeric value controlling point/label size (currently not used in plotting directly).
#'
#' @return A ggplot2 object containing the faceted boxplot
#' @export
plot_stratified_feature <- function(
  X, 
  Y, 
  R, 
  MX, 
  MY, 
  MR,
  fill_var, 
  features = NULL,
  title = NULL, 
  x_label = NULL,
  y_label = NULL,
  point_size = 0.5
) {

  # Infer features if not provided
  if (is.null(features)) {
    all_feats <- Reduce(intersect, list(colnames(X), colnames(Y), colnames(R)))
    if (length(all_feats) < 1) {
      stop("No common features found across X, Y, R.")
    }
    features <- sort(all_feats)
    features <- head(features, 9)
  }

  # Subset matrices to selected features
  X <- X[, features, drop = FALSE]
  Y <- Y[, features, drop = FALSE]
  R <- R[, features, drop = FALSE]

  # Combine and scale all rows together
  all_mat <- rbind(R, X, Y)
  all_scaled <- scale(all_mat)

  # Assign split labels
  split_labels <- c(rep("R", nrow(R)), rep("X", nrow(X)), rep("Y", nrow(Y)))

  all_conditions <- c(MR[[fill_var]], MX[[fill_var]], MY[[fill_var]])

  df_all <- data.frame(
    feature = rep(colnames(all_scaled), each = nrow(all_scaled)),
    value = as.vector(all_scaled),
    condition = rep(all_conditions, ncol(all_scaled)),
    split = rep(split_labels, ncol(all_scaled)),
    stringsAsFactors = FALSE
  )

  df_all$feature <- factor(df_all$feature, levels = features)
  df_all$split <- factor(df_all$split, levels = c("R", "X", "Y"))

  p <- ggplot(
    data = df_all, 
    aes(
      x = split, 
      y = value, 
      fill = condition
      )
    ) +
    geom_boxplot(
      outlier.size = 0.25, 
      width = 0.7, 
      position = position_dodge2(
        preserve = "single"
      )
    ) +
    facet_wrap(
      ~feature, 
      scales = "free"
    ) +
    labs(
      title = title,
      x = ifelse(is.null(x_label), x_var, x_label),
      y = ifelse(is.null(y_label), "z-score", y_label),,
      fill = fill_var
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend()

  return(p)
}
