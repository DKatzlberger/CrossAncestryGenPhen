#' Plot stratified feature distributions across splits
#'
#' Generates boxplots of scaled feature values across data splits
#' (X, Y, R) stratified by a grouping column, optionally for
#' specified features.
#'
#' @param X Numeric matrix or data frame for the test split.
#' @param Y Numeric matrix or data frame for the inference split.
#' @param R Numeric matrix or data frame for the train split.
#' @param MX Metadata for X with IDs and group information.
#' @param MY Metadata for Y with IDs and group information.
#' @param MR Metadata for R with IDs and group information.
#' @param fill_var Metadata column used for fill grouping.
#' @param features Features to plot; uses first nine common ones if NULL.
#' @param id_col Column with sample identifiers. Uses rownames when NULL.
#' @param title Optional plot title.
#' @param x_label Label for the x-axis.
#' @param y_label Label for the y-axis.
#' @param point_size Numeric value controlling point size.
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
  id_col = NULL,
  title = NULL, 
  x_label = NULL,
  y_label = NULL,
  point_size = 0.5
) {
  
  if (is.null(rownames(X)) || is.null(rownames(Y)) || is.null(rownames(R))) {
    stop("All input matrices (X, Y, R) must have rownames corresponding to sample IDs.")
  }

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

  # Combine metadata and ensure proper ordering
  if (!is.null(id_col)) {
    if (!(id_col %in% colnames(MR)) || !(id_col %in% colnames(MX)) || !(id_col %in% colnames(MY))) {
      stop(sprintf("id_col '%s' not found in one or more metadata sets.", id_col))
    }
    all_ids <- c(MR[[id_col]], MX[[id_col]], MY[[id_col]])
  } else {
    all_ids <- c(rownames(R), rownames(X), rownames(Y))
  }

  all_conditions <- c(MR[[fill_var]], MX[[fill_var]], MY[[fill_var]])

  df_all <- data.frame(
    sample_id = rep(all_ids, ncol(all_scaled)),
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
