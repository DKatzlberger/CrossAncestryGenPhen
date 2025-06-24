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
#' @param g_col Column name in metadata indicating the grouping (e.g., condition)
#' @param id_col Column name in metadata for sample identifiers (assumes rownames of matrix are IDs)
#' @param title Optional title for the plot
#' @param point_size Numeric value controlling point/label size (currently not used in plotting directly)
#'
#' @return A ggplot2 object containing the faceted boxplot
#' @export
plot_feature <- function(
  X, 
  Y, 
  MX, 
  MY, 
  g_col, 
  features = NULL,
  id_col = NULL,
  title = NULL,
  point_size = 0.5
) {
  if (is.null(rownames(X)) || is.null(rownames(Y))) {
      stop("X and Y must have rownames corresponding to sample IDs.")
    }

  # Infer features if not provided
  if (is.null(features)) {
    all_feats <- intersect(colnames(X), colnames(Y))
    if (length(all_feats) < 1) {
      stop("No common features found across X and Y.")
    }
    features <- sort(all_feats)
    features <- head(features, 9)
  }

  # Subset to common features
  X <- X[, features, drop = FALSE]
  Y <- Y[, features, drop = FALSE]

  # Combine matrices and scale together
  combined <- rbind(X, Y)
  combined_scaled <- scale(combined)

  # Assign split labels
  split_labels <- c(rep("X", nrow(X)), rep("Y", nrow(Y)))

  # Combine metadata
  if (!is.null(id_col)) {
    if (!(id_col %in% colnames(MX)) || !(id_col %in% colnames(MY))) {
      stop(sprintf("id_col '%s' not found in metadata.", id_col))
    }
    ids <- c(MX[[id_col]], MY[[id_col]])
  } else {
    ids <- c(rownames(X), rownames(Y))
  }

  conditions <- c(MX[[g_col]], MY[[g_col]])

  df_all <- data.frame(
    sample_id = rep(ids, ncol(combined_scaled)),
    feature = rep(colnames(combined_scaled), each = nrow(combined_scaled)),
    value = as.vector(combined_scaled),
    condition = rep(conditions, ncol(combined_scaled)),
    split = rep(split_labels, ncol(combined_scaled)),
    stringsAsFactors = FALSE
  )

  df_all$feature <- factor(df_all$feature, levels = features)
  df_all$split <- factor(df_all$split, levels = c("X", "Y"))

  p <- ggplot(
    data = df_all, 
    aes(
      x = split, 
      y = value, 
      fill = condition)
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
      x = NULL,
      y = "Z-score",
      fill = g_col,
      title = title
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend()
}
