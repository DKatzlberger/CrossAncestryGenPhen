#' Plot feature distributions across two splits
#'
#' Generates boxplots of scaled feature values across two data splits
#' (X and Y) stratified by a grouping column.
#'
#' @param X Numeric matrix or data frame for the test split.
#' @param Y Numeric matrix or data frame for the inference split.
#' @param MX Metadata for X with IDs and group information.
#' @param MY Metadata for Y with IDs and group information.
#' @param g_col Metadata column indicating the grouping.
#' @param features Features to plot. If NULL, first nine common features
#'   are used.
#' @param title Optional plot title.
#' @param x_label Label for the x-axis.
#' @param y_label Label for the y-axis.
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
  title = NULL,
  x_label = NULL,
  y_label = NULL
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

  conditions <- c(MX[[g_col]], MY[[g_col]])

  df_all <- data.frame(
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
      x = ifelse(is.null(x_label), "", x_label),
      y = ifelse(is.null(y_label), "z-score", y_label),
      fill = g_col,
      title = title
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend()
}
