#' Plot P-value Distribution Colored by Aggregated T_obs
#'
#' Visualizes empirical p-value distributions for selected features using a diverging color scale.
#' The fill color represents the aggregated T-statistic (`T_obs`) per feature, using a user-defined
#' aggregation function (e.g., mean, median).
#'
#' @param x A data frame containing at least the columns \code{feature}, \code{p_value}, and \code{T_obs}.
#' @param features Character vector of feature names to include in the plot. If \code{NULL}, the first
#'   9 unique features in \code{x} are used.
#' @param aggregation_fun A function used to summarize \code{T_obs} values per feature. Default is \code{mean}.
#'   Other options include \code{median}, or any custom summary function.
#' @param title Optional character string to set the plot title.
#' @param point_size Numeric. Controls base font size (currently not actively used, retained for consistency).
#'
#' @return A \code{ggplot2} object with facetted histograms showing p-value distributions, colored by the
#'   aggregated T-statistic value per feature.
#'
#' @examples
#' # Basic usage with default mean aggregation
#' plot_pvalue_distribution(combined_results)
#'
#' # Using median aggregation
#' plot_pvalue_distribution(combined_results, aggregation_fun = median)
#'
#' @importFrom ggplot2 ggplot aes geom_histogram facet_wrap scale_fill_gradient2
#' @importFrom ggplot2 labs theme
#' @export
plot_pvalue_distribution <- function(
  x,
  features = NULL,
  aggregation_fun = mean,
  title = NULL,
  point_size = 0.5
) {

  required_cols <- c("feature", "p_value", "T_obs")
  if (!all(required_cols %in% colnames(x))) {
    stop("Input data frame must include columns: feature, p_value, T_obs")
  }

  # Default to first 9 features if none specified
  if (is.null(features)) {
    features <- unique(x$feature)[1:min(9, length(unique(x$feature)))]
  }

  # Validate and subset features
  features <- intersect(features, unique(x$feature))
  if (length(features) == 0) stop("None of the selected features are available.")
  df <- x[x$feature %in% features, , drop = FALSE]

  # Compute aggregated T_obs per feature
  agg_name <- deparse(substitute(aggregation_fun))
  mean_df <- aggregate(T_obs ~ feature, data = df, FUN = aggregation_fun)
  colnames(mean_df)[2] <- paste0(agg_name, "_T_obs")

  # Merge back into main df
  df <- merge(df, mean_df, by = "feature")
  colnames(df)[colnames(df) == paste0(agg_name, "_T_obs")] <- "agg_T_obs"
  df$feature <- factor(df$feature, levels = features)

  # Plot
  p <- ggplot(
    data = df, 
    aes(
      x = p_value, 
      fill = agg_T_obs
    )
  ) +
    geom_histogram(
      bins = 50, 
      color = "black"
    ) +
    facet_wrap(
      ~feature, 
      scales = "free_y"
    ) +
    scale_fill_gradient2(
      low = "blue", 
      mid = "white", 
      high = "red", 
      midpoint = 0,
      name = paste0(agg_name, " T_obs")
    ) +
    labs(
      title = title,
      x = "Empirical p-value",
      y = "Count"
    ) +
    theme_nature_fonts() +
    theme_small_legend() +
    theme_white_background()

  return(p)
}
