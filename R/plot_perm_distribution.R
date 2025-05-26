#' Plot Empirical Null Distributions for Selected Features
#'
#' Visualizes the permutation null distribution for selected features,
#' showing how empirical p-values were derived from permuted test statistics.
#'
#' @param result A result list from \code{perm_diff_interaction()}.
#' @param features Character vector of feature (gene) names to plot.
#' @param point_size Numeric value to control fontsize.
#'
#' @return A \code{ggplot} object with histograms and observed statistic overlay.
#'
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline facet_wrap labs theme_minimal
#' @export

plot_perm_distribution <- function(
  result,
  features,
  point_size = 0.5
) {
  T_boot <- result$T_boot
  stats <- result$summary_stats

  features <- intersect(features, colnames(T_boot))
  if (length(features) == 0) stop("None of the selected features are present in the result.")

  # Create long-format data using base R
  boot_long <- data.frame(
    T_stat  = as.vector(T_boot[, features, drop = FALSE]),
    feature = rep(features, each = nrow(T_boot))
  )

  obs_idx <- match(features, stats$feature)
  obs_df <- data.frame(
    feature = features,
    T_obs   = stats$T_obs[obs_idx]
  )

  # Create and return ggplot object
  p <- ggplot(
    data = boot_long, 
    aes(
      x = T_stat
      )
      ) +
    geom_histogram(
      bins = 50, 
      fill = "gray80", 
      color = "black"
      ) +
    geom_vline(
      data = obs_df, 
      aes(
        xintercept = T_obs
        ), 
      color = "blue", 
      linewidth = 0.5,
      linetype = "dashed"
      ) +
    facet_wrap(
      ~feature
      ) +
    labs(
      title = "Permutation Null Distribution of Test Statistics",
      x = "T-statistic",
      y = "Frequency"
      ) +
    theme_minimal(
      base_size = point_size * 10
      )

  return(p)
}
