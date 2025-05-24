#' Plot Empirical Null Distributions for Selected Features
#'
#' Visualizes the permutation null distribution for selected features,
#' showing how empirical p-values were derived from permuted test statistics.
#'
#' @param result A result list from \code{perm_diff_interaction()}.
#' @param features Character vector of feature (gene) names to plot.
#'
#' @return A \code{ggplot} object with histograms and observed statistic overlay.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20)
#' Y <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20)
#' colnames(X) <- colnames(Y) <- paste0("Gene_", seq_len(ncol(X)))
#' MX <- data.frame(condition = factor(rep(c("A", "B"), each = 50)), ancestry = "EUR")
#' MY <- data.frame(condition = factor(rep(c("A", "B"), each = 50)), ancestry = "AFR")
#' result <- perm_diff_interaction(X, Y, MX, MY, "condition", "ancestry", B = 100)
#' p <- plot_perm_distribution(result, features = c("Gene_1", "Gene_2"))
#' ggsave("perm_plot.png", p, width = 8, height = 4)
#'
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline facet_wrap labs theme_minimal
#' @export

plot_perm_distribution <- function(
  result,
  features
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
  p <- ggplot(boot_long, aes(x = T_stat)) +
    geom_histogram(
      bins = 50, 
      fill = "gray80", 
      color = "black"
      ) +
    geom_vline(
      data = obs_df, 
      aes(xintercept = T_obs), 
      color = "red", 
      linewidth = 0.8
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
      base_size = 12
      )

  return(p)
}
