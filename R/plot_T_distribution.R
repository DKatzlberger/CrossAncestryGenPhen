#' Plot T-Statistic Distributions with Observed Values and Annotations
#'
#' This function visualizes empirical null distributions of T-statistics for selected features (e.g., genes)
#' using either permutation- or bootstrap-based inference. It overlays observed T-statistics as dashed lines,
#' and when bootstrapping is used, adds 95% confidence intervals as dotted red lines. For permutation-based
#' inference, empirical p-values are shown in the facet titles.
#'
#' @param x A named list output from \code{perm_interaction()} or similar, containing null distributions
#'        and summary statistics. Must include a \code{summary_stats} data frame and either \code{T_perm} or
#'        \code{T_boot} depending on the statistic type.
#' @param features Character vector of feature names to plot. If NULL, the first 9 features (columns) from the
#'        selected null matrix are used.
#' @param statistic Character string indicating which null distribution to plot. Must be either \code{"T_perm"}
#'        (for permutation) or \code{"T_boot"} (for bootstrap).
#' @param title Optional title for the entire plot.
#' @param point_size Numeric value to control relative font or point size (not actively used in this version).
#'
#' @return A \code{ggplot} object consisting of histograms of test statistics per feature, with vertical
#'         lines showing observed values and either p-values (for permutations) or confidence intervals
#'         (for bootstraps).
#'
#' @details
#' The function automatically creates facet plots for each feature selected. If \code{statistic = "T_perm"},
#' the observed test statistic is plotted with a dashed blue line, and the empirical p-value is included
#' in the facet label. If \code{statistic = "T_boot"}, the observed value is plotted with a dashed blue line
#' and the 95% CI is shown using dotted red lines.
#'
#' If the required null distribution is missing from the input object, the function will throw an informative error.
#'
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline facet_wrap labs scale_linetype_manual scale_color_manual
#' @export

plot_T_distribution <- function(
  x,
  features = NULL,
  statistic = c("T_perm", "T_boot"),
  title = NULL,
  point_size = 0.5
) {
  statistic <- match.arg(statistic)

  if (!statistic %in% names(x)) {
    msg <- switch(
      statistic,
      "T_perm" = "Can't find permuted null distribution — maybe you forgot to run permutation?",
      "Unknown test statistic."
    )
    stop(msg)
  }

  T_matrix <- x[[statistic]]
  stats <- x$summary_stats

  if (is.null(features)) {
    features <- colnames(T_matrix)[1:min(9, ncol(T_matrix))]
  }

  features <- intersect(features, colnames(T_matrix))
  if (length(features) == 0) stop("None of the selected features are present in the result.")

  long_df <- data.frame(
    T_stat  = as.vector(T_matrix[, features, drop = FALSE]),
    feature = rep(features, each = nrow(T_matrix))
  )

  obs_idx <- match(features, stats$feature)
  obs_df <- data.frame(
    feature = features,
    T_obs   = stats$T_obs[obs_idx],
    linetype = "Observed T",
    color = "Observed T"
  )

  if (statistic == "T_boot") {
    obs_df$CI_lower <- stats$CI_lower[obs_idx]
    obs_df$CI_upper <- stats$CI_upper[obs_idx]
    ci_df <- data.frame(
      feature = rep(features, 2),
      xintercept = c(obs_df$CI_lower, obs_df$CI_upper),
      linetype = "95% CI",
      color = "95% CI"
    )
  }

  if (statistic == "T_perm") {
    pvals <- stats$p_value[obs_idx]
    formatted <- paste0(features, "\n(p = ", signif(pvals, 3), ")")
    long_df$feature <- factor(long_df$feature, levels = features, labels = formatted)
    obs_df$feature  <- factor(obs_df$feature, levels = features, labels = formatted)
    ci_df <- NULL
  }

  p <- ggplot(long_df, aes(x = T_stat)) +
    geom_histogram(
      bins = 50,
      fill = "gray80",
      color = "black"
    ) +
    geom_vline(
      data = obs_df,
      aes(xintercept = T_obs, linetype = linetype, color = color),
      linewidth = 0.6
    )

  if (statistic == "T_boot") {
    p <- p + geom_vline(
      data = ci_df,
      aes(xintercept = xintercept, linetype = linetype, color = color),
      linewidth = 0.6
    )
  }

  p <- p +
    scale_linetype_manual(
      name = "Line type",
      values = c("Observed T" = "dashed", "95% CI" = "dotted")
    ) +
    scale_color_manual(
      name = "Line type",
      values = c("Observed T" = "blue", "95% CI" = "red")
    ) +
    facet_wrap(~feature) +
    labs(
      title = title,
      x = paste0(statistic, " (T-statistic)"),
      y = "Count"
    ) +
    theme_nature_fonts() +
    theme_small_legend() +
    theme_white_background()

  return(p)
}
