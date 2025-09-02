#' Plot permutation-based T-statistics with observed values and p-values
#'
#' Plots null distributions of T-statistics from \code{perm_interaction_effect()}.
#' Adds dashed lines for observed T-statistics, optional normal fit overlay, and optional density overlay.
#'
#' @param data A list from \code{perm_interaction_effect()}, must include
#'   \code{T_null} and \code{summary_stats}.
#' @param features Character vector of features to plot. Defaults to first 9 if NULL.
#' @param show_p Logical. Show parametric p-values? (default TRUE)
#' @param show_normal Logical. Overlay fitted normal curve? (default TRUE)
#' @param title Optional plot title.
#' @param x_label Optional string for the x-axis label. 
#' @param y_label Optional string for the y-axis label. 
#' @param bins Integer number of histogram bins (default 50)
#'
#' @return A ggplot2 object.
#' @export
plot_null_distribution <- function(
  data,
  features = NULL,
  show_p = TRUE,
  show_normal = TRUE,
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  bins = 50
) {
  if (!"T_null" %in% names(data)) stop("List must include 'T_null'.")
  if (!"summary_stats" %in% names(data)) stop("List must include 'summary_stats'.")

  T_matrix <- data$T_null
  stats <- data$summary_stats

  if (is.null(features)) {
    features <- colnames(T_matrix)
    features <- sort(features)
    features <- features[1:min(9, ncol(T_matrix))]
  }
  features <- intersect(features, colnames(T_matrix))
  if (length(features) == 0) stop("None of the selected features are found.")

  # Long format
  long_df <- do.call(rbind, lapply(features, function(f) {
    data.frame(T_stat = T_matrix[, f], feature = f, stringsAsFactors = FALSE)
  }))

  obs_idx <- match(features, stats$feature)
  obs_df <- data.frame(feature = features, T_obs = stats$T_obs[obs_idx])

  label_text <- mapply(function(i) {
    vals <- c()
    if (show_p && "p_param_value" %in% names(stats)) {
      vals <- c(vals, paste0("p_param = ", signif(stats$p_param_value[i], 3)))
    }
    paste(vals, collapse = "\n")
  }, i = obs_idx, SIMPLIFY = TRUE)

  pval_df <- data.frame(feature = features, label = label_text)

  # Determine full x-range across all features
  T_range <- range(unlist(T_matrix[, features]), na.rm = TRUE)
  x_buffer <- diff(T_range) * 0.05  # 5% padding
  x_min <- T_range[1] - x_buffer
  x_max <- T_range[2] + x_buffer
  x_vals_full <- seq(x_min, x_max, length.out = 300)

  # Normal curves (scaled to histogram bin height)
  norm_curves <- do.call(rbind, lapply(features, function(f) {
    vals <- T_matrix[, f]
    mu <- mean(vals)
    sd <- sd(vals)
    hist_data <- hist(vals, plot = FALSE, breaks = bins)
    bin_width <- hist_data$breaks[2] - hist_data$breaks[1]
    y_vals <- dnorm(x_vals_full, mean = mu, sd = sd) * length(vals) * bin_width
    data.frame(x = x_vals_full, y = y_vals, feature = f)
  }))

  # Plot
  p <- ggplot(
    data = long_df, 
    aes(
      x = T_stat
      )
    ) +
    geom_histogram(
      bins = bins,
      fill = "gray80",
      color = "black",
      linewidth = 0.1
    ) +
    geom_vline(
      data = obs_df,
      aes(
        xintercept = T_obs, 
        color = "Obs. T-stat"
      ),
      linetype = "dashed",
      linewidth = 0.5
    ) +
    geom_text(
      data = pval_df,
      aes(
        x = -Inf, 
        y = Inf, 
        label = label
      ),
      hjust = -0.1, 
      vjust = 1.1,
      size = 2,
      inherit.aes = FALSE
    )

  if (show_normal) {
    p <- p + geom_line(
      data = norm_curves,
      aes(
        x = x, 
        y = y, 
        color = "Std. normal dist."
      ),
      linewidth = 0.5,
      inherit.aes = FALSE
    )
  }

  p <- p +
    scale_color_manual(
      name = NULL,
      values = c(
        "Obs. T-stat" = "blue",
        "Std. normal dist." = "red"
      ),
      breaks = c(
        "Obs. T-stat",
        "Std. normal dist."
      )
    ) +
    labs(
      title = title,
      x = ifelse(is.null(x_label), "T_obs (under null)", x_label),
      y = ifelse(is.null(y_label), "Count", y_label),
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend()

  if (length(features) > 1) {
    p <- p + facet_wrap(
      ~feature,
      )
  }

  return(p)
}
