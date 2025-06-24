#' Plot P-value Distribution Colored by a Binned Fill Variable
#'
#' Visualizes p-value distributions for selected features using a stacked histogram.
#' The fill color of each bin segment reflects the binned version of the chosen fill variable.
#'
#' @param data A data frame containing at least the columns \code{feature}, and user-specified \code{x_var} and \code{fill_var}.
#' @param x_var Name of the column to use for the x-axis (e.g., "p_value").
#' @param fill_var Name of the column to use for fill coloring (e.g., "T_obs").
#' @param features Character vector of feature names to include in the plot. If \code{NULL}, the first
#'   9 unique features in \code{feature} are used.
#' @param title Optional character string to set the plot title.
#' @param x_label Label for the x-axis.
#' @param y_label Label for the y-axis. Defaults to \code{"Count"}.
#' @param bins Integer. Number of bins to use for the histogram (default is 50).
#' @param bin_width Numeric. Width used to bin the fill variable (default is 0.5).
#'
#' @return A \code{ggplot2} object showing facetted histograms.
#' 
#' @import ggplot2
#' @export
plot_pvalue_distribution <- function(
  data,
  x_var,
  fill_var = NULL,
  features = NULL,
  title = NULL,
  x_label = x_var,
  y_label = "Count",
  bins = 50,
  bin_width = 0.5
) {
  if (!all(c("feature", x_var, fill_var) %in% names(data))) {
    stop("Input data must contain 'feature', the specified x_var, and fill_var columns.")
  }

  if (is.null(features)) {
    features <- unique(data$feature)
    features <- sort(features)
    features <- features[1:min(9, length(unique(data$feature)))]
  }

  features <- intersect(features, unique(data$feature))
  if (length(features) == 0) stop("None of the selected features are available.")

  df <- data[data$feature %in% features, , drop = FALSE]
  df$feature <- factor(df$feature, levels = features)

  # Ensure bin_levels always includes 0
  t_vals <- df[[fill_var]]
  t_min <- -3
  t_max <- 3

  bin_levels <- seq(t_min, t_max, by = bin_width)
  if (!any(bin_levels == 0)) {
    bin_levels <- sort(c(bin_levels, 0))  # ensure 0 is in the middle
  }
  bin_labels <- format(bin_levels, nsmall = 1)
  bin_labels[1] <- paste0("≤ ", bin_labels[1])
  bin_labels[length(bin_labels)] <- paste0("≥ ", bin_labels[length(bin_labels)])

  # Bin the data
  bin_raw <- floor(t_vals / bin_width) * bin_width
  bin_raw[bin_raw <= t_min] <- t_min
  bin_raw[bin_raw >= t_max] <- t_max
  df$fill_bin <- factor(bin_raw, levels = bin_levels, labels = bin_labels)

  # Color palette: ensure white is at 0
  zero_idx <- which(bin_levels == 0)
  n_below <- zero_idx - 1
  n_above <- length(bin_levels) - zero_idx

  colors <- c(
    colorRampPalette(c("blue", "white"))(n_below + 1),
    colorRampPalette(c("white", "red"))(n_above + 1)[-1]
  )
  names(colors) <- bin_labels


  # Plot
  p <- ggplot(
    df, 
    aes(
      x = .data[[x_var]], 
      fill = fill_bin
      )
    ) +
    geom_histogram(
      bins = bins,
      position = "stack",
      color = "black",
      linewidth = 0.1
    ) +
    facet_wrap(
      ~feature, 
      # scales = "free"
    ) +
    scale_fill_manual(
      values = colors, 
      name = fill_var
    ) +
    labs(
      title = title,
      x = x_label,
      y = y_label
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend()

  return(p)
}
