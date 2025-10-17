#' Plot P-value Distribution Colored by a Binned Fill Variable
#'
#' Visualizes the distribution of a variable (typically p-values) using a stacked histogram.
#' Bars are colored based on binned values of a second variable (e.g., a test statistic), centered around zero.
#' Out-of-range values are grouped into outermost bins labeled with `<=` and `>=`.
#' The plot can be optionally faceted by a categorical column.
#'
#' @param data A data frame containing at least the columns specified in \code{x_var} and \code{fill_var}.
#' @param x_var String. Name of the column to use for the x-axis (e.g., "p_value").
#' @param fill_var String. Name of the column used for fill coloring (e.g., "T_obs").
#' @param facet_col Optional. String. Name of a column to use for faceting the plot (e.g., "feature", "cluster").
#' @param facet_levels Optional. A character vector of levels to include in faceting.
#' @param title Optional Plot title.
#' @param x_label Optional x-axis label.
#' @param y_label Optional y-axis label.
#' @param bins Number of bins for the x-axis histogram (default 50).
#' @param fill_bins Number of bins for the fill variable (must be odd; default 9).
#' @param fill_limits Optional numeric vector of length 2 giving the global limits for the fill variable.
#'
#' @return A ggplot2 object.
#' @import ggplot2
#' @export
plot_pvalue_distribution <- function(
  data,
  x_var,
  fill_var = NULL,
  facet_col = NULL,
  facet_levels = NULL,
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  fill_label = NULL,
  bins = 50,
  fill_bins = 9,
  fill_limits = NULL   
) {

  ## --- Convert to dataframe ---
  df <- as.data.frame(data)

  ## --- Validate x_var ---
  if (!(x_var %in% names(df))) {
    stop("Data must contain column: ", x_var)
  }

  ## --- Faceting ---
  do_facet <- !is.null(facet_col)
  if (do_facet) {
    if (!(facet_col %in% names(df))) {
      stop("Data must contain column: ", facet_col)
    }
    if (is.null(facet_levels)) {
      facet_levels <- head(unique(df[[facet_col]]), 9)
    }
    df <- df[df[[facet_col]] %in% facet_levels, , drop = FALSE]
    df[[facet_col]] <- factor(df[[facet_col]], levels = facet_levels)
  }

  ## --- Fill logic ---
  if (!is.null(fill_var)) {
    if (!(fill_var %in% names(df))) {
      stop("Data must contain column: ", fill_var)
    }

    t_vals <- df[[fill_var]]

    # --- Determine range for fill scaling ---
    if (is.null(fill_limits)) {
      max_abs <- max(abs(t_vals), na.rm = TRUE)
      if (max_abs < 1e-6) max_abs <- 1
      neg_limit <- -max_abs
      pos_limit <- max_abs
    } else {
      neg_limit <- fill_limits[1]
      pos_limit <- fill_limits[2]
    }

    # Ensure fill_bins is odd
    if (fill_bins %% 2 == 0) {
      stop("fill_bins must be odd to ensure a center bin at 0.")
    }

    # Define bin centers & breaks
    bin_centers <- seq(neg_limit, pos_limit, length.out = fill_bins)
    bin_width <- diff(bin_centers)[1]
    bin_breaks <- c(-Inf, head(bin_centers, -1) + bin_width / 2, Inf)
    bin_labels <- format(round(bin_centers, 1), nsmall = 1)
    bin_labels[1] <- paste0("<= ", round(neg_limit, 1))
    bin_labels[length(bin_labels)] <- paste0(">= ", round(pos_limit, 1))

    df$fill_bin <- cut(
      t_vals,
      breaks = bin_breaks,
      labels = bin_labels,
      include.lowest = TRUE
    )

    # Number of bins below/above zero for color gradients
    zero_idx <- which.min(abs(bin_centers))
    n_below <- zero_idx - 1
    n_above <- length(bin_centers) - zero_idx

    colors <- c(
      colorRampPalette(c("blue", "white"))(n_below + 1),
      colorRampPalette(c("white", "red"))(n_above + 1)[-1]
    )
    names(colors) <- bin_labels

    fill_aes <- aes(fill = fill_bin)
  } else {
    fill_aes <- NULL
  }

  ## --- Plot ---
  p <- ggplot(
      data = df, 
      mapping = aes(
        x = .data[[x_var]]
      )
    ) +
    (
      if (!is.null(fill_aes)) {
        geom_histogram(
          fill_aes,
          bins = bins,
          position = "stack",
          color = "black",
          linewidth = 0.1
        )
      } else {
        geom_histogram(
          bins = bins,
          fill = "grey80",
          color = "black",
          linewidth = 0.1
        )
      }
    )

  # Fill scale
  if (!is.null(fill_var)) {
    p <- p + scale_fill_manual(values = colors)
  }

  # Facets
  if (do_facet) {
    p <- p + facet_wrap(as.formula(paste("~", facet_col)))
  }

  # Labels and theme
  p <- p +
    labs(
      title = title,
      x = ifelse(is.null(x_label), x_var, x_label),
      y = ifelse(is.null(y_label), "Count", y_label),
      fill = fill_label
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend()

  return(p)
}
