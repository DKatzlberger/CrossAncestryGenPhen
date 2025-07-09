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
#' @param facet_col Optional. String. Name of a column to use for faceting the plot (e.g., "feature", "cluster"). If NULL, no faceting is applied.
#' @param facet_levels Optional. A character vector of levels to include in faceting. If NULL and \code{facet_col} is provided, the first 9 unique values are used.
#' @param title Optional. Character string to set the plot title.
#' @param x_label Label for the x-axis. Defaults to the value of \code{x_var}.
#' @param y_label Label for the y-axis. Defaults to \code{"Count"}.
#' @param bins Integer. Number of bins to use for the x-axis histogram (default is 50).
#' @param n_fill_bins Integer. Number of bins to use for binning the fill variable (must be odd; default is 9). Bins are centered at 0 and symmetric. The outermost bins are labeled with \code{<=} and \code{>=}.
#'
#' @return A \code{ggplot2} object showing a histogram of \code{x_var}, with bars stacked and colored by binned values of \code{fill_var}. Faceting is applied if specified.
#'
#' @examples
#' # Basic histogram
#' plot_pvalue_distribution(data = agg_perm, x_var = "p_value", fill_var = "T_obs")
#'
#' # Faceted by cluster
#' plot_pvalue_distribution(
#'   data = agg_perm,
#'   x_var = "p_value",
#'   fill_var = "T_obs",
#'   facet_col = "cluster"
#' )
#'
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
  bins = 50,
  n_fill_bins = 9
) {

  # Convert to dataframe
  df <- as.data.frame(data)

  # Validate x_var
  if (!(x_var %in% names(df))) {
    stop("Input data must contain column: ", x_var)
  }

  # Faceting
  do_facet <- !is.null(facet_col)
  if (do_facet) {
    if (!(facet_col %in% names(df))) {
      stop(sprintf("Column '%s' specified as facet_col not found in data.", facet_col))
    }
    if (is.null(facet_levels)) {
      facet_levels <- head(unique(df[[facet_col]]), 9)
    }
    df <- df[df[[facet_col]] %in% facet_levels, , drop = FALSE]
    df[[facet_col]] <- factor(df[[facet_col]], levels = facet_levels)
  }

  # Fill logic
  if (!is.null(fill_var)) {
    if (!(fill_var %in% names(df))) {
      stop("Input data must contain column: ", fill_var)
    }

    t_vals <- df[[fill_var]]
    max_abs <- max(abs(t_vals), na.rm = TRUE)

    pos_limit <- floor(max_abs)
    neg_limit <- ceiling(-max_abs)

    if (n_fill_bins %% 2 == 0) {
      stop("n_fill_bins must be odd to ensure a center bin at 0.")
    }

    bin_centers <- seq(neg_limit, pos_limit, length.out = n_fill_bins)
    bin_width <- diff(bin_centers)[1]

    bin_breaks <- c(-Inf, head(bin_centers, -1) + bin_width / 2, Inf)
    bin_labels <- format(round(bin_centers, 1), nsmall = 1)
    bin_labels[1] <- paste0("<= ", abs(neg_limit))
    bin_labels[length(bin_labels)] <- paste0(">= ", pos_limit)

    df$fill_bin <- cut(
      t_vals,
      breaks = bin_breaks,
      labels = bin_labels,
      include.lowest = TRUE
    )

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

  # Plot
  p <- ggplot(
    data = df, 
    aes(
      x = .data[[x_var]]
      )
    ) +
    (
    if (!is.null(fill_aes)) 
      geom_histogram(
        fill_aes, 
        bins = bins, 
        position = "stack", 
        color = "black", 
        linewidth = 0.1
      )
     else 
      geom_histogram(
        bins = bins, 
        fill = "grey", 
        color = "black", 
        linewidth = 0.1
      )
    ) 

  # Final styling
  p <- p +
    labs(
      title = title,
      x = ifelse(is.null(x_label), x_var, x_label),
      y = ifelse(is.null(y_label), "count", y_label)
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend()


  # Facetting
  if (!is.null(fill_var)) {
    p <- p + scale_fill_manual(
      values = colors, 
      name = fill_var
    )
  }

  if (do_facet) {
    p <- p + facet_wrap(as.formula(paste("~", facet_col)))
  }

  return(p)
}