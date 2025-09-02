#' Plot Mean-Variance Trend Using log2-CPM (voom-style)
#'
#' This function plots the mean-variance relationship across genes based on 
#' log2-counts-per-million (log2-CPM) values derived from raw count data.
#' Each gene is plotted by its average log2-CPM (mean expression) and 
#' standard deviation across samples. A LOWESS curve is overlaid to show 
#' the global trend, similar to the diagnostic plot produced by limma::voom.
#'
#' @param X A numeric matrix or data.frame of **raw count data** (samples in rows, genes in columns).
#' @param title An optional character string specifying the plot title.
#' @param x_label An optional label for the x-axis. Defaults to "Mean log2-CPM".
#' @param y_label An optional label for the y-axis. Defaults to "Standard deviation".
#' @param point_size A numeric value for point size in the scatter plot (default = 1).
#'
#' @return A ggplot2 object showing the mean-standard deviation trend of log2-CPM values.
#'
#' @import ggplot2
#' @export
plot_mean_variance_trend <- function(
  X,
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  point_size = 1
) {

  # Run voom (requires samples in rows, genes in columns)
  v <- limma::voom(t(X), plot = FALSE, save.plot = TRUE)

  df_points <- data.frame(
    x = v$voom.xy$x,
    y = v$voom.xy$y
  )

  df_line <- data.frame(
    x = v$voom.line$x,
    y = v$voom.line$y
  )

  p <- ggplot() +
    geom_point(
      data = df_points, 
      mapping = aes(
        x = x, 
        y = y
      ),
      size = point_size
    ) +
    geom_line(
      data = df_line, 
      mapping = aes(
        x = x, 
        y = y
      ),
      color = "red", 
      linewidth = 0.3
    ) 

    # Final styling
    p <- p + labs(
      title = title,
      x = ifelse(is.null(x_label), "Log2 CPM", x_label),
      y = ifelse(is.null(y_label), "sqrt(Standard deviation)", y_label)
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend() +
    theme(legend.position = "none")

    return(p)
}
