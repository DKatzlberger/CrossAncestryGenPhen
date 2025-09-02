#' Plot Estimated Gene Dispersions Between Two Datasets
#'
#' This function creates a scatter plot comparing gene-wise dispersion estimates
#' between two datasets (`estimates_X` and `estimates_Y`) for a specified dispersion method.
#' A diagonal reference line (slope = 1) is included. Points can optionally be colored
#' by the number of random outliers introduced per gene (continuous scale).
#'
#' @param estimates_X A list returned by `estimate_params()` for dataset X. Must contain a `disps` sublist.
#' @param estimates_Y A list returned by `estimate_params()` for dataset Y. Must contain a `disps` sublist.
#' @param method Character. Dispersion type to compare. One of `"mle"` or `"map"`.
#' @param log Logical. If `TRUE` (default), applies `log10` transformation to dispersion values.
#' @param ancestry_X Label for dataset X (used in axis titles).
#' @param ancestry_Y Label for dataset Y (used in axis titles).
#' @param title Plot title.
#' @param x_label X-axis label. If `NULL`, auto-generated.
#' @param y_label Y-axis label. If `NULL`, auto-generated.
#' @param point_size Size of points. Default is 1.
#'
#' @return A ggplot2 scatter plot object.
#' @export
plot_estimated_dispersions <- function(
  estimates_X,
  estimates_Y,
  method = c("mle", "map"),
  log = TRUE,
  ancestry_X = "Dataset X", 
  ancestry_Y = "Dataset Y",
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  point_size = 1
) {
  method <- match.arg(method)

  # Extract dispersion vectors
  x_vals <- estimates_X$disps[[method]]
  y_vals <- estimates_Y$disps[[method]]

  if (length(x_vals) != length(y_vals)) {
    stop("x_disps and y_disps must have the same number of genes.")
  }

  # Optionally log-transform
  if (log) {
    x_vals <- log2(x_vals)
    y_vals <- log2(y_vals)
  }

  # Assemble data frame
  df <- data.frame(
    gene = seq_along(x_vals),
    x = x_vals,
    y = y_vals
  )

  # Build plot
  p <- ggplot(
    data = df, 
    mapping = aes(
      x = x, 
      y = y
      )
    ) +
    geom_point(
      size = point_size
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      color = "gray50",
      linewidth = 0.3
    )

  # Labels and theme
  p <- p + labs(
    title = title,
    x = ifelse(
      is.null(x_label),
      paste0(if (log) "Log2 " else "", toupper(method), " dispersions (", ancestry_X, ")"),
      x_label
    ),
    y = ifelse(
      is.null(y_label),
      paste0(if (log) "Log2 " else "", toupper(method), " dispersions (", ancestry_Y, ")"),
      y_label
    )
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_small_legend()

  return(p)
}
