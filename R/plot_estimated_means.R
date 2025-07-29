#' Plot Estimated Gene Means Between Two Datasets
#'
#' This function generates a scatter plot comparing gene-wise expression means
#' between two datasets (`estimates_X` and `estimates_Y`) for a specified normalization method.
#' A reference line with slope = 1 is included to help visualize agreement between datasets.
#'
#' @param estimates_X A list returned by `estimate_params()` for dataset X. Must contain a `means` sublist.
#' @param estimates_Y A list returned by `estimate_params()` for dataset Y. Must contain a `means` sublist.
#' @param method Character. Type of gene-level mean to compare. One of `"raw"` (unadjusted counts),
#'   `"libnorm"` (library-normalized counts), or `"logcpm"` (log2 counts per million).
#' @param ancestry_X Character. Label used for dataset X in axis titles. Default is `"Dataset X"`.
#' @param ancestry_Y Character. Label used for dataset Y in axis titles. Default is `"Dataset Y"`.
#' @param title Plot title.
#' @param x_label X-axis label.
#' @param y_label Y-axis label.
#' @param point_size Numeric. Size of the scatter plot points. Default is `1`.
#'
#' @return A [ggplot2::ggplot()] object showing gene-wise mean comparison.
#'
#' @examples
#' \dontrun{
#' estimates_X <- estimate_params(count_matrix_X)
#' estimates_Y <- estimate_params(count_matrix_Y)
#' plot_estimated_mean(estimates_X, estimates_Y, method = "logcpm")
#' }
#'
#' @export
plot_estimated_mean <- function(
  estimates_X,
  estimates_Y,
  method = c("raw", "libnorm", "logcpm"),
  ancestry_X = "Dataset X", 
  ancestry_Y = "Dataset Y",
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  point_size = 1
){

  method <- match.arg(method)

  # Extract selected mean type
  x_vals <- estimates_X$means[[method]]
  y_vals <- estimates_Y$means[[method]]

  if (length(x_vals) != length(y_vals)) {
    stop("x_means and y_means must have the same number of genes.")
  }

  df <- data.frame(
    x = x_vals,
    y = y_vals
  )

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
    color = "gray50"
  )

  # Final style
  p <- p + labs(
    title = title,
    x = ifelse(is.null(x_label), paste0(ancestry_X, " means ", "(", toupper(method), ")"), x_label),
    y = ifelse(is.null(y_label), paste0(ancestry_Y, " means ", "(", toupper(method), ")"), y_label)
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_small_legend()

  return(p)
}