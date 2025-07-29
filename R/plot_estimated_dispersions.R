#' Plot Estimated Gene Dispersions Between Two Datasets
#'
#' This function creates a scatter plot comparing gene-wise dispersion estimates
#' between two datasets (`estimates_X` and `estimates_Y`) for a specified dispersion method.
#' It is useful for assessing agreement or shifts in variability between conditions or sample groups.
#' A diagonal reference line (slope = 1) is included.
#'
#' @param estimates_X A list returned by `estimate_params()` for dataset X. Must contain a `disps` sublist.
#' @param estimates_Y A list returned by `estimate_params()` for dataset Y. Must contain a `disps` sublist.
#' @param method Character. Dispersion type to compare. One of `"mle"` (maximum likelihood estimates) or `"map"` (posterior mode estimates).
#' @param ancestry_X Character. Label used for dataset X in axis titles. Default is `"Dataset X"`.
#' @param ancestry_Y Character. Label used for dataset Y in axis titles. Default is `"Dataset Y"`.
#' @param title Plot title.
#' @param x_label X-axis label.
#' @param y_label Y-axis label.
#' @param point_size Numeric. Size of the points in the scatter plot. Default is `1`.
#'
#' @return A [ggplot2::ggplot()] object comparing gene-wise dispersion estimates.
#'
#' @examples
#' \dontrun{
#' estimates_X <- estimate_params(count_matrix_X)
#' estimates_Y <- estimate_params(count_matrix_Y)
#' plot_estimated_dispersion(estimates_X, estimates_Y, method = "mle")
#' }
#'
#' @export
plot_estimated_dispersion <- function(
  estimates_X,
  estimates_Y,
  method = c("mle", "map"),
  ancestry_X = "Dataset X", 
  ancestry_Y = "Dataset Y",
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  point_size = 1
){

  method <- match.arg(method)

  # Extract selected mean type
  x_vals <- estimates_X$disps[[method]]
  y_vals <- estimates_Y$disps[[method]]

  if (length(x_vals) != length(y_vals)) {
    stop("x_disps and y_disps must have the same number of genes.")
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
    x = ifelse(is.null(x_label), paste0(ancestry_X, " dispersions ", "(", toupper(method), ")"), x_label),
    y = ifelse(is.null(y_label), paste0(ancestry_Y, " dispersions ", "(", toupper(method), ")"), y_label)
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_small_legend()

  return(p)
}