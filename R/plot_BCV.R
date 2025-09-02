#' Mean–Dispersion Scatter Plot with Optional Overlay
#'
#' Plots (optionally log-transformed) means vs dispersions for one dataset,
#' with an optional second dataset overlaid. Points from the first dataset
#' are drawn in black; the optional overlay is drawn in red.
#'
#' If `log = TRUE`, values are transformed with `log2()` **before** plotting,
#' so the axes in ggplot are linear (no log scales). Non-finite points
#' created by the transform (e.g., non-positive inputs) are dropped.
#'
#' @param estimates_X A list with components `means` and `disps`. Each of these
#'   should be a list keyed by method (e.g., `list(mle = <num>, map = <num>)`),
#'   where each entry is a numeric vector of the same length.
#' @param estimates_Y Optional overlay in the same structure as `estimates_X`
#'   (i.e., a list with `means` and `disps`, each a list keyed by `method`).
#' @param method Character; which estimate to plot. One of `"mle"` or `"map"`.
#' @param log Logical; if `TRUE`, apply `log2()` to means and dispersions
#'   **before** plotting (recommended if values span orders of magnitude).
#' @param ancestry_X Character label for the primary dataset (used in legend/title/labels).
#' @param ancestry_Y Character label for the overlay dataset (used in legend/title/labels).
#' @param title Optional plot title. If `NULL`, a sensible default is created,
#'   mentioning the ancestry labels and colors (black/red).
#' @param x_label Optional x-axis label. If `NULL`, a dynamic default is constructed:
#'   `"Log2 MLE means (Dataset X)"` (or similar), depending on `log`, `method`,
#'   and whether an overlay is present.
#' @param y_label Optional y-axis label. If `NULL`, a dynamic default is constructed:
#'   `"Log2 MLE dispersions (Dataset X)"` (or similar), depending on `log`, `method`,
#'   and whether an overlay is present.
#' @param point_size Numeric; point size passed to `geom_point()`. Default `1`.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' # Suppose you have lists like:
#' # estimates_X <- list(
#' #   means = list(mle = runif(100, 1, 100), map = runif(100, 1, 100)),
#' #   disps = list(mle = runif(100, 1, 10),  map = runif(100, 1, 10))
#' # )
#' # estimates_Y <- similar structure...
#'
#' # Single dataset (black), pre-log2:
#' plot_BCV(estimates_X, method = "mle", log = TRUE, ancestry_X = "EUR")
#'
#' # With overlay (EUR black vs AFR red), no log transform:
#' plot_BCV(
#'   estimates_X, estimates_Y,
#'   method = "map", log = FALSE,
#'   ancestry_X = "EUR", ancestry_Y = "AFR"
#' )
#' }
#'
#' @import ggplot2
#' @export
plot_BCV <- function(
  estimates_X,
  estimates_Y = NULL,
  method = c("mle", "map"),
  log = TRUE,
  ancestry_X = "Dataset X",
  ancestry_Y = "Dataset Y",
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  point_size = 1
){
  method <- match.arg(method)

  X_means <- estimates_X$means[[method]]
  X_disps <- estimates_X$disps[[method]]

  if (!is.null(estimates_Y)) {
    Y_means <- estimates_Y$means[[method]]
    Y_disps <- estimates_Y$disps[[method]]
  } else {
    Y_means <- NULL
    Y_disps <- NULL
  }

  if (isTRUE(log)) {
    X_means <- log2(X_means)
    X_disps <- log2(X_disps)
    if (!is.null(Y_means)) {
      Y_means <- log2(Y_means)
      Y_disps <- log2(Y_disps)
    }
  }

  df_X <- data.frame(mean = X_means, disp = X_disps, group = ancestry_X)
  df <- df_X

  has_overlay <- !is.null(Y_means) && !is.null(Y_disps)
  if (has_overlay) {
    df_Y <- data.frame(mean = Y_means, disp = Y_disps, group = ancestry_Y)
    df <- rbind(df_X, df_Y)
  }

  # Drop non-finite (e.g., log2 of <= 0)
  df <- subset(df, is.finite(mean) & is.finite(disp))

  if (has_overlay) {
    p <- ggplot(
      data = df, 
      mapping = aes(
        x = mean, 
        y = disp, 
        color = group
        )
      ) +
      geom_point(
        size = point_size
      )

  } else {
    p <- ggplot(
      data = df, 
      mapping = aes(
        x = mean, 
        y = disp
        )
      ) +
      geom_point(
        size = point_size
      )
  }

  meth_up <- toupper(method)
  prefix   <- if (isTRUE(log)) "Log2 " else ""

  default_x <- if (has_overlay) {
    paste0(prefix, meth_up, " means")
  } else {
    paste0(prefix, meth_up, " means (", ancestry_X, ")")
  }

  default_y <- if (has_overlay) {
    paste0(prefix, meth_up, " dispersions")
  } else {
    paste0(prefix, meth_up, " dispersions (", ancestry_X, ")")
  }

  p <- p +
    labs(
      title = title,
      x = ifelse(is.null(x_label), default_x, x_label),
      y = ifelse(is.null(y_label), default_y, y_label),
      color = NULL
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend()

  return(p)
}
