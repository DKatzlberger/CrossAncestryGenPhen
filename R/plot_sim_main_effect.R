#' Plot Simulated Main Effects
#'
#' Create a scatter plot comparing the true log2 fold-change (log2FC) values
#' between two simulated datasets \code{fX} and \code{fY}.  
#' The plot includes a 1:1 reference line to visualize agreement between datasets.
#'
#' @param fX Data frame containing a column named \code{true_log2FC} for dataset X.
#' @param fY Data frame containing a column named \code{true_log2FC} for dataset Y.
#' @param title Optional character string for the plot title.
#' @param x_label Optional x-axis label. Defaults to \code{"Main effect in fX"}.
#' @param y_label Optional y-axis label. Defaults to \code{"Main effect in fY"}.
#' @param point_size Numeric value controlling point size in the plot.
#'
#' @return A \code{ggplot} object showing main effect comparison.
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_abline labs
plot_sim_main_effect <- function(
  fX,
  fY,
  ancestry_X = "Dataset X", 
  ancestry_Y = "Dataset Y",
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  point_size = 1
){

  stopifnot(dim(fX) == dim(fY))

  df <- data.frame(
    x = fX$T_obs, 
    y = fY$T_obs
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
    color = "gray50",
    linewidth = 0.3
  ) 

  # Final styling
  p <- p + labs(
    title = title,
    x = ifelse(is.null(x_label), paste("Main effect in", ancestry_X, "(log2FC)"), x_label),
    y = ifelse(is.null(y_label), paste("Main effect in", ancestry_Y, "(log2FC)"), y_label)
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_small_legend()

  return(p)
}