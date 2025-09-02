#' Plot Simulated Interaction Effects
#'
#' Create a histogram of the interaction effects (\code{true_log2FC}) from a
#' simulated interaction effects table \code{fI}.
#'
#' @param fI Data frame containing a column named \code{true_log2FC} for interaction effects.
#' @param exclude_zeros Logical, whether to exclude exact zeros from the histogram.
#' @param title Optional character string for the plot title.
#' @param x_label Optional x-axis label. Defaults to \code{"Interaction effect (log2FC)"}.
#' @param y_label Optional y-axis label. Defaults to \code{"Count"}.
#' @param bins Integer specifying the number of histogram bins. Default is \code{50}.
#'
#' @return A \code{ggplot} object showing the distribution of interaction effects.
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_histogram labs
plot_sim_interaction_effects <- function(
  fI,
  exclude_zeros = TRUE,
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  bins = 50
){

  vals <- fI$true_log2FC
  if (exclude_zeros) {
    vals <- vals[vals != 0]
  }
  
  df <- data.frame(x = vals)
  
   p <- ggplot(
    data = df,
    mapping = aes(
      x = x
    )
  ) +
  geom_histogram(
    bins = bins,
    fill = "gray80",
    color = "black",
    linewidth = 0.1
  )

  # Final styling
  p <- p + labs(
    title = title,
    x = ifelse(is.null(x_label), "Interaction effect (log2FC)", x_label),
    y = ifelse(is.null(y_label), "Count", y_label)
  ) +
  theme_nature_fonts() +
  theme_white_background() +
  theme_small_legend()
  
  return(p)
}
