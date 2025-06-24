#' Volcano Plot
#'
#' Creates a volcano plot using ggplot2 based on specified x and y variables, 
#' highlighting points that exceed significance and effect size thresholds.
#'
#' @param data A data frame containing the variables to be plotted.
#' @param x_var A string specifying the name of the column to be used as the x-axis (e.g., effect size or test statistic).
#' @param y_var A string specifying the name of the column to be used for computing \code{-log10(p)} on the y-axis (typically a p-value or adjusted p-value).
#' @param sig_thr Optional numeric value specifying the significance threshold for p-values. If \code{NULL}, no significance filtering is applied.
#' @param effect_thr Optional numeric value specifying the minimum absolute effect size threshold. If \code{NULL}, no effect size filtering is applied.
#' @param x_label Optional string for the x-axis label. If \code{NULL}, uses \code{x_var} as default.
#' @param y_label Optional string for the y-axis label. Default is \code{paste0("-log10(", y_var, ")")}.
#' @param title Optional string for the plot title.
#' @param point_size Numeric value for the size of the points in the scatter plot. Default is \code{0.5}.
#'
#' @return A \code{ggplot} object representing the volcano plot.
#'
#' @examples
#' df <- data.frame(
#'   feature = LETTERS[1:10],
#'   T_obs = rnorm(10),
#'   p_value = runif(10),
#'   p_adj = p.adjust(runif(10))
#' )
#' plot_volcano(df, x_var = "T_obs", y_var = "p_adj", sig_thr = 0.05, effect_thr = 1)
#'
#' @import ggplot2
#' @export

plot_volcano <- function(
  data, 
  x_var, 
  y_var, 
  sig_thr = NULL, 
  effect_thr = NULL, 
  x_label = NULL,
  y_label = paste0("-log10(", y_var, ")"),
  title = NULL, 
  point_size = 0.5
) {
  
  df <- as.data.frame(data)

  # Extract relevant columns
  x_col <- df[[x_var]]
  y_col <- df[[y_var]]
  
  # Compute -log10 p-values
  log_p <- -log10(y_col)
  
  # Define significance color (base R way)
  sig <- rep("Not Significant", length(x_col))
  if (!is.null(sig_thr) && !is.null(effect_thr)) {
    sig[y_col < sig_thr & abs(x_col) > effect_thr] <- "Significant"
  }
  
  # Create a plotting data frame
  plot_data <- data.frame(
    x = x_col, 
    y = log_p, 
    sig = factor(
      sig, 
      levels = c(
        "Not Significant", 
        "Significant"
        )
      )
  )
  
  # Base plot
  p <- ggplot(
    data = plot_data, 
    aes(
      x = x, 
      y = y, 
      color = sig
      )
    ) +
    geom_point() +
    scale_color_manual(
      values = c(
        "Not Significant" = "gray", 
        "Significant" = "red"
        )
    ) 
  
  # Add thresholds if specified
  if (!is.null(sig_thr)) {
    p <- p + geom_hline(
      yintercept = -log10(sig_thr), 
      linetype = "dashed", 
      color = "blue"
      )
  }
  if (!is.null(effect_thr)) {
    p <- p + geom_vline(
      xintercept = c(-effect_thr, effect_thr), 
      linetype = "dashed", 
      color = "blue"
      )
  }

  # Final plot decoration
  p <- p +
    labs(
      title = title,
      x = x_label,
      y = y_label
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend() +
    theme(legend.position = "none")

  return(p)
}
