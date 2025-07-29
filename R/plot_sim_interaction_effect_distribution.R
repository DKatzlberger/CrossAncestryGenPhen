#' Plot Interaction Effect Distribution with Optional Split
#'
#' This function plots the distribution of an interaction effect variable 
#' using `ggplot2`. If `effect_thr` is provided, the distribution is split 
#' into three symmetric panels: 
#' - Below region: values less than `-effect_thr`
#' - Inside region: values between `-effect_thr` and `+effect_thr`
#' - Above region: values greater than `+effect_thr`
#'
#' Panels are combined side-by-side using `{patchwork}`. Each panel shows 
#' symmetric x-axis tick marks and optional axis labels.
#'
#' @param data A `data.frame` or similar object containing the variable to plot.
#' @param x_var A string. The name of the column in `data` to plot on the x-axis.
#' @param effect_thr Numeric. Threshold for splitting the distribution 
#'   (symmetric around zero). If `NULL`, no split is done and a single 
#'   histogram is plotted.
#' @param title Optional plot title.
#' @param x_label Optional x-axis label for the middle panel.
#' @param y_label Optional y-axis label for the left panel.
#' @param fill_label Optional fill legend label (not used by default).
#' @param bins Number of bins for the histogram.
#'
#' @return A `ggplot` object or a `{patchwork}` combined plot if 
#'   `effect_thr` is provided.
#'
#'
#' @import ggplot2
#' @import patchwork
#' 
#' @export
plot_sim_interaction_effect_distribution <- function(
  data,
  x_var,
  effect_thr = NULL,
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  bins = 50
) {
  df <- as.data.frame(data)

  if (is.null(effect_thr)) {
    # No split, regular plot
    p <- ggplot(
      data = df,
      mapping = aes(
        x = .data[[x_var]]
      )
    ) +
      geom_histogram(
        bins = bins,
        fill = "grey80",
        color = "black",
        linewidth = 0.1
      ) +
      labs(
        title = title,
        x = x_label,
        y = y_label
      )
    return(p)
  }

  # Symmetric split
  x_min <- min(df[[x_var]], na.rm = TRUE)
  x_max <- max(df[[x_var]], na.rm = TRUE)
  max_abs <- max(abs(x_min), abs(x_max), effect_thr)

  # Below region: [-max_abs, -effect_thr]
  df_below <- df[df[[x_var]] < -effect_thr, , drop = FALSE]
  # Inside region: [-effect_thr, +effect_thr]
  df_inside <- df[df[[x_var]] >= -effect_thr & df[[x_var]] <= effect_thr, , drop = FALSE]
  # Above region: [+effect_thr, max_abs]
  df_above <- df[df[[x_var]] > effect_thr, , drop = FALSE]

  # Below plot
  p_below <- ggplot(
    data = df_below, 
    mapping = aes(
      x = .data[[x_var]]
      )
    ) +
    geom_histogram(
      bins = bins,
      fill = "grey80",
      color = "black",
      linewidth = 0.1
    ) +
    coord_cartesian(
      xlim = c(-max_abs, -effect_thr)
    ) +
    scale_x_continuous(
      breaks = pretty(c(-max_abs, -effect_thr), n = 2)
    ) +
    labs(
      subtitle = paste0("Effect size: [-∞", ", -", effect_thr,"]"),
      y = y_label,
      x = NULL
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend() +
    theme(
      axis.title.x = element_blank()
    )

  # Inside plot
  p_inside <- ggplot(
    data = df_inside, 
    mapping = aes(
      x = .data[[x_var]]
      )
    ) +
    geom_histogram(
      bins = bins,
      fill = "grey80",
      color = "black",
      linewidth = 0.1
    ) +
    coord_cartesian(
      xlim = c(-effect_thr, effect_thr)
    ) +
    scale_x_continuous(
      breaks = pretty(c(-effect_thr, effect_thr), n = 2)
    ) +
    labs(
      subtitle = paste0("Effect size: [-", effect_thr, ", +", effect_thr,"]"),
      x = x_label,
      y = NULL
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend() +
    theme(
      axis.title.y = element_blank()
    )

  # Above plot
  p_above <- ggplot(
    data = df_above, 
    mapping = aes(
      x = .data[[x_var]]
      )
    ) +
    geom_histogram(
      bins = bins,
      fill = "grey80",
      color = "black",
      linewidth = 0.1
    ) +
    coord_cartesian(
      xlim = c(effect_thr, max_abs)
    ) +
    scale_x_continuous(
      breaks = pretty(c(effect_thr, max_abs), n = 2)
    ) +
    labs(
      subtitle = paste0("Effect size: [+", effect_thr, ", +∞]"),
      x = NULL,
      y = NULL
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend() + 
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )

  # Combine with patchwork
  combined_theme <- theme_nature_fonts() + theme_white_background() + theme_small_legend() 
  combined_plot <- (p_below | p_inside | p_above) +
    plot_annotation(
        title = title,
        theme = combined_theme
        )

  return(combined_plot)
}

