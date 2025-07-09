#' Volcano Plot
#'
#' Creates a volcano plot using ggplot2 based on specified x and y variables, 
#' highlighting points that exceed significance and effect size thresholds.
#' Supports p-value capping, customizable point size, and optional fill aesthetics.
#'
#' @param data A data frame containing the variables to be plotted.
#' @param x_var A string specifying the name of the column to be used as the x-axis (e.g., effect size or test statistic).
#' @param y_var A string specifying the name of the column to be used for computing \code{-log10(p)} on the y-axis (typically a p-value or adjusted p-value).
#' @param fill_var Optional string. Name of a column to map to point fill color (used only with shapes that support fill).
#' @param features Features to annotate in the plot.
#' @param sig_thr Optional numeric value specifying the significance threshold for p-values. If \code{NULL}, no significance filtering is applied.
#' @param effect_thr Optional numeric value specifying the minimum absolute effect size threshold. If \code{NULL}, no effect size filtering is applied.
#' @param x_label Optional string for the x-axis label. 
#' @param y_label Optional string for the y-axis label. 
#' @param title Optional string for the plot title.
#' @param log_cap Numeric value specifying the maximum \code{-log10(p)} to display. Values above this are capped. Default is \code{5}.
#' @param epsilon Numeric constant added to p-values before log transformation to avoid \code{log10(0)}. Default is \code{1e-16}.
#'
#' @return A \code{ggplot} object representing the volcano plot.
#'
#' @import ggplot2
#' @export
plot_volcano <- function(
  data, 
  x_var, 
  y_var, 
  fill_var = NULL,
  features = NULL,
  sig_thr = NULL, 
  effect_thr = NULL, 
  title = NULL, 
  x_label = NULL,
  y_label = NULL,
  point_size = 1,
  log_cap = 5,
  epsilon = 1e-15
) {
  
  df <- as.data.frame(data)

  # Extract relevant columns
  x_col <- df[[x_var]]
  y_col <- df[[y_var]]
  
  # Compute -log10(p-values) and cap them
  log_p <- -log10(y_col + epsilon)
  log_p_capped <- pmin(log_p, log_cap)
  capped <- factor(ifelse(log_p > log_cap, "Capped", "Uncapped"), levels = c("Uncapped", "Capped"))

  # Create plotting data frame
  plot_data <- data.frame(
    x = x_col, 
    y = log_p_capped, 
    capped = capped
  )

  # Add fill_var if provided
  if (!is.null(fill_var)) {
    plot_data$fill_var <- df[[fill_var]]
  }

  # Add feature names if available
  if ("feature" %in% colnames(df)) {
    plot_data$feature <- df$feature
  } else {
    plot_data$feature <- rownames(df)
  }

  # Base aesthetics
  aes_mapping <- if (!is.null(fill_var)) {
    aes(
      x = x, 
      y = y, 
      shape = capped, 
      fill = fill_var,
      color = fill_var
    )
  } else {
    aes(
      x = x, 
      y = y, 
      shape = capped
    )
  }

  # Base plot
  p <- ggplot(
      data = plot_data, 
      mapping = aes_mapping
    ) +
    geom_point(
      size = point_size
    ) +  
    scale_shape_manual(
      name = "p-value capping",
      values = c(
        "Uncapped" = 21, 
        "Capped" = 23
      )
    ) 

  if (!is.null(fill_var)) {
    p <- p + scale_fill_viridis_c( 
      name = fill_var,
      option = "viridis"
    ) +
    scale_color_viridis_c(
      name = fill_var,
      option = "viridis"
    )
  }

  # Threshold lines
  if (!is.null(sig_thr)) {
    p <- p + geom_hline(
      yintercept = -log10(
        sig_thr + epsilon
      ), 
      linetype = "dashed", 
      linewidth = 0.3,
      color = "blue"
    )
  }
  if (!is.null(effect_thr)) {
    p <- p + geom_vline(
      xintercept = c(
        -effect_thr, 
        effect_thr
      ), 
      linetype = "dashed", 
      linewidth = 0.3,
      color = "blue"
    )
  }

  # Feature annotation
  if (!is.null(features)) {
    plot_data$annotate <- ifelse(plot_data$feature %in% features, plot_data$feature, NA)

    p <- p + ggrepel::geom_text_repel(
      data = subset(plot_data, !is.na(annotate)),
      aes(
        x = x,
        y = y,
        label = annotate
      ),
      size = 2,
      min.segment.length = 0,
      segment.size = 0.3,
      max.overlaps = Inf,
      force = 2,
      box.padding = 0.5,
      bg.color = "white",
      bg.r = 0.2,
      inherit.aes = FALSE
    )
  }

  # Final styling
  p <- p +
    labs(
      title = title,
      x = ifelse(is.null(x_label), x_var, x_label),
      y = ifelse(is.null(y_label), paste0("-log10(", y_var, ")"), y_label)
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend() +
    theme(legend.position = "right")

  return(p)
}
