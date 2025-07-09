#' Plot observed vs. true interaction effects
#'
#' Visualize how your difference-of-differences estimate captures 
#' true interaction effects in a simulation.
#'
#' True interactions should align on the 45° line; others scatter vertically.
#'
#' @param data Data frame with columns `x_var`, `y_var`, and `color_var`.
#' @param x_var Name of column for true interaction effect.
#' @param y_var Name of column for observed interaction effect.
#' @param color_var Name of column indicating true interaction status.
#' @param features Optional vector of feature names to annotate.
#' @param title Plot title.
#' @param x_label X-axis label.
#' @param y_label Y-axis label.
#' @param point_size Point size (default 1).
#'
#' @return A ggplot object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_point geom_abline scale_color_manual
#' @importFrom ggplot2 labs
#' @importFrom ggrepel geom_text_repel
plot_sim_interaction_effect <- function(
  data,
  x_var,
  y_var,
  color_var,
  features = NULL,
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  point_size = 1
) {
  stopifnot(all(c(x_var, y_var, color_var) %in% names(data)))

  # Add annotate column if features are given
  data$annotate <- if (!is.null(features)) {
    ifelse(data$feature %in% features, data$feature, NA)
  } else {
    NA
  }

  p <- ggplot(
      data = data,
      mapping = aes(
        x = .data[[x_var]],
        y = .data[[y_var]],
        color = factor(.data[[color_var]])
      )
    ) +
    geom_point(
      size = point_size
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      linewidth = 0.3,
      color = "grey50"
    ) +
    scale_color_manual(
      values = c("1" = "#1f77b4", "0" = "grey80"),
      labels = c("1" = "Interaction", "0" = "No interaction")
    ) +
    labs(
      title = title,
      x = ifelse(is.null(x_label), "True interaction effect (log2FC)", x_label),
      y = ifelse(is.null(y_label), "Observed interaction effect (log2FC)", y_label),
      color = "Intended interaction"
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend()

  # Add labels if features given
  if (!is.null(features)) {
    p <- p + ggrepel::geom_text_repel(
      data = subset(data, !is.na(annotate)),
      aes(
        x = .data[[x_var]],
        y = .data[[y_var]],
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

  return(p)
}
