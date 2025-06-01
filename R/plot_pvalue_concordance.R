#' Plot Concordance of P-values Between Two Methods
#'
#' Generates a scatter plot comparing -log10(p-values) from two statistical methods,
#' highlighting concordance and significance patterns. Points are optionally capped
#' to a maximum -log10(p) value for visual clarity, and a significance threshold is
#' visualized with intersecting dashed lines. Faceting can be used to compare across iterations.
#'
#' @param data A `data.frame` containing raw p-values and grouping information.
#' @param p_col_x Character. Column name of p-values from method X (e.g., "p_limma").
#' @param p_col_y Character. Column name of p-values from method Y (e.g., "p_bootstrap").
#' @param signif_source Character. Column name indicating the source of significance (e.g., "Both", "Limma only").
#' @param x_label Character. Label for the x-axis (usually the name of method X).
#' @param y_label Character. Label for the y-axis (usually the name of method Y).
#' @param facet_col Character (optional). Column name used to facet the plot (e.g., "iteration").
#' @param facet_levels Character vector (optional). Specific levels of `facet_col` to include in the plot.
#' @param log_cap Numeric. Maximum -log10(p) value to display (default = 5); values beyond this are capped.
#' @param point_size Numeric. Size of points in the scatter plot (default = 0.5).
#' @param title Character (optional). Title of the plot.
#'
#' @return A `ggplot2` object representing the p-value concordance plot.
#'
#' @keywords internal
#' @noRd
plot_pvalue_concordance <- function(
  data,
  p_col_x,
  p_col_y,
  signif_source,
  x_label,
  y_label,
  facet_col,
  facet_levels,
  log_cap = 5,
  point_size = 0.5,
  title = NULL
) {

  df <- data

  # Filter by specified facet levels if faceting is used
  if (!missing(facet_col)) {
    df <- df[df[[facet_col]] %in% facet_levels, , drop = FALSE]
    df[[facet_col]] <- factor(df[[facet_col]], levels = facet_levels)
  }

  # Compute log-transformed p-values
  df$logp_x <- -log10(df[[p_col_x]])
  df$logp_y <- -log10(df[[p_col_y]])

  # Cap extremes
  df$logp_x_capped <- pmin(df$logp_x, log_cap)
  df$logp_y_capped <- pmin(df$logp_y, log_cap)
  df$capped <- factor(ifelse(df$logp_x > log_cap | df$logp_y > log_cap, "Capped", "Uncapped"))

  # Plot
  p <- ggplot(
    df,
    aes(
      x = logp_x_capped,
      y = logp_y_capped,
      color = !!sym(signif_source),
      shape = capped
      )
    ) +
    geom_point() +
    geom_abline(
      slope = 1, 
      intercept = 0, 
      linetype = "solid", 
      color = "blue",
      linewidth = 0.3
    ) +
    scale_color_manual(
      values = c(
        "None"           = "gray50",
        "Limma only"     = "#1f77b4",
        "Bootstrap only" = "#ff7f0e",
        "Both"           = "purple"
      ),
      name = "Significances source"
    ) +
    scale_shape_manual(
      values = c("Uncapped" = 1, "Capped" = 25),
      name = "Point status"
    ) +
    labs(
      title = title,
      x = x_label,
      y = y_label
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend() +
    theme(legend.position = "right")

  # Apply facetting if specified
  if (!missing(facet_col)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_col)))
  }

  return(p)
}
