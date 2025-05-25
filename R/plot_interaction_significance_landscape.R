#' Plot Interaction Significance Landscape
#'
#' Creates a scatter plot summarizing the reproducibility and statistical significance
#' of interaction effects across multiple stratified resampling iterations.
#'
#' @param aggregated_df A data frame (e.g., from `repeated_perm_diff_interaction()$aggregated`)
#'   containing at least the feature, effect size, variability, and p-value columns.
#' @param x_var Name of the column to use on the x-axis. Default: "mean_T_obs"
#' @param y_var Name of the column to use on the y-axis. Will be plotted as -log10(y). Default: "median_p"
#' @param color_var Name of the column to map to color. Default: "prob_signif"
#' @param label_var Name of the column used for labeling points. Default: "feature"
#' @param label_threshold Optional numeric. If provided, labels features where
#'   `color_var >= threshold`. If NULL (default), no labels are shown.
#' @param point_size Numeric value to control base font size. Default: 0.5
#'
#' @return A ggplot2 object.
#' @importFrom ggplot2 ggplot aes geom_point geom_text scale_color_viridis_c
#' @importFrom ggplot2 labs theme_minimal theme element_text
#' @importFrom rlang sym !!
#' @importFrom scales trans_new
#' @export
plot_interaction_significance_landscape <- function(
  aggregated_df,
  x_var = "mean_T_obs",
  y_var = "median_p",
  color_var = "prob_signif",
  label_var = "feature",
  label_threshold = NULL,
  point_size = 0.5
) {
  stopifnot(all(c(x_var, y_var, color_var, label_var) %in% colnames(aggregated_df)))

  # Set up plot
  p <- ggplot(
    data = aggregated_df,
    aes(
      x = !!sym(x_var),
      y = -log10(!!sym(y_var)),
      color = !!sym(color_var),
      label = !!sym(label_var)
    )
  ) +
    geom_point() +
    scale_color_viridis_c(
      option = "D",
      direction = 1,
      end = 0.95,
      name = paste("Color:", color_var)
    ) +
    labs(
      title = "Interaction Significance Landscape",
      x = x_var,
      y = paste0("-log10(", y_var, ")")
    ) +
    theme_minimal(base_size = point_size * 10)

  # Conditional labeling
  if (!is.null(label_threshold)) {
    p <- p +
      geom_text(
        aes(label = ifelse(
          !!sym(color_var) >= label_threshold,
          !!sym(label_var),
          ""
        )),
        hjust = 0.5,
        vjust = -1,
        size = 3.5,
        check_overlap = TRUE
      )
  }

  return(p)
}
