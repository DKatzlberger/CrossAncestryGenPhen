#' Plot P-value Distribution Colored by Continuous Mean log2FC (Diverging)
#'
#' Shows empirical p-value distributions across repeated permutation iterations.
#' Bars are filled using a diverging color scale: blue (low logFC), white (0), red (high logFC).
#'
#' @param x Output from `repeated_perm_diff_interaction()`.
#' @param features Character vector of gene names to include.
#' @param threshold Optional numeric. Adds a vertical dashed line at a p-value threshold.
#'   Use e.g., 0.05 for significance threshold, or set to NULL to skip it (default = NULL).
#' @param point_size Numeric. Controls base font size (default = 0.5).
#' 
#' @return A ggplot2 object.
#'
#' @importFrom ggplot2 ggplot aes geom_histogram facet_wrap geom_vline scale_fill_gradient2
#' @importFrom ggplot2 labs theme_minimal theme element_text
#' @export
plot_pvalue_distribution <- function(
  x,
  features,
  threshold = NULL,
  point_size = 0.5
) {
  df <- x$all_iterations
  agg <- x$aggregated

  # Validate and filter features
  features <- intersect(features, df$feature)
  if (length(features) == 0) stop("None of the selected features are available.")

  df <- df[df$feature %in% features, , drop = FALSE]
  agg <- agg[agg$feature %in% features, c("feature", "mean_T_obs")]

  # Merge logFC into per-iteration results
  df <- merge(df, agg, by = "feature", all.x = TRUE)

  # Start plot
  p <- ggplot(
    data = df, 
    aes(
      x = p_value, 
      fill = mean_T_obs
      )
    ) +
    geom_histogram(
      bins = 50,
      color = "black"
    ) +
    facet_wrap(
        ~feature, 
        scales = "free_y"
    ) +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      name = "Mean log2FC"
    ) +
    labs(
      title = "P-value Distribution Colored by Mean Interaction Effect",
      x = "Empirical P-value",
      y = "Count"
    ) +
    theme_nature_fonts() +
    theme_small_legend() +
    theme_white_background() 


  # Optional vertical threshold line
  if (!is.null(threshold)) {
    p <- p + geom_vline(
      xintercept = threshold,
      color = "blue",
      linewidth = 0.5,
      linetype = "dashed"
    )
  }

  return(p)
}
