#' Plot Estimated Expression Outlier Probabilities Between Two Datasets
#'
#' This function compares expression outlier rates (random and single-sample outliers)
#' between two datasets. It uses results from `estimate_params()` and visualizes the
#' estimated probabilities for each type of outlier using grouped bar plots and faceting.
#'
#' @param estimates_X A list returned by `estimate_params()` for dataset X. Must contain `random_outliers` and `single_outliers` sublists.
#' @param estimates_Y A list returned by `estimate_params()` for dataset Y. Must contain `random_outliers` and `single_outliers` sublists.
#' @param ancestry_X Character. Label used for dataset X in the plot. Default is `"Dataset X"`.
#' @param ancestry_Y Character. Label used for dataset Y in the plot. Default is `"Dataset Y"`.
#' @param title Plot title.
#' @param x_label X-axis label.
#' @param y_label Y-axis label.
#' 
#' @return A [ggplot2::ggplot()] object showing bar plots of outlier probabilities across metrics.
#'
#' @examples
#' \dontrun{
#' estimates_X <- estimate_params(count_matrix_X)
#' estimates_Y <- estimate_params(count_matrix_Y)
#' plot_estimated_outliers(estimates_X, estimates_Y)
#' }
#'
#' @export
plot_estimated_outliers <- function(
  estimates_X, 
  estimates_Y, 
  ancestry_X = "Dataset X", 
  ancestry_Y = "Dataset Y",
  title = NULL,
  x_label = NULL,
  y_label = NULL,
) {
  # Extract the key values from each estimate object
  extract_summary <- function(estimates, label) {
    data.frame(
      dataset = label,
      metric = c("Random High", "Random Low", "Single High", "Single Low"),
      probability = c(
        estimates$random_outliers$high$prob,
        estimates$random_outliers$low$prob,
        estimates$single_outliers$high$prob,
        estimates$single_outliers$low$prob
      ),
      stringsAsFactors = FALSE
    )
  }

  # Prepare data.frame
  df1 <- extract_summary(estimates_X, ancestry_X)
  df2 <- extract_summary(estimates_Y, ancestry_Y)
  plot_df <- rbind(df1, df2)

  # Plot with ggplot2
  p <- ggplot(
    data = plot_df, 
    mapping= aes(
      x = metric, 
      y = probability,
      fill = dataset
      )
    ) +
    geom_bar(
      stat = "identity", 
      position = "dodge"
    ) +
    facet_grid(
      cols = vars(metric),
      scales = "free",
      space = "free"
    )

    # Final style
    p <- p + labs(
      title = title,
      x = ifelse(is.null(x_label), "Metric", x_label),
      y = ifelse(is.null(y_label), "Estimated probability", y_label),
      fill = "Ancestry group"
    ) +
    theme_nature_fonts() +
    theme_white_background(show_facets = FALSE) +
    theme_small_legend()

  return(p)
}