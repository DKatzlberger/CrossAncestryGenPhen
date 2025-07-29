#' Plot Estimated Library Sizes Between Two Datasets
#'
#' This function creates a violin and boxplot comparison of estimated library sizes
#' between two datasets, as returned by `estimate_params()`. This can be used to
#' check for sequencing depth differences or normalization bias between groups.
#'
#' @param estimates_X A list returned by `estimate_params()` for dataset X. Must contain a `libsizes` sublist.
#' @param estimates_Y A list returned by `estimate_params()` for dataset Y. Must contain a `libsizes` sublist.
#' @param ancestry_X Character. Label used for dataset X in the plot. Default is `"Dataset X"`.
#' @param ancestry_Y Character. Label used for dataset Y in the plot. Default is `"Dataset Y"`.
#' @param title Plot title.
#' @param x_label X-axis label.
#' @param y_label Y-axis label.
#' 
#' @return A [ggplot2::ggplot()] object comparing the distribution of library sizes between datasets.
#'
#' @examples
#' \dontrun{
#' estimates_X <- estimate_params(count_matrix_X)
#' estimates_Y <- estimate_params(count_matrix_Y)
#' plot_estimated_library_size(estimates_X, estimates_Y)
#' }
#'
#' @export
plot_estimated_library_size <- function(
  estimates_X,
  estimates_Y,
  ancestry_X = "Dataset X",
  ancestry_Y = "Dataset Y",
  title = NULL,
  x_label = NULL,
  y_label = NULL,
){

  df <- data.frame(
    libsize = c(
      estimates_X$libsizes$libsizes, 
      estimates_Y$libsizes$libsizes
    ),
    dataset = rep(
      x = c(
        ancestry_X, 
        ancestry_Y
      ),
      times = c(
        length(estimates_X$libsizes$libsizes), 
        length(estimates_Y$libsizes$libsizes)
        )
    )
  )

  p <- ggplot(
    data = df, 
    mapping = aes(
        x = dataset, 
        y = libsize
      )
    ) +
    geom_violin(
      trim = FALSE
    ) +
    geom_boxplot(
      width = 0.15, 
      outlier.shape = NA
    )

    # Final style
    p <- p + labs(
      title = title,
      x = ifelse(is.null(x_label), "Ancestry group", x_label),
      y = ifelse(is.null(y_label), "Estimated library size", y_label)
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend()

  return(p)
}