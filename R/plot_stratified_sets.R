#' Plot Sample Size by Stratified Variable and Dataset Split
#'
#' Visualizes the number of samples in train, test, and inference sets,
#' with bars grouped by dataset split and filled by stratum values. Each
#' stratification column gets its own facet for clearer interpretation.
#'
#' @param result The output from `split_stratified_ancestry_sets()`.
#' @param stratify_cols Character vector of metadata column names used for stratification.
#' @param point_size Numeric value to control fontsize.
#'
#' @return A ggplot2 object with dodged bars grouped by set, filled by strata values, and faceted by stratify column.
#' @export
plot_stratified_sets <- function(
  x,
  stratify_cols,
  title = NULL,
  point_size = 0.5
) {
  # Add dataset split labels
  train_M <- x$train$M
  test_M  <- x$test$M
  infer_M <- x$inference$M

  train_M$set <- "Train"
  test_M$set  <- "Test"
  infer_M$set <- "Inference"

  all_M <- rbind(train_M, test_M, infer_M)
  all_M$set <- factor(all_M$set, levels = c("Train", "Test", "Inference"))

  # Combine counts for each stratification column separately
  plot_data_list <- lapply(stratify_cols, function(col) {
    tab <- as.data.frame(table(set = all_M$set, stratum = all_M[[col]]))
    tab$variable <- col
    names(tab)[2] <- "value"
    tab
  })

  plot_data <- do.call(rbind, plot_data_list)

  # Plot: dodged bars per split, fill by stratum value, facet by stratification variable
  p <- ggplot(
    data = plot_data, 
    aes(
      x = set, 
      y = Freq, 
      fill = value
      )
    ) +
    geom_col(
      position = "dodge",
      color = "black"
    ) +
    facet_wrap(
      ~ variable
    ) +
    labs(
      title = title,
      x = NULL,
      y = "Count",
      fill = "Stratum value"
    ) +
    theme_nature_fonts() +
    theme_small_legend() +
    theme_white_background() 

  return(p)
}

