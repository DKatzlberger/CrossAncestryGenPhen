#' Plot Sample Counts by Stratified Variable and Dataset Split
#'
#' Creates a bar plot showing the number of samples in each dataset split
#' (train, test, inference), grouped and colored by stratification variables.
#' Each stratification column is shown in a separate facet for clear comparison.
#'
#' This is useful for visually checking balance or representation across
#' strata (e.g., ancestry or condition) within each split.
#'
#' @param MX A data.frame containing metadata for the test set (usually output from `split_stratified_ancestry_sets()`).
#' @param MY A data.frame containing metadata for the inference set.
#' @param MR A data.frame containing metadata for the train set.
#' @param g_col Character vector of column names in the metadata to stratify and plot by (e.g., c("ancestry", "sex")).
#' @param title Optional character string to use as the plot title.
#'
#' @return A ggplot2 object showing counts per stratum and dataset split, faceted by stratification variable.
#' @export
plot_stratified_sets <- function(
  MX, 
  MY, 
  MR, 
  g_col,
  title = NULL
) {
  # Add dataset split labels
  train_M <- MR
  test_M  <- MX
  infer_M <- MY

  train_M$set <- "R"
  test_M$set  <- "X"
  infer_M$set <- "Y"

  all_M <- rbind(train_M, test_M, infer_M)
  all_M$set <- factor(all_M$set, levels = c("R", "X", "Y"))

  # Combine counts for each stratification column separately
  plot_data_list <- lapply(g_col, function(col) {
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

