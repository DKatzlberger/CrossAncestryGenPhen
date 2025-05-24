#' Plot Sample Size by Set Colored by Stratum
#'
#' Visualizes the number of samples in train, test, and inference sets, colored by stratum.
#'
#' @param result The output from `split_stratified_ancestry_sets()`.
#' @param stratify_cols Character vector of metadata column names used for stratification.
#'
#' @return A ggplot object showing sample count per set, stacked by stratum.
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_brewer labs theme_minimal theme element_text
#' @export
plot_stratified_sets <- function(
  result,
  stratify_cols
) {
  # Extract metadata and add set labels
  train_M <- result$train$M
  test_M  <- result$test$M
  infer_M <- result$inference$M

  train_M$set <- "Train"
  test_M$set  <- "Test"
  infer_M$set <- "Inference"

  # Combine all metadata
  all_M <- rbind(train_M, test_M, infer_M)

  # Create stratum variable using interaction
  all_M$stratum <- interaction(all_M[, stratify_cols], drop = TRUE)

  # Count samples by set and stratum
  tab <- table(all_M$set, all_M$stratum)
  plot_data <- as.data.frame(tab)
  names(plot_data) <- c("set", "stratum", "count")

  # Plot stacked bar chart
  p <- ggplot(
    data = plot_data, 
    aes(
      x = set, 
      y = count, 
      fill = stratum
      )
      ) +
    geom_col(
      position = "stack"
      ) +
    scale_fill_brewer(
      palette = "Set3"
      ) +
    labs(
      title = "Sample Sizes by Set, Colored by Stratum",
      x = "Dataset Split",
      y = "Sample Count",
      fill = "Stratum"
      ) +
    theme_minimal(
      base_size = 12
      ) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      strip.text = element_text(face = "bold")
    )

  return(p)
}


