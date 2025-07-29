#' Plot sensitivity and specificity as a bar plot.
#'
#' Computes sensitivity and specificity for binary classification tasks.
#' Takes binary vectors of true and predicted labels.
#' Uses the `target` value to define the positive class.
#' Sensitivity and specificity are computed from TP, FP, FN, TN counts.
#' Plots the results as two bars with percentage labels.
#'
#' @param y_true Vector of true labels (numeric, factor, or logical).
#' @param y_selected Vector of predicted labels (numeric, factor, or logical).
#' @param target The value that should be treated as the positive class.
#' @param title Optional plot title.
#' @param x_label Optional x-axis label.
#' @param y_label Optional y-axis label (default: "Value label").
#'
#' @return A ggplot2 object showing sensitivity and specificity as bars.
#' @import ggplot2
#' @export
plot_sensitivity_specificity <- function(
  y_true,
  y_selected,
  target = 1,
  title = NULL,
  x_label = NULL,
  y_label = NULL
) {
  # Recode true and predicted labels relative to target
  y <- ifelse(y_true == target, 1, 0)
  y_hat <- ifelse(y_selected == target, 1, 0)
  
  # Convert to factors with consistent levels and labels
  y <- factor(y, levels = c(1, 0), labels = c("1", "0"))
  y_hat <- factor(y_hat, levels = c(1, 0), labels = c("1", "0"))
  
  # Compute confusion matrix
  cm <- table(Predicted = y_hat, True = y)
  
  # Extract counts safely
  TP <- ifelse(!is.na(cm["1", "1"]), cm["1", "1"], 0)
  FP <- ifelse(!is.na(cm["1", "0"]), cm["1", "0"], 0)
  FN <- ifelse(!is.na(cm["0", "1"]), cm["0", "1"], 0)
  TN <- ifelse(!is.na(cm["0", "0"]), cm["0", "0"], 0)
  
  # Compute metrics
  sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)
  specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), NA)
  
  # Build tidy DF for ggplot
  metrics_df <- data.frame(
    Metric = c("Sensitivity", "Specificity"),
    Value = c(sensitivity, specificity)
  )
  
  # Plot
  p <- ggplot(
    data = metrics_df, 
    mapping = aes(
        x = Metric, 
        y = Value
      )
    ) +
    geom_col(
      width = 0.5
    ) +
    geom_text(
      mapping = aes(
        label = scales::percent(Value, accuracy = 0.1)
      ),
      size = 2,
      vjust = -0.5
    ) +
    ylim(0, 1) +
    labs(
      title = title,
      y = ifelse(is.null(y_label), "Value", y_label),
      x = x_label
    ) +
    theme_nature_fonts() +
    theme_white_background() +
    theme_small_legend()
  
  return(p)
}
