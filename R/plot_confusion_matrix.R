#' Plot a binary confusion matrix as a heatmap with TP, FP, FN, TN labels.
#'
#' Computes a confusion matrix for binary classification tasks.
#' Takes binary vectors of true and predicted labels.
#' Uses the `target` value to define the positive class.
#' Cells are labeled with TP, FP, FN, TN and counts.
#' Class labels shown as "1" (positive) and "0" (negative).
#' Uses ggplot2 with a simple blue color scale.
#'
#' @param y_true Vector of true labels (numeric, factor, or logical).
#' @param y_selected Vector of predicted labels (numeric, factor, or logical).
#' @param target The value that should be treated as the positive class.
#' @param title Optional plot title.
#' @param x_label Optional x-axis label (default: "True Label").
#' @param y_label Optional y-axis label (default: "Predicted Label").
#'
#' @return A ggplot2 object showing the confusion matrix heatmap.
#' @import ggplot2
#' @export
plot_confusion_matrix <- function(
  y_true,
  y_selected,
  target = 1,
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  caption = NULL
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
  
  # Build tidy DF for ggplot
  cm_df <- data.frame(
    Predicted = factor(rep(c("1", "0"), each = 2), levels = c("1", "0")),
    True = factor(rep(c("1", "0"), times = 2), levels = c("1", "0")),
    Count = c(TP, FP, FN, TN),
    Label = c("TP", "FP", "FN", "TN")
  )
  cm_df$Text <- paste0(cm_df$Label, "\n", cm_df$Count)
  
  # Plot
  p <- ggplot(
    data = cm_df, 
    mapping = aes(
      x = True, 
      y = Predicted, 
      fill = Count
      )
    ) +
    geom_tile(
      color = "white"
    ) +
    geom_text(
      aes(
        label = Text
      ),
      size = 2, 
      color = "white"
    ) +
    coord_fixed(
      ratio = 1
    ) +
    labs(
      title = title,
      x = ifelse(is.null(x_label), "True label", x_label),
      y = ifelse(is.null(y_label), "Predicted label", y_label),
      caption = ifelse(is.null(caption), paste0("\n", target, ": Interaction, 0: No interaction"), caption),
      fill = NULL
    ) +
    theme_CrossAncestryGenPhen()
    theme(
      panel.grid = element_blank()
    )
  
  return(p)
}
