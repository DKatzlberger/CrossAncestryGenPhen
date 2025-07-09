#' Compute sensitivity and specificity for feature selection
#'
#' Compares selected features to known true features and calculates 
#' TP, FP, FN, TN, sensitivity, and specificity.
#'
#' @param true Character vector of true significant features (ground truth).
#' @param selected Character vector of features selected by your method.
#' @param all Character vector of all features tested or eligible to select.
#'
#' @return A data frame with TP, FP, FN, TN, sensitivity, and specificity.
#'
#' @details Sensitivity = TP / (TP + FN). Specificity = TN / (TN + FP).
#' Use only features actually tested in `all` to get meaningful TNs.
#'
#' @examples
#' true <- c("gene1", "gene2", "gene3")
#' selected <- c("gene2", "gene3", "gene4")
#' all <- paste0("gene", 1:10)
#' compute_sens_spec(true, selected, all)
#'
#' @export
compute_sens_spec <- function(
    true, 
    selected, 
    all
) {

  true_features <- unique(true)
  selected_features <- unique(selected)
  all_features <- unique(all)

  TP <- length(intersect(true_features, selected_features))
  FP <- length(setdiff(selected_features, true_features))
  FN <- length(setdiff(true_features, selected_features))
  TN <- length(setdiff(all_features, union(true_features, selected_features)))

  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)

  return(
    data.frame(
      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN,
      Sensitivity = sensitivity,
      Specificity = specificity
    )
  )
}
