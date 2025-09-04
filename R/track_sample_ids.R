#' Track Sample Roles and IDs from a Split
#'
#' Extracts the sample IDs from a split object and labels their roles (test, train, inference),
#' associating each sample with its assigned role in a given iteration.
#'
#' @param split A split object returned from \code{split_stratified_ancestry_sets()}.
#' @param iteration An integer iteration number to annotate the split with.
#'
#' @return A data.frame with columns: \code{ids}, \code{role}, and \code{iteration}.
#'
#' @export
track_sample_ids <- function(
  split, 
  iteration
) {
  
  data.frame(
    ids = c(
      split$R$ids, 
      split$X$ids, 
      split$Y$ids
    ),
    role = c(
      rep("R", length(split$R$ids)),
      rep("X", length(split$X$ids)),
      rep("Y", length(split$Y$ids))
    ),
    iteration = iteration
  )
}
