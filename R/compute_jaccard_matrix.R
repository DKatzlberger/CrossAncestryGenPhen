#' Compute Jaccard Similarity Matrix Between Iterations
#'
#' Calculates pairwise Jaccard indices between iterations based on sample usage.
#' Useful for assessing overlap and dependence in train/test/inference sets across iterations.
#'
#' The Jaccard index is computed as:
#' \deqn{J(i, j) = |A_i ∩ A_j| / |A_i ∪ A_j|}
#' where \eqn{A_i} and \eqn{A_j} are the sets of sample IDs used in the specified role for iterations i and j.
#'
#' @param id_usage A data.frame with columns \code{ids}, \code{role}, and \code{iteration},
#'        such as the output from \code{\link{track_sample_ids}}.
#' @param role A character string specifying which sample role to evaluate ("test", "train", or "inference").
#'
#' @return A square matrix of Jaccard indices with dimensions \code{iterations × iterations}.
#'
#' @examples
#' # Assume `id_usage` is created using multiple splits and track_sample_ids()
#' jaccard_test <- compute_jaccard_matrix(id_usage, role = "test")
#' jaccard_train <- compute_jaccard_matrix(id_usage, role = "train")
#'
#' @export
compute_jaccard_matrix <- function(
  id_usage, 
  role
) {

  usage <- id_usage[id_usage$role == role, ]
  sample_matrix <- table(usage$ids, usage$iteration)
  n_iter <- ncol(sample_matrix)
  jaccard_matrix <- matrix(NA_real_, n_iter, n_iter)
  
  # Jaccard index computation
  for (i in seq_len(n_iter)) {
    for (j in seq_len(n_iter)) {
      a <- sample_matrix[, i]
      b <- sample_matrix[, j]
      jaccard_matrix[i, j] <- sum(a & b) / sum(a | b)
    }
  }

  colnames(jaccard_matrix) <- rownames(jaccard_matrix) <- as.character(seq_len(n_iter))
  return(jaccard_matrix)
}
