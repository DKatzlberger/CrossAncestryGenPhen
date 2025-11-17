#' Collapse duplicate genes with sum and remove NAs from count matrix
#'
#' @param matr Numeric matrix (samples × features).
#' @param verbose Logical; show summary (default TRUE).
#'
#' @return Matrix with NAs removed and duplicate genes summed.
#' @export
gdc_gene_level_count_matrix <- function(
  matr,
  verbose = TRUE
){

  ## --- Check for matrix ---
  if (!is.matrix(matr)) stop("Input must be a matrix (samples × features).")


  ## --- Record start dimensions ---
  n_start_samples  <- nrow(matr)
  n_start_features <- ncol(matr)


  ## --- Remove feature with any NA ---
  na_mask <- colSums(is.na(matr)) == 0
  matr    <- matr[, na_mask, drop = FALSE]
  n_removed_na <- ncol(matr)

  ## --- Aggregate duplicate features by SUM ---
  feature_names <- colnames(matr)
  agg_matr <- t(rowsum(t(matr), group = feature_names))


  ## --- Verbose message ---
  if (verbose){
    message("\nGDC gene-level count matrix summary:")
    message(sprintf("Features:   %d -> %d (no NAs) -> %d (sum)", n_start_features, n_removed_na, ncol(agg_matr)))
    message(sprintf("Dimensions: %d (samples) x %d (features)", nrow(agg_matr), ncol(agg_matr)))
  }


  ## --- Return ---
  return(agg_matr)
}
