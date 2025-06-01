#' Create Stratification Factor
#'
#' Combines metadata columns into a factor for stratified sampling.
#'
#' @param M A data.frame containing metadata.
#' @param stratify_cols A character vector of column names to stratify by.
#'
#' @return A factor indicating stratification strata.
#' @keywords internal
get_strata <- function(M, stratify_cols) {
  if (is.null(stratify_cols)) {
    factor(rep("all", nrow(M)))
  } else {
    interaction(M[, stratify_cols, drop = FALSE], drop = TRUE)
  }
}

#' Check Stratification Feasibility
#'
#' Identifies usable, missing, or insufficient strata based on availability in each group.
#'
#' @param strata_X A factor for overrepresented group strata (e.g., EUR).
#' @param strata_Y A factor for underrepresented group strata (e.g., AFR).
#'
#' @return A list with `usable`, `missing`, and `insufficient` strata names.
#' @keywords internal
check_strata_feasibility <- function(strata_X, strata_Y) {
  count_X <- table(strata_X)
  count_Y <- table(strata_Y)

  matched <- intersect(names(count_X), names(count_Y))
  missing <- setdiff(names(count_Y), names(count_X))
  insufficient <- matched[count_Y[matched] > count_X[matched]]

  list(
    usable = setdiff(matched, insufficient),
    missing = missing,
    insufficient = insufficient
  )
}

#' Stratified Sampling by Sample IDs
#'
#' Performs stratified sampling from a vector of sample IDs grouped by strata.
#'
#' @param strata A factor indicating strata for each sample (must match `ids` length).
#' @param target_counts A named vector of how many samples to draw per stratum.
#' @param ids A character or numeric vector of sample IDs, same length as `strata`.
#'
#' @return A character or numeric vector of sampled IDs.
#' @keywords internal
stratified_sample_ids <- function(strata, target_counts, ids) {
  unlist(lapply(names(target_counts), function(stratum) {
    stratum_ids <- ids[strata == stratum]
    n <- target_counts[stratum]
    sample(stratum_ids, size = min(length(stratum_ids), n))
  }))
}

#' Apply Row Mask to Expression and Metadata
#'
#' Subsets expression matrix and metadata using the same row mask.
#'
#' @param X A numeric matrix (samples x features).
#' @param M A data.frame of metadata (samples x variables).
#' @param mask A logical or integer vector of rows to keep.
#'
#' @return A list with `X` and `M` subsets.
#' @keywords internal
apply_mask <- function(X, M, mask) {
  list(
    X = X[mask, , drop = FALSE],
    M = M[mask, , drop = FALSE]
  )
}
