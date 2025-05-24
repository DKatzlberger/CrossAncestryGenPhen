#' Create Stratification Factor
#'
#' Combines metadata columns into a factor for stratified sampling.
#'
#' @param M A data.frame containing metadata.
#' @param stratify_cols A character vector of column names to stratify by.
#'
#' @return A factor indicating stratification strata.
#' @keywords internal
get_strata <- function(
  M,
  stratify_cols
) {
  if (is.null(stratify_cols)) {
    factor(rep("all", nrow(M)))
  } else {
    interaction(M[, stratify_cols], drop = TRUE)
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
check_strata_feasibility <- function(
  strata_X,
  strata_Y
) {
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

#' Stratified Sampling Indices
#'
#' Performs stratified random sampling from a vector of strata.
#'
#' @param strata A factor representing sample strata (from overrepresented group).
#' @param target_counts A named vector of counts to sample per stratum.
#'
#' @return An integer vector of row indices to sample.
#' @keywords internal
stratified_sample_indices <- function(
  strata,
  target_counts
) {
  unlist(lapply(names(target_counts), function(stratum) {
    candidates <- which(strata == stratum)

    if (length(candidates) == 0) {
      return(integer(0))
    }

    n <- target_counts[stratum]

    sample(candidates, size = min(length(candidates), n))
  }))
}

#' Apply Row Mask to Expression and Metadata
#'
#' Subsets expression matrix and metadata data.frame using the same row filter.
#'
#' @param X A numeric matrix of gene expression (samples x genes).
#' @param M A data.frame of sample-level metadata.
#' @param mask A logical or integer vector indicating which rows to keep.
#'
#' @return A list containing `expr` and `meta` components.
#' @keywords internal
apply_mask <- function(
  X,
  M,
  mask
) {
  list(
    X = X[mask, , drop = FALSE],
    M = M[mask, , drop = FALSE]
  )
}
