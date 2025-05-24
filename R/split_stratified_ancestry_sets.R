#' Split Expression and Metadata into Train, Test, and Inference Sets
#'
#' Subsamples the overrepresented ancestry to match the underrepresented ancestry,
#' stratifying by metadata columns. Returns train, test, and inference sets.
#' Issues a warning when any strata from the underrepresented group
#' cannot be matched due to missing or insufficient samples in the overrepresented group.
#'
#' @param X A numeric matrix (samples x genes) for the overrepresented ancestry group.
#' @param Y A numeric matrix for the underrepresented ancestry group.
#' @param MX A data.frame of metadata corresponding to X (same number of rows).
#' @param MY A data.frame of metadata corresponding to Y (same number of rows).
#' @param stratify_cols A character vector of column names in metadata to stratify by.
#' @param seed An optional integer used to set the random seed.
#'
#' @return A list with named elements:
#' \describe{
#'   \item{train}{List with `X` and `M` for subsampled overrepresented group.}
#'   \item{test}{List with `X` and `M` for the remaining overrepresented samples.}
#'   \item{inference}{List with `X` and `M` for the underrepresented ancestry.}
#'   \item{strata_info}{Stratum match summary with `usable`, `missing`, and `insufficient`.}
#' }
#' @export
split_stratified_ancestry_sets <- function(
  X,
  Y,
  MX,
  MY,
  stratify_cols = NULL,
  seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Generate strata identifiers
  strata_X <- get_strata(
    M = MX,
    stratify_cols = stratify_cols
  )

  strata_Y <- get_strata(
    M = MY,
    stratify_cols = stratify_cols
  )

  # Determine feasible strata
  strata_info <- check_strata_feasibility(
    strata_X = strata_X,
    strata_Y = strata_Y
  )

  if (length(strata_info$missing) > 0 || length(strata_info$insufficient) > 0) {
    warning(
      "Some strata in the underrepresented ancestry (Y) could not be matched:\n",
      if (length(strata_info$missing) > 0) {
        paste0("- Missing in X: ", paste(strata_info$missing, collapse = ", "), "\n")
      },
      if (length(strata_info$insufficient) > 0) {
        paste0("- Too few in X: ", paste(strata_info$insufficient, collapse = ", "), "\n")
      }
    )
  }

  # Filter to usable strata
  mask_X <- strata_X %in% strata_info$usable
  mask_Y <- strata_Y %in% strata_info$usable

  X_sub <- apply_mask(
    X = X,
    M = MX,
    mask = mask_X
  )

  Y_sub <- apply_mask(
    X = Y,
    M = MY,
    mask = mask_Y
  )

  # Sample from overrepresented to match underrepresented
  sampled_indices <- stratified_sample_indices(
    strata = strata_X[mask_X],
    target_counts = table(strata_Y[mask_Y])
  )

  train <- list(
    X = X_sub$X[sampled_indices, , drop = FALSE],
    M = X_sub$M[sampled_indices, , drop = FALSE]
  )

  test <- list(
    X = X_sub$X[-sampled_indices, , drop = FALSE],
    M = X_sub$M[-sampled_indices, , drop = FALSE]
  )

  inference <- list(
    X = Y_sub$X,
    M = Y_sub$M
  )

  list(
    train = train,
    test = test,
    inference = inference,
    strata_info = strata_info
  )
}

