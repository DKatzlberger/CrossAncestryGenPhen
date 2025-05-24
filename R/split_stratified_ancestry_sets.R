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
#' @return A list with three named elements:
#' \describe{
#'   \item{train}{List with `expr` and `meta` for subsampled overrepresented group.}
#'   \item{test}{List with `expr` and `meta` for the remaining overrepresented samples.}
#'   \item{inference}{List with `expr` and `meta` for the underrepresented ancestry.}
#' }
#' @importFrom stats interaction
#' @export
split_stratified_ancestry_sets <- function(
  X,
  Y,
  MX,
  MY,
  stratify_cols = NULL,
  seed = NULL
) {
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Generate strata for both groups based on metadata
  strata_X <- get_strata(
    meta = MX,
    stratify_cols = stratify_cols
  )

  strata_Y <- get_strata(
    meta = MY,
    stratify_cols = stratify_cols
  )

  # Check which strata are matchable
  strata_info <- check_strata_feasibility(
    strata_X = strata_X,
    strata_Y = strata_Y
  )

  # Warn user about any unmatched strata
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

  # Filter samples to only usable strata
  mask_X <- strata_X %in% strata_info$usable
  mask_Y <- strata_Y %in% strata_info$usable

  X_sub <- apply_mask(
    expr = X,
    meta = MX,
    mask = mask_X
  )

  Y_sub <- apply_mask(
    expr = Y,
    meta = MY,
    mask = mask_Y
  )

  # Sample from overrepresented group to match counts in underrepresented group
  sampled_indices <- stratified_sample_indices(
    strata = strata_X[mask_X],
    target_counts = table(strata_Y[mask_Y])
  )

  # Define train and test sets
  train <- apply_mask(
    expr = X_sub$expr,
    meta = X_sub$meta,
    mask = sampled_indices
  )

  test <- apply_mask(
    expr = X_sub$expr,
    meta = X_sub$meta,
    mask = setdiff(seq_len(nrow(X_sub$expr)), sampled_indices)
  )

  # Inference set is full Y subset
  inference <- Y_sub

  # Return all splits
  list(
    train = train,
    test = test,
    inference = inference
  )
}
