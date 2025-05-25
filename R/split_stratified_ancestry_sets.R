#' Split Expression and Metadata into Train, Test, and Inference Sets
#'
#' Subsamples the overrepresented ancestry to match the underrepresented ancestry
#' for a test set, stratifying by metadata columns. The remaining overrepresented
#' samples are used for training. The underrepresented ancestry group is returned
#' as the inference set. This function is useful when training a model on the
#' majority ancestry and testing fairness or generalizability using both a balanced
#' subset and the underrepresented group.
#'
#' @param X A numeric matrix (samples × features) for the overrepresented ancestry group.
#' @param Y A numeric matrix (samples × features) for the underrepresented ancestry group.
#' @param MX A data.frame of metadata for X (same number of rows as X).
#' @param MY A data.frame of metadata for Y (same number of rows as Y).
#' @param stratify_cols A character vector of metadata column names to stratify on.
#' @param seed Optional integer to set random seed for reproducibility.
#'
#' @return A named list with:
#' \describe{
#'   \item{train}{List with `X` and `M` for the **remaining overrepresented samples**.}
#'   \item{test}{List with `X` and `M` for the **subsampled overrepresented group**, stratified to match the underrepresented group.}
#'   \item{inference}{List with `X` and `M` for the **underrepresented ancestry**.}
#'   \item{strata_info}{A summary list with usable, missing, and insufficient strata.}
#' }
#'
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
  strata_X <- get_strata(M = MX, stratify_cols = stratify_cols)
  strata_Y <- get_strata(M = MY, stratify_cols = stratify_cols)

  # Determine feasible strata
  strata_info <- check_strata_feasibility(strata_X = strata_X, strata_Y = strata_Y)

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

  X_sub <- apply_mask(X = X, M = MX, mask = mask_X)
  Y_sub <- apply_mask(X = Y, M = MY, mask = mask_Y)

  # Sample from overrepresented group to match underrepresented strata counts
  sampled_indices <- stratified_sample_indices(
    strata = strata_X[mask_X],
    target_counts = table(strata_Y[mask_Y])
  )

  # Subset = test; remaining = train
  test <- list(
    X = X_sub$X[sampled_indices, , drop = FALSE],
    M = X_sub$M[sampled_indices, , drop = FALSE]
  )

  train <- list(
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
