#' Split Expression and Metadata into Train, Test, and Inference Sets
#'
#' Subsamples the overrepresented ancestry (e.g., EUR) to match the underrepresented ancestry
#' (e.g., AFR) for a test set, stratifying by metadata columns. The remaining overrepresented
#' samples are used for training. The underrepresented ancestry group is returned as the inference set.
#'
#' This version supports tracking and subsetting using either rownames or a specified sample ID column.
#'
#' @param X A numeric matrix (samples × features) for the overrepresented ancestry group.
#' @param Y A numeric matrix (samples × features) for the underrepresented ancestry group.
#' @param MX A data.frame of metadata for X (same number of rows as X).
#' @param MY A data.frame of metadata for Y (same number of rows as Y).
#' @param id_col Optional string specifying the column name in metadata that holds sample IDs. If NULL, rownames are used.
#' @param stratify_cols A character vector of metadata column names to stratify on.
#' @param seed Optional integer to set random seed for reproducibility.
#'
#' @return A named list with:
#' \describe{
#'   \item{train}{List with `X`, `M`, and `ids` for the remaining overrepresented samples.}
#'   \item{test}{List with `X`, `M`, and `ids` for the subsampled overrepresented group.}
#'   \item{inference}{List with `X`, `M`, and `ids` for the underrepresented ancestry group.}
#'   \item{strata_info}{A summary list with usable, missing, and insufficient strata.}
#' }
#'
#' @export
split_stratified_ancestry_sets <- function(
  X,
  Y,
  MX,
  MY,
  id_col = NULL,
  stratify_cols = NULL,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  # Get sample IDs
  ids_X <- if (!is.null(id_col)) MX[[id_col]] else rownames(MX)
  ids_Y <- if (!is.null(id_col)) MY[[id_col]] else rownames(MY)
  if (is.null(ids_X) || is.null(ids_Y)) stop("Sample IDs must be in rownames or specified by `id_col`.")

  # Generate strata
  strata_X <- get_strata(MX, stratify_cols)
  strata_Y <- get_strata(MY, stratify_cols)

  # Identify usable strata
  strata_info <- check_strata_feasibility(strata_X, strata_Y)

  # Print summary of strata
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

  X_sub <- apply_mask(X, MX, mask_X)
  Y_sub <- apply_mask(Y, MY, mask_Y)
  ids_X_sub <- ids_X[mask_X]
  ids_Y_sub <- ids_Y[mask_Y]

  # Sample based on sample IDs
  sampled_ids <- stratified_sample_ids(
    strata = strata_X[mask_X],
    target_counts = table(strata_Y[mask_Y]),
    ids = ids_X_sub
  )

  # Masks for test/train
  mask_test <- ids_X_sub %in% sampled_ids
  mask_train <- !mask_test

  # Subsets
  test <- list(
    X = X_sub$X[mask_test, , drop = FALSE],
    M = X_sub$M[mask_test, , drop = FALSE],
    ids = ids_X_sub[mask_test]
  )

  train <- list(
    X = X_sub$X[mask_train, , drop = FALSE],
    M = X_sub$M[mask_train, , drop = FALSE],
    ids = ids_X_sub[mask_train]
  )

  inference <- list(
    X = Y_sub$X,
    M = Y_sub$M,
    ids = ids_Y_sub
  )

  list(
    train = train,
    test = test,
    inference = inference,
    strata_info = strata_info
  )
}


