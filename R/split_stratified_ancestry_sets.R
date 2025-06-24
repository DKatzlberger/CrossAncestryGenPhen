#' Split Expression and Metadata into Train, Test, and Inference Sets
#'
#' Performs stratified sampling of an overrepresented group (X) to match the distribution of
#' an underrepresented group (Y) based on one or more grouping variables (e.g., ancestry).
#'
#' Sample IDs must be provided as rownames in both expression matrices. If `id_col` is specified,
#' it must match those rownames exactly in the metadata.
#'
#' @param X A numeric matrix (samples × features) for the overrepresented group. Row names must be sample IDs.
#' @param Y A numeric matrix (samples × features) for the underrepresented group. Row names must be sample IDs.
#' @param MX A data frame of metadata for `X`, with the same number of rows. Must align with `X` either by row order or `id_col`.
#' @param MY A data frame of metadata for `Y`, same format and row count as `Y`.
#' @param g_col Optional character vector of column names in metadata to stratify by (e.g., ancestry, sex). If `NULL`, all samples are grouped together.
#' @param id_col Optional column name in `MX` and `MY` that holds sample IDs. If `NULL`, metadata is matched to expression matrices by rownames.
#' @param seed Optional numeric seed for reproducibility of sampling.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{train}{A list with matrix, metadata, and sample IDs for the training set.}
#'   \item{test}{A list with matrix, metadata, and sample IDs for the testing set.}
#'   \item{inference}{A list with matrix, metadata, and sample IDs for the inference set (group Y).}
#'   \item{strata_info}{A summary list of usable, missing, and insufficient strata.}
#' }
#' @export
split_stratified_ancestry_sets <- function(
  X, 
  Y, 
  MX, 
  MY,
  g_col = NULL,
  id_col = NULL,
  seed = NULL
) {
  
  if (!is.matrix(X) || !is.matrix(Y)) stop("X and Y must be matrices.")
  if (nrow(X) != nrow(MX) || nrow(Y) != nrow(MY)) stop("X/MX and Y/MY row counts must match.")
  if (is.null(rownames(X)) || is.null(rownames(Y))) {
    stop("Both X and Y must have rownames corresponding to sample IDs.")
  }
  if (!is.null(seed)) set.seed(seed)

  # Get sample IDs from rownames of matrices (now a strict requirement)
  ids_X <- rownames(X)
  ids_Y <- rownames(Y)

  # If id_col is provided, check consistency with metadata
  if (!is.null(id_col)) {
    if (!id_col %in% colnames(MX) || !id_col %in% colnames(MY)) {
      stop("id_col not found in metadata.")
    }
    if (!all(rownames(X) == MX[[id_col]])) {
      stop("rownames of X do not match MX[[id_col]]; ensure alignment.")
    }
    if (!all(rownames(Y) == MY[[id_col]])) {
      stop("rownames of Y do not match MY[[id_col]]; ensure alignment.")
    }
  } else {
    if (nrow(X) != nrow(MX) || nrow(Y) != nrow(MY)) {
      stop("When id_col is NULL, metadata must align with expression matrices by row order.")
    }
    # Use row order only; no further ID checking
  }


  # Create stratification factors
  get_strata <- function(M, stratify_cols) {
    M <- as.data.frame(M)
    if (is.null(stratify_cols)) {
      return(factor(rep("all", nrow(M))))
    }
    if (!all(stratify_cols %in% colnames(M))) {
      stop("Some stratify_cols not found in metadata.")
    }
    interaction(M[, stratify_cols, drop = FALSE], drop = TRUE)
  }

  strata_X <- get_strata(MX, g_col)
  strata_Y <- get_strata(MY, g_col)

  # Check strata feasibility
  count_X <- table(strata_X)
  count_Y <- table(strata_Y)
  matched <- intersect(names(count_X), names(count_Y))
  missing <- setdiff(names(count_Y), names(count_X))
  insufficient <- matched[count_Y[matched] > count_X[matched]]
  usable <- setdiff(matched, insufficient)

  strata_info <- list(
    usable = usable,
    missing = missing,
    insufficient = insufficient
  )

  # Mask usable strata
  mask_X <- strata_X %in% usable
  mask_Y <- strata_Y %in% usable

  X_sub <- X[mask_X, , drop = FALSE]
  MX_sub <- as.data.frame(MX)[mask_X, , drop = FALSE]
  Y_sub <- Y[mask_Y, , drop = FALSE]
  MY_sub <- as.data.frame(MY)[mask_Y, , drop = FALSE]
  ids_X_sub <- ids_X[mask_X]
  ids_Y_sub <- ids_Y[mask_Y]

  # Stratified sampling
  sampled_ids <- unlist(lapply(names(table(strata_Y[mask_Y])), function(stratum) {
    stratum_ids <- ids_X_sub[strata_X[mask_X] == stratum]
    n <- table(strata_Y[mask_Y])[stratum]
    sample(stratum_ids, size = min(length(stratum_ids), n))
  }))

  # Split into train/test
  mask_test <- ids_X_sub %in% sampled_ids
  mask_train <- !mask_test

  list(
    train = list(
      X = X_sub[mask_train, , drop = FALSE],
      M = MX_sub[mask_train, , drop = FALSE],
      ids = ids_X_sub[mask_train]
    ),
    test = list(
      X = X_sub[mask_test, , drop = FALSE],
      M = MX_sub[mask_test, , drop = FALSE],
      ids = ids_X_sub[mask_test]
    ),
    inference = list(
      X = Y_sub,
      M = MY_sub,
      ids = ids_Y_sub
    ),
    strata_info = strata_info
  )
}
