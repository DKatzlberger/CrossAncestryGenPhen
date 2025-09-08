#' Split Expression and Metadata into Reference (R), Subset (X), and Inference (Y) Sets 
#'
#' Performs stratified sampling of an overrepresented group (X, e.g. EUR) to match the distribution of 
#' an underrepresented group (Y, e.g. AFR) based on a grouping variable (e.g., condition).
#'
#' @param X Numeric matrix or data.frame of features for cohort X; rows are samples and must align with MX.
#' @param Y Numeric matrix or data.frame of features for cohort Y; rows are samples and must align with MY.
#' @param MX Data.frame with metadata for X.
#' @param MY Data.frame with metadata for Y.
#' @param g_col Name of the metadata column holding the stratification label.
#' @param a_col Name of the metadata column holding the ancestry label.
#' @param seed Optional numeric seed for reproducibility of sampling.
#' @param verbose Logical, whether to print messages.
#'
#' @return A list with the following elements (all matrices with rownames):
#' \describe{
#'   \item{R}{Reference set: remaining X after subsampling}
#'   \item{X}{Subset set: subsampled X matching Y}
#'   \item{Y}{Inference set: full Y, untouched}
#'   \item{strata_info}{list with usable/missing/insufficient strata}
#' }
#' @export
split_stratified_ancestry_sets <- function(
  X, 
  Y, 
  MX, 
  MY,
  g_col, 
  a_col,
  seed = NULL,
  verbose = TRUE
) {
  ## --- Seed ---
  if (!is.null(seed)) set.seed(seed)

  ## --- Input checks ---
  assert_input(
    X = X, 
    Y = Y,
    MX = MX, 
    MY = MY,
    g_col = g_col, 
    a_col = a_col
  )


  ## --- Factor setup ---
  ancestry_X <- unique(MX[[a_col]])
  ancestry_Y <- unique(MY[[a_col]])

  g_levels <- levels(MX[[g_col]])
  a_levels <- c(ancestry_X, ancestry_Y)

  if (length(g_levels) != 2 || length(a_levels) != 2) {
    stop("[split_stratified_ancestry_sets] Function supports only 2x2 designs (two levels in g_col × two levels a_col).")
  }

  ## --- Define strata (grouping) ---
  get_strata <- function(M, stratify_col) {
    interaction(M[[stratify_col]], drop = TRUE)
  }
  strata_X <- get_strata(MX, g_col)
  strata_Y <- get_strata(MY, g_col)

  ## --- Feasibility check ---
  count_X <- table(strata_X)
  count_Y <- table(strata_Y)

  matched      <- intersect(names(count_X), names(count_Y))
  missing      <- setdiff(names(count_Y), names(count_X))
  insufficient <- matched[count_Y[matched] > count_X[matched]]
  usable       <- setdiff(matched, insufficient)

  # Error if Y has strata that cannot be matched
  if (length(missing) > 0 || length(insufficient) > 0) {
    stop("[split_stratified_ancestry_sets] Y contains strata that X cannot match.\n",
         "Missing strata: ", paste(missing, collapse = ", "), "\n",
         "Insufficient strata: ", paste(names(insufficient), collapse = ", "))
  }

  strata_info <- list(
    usable = usable,
    missing = missing,
    insufficient = insufficient
  )

  ## --- Subset usable X ---
  mask_X <- strata_X %in% usable
  X_sub  <- as.matrix(X[mask_X, , drop = FALSE])
  MX_sub <- as.data.frame(MX[mask_X, , drop = FALSE])

  ids_X_sub <- rownames(X_sub)
  strata_X_sub <- strata_X[mask_X]

  ## --- Stratified sampling of X to match Y ---
  sampled_ids <- unlist(lapply(names(count_Y), function(stratum) {
    stratum_ids <- ids_X_sub[strata_X_sub == stratum]
    n <- count_Y[stratum]
    sample(stratum_ids, size = n, replace = FALSE)
  }), use.names = FALSE)

  mask_subset <- ids_X_sub %in% sampled_ids
  mask_ref    <- !mask_subset

  ## --- Reference set (R = remaining X) ---
  R_counts <- X_sub[mask_ref, , drop = FALSE]
  R_meta   <- MX_sub[mask_ref, , drop = FALSE]

  ## --- Subset set (X = sampled X) ---
  X_counts <- X_sub[mask_subset, , drop = FALSE]
  X_meta   <- MX_sub[mask_subset, , drop = FALSE]

  ## --- Inference set (Y = full Y) ---
  Y_counts <- Y
  Y_meta   <- MY

  ## --- Verbose summary ---
  if (verbose) {
    fmt_counts <- function(ids, M, g_col) {
      if (length(ids) == 0) return("")
      tab <- table(M[ids, g_col, drop = TRUE])
      paste(sprintf("%s: %-4d", names(tab), as.integer(tab)), collapse = " ")
    }

    message("\nStratified split summary:")
    message(sprintf("%s (Reference, R):    N: %-4d %s", ancestry_X, nrow(R_counts), fmt_counts(rownames(R_counts), R_meta, g_col)))
    message(sprintf("%s (Subset,    X):    N: %-4d %s", ancestry_X, nrow(X_counts), fmt_counts(rownames(X_counts), X_meta, g_col)))
    message(sprintf("%s (Inference, Y):    N: %-4d %s", ancestry_Y, nrow(Y_counts), fmt_counts(rownames(Y_counts), Y_meta, g_col)))
  }

  ## --- Return ---
  return(
    list(
      R = list(counts = R_counts, meta = R_meta, ids = rownames(R_counts)),
      X = list(counts = X_counts, meta = X_meta, ids = rownames(X_counts)),
      Y = list(counts = Y_counts, meta = Y_meta, ids = rownames(Y_counts)),
      strata_info = strata_info
    )
  )
}
