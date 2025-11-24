#' Filter samples using a metadata column and selected levels
#'
#' Filters the sample matrix and metadata to only include rows whose metadata
#' column matches specified levels.
#'
#' @param X Matrix with samples as rows.
#' @param M Data frame of metadata with matching row names.
#' @param s_col Metadata column name to filter by.
#' @param s_levels Character vector of allowed levels.
#' @param verbose Logical; print summary of filtering.
#'
#' @return A list containing filtered X and M.
#' @export
filter_sample_type <- function(
  X,
  M,
  s_col,
  s_levels,
  verbose = TRUE
){

  ## --- Validate input ---
  stopifnot(is.matrix(X), !is.null(rownames(X)))
  stopifnot(is.data.frame(M), !is.null(rownames(M)))
  if (!identical(rownames(X), rownames(M))) stop("Row names of `M` and `X` must match exactly (same sampleIDs).")


  ## --- Filter meta by columns ---
  meta_cols <- unique(c(s_col))
  meta_cols <- meta_cols[meta_cols %in% colnames(M)]

  if (length(meta_cols) == 0) {
    stop("[filter_sample_type] No valid metadata columns found (check s_col).")
  }

  ## --- Validate requested levels ---
  missing <- setdiff(s_levels, unique(M[[s_col]]))
  if (length(missing) > 0) {
    stop(sprintf("[filter_sample_type] Requested sample level(s) NOT found: %s", paste(missing, collapse = ", ")))
  }

  ## --- Filter by s_col, s_levels ---
  m_sub <- M[M[[s_col]] %in% s_levels, , drop = FALSE]
  if (nrow(m_sub) == 0L) stop(sprintf("No samples found for sample levels '%s'.", patse(s_levels, collapse = ", ")))

  # Matrix filter
  X_sub <- X[rownames(m_sub), , drop = FALSE]

  ## --- Verbose message ---
if (verbose) {
  b_counts <- table(M[[s_col]])
  a_counts <- table(m_sub[[s_col]])
  message("\nFilter sample type summary:")
  message(sprintf("Before filtering: %d samples (%s)", nrow(M), paste(names(b_counts), b_counts, sep = " = ", collapse = ", ")))
  message(sprintf("After filtering:  %d samples (%s)", nrow(m_sub), paste(names(a_counts), a_counts, sep = " = ", collapse = ", ")))
}


  ## --- Return ---
  return(
    list(
      matr = X_sub,
      meta = m_sub
    )
  )
}
