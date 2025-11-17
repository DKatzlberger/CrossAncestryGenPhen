#' Align metadata and matrix by IDs
#'
#' Reorders matrix rows to match metadata.
#' Checks ID sets for consistency.
#'
#' @param meta Metadata table.
#' @param matr Omics matrix.
#' @param id ID column name.
#' @param tech Technology name.
#' @param verbose Show messages.
#'
#' @return Named list with meta and matr.
#' @export
gdc_omics_layer <- function(
  meta,
  matr,
  id,
  tech,
  verbose = TRUE
){

  ## --- Check ID column ---
  if (!id %in% colnames(meta)) {
    stop(sprintf("ID column '%s' not found in metadata.", id))
  }

  ## --- Extract IDs ---
  meta_id <- meta[[id]]
  matr_id <- rownames(matr)


  ## --- Check ID agreement (ignoring order) ---
  if (!setequal(meta_id, matr_id)) {
    stop(sprintf("IDs do not match between metadata and matrix."))
  }


  ## --- Reorder matrix if needed ---
  if (!identical(meta_id, matr_id)) {
    matr <- matr[meta_id, , drop = FALSE]
    message(sprintf("Reordered matrix rows to match metadata."))
  }


  ## --- Assign rownames to metadata ---
  rownames(meta) <- meta_id
  stopifnot(identical(rownames(meta), rownames(matr)))


  ## --- Verbose message ---
  message(sprintf("\nOmics %s summary:", tech))
  message(sprintf("Number of samples: %d", length(rownames(meta))))


  ## --- Create named list output ---
  out <- list(
    meta = meta,
    matr = matr
  )
  

  ## --- Return ---
  return(structure(list(out), names = tech))
}