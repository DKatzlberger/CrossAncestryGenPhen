#' Map metadata columns and values
#'
#' Applies a user-defined mapping table to harmonize column names and values in a metadata data frame.
#'
#' @param file_map A data frame containing the metadata to be mapped.
#' @param map A data frame with columns \code{old_colname}, \code{new_colname}, \code{old_colvalue}, and \code{new_colvalue} specifying how to map columns and values.
#' @param verbose Logical; if \code{TRUE}, prints progress messages.
#'
#' @return A data frame with updated columns and values according to the mapping table.
#' @export
gdc_map_metadata <- function(
  file_map,
  map,
  verbose = TRUE
){
  ## --- Input validation ---
  stopifnot(all(c("old_colname", "new_colname", "old_colvalue", "new_colvalue") %in% names(map)))

  ## --- Trim ALL whitespace in metadata BEFORE mapping ---
  file_map <- as.data.frame(lapply(file_map, function(x) {
    if (is.character(x)) trimws(x) else x
  }), stringsAsFactors = FALSE)


  ## --- Trim whitespace in the mapping table itself ---
  map$old_colvalue <- trimws(map$old_colvalue)
  map$new_colvalue <- trimws(map$new_colvalue)
  map$old_colname  <- trimws(map$old_colname)
  map$new_colname  <- trimws(map$new_colname)

  
  ## --- Apply mapping ---
  for (i in seq_len(nrow(map))) {
    old_col  <- map$old_colname[i]
    new_col  <- map$new_colname[i]
    old_val  <- map$old_colvalue[i]
    new_val  <- map$new_colvalue[i]
    
    if (!old_col %in% names(file_map)) {
      message(sprintf("Skipping: column '%s' not found in file_map", old_col))
      next
    }

    # Trim the active column BEFORE matching ---
    if (is.character(file_map[[old_col]])) {
      file_map[[old_col]] <- trimws(file_map[[old_col]])
    }
    
    # Create the new column if it doesn’t exist yet
    if (!new_col %in% names(file_map)) {
      file_map[[new_col]] <- file_map[[old_col]]
      if (verbose) message(sprintf("\nCreated new column %-25s from %-25s", paste0("'", new_col, "'"), paste0("'", old_col, "'")))
    }
    
    # Replace values
    matches <- file_map[[old_col]] == old_val
    if (any(matches, na.rm = TRUE)) {
      file_map[[new_col]][matches] <- new_val
      if (verbose) message(sprintf("Mapped %-25s == %-15s → %-15s in %-25s", paste0("'", old_col, "'"), paste0("'", old_val, "'"), paste0("'", new_val, "'"), paste0("'", new_col, "'")))
    }
  }
  
  ## --- Return ---
  return(file_map)
}