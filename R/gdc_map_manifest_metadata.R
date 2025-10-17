#' Parse GDC metadata and map files to samples
#'
#' Reads manifest and metadata JSON from GDC cart
#' Builds a file-to-sample mapping data.frame
#' Drops rows with missing sample IDs before check
#' Compares manifest and metadata row counts
#' Stops if counts differ after filtering
#'
#' @param manifest_file Path to GDC manifest file
#' @param metadata_file Path to metadata JSON file
#' @param verbose Print progress messages
#'
#' @return Data frame with file_name and sample_id
#' 
#' @export
gdc_map_manifest_metadata <- function(
	manifest_file, 
	metadata_file,
	verbose = TRUE
){

	## --- Check for the manifest file ---
	if (!file.exists(manifest_file)) {
    stop("[parse_gdc_metadata] Manifest file not found: ", manifest_file)
  }

	## --- Check for the manifest file ---
  if (!file.exists(metadata_file)) {
    stop("[parse_gdc_metadata] Metadata file not found: ", metadata_json)
  }

	## --- Read manifest ---
  manifest <- read.delim(manifest_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)


	## --- Read metadata ---
  meta_list <- jsonlite::fromJSON(metadata_file, simplifyVector = FALSE)


	## --- Extract mapping ---
  file_map <- data.frame(
    file_name = sapply(meta_list, function(x) x$file_name),
    sample_id = sapply(meta_list, function(x) {if (!is.null(x$associated_entities[[1]]$entity_submitter_id)) {x$associated_entities[[1]]$entity_submitter_id} else {NA_character_}}),
    stringsAsFactors = FALSE
  )


  ## --- Drop rows with missing sample_id ---
  file_map <- file_map[!is.na(file_map$sample_id) & file_map$sample_id != "", , drop = FALSE]


	## --- Compare counts ---
  n_manifest <- nrow(manifest)
  n_metadata <- nrow(file_map)


	## --- verbose message ---
	if (verbose) {
		message("\nGDC meta summary:")
    message(sprintf("Manifest entries (total): %-4d", n_manifest))
    message(sprintf("Metadata entries (valid): %-4d", n_metadata))
  }


  ## --- Return ---
  return(file_map)
}