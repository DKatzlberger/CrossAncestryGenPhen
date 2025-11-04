#' Download data from GDC via gdc-client
#'
#' Reads a manifest file and downloads all files
#' Uses current dir if out_dir is not specified
#' Supports optional token for protected data
#' Prints file and download counts when verbose
#'
#' @param manifest_file Path to GDC manifest txt file
#' @param metadata_file Path to GDC metadata json file
#' @param token_file Optional path to token file
#' @param download Wether to download (default: True)
#' @param file_dir Output dir, defaults to getwd()
#' @param client_path Path to gdc-client binary
#' @param verbose Print progress messages
#'
#' @return ID to file
#' 
#' @export
gdc_download_manifest <- function(
  manifest_file,
  metadata_file,
  token_file = NULL,
  download = TRUE,
  file_dir = NULL,
  client_path,
  verbose = TRUE
){


	## --- Check file existence ---
	if (!file.exists(manifest_file)) stop("[parse_gdc_metadata] Manifest file not found: ", manifest_file)
  if (!file.exists(metadata_file)) stop("[parse_gdc_metadata] Metadata file not found: ", metadata_file)


  ## --- Check correct file formats ---
  if (!grepl("\\.txt$", manifest_file, ignore.case = TRUE)) stop("[download_from_gdc] Manifest file must have extension '.txt': ", manifest_file)
  if (!grepl("\\.json$", metadata_file, ignore.case = TRUE)) stop("[download_from_gdc] Metadata file must have extension '.json': ", metadata_file)
  if (!is.null(token_file) && !grepl("\\.txt$", token_file, ignore.case = TRUE)) stop("[download_from_gdc] Token file must have extension '.txt': ", token_file)
  

  ## --- Set output directory ---
  if (is.null(file_dir)) {
    out_dir <- getwd()
  } else if (!dir.exists(file_dir)) {
    dir.create(file_dir, recursive = TRUE)
  }


  ## --- Read manifest ---
  manifest_data <- read.delim(manifest_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  n_manifest    <- nrow(manifest_data)


  ## --- Map manifest with meta ---
  file_map <- gdc_map_manifest(
    manifest_file = manifest_file, 
    metadata_file = metadata_file,
    verbose = TRUE
  )
  n_file_map <- nrow(file_map)


  ## --- Check existing folders ---
  expected_folders <- unique(file_map$FILE_ID)
  existing_folders <- list.dirs(file_dir, full.names = FALSE, recursive = FALSE)
  existing_folders <- basename(existing_folders)

  exists_flag <- expected_folders %in% existing_folders
  n_expected  <- length(expected_folders)
  n_existing  <- length(existing_folders)
  n_correct   <- sum(exists_flag)
  n_missing   <- n_expected - n_correct

  # Detect extra folders that don’t match the manifest
  n_extra_folders <- length(setdiff(existing_folders, expected_folders))


  ## --- Verbose message ---
  if (verbose) {
    message("\nGDC download summary:")
    message(sprintf("Entries:  %-9s  %-9s  %-9s", paste("Manifest:", n_manifest), paste("Metadata:", n_file_map), paste("directory:", file_dir)))
    message(sprintf("Existing: %-5d | Expected: %-5d | Missing: %-5d", n_existing, n_expected, n_missing))
  }

  # Warn if extra / mismatched folders are found
  if (n_extra_folders > 0) {
    stop(sprintf("Warning: %d entrie(s) in specified directory not match any FILE ID in the manifest.", n_extra_folders))
  }


  ## --- Base command ---
  cmd <- sprintf('%s download -m "%s" -d "%s"', client_path, manifest_file, file_dir)


  ## --- Token (for protected data) ---
  if (!is.null(token_file)) {
    if (!file.exists(token_file)) {
      stop("[download_from_gdc] Token file not found: ", token_file)
    }
    cmd <- paste(cmd, "-t", token_file)
  }


  ## --- Run command ---
  if (download) {
    message(sprintf("\nDownload: TRUE -> fetch all manifest entries."))
    status <- system(cmd, intern = verbose)

    ## --- Count final downloaded folders ---
    final_folders    <- list.dirs(file_dir, full.names = FALSE, recursive = FALSE)
    final_folders    <- basename(final_folders)
    n_final_folders  <- length(final_folders)

    message(sprintf("Final entrie(s) in directory: %-4d", n_final_folders))
  }


  ## --- Return ---
  return(file_map)
}