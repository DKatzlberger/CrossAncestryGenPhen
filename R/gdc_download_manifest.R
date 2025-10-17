#' Download data from GDC via gdc-client
#'
#' Reads a manifest file and downloads all files
#' Uses current dir if out_dir is not specified
#' Supports optional token for protected data
#' Prints file and download counts when verbose
#'
#' @param manifest_file Path to GDC manifest file
#' @param toke_file Optional path to token file
#' @param out_dir Output dir, defaults to getwd()
#' @param client_path Path to gdc-client binary
#' @param verbose Print progress messages
#'
#' @return Invisibly returns system status
#' 
#' @export
gdc_download_manifest <- function(
  manifest_file,
  toke_file = NULL,
  out_dir = NULL,
  client_path,
  verbose = TRUE
){


  ## --- Check for the manifest file ---
  if (!file.exists(manifest_file)) {
    stop("[download_from_gdc] Manifest file not found: ", manifest_file)
  }

  ## --- Set output directory ---
  if (is.null(out_dir)) {
    out_dir <- getwd()
  } else if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }


  ## --- Count files in manifest ---
  manifest_data <- tryCatch(
    read.delim(manifest_file, header = TRUE, stringsAsFactors = FALSE),
    error = function(e) stop("[download_from_gdc] Failed to read manifest file.")
  )
  n_manifest <- nrow(manifest_data)


  ## --- Base command ---
  cmd <- sprintf('%s download -m "%s" -d "%s"', client_path, manifest_file, out_dir)


  ## --- Token (for protected data) ---
  if (!is.null(token_file)) {
    if (!file.exists(token_file)) {
      stop("[download_from_gdc] Token file not found: ", token_file)
    }
    cmd <- paste(cmd, "-t", token_file)
  }


  ## --- Run command ---
  status <- system(cmd, intern = verbose)


  ## --- Count downloaded files ---
  downloaded_folders <- list.dirs(out_dir, full.names = FALSE, recursive = FALSE)
  n_downloaded <- length(downloaded_folders)


  ## --- Verbose output ---
  if (verbose) {
    message("\nGDC dowload summary:")
    message(sprintf("Manifest entries: %-4d",    n_manifest))
    message(sprintf("Manifest downloaded: %-4d", n_downloaded))
  }


  ## --- Return ---
  return(invisible(status))
}