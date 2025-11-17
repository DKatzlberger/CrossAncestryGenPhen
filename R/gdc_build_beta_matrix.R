#' Build GDC Beta-values count matrix
#'
#' Reads GDC Illumina 450k files and builds a numeric matrix with
#' samples as rows and features as columns.
#'
#' @param file_map Data frame with FILE_ID, FILE_NAME, SAMPLE_ID.
#' @param file_dir Directory containing GDC subfolders.
#' @param verbose Print progress messages.
#'
#' @return Numeric matrix (samples × features).
#' @export
gdc_build_beta_matrix <- function(
  file_map,
  file_dir,
  verbose = TRUE
){
  
  ## --- FILE and SAMPLE_ID ---
  file_ids   <- file_map$FILE_ID
  file_names <- file_map$FILE_NAME
  sample_ids <- file_map$SAMPLE_ID


  ## --- Initialize ---
  feature_info  <- NULL
  matr       <- NULL

  ## --- Loop over files ---
  for (i in seq_along(file_ids)) {

    fpath <- file.path(file_dir, file_ids[i], file_names[i])
    sid   <- sample_ids[i]

    if (!file.exists(fpath)) {
      warning(sprintf("File not found: %s", fpath))
      next
    }

    # Check data format
    fdat <- read.table(fpath, sep = "\t", colClasses = c("character", "numeric"))
    if (ncol(fdat) < 2) stop("Invalid format in ", fpath)

    # Initialize matrix and gene info from the first file
    if (is.null(feature_info)) {
      feature_info <- fdat[[1]]
      n_feat <- length(feature_info)
      n_samp <- nrow(file_map)
      # pre-allocate matrix (samples × features)
      matr <- matrix(NA_real_, nrow = n_samp, ncol = n_feat)
      colnames(matr) <- feature_info
      rownames(matr) <- sample_ids
    }

    # Expand matrix
    matr[sid, ] <- fdat[[2]]
    rm(fdat); gc(verbose = FALSE)

    # Verboe message
    if (verbose && (i %% 50 == 1 || i == nrow(file_map))) {
      pct    <- round(100 * i / nrow(file_map), 1)
      mem_gb <- as.numeric(object.size(matr)) / (1024^3)
      message(sprintf("[gdc_build_beta_matrix] Reading %-4d / %d: %s %-9s | approx. mem: %.1f GB", i, nrow(file_map), sid, sprintf("(%.1f%%)", pct), mem_gb))
    }
  }


  ## --- Check final matrix ---
  if (is.null(matr)) stop("No valid files read successfully.")


  ## --- Verbose message ---
  if (verbose) {
    message("\nGDC beta matrix summary:")
    message(sprintf("Dimensions: %d (samples) x %d (features)", nrow(matr), ncol(matr)))
    message(sprintf("Unique features: %d", length(unique(colnames(matr)))))
  }

  ## --- Return ---
  return(matr)
}