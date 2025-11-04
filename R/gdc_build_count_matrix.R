#' Build GDC RNA-seq count matrix
#'
#' Reads GDC STAR gene-count files and builds a numeric matrix with
#' samples as rows and genes as columns. Skips QC lines and extracts
#' a specified expression column (e.g. "unstranded", "tpm_unstranded").
#'
#' @param file_map Data frame with FILE_ID, FILE_NAME, SAMPLE_ID.
#' @param file_dir Directory containing GDC subfolders.
#' @param gene_id Column to use as gene identifier (default "gene_name").
#' @param gene_value Column with expression values (default "unstranded").
#' @param verbose Print progress messages.
#'
#' @return Numeric matrix (samples × genes).
#' @export
gdc_build_count_matrix <- function(
  file_map,
  file_dir,
  gene_id = "gene_name",
  gene_value = "unstranded",
  verbose = TRUE
){

  ## --- FILE and SAMPLE_ID ---
  file_ids   <- file_map$FILE_ID
  file_names <- file_map$FILE_NAME
  sample_ids <- file_map$SAMPLE_ID


  ## --- Read and extract ---
  value_list <- list()
  missing    <- character(0)
  gene_info  <- NULL

  for (i in seq_along(file_ids)){

    # Create path
    fpath <- file.path(file_dir, file_ids[i], file_names[i])
    sid   <- sample_ids[i]

    if (!file.exists(fpath)) {
      warning(sprintf("File not found: %s", fpath))
      missing <- c(missing, file_map$FILE_ID[i])
      next
    }

    # Read content
    if (verbose && (i %% 50 == 1 || i == nrow(file_map))) {
      pct <- round(100 * i / nrow(file_map), 1)
      message(sprintf("[gdc_build_count_matrix] Reading %-5d / %-5d: %s (%.1f%%)", i, nrow(file_map), sid, pct))
    }
    fdat <- read.delim(fpath, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, skip = 1)
    fdat <- fdat[grepl("^ENSG", fdat$gene_id), , drop = FALSE]

    # Check data format
    if (ncol(fdat) < 9) stop("Invalid format in ", fpath)
    if (!gene_value %in% names(fdat)) stop(sprintf("Column '%s' not found in %s", gene_value, fpath))

    # Store gene info (once) and appand value col
    if (is.null(gene_info)) gene_info <- fdat[, c("gene_id", "gene_name", "gene_type")]
    value_list[[sid]] <- fdat[[gene_value]]
  }

  ## --- Assemble matrix ---
  if (length(value_list) == 0) stop("No valid files were read successfully.")

  matr <- do.call(cbind, value_list)
  rownames(matr) <- gene_info[[gene_id]]
  colnames(matr) <- names(value_list)
  storage.mode(matr) <- "numeric"
  matr <- t(matr)


  ## --- Verbose message ---
  if (verbose){
    message("\nGDC count matrix summary:")
    message(sprintf("Dimensions: %d (samples) x %d (features)", nrow(matr), ncol(matr)))
    message(sprintf("Unique features: %d", length(unique(colnames(matr)))))
  }

  ## --- Return ---
  return(matr)
}