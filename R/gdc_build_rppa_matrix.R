#' Build GDC RPPA protein expression matrix (robust to inconsistent probe sets)
#'
#' @param file_map Data frame with FILE_ID, FILE_NAME, SAMPLE_ID
#' @param file_dir Directory with GDC subfolders
#' @param feature_id Column used as protein ID (default "peptide_target")
#' @param feature_value Column with expression values (default "protein_expression")
#' @param verbose Print progress messages
#'
#' @return Numeric matrix (samples × proteins)
#' 
#' @export
gdc_build_rppa_matrix <- function(
  file_map,
  file_dir,
  feature_id    = "peptide_target",
  feature_value = "protein_expression",
  verbose = TRUE
){

  file_ids   <- file_map$FILE_ID
  file_names <- file_map$FILE_NAME
  sample_ids <- file_map$SAMPLE_ID

  ## --- First pass: collect ALL features across all samples ---
  all_features <- character()

  if (verbose) message("[gdc_build_rppa_matrix] Collecting RPPA feature list...")

  for (i in seq_along(file_ids)) {
    fpath <- file.path(file_dir, file_ids[i], file_names[i])

    if (!file.exists(fpath)) next

    fdat <- read.delim(fpath, sep="\t", header=TRUE, stringsAsFactors=FALSE)

    if (!(feature_id %in% colnames(fdat)))
      stop(sprintf("Column '%s' missing in file %s", feature_id, fpath))

    all_features <- union(all_features, fdat[[feature_id]])
  }

  all_features <- sort(all_features)
  n_feat <- length(all_features)
  n_samp <- length(sample_ids)

  if (verbose) {
    message(sprintf("[gdc_build_rppa_matrix] Total unique RPPA features found: %d\n", n_feat))
  }

  ## --- Allocate matrix ---
  matr <- matrix(
    NA_real_,
    nrow = n_samp,
    ncol = n_feat,
    dimnames = list(sample_ids, all_features)
  )

  ## --- Fill matrix sample by sample ---
  for (i in seq_along(file_ids)) {

    fpath <- file.path(file_dir, file_ids[i], file_names[i])
    sid <- sample_ids[i]

    if (!file.exists(fpath)) {
      warning(sprintf("[gdc_build_rppa_matrix] Missing file: %s", fpath))
      next
    }

    fdat <- read.delim(
      fpath, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE
    )

    if (!all(c(feature_id, feature_value) %in% colnames(fdat))) {
      stop(sprintf("[gdc_build_rppa_matrix] Missing RPPA columns in %s", fpath))
    }

    # Map features to matrix positions
    idx <- match(fdat[[feature_id]], all_features)

    matr[i, idx] <- fdat[[feature_value]]

    if (verbose && (i %% 50 == 1 || i == n_samp)) {
      pct <- round(100 * i / n_samp, 1)
      message(sprintf("[gdc_build_rppa_matrix] %d / %d (%.1f%%): %s", i, n_samp, pct, sid))
    }
  }

  ## --- Verbose message ---
  if (verbose) {
    message("\nGDC RPPA matrix summary:")
    message(sprintf("Dimensions: %d (samples) x %d (features)", nrow(matr), ncol(matr)))
    message(sprintf("Unique features: %d", length(unique(colnames(matr)))))
  }

  return(matr)
}