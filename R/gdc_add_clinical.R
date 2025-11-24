#' Merge GDC manifest with clinical data
#'
#' Combines file map from GDC with cBioPortal
#' sample and patient metadata. Keeps one row
#' per SAMPLE_ID after optional filtering by vial/plate.
#'
#' @param file_map Data frame from GDC manifest
#' @param clinical_data List with SAMPLE, PATIENT attributes
#' @param vial Optional vial code to keep (default NULL = no filtering)
#' @param plate Optional plate filter, e.g. "first" or NULL for no filtering
#' @param verbose Print progress messages
#'
#' @return Data frame merged on SAMPLE_ID,
#'   containing file, vial, plate, and clinical
#'   annotations. Duplicate columns removed and
#'   SAMPLE_ID made unique if filtering applied.
#'
#' @export
gdc_add_clinical <- function(
  file_map,
  clinical_data,
  vial = NULL,
  plate = NULL,
  verbose = TRUE
){

  ### SAMPLE
  merged <- merge(
    file_map,
    clinical_data$SAMPLE,
    by = "SAMPLE_ID",
    all.x = TRUE,
    suffixes = c(".gdc", ".cbio")
  )

  ## --- Remove duplicated columns (keep .gdc versions) ---
  prefixes <- sub("\\.gdc$", "", grep("\\.gdc$", names(merged), value = TRUE))
  keep_cols <- names(merged)[!names(merged) %in% paste0(prefixes, ".cbio")]
  merged <- merged[, keep_cols]
  names(merged) <- sub("\\.gdc$", "", names(merged))

  ### PATIENT
  merged <- merge(
    merged,
    clinical_data$PATIENT,
    by = "PATIENT_ID",
    all.x = TRUE,
    suffixes = c(".gdc", ".cbio")
  )

  ## --- Remove duplicated columns (keep .gdc versions) ---
  prefixes <- sub("\\.gdc$", "", grep("\\.gdc$", names(merged), value = TRUE))
  keep_cols <- names(merged)[!names(merged) %in% paste0(prefixes, ".cbio")]
  merged <- merged[, keep_cols]
  names(merged) <- sub("\\.gdc$", "", names(merged))


  ## --- Set order ---
  ordered_cols <- c(
    "RELEASE", "FILE_ID", "FILE_NAME", 
    "ID", "SAMPLE_ID", "PATIENT_ID", 
    "VIAL", "PLATE", "SAMPLE_TYPE", 
    "CANCER_TYPE", "CANCER_TYPE_DETAILED",
    "SEX", "AGE", "SUBTYPE", "GENETIC_ANCESTRY_LABEL"
  )
  ordered_cols <- ordered_cols[ordered_cols %in% names(merged)]
  remaining_cols <- setdiff(names(merged), ordered_cols)
  merged <- merged[, c(ordered_cols, remaining_cols), drop = FALSE]


  ## --- Validate merge ---
  if (nrow(file_map) != nrow(merged)) message("[gdc_add_clinical] Merging introduced new rows.")
  if (anyDuplicated(merged$ID))       message("[gdc_add_clinical] ID is not unique anymore.")


  ## --- N unique ---
  n_IDs         <- length(unique(merged$ID))
  n_SAMPLE_IDs  <- length(unique(merged$SAMPLE_ID))
  n_PATIENT_IDs <- length(unique(merged$PATIENT_ID))


  ## --- Optional filtering by vial ---
  if (!is.null(vial)) {
    vial_arg <- match.arg(vial, choices = c("first", LETTERS))  # support both
    if (vial_arg == "first") {
      merged <- merged[order(merged$SAMPLE_ID, merged$VIAL), , drop = FALSE]
      merged <- merged[!duplicated(merged$SAMPLE_ID), , drop = FALSE]
    } else {
      merged <- merged[merged$VIAL == vial_arg, , drop = FALSE]
    }
  }


  ## --- Optional filtering by plate ---
  if (!is.null(plate)) {
    plate_arg <- match.arg(plate, choices = c("first"))
    if (plate_arg == "first") {
      merged <- merged[order(merged$SAMPLE_ID, merged$PLATE), , drop = FALSE]
      merged <- merged[!duplicated(merged$SAMPLE_ID), , drop = FALSE]
    }
  }


  ## --- Check unique SAMPLE_ID ---
  dup_ids <- merged$SAMPLE_ID[duplicated(merged$SAMPLE_ID)]
  n_dup <- length(dup_ids)


  ## --- Verbose message ---
  if (verbose){

    fmt_counts <- function(M, g_col) {
      if (is.null(M) || !nrow(M)) return("")
      tab <- table(M[[g_col]], dnn = NULL)
      paste(sprintf("%s: %-4d", names(tab), as.integer(tab)), collapse = "  ")
    }

    message("\nGDC add clinical data summary:")
    message(sprintf("%-9s  %-9s  %-9s  ->  %-9s", paste("IDs:", n_IDs), paste("SAMPLE IDs:", n_SAMPLE_IDs), paste("PATIENT IDs:", n_PATIENT_IDs), paste("output:", nrow(merged))))
    message(sprintf("Vials:   %-18s", fmt_counts(merged, "VIAL")))
    message(sprintf("Types:   %-18s", fmt_counts(merged, "SAMPLE_TYPE")))
  }


  ## --- Print duplication warning ---
  if (n_dup > 0) {
    message(sprintf("Warning: Found %d duplicated SAMPLE_ID(s)!", n_dup))
    dup_unique <- unique(dup_ids)
    dup_df     <- merged[merged$SAMPLE_ID %in% dup_unique, , drop = FALSE]

    # Identify which columns differ among duplicates
    diff_cols <- sapply(names(dup_df), function(col) {any(tapply(dup_df[[col]], dup_df$SAMPLE_ID, function(x) length(unique(x)) > 1))})
    diff_cols <- names(diff_cols[diff_cols])
    diff_cols <- setdiff(diff_cols, c("RELEASE", "FILE_ID", "FILE_NAME", "ID"))

    message(sprintf("Columns with differing values among duplicates: %s", paste(diff_cols, collapse = ", ")))
  }
  

  ## --- Return ---
  return(merged)
}
