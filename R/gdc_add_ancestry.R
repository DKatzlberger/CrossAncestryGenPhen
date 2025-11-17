#' Add and compare GDC ancestry information
#'
#' Merges a GDC ancestry dataset with an existing metadata table (e.g., from cBioPortal)
#' by patient ID, compares missingness between ancestry sources, and retains
#' the GDC ancestry labels by default.
#'
#' @param file_map A data frame containing metadata from cBioPortal.
#' @param ancestry_data A data frame containing ancestry information from GDC (with columns \code{PATIENT_ID} and \code{GENETIC_ANCESTRY_LABEL})
#' @param verbose Logical; if \code{TRUE}, prints summary messages
#'
#' @return A merged data frame including ancestry information, prioritizing GDC ancestry values when duplicates occur.
#'
#' @export
#'
gdc_add_ancestry <- function(
  file_map,
  ancestry_data,
  verbose = TRUE
){

  ### PATIENT
  merged <- merge(
    file_map,
    ancestry_data,
    by = "PATIENT_ID",
    all.x = TRUE,
    suffixes = c(".cbio", ".gdc")
  )

  ## --- Count which has more NAs ---
  cbio_col <- "GENETIC_ANCESTRY_LABEL.cbio"
  gdc_col  <- "GENETIC_ANCESTRY_LABEL.gdc"

  if (cbio_col %in% names(merged) && gdc_col %in% names(merged)) {
    n_na_cbio <- sum(is.na(merged[[cbio_col]]))
    n_na_gdc  <- sum(is.na(merged[[gdc_col]]))
    total     <- nrow(merged)
    
    if (verbose) {
      message("\nAncestry coverage check")
      message(sprintf("cBioPortal ancestry missing: %d / %d (%.1f%%)", n_na_cbio, total, 100 * n_na_cbio / total))
      message(sprintf("GDCPortal ancestry missing:  %d / %d (%.1f%%)", n_na_gdc, total, 100 * n_na_gdc / total))
      
      if (n_na_gdc < n_na_cbio) {
        message("GDC ancestry has *better coverage* (fewer missing values).")
      } else if (n_na_gdc > n_na_cbio) {
        message("cBioPortal ancestry has *better coverage* (fewer missing values).")
      } else {
        message("Both sources have equal missing ancestry values.")
      }
    }
  }


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
  if (nrow(file_map) != nrow(merged))   message("[gdc_add_ancestry] Merging introduced new rows.")


  ## --- N unique ---
  n_IDs         <- length(unique(merged$ID))
  n_SAMPLE_IDs  <- length(unique(merged$SAMPLE_ID))
  n_PATIENT_IDs <- length(unique(merged$PATIENT_ID))


  ## --- Verbose message ---
  if (verbose){

    fmt_counts <- function(M, g_col) {
      if (is.null(M) || !nrow(M)) return("")
      tab <- table(M[[g_col]], dnn = NULL)
      paste(sprintf("%s: %-4d", names(tab), as.integer(tab)), collapse = "  ")
    }

    message("\nGDC add clinical data summary:")
    message(sprintf("%-9s  %-9s  %-9s  ->  %-9s", paste("IDs:", n_IDs), paste("SAMPLE IDs:", n_SAMPLE_IDs), paste("PATIENT IDs:", n_PATIENT_IDs), paste("output:", nrow(merged))))
    message(sprintf("Vials:    %-18s", fmt_counts(merged, "VIAL")))
    message(sprintf("Types:    %-18s", fmt_counts(merged, "SAMPLE_TYPE")))
    message(sprintf("Ancestry: %-18s", fmt_counts(merged, "GENETIC_ANCESTRY_LABEL")))
  }


  ## --- Return ---
  return(merged)

}
