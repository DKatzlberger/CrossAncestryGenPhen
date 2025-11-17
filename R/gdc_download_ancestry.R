#' Download and process TCGA ancestry assignments
#'
#' Retrieves the TCGA genetic ancestry assignment data from
#' Carrot-Zhang et al. (2020) and extracts the patient ID and
#' consensus ancestry label.
#'
#' @param url Character string; URL of the Excel file to download.
#' @param sheet Character string; name of the sheet containing the ancestry data.
#' @param verbose Logical; if \code{TRUE}, prints summary messages.
#'
#' @return A data frame with two columns: \code{PATIENT_ID} and
#'   \code{GENETIC_ANCESTRY_LABEL}, where ancestry labels are in uppercase.
#'
#' @importFrom utils download.file
#' @importFrom openxlsx read.xlsx
#' @export
#'
gdc_download_ancestry <- function(
  url = "https://ars.els-cdn.com/content/image/1-s2.0-S1535610820302117-mmc2.xlsx",
  sheet = "S1 Calls per Patient",
  verbose = TRUE
){


  ## --- Create tmp file ---
  tmp <- tempfile(fileext = ".xlsx")
  download.file(url, tmp, mode = "wb", quiet = TRUE)


  ## --- Extract sheet ---
  raw_sheet <- openxlsx::read.xlsx(tmp, sheet = sheet)
  headers   <- as.character(unlist(raw_sheet[1, ]))
  raw_sheet <- raw_sheet[-1, ]
  names(raw_sheet) <- headers

  ## --- Rename ---
  ancestry_df <- raw_sheet[, c("patient", "consensus_ancestry")]
  names(ancestry_df) <- c("PATIENT_ID", "GENETIC_ANCESTRY_LABEL")
  ancestry_df$GENETIC_ANCESTRY_LABEL <- toupper(trimws(ancestry_df$GENETIC_ANCESTRY_LABEL))
  


  ## --- Verbose message ---
  if (verbose) {
    message("\nGDC ancestry call summary:")
    message("Loaded ancestry data with ", nrow(ancestry_df), " patients.")
    message("Columns: ", paste(names(ancestry_df), collapse = ", "))
  }

  ## --- Return ---
  return(ancestry_df)
}
