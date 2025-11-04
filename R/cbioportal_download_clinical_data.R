#' Download and Combine Clinical Data from cBioPortal
#'
#' This function retrieves and merges patient- and sample-level clinical data 
#' for a given study from the [cBioPortal API](https://www.cbioportal.org/api). 
#' It fetches both *SAMPLE* and *PATIENT* clinical attributes, converts them 
#' into a wide format, and merges them by `PATIENT_ID`.
#'
#' @param study_id The cBioPortal study ID (e.g., `"brca_tcga"`).
#' @param base_url The base URL of the cBioPortal API.  Defaults to `"https://www.cbioportal.org/api"`.
#' @param verbose Whether to print progress and summary messages. Defaults to `TRUE`.
#'
#' @section API Reference:
#' - Endpoint: `/studies/{studyId}/clinical-data`
#' - Projection used: `DETAILED`
#' - Example: `https://www.cbioportal.org/api/studies/brca_tcga/clinical-data?projection=DETAILED&clinicalDataType=PATIENT`
#'
#' @import httr
#' @import jsonlite
#' @export
cbioportal_download_clinical_data <- function(
  study_id,
  base_url = "https://www.cbioportal.org/api",
  verbose = TRUE
){

  ## --- Helper: Fetches one type of attribute ---
  fetch_clinical_type <- function(
    type
  ){

    ## --- Define keys ---
    if (type == "SAMPLE"){
      fetch_key <- c("sampleId")
      fetch_col <- c("sampleId", "patientId", "clinicalAttributeId", "value")
    } else if (type == "PATIENT") {
      fetch_key  <- c("patientId")
      fetch_col  <- c("patientId", "clinicalAttributeId", "value")
    }


    ## --- Build API endpoint ---
    url <- sprintf("%s/studies/%s/clinical-data?projection=DETAILED&clinicalDataType=%s", base_url, study_id, type)
    res <- httr::GET(url, httr::add_headers(Accept = "application/json"))


    ## --- Error handling --- 
    if (httr::status_code(res) != 200) {
      stop("Request failed: ", httr::status_code(res), " - ", httr::content(res, "text"))
    }

    ## --- Parse output ---
    json_txt <- rawToChar(httr::content(res, "raw"))
    json_res <- jsonlite::fromJSON(json_txt)


    ## --- Check if keys exist ---
    missing_cols <- setdiff(fetch_col, colnames(json_res))
    if (length(missing_cols) > 0) {
      stop(sprintf("[cbioportal_download_clinical_data] Missing expected columns '%s' for %s attributes.", paste(missing_cols, collapse = ", "), type))
    }


    ## --- Format output to wide ---
    res <- json_res[, fetch_col]
    key <- unique(json_res[, fetch_key])

    if (is.null(key)){
      stop("[cbioportal_download_clinical_data] Specified ID not found! There is no fallback...")
    }

    # Create wide data.frame
    if (type ==  "SAMPLE"){
      wide <- data.frame(SAMPLE_ID = key, stringsAsFactors = FALSE)
      # Add a patient_id as second key
      map  <- unique(json_res[, c("sampleId", "patientId")])
      wide <- merge(wide, map, by.x = "SAMPLE_ID", by.y = "sampleId", all.y = TRUE)
      if ("patientId" %in% names(wide)) {
        names(wide)[names(wide) == "patientId"] <- "PATIENT_ID"
      }
    } else if (type == "PATIENT") {
      wide <- data.frame(PATIENT_ID = key, stringsAsFactors = FALSE)
    }

    # Add attributes
    attributes <- unique(json_res$clinicalAttributeId)
    for (attr in attributes){
      attribute_rows <- json_res[json_res$clinicalAttributeId == attr, ]
      # Choose the right matching key
      if (type == "SAMPLE") {
        match_idx <- match(wide$SAMPLE_ID, attribute_rows$sampleId)
      } else if (type == "PATIENT") {
        match_idx <- match(wide$PATIENT_ID, attribute_rows$patientId)
      }
      # Add matched attribute
      wide[[attr]] <- attribute_rows$value[match_idx]
    }


    ## --- Return ---
    return(wide)
  }

  ## --- Call helper function ---
  sample_attr  <- fetch_clinical_type(type = "SAMPLE")
  patient_attr <- fetch_clinical_type(type = "PATIENT")


  ## --- Add study id ---
  sample_attr$STUDY_ID <- study_id
  patient_attr$STUDY_ID <- study_id


  ## --- Verbose message ---
  if (verbose){
    message("\ncBioportal clinical data summary:")
    message(sprintf("Study: %-10s", study_id))

    message(sprintf(
      "Sample-level   N: %-5d  features: %-5d",
      nrow(sample_attr),
      length(setdiff(colnames(sample_attr), c("PATIENT_ID", "SAMPLE_ID")))
    ))

    message(sprintf(
      "Patient-level  N: %-5d  features: %-5d",
      nrow(patient_attr),
      length(setdiff(colnames(patient_attr), c("PATIENT_ID", "SAMPLE_ID")))
    ))
  }

  ## --- Return ---
  return(
    list(
      SAMPLE  = sample_attr,
      PATIENT = patient_attr
    )
  )
}
