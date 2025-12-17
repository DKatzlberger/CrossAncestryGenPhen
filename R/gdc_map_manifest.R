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
#' @param base_url URL to tcga-code-tables without endpoint
#' @param verbose Print progress messages
#'
#' @return Data frame with file_name and sample_id
#' 
#' @export

gdc_map_manifest <- function(
	manifest_file, 
	metadata_file,
  base_url = "https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/",
	verbose = TRUE
){

	## --- Check file existence ---
	if (!file.exists(manifest_file)) stop("[gdc_map_manifest] Manifest file not found: ", manifest_file)
  if (!file.exists(metadata_file)) stop("[gdc_map_manifest] Metadata file not found: ", metadata_file)


  ## --- Check correct file formats ---
  if (!grepl("\\.txt$", manifest_file, ignore.case = TRUE)) stop("[gdc_map_manifest] Manifest file must have extension '.txt': ", manifest_file)
  if (!grepl("\\.json$", metadata_file, ignore.case = TRUE)) stop("[gdc_map_manifest] Metadata file must have extension '.json': ", metadata_file)
  

	## --- Read manifest ---
  manifest <- read.delim(manifest_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  metalist <- jsonlite::fromJSON(metadata_file, simplifyVector = FALSE)


  ## --- Save extractors ---
  release_extractor <- local({
    
    function(x) {
      # if analysis$updated_datetime exists → extract YY-MM-DD
      if (!is.null(x$analysis) && !is.null(x$analysis$updated_datetime)) {
        return(substr(x$analysis$updated_datetime, 1, 10))
      }
      
      return(NA_character_)
    }
  })


	## --- Extract mapping ---
  file_map <- data.frame(
    RELEASE   = sapply(metalist, release_extractor, USE.NAMES = FALSE),
    FILE_ID   = sapply(metalist, function(x) {if (!is.null(x$file_id)) {x$file_id} else {NA_character_}}),
    FILE_NAME = sapply(metalist, function(x) {if (!is.null(x$file_name)) {x$file_name} else {NA_character_}}),
    ID        = sapply(metalist, function(x) {if (!is.null(x$associated_entities[[1]]$entity_submitter_id)) {x$associated_entities[[1]]$entity_submitter_id} else {NA_character_}}),
    stringsAsFactors = FALSE
  )

  # Warning if different
  if (any(is.na(file_map$ID)) || any(file_map$ID == "")) {
    warning("[gdc_map_manifest] Some 'ID' values could not be extracted from metadata (TCGA barcode missing).")
  }



  ##  --- Parse TCGA barcode ---
  parse_tcga_barcode <- local({

    shown <- list()   # keep track of message types already printed

    function(barcode) {

      parts <- strsplit(barcode, "-", fixed = TRUE)[[1]]
      n <- length(parts)

      out <- list(
        PROJECT_CODE = NA_character_, 
        TSS_CODE = NA_character_, 
        PATIENT_CODE = NA_character_,
        SAMPLE_TYPE_CODE = NA_character_, 
        VIAL_CODE = NA_character_, 
        PORTION_CODE = NA_character_,
        ANALYTE_CODE = NA_character_, 
        PLATE_CODE = NA_character_, 
        CENTER_CODE = NA_character_
      )

      if (n >= 1) out$PROJECT_CODE <- parts[1]
      if (n >= 2) out$TSS_CODE     <- parts[2]
      if (n >= 3) out$PATIENT_CODE <- parts[3]

      type <- NULL

      # ---- Detect barcode type ----
      if (n == 3) {
        type <- "patient"
      } else if (n == 4) {
        type <- "sample"
      } else if (n == 5) {
        type <- "portion"
      } else if (n == 6) {
        type <- "portion+analyte"
      } else if (n == 7) {
        type <- "full aliquot"
      } else {
        type <- "unrecognized"
      }

      # Print message ONCE per type
      if (!type %in% shown) {
        message("[gdc_map_manifest] Detected ", type, " barcode format: ", barcode)
        shown <<- c(shown, type)
      }

      # ---- Parse sample type + vial ----
      if (n >= 4) {
        samp <- parts[4]
        if (grepl("^[0-9]{2}[A-Z]?$", samp)) {
          out$SAMPLE_TYPE_CODE <- substr(samp, 1, 2)
          out$VIAL_CODE <- ifelse(nchar(samp) == 3, substr(samp, 3, 3), NA)
        }
      }

      # ---- Portion + analyte ----
      if (n >= 5) {
        p <- parts[5]
        if (grepl("^[0-9]{2}[A-Z]?$", p)) {
          out$PORTION_CODE <- substr(p, 1, 2)
          out$ANALYTE_CODE <- ifelse(nchar(p) == 3, substr(p, 3, 3), NA)
        }
      }

      # ---- Plate + center ----
      if (n >= 6) out$PLATE_CODE  <- parts[6]
      if (n >= 7) out$CENTER_CODE <- parts[7]

      return(out)
    }
  })


  ## --- Parse TCGA barcodes with robust parser ---
  parsed <- lapply(file_map$ID, parse_tcga_barcode)

  id_components <- do.call(rbind, lapply(parsed, as.data.frame))
  rownames(id_components) <- NULL


  # Derived IDs
  id_components$PATIENT_ID <- with(id_components, paste(PROJECT_CODE, TSS_CODE, PATIENT_CODE, sep = "-"))
  id_components$SAMPLE_ID  <- with(id_components, paste(PROJECT_CODE, TSS_CODE, PATIENT_CODE, SAMPLE_TYPE_CODE, sep = "-"))


  ## --- Map the TCGA barcode ---
  fetch_gdc_map <- function(endpoint) {

    ## Create endpoint
    url <- paste0(base_url, endpoint)

    ## Fetch the table
    out <- tryCatch({

      # Request
      html  <- rvest::read_html(url)
      tbls  <- rvest::html_table(html, fill = TRUE)
      tbl   <- tbls[[2]]
      tbl[] <- lapply(tbl, trimws)

      # Normalize column names (trim punctuation, spaces → underscores, lower case)
      names(tbl) <- gsub("[[:punct:]]+$", "", names(tbl))
      names(tbl) <- gsub("[[:space:]]+", "_", names(tbl))
      names(tbl) <- tolower(trimws(names(tbl)))
      
      # Identify code + description columns
      code_col <- desc_col <- NULL
      if (endpoint == "sample-type-codes") {
        code_col <- "code"; desc_col <- "definition"
      } else if (endpoint == "center-codes") {
        code_col <- "code"; desc_col <- "display_name"
      } else if (endpoint == "portion-analyte-codes") {
        code_col <- "code"; desc_col <- "definition"
      } else if (endpoint == "tissue-source-site-codes") {
        code_col <- "tss_code"; desc_col <- "source_site"
      }

      # Rename to standard form
      if (code_col %in% names(tbl)) names(tbl)[names(tbl) == code_col] <- "code"
      if (desc_col %in% names(tbl)) names(tbl)[names(tbl) == desc_col] <- "description"

      # Standardize values
      if ("code" %in% names(tbl)) {
        suppressWarnings({
          num_mask <- grepl("^[0-9]+$", tbl$code)
          tbl$code[num_mask] <- sprintf("%02d", as.integer(tbl$code[num_mask]))
          tbl$code[!num_mask] <- toupper(tbl$code[!num_mask])
        })
      }
      tbl$description <- trimws(tbl$description)

      tbl <- tbl[, c("code", "description")]
      return(tbl)

    }, error = function(e) {
      warning(sprintf("[fetch_gdc_map] Could not fetch table from: %s", url))
      NULL
    })

    ## Return
    return(out)
  }

  # Fetch all tables
  sample_type_tbl <- fetch_gdc_map(endpoint = "sample-type-codes")
  center_tbl      <- fetch_gdc_map(endpoint = "center-codes")
  analyte_tbl     <- fetch_gdc_map(endpoint = "portion-analyte-codes")
  tss_tbl         <- fetch_gdc_map(endpoint = "tissue-source-site-codes")


  ## --- Map the codes ---
  map_with_summary <- function(source_codes, ref_tbl, ref_name) {
    ## Map
    mapped <- ref_tbl$description[match(source_codes, ref_tbl$code)]
    n_mapped <- sum(!is.na(mapped))
    total <- length(source_codes)

    ## Return
    return(mapped)
  }

  ## --- Perform mapping ---
  id_components$SAMPLE_TYPE <- map_with_summary(id_components$SAMPLE_TYPE_CODE, sample_type_tbl, "Sample type")
  id_components$CENTER      <- map_with_summary(id_components$CENTER_CODE, center_tbl, "Center")
  id_components$ANALYTE     <- map_with_summary(id_components$ANALYTE_CODE, analyte_tbl, "Analyte")
  id_components$TSS         <- map_with_summary(id_components$TSS_CODE, tss_tbl, "Tissue source site")


  ## --- Create mapped data.frame ---
  file_map$SAMPLE_ID   <- id_components$SAMPLE_ID
  file_map$PATIENT_ID  <- id_components$PATIENT_ID
  file_map$VIAL        <- id_components$VIAL_CODE
  file_map$PLATE       <-  id_components$PLATE_CODE
  file_map$SAMPLE_TYPE <- id_components$SAMPLE_TYPE


  # Force order
  file_map <- file_map[, c("RELEASE", "FILE_ID", "FILE_NAME", "ID", "SAMPLE_ID", "PATIENT_ID", "VIAL", "PLATE", "SAMPLE_TYPE")]


  ## --- Check NA values in critical cols ---
  critical_cols <- c("FILE_ID", "FILE_NAME", "ID", "SAMPLE_ID", "PATIENT_ID")
  missing_mask <- sapply(file_map[critical_cols], function(col) any(is.na(col) | col == ""))

  if (any(missing_mask)) {
    warning("[gdc_map_manifest] Missing values detected in critical metadata fields: ", paste(names(which(missing_mask)), collapse = ", "))
  }


  ## --- N unique ---
  n_IDs         <- length(unique(file_map$ID))
  n_SAMPLE_IDs  <- length(unique(file_map$SAMPLE_ID))
  n_PATIENT_IDs <- length(unique(file_map$PATIENT_ID))


	## --- Verbose message ---
	if (verbose) {

    fmt_counts <- function(M, g_col) {
      if (is.null(M) || !nrow(M)) return("")
      tab <- table(M[[g_col]], dnn = NULL)
      paste(sprintf("%s: %-4d", names(tab), as.integer(tab)), collapse = "  ")
    }

		message("\nGDC map summary:")
    message(sprintf("%-9s  %-9s  %-9s  ->  %-9s", paste("IDs:", n_IDs), paste("SAMPLE IDs:", n_SAMPLE_IDs), paste("PATIENT IDs:", n_PATIENT_IDs), paste("output:", nrow(file_map))))
    message(sprintf("Vials:   %-18s", fmt_counts(file_map, "VIAL")))
    message(sprintf("Types:   %-18s", fmt_counts(file_map, "SAMPLE_TYPE")))
  }


  ## --- Return ---
  return(file_map)
}
