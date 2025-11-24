#' Download 450k probe map from cBioPortal
#'
#' Downloads a methylation archive from cBioPortal,
#' extracts the probe-to-gene mapping, summarizes probe
#' coverage, removes unmapped probes, and returns a
#' clean probe–gene table.
#'
#' @param study_id cBioPortal study identifier
#' @param profile_id Methylation profile filename
#' @param base_url Base URL for data download
#' @param verbose Print progress messages
#'
#' @return Data frame with probe–gene mapping
#'
#' @importFrom data.table fread
#' @export
#'
cbioportal_download_450kprobes <- function(
  study_id,
  profile_id = "data_methylation_hm27_hm450_merged.txt",
  base_url = "https://cbioportal-datahub.s3.amazonaws.com",
  verbose = TRUE
){

  ## --- Download ---
  options(timeout = 1000000)
  url  <- sprintf("%s/%s.tar.gz", base_url, study_id)
  tar  <- paste0(study_id, ".tar.gz")
  if (!file.exists(tar)) download.file(url, tar, mode = "wb",  method = "curl")


  ## --- Extract file ---
  file <- file.path(study_id, profile_id)
  if (!file %in% untar(tar, list = TRUE)) stop("File NOT found in the tar archive")
  untar(tar, files = file)


  ## --- Probe -> gene mapping ---
  data <- fread(file)
  probe_map <- data.frame(
    probe = data$ENTITY_STABLE_ID,
    gene  = data$NAME,
    stringsAsFactors = FALSE
  )

  ## --- Summary counts ---
  total_probes  <- nrow(probe_map)
  mapped_probes <- sum(!is.na(probe_map$gene) & probe_map$gene != "")
  pct_mapped    <- round(mapped_probes / total_probes * 100, 2)


  ## --- Verbose message ---
  if (verbose){
    message("\nProbe summary:")
    message("Probes returned from cBioPortal: ", total_probes)
    message("Probes with geneSymbol: ", mapped_probes, " (", pct_mapped, "%)")
  }


  ## --- Remove empty geneSymbols ---
  probe_map <- probe_map[probe_map$gene != "" & !is.na(probe_map$gene), ]


  ## --- Remove links ---
  unlink(tar); unlink("brca_tcga_pan_can_atlas_2018", recursive = TRUE)


  ## --- Return ---
  return(probe_map)
}