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
  matr          <- NULL

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


  ## --- QC step ---
  # ann450k <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19)

  # # Remove feature with any NA
  # na_mask <- colSums(is.na(matr)) == 0
  # matr    <- matr[, na_mask, drop = FALSE]

  # # Remove probes that match chromosomes X and Y
  # sex_mask <- !(colnames(matr) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
  # matr     <- matr[, sex_mask, drop = FALSE]
  # rm(sex_mask); gc(verbose = FALSE)

  # # Remove SNPs overlapped probe
  # no_snp_probe <- ann450k$Name[is.na(ann450k$Probe_rs)]
  # snp_probe    <- ann450k[!is.na(ann450k$Probe_rs), ]
  # # Snps with maf <= 0.05
  # snp5_probe   <- snp_probe$Name[snp_probe$Probe_maf <= 0.05]
  # matr         <- matr[colnames(matr) %in% c(no_snp_probe, snp5_probe), ]
  # rm(no_snp_probe, probe, probe_na, snp_probe, snp5_probe); gc(verbose = FALSE)

  # # Removing probes that have been demonstrated to map to multiple places in the genome
  # data("crossReactiveProbes", package = "IlluminaHumanMethylation450kanno.ilmn12.hg19", envir = environment())
  # cr_mask <- colnames(matr) %in% crossReactiveProbes
  # matr    <- matr[, !cr_mask, drop = FALSE]


  ## --- Verbose message ---
  if (verbose) {
    message("\nGDC beta matrix summary:")
    message(sprintf("Dimensions: %d (samples) x %d (features)", nrow(matr), ncol(matr)))
    message(sprintf("Unique features: %d", length(unique(colnames(matr)))))
  }

  ## --- Return ---
  return(matr)
}