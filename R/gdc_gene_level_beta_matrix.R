#' Collapse gene-level probes with mean and remove NAs from beta matrix
#'
#' @param matr Numeric matrix (samples × features).
#' @param verbose Logical; show summary (default TRUE).
#'
#' @return Matrix with NAs removed and duplicate genes summed.
#' 
#' @importFrom minfi getAnnotation
#' @import IlluminaHumanMethylation450kanno.ilmn12.hg19
#' @export
gdc_gene_level_beta_matrix <- function(
  matr,
  verbose = TRUE
){

  ## --- Check for matrix ---
  if (!is.matrix(matr)) stop("Input must be a matrix (samples × features).")


  ## --- Record start dimensions ---
  n_start_samples  <- nrow(matr)
  n_start_features <- ncol(matr)


  ## --- Remove feature with any NA ---
  na_mask <- colSums(is.na(matr)) == 0
  matr    <- matr[, na_mask, drop = FALSE]
  n_removed_na <- ncol(matr)


  ## --- Create probe2gene mapping ---
  ann450k <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19)
  probe2gene <- data.frame(
    probe_id = ann450k$Name,
    gene_id  = ann450k$UCSC_RefGene_Name,
    stringsAsFactors = FALSE
  )

  # Remove probes with no gene mapping
  probe2gene <- probe2gene[probe2gene$gene_id != "", , drop = FALSE]
  gene_split <- strsplit(probe2gene$gene_id, ";")

  probe2gene <- data.frame(
    probe_id = rep(probe2gene$probe_id, lengths(gene_split)),
    gene_id  = unlist(gene_split),
    stringsAsFactors = FALSE
  )


  ## --- Filter mapping to probes in matrix ---
  common <- intersect(colnames(matr), probe2gene$probe_id)
  matr <- matr[, common, drop = FALSE]
  probe2gene <- probe2gene[match(common, probe2gene$probe_id), , drop = FALSE]


  ## --- Aggregate by mean beta per gene ---
  feature_names <- probe2gene$gene_id
  agg_matr <- t(rowsum(t(matr), group = feature_names)) / as.numeric(table(feature_names))


  ## --- Verbose message ---
  if (verbose){
    message("\nGDC gene-level count matrix summary:")
    message(sprintf("Features:   %d -> %d (no NAs) -> %d (mean)", n_start_features, n_removed_na, ncol(agg_matr)))
    message(sprintf("Dimensions: %d (samples) x %d (features)", nrow(agg_matr), ncol(agg_matr)))
  }


  ## --- Return ---
  return(agg_matr)
}
